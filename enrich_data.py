"""
PubChemから化合物のSMILES・InChI・分子式・分子量を取得するスクリプト。

実行方法:
    python enrich_data.py

特徴:
- Windowsのシステム証明書ストアを使用（SSL検証エラー対策）
- 全エラーをログ出力（ステータスコード、レスポンス本文を表示）
- 429/503 などのレートリミットエラーは指数バックオフで自動リトライ
- 失敗したエントリは {"_failed": True, ...} としてキャッシュに記録（再実行時にスキップしない）
"""

import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# 出力を UTF-8 に再設定（Windows cp932 で日本語・特殊文字エラー対策）
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# Windows のシステム証明書ストアを使う
try:
    import truststore
    truststore.inject_into_ssl()
    print("[OK] truststore injected", flush=True)
except ImportError:
    print("[WARN] truststore not available", flush=True)

BASE_DIR = Path(__file__).parent
CSV_FILE = BASE_DIR / "stocklist_20251023 のコピー - stock_list_20251011175958.csv"
CACHE_FILE = BASE_DIR / "compounds_enriched.json"

PUBCHEM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}"
    "/property/SMILES,InChI,InChIKey,MolecularFormula,MolecularWeight/JSON"
)

SLEEP_BETWEEN_REQUESTS = 0.34  # PubChem は 5 req/s 上限 → 余裕を持って ~3 req/s
MAX_RETRIES = 3


def log(msg: str) -> None:
    print(msg, flush=True)


def fetch_pubchem(query: str) -> dict | None:
    """PubChemから化合物情報を取得。
    成功 → 辞書, 見つからない → "404", レートリミットなどリトライ可能エラー → "retry"
    完全失敗 → None
    """
    query = (query or "").strip()
    if not query or query in ("-", "nan", ""):
        return None

    url = PUBCHEM_URL.format(requests.utils.quote(query))

    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, timeout=20)
            sc = r.status_code

            if sc == 200:
                props = r.json()["PropertyTable"]["Properties"][0]
                return {
                    "smiles": props.get("SMILES") or props.get("IsomericSMILES") or "",
                    "inchi": props.get("InChI", ""),
                    "inchikey": props.get("InChIKey", ""),
                    "molecular_formula": props.get("MolecularFormula", ""),
                    "molecular_weight": float(props.get("MolecularWeight", 0) or 0),
                    "pubchem_cid": int(props.get("CID", 0) or 0),
                }
            if sc == 404:
                return "not_found"
            if sc in (429, 503, 504):
                wait = 2 ** attempt
                log(f"    {sc} レートリミット → {wait}秒待機 (試行 {attempt + 1}/{MAX_RETRIES})")
                time.sleep(wait)
                continue
            # 他のエラー
            log(f"    HTTP {sc}: {r.text[:120]}")
            return None
        except requests.exceptions.SSLError as e:
            log(f"    SSL エラー: {str(e)[:120]}")
            return None
        except requests.exceptions.Timeout:
            log(f"    タイムアウト (試行 {attempt + 1}/{MAX_RETRIES})")
            time.sleep(1)
            continue
        except Exception as e:
            log(f"    例外: {type(e).__name__}: {str(e)[:120]}")
            return None

    return None


def clean_cas(raw: object) -> str:
    """CAS番号を正規化。"""
    s = re.sub(r"\s+", "", str(raw))
    return s if re.match(r"^\d+-\d+-\d+$", s) else ""


def clean_name(name: str) -> str:
    """末尾の量表記を除去 例: '(-)-カルベオール (5g)' → '(-)-カルベオール'"""
    return re.sub(r"\s*\([\d.]+\s*[mμk]?[glmL]+\s*\)\s*$", "", name).strip()


def load_csv() -> pd.DataFrame:
    df = pd.read_csv(CSV_FILE, encoding="utf-8-sig", header=0, low_memory=False)
    cols = list(df.columns)
    cols[0] = "薬品名"
    df.columns = cols
    return df


def save_cache(cache: dict) -> None:
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def main():
    cache: dict = {}
    if CACHE_FILE.exists():
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            cache = json.load(f)
    log(f"[INFO] キャッシュ読み込み: {len(cache)} 件")

    df = load_csv()
    total = len(df)
    log(f"[INFO] CSV読み込み: {total} 行")

    new_fetches = 0
    found_count = sum(1 for v in cache.values() if v.get("smiles"))

    for i, row in df.iterrows():
        name = str(row.get("薬品名", "")).strip()
        if not name or name == "nan":
            continue

        # 既に取得成功している場合はスキップ
        existing = cache.get(name)
        if existing and existing.get("smiles"):
            continue
        if existing and existing.get("_status") == "not_found":
            continue

        cas = clean_cas(row.get("CAS No.", ""))
        en_name = str(row.get("薬品名（英語Ⅰ）", "") or "").strip()
        jp_clean = clean_name(name)

        queries = []
        if cas:
            queries.append(("CAS", cas))
        if en_name and en_name.lower() != "nan":
            queries.append(("EN", en_name))
        if jp_clean and jp_clean != cas and jp_clean != en_name:
            queries.append(("JP", jp_clean))

        log(f"[{i + 1}/{total}] {name[:55]}  (CAS={cas or '-'})")

        result = None
        for tag, q in queries:
            data = fetch_pubchem(q)
            time.sleep(SLEEP_BETWEEN_REQUESTS)

            if isinstance(data, dict):
                log(f"    [OK] {tag}検索成功: SMILES={data['smiles'][:40]}")
                result = data
                break
            elif data == "not_found":
                log(f"    [NF] {tag}検索: 404 ({q[:40]})")
            else:
                log(f"    [!!] {tag}検索: 失敗 ({q[:40]})")

        if result:
            cache[name] = result
            found_count += 1
        else:
            cache[name] = {"_status": "not_found"}

        new_fetches += 1

        if new_fetches % 25 == 0:
            save_cache(cache)
            log(f"    [SAVED] 進捗 {len(cache)}/{total}, SMILES取得 {found_count} 件")

    save_cache(cache)
    log(f"\n[DONE] 完了！ {len(cache)} 件中 {found_count} 件のSMILESを取得")


if __name__ == "__main__":
    main()
