"""
化合物の pKa を取得する。

第1優先: PubChem PUG-View の "Dissociation Constants"（文献由来の実験値）
第2優先: 実験値がない化合物は RDKit SMARTS パターンマッチで官能基ベース予測

実行方法:
    python enrich_pka.py

- 既に "pka" フィールドが入っている化合物はスキップ
- 結果は compounds_enriched.json に追記保存
- 50件ごとに途中保存
"""

import json
import re
import sys
import time
from pathlib import Path

import requests

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

try:
    import truststore
    truststore.inject_into_ssl()
    print("[OK] truststore injected", flush=True)
except ImportError:
    pass

from rdkit import Chem

BASE_DIR = Path(__file__).parent
CACHE_FILE = BASE_DIR / "compounds_enriched.json"

PUGVIEW_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
)
SLEEP = 0.34
MAX_RETRIES = 3


# ── 予測ルール（官能基 → 代表的 pKa） ────────────────────────────────────────
# (group_name, SMARTS, pKa, type)
# type は表示用: "acid" (HA→H+ + A-) / "base" (B + H+ → BH+ の共役酸 pKa)
PKA_RULES: list[tuple[str, str, float, str]] = [
    # 強酸〜中強度の酸（具体的なものを先に）
    ("スルホン酸",         "[SX4](=O)(=O)[OX2H1]",                   -2.5, "acid"),
    ("ホスホン酸",         "[PX4](=O)([OX2H1])[OX2H1]",               2.0, "acid"),
    ("ヒドロキサム酸",     "[CX3](=O)[NX3H1][OX2H1]",                 9.0, "acid"),
    ("テトラゾール (NH)",  "c1nnn[nH]1",                              4.9, "acid"),
    ("カルボン酸",         "[CX3](=O)[OX2H1]",                        4.5, "acid"),
    ("フェノール",         "[c][OX2H1]",                             10.0, "acid"),
    ("チオール",           "[#6][SX2H1]",                            10.5, "acid"),
    # 塩基（ pKa は共役酸の値）
    ("グアニジン",         "[NX3H2][CX3](=[NX2H1,NX3H2])[NX3H1,NX3H2]", 13.0, "base"),
    ("アミジン",           "[NX3H2,NX3H1][CX3]=[NX2H1,NX3H2]",       11.5, "base"),
    ("イミダゾール",       "[nH0]1cnc[cH]1",                          7.0, "base"),
    ("ピリジン",           "[nX2H0;R1]",                              5.2, "base"),
    ("芳香族アミン",       "[cX3][NX3H2]",                            4.5, "base"),
    ("第一級アミン (脂肪)", "[NX3H2;!$(N=*);!$(Nc);!$(N[C,S,P]=O)]",  10.6, "base"),
    ("第二級アミン (脂肪)", "[NX3H1;!$(N=*);!$(Nc);!$(N[C,S,P]=O)]",  10.8, "base"),
    ("第三級アミン (脂肪)", "[NX3H0;!$(N=*);!$(Nc);!$(N[C,S,P]=O);!R]", 10.0, "base"),
]
_COMPILED_RULES = [
    (name, Chem.MolFromSmarts(smarts), pka, ptype)
    for name, smarts, pka, ptype in PKA_RULES
]
_COMPILED_RULES = [r for r in _COMPILED_RULES if r[1] is not None]


# ── PubChem PUG-View 実験値 ─────────────────────────────────────────────────
def _walk_sections(node, target_heading: str):
    """JSON ツリーを再帰走査して、指定 TOCHeading のセクションを yield。"""
    if isinstance(node, dict):
        if node.get("TOCHeading") == target_heading:
            yield node
        for v in node.values():
            yield from _walk_sections(v, target_heading)
    elif isinstance(node, list):
        for v in node:
            yield from _walk_sections(v, target_heading)


_PKA_LABEL_RE = re.compile(r"p[Kk]a\d?", re.IGNORECASE)
_NUM_RE = re.compile(r"-?\d+\.\d+|-?\d+")


def _parse_pka_strings(strings: list[str]) -> list[float]:
    """Dissociation Constants セクション内のテキストから pKa 数値を抽出。

    PubChem の "Dissociation Constants" セクション内では：
    - "3.47"  (数値だけ)
    - "pKa = 3.49 at 25 °C"
    - "pKa1 = 3.13, pKa2 = 4.76"
    などの形式で記載されている。

    手順: pKaラベル(pKa, pKa1 など)を除去 → 全ての小数/整数を抽出
    → pKa 妥当範囲 [-3, 18] でフィルタ → 0.1刻みで重複除去
    """
    values: list[float] = []
    seen = set()
    for s in strings:
        s_clean = _PKA_LABEL_RE.sub(" ", s)
        for m in _NUM_RE.findall(s_clean):
            try:
                v = float(m)
            except ValueError:
                continue
            # 温度・参照番号・室温(25)などを除外
            if v < -3 or v > 18:
                continue
            key = round(v, 1)  # 0.1 刻みで重複除去（3.47 と 3.5 は同じ pKa）
            if key in seen:
                continue
            seen.add(key)
            values.append(round(v, 2))
    return values


def fetch_experimental_pka(cid: int) -> list[dict]:
    """PubChem PUG-View から実験 pKa を取得。"""
    if not cid:
        return []
    url = PUGVIEW_URL.format(cid=cid)
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, timeout=20)
            sc = r.status_code
            if sc == 200:
                data = r.json()
                break
            if sc == 404:
                return []
            if sc in (429, 503, 504):
                time.sleep(2 ** attempt)
                continue
            return []
        except requests.exceptions.Timeout:
            time.sleep(1)
            continue
        except Exception:
            return []
    else:
        return []

    # 全 Information の文字列を集めてから一度にパース（化合物単位で重複除去するため）
    all_strings: list[str] = []
    for section in _walk_sections(data, "Dissociation Constants"):
        for info in section.get("Information", []):
            value = info.get("Value", {})
            for sm in value.get("StringWithMarkup", []):
                if "String" in sm:
                    all_strings.append(sm["String"])
    return [{"value": v, "type": "exp"} for v in _parse_pka_strings(all_strings)]


# ── 官能基ベース予測 ────────────────────────────────────────────────────────
def predict_pka(smiles: str) -> list[dict]:
    """SMARTS マッチで pKa を予測。マッチした官能基1つにつき1値を返す。"""
    if not smiles:
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    # 同じ原子が複数ルールにマッチする場合の重複防止
    used_atoms: set[int] = set()
    predictions = []
    for name, pat, pka, ptype in _COMPILED_RULES:
        matches = mol.GetSubstructMatches(pat)
        for atoms in matches:
            atom_set = frozenset(atoms)
            if atom_set & used_atoms:
                continue  # 既にこの原子は別ルールでマッチ済み
            used_atoms |= atom_set
            predictions.append({
                "value": pka,
                "type": "pred",
                "group": name,
                "kind": ptype,
            })
    return predictions


# ── メイン ────────────────────────────────────────────────────────────────
def main() -> None:
    with open(CACHE_FILE, "r", encoding="utf-8") as f:
        cache = json.load(f)

    total = len(cache)
    print(f"[INFO] {total} 件をチェック", flush=True)

    new_processed = 0
    exp_count = 0
    pred_count = 0

    for i, (name, data) in enumerate(cache.items()):
        if not isinstance(data, dict):
            continue
        if "pka" in data:
            # 既に処理済みはスキップ（再実行時の負荷軽減）
            continue

        smiles = data.get("smiles", "")
        cid = int(data.get("pubchem_cid", 0) or 0)

        if not smiles:
            data["pka"] = []
            continue

        print(f"[{i + 1}/{total}] {name[:55]}", flush=True)

        # 第1優先: 実験値
        exp = fetch_experimental_pka(cid)
        time.sleep(SLEEP)

        if exp:
            data["pka"] = exp
            exp_count += 1
            print(f"    [EXP] {len(exp)}個: {[v['value'] for v in exp]}", flush=True)
        else:
            pred = predict_pka(smiles)
            data["pka"] = pred
            if pred:
                pred_count += 1
                print(f"    [PRED] {len(pred)}個: "
                      f"{[(v['group'], v['value']) for v in pred]}", flush=True)
            else:
                print(f"    [NONE] イオン化基なし", flush=True)

        new_processed += 1
        if new_processed % 50 == 0:
            _save(cache)
            print(f"    [SAVED] 新規処理 {new_processed} 件 "
                  f"(実験 {exp_count} / 予測 {pred_count})", flush=True)

    _save(cache)
    print(f"\n[DONE] 新規処理 {new_processed} 件 "
          f"(実験 {exp_count} / 予測 {pred_count})", flush=True)


def _save(cache: dict) -> None:
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
