"""
植草研 薬品ストックリスト検索システム

実行方法:
    streamlit run app.py
"""

import base64
from io import BytesIO
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="植草研 薬品リスト",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "About": "植草研 薬品ストックリスト検索システム\n\nPubChemから取得した構造式・官能基情報で薬品を検索できます。",
        "Get help": None,
        "Report a bug": None,
    },
)

# Streamlit標準UI要素を非表示にするCSS
st.markdown(
    """
    <style>
    /* 右上のヘッダー（Deploy・三点リーダーを含む領域）を非表示 */
    header[data-testid="stHeader"] { display: none; }
    /* "Made with Streamlit" フッターを非表示 */
    footer { display: none; }
    /* メインコンテナの上余白を縮める */
    .block-container { padding-top: 1.5rem; }
    /* サイドバーの上余白も縮める */
    section[data-testid="stSidebar"] .block-container { padding-top: 1.5rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

BASE_DIR = Path(__file__).parent
CSV_FILE = BASE_DIR / "stocklist_20251023 のコピー - stock_list_20251011175958.csv"
CACHE_FILE = BASE_DIR / "compounds_enriched.json"

# ── 官能基 SMARTS 定義 ─────────────────────────────────────────────────────────
FUNCTIONAL_GROUPS: dict[str, str] = {
    "アルコール":       "[OX2H][CX4]",
    "フェノール":       "[OX2H]c",
    "エーテル":         "[OD2]([#6])[#6]",
    "エポキシ":         "[C]1[O][C]1",
    "アルデヒド":       "[CX3H1](=O)[#6]",
    "ケトン":           "[#6][CX3](=O)[#6]",
    "カルボン酸":       "[CX3](=O)[OX2H1]",
    "エステル":         "[#6][CX3](=O)[OX2][#6]",
    "アミド":           "[NX3][CX3](=O)[#6]",
    "第一級アミン":     "[NX3H2;!$(NC=O)]",
    "第二級アミン":     "[NX3H1;!$(NC=O)]([#6])[#6]",
    "第三級アミン":     "[NX3H0;!$(NC=O)]([#6])([#6])[#6]",
    "ニトリル":         "[NX1]#[CX2]",
    "ニトロ基":         "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "アルケン":         "[CX3]=[CX3]",
    "アルキン":         "[CX2]#[CX2]",
    "芳香環":           "c1ccccc1",
    "塩素 (Cl)":        "[Cl]",
    "臭素 (Br)":        "[Br]",
    "フッ素 (F)":       "[F]",
    "ヨウ素 (I)":       "[I]",
    "チオール":         "[SX2H]",
    "チオエーテル":     "[#6][SX2][#6]",
    "スルホニル":       "[#6][SX4](=O)(=O)[#6]",
    "イミン":           "[CX3]=[NX2H1]",
    "ヒドラジン":       "[NX3H1][NX3H2]",
    "ボロン酸":         "[BX3]([OX2H])[OX2H]",
}

STORAGE_COLS = [
    "A薬品庫右", "A薬品庫左", "B薬品庫",
    "冷蔵①", "冷蔵②", "冷蔵③", "冷蔵④", "冷蔵⑤", "冷蔵⑥",
    "冷凍", "D薬品庫", "無水箱", "個人所有", "その他", "3階試薬庫",
]


# ── データ読み込み ──────────────────────────────────────────────────────────────
@st.cache_data
def load_data() -> pd.DataFrame:
    import json

    df = pd.read_csv(CSV_FILE, encoding="utf-8-sig", header=0, low_memory=False)
    cols = list(df.columns)
    cols[0] = "薬品名"
    df.columns = cols

    enrichment: dict = {}
    if CACHE_FILE.exists():
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            enrichment = json.load(f)

    def _get(name: str, key: str, default=""):
        return enrichment.get(str(name).strip(), {}).get(key, default)

    df["_smiles_fetched"] = df["薬品名"].map(lambda n: _get(n, "smiles"))
    df["_inchikey"]       = df["薬品名"].map(lambda n: _get(n, "inchikey"))
    df["_mol_formula"]    = df["薬品名"].map(lambda n: _get(n, "molecular_formula"))
    df["_mol_weight"]     = pd.to_numeric(
        df["薬品名"].map(lambda n: _get(n, "molecular_weight", None)), errors="coerce"
    )
    df["_pubchem_cid"]    = df["薬品名"].map(lambda n: _get(n, "pubchem_cid", 0))

    # CSV既存SMILESを優先、なければ取得分を使用
    csv_smiles = df.get("SMILES", pd.Series([""] * len(df), index=df.index))
    valid = csv_smiles.notna() & (csv_smiles.astype(str).str.strip() != "") & (csv_smiles.astype(str) != "nan")
    df["_smiles"] = csv_smiles.where(valid, df["_smiles_fetched"]).fillna("")

    # 保管場所
    def _storage(row):
        locs = [c for c in STORAGE_COLS if c in df.columns and str(row.get(c, "")).upper() == "TRUE"]
        return "、".join(locs) if locs else "不明"

    df["_storage"] = df.apply(_storage, axis=1)

    return df


# ── 官能基検出 (RDKit) ─────────────────────────────────────────────────────────
@st.cache_data
def compute_functional_groups(smiles_series: pd.Series) -> list[list[str]]:
    if not RDKIT_AVAILABLE:
        return [[] for _ in range(len(smiles_series))]

    patterns = {}
    for name, smarts in FUNCTIONAL_GROUPS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat:
            patterns[name] = pat

    results = []
    for smi in smiles_series:
        if not smi or pd.isna(smi) or str(smi).strip() == "":
            results.append([])
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            results.append([])
        else:
            results.append([n for n, p in patterns.items() if mol.HasSubstructMatch(p)])
    return results


# ── 構造式画像 ──────────────────────────────────────────────────────────────────
@st.cache_data
def smiles_to_png_b64(smiles: str, width: int = 280, height: int = 200) -> str:
    """SMILESをbase64 PNG文字列に変換。失敗時は空文字。"""
    if not RDKIT_AVAILABLE or not smiles or str(smiles).strip() == "":
        return ""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    img = Draw.MolToImage(mol, size=(width, height))
    buf = BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()


# ── カード HTML ─────────────────────────────────────────────────────────────────
def render_card(row: pd.Series, img_b64: str) -> str:
    name     = str(row.get("薬品名", ""))
    en_name  = str(row.get("薬品名（英語Ⅰ）", "") or "")
    cas      = str(row.get("CAS No.", "") or "").strip()
    formula  = str(row.get("_mol_formula", "") or "")
    mw       = row.get("_mol_weight", "")
    storage  = str(row.get("_storage", ""))
    unopened = int(float(str(row.get("未開封", 0) or 0)))
    opened   = int(float(str(row.get("開封", 0) or 0)))
    fg_list  = row.get("_functional_groups", [])
    cid      = row.get("_pubchem_cid", 0)

    mw_str = f"{mw:.1f}" if isinstance(mw, float) and mw > 0 else ""

    img_html = (
        f'<img src="data:image/png;base64,{img_b64}" '
        'style="width:100%;border-radius:6px;background:#fff;" />'
        if img_b64
        else '<div style="height:160px;background:#f5f5f5;display:flex;'
             'align-items:center;justify-content:center;color:#aaa;'
             'font-size:12px;border-radius:6px;">構造式なし</div>'
    )
    fg_html = "".join(
        f'<span style="background:#e3f2fd;color:#1565c0;padding:1px 7px;'
        f'border-radius:10px;font-size:11px;margin:2px;display:inline-block;">'
        f'{fg}</span>'
        for fg in fg_list
    )
    stock_color = "#2e7d32" if (unopened + opened) > 0 else "#c62828"
    stock_label = f"未開封 {unopened} / 開封 {opened}"

    pubchem_link = (
        f'<a href="https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" '
        'target="_blank" style="font-size:11px;color:#1976d2;">PubChem ↗</a>'
        if cid else ""
    )

    return f"""
<div style="border:1px solid #e0e0e0;border-radius:10px;padding:12px;
            background:white;height:100%;box-sizing:border-box;">
  {img_html}
  <div style="margin-top:8px;">
    <div style="font-size:12px;font-weight:600;line-height:1.4;
                margin-bottom:3px;word-break:break-word;">{name}</div>
    <div style="font-size:11px;color:#555;margin-bottom:4px;
                font-style:italic;">{en_name}</div>
    <div style="font-size:11px;color:#777;">
      CAS: {cas} &nbsp;|&nbsp; {formula}
      {"&nbsp;|&nbsp; MW: " + mw_str if mw_str else ""}
    </div>
    <div style="font-size:11px;margin-top:4px;">📦 {storage}</div>
    <div style="font-size:11px;color:{stock_color};margin-top:2px;">
      {stock_label}
    </div>
    <div style="margin-top:6px;">{fg_html}</div>
    <div style="margin-top:4px;">{pubchem_link}</div>
  </div>
</div>"""


# ── メインアプリ ────────────────────────────────────────────────────────────────
def main():
    st.title("⚗️ 植草研 薬品ストックリスト")

    with st.spinner("データ読み込み中…"):
        df = load_data()

    if not RDKIT_AVAILABLE:
        st.warning(
            "⚠️ RDKitが未インストールのため構造式表示・官能基検索が使えません。  \n"
            "`pip install rdkit` を実行後アプリを再起動してください。"
        )

    # 官能基を計算してDataFrameに追加
    if RDKIT_AVAILABLE:
        with st.spinner("官能基を解析中…（初回のみ）"):
            fg_lists = compute_functional_groups(df["_smiles"])
        df["_functional_groups"] = fg_lists
    else:
        df["_functional_groups"] = [[] for _ in range(len(df))]

    # ── サイドバー ──────────────────────────────────────────────────────────────
    st.sidebar.header("🔍 検索・フィルター")

    search = st.sidebar.text_input(
        "化合物名で検索", placeholder="例: エタノール、benzene…"
    )

    selected_fgs = st.sidebar.multiselect(
        "官能基で絞り込み（AND検索）",
        options=list(FUNCTIONAL_GROUPS.keys()),
        disabled=not RDKIT_AVAILABLE,
        help="選択した官能基をすべて含む化合物を表示します",
    )

    storage_options = ["すべて"] + [c for c in STORAGE_COLS if c in df.columns]
    selected_storage = st.sidebar.selectbox("保管場所で絞り込み", storage_options)

    only_in_stock = st.sidebar.checkbox("在庫あり（未開封＋開封 > 0）のみ", value=False)

    st.sidebar.divider()
    st.sidebar.subheader("並び替え")
    sort_key = st.sidebar.selectbox(
        "基準",
        ["薬品名（昇順）", "薬品名（降順）", "分子量（昇順）", "分子量（降順）", "分子式"],
    )

    st.sidebar.divider()
    view_mode = st.sidebar.radio("表示形式", ["カード", "テーブル"])
    per_page = st.sidebar.select_slider(
        "1ページの表示件数", options=[12, 24, 48, 100], value=24
    )

    # ── フィルター適用 ──────────────────────────────────────────────────────────
    mask = pd.Series(True, index=df.index)

    if search:
        q = search.lower()
        name_hit    = df["薬品名"].str.lower().str.contains(q, na=False)
        en_name_hit = df.get("薬品名（英語Ⅰ）", pd.Series("", index=df.index)).str.lower().str.contains(q, na=False)
        cas_hit     = df.get("CAS No.", pd.Series("", index=df.index)).astype(str).str.contains(q, na=False)
        mask &= name_hit | en_name_hit | cas_hit

    if selected_fgs and RDKIT_AVAILABLE:
        mask &= df["_functional_groups"].apply(
            lambda fgs: all(fg in fgs for fg in selected_fgs)
        )

    if selected_storage != "すべて" and selected_storage in df.columns:
        mask &= df[selected_storage].astype(str).str.upper() == "TRUE"

    if only_in_stock:
        unopened = pd.to_numeric(df.get("未開封", 0), errors="coerce").fillna(0)
        opened   = pd.to_numeric(df.get("開封", 0), errors="coerce").fillna(0)
        mask &= (unopened + opened) > 0

    filtered = df[mask].copy()

    # ── 並び替え ────────────────────────────────────────────────────────────────
    sort_map = {
        "薬品名（昇順）":  ("薬品名", True),
        "薬品名（降順）":  ("薬品名", False),
        "分子量（昇順）":  ("_mol_weight", True),
        "分子量（降順）":  ("_mol_weight", False),
        "分子式":          ("_mol_formula", True),
    }
    sort_col, sort_asc = sort_map[sort_key]
    filtered = filtered.sort_values(sort_col, ascending=sort_asc, na_position="last")

    # ── 統計 ────────────────────────────────────────────────────────────────────
    total_rows = len(df)
    smiles_count = (df["_smiles"].str.strip() != "").sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("総エントリ数", total_rows)
    c2.metric("フィルター後", len(filtered))
    c3.metric("SMILES取得済み", int(smiles_count))
    c4.metric(
        "SMILES未取得",
        int(total_rows - smiles_count),
        delta_color="inverse",
        help="enrich_data.py を実行すると取得できます",
    )

    # PubChem未実行の案内
    if smiles_count == 0:
        st.info(
            "💡 構造式データがありません。以下のコマンドでPubChemから一括取得できます：\n\n"
            "```\npython enrich_data.py\n```\n\n"
            "完了後アプリを再起動してください（約5〜10分）。"
        )

    # ── ページネーション ─────────────────────────────────────────────────────────
    total_pages = max(1, (len(filtered) + per_page - 1) // per_page)
    page_col, _ = st.columns([1, 5])
    with page_col:
        page = st.number_input("ページ", min_value=1, max_value=total_pages, value=1, step=1)
    page -= 1
    start = page * per_page
    end   = start + per_page
    page_df = filtered.iloc[start:end]

    st.caption(
        f"{start + 1}–{min(end, len(filtered))} 件目 / 全 {len(filtered)} 件"
        f"（{total_pages} ページ）"
    )

    # ── 表示 ────────────────────────────────────────────────────────────────────
    if view_mode == "カード":
        n_cols = 4
        rows = [page_df.iloc[i : i + n_cols] for i in range(0, len(page_df), n_cols)]
        for row_group in rows:
            cols = st.columns(n_cols)
            for j, (_, row) in enumerate(row_group.iterrows()):
                with cols[j]:
                    smi = str(row.get("_smiles", "")).strip()
                    img = smiles_to_png_b64(smi)
                    st.markdown(render_card(row, img), unsafe_allow_html=True)
                    st.write("")  # spacing

    else:
        display_cols = {
            "薬品名": "薬品名",
            "薬品名（英語Ⅰ）": "英語名",
            "CAS No.": "CAS No.",
            "_mol_formula": "分子式",
            "_mol_weight": "分子量",
            "_storage": "保管場所",
            "未開封": "未開封",
            "開封": "開封",
            "メーカー名": "メーカー",
        }
        show_cols = [c for c in display_cols if c in filtered.columns]
        renamed = {c: display_cols[c] for c in show_cols}
        st.dataframe(
            page_df[show_cols].rename(columns=renamed),
            use_container_width=True,
            hide_index=True,
        )

    # ── CSVエクスポート ──────────────────────────────────────────────────────────
    st.divider()
    if st.button("📥 フィルター結果をCSVダウンロード"):
        export_cols = [
            "薬品名", "薬品名（英語Ⅰ）", "CAS No.",
            "_smiles", "_mol_formula", "_mol_weight",
            "_storage", "未開封", "開封", "メーカー名",
        ]
        export_cols = [c for c in export_cols if c in filtered.columns]
        col_rename  = {
            "_smiles": "SMILES",
            "_mol_formula": "分子式",
            "_mol_weight": "分子量",
            "_storage": "保管場所",
        }
        csv_bytes = (
            filtered[export_cols]
            .rename(columns=col_rename)
            .to_csv(index=False, encoding="utf-8-sig")
            .encode("utf-8-sig")
        )
        st.download_button(
            label="ダウンロード",
            data=csv_bytes,
            file_name="stocklist_filtered.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    main()
