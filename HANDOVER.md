# 植草研 薬品ストックリスト 引き継ぎマニュアル

このシステムを引き継ぐ後任者向けのドキュメントです。「何ができるか」「どこに何があるか」「どうやって維持するか」を一通り書いています。

---

## 1. システム概要

研究室の薬品在庫を、化合物名・CAS番号・官能基・保管場所などで検索・絞り込みできるWebアプリです。

- **公開URL**: <https://uekusa-stocklist-new.streamlit.app>
- **データ更新**: Google Sheetsを編集すると最大5分で自動反映
- **アクセス制限**: Streamlitサインイン（招待制）+ パスワードの二重ロック
- **当該PC不要**: クラウドホスティング、初代運用者のPCを介さず動作

---

## 2. 構成図

```
[利用者のブラウザ]
       ↓
[Streamlit Community Cloud]  ← ここでアプリが24時間稼働（無料プラン）
       ↓ コード読み込み
[GitHub: amkent-lgtm/stocklist-new]  ← Privateリポジトリ
       ↓ アプリ実行時にデータ取得
[Google Sheets 公開CSV]  ← 試薬データの真の出典
```

| 層 | 役割 | サービス |
|---|---|---|
| データ | 試薬リスト本体（追加・編集・削除はここ） | Google Sheets |
| ロジック | アプリのコード（変更したい時はここ） | GitHub Private repo |
| 実行環境 | アプリを24時間動かす場所 | Streamlit Community Cloud |
| 補助データ | PubChemから取得したSMILES/構造式キャッシュ | リポジトリ内の `compounds_enriched.json` |

---

## 3. アカウント一覧と引き継ぎ事項

| サービス | 現所有者 | アカウント | 備考 |
|---|---|---|---|
| GitHub | 〔初代運用者〕 | `amkent-lgtm` | リポジトリ `stocklist-new` のオーナー |
| Streamlit Cloud | 〔初代運用者〕 | `am.kent0617@gmail.com`（GitHub認証） | アプリのデプロイ管理 |
| Google アカウント | 〔初代運用者〕 | `am.kent0617@gmail.com` | Sheets のオーナー |
| 共有パスワード | 研究室メンバー共有 | — | Streamlit Cloud の Secrets に保管 |

⚠️ **重要**: いずれかのアカウントが削除されるとシステムが停止します。**継続運用するには引き継ぎが必須**です。後任者は本ドキュメント末尾の「7. 引き継ぎ手順」を実施してください。

---

## 4. 日常運用 — よくある作業マニュアル

### 4-1. 試薬を追加・編集・削除する

Google Sheetsを開いて直接編集してください。アプリは自動で最大5分以内に反映します。即時反映したい場合はサイドバーの **「🔄 データを再取得」** をクリック。

**Google Sheets URL（限定共有）**: 初代運用者から引き継ぎ時に伝達してください。

⚠️ Sheets の **列の追加・削除・順番変更・列名変更** は絶対にしないでください。アプリのコードが特定の列名（`CAS No.`、`A薬品庫右` など）に依存しています。**行の追加・編集・削除のみ** 安全です。

### 4-2. メンバーを追加する（アプリにアクセスさせる）

1. <https://share.streamlit.io> にサインイン（オーナーアカウントで）
2. 該当アプリの右下「**⋮**」→「**Settings**」をクリック
3. 左メニュー **「Sharing」** タブ
4. **「Invite viewers by email」** 欄に追加したいメンバーのメールアドレスを追記
5. **「Save changes」** をクリック
6. 新メンバーに以下を伝える：
   - URL: <https://uekusa-stocklist-new.streamlit.app>
   - 招待されたメールアドレスのGoogleアカウントでサインインすること
   - アプリのパスワード（合言葉）

### 4-3. メンバーから外す

同じ画面でメールアドレスを削除して **「Save changes」**。

### 4-4. アプリのパスワード（合言葉）を変更する

1. <https://share.streamlit.io> → 該当アプリの **「Settings」**
2. **「Secrets」** タブ
3. `password = "新しいパスワード"` に書き換え
4. **「Save」** をクリック（数秒で反映）
5. 全メンバーに新パスワードを通知

### 4-5. アプリの機能を変更したい（コード変更）

1. このリポジトリ（`amkent-lgtm/stocklist-new`）をローカルにクローン or ブラウザでGitHubから直接編集
2. `app.py` を修正
3. `git push` すれば Streamlit Cloud が自動再デプロイ（1〜2分）

ローカル開発する場合の事前準備：
```powershell
pip install -r requirements.txt
streamlit run app.py
```

### 4-6. 構造式データを最新化する（PubChem再取得）

Google Sheetsに新しい試薬を大量追加した場合、構造式キャッシュが古くなります。最新化手順：

```powershell
python enrich_data.py
```

ローカルで実行 → `compounds_enriched.json` が更新される → `git add compounds_enriched.json && git commit && git push` で反映。

5〜10分かかります。詳細は `enrich_data.py` 冒頭のコメント参照。

---

## 5. ファイル構成

```
植草研_stocklist/
├── app.py                            # メインアプリのコード
├── enrich_data.py                    # PubChem から SMILES/InChI を取得するスクリプト
├── requirements.txt                  # Python依存パッケージ
├── packages.txt                      # Streamlit Cloud に追加するシステムライブラリ
├── compounds_enriched.json           # PubChemから取得した構造情報キャッシュ
├── stocklist_20251023 のコピー - ....csv  # ローカルCSVフォールバック
├── .streamlit/
│   └── config.toml                   # Streamlit設定（テーマ・ツールバー非表示など）
├── .gitignore
└── HANDOVER.md                       # この文書
```

---

## 6. トラブルシューティング

| 症状 | 原因の可能性 | 対処 |
|---|---|---|
| アプリが「Sleep」している | 何日もアクセスがないと Streamlit Cloud がスリープ | URLにアクセス → 10秒程度で復帰 |
| データが古い | キャッシュ TTL（5分） | サイドバーの「🔄 データを再取得」をクリック |
| 構造式が表示されない | SMILES未取得（22件） or PubChem変更 | `enrich_data.py` を再実行 |
| 「Access denied」と出る | 招待リストに無い、または別のGoogleアカウントでサインイン中 | 4-2の手順で招待 / 招待されたアカウントでサインインし直す |
| メンバーがGoogleサインインから先に進まない | 招待リスト未追加 | 4-2の手順で追加 |
| Streamlit Cloud から「ビルド失敗」の通知 | requirements.txt のパッケージバージョン不一致など | Streamlit Cloud の **Manage app → Logs** で詳細確認 |
| Google Sheets が読めなくなった | 公開設定がOFFになった | Sheets を開いて「ファイル→共有→ウェブに公開」を再設定 |
| パスワードを入れても通らない | Secretsの書式エラー | `password = "値"` 形式（クォート必須）になっているか確認 |

### ログ確認方法
1. <https://share.streamlit.io>
2. 該当アプリの「**Manage app**」をクリック
3. 右側にログが表示。エラーメッセージを確認

---

## 7. 引き継ぎ手順（重要）

**初代運用者が異動・卒業する前に必ず実施してください**。

### 推奨：研究室共有アカウント方式

1. **研究室共有のGoogleアカウントを作成**（例: `uesa.lab.uoa@gmail.com`）
2. そのアカウントで **GitHubアカウント** を新規作成
3. リポジトリを移転：
   - GitHubで `amkent-lgtm/stocklist-new` の **Settings** → **Transfer ownership**
   - 移転先に研究室共有GitHubアカウントを指定
4. **Streamlit Cloud に研究室共有Googleアカウントでサインイン** → 新アプリとして再デプロイ
   - 旧URLが変わるので、新URLをメンバーに再周知
5. **Google Sheets のオーナー権限を共有アカウントに移譲**：
   - Sheets を開いて **「共有」** → 共有アカウントを「編集者」追加
   - 右上3点メニュー → **「オーナーシップを移行」**
6. 旧アプリ（初代運用者管理）を停止：
   - Streamlit Cloud の旧アプリ「Manage app」→「Delete」
7. 共有アカウントの **ID/パスワードを教員 or 後任者に文書で引き継ぐ**

### 代替：個人後任者への譲渡

1. GitHubで `amkent-lgtm/stocklist-new` を後任者の個人アカウントへ Transfer
2. 後任者が Streamlit Cloud に自分のアカウントでログイン → リポジトリを再デプロイ
3. Google Sheets のオーナーシップを後任者に移譲
4. メンバー招待リスト・パスワードは新環境で再設定

---

## 8. 緊急時の連絡先

| 用件 | 連絡先 |
|---|---|
| 初代運用者（システム構築者） | am.kent0617@gmail.com |
| Streamlit Cloud サポート | <https://discuss.streamlit.io> |
| Google Sheets ヘルプ | <https://support.google.com/docs> |

---

## 9. 技術スタック

このアプリで使用している主な技術：

- **Python 3.14**（Streamlit Cloud デフォルト）
- **Streamlit 1.57** — Web アプリフレームワーク
- **pandas 3.0** — CSV/Sheets データ処理
- **RDKit 2026.3** — 化学構造式の描画・官能基検出
- **requests** — PubChem APIへの問い合わせ（`enrich_data.py`内）
- **truststore** — Windows SSL証明書ストア利用（`enrich_data.py`内）

---

## 10. ライセンス・引用について

このアプリのデータ源：
- **PubChem**（米国国立医学図書館 NCBI） — SMILES / InChI / 分子式 / 分子量
  引用が必要な場合: <https://pubchem.ncbi.nlm.nih.gov/about>
- **薬品ストックリスト** — 植草研究室の在庫情報

---

*最終更新: 2026-05-27*  
*作成: 初代運用者*
