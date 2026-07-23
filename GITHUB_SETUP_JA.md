# GitHub公開・pkgdownサイト設定手順

このリポジトリ一式は、CRANへ提出した最終ソース
`TTE_1.1.1.tar.gz`を基礎として作成されています。

## 1. GitHubで空のリポジトリを作る

1. <https://github.com/new> を開く
2. Owner：`nomahi`
3. Repository name：`TTE`
4. Publicを選択
5. GitHub側ではREADME、LICENSE、`.gitignore`を追加しない
6. Create repository

## 2. 最初のpush

展開した `TTE` フォルダで実行します。

```bash
git init
git add .
git commit -m "Initial public release of TTE 1.1.1"
git branch -M main
git remote add origin https://github.com/nomahi/TTE.git
git push -u origin main
```

## 3. GitHub Pagesを有効化

1. `nomahi/TTE` の **Settings**
2. 左側の **Pages**
3. **Build and deployment**
4. Sourceとして **GitHub Actions** を選択

`.github/workflows/pkgdown.yaml` が、pkgdownサイトを作成して
GitHub Pagesへ配信します。

公開予定URL：

<https://nomahi.github.io/TTE/>

## 4. Actionsの確認

**Actions** タブで次のworkflowを確認します。

- `R-CMD-check.yaml`
- `pkgdown`

初回はworkflowの実行許可が必要な場合があります。

## 5. About欄の推奨設定

- Description:
  `Design and analysis tools for target trial emulation in R`
- Website:
  `https://nomahi.github.io/TTE/`
- Topics:
  `r-package`, `causal-inference`, `target-trial-emulation`,
  `inverse-probability-weighting`, `survival-analysis`,
  `competing-risks`, `biostatistics`, `epidemiology`

## 6. CRAN掲載後のrelease

```bash
git tag -a v1.1.1 -m "TTE 1.1.1"
git push origin v1.1.1
```

GitHubの **Releases → Draft a new release** で `v1.1.1` を選びます。
CRANへ提出したtar.gzをRelease assetとして添付しても構いません。

## 7. 通常の更新

```bash
git add .
git commit -m "Describe the change"
git push
```

push後、R CMD checkとpkgdownサイト構築が自動実行されます。

## 8. 将来のCRAN用tar.gz

GitHub専用ファイルは `.Rbuildignore` によって除外されます。

```r
devtools::test()
devtools::check()
devtools::build()
```

CRANへ送るのは `devtools::build()` が生成した `.tar.gz` です。
