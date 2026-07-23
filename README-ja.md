# TTE

**Design and Analysis Tools for Target Trial Emulation**

`TTE`は、縦断的観察データを用いた標的試験エミュレーションの
設計・解析を支援するRパッケージです。

このリポジトリは、CRAN初回リリースである **TTE 1.1.1** に対応します。

- パッケージサイト：<https://nomahi.github.io/TTE/>
- CRAN：<https://CRAN.R-project.org/package=TTE>
- PDFマニュアル：[`manual/TTE_1.1.1.pdf`](manual/TTE_1.1.1.pdf)

## インストール

```r
install.packages("TTE")
```

GitHub上の開発版：

```r
install.packages("remotes")
remotes::install_github("nomahi/TTE")
```

## 主な関数

- `check_tte()`：person-periodデータの構造確認
- `seqdesign_tte()`：sequentially nested trialsの構築
- `est_wt()`：treatment、censoring、adherence weightの推定
- `combine_wt()`：複数のweight componentの結合
- `balance_wt()`：重み付け前後の共変量バランス
- `diagnose_wt()`：weight、ESS、weighted risk setの診断
- `discsurvreg()`：weighted pooled discrete-time outcome model
- `std_tte()`：標準化リスク、CIF、RD、RR、NNT/NNH
- `curve_tte()`：weighted Kaplan-Meier／Aalen-Johansen
- `boot_tte()`：元の個人単位のcluster bootstrap

## 公開チュートリアル

### SGLT2阻害薬 versus DPP-4阻害薬：全死亡

- [日本語版](tutorials/ja/SGLT2_TTE_tutorial_ja.R)
- [英語版](tutorials/en/SGLT2_TTE_tutorial_en.R)

### ARB versus CCB：心不全入院、死亡を競合イベントとして考慮

- [日本語版](tutorials/ja/ARB_TTE_tutorial_ja.R)
- [英語版](tutorials/en/ARB_TTE_tutorial_en.R)

両チュートリアルとも、収載済みの参照weightを主解析に使わず、
baseline IPTW、LTFU-IPCW、PP解析用joint censoring weightを
合成データから推定します。

## 合成データ

`SGLT2`、`SGLT2_baseline`、`ARB`、`ARB_baseline`は、すべて
教育・ソフトウェア検証・再現可能な事例のための完全な合成データです。
実在する患者記録や元研究の実データは含まれていません。

## ライセンス

GPL-3.
