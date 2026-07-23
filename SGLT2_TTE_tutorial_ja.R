# =============================================================================
# TTE package：日本語ハンズオン・チュートリアル
#
# 1. Target Trial Emulation for SGLT2i versus DPP-4i:
#    All-cause death
# =============================================================================
#
# モデルとした研究：
# Noma, H., Goto, A., Sugimoto, T., Sunada, H., Oda, F., Maeda, M.,
# and Fukuda, H. (2026).
# Real-world effectiveness of SGLT2 inhibitors in adults aged 75 years
# or older: a target trial emulation.
# Age and Ageing, in press.
#
# どんな研究？
# - active-comparator new-user design
# - 75歳以上の2型糖尿病患者を対象
# - SGLT2阻害薬開始とDPP-4阻害薬開始を比較するRCTをemulation
# - primary outcomeは全死亡
# - 1か月のinduction periodを設定
# - intention-to-treat（ITT）effectとper-protocol（PP）effectを推定
#
# コンピュータ乱数により、この研究を模した小規模な合成データが
# TTEパッケージに収載されています。
#
#   DPP-4i群：497人
#   SGLT2i群：203人
#   合計      ：700人
#
# 【重要】
# - 収載データはすべて教育目的の合成データです。
# - 実在する患者記録や、元研究の実データは含まれていません。
# - このスクリプトから得られる数値は、実際の薬剤効果を示すものでは
#   ありません。
#
# この事例で学ぶ内容：
# - baseline treatment weight（IPTW）
# - loss-to-follow-upに対する累積IPCW
# - ITT解析用の最終weight
# - weighted pooled discrete-time outcome model
# - standardizationによる絶対リスク、risk difference、risk ratio
# - weighted Kaplan-Meier curve
# - protocol deviationを考慮したPP解析
#
# TTE version 1.1.1以降を想定しています。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 準備
# -----------------------------------------------------------------------------

# 未インストールの場合は、CRAN収載後に一度だけ実行します。
# install.packages("TTE")

library(TTE)

options(width = 110)

packageVersion("TTE")


# -----------------------------------------------------------------------------
# 1-1. データの読み込みと構造の確認
# -----------------------------------------------------------------------------

data(SGLT2_baseline)
data(SGLT2)

# Baselineデータ：
#   1人につき1行の治療開始時データです。
#   主として、baseline treatment model、IPTW、
#   ベースライン共変量の分布評価に使います。
dim(SGLT2_baseline)
head(SGLT2_baseline, 20)

# Longデータ：
#   1人が複数行のperson-monthデータとして記録されています。
#   pooled discrete-time modelやweighted survival curveに使います。
dim(SGLT2)
head(SGLT2, 6)

# 主要変数
#   id                    ：元の個人ID
#   trial                 ：emulated trialの識別子
#   time                  ：trial baselineからの経過月（0始まり）
#   A                     ：治療指標（0 = DPP-4i、1 = SGLT2i）
#   treatment             ：治療名を表すfactor
#   Y_death               ：その区間に全死亡が発生したか（0/1）
#   stay_ltfu             ：次区間まで追跡下に残ったか（0/1）
#   stay_adherent         ：次区間までassigned strategyに従ったか（0/1）
#   pp_at_risk            ：PP解析のrisk setに含まれる区間か（0/1）
#
# データセットのヘルプファイル：
# ?SGLT2

# 治療群別の人数
table(SGLT2_baseline$treatment)

# person-monthデータ上の死亡イベント数
with(SGLT2, tapply(Y_death, treatment, sum))

# 必要な変数が存在することを確認します。
required_baseline <- c(
  "id", "trial", "A", "treatment", "age", "female", "bmi", "hba1c",
  "egfr", "proteinuria", "prior_heart_failure", "prior_stroke",
  "recent_hospitalization", "trial_period"
)

required_long <- c(
  "id", "trial", "time", "A", "treatment", "Y_death",
  "stay_ltfu", "stay_adherent", "pp_at_risk",
  "age", "female", "bmi", "hba1c", "egfr", "proteinuria",
  "prior_heart_failure", "prior_stroke", "recent_hospitalization",
  "trial_period"
)

setdiff(required_baseline, names(SGLT2_baseline))
setdiff(required_long, names(SGLT2))

stopifnot(all(required_baseline %in% names(SGLT2_baseline)))
stopifnot(all(required_long %in% names(SGLT2)))

# パッケージ内部のデータを直接変更しないよう、作業用にコピーします。
sglt2_baseline <- SGLT2_baseline
sglt2_long <- SGLT2

# Aの符号と表示順を確実に一致させます。
# A = 0をDPP-4i、A = 1をSGLT2iとして表示します。
sglt2_baseline$treatment <- factor(
  sglt2_baseline$A,
  levels = c(0, 1),
  labels = c("DPP-4i", "SGLT2i")
)

sglt2_long$treatment <- factor(
  sglt2_long$A,
  levels = c(0, 1),
  labels = c("DPP-4i", "SGLT2i")
)

# Person-periodデータの構造を確認します。
# この段階では、まだ解析weightを推定していないためweightsは指定しません。
check_sglt2 <- check_tte(
  sglt2_long,
  id = id,
  trial = trial,
  time = time,
  event = Y_death,
  treatment = A
)

print(check_sglt2)


# -----------------------------------------------------------------------------
# 1-2. Baseline treatment weight
# -----------------------------------------------------------------------------
#
# est_wt()は、デフォルトではstabilized weightを推定します。
#
# このbaseline treatment weightは、各人について、
#
#   分子：観察された治療を受ける周辺確率
#   分母：baseline covariatesを条件とした観察治療の確率
#
# の比として計算されます。
#
# この例では、
#   - 年齢、性別
#   - BMI、HbA1c、eGFR、proteinuria
#   - 心不全・脳卒中の既往
#   - 最近の入院
#   - calendar period
# をbaseline confoundersとして使用します。
#
# Componentごとにはtruncationを行わず、最終weightを作った後に
# 1・99 percentileでtruncationします。
#
# est_wt()でtruncate = c(0, 1)とすると、実質的にtruncationなしです。

wt_treatment_sglt2 <- est_wt(
  A ~
    age + I(age^2) +
    female +
    bmi +
    hba1c + I(hba1c^2) +
    egfr +
    proteinuria +
    prior_heart_failure +
    prior_stroke +
    recent_hospitalization +
    trial_period,

  data = sglt2_baseline,
  type = "treatment",
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_treatment_sglt2)
plot(wt_treatment_sglt2)


# -----------------------------------------------------------------------------
# 1-3. Baseline covariate balance
# -----------------------------------------------------------------------------
#
# Baseline IPTWによって、治療群間の共変量バランスが改善したかを
# standardized mean difference（SMD）で確認します。
#
# |SMD| <= 0.10はよく用いられる目安ですが、
# 0.10未満であれば未測定交絡まで解消された、という意味ではありません。

balance_sglt2 <- balance_wt(
  A ~
    age +
    female +
    bmi +
    hba1c +
    egfr +
    proteinuria +
    prior_heart_failure +
    prior_stroke +
    recent_hospitalization +
    trial_period,

  data = sglt2_baseline,
  weights = wt_treatment_sglt2
)

print(balance_sglt2)
summary(balance_sglt2, threshold = 0.10)
plot(balance_sglt2, threshold = 0.10)


# -----------------------------------------------------------------------------
# 1-4. Baseline weightをperson-periodデータへ対応づける
# -----------------------------------------------------------------------------
#
# Baselineデータは1人1行、longデータは1人複数行です。
# 各人のbaseline treatment weightを、その人のすべてのperson-monthへ
# 対応づけます。
#
# この事例ではtrialは実質的に1つですが、一般化可能な書き方として
# idとtrialの組で対応づけます。

key_baseline <- paste(
  sglt2_baseline$id,
  sglt2_baseline$trial,
  sep = ":"
)

key_long <- paste(
  sglt2_long$id,
  sglt2_long$trial,
  sep = ":"
)

baseline_index <- match(key_long, key_baseline)

stopifnot(!anyNA(baseline_index))

sglt2_long$w_treatment_est <- weights(
  wt_treatment_sglt2,
  which = "untruncated"
)[baseline_index]

summary(sglt2_long$w_treatment_est)


# -----------------------------------------------------------------------------
# 1-5. ITT解析用のloss-to-follow-up weight
# -----------------------------------------------------------------------------
#
# stay_ltfu = 1：
#   次のintervalまで追跡下に残る
#
# cumulative = TRUE：
#   interval-specific weight factorを、id・trial内で時間順に累積します。
#
# lag = 1：
#   interval tのoutcome contributionには、interval tの開始時点までに
#   累積したweightを使います。
#
# 分母モデルには、追跡不能とoutcomeの両方に関連し得る変数を含めます。
#
# 分子モデルはweightを安定化するため、
# treatment、follow-up time、calendar periodに限定します。
#
# 【実データへの注意】
# この教育用SGLT2データのclinical covariatesは、主としてbaseline値です。
# 実研究でtime-varying covariatesをIPCWモデルに入れる場合は、
# 必ずその区間の打切り判断より前に測定された値を使用します。

wt_ltfu_itt_sglt2 <- est_wt(
  stay_ltfu ~
    A +
    splines::ns(time, df = 3) +
    age +
    female +
    bmi +
    hba1c +
    egfr +
    proteinuria +
    prior_heart_failure +
    prior_stroke +
    recent_hospitalization +
    trial_period,

  numerator =
    stay_ltfu ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = sglt2_long,
  type = "censoring",
  id = id,
  trial = trial,
  time = time,
  cumulative = TRUE,
  lag = 1,
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_ltfu_itt_sglt2)

# Interval-specific factor、累積weight、解析に用いるlagged weightを確認
head(as.data.frame(wt_ltfu_itt_sglt2), 10)


# -----------------------------------------------------------------------------
# 1-6. Final weight for ITT analysis:
#      Baseline IPTW × LTFU-IPCW
# -----------------------------------------------------------------------------
#
# Baseline treatment weightとloss-to-follow-up weightを掛け合わせた後、
# 最終weightを1・99 percentileでtruncationします。

wt_itt_sglt2 <- combine_wt(
  sglt2_long$w_treatment_est,
  wt_ltfu_itt_sglt2,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

sglt2_long$w_itt_est <- weights(wt_itt_sglt2)

summary(sglt2_long$w_itt_est)

diag_itt_sglt2 <- diagnose_wt(
  w_itt_est,
  data = sglt2_long,
  treatment = treatment,
  time = time,
  id = id,
  trial = trial
)

print(diag_itt_sglt2)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(diag_itt_sglt2, type = "weights")
plot(diag_itt_sglt2, type = "risk_set")
par(old_par)


# -----------------------------------------------------------------------------
# 1-7. ITT effect：weighted pooled discrete-time outcome model
# -----------------------------------------------------------------------------
#
# 1か月単位の区間ハザードを、weighted pooled GLMとしてモデル化します。
#
# ここではquasibinomial familyとcomplementary log-log linkを用います。
# 平均モデルは離散時間比例ハザードモデルに対応し、
# exp(Aの回帰係数)はdiscrete-time hazard ratioとして解釈できます。
#
# splines::ns(time, df = 3)により、baseline hazardの時間変化を
# 自然スプラインで柔軟にモデル化します。
#
# 同一人物が複数のperson-monthに寄与するため、
# 元の個人IDでcluster化したロバスト分散を用います。

fit_death <- discsurvreg(
  Y_death ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = sglt2_long,
  id = id,
  weights = w_itt_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_death)
summary(fit_death)

# SGLT2i versus DPP-4iのhazard ratioと95%信頼区間
confint(fit_death, parm = "A", eform = TRUE)

# 【解釈上の注意】
# Hazard ratioは相対的な区間ハザードの指標です。
# 3年・5年の絶対リスク、risk difference、risk ratioを得るには、
# 次のstandardizationが必要です。


# -----------------------------------------------------------------------------
# 1-8. ITT effect：standardizationによる絶対リスク
# -----------------------------------------------------------------------------
#
# std_tte()は、baseline populationの各人について、
#
#   全員がDPP-4iを開始した世界（A = 0）
#   全員がSGLT2iを開始した世界（A = 1）
#
# のcounterfactual predictionを作り、それを対象集団で平均します。
#
# これにより、
#   - 治療別のsurvival probability
#   - 治療別のcumulative risk
#   - risk difference
#   - risk ratio
#   - time-specific NNT/NNH
# を推定できます。

std_death <- std_tte(
  fit_death,

  # 標準化の対象集団
  data = sglt2_baseline,

  treatment = A,
  time = time,

  # 0～59の60区間＝60か月
  times = 0:59,

  values = c(0, 1),
  labels = c("DPP-4i", "SGLT2i")
)

# 36か月時点
summary(std_death, horizon = 36)

# 60か月時点
summary(std_death, horizon = 60)

# モデルベースの標準化累積リスク曲線
plot(std_death, measure = "risk")

# Survival curveを描く場合
# plot(std_death, measure = "survival")


# -----------------------------------------------------------------------------
# 1-9. ITT effect：weighted Kaplan-Meier curve
# -----------------------------------------------------------------------------
#
# curve_tte(type = "km")は、IPW pseudo-populationの各risk setから、
# Kaplan-Meier型のsurvival/risk curveを直接推定します。
#
# 1-8の標準化曲線：
#   outcome modelを用いるmodel-based estimate
#
# Weighted Kaplan-Meier curve：
#   各risk setの重み付きイベント数から求める、よりnonparametricな推定
#
# 両者が完全に一致する必要はありませんが、大きく異なる場合には、
# outcome model、weight model、positivity、time function、
# follow-up後半のrisk-set sizeなどを再検討します。

km_death <- curve_tte(
  sglt2_long,
  time = time,
  event = Y_death,
  treatment = treatment,
  weights = w_itt_est,
  type = "km",
  id = id,
  trial = trial
)

# 利用者向けには、経過月をそのまま指定します。
summary(km_death, time = 36)
summary(km_death, time = 60)

# Kaplan-Meier推定量は階段関数として描画されます。
plot(km_death, measure = "risk")

# Survival curveを描く場合
# plot(km_death, measure = "survival")


# -----------------------------------------------------------------------------
# 1-10. Per-protocol解析用データ
# -----------------------------------------------------------------------------
#
# ITT estimand：
#   baselineで開始したstrategyに割り付けたままとみなし、
#   その後の中止・switchingにかかわらず比較します。
#
# PP estimand：
#   assigned strategyから逸脱した時点でartificial censoringし、
#   その選択をIPCWで補正します。
#
# pp_at_risk == 1のintervalだけをPP解析に使用します。

sglt2_pp <- subset(
  sglt2_long,
  pp_at_risk == 1
)

# 次区間まで、
#   - 自然な追跡不能がなく
#   - assigned strategyからの逸脱もない
# 場合を1とするjoint censoring indicatorを作ります。
sglt2_pp$stay_pp <- as.integer(
  sglt2_pp$stay_ltfu == 1 &
  sglt2_pp$stay_adherent == 1
)

table(sglt2_pp$stay_pp)


# -----------------------------------------------------------------------------
# 1-11. PP解析用のjoint censoring weight
# -----------------------------------------------------------------------------
#
# Loss to follow-upとprotocol deviationをまとめて、
# joint censoring processとしてモデル化します。
#
# 実研究では、治療継続、中止、switchingを決めるtime-varying predictorsを
# 分母モデルに十分含めることが重要です。
#
# この教育用データでは利用できるtime-varying clinical covariatesが
# 限られているため、PP weight modelは簡略化された教材上の例です。
#
# 実研究では、HbA1c、eGFR、入院、併用薬、frailtyなどの更新値が
# protocol deviationとoutcomeに影響する場合、それらを適切なlagで
# 分母モデルに含める必要があります。

wt_censor_pp_sglt2 <- est_wt(
  stay_pp ~
    A +
    splines::ns(time, df = 3) +
    age +
    female +
    bmi +
    hba1c +
    egfr +
    proteinuria +
    prior_heart_failure +
    prior_stroke +
    recent_hospitalization +
    trial_period,

  numerator =
    stay_pp ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = sglt2_pp,
  type = "censoring",
  id = id,
  trial = trial,
  time = time,
  cumulative = TRUE,
  lag = 1,
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_censor_pp_sglt2)


# -----------------------------------------------------------------------------
# 1-12. Final weight for PP analysis:
#       Baseline IPTW × IPCW for joint censoring
# -----------------------------------------------------------------------------

wt_pp_sglt2 <- combine_wt(
  sglt2_pp$w_treatment_est,
  wt_censor_pp_sglt2,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

sglt2_pp$w_pp_est <- weights(wt_pp_sglt2)

summary(sglt2_pp$w_pp_est)

diag_pp_sglt2 <- diagnose_wt(
  w_pp_est,
  data = sglt2_pp,
  treatment = treatment,
  time = time,
  id = id,
  trial = trial
)

print(diag_pp_sglt2)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(diag_pp_sglt2, type = "weights")
plot(diag_pp_sglt2, type = "risk_set")
par(old_par)


# -----------------------------------------------------------------------------
# 1-13. PP effect：weighted pooled discrete-time outcome model
# -----------------------------------------------------------------------------

fit_death_pp <- discsurvreg(
  Y_death ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = sglt2_pp,
  id = id,
  weights = w_pp_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_death_pp)
summary(fit_death_pp)

# SGLT2i versus DPP-4iのPP hazard ratioと95%信頼区間
confint(fit_death_pp, parm = "A", eform = TRUE)


# -----------------------------------------------------------------------------
# 1-14. PP effect：standardizationによる絶対リスク
# -----------------------------------------------------------------------------
#
# PP outcome modelを用いて、全員が各strategyを継続した場合の
# 標準化累積リスクを推定します。

std_death_pp <- std_tte(
  fit_death_pp,
  data = sglt2_baseline,
  treatment = A,
  time = time,
  times = 0:59,
  values = c(0, 1),
  labels = c("DPP-4i", "SGLT2i")
)

summary(std_death_pp, horizon = 36)
summary(std_death_pp, horizon = 60)

plot(std_death_pp, measure = "risk")


# -----------------------------------------------------------------------------
# 1-15. PP effect：weighted Kaplan-Meier curve
# -----------------------------------------------------------------------------
#
# PP risk setと推定したPP weightを使い、
# PP estimandに対応するweighted Kaplan-Meier curveを直接推定します。

km_death_pp <- curve_tte(
  sglt2_pp,
  time = time,
  event = Y_death,
  treatment = treatment,
  weights = w_pp_est,
  type = "km",
  id = id,
  trial = trial
)

summary(km_death_pp, time = 36)
summary(km_death_pp, time = 60)

plot(km_death_pp, measure = "risk")


# -----------------------------------------------------------------------------
# 1-16. 収載済み参照weightとの比較（任意）
# -----------------------------------------------------------------------------
#
# SGLT2データには、合成データ生成過程で用意した参照用のw_ittとw_ppも
# 含まれています。
#
# 今回推定したw_itt_estとw_pp_estは、観測された合成データに
# probability modelを当てはめて推定したweightなので、
# 参照weightと完全に一致する必要はありません。
#
# 必要に応じて、以下を実行して比較できます。

# cor(
#   log(sglt2_long$w_itt),
#   log(sglt2_long$w_itt_est)
# )

# cor(
#   log(sglt2_pp$w_pp),
#   log(sglt2_pp$w_pp_est)
# )


# =============================================================================
# 解析上の重要な整理
# =============================================================================
#
# 1. fit_death
#    Baseline treatmentとLTFUを補正したITT hazard ratioを推定
#
# 2. std_death
#    ITT outcome modelをstandardizeし、
#    治療別絶対リスク、risk difference、risk ratioを推定
#
# 3. km_death
#    ITT pseudo-population上のweighted Kaplan-Meier curveを直接推定
#
# 4. fit_death_pp
#    Treatment strategyからの逸脱を人工的に打ち切り、
#    joint censoring IPCWで補正したPP hazard ratioを推定
#
# 5. std_death_pp
#    PP outcome modelをstandardizeし、PP絶対リスクを推定
#
# 6. km_death_pp
#    PP risk setとPP weightによるweighted Kaplan-Meier curveを直接推定
#
# この全死亡事例にはcompeting eventがないため、
# ARB/CCB事例のようなAalen-Johansen推定は必要ありません。
# =============================================================================


# -----------------------------------------------------------------------------
# 解析環境の記録
# -----------------------------------------------------------------------------

sessionInfo()

# =============================================================================
# End of SGLT2i versus DPP-4i tutorial
# =============================================================================
