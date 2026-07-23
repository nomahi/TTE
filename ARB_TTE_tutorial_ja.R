# =============================================================================
# 2. Target Trial Emulation for ARB versus CCB:
#    Heart-failure hospitalization with death as a competing event
# =============================================================================
#
# モデルとした研究：
# Noma, H., Kurita, N., Fukuma, S., Fujisawa, T., Oda, F., Maeda, M.,
# and Fukuda, H. (2026).
# Heart failure and renal outcomes with angiotensin receptor blockers
# compared with calcium channel blockers in patients with chronic kidney
# disease: a target trial emulation.
# Heart. DOI: 10.1136/heartjnl-2026-328193.
#
# どんな研究？
# - CKD患者において、ARB-based strategyとCCB-based strategyを比較
# - 複数のcalendar timeをtrial baselineとするsequentially nested trials
# - primary outcomeは心不全入院
# - 全死亡はcompeting event
# - 1か月のinduction periodを設定
# - ITT effectとper-protocol（PP）effectの両方を推定
#
# コンピュータ乱数により、この研究を模した小規模な合成データが
# TTEパッケージに収載されています。
#
#   750人から構成される900件のperson-trial entry
#   ARB：370 person-trial entries
#   CCB：530 person-trial entries
#
# 【重要】
# - 同じ人が複数のsequential trialに参加することがあります。
# - したがって、ARB_baselineの1行は「1人」ではなく、
#   「1つのperson-trial entry」を表します。
# - 収載データはすべて教育用の合成データです。
# - 実在する患者記録や元研究の実データは含まれていません。
#
# この事例で学ぶ内容：
# - baseline IPTW
# - loss-to-follow-upに対する累積IPCW
# - ITT解析用の最終weight
# - cause-specific discrete-time hazard model
# - competing riskを考慮したstandardized cumulative incidence function
# - weighted Aalen-Johansen curve
# - protocol deviationを考慮したPP解析
# =============================================================================


# -----------------------------------------------------------------------------
# 2-1. データの読み込みと構造の確認
# -----------------------------------------------------------------------------

data(ARB_baseline)
data(ARB)

# Baselineデータ：
#   1行は1つのperson-trial entryです。
#   同一人物が異なるtrial baselineで複数回適格となることがあります。
#   主として、baseline treatment weightと共変量バランスの評価に使います。
dim(ARB_baseline)
length(unique(ARB_baseline$id))
head(ARB_baseline, 20)

# Longデータ：
#   1つのperson-trial entryが複数のperson-monthに展開されています。
#   pooled discrete-time modelやweighted survival curveに使います。
dim(ARB)
head(ARB, 6)

# 主要変数
#   id          ：元の個人ID
#   trial       ：emulated trialの識別子
#   time        ：trial baselineからの経過月（0始まり）
#   A           ：治療指標（0 = CCB、1 = ARB）
#   treatment   ：治療名を表すfactor
#   Y_hf        ：その区間に心不全入院が発生したか（0/1）
#   Y_death     ：その区間に死亡が発生したか（0/1）
#   event_code  ：0 = イベントなし、1 = 心不全入院、2 = 死亡
#   stay_ltfu   ：次区間まで追跡下に残ったか（0/1）
#   stay_adherent：次区間までassigned strategyに従ったか（0/1）
#   pp_at_risk  ：PP解析のrisk setに含まれる区間か（0/1）
#
# データセットのヘルプファイル：
# ?ARB

# 治療群別のperson-trial entry数
table(ARB_baseline$treatment)

# 元の個人数とperson-trial entry数
length(unique(ARB_baseline$id))
nrow(ARB_baseline)

# 治療群別のprimary eventとcompeting eventの件数
with(ARB, table(treatment, event_code))

# パッケージ内部のデータを直接変更しないよう、作業用にコピーします。
arb_baseline <- ARB_baseline
arb_long <- ARB

# Aの符号と表示順を確実に一致させます。
# A = 0をCCB、A = 1をARBとして表示します。
arb_baseline$treatment <- factor(
  arb_baseline$A,
  levels = c(0, 1),
  labels = c("CCB", "ARB")
)

arb_long$treatment <- factor(
  arb_long$A,
  levels = c(0, 1),
  labels = c("CCB", "ARB")
)

# Person-periodデータの構造を確認します。
# この段階では、まだ解析weightを推定していないためweightsは指定しません。
check_arb <- check_tte(
  arb_long,
  id = id,
  trial = trial,
  time = time,
  event = event_code,
  treatment = A
)

print(check_arb)


# -----------------------------------------------------------------------------
# 2-2. Baseline treatment weight
# -----------------------------------------------------------------------------
#
# est_wt()は、デフォルトではstabilized weightを推定します。
#
# このbaseline treatment weightは、各person-trial entryについて、
#
#   分子：観察された治療を受ける周辺確率
#   分母：baseline covariatesを条件とした観察治療の確率
#
# の比として計算されます。
#
# この例では、
#   - 年齢、性別
#   - BMI、血圧
#   - eGFR、proteinuria、diabetes
#   - 既往歴
#   - calendar period
# をbaseline confoundersとして使用します。
#
# Componentごとにはtruncationを行わず、最終weightを作った後に
# 1・99 percentileでtruncationします。
#
# est_wt()でtruncate = c(0, 1)とすると、実質的にtruncationなしです。

wt_treatment_arb <- est_wt(
  A ~
    age + I(age^2) +
    female +
    bmi +
    sbp +
    dbp +
    egfr +
    proteinuria +
    diabetes +
    prior_heart_failure +
    prior_stroke +
    trial_period,

  data = arb_baseline,
  type = "treatment",
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_treatment_arb)
plot(wt_treatment_arb)


# -----------------------------------------------------------------------------
# 2-3. Baseline covariate balance
# -----------------------------------------------------------------------------
#
# Baseline IPTWによって治療群間の共変量バランスが改善したかを、
# standardized mean difference（SMD）で確認します。
#
# |SMD| <= 0.10はよく使われる目安ですが、
# 0.10未満であれば未測定交絡まで解消された、という意味ではありません。

balance_arb <- balance_wt(
  A ~
    age +
    female +
    bmi +
    sbp +
    dbp +
    egfr +
    proteinuria +
    diabetes +
    prior_heart_failure +
    prior_stroke +
    trial_period,

  data = arb_baseline,
  weights = wt_treatment_arb
)

print(balance_arb)
summary(balance_arb, threshold = 0.10)
plot(balance_arb, threshold = 0.10)


# -----------------------------------------------------------------------------
# 2-4. Baseline weightをperson-periodデータへ対応づける
# -----------------------------------------------------------------------------
#
# Sequential trial designでは同一人物が複数trialに参加するため、
# idだけではなく、idとtrialの組でbaseline rowとlong rowを対応づけます。

key_baseline <- paste(
  arb_baseline$id,
  arb_baseline$trial,
  sep = ":"
)

key_long <- paste(
  arb_long$id,
  arb_long$trial,
  sep = ":"
)

baseline_index <- match(key_long, key_baseline)

stopifnot(!anyNA(baseline_index))

arb_long$w_treatment_est <- weights(
  wt_treatment_arb,
  which = "untruncated"
)[baseline_index]

summary(arb_long$w_treatment_est)


# -----------------------------------------------------------------------------
# 2-5. ITT解析用のloss-to-follow-up weight
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
# ARBデータでは、sbpとegfrはlong data内で時間変化します。
# 実データでtime-varying variablesを使用するときは、必ず
# 「その区間の打切り判断より前に測定された値」を使う必要があります。
#
# 分子モデルはweightを安定化するため、
# treatment、follow-up time、calendar periodに限定します。

wt_ltfu_itt_arb <- est_wt(
  stay_ltfu ~
    A +
    splines::ns(time, df = 3) +
    age +
    female +
    bmi +
    sbp +
    dbp +
    egfr +
    proteinuria +
    diabetes +
    prior_heart_failure +
    prior_stroke +
    trial_period,

  numerator =
    stay_ltfu ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_long,
  type = "censoring",
  id = id,
  trial = trial,
  time = time,
  cumulative = TRUE,
  lag = 1,
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_ltfu_itt_arb)

# Interval-specific factor、累積weight、解析に使用するlagged weightを確認
head(as.data.frame(wt_ltfu_itt_arb), 10)


# -----------------------------------------------------------------------------
# 2-6. Final weight for ITT analysis:
#      Baseline IPTW × LTFU-IPCW
# -----------------------------------------------------------------------------
#
# Baseline treatment weightとloss-to-follow-up weightを掛け合わせた後、
# 最終weightを1・99 percentileでtruncationします。

wt_itt_arb <- combine_wt(
  arb_long$w_treatment_est,
  wt_ltfu_itt_arb,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

arb_long$w_itt_est <- weights(wt_itt_arb)

summary(arb_long$w_itt_est)

diag_itt_arb <- diagnose_wt(
  w_itt_est,
  data = arb_long,
  treatment = treatment,
  time = time,
  id = id,
  trial = trial
)

print(diag_itt_arb)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(diag_itt_arb, type = "weights")
plot(diag_itt_arb, type = "risk_set")
par(old_par)


# -----------------------------------------------------------------------------
# 2-7. ITT effect：心不全入院のcause-specific outcome model
# -----------------------------------------------------------------------------
#
# 1か月単位のcause-specific hazardを、weighted pooled GLM
# （binomial、complementary log-log link）でモデル化します。
#
# 心不全入院モデルでは、死亡した時点で心不全入院のrisk setから外れます。
#
# exp(Aの回帰係数)は、心不全入院に対するcause-specific
# discrete-time hazard ratioです。
#
# これはFine-Gray modelのsubdistribution hazard ratioではありません。
# Cause-specific HRとcumulative incidenceは異なる量なので、
# 絶対リスクは後のstandardizationによって推定します。
#
# 同一人物が複数のtrialとperson-monthに寄与するため、
# 元の個人IDでcluster化したロバスト分散を用います。

fit_hf <- discsurvreg(
  Y_hf ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_long,
  id = id,
  weights = w_itt_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_hf)
summary(fit_hf)
confint(fit_hf, parm = "A", eform = TRUE)


# -----------------------------------------------------------------------------
# 2-8. ITT effect：競合死亡のcause-specific outcome model
# -----------------------------------------------------------------------------
#
# 心不全入院のcumulative incidenceを推定するには、
# 心不全入院モデルだけでなく、競合死亡のhazard modelも必要です。
#
# exp(Aの回帰係数)は、死亡に対するcause-specific
# discrete-time hazard ratioです。

fit_death_competing <- discsurvreg(
  Y_death ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_long,
  id = id,
  weights = w_itt_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_death_competing)
summary(fit_death_competing)
confint(fit_death_competing, parm = "A", eform = TRUE)


# -----------------------------------------------------------------------------
# 2-9. Competing riskを考慮したstandardized CIF
# -----------------------------------------------------------------------------
#
# std_tte()に、
#
#   fit             ：心不全入院のcause-specific model
#   competing_fit   ：死亡のcause-specific model
#
# を渡します。
#
# 各baseline person-trial entryについて、
#
#   全員がCCBを開始した世界（A = 0）
#   全員がARBを開始した世界（A = 1）
#
# の区間別心不全入院確率と死亡確率を予測し、
# Aalen-Johansen型の再帰計算で心不全入院の
# cumulative incidence function（CIF）を求めます。
#
# competing riskがある場合、
# 心不全入院モデルだけから1 - survivalを計算してはいけません。
# それでは死亡を単純なcensoringとして扱い、一般に心不全入院の
# cumulative incidenceを過大評価します。

std_hf <- std_tte(
  fit_hf,

  # 標準化の対象集団
  data = arb_baseline,

  treatment = A,
  time = time,

  # 0～59の60区間＝60か月
  times = 0:59,

  competing_fit = fit_death_competing,
  values = c(0, 1),
  labels = c("CCB", "ARB")
)

# 36か月時点の標準化CIF、risk difference、risk ratio、NNT/NNH
summary(std_hf, horizon = 36)

# 60か月時点
summary(std_hf, horizon = 60)

# モデルベースの標準化CIF
plot(std_hf, measure = "risk")


# -----------------------------------------------------------------------------
# 2-10. Weighted Aalen-Johansen curve
# -----------------------------------------------------------------------------
#
# event_code:
#   0 = eventなし
#   1 = 心不全入院
#   2 = competing death
#
# type = "aj"、cause = 1を指定すると、
# IPW pseudo-population上で心不全入院のweighted CIFを直接推定します。
#
# 2-9のstandardized CIFはcause-specific outcome modelsに基づく
# model-based estimateです。
#
# Weighted Aalen-Johansen curveは、各risk setの重み付きイベント数から
# 直接求める、よりnonparametricな推定です。
#
# 両者が完全に一致する必要はありませんが、大きく異なる場合は、
# outcome model、weight model、positivity、time function、
# follow-up後半のrisk-set sizeなどを再検討します。

aj_hf <- curve_tte(
  arb_long,
  time = time,
  event = event_code,
  treatment = treatment,
  weights = w_itt_est,
  type = "aj",
  cause = 1,
  id = id,
  trial = trial
)

summary(aj_hf, time = 60)
plot(aj_hf, measure = "risk")


# -----------------------------------------------------------------------------
# 2-11. Per-protocol解析用データ
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

arb_pp <- subset(
  arb_long,
  pp_at_risk == 1
)

# 次区間まで、
#   - 自然な追跡不能がなく
#   - assigned strategyからの逸脱もない
# 場合を1とするjoint censoring indicatorを作ります。
arb_pp$stay_pp <- as.integer(
  arb_pp$stay_ltfu == 1 &
  arb_pp$stay_adherent == 1
)

table(arb_pp$stay_pp)


# -----------------------------------------------------------------------------
# 2-12. PP解析用のjoint censoring weight
# -----------------------------------------------------------------------------
#
# loss to follow-upとprotocol deviationをまとめて、
# joint censoring processとしてモデル化します。
#
# 実研究では、治療継続、中止、switchingを決めるtime-varying predictorsを
# 分母モデルに十分含めることが重要です。
#
# この合成データでは、time-varying sbpとegfrを使用します。

wt_censor_pp_arb <- est_wt(
  stay_pp ~
    A +
    splines::ns(time, df = 3) +
    age +
    female +
    bmi +
    sbp +
    dbp +
    egfr +
    proteinuria +
    diabetes +
    prior_heart_failure +
    prior_stroke +
    trial_period,

  numerator =
    stay_pp ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_pp,
  type = "censoring",
  id = id,
  trial = trial,
  time = time,
  cumulative = TRUE,
  lag = 1,
  stabilize = TRUE,
  truncate = c(0, 1)
)

summary(wt_censor_pp_arb)


# -----------------------------------------------------------------------------
# 2-13. Final weight for PP analysis:
#       Baseline IPTW × IPCW for joint censoring
# -----------------------------------------------------------------------------

wt_pp_arb <- combine_wt(
  arb_pp$w_treatment_est,
  wt_censor_pp_arb,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

arb_pp$w_pp_est <- weights(wt_pp_arb)

summary(arb_pp$w_pp_est)

diag_pp_arb <- diagnose_wt(
  w_pp_est,
  data = arb_pp,
  treatment = treatment,
  time = time,
  id = id,
  trial = trial
)

print(diag_pp_arb)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(diag_pp_arb, type = "weights")
plot(diag_pp_arb, type = "risk_set")
par(old_par)


# -----------------------------------------------------------------------------
# 2-14. PP effect：心不全入院のcause-specific outcome model
# -----------------------------------------------------------------------------

fit_hf_pp <- discsurvreg(
  Y_hf ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_pp,
  id = id,
  weights = w_pp_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_hf_pp)
summary(fit_hf_pp)
confint(fit_hf_pp, parm = "A", eform = TRUE)


# -----------------------------------------------------------------------------
# 2-15. PP effect：競合死亡のcause-specific outcome model
# -----------------------------------------------------------------------------
#
# PPの標準化CIFを推定する場合も、心不全入院と死亡の両方について、
# 同じPP risk setと同じw_pp_estを用いてモデルを推定します。

fit_death_pp <- discsurvreg(
  Y_death ~
    A +
    splines::ns(time, df = 3) +
    trial_period,

  data = arb_pp,
  id = id,
  weights = w_pp_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

print(fit_death_pp)
summary(fit_death_pp)
confint(fit_death_pp, parm = "A", eform = TRUE)


# -----------------------------------------------------------------------------
# 2-16. PP effect：competing riskを考慮したstandardized CIF
# -----------------------------------------------------------------------------

std_hf_pp <- std_tte(
  fit_hf_pp,
  data = arb_baseline,
  treatment = A,
  time = time,
  times = 0:59,
  competing_fit = fit_death_pp,
  values = c(0, 1),
  labels = c("CCB", "ARB")
)

summary(std_hf_pp, horizon = 36)
summary(std_hf_pp, horizon = 60)

plot(std_hf_pp, measure = "risk")


# -----------------------------------------------------------------------------
# 2-17. PP effect：Weighted Aalen-Johansen curve
# -----------------------------------------------------------------------------
#
# PP risk setと推定したPP weightを用いて、
# 心不全入院のweighted CIFを直接推定します。

aj_hf_pp <- curve_tte(
  arb_pp,
  time = time,
  event = event_code,
  treatment = treatment,
  weights = w_pp_est,
  type = "aj",
  cause = 1,
  id = id,
  trial = trial
)

summary(aj_hf_pp, time = 60)
plot(aj_hf_pp, measure = "risk")


# -----------------------------------------------------------------------------
# 2-18. 収載済み参照weightとの比較（任意）
# -----------------------------------------------------------------------------
#
# ARBデータには、合成データ生成過程で用意した参照用のw_ittとw_ppも
# 含まれています。
#
# 今回推定したw_itt_estとw_pp_estは、観測された合成データに
# probability modelを当てはめて推定したweightなので、
# 参照weightと完全に一致する必要はありません。
#
# 必要に応じて、以下を実行して比較できます。

# cor(
#   log(arb_long$w_itt),
#   log(arb_long$w_itt_est)
# )

# cor(
#   log(arb_pp$w_pp),
#   log(arb_pp$w_pp_est)
# )


# =============================================================================
# 解析上の重要な整理
# =============================================================================
#
# 1. fit_hf / fit_hf_pp
#    心不全入院に対するcause-specific hazard ratioを推定
#
# 2. fit_death_competing / fit_death_pp
#    競合死亡に対するcause-specific hazard ratioを推定
#
# 3. std_hf / std_hf_pp
#    2つのcause-specific modelsを組み合わせて、
#    心不全入院の標準化CIF、risk difference、risk ratioを推定
#
# 4. aj_hf / aj_hf_pp
#    weighted Aalen-Johansen estimatorによってCIFを直接推定
#
# 競合リスク下では、
#   1 - Kaplan-Meier
# を心不全入院のcumulative incidenceとして用いてはいけません。
# =============================================================================
