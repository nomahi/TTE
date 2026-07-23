# =============================================================================
# TTE package: hands-on tutorial
#
# 2. Target trial emulation for ARB versus CCB strategies
#    Outcome: heart-failure hospitalization
#    Competing event: all-cause death
# =============================================================================
#
# Motivating study:
# Noma H, Kurita N, Fukuma S, Fujisawa T, Oda F, Maeda M, Fukuda H.
# Heart failure and renal outcomes with angiotensin receptor blockers
# compared with calcium channel blockers in patients with chronic kidney
# disease: a target trial emulation.
# Heart. 2026. doi:10.1136/heartjnl-2026-328193.
#
# Design represented by the synthetic example:
# - patients with chronic kidney disease
# - ARB-based versus CCB-based treatment strategies
# - sequentially nested new-user trials across calendar time
# - primary outcome: heart-failure hospitalization
# - competing event: all-cause death
# - one-month induction period
# - intention-to-treat (ITT) and per-protocol (PP) estimands
#
# The package contains a small synthetic dataset:
#   750 unique individuals
#   900 person-trial entries
#   ARB: 370 person-trial entries
#   CCB: 530 person-trial entries
#
# IMPORTANT:
# - A person may enter more than one sequential trial.
# - One row of `ARB_baseline` therefore represents a person-trial entry, not
#   necessarily a unique person.
# - These data are entirely synthetic and contain no actual patient records.
# - Numerical results must not be interpreted as estimates from the motivating
#   clinical study.
#
# This tutorial estimates all analysis weights from the observed synthetic data.
# It emphasizes cause-specific hazards, standardized cumulative incidence
# functions, and weighted Aalen-Johansen curves.
#
# TTE version 1.1.1 or later is assumed.
# =============================================================================


# Install once after the package becomes available on CRAN:
# install.packages("TTE")

library("TTE")

options(width = 110)

packageVersion("TTE")


# -----------------------------------------------------------------------------
# 2-1. Load and inspect the sequential-trial data
# -----------------------------------------------------------------------------
#
# `ARB_baseline` contains one row per person-trial entry. `ARB` expands each
# entry into person-month records. The same original person may contribute to
# multiple emulated trials.
# The code verifies the expected variables, creates working copies, fixes the
# display order as CCB followed by ARB, and checks the person-period structure
# before estimating any weights.

data(ARB_baseline)
data(ARB)

dim(ARB_baseline)
length(unique(ARB_baseline$id))
head(ARB_baseline, 20)

dim(ARB)
head(ARB, 6)


table(ARB_baseline$treatment)

length(unique(ARB_baseline$id))
nrow(ARB_baseline)

with(ARB, table(treatment, event_code))

arb_baseline <- ARB_baseline
arb_long <- ARB

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
# 2-2. Estimate the baseline treatment weight
# -----------------------------------------------------------------------------
#
# The stabilized baseline treatment weight compares the marginal
# probability of the observed strategy with the probability conditional on
# measured baseline covariates.
# The model includes demographics, BMI, blood pressure, kidney function,
# proteinuria, diabetes, prior cardiovascular disease, and calendar period.
# Component weights are left untruncated until the final ITT or PP weight is
# constructed.

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
# 2-3. Assess baseline covariate balance
# -----------------------------------------------------------------------------
#
# `balance_wt()` reports unweighted and weighted standardized mean
# differences. Balance is assessed across person-trial entries, which is the
# unit randomized in each emulated sequential trial.

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
# 2-4. Map the baseline weight to person-month records
# -----------------------------------------------------------------------------
#
# Because an individual may enter several trials, matching must use the
# combination of `id` and `trial`. Each baseline treatment weight is then copied
# to all person-months belonging to that person-trial entry.

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
# 2-5. Estimate the loss-to-follow-up weight for the ITT analysis
# -----------------------------------------------------------------------------
#
# Interval-specific censoring probabilities are estimated and accumulated
# within each `id`-`trial` trajectory. `lag = 1` assigns to each outcome interval
# the weight available at the start of that interval.
# In this synthetic dataset, systolic blood pressure and eGFR vary over time.
# In real data, any time-varying covariate must be measured before the censoring
# decision it predicts.

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

head(as.data.frame(wt_ltfu_itt_arb), 10)


# -----------------------------------------------------------------------------
# 2-6. Construct and diagnose the final ITT weight
# -----------------------------------------------------------------------------
#
# The baseline treatment weight and cumulative loss-to-follow-up weight are
# multiplied, then truncated at the 1st and 99th percentiles.
# Diagnostics summarize the weight distribution, effective sample size,
# treatment-specific behavior, and weighted risk sets over follow-up.

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
# 2-7. Fit the cause-specific model for heart-failure hospitalization
# -----------------------------------------------------------------------------
#
# The weighted complementary log-log model estimates the cause-specific
# discrete-time hazard ratio for heart-failure hospitalization. Death removes an
# individual from the heart-failure risk set.
# This is not a Fine-Gray subdistribution hazard model. The cause-specific hazard
# ratio and cumulative incidence answer different questions.

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
# 2-8. Fit the cause-specific model for competing death
# -----------------------------------------------------------------------------
#
# A second weighted outcome model estimates the cause-specific hazard of
# death. Both cause-specific models are required to obtain a model-based
# cumulative incidence function for heart-failure hospitalization.

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
# 2-9. Standardize the cause-specific models to a cumulative incidence function
# -----------------------------------------------------------------------------
#
# For every baseline person-trial entry, `std_tte()` predicts the
# interval-specific probabilities of heart-failure hospitalization and death
# under CCB and under ARB. An Aalen-Johansen-type recursion combines the two
# causes to obtain the standardized cumulative incidence function (CIF).
# Using one minus Kaplan-Meier while treating death as ordinary censoring would
# generally overestimate the cumulative incidence of heart failure.

std_hf <- std_tte(
  fit_hf,

  data = arb_baseline,

  treatment = A,
  time = time,

  times = 0:59,

  competing_fit = fit_death_competing,
  values = c(0, 1),
  labels = c("CCB", "ARB")
)

summary(std_hf, horizon = 36)

summary(std_hf, horizon = 60)

plot(std_hf, measure = "risk")


# -----------------------------------------------------------------------------
# 2-10. Estimate an IPW-weighted Aalen-Johansen curve
# -----------------------------------------------------------------------------
#
# The standardized CIF is model based. The weighted Aalen-Johansen curve is
# estimated directly from weighted risk sets and cause-specific event counts.
# The two estimates need not be identical. Large differences can point to
# outcome-model misspecification, extreme weights, limited positivity, or a
# small weighted risk set late in follow-up.

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
# 2-11. Create the per-protocol risk set
# -----------------------------------------------------------------------------
#
# For the PP estimand, follow-up is artificially censored at protocol
# deviation. Intervals with `pp_at_risk == 1` are retained, and `stay_pp`
# indicates remaining both observed and adherent through the next interval.

arb_pp <- subset(
  arb_long,
  pp_at_risk == 1
)

arb_pp$stay_pp <- as.integer(
  arb_pp$stay_ltfu == 1 &
  arb_pp$stay_adherent == 1
)

table(arb_pp$stay_pp)


# -----------------------------------------------------------------------------
# 2-12. Estimate the joint censoring weight for the PP analysis
# -----------------------------------------------------------------------------
#
# Loss to follow-up and protocol deviation are modeled jointly.
# The denominator includes measured predictors of treatment continuation,
# switching, discontinuation, and the outcome. The time-varying systolic blood
# pressure and eGFR values in this synthetic dataset illustrate how updated
# predictors can enter the censoring model.

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
# 2-13. Construct and diagnose the final PP weight
# -----------------------------------------------------------------------------
#
# The baseline treatment weight is multiplied by the cumulative joint
# censoring weight, and the final result is truncated at the 1st and 99th
# percentiles.

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
# 2-14. Fit the PP cause-specific model for heart-failure hospitalization
# -----------------------------------------------------------------------------
#
# The heart-failure model is fitted in the PP risk set using `w_pp_est`.
# Its treatment coefficient estimates the cause-specific hazard ratio under
# continued adherence, subject to the required exchangeability and positivity
# assumptions.

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
# 2-15. Fit the PP cause-specific model for competing death
# -----------------------------------------------------------------------------
#
# The competing-death model must use the same PP risk set and PP weights as
# the heart-failure model. This is necessary for a coherent PP cumulative
# incidence calculation.

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
# 2-16. Standardize the PP cause-specific models to a cumulative incidence function
# -----------------------------------------------------------------------------
#
# The two PP cause-specific models are combined to estimate the
# heart-failure CIF under sustained CCB and sustained ARB strategies.

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
# 2-17. Estimate a weighted Aalen-Johansen curve for the PP effect
# -----------------------------------------------------------------------------
#
# The PP Aalen-Johansen curve is estimated directly from the PP risk set and
# the estimated PP weights, providing a less outcome-model-dependent comparison
# with the standardized PP CIF.

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
# 2-18. Optionally compare estimated and bundled reference weights
# -----------------------------------------------------------------------------
#
# The synthetic data also contain reference weights produced during data
# generation. They are not used in the main analysis. Estimated finite-sample
# weights are not expected to equal those reference values exactly.

# Optional comparisons:
# cor(
#   log(arb_long$w_itt),
#   log(arb_long$w_itt_est)
# )
#
# cor(
#   log(arb_pp$w_pp),
#   log(arb_pp$w_pp_est)
# )


# =============================================================================
# Interpretation map
# =============================================================================
#
# fit_hf / fit_hf_pp:
#   Cause-specific hazard ratios for heart-failure hospitalization.
#
# fit_death_competing / fit_death_pp:
#   Cause-specific hazard ratios for the competing event of death.
#
# std_hf / std_hf_pp:
#   Model-based standardized heart-failure CIFs, risk differences, and risk
#   ratios obtained by combining the two cause-specific models.
#
# aj_hf / aj_hf_pp:
#   Direct weighted Aalen-Johansen estimates of the heart-failure CIF.
#
# In the presence of competing death, one minus Kaplan-Meier must not be used
# as the cumulative incidence of heart-failure hospitalization.
# =============================================================================

sessionInfo()

# =============================================================================
# End of the ARB versus CCB tutorial
# =============================================================================
