# TTE

**Transparent Tools for Target Trial Emulation**

[![R-CMD-check](https://github.com/nomahi/TTE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nomahi/TTE/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/TTE)](https://CRAN.R-project.org/package=TTE)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

`TTE` is an R package for transparent, teaching-oriented target trial
emulation with longitudinal observational data. It exposes the major analytic
steps—data checking, trial construction, inverse probability weighting,
diagnostics, weighted discrete-time outcome modeling, standardization, and
weighted survival or cumulative-incidence estimation—rather than hiding the
analysis inside a single black-box function.

Package website: <https://nomahi.github.io/TTE/>

## Installation

Install the CRAN release with:

```r
install.packages("TTE")
```

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("nomahi/TTE")
```

## Main workflow

| Task | Function |
|---|---|
| Check a person-period dataset | `check_tte()` |
| Construct sequentially nested trials | `seqdesign_tte()` |
| Estimate treatment, censoring, or adherence weights | `est_wt()` |
| Combine weight components | `combine_wt()` |
| Assess baseline covariate balance | `balance_wt()` |
| Diagnose weights and effective sample size | `diagnose_wt()` |
| Fit a weighted discrete-time outcome model | `discsurvreg()` |
| Standardize to marginal risks and contrasts | `std_tte()` |
| Estimate weighted Kaplan-Meier or Aalen-Johansen curves | `curve_tte()` |
| Run an individual-cluster bootstrap | `boot_tte()` |

## Synthetic teaching datasets

The package contains two fully synthetic examples. They include no actual
study records.

- `SGLT2_baseline` and `SGLT2`: an active-comparator new-user example
  comparing initiation of an SGLT2 inhibitor with initiation of a DPP-4
  inhibitor among older adults with type 2 diabetes. The outcome is all-cause
  death.
- `ARB_baseline` and `ARB`: a sequentially nested trial example comparing
  ARB-based and CCB-based strategies in chronic kidney disease. The outcome is
  heart-failure hospitalization and death is a competing event.

## Quick start: estimate the ITT weight from the data

The complete tutorials below are recommended for applied use. The compact
example here illustrates the main logic.

```r
library(TTE)

data(SGLT2_baseline)
data(SGLT2)

baseline <- SGLT2_baseline
long <- SGLT2

# 1. Baseline treatment weight
w_a <- est_wt(
  A ~ age + I(age^2) + female + bmi +
    hba1c + I(hba1c^2) + egfr + proteinuria +
    prior_heart_failure + prior_stroke +
    recent_hospitalization + trial_period,
  data = baseline,
  type = "treatment",
  stabilize = TRUE,
  truncate = c(0, 1)
)

# 2. Copy the baseline weight to each person-month
key_baseline <- paste(baseline$id, baseline$trial, sep = ":")
key_long <- paste(long$id, long$trial, sep = ":")

long$w_a <- weights(w_a, which = "untruncated")[
  match(key_long, key_baseline)
]

# 3. Loss-to-follow-up weight
w_c <- est_wt(
  stay_ltfu ~ A + splines::ns(time, df = 3) +
    age + female + bmi + hba1c + egfr + proteinuria +
    prior_heart_failure + prior_stroke +
    recent_hospitalization + trial_period,
  numerator = stay_ltfu ~
    A + splines::ns(time, df = 3) + trial_period,
  data = long,
  type = "censoring",
  id = id,
  trial = trial,
  time = time,
  cumulative = TRUE,
  lag = 1,
  stabilize = TRUE,
  truncate = c(0, 1)
)

# 4. Final ITT weight
w_itt <- combine_wt(
  long$w_a,
  w_c,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

long$w_itt_est <- weights(w_itt)

# 5. Weighted outcome model
fit <- discsurvreg(
  Y_death ~ A + splines::ns(time, df = 3) + trial_period,
  data = long,
  id = id,
  weights = w_itt_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

confint(fit, parm = "A", eform = TRUE)

# 6. Standardized absolute risks
std <- std_tte(
  fit,
  data = baseline,
  treatment = A,
  time = time,
  times = 0:59,
  values = c(0, 1),
  labels = c("DPP-4i", "SGLT2i")
)

summary(std, horizon = 60)
plot(std, measure = "risk")
```

## Hands-on tutorials

Each tutorial estimates the analysis weights from the observed synthetic
data rather than relying on the reference weights bundled with the datasets.

### SGLT2i versus DPP-4i: all-cause death

- [Japanese tutorial](https://github.com/nomahi/TTE/blob/main/tutorials/ja/SGLT2_TTE_tutorial_ja.R)
- [English tutorial](https://github.com/nomahi/TTE/blob/main/tutorials/en/SGLT2_TTE_tutorial_en.R)

### ARB versus CCB: heart-failure hospitalization with competing death

- [Japanese tutorial](https://github.com/nomahi/TTE/blob/main/tutorials/ja/ARB_TTE_tutorial_ja.R)
- [English tutorial](https://github.com/nomahi/TTE/blob/main/tutorials/en/ARB_TTE_tutorial_en.R)

The SGLT2i example introduces IPTW, IPCW, standardized risks, weighted
Kaplan-Meier estimation, and ITT/PP estimands. The ARB example extends the
workflow to sequentially nested trials and competing-risk estimation using
cause-specific models and weighted Aalen-Johansen curves.

## Time scale and curve display

Person-period indices are stored from zero internally. Public summaries use
ordinary elapsed follow-up time, so `time = 60` and `horizon = 60` both refer
to 60 months.

Weighted Kaplan-Meier and Aalen-Johansen curves are displayed as step
functions.

## Scope and assumptions

`TTE` provides computational tools, not automatic identification of causal
effects. Applied analyses still require defensible definitions of eligibility,
treatment strategies, time zero, follow-up, outcomes, censoring, confounders,
and estimands. Users should assess consistency, exchangeability, positivity,
model specification, missing data, and measurement timing in the context of
their own study.

For publication-quality uncertainty intervals for standardized risks or
weighted curves, a full individual-cluster bootstrap that re-estimates all
weight and outcome models within each bootstrap sample is generally preferable
to treating estimated weights as fixed.

## Motivating studies

- Noma H, Goto A, Sugimoto T, Sunada H, Oda F, Maeda M, Fukuda H.
  *Real-world effectiveness of SGLT2 inhibitors in adults aged 75 years or
  older: a target trial emulation.* Age and Ageing. 2026; in press.
- Noma H, Kurita N, Fukuma S, Fujisawa T, Oda F, Maeda M, Fukuda H.
  *Heart failure and renal outcomes with angiotensin receptor blockers
  compared with calcium channel blockers in patients with chronic kidney
  disease: a target trial emulation.* Heart. 2026.
  <https://doi.org/10.1136/heartjnl-2026-328193>

## Citation

After installation, run:

```r
citation("TTE")
```

## License

GPL-3.

## Author

Hisashi Noma  
Institute of Statistical Mathematics  
ORCID: <https://orcid.org/0000-0002-2520-9949>
