# TTE

**Design and Analysis Tools for Target Trial Emulation**

`TTE` is an R package providing design and analysis tools for target trial
emulation with longitudinal observational data. Its functions cover
person-period data checks, sequentially nested trial construction, inverse
probability weighting, diagnostics, weighted pooled discrete-time outcome
models, standardization, and weighted Kaplan-Meier and Aalen-Johansen
estimation.

The repository accompanies **TTE version 1.1.1**, the first CRAN release.

- Package website: <https://github.com/nomahi/TTE/>
- CRAN: <https://CRAN.R-project.org/package=TTE>
- PDF manual: <https://cran.r-project.org/web/packages/TTE/TTE.pdf>
- Vignette: <[https://github.com/nomahi/TTE/vignette.pdf](https://github.com/nomahi/TTE/blob/main/vignette.pdf)>

## Installation

Install the CRAN release with:

```r
install.packages("TTE")
```

## Main functions

| Task | Function |
|---|---|
| Check a target-trial person-period dataset | `check_tte()` |
| Expand longitudinal data into sequentially nested trials | `seqdesign_tte()` |
| Estimate treatment, censoring, or adherence weights | `est_wt()` |
| Combine inverse probability weight components | `combine_wt()` |
| Assess covariate balance before and after weighting | `balance_wt()` |
| Diagnose weights, effective sample size, and weighted risk sets | `diagnose_wt()` |
| Fit a weighted pooled discrete-time outcome model | `discsurvreg()` |
| Standardize fitted models to marginal risks and treatment contrasts | `std_tte()` |
| Estimate weighted Kaplan-Meier or Aalen-Johansen curves | `curve_tte()` |
| Run an individual-cluster bootstrap | `boot_tte()` |

## Synthetic example datasets

The package includes two pairs of fully synthetic datasets. They contain no
actual participant records from the motivating studies.

### `SGLT2_baseline` and `SGLT2`

An active-comparator new-user target trial emulation comparing initiation of
an SGLT2 inhibitor with initiation of a DPP-4 inhibitor among older adults
with type 2 diabetes. The primary outcome is all-cause death.

### `ARB_baseline` and `ARB`

Sequentially nested new-user trials comparing ARB and CCB treatment strategies
among people with chronic kidney disease. The primary outcome is heart-failure
hospitalization, and death is a competing event.

## Quick start

The full tutorials estimate all analysis weights from the observed synthetic
data. This compact example shows the principal ITT workflow.

```r
library("TTE")

data(SGLT2_baseline)
data(SGLT2)

baseline <- SGLT2_baseline
long <- SGLT2

# Baseline stabilized treatment weight
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

# Map the baseline treatment weight to all person-month records
key_baseline <- paste(baseline$id, baseline$trial, sep = ":")
key_long <- paste(long$id, long$trial, sep = ":")

long$w_a <- weights(w_a, which = "untruncated")[
  match(key_long, key_baseline)
]

# Stabilized loss-to-follow-up weight
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

# Final ITT weight, truncated after combining components
w_itt <- combine_wt(
  long$w_a,
  w_c,
  truncate = c(0.01, 0.99),
  normalize = "none"
)

long$w_itt_est <- weights(w_itt)

# Weighted pooled discrete-time outcome model
fit <- discsurvreg(
  Y_death ~ A + splines::ns(time, df = 3) + trial_period,
  data = long,
  id = id,
  weights = w_itt_est,
  family = quasibinomial(link = "cloglog"),
  var_method = "standard"
)

confint(fit, parm = "A", eform = TRUE)

# Standardized absolute risks
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

The public scripts provide complete ITT and PP workflows. The main analyses
estimate treatment and censoring weights from the observed synthetic data
rather than relying on the reference weights bundled with the datasets.

### SGLT2i versus DPP-4i: all-cause death

- [Tutorial](SGLT2_TTE_tutorial_en.R)

This example covers baseline IPTW, loss-to-follow-up IPCW, weighted pooled
discrete-time modeling, standardized absolute risks, weighted Kaplan-Meier
curves, and per-protocol artificial-censoring analyses.

### ARB versus CCB: heart-failure hospitalization with competing death

- [Tutorial](ARB_TTE_tutorial_en.R)

This example extends the workflow to sequentially nested trials,
cause-specific outcome models, standardized cumulative incidence functions,
weighted Aalen-Johansen curves, and per-protocol estimation.

## Follow-up time and curve display

Person-period interval indices are stored from zero internally. Public
summaries use elapsed follow-up time. Thus, with `times = 0:59`:

```r
summary(std, horizon = 60)
summary(km, time = 60)
```

both report results after 60 follow-up intervals.

Weighted Kaplan-Meier and Aalen-Johansen curves are displayed as
right-continuous step functions.

## Interpretation and assumptions

`TTE` supplies computational tools; it does not automatically identify a
causal effect. Applied analyses still require defensible definitions of:

- eligibility criteria;
- treatment strategies;
- assignment procedures;
- time zero;
- induction or grace periods;
- outcomes and competing events;
- follow-up and censoring;
- confounders;
- estimands; and
- analysis populations.

Users must assess consistency, exchangeability, positivity, model
specification, missing data, measurement timing, and the clinical meaning of
the intervention strategies.

For publication-quality uncertainty intervals for standardized risks or
weighted curves, an individual-cluster bootstrap that re-estimates all
treatment, censoring, adherence, and outcome models within each bootstrap
sample is generally preferable to treating estimated weights as fixed.

## Motivating studies

- Noma H, Goto A, Sugimoto T, Sunada H, Oda F, Maeda M, Fukuda H.
  *Real-world effectiveness of SGLT2 inhibitors in adults aged 75 years or
  older: a target trial emulation.* Age and Ageing. 2026; in press.
- Noma H, Kurita N, Fukuma S, Fujisawa T, Oda F, Maeda M, Fukuda H.
  *Heart failure and renal outcomes with angiotensin receptor blockers
  compared with calcium channel blockers in patients with chronic kidney
  disease: a target trial emulation.* Heart. 2026.
  <https://doi.org/10.1136/heartjnl-2026-328193>

## Author

Hisashi Noma  
The Institute of Statistical Mathematics  
ORCID: <https://orcid.org/0000-0002-2520-9949>
