# copyfair

An R package for item- and ability-conditioned answer-copying
detection in multiple-choice examinations.

Classical answer-copying indices (K-index, S2, omega, GBT) compare an
examinee pair's error similarity against a single global or
ability-marginal baseline. This can conflate genuine collusion with
population-wide "attractive distractor" effects, where unrelated
examinees independently select the same wrong option because it
reflects a common misconception rather than copying.

`copyfair` extends the ability-conditioned diagnostic logic of
[condfair](https://github.com/tzningsih/condfair) (Ningsih, Aman, &
Nasrulloh, 2026, *Applied Psychological Measurement*) to this problem:
it models the expected number of matching incorrect responses jointly
as a function of examinee ability **and** item-specific distractor
popularity, tests localized disparity with a wild-bootstrap omnibus
procedure, and returns Benjamini-Hochberg-corrected pairwise flags.

## Installation

```r
# install.packages("remotes")
remotes::install_github("tzningsih/copyfair")
```

## Quick start

```r
library(copyfair)

# simulate a 40-item exam with an ability-concentrated copying signal
exam <- build_copyfair_demo(n_examinees = 300, n_items = 40,
                             n_copy_pairs = 10, ability_band = c(-0.5, 0.5))

key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
result <- copyfair_analyze(exam$responses, key)

result$omnibus_test$p_value      # test for localized disparity
head(result$pairs[order(-result$pairs$z), ])  # most suspicious pairs
```

## Core functions

| Function | Purpose |
|---|---|
| `copyfair_pairs()` | Computes observed, expected, and standardized-residual match statistics for every examinee pair sharing a minimum number of jointly administered items |
| `copyfair_flag()` | One-sided p-values with Benjamini-Hochberg FDR correction |
| `copyfair_omnibus_test()` | Wild-bootstrap test for localized (ability-concentrated) disparity |
| `copyfair_analyze()` | Wrapper chaining the three functions above |
| `simulate_exam()`, `inject_copying()`, `build_copyfair_demo()` | Simulate exam data with a known, injected collusion signal for validation |
| `estimate_ability()` | Lightweight 2PL MLE ability estimator (no external IRT dependency) |

## Reproducing the paper's analyses

The `analysis/` directory contains scripts reproducing the simulation
study and real-data applications described in the accompanying
manuscript:

- `analysis/simulation_study.R` — self-contained; generates its own
  data.
- `analysis/pirls_reproduce.R` — requires the PIRLS 2011 International
  Database, obtained separately from the
  [IEA TIMSS & PIRLS International Study
  Center](https://timssandpirls.bc.edu/pirls2011/international-database.html)
  under its non-commercial research-use terms (not redistributed
  here).

## Citation

If you use `copyfair`, please cite:

> Ningsih, T. Z., Aman, A., & Nasrulloh, A. (2026). copyfair: An R
> package for item- and ability-conditioned answer-copying detection
> in multiple-choice examinations. *Applied Psychological
> Measurement*. [manuscript in preparation]

and the companion package:

> Ningsih, T. Z., Aman, A., & Nasrulloh, A. (2026). condfair: An R
> package for ability-conditioned fairness and explanation diagnostics
> in automated scoring. *Applied Psychological Measurement*.
> https://doi.org/10.1177/01466216261484166

## License

MIT © Tri Zahra Ningsih, Aman Aman, Ahmad Nasrulloh
