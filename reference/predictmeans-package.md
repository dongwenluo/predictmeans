# Predicted Means for Linear and Semiparametric Models

This package provides functions to diagnose and make inferences from
linear, generalized linear, mixed and semiparametric models, such as
those obtained from 'aov', 'lm', 'glm', 'gls', 'lme', 'lmer', 'glmer',
'glmmTMB', 'gam' and 'semireg'. Inferences include predicted means and
standard errors, contrasts, multiple comparisons, permutation tests,
adjusted R-square and graphs (`autoplot()`). Vignettes
(`browseVignettes("predictmeans")`) show worked examples. For AI
assistants and scripted use,
`system.file("llms-context.md", package = "predictmeans")` gives a
compact description of the result tables, defaults and common recipes.

## Details

|          |              |
|----------|--------------|
| Package: | predictmeans |
| Type:    | Package      |
| License: | GPL (\>= 2)  |

## Workflow

![From a fitted model to the report: diagnostics, tests, means and
contrasts, variance components, then tables and
plots](figures/workflow.png)

## What predictmeans aims for

- Defaults chosen on evidence (simulation, literature): Kenward-Roger df
  and adjusted SEs for REML `lmer` fits; permutation p-values \\(b +
  1)/(B + 1)\\, never 0; Genstat-style stratum SEs for `aov` models with
  `Error()` strata.

- Checked against outside references: SAS PROC MIXED and Genstat runs
  (`dev/sas`, `dev/genstat` in the package sources), pbkrtest, nlme and
  emmeans, in tests with stated tolerances.

- Works across `lm`, `glm`, `aov`, `lmer`, `glmer`, `lme`, `gls`,
  `glmmTMB`, `gam` and `semireg` fits; when a result cannot be computed
  correctly it stops with a clear message.

- Advice triggered by this fit and data, read with
  [`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md).

- Results are the same for any `ncore`.

## References

Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
*Prediction in linear mixed models*, Australian and New Zealand Journal
of Statistics, 46(3), 325-347.

## See also

[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`autoplot.predictmeans`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md),
[`tidy.predictmeans`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md),
[`pm_check`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md)

## Author

Dongwen Luo, Siva Ganesh, John Koolaard and James Curran

Maintainer: Dongwen Luo \<dongwen.luo@agresearch.co.nz\>
