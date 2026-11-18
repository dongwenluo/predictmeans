# Functions removed in predictmeans 2.0

These functions were removed in 2.0. Each stops with what to use
instead;
[`vignette("migrating")`](https://dongwenluo.github.io/predictmeans/articles/migrating.md)
shows the 1.1.1 code beside the 2.0 code.

|  |  |
|----|----|
| **Removed** | **Use** |
| `smZ()`, `smoothZ =` | smooths in the formula of [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) or [`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md), e.g. `smoothZ = list(sm = smZ(x))` is `s(x, bs = "os")` |
| `permglmer()`, `permglmmTMB()` | [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md) (they were the same function) |
| `ci_mcp()` | [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)`(..., letterCI = TRUE)`: letters from comparison intervals, whose overlap agrees with the pairwise tests (overlap of confidence intervals does not) |
| `R2_glmm()` | [`R2_model`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)`(model)`: the same values, as numbers |
| `permanova.lmer()` | [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)`(model, type = )`: Freedman-Lane permutations, so the p-values differ from the refits of `permanova.lmer()` |

## Arguments

- ...:

  Not used.
