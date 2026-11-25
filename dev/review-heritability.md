# Heritability in predictmeans (review note, 2026-11-25)

For the maintainer's decisions before any code. Sources: Falconer and Mackay
(1996); Cullis, Smith and Coombes (2006); Oakey et al. (2006); Piepho and Möhring
(2007); Schmidt, Hartung, Bennewitz and Piepho (2019, the review of all
line-mean forms); Mrode (2014); the maintainer's summary of 2026-11-25.

## What to compute

Line-mean (entry-mean) heritability: the share of the variance of genotype
means that is genetic. It sets the expected response to selection on means and
depends on the design (replicates, environments, blocks, spatial errors).

| method | formula | needs |
|---|---|---|
| `"standard"` | s2g / (s2g + s2ge/e + s2e/(r e)) | balanced data; r, e from the design |
| `"cullis"` | 1 - vBLUP / (2 s2g) | genotype random; vBLUP = mean variance of a difference of BLUPs (PEV of the differences) |
| `"piepho"` | s2g / (s2g + vBLUE/2) | s2g from the random-genotype fit; vBLUE = mean variance of a difference of BLUEs, from the same model with genotype fixed (one refit) |
| `"oakey"` | 1 - tr(G^-1 C_gg)/m (non-zero eigenvalues) | genotype random; the reliability form; = Cullis with independent genotypes |
| `"family"` | (s2f / rho) / (s2f + s2plot + s2e) on a single-tree basis; family-mean form with r, e as above | half-sib / full-sib families: rho = 1/4 (half-sibs), 1/2 (full-sibs), the user's value for open-pollinated progeny (often 1/3) |
| `"repeatability"` | (s2g + s2pe) / (s2g + s2pe + s2e) | repeated records on the same subject |

Correction to the summary received: Piepho and Möhring's ad hoc measure is
s2g / (s2g + vBLUE/2), not 1 - vBLUE/(2 s2g) (that is the Cullis form with BLUEs,
and can be negative). All forms agree for balanced data with independent
genotypes (Schmidt et al. 2019); they differ with imbalance, incomplete blocks,
spatial errors and kinship.

Not possible: narrow-sense h2 or reliability from a pedigree or genomic
relationship matrix (animal models, GBLUP). lme4 has no known relationship
matrix; `heritability()` says so and names software that fits it.

## How it fits predictmeans

New function `heritability(model, genotype, method = , ...)`; no options added
to `varcomp()`. The pieces exist, each computed in one place:
- s2g, s2ge, s2e and their covariance: `varcomp()` (observed information).
- vBLUP: the prediction-error covariance of the BLUPs, `pm_ranef_cov()`
  (`se_ranef(type = "prediction")`, checked against SAS Std Err Pred).
- vBLUE: SEDs of genotype means from `predictmeans()` on the fixed-genotype
  refit (`pm_refit()`).
- r and e for `"standard"`: from the data, the harmonic mean when unequal, with
  a caution that the standard form then is approximate.

Output: one table in the package layout (`estimate`, `std.error`, `conf.low`,
`conf.high`, method, components used), and advice lines only when the fit calls
for them (unbalanced data with `"standard"`; s2g at zero; methods that disagree
by more than 0.1).

## Decisions

**H1. Default method.** Proposal: `"cullis"` for a fit with random genotypes
(any design, the most used; = standard when balanced); `"standard"` only on
request. `method = "all"` prints the forms side by side, as Schmidt et al.
recommend for reporting.

**H2. SE and interval.** Proposal: delta method from `varcomp()`'s covariance
of the components for `"standard"`, `"family"` and `"repeatability"`; for
`"cullis"`/`"piepho"`/`"oakey"`, where vBLUP and vBLUE also depend on the
variance parameters, a parametric bootstrap (`bootmodels()` engine, seeds per
draw) with the delta method as the quick option. Intervals on the logit scale
of h2, so they stay in [0, 1].

**H3. Models.** lmer first (one or several environments, G x E, blocks);
then lme/gls (spatial residuals), glmmTMB, semireg; aov with Error() for
balanced `"standard"`/`"family"`. GLMMs (binary traits): latent-scale h2 with
the link variance (pi^2/3 for logit), as the R2 code does, in a later step.

**H4. Several environments.** G x E as a random genotype:environment term
(`"standard"` uses e; the mixed-model forms read it from the fit). Per
environment heritability with `by = "env"` later.

**H5. References.** Genstat (its heritability procedure, `VHERITABILITY` if the version has it; run by the maintainer); the R
package inti (`H2cal()`: standard, Cullis, Piepho; CRAN, installable here);
the family examples of the tree-improvement book (eda4treeR data). Tolerances:
mixed-model 2%; exact formulas at given components 1e-6.

**H6. Name and place.** `heritability()` (also `repeatability`, `family` as
methods); a section in the tree-improvement article with the book's examples.
