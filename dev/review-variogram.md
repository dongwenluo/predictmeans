# Residual variogram for field trials (review note, 2026-11-25)

For the maintainer's decisions before any code.

## Why

Field trials laid out in rows and columns often have spatial trend: neighbouring
plots more alike than the blocks allow. A residual-versus-fitted plot does not
show it, and the means and their SEDs are then wrong (too optimistic for some
comparisons, too pessimistic for others). The sample variogram of the residuals
is the standard check (Gilmour, Cullis and Verbyla 1997; Stefanova, Smith and
Cullis 2009). predictmeans has no spatial diagnostic now.

What exists elsewhere: `nlme::Variogram()` (lme/gls fits only; one isotropic
curve by distance, not the row x column surface); gstat and similar
(geostatistics on coordinate objects, no model fits); ASReml and Genstat (the
row x column sample variogram for their own fits). Nothing reads lmer, glmer,
glmmTMB or semireg fits.

## What it would compute

`residvariogram(model, row = , col = , coords = , by = , ...)`: a classed
result (table of lags, semivariances and numbers of pairs) with `autoplot()`.

- **Residuals**: the ones `residplot()` uses (one place): conditional residuals
  y - X b - Z u for mixed models (the blocks removed, as ASReml), randomized
  quantile residuals for Poisson, binomial and negative binomial, the fit's
  residuals for semireg and gam.
- **Row x column layout** (`row`, `col`, integer positions): the sample
  variogram gamma(dr, dc) = sum (e_i - e_j)^2 / (2 N(dr, dc)) over pairs of
  plots dr rows and dc columns apart (absolute lags), computed by shifting the
  residual matrix of the grid (no list of n^2 pairs; fast for large trials).
  Shown as a surface (tile plot) with the two edge curves gamma(dr, 0) and
  gamma(0, dc), the directions an AR1 x AR1 model would describe.
- **Coordinates** (`coords = c("x", "y")`, or `c("lat", "long")` converted to km
  by the haversine formula): the isotropic variogram by distance class
  (Matheron's estimator; `robust = TRUE` for Cressie and Hawkins'), up to half
  the largest distance, 15 classes by default; `directions` (degrees from east)
  and `segment` (the angle each covers) for directional variograms, to see
  anisotropy, as Genstat's `FVARIOGRAM` (DIRECTIONS, SEGMENTS). Genstat's
  guide (8.3.1, Webster and Oliver) advises some 300 data for directions.
- **Several trials** (`by = "trial"`): one variogram per trial, never pooled
  (lags across trials have no meaning); facets in the plot.
- **Reference band**: the variogram expected without spatial correlation, from
  permutations of the residuals within each trial (99 by default, seeded).
  A global envelope (Myllymäki et al. 2017), so that looking at many lags does
  not give false alarms; the plot shades it.

## Advice (only when the user gives a layout)

In `residvariogram()` and, with `row`/`col` or `coords` given, in
`pm_check()`: when the variogram leaves the band at short lags, one line saying
in which direction and what to fit: an AR1 x AR1 covariance (glmmTMB `ar1()`
per direction; nlme has no separable two-dimensional AR1, `corExp`/`corGaus` are
isotropic), or a smooth surface `s(row, col)` / `t2(row, col)` in `semireg()`
(the SpATS approach, Rodríguez-Álvarez et al. 2018). A note when the layout is
not a full grid (missing plots: fewer pairs at some lags) or when a cell holds
several plots.

Not done: fitting spatial models, or the theoretical variogram of a fitted
spatial model (possible later, for glmmTMB `ar1()`/`exp()` and nlme
`corSpatial`).

## Decisions

**V1. Name and place.** Proposal: a new function `residvariogram()` (a
diagnostic like `residplot()`), and `pm_check(row = , col = , coords = )`
calling it. Alternative: `residplot(type = "variogram")`.

**V2. Residuals.** Proposal: conditional residuals (as ASReml and
`residplot()`), standardised by their SD so that trials with different error
variances share a scale.

**V3. Reference band.** Proposal: a global envelope from permutations within
each trial (99), with the advice based on it; pointwise bands only as an
option.

**V4. Coordinates.** Proposal: grid lags for row/column data, distance classes
for coordinates; latitude/longitude in km.

**V5. References.** Genstat Guide, section 8.3.1 (extract from the
maintainer, 2026-11-25): `FVARIOGRAM` gives Matheron's estimator (Eq. 8.3.12)
and Cressie-Hawkins (8.3.13; also Dowd, Genton); its Example 8.3.1 prints the
semivariances and pair counts of log10 potassium at Brooms Barn (434 points,
4 directions, lags 1-12; Webster and Oliver 1990) - an exact reference for the
estimator if the data file Broomesb.dat (shipped with Genstat) can be shared.
That chapter is geostatistics (model fitting `MVARIOGRAM`, `KRIGE`), not the
variogram of REML residuals of a field trial, which is the main use here.
Other references: `nlme::Variogram()` for the isotropic estimator on lme/gls
fits (to 1e-8); Gilmour et al. (1997) wheat trial (agridat
`gilmour.serpentine`, 330 plots, 22 rows x 15 columns, 107 genotypes, 3
replicates) for the row x column variogram and its AR1 x AR1 interpretation; a
Genstat or ASReml run by the maintainer on the same trial, if possible.

**V6. Dependencies.** None new (ggplot2, stats; the permutations through
`pm_perm_map()`).

## References

Cressie, N. (1993). *Statistics for Spatial Data*, revised edition. Wiley.

Gilmour, A. R., Cullis, B. R. and Verbyla, A. P. (1997). Accounting for natural
and extraneous variation in the analysis of field experiments. *Journal of
Agricultural, Biological, and Environmental Statistics*, 2, 269-293.

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017).
Global envelope tests for spatial processes. *Journal of the Royal Statistical
Society B*, 79, 381-404.

Rodríguez-Álvarez, M. X., Boer, M. P., van Eeuwijk, F. A. and Eilers, P. H. C.
(2018). Correcting for spatial heterogeneity in plant breeding experiments with
P-splines. *Spatial Statistics*, 23, 52-71.

Stefanova, K. T., Smith, A. B. and Cullis, B. R. (2009). Enhanced diagnostics
for the spatial analysis of field trials. *Journal of Agricultural, Biological,
and Environmental Statistics*, 14, 392-410.
