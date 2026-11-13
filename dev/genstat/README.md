# Genstat reference runs

`pm-reference.gen` asks Genstat's REML for:

1. **Oats split plot with a whole plot missing** (Block I, Golden Rain).
   pbkrtest's Kenward–Roger gives 9.25 df for the Variety differences, SAS's
   `DDFM=KR` 8.35 (it uses the observed information; `df_method = "Kenward-Roger"`).
   Genstat fits REML by average information (`METHOD=AI`), so its df may
   match either, or neither.
2. **AR(1) repeated measures with 5 of 48 values missing** — the same data
   as `dev/sas/pm-reference-ar1.sas` (`dev/sas/ar1.txt`). For AR(1), SAS's
   default `DDFM=KR` adds second-derivative terms; `KR(FIRSTORDER)` should
   match predictmeans.

The script was written without Genstat at hand; the first run may need small
syntax fixes. Reference values come only from the output Genstat writes
(the output window saved as a file, or the `.out` file), not from values
typed into chat.

## The first run (Genstat 23.1, 2026-11-13; output window pasted into chat)

1. **Oats, whole plot missing.** Variance components as SAS (Block 179.2,
   Block.Variety 131.3, residual 169.4). The Variety F test has **8.7** df,
   between SAS's observed-information Kenward–Roger (8.26) and pbkrtest's
   expected information (9.16): predictmeans' Kenward–Roger computed with
   Genstat's average-information matrix gives 8.75 (see the second run). Genstat's SEDs are the **unadjusted** ones
   (Variety 8.101, 7.609, average 7.937; the Kenward–Roger adjusted are
   8.18, 7.609), and its F statistic is Wald/df (1.28). Nitro and
   Variety.Nitro: 42.0 df.
2. **AR(1)** did not run as intended: Genstat added a residual term next to
   the AR(1) term on `Subj.Time` ("aliased random term", no convergence, no
   F-test df), because the 5 units with a missing response were dropped and
   `Subj.Time` was no longer the unit term. The script now asks REML to keep
   them (`MVINCLUDE=yvariate`), so that `Subj.Time` is the residual term with
   the AR(1) structure.

## The second run (AR(1) part, 2026-11-13; output window pasted into chat)

With `MVINCLUDE=yvariate` the AR(1) model fits:

- Variance parameters as SAS and `nlme::gls()`: residual variance 1.194,
  AR(1) correlation 0.6157.
- SEDs are the **unadjusted** ones, equal to predictmeans with
  `vcov_adjust = FALSE`: Group 0.4970; Time max 0.4139, min 0.2940;
  Time.Group max 0.6692, min 0.3912.
- Denominator df: Time.Group F 0.77 (Wald 2.32, as predictmeans) on **25.8**
  df. Kenward–Roger with the observed information (SAS, predictmeans) gives
  26.62, with the expected information (pbkrtest) 25.04. Kenward–Roger with
  Genstat's average-information matrix, A_ij = ½ y'P G_i P G_j P y (the
  matrix its REML iterations use), gives **25.76**, and **8.75** for the
  oats Variety test (Genstat 8.7). So Genstat's df are Kenward–Roger with
  the average information, while its SEDs and F statistics are not adjusted.
  `df_method = "KR-average"` gives these df (with `vcov_adjust = FALSE`
  also Genstat's SEDs).
- Time (3.71 on 25.6 df) and Group (5.62 on 11.7) are sequential (added in
  turn) tests, not SAS's Type 3, so they are not compared.

These are notes, not test references: a reference needs the output file.

## The second program (`pm-reference-2.gen`, to be run)

The oats split plot with four patterns of missing data, each testing a
different rule of predictmeans. What predictmeans gives (lmer REML), to
compare with Genstat's output once it is run (values from the output file
only):

| pattern | n | what it tests | predictmeans |
|---|---|---|---|
| A: whole plot I x GoldenRain + sub-plot II x Victory x 0.6 (= SAS run 3) | 67 | which Kenward-Roger information Genstat's df use | Variety F 1.21 on 8.30 (observed, = SAS), 9.14 (expected), 8.75 (average); unadjusted SEDs 7.813, 7.846, 7.364 |
| B: empty cell Victory x 0.6 (all blocks) | 66 | means over an empty cell (Genstat `PRESENT`) | Nitro means 79.39, 98.89, 114.22, 125.83; Variety 104.50, 109.79, 90.67; 0.6 - 0 compared on GoldenRain and Marvellous: 42.5 (SED 5.474) |
| C: sub-plot II x Victory x 0.6 | 71 | `aov(... Error())` by strata (predictmeans combines the strata) vs Genstat ANOVA (missing value estimated) and REML | REML: Variety F 1.38 on 9.91 / 9.99 / 9.95; aov strata: GoldenRain x 0 mean 79.94 |
| D: five scattered sub-plots | 67 | REML df with an unbalanced sub-plot stratum | Variety F 1.24 on 9.96 / 9.95 / 9.96; Nitro F 36.86 on 40.55 / 40.52 / 40.53 |

## Runs of 2026-11-14 (output files): reference tests

The maintainer ran both programs in Genstat 23.1 and saved the Output windows
as RTF. `rtf2txt.py` turned them into `genstat-results-2026-11-14-run1.txt`
(`pm-reference.gen`: oats with a whole plot missing, AR(1)) and `-run2.txt`
(`pm-reference-2.gen`: patterns A-D); the line with the user's working
directory is dropped. `parse-genstat.py` writes `genstat-results.csv` (383
printed values), tested in `tests/testthat/test-genstat-reference.R`.

What agrees (within one unit of Genstat's last printed digit), for all six
REML analyses:

- variance components (and the AR(1) correlation);
- predicted means (random terms excluded) and their SEDs (average, maximum,
  minimum, and within the same level of each factor): those of
  `vcov_adjust = FALSE`;
- the Wald tests that add the terms in turn (each computed from the full fit:
  the rows of the Cholesky factor of X'V^-1X for the term): Wald statistics,
  n.d.f.;
- their denominator df: Kenward-Roger with the average information
  (`df_method = "KR-average"`) in 17 of 18 tests;
- the F statistic: Kenward-Roger's scaling factor (average information) times
  Wald / n.d.f.; the factor is about 1 except for AR(1) Time (0.9959:
  3.7264 -> 3.71 as Genstat prints);
- pattern B (empty cell): Genstat's predictions averaging over the
  combinations present (`VPREDICT [PRESENTCOMBINATIONS=Variety, Nitro]`) equal
  predictmeans' means over empty cells (Nitro 79.39, 98.89, 114.22, 125.83;
  Variety 104.5, 109.8, 90.7), with the same SEs for Nitro.

What differs, and why:

- pattern B, Variety: d.d.f. 10.4 against Kenward-Roger 10.10 (average),
  10.19 (observed), 10.01 (expected); the Wald statistic agrees. It is the only
  rank-deficient model (the empty cell aliases one interaction coefficient),
  so Genstat's df probably treat the aliased parameter differently.
- pattern B, SEDs of the PRESENT predictions: Genstat gives those of the
  differences of the averages (0.6 against the others 5.214); predictmeans
  compares on the cells both means have (5.474), as decided after
  `dev/sim-empty-cells.R`.
- pattern B, `VPREDICT` for Variety: Genstat included the random term
  Block.Variety in the prediction (SEs 6.457, SEDs 3.871-4.181), a different
  prediction from the population means; not compared.
- pattern C, ANOVA by strata: Genstat estimates the missing value (113.0) and
  analyses the balanced design (Victory x 0.6 mean 120.7; SEDs Variety 6.88 on
  10 df, Nitro 4.46 on 44, Variety.Nitro 9.59 on 31.42, "not adjusted for
  missing values"). predictmeans' `aov(... Error())` combines the strata
  instead (Victory x 0.6 121.8; Variety SED about 6.63 on 9 df); Genstat's
  REML gives 121.5. A decision for the maintainer (see the plan).
