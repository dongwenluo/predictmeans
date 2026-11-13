# Genstat reference runs

`pm-reference.gen` asks Genstat's REML for:

1. **Oats split plot with a whole plot missing** (Block I, Golden Rain).
   pbkrtest's Kenward–Roger gives 9.25 df for the Variety differences, SAS's
   `DDFM=KR` 8.35 (it uses the observed information; `df_method = "KR-SAS"`).
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
   the *average* information (the half-sum of the two, which Genstat's AI
   algorithm uses) gives 8.75. Genstat's SEDs are the **unadjusted** ones
   (Variety 8.101, 7.609, average 7.937; the Kenward–Roger adjusted are
   8.18, 7.609), and its F statistic is Wald/df (1.28). Nitro and
   Variety.Nitro: 42.0 df.
2. **AR(1)** did not run as intended: Genstat added a residual term next to
   the AR(1) term on `Subj.Time` ("aliased random term", no convergence, no
   F-test df), because the 5 units with a missing response were dropped and
   `Subj.Time` was no longer the unit term. The script now asks REML to keep
   them (`MVINCLUDE=yvariate`), so that `Subj.Time` is the residual term with
   the AR(1) structure; to be rerun.

These are notes, not test references: a reference needs the output file.
