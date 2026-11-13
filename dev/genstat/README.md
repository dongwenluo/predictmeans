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
