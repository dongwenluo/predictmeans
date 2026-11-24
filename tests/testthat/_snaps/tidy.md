# with a permlist the limits of the means are the model's t limits at conf_level (not +- 2 SE)

    Code
      print(summary(predictmeans(fm, "nitro", permlist = pl)))
    Output
      Predicted means for 'nitro' (lm model; df: permutation)
      Confidence limits: 95% (the model's t limits of each mean; the permutation tests are of the comparisons)
      
       nitro estimate std.error df conf.low conf.high group
           0    79.39     3.767 NA    71.86     86.92     A
         0.2    98.89     3.767 NA    91.36    106.42     A
         0.4   114.22     3.767 NA   106.69    121.75     A
         0.6   123.39     3.767 NA   115.86    130.92     A
      
      6 pairwise comparisons.
      SED: min 5.327, mean 5.327, max 5.327
      LSD (5%): min 7.066, mean 8.75, max 10.97 (they differ: use each pair's p-value or LSD)
      p-values: permutation (19 permutations), adjusted by 'none'.
      Means sharing a letter are not significantly different at the 5% level.
      Comparison intervals: mean +/- 1.161 SE (75.5% intervals); means differ at the 5% level when these do not overlap (the SEDs differ: 1 of 6 pairs are judged differently by their own test). The 95% confidence limits are not for comparing means.
      
      Advice:
      i 6 pairs tested without adjustment: some will be significant by chance. If
        they are a family, use adj = "holm" (any design) or "tukey"/"mvt" (all
        pairs).

