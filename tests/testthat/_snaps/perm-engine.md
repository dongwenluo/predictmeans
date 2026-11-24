# exact = TRUE counts the observed statistic unless perms enumerate every permutation

    Code
      r <- permmodels(lm(Reaction ~ Days, sleepstudy), nperm = 19, seed = 1, exact = TRUE)
    Condition
      Warning:
      'exact = TRUE' needs 'perms' holding all the permutations the design allows; with random permutations the observed statistic is counted as one of them (p = (b + 1)/(B + 1)).

---

    Code
      e <- permmodels(lm(y ~ x, d4), perms = cbind(P, P[, 1]), exact = TRUE, method = "raw")
    Condition
      Warning:
      'exact = TRUE' ignored: the columns of 'perms' are not distinct, so they are not a complete enumeration.

# permmodels() prints tiny model p-values as <1e-04 and notes the smallest attainable permutation p-value

    Code
      print(r)
    Output
      
      Coefficients of (fixed) effects:
                  Estimate Std. Error t value Pr(>|t|) Perm_p_value
      (Intercept) 251.4051     6.6102 38.0332   <1e-04           NA
      Days         10.4673     1.2382  8.4537   <1e-04         0.05
      
      Note: Perm_p_value of t test is obtained using '19' permutations. The smallest attainable p-value is 1/(19 + 1) = 0.05.
      
      ANOVA:
                 Df   Sum Sq    Mean Sq F value Pr(>F) Perm_p_value
      Days        1 162702.7 162702.652 71.4644 <1e-04         0.05
      Residuals 178 405251.6   2276.695                            
      
      Note: Perm_p_value of F test is obtained using '19' permutations. The smallest attainable p-value is 1/(19 + 1) = 0.05.
      
      Advice:
      ! Only 19 resamples: p-values at the floor (0.05) only say 'as small as can
        be shown' (Monte Carlo SE about 0.05 at p = 0.05). Use nperm = 999 or
        more (4999 for reporting); 'seed' makes it reproducible.
      

