# an invalid Kenward-Roger df stops, and a failed adjustment falls back to Satterthwaite with a message

    Code
      predictmeans:::pm_kr_ddf_formula(A1 = 0.85, A2 = 0.1, q = 1)
    Condition
      Error:
      ! Kenward-Roger df could not be computed (-1.669): the approximation breaks down for this model, e.g. a variance parameter near its boundary or a poorly identified one. Use df_method = "Satterthwaite" or "KR-expected".

---

    Code
      dfm <- predictmeans:::pm_df_method(fm, "auto")
    Message
      Kenward-Roger adjustment failed (simulated failure); using Satterthwaite df and unadjusted SEs.

# KR: observed information that is not positive definite falls back to the expected one, with a message

    Code
      ok <- predictmeans:::pm_kr_info_check(matrix(c(1, 2, 2, 1), 2), "observed")
    Message
      Kenward-Roger: the observed information of the variance parameters is not positive definite (smallest eigenvalue -1); using the expected information (as "KR-expected").

