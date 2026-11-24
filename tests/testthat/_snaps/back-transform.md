# a user's trans is used for models with a link (as covariatemeans()), with a message; trans = I keeps the link scale

    Code
      r1 <- predictmeans(fp, "g", trans = function(x) x^2)
    Message
      'trans' is used to back-transform the means instead of the model's inverse link.

# bias_adjust: one rule for predictmeans() and covariatemeans() (pm_bias_setup())

    Code
      predictmeans:::pm_bias_setup("yes", NULL, NULL, TRUE)
    Condition
      Error:
      ! 'bias_adjust' must be TRUE, FALSE or "taylor".

