## The stages of predictmeans(), in this order. Each is an ordinary function
## that takes the results of the stages before it and returns a named list:
##   opt  pm_setup()     the arguments as used: checked and resolved options,
##                       link or response back-transformation
##   est  pm_estimate()  pm_inference(), Kmatrix(), pm_estimate_means(): K
##                       matrix, predicted means, SEs, legacy tables
##   cmp  pm_compare()   pm_compare_sed() or pm_sed_summary(), pm_compare_df(),
##                       pm_pair_table() (every statistic of the pairs, once),
##                       pm_lsd_crit(), pm_comparison_intervals(),
##                       pm_pairs_overall() / pm_pairs_atvar() (the printed
##                       matrices and letters, from the pair table)
##   tab  pm_tables()    pm_mean_table(), pm_back_transform(),
##                       pm_final_mean_table(): the mean table, the response
##                       scale, letters
## pm_flat() puts the four lists together (a later stage's value replaces an
## earlier one: 'atvar', 'mt' and 'bkmt' are updated by pm_compare()) for the
## reader of everything: pm_assemble() (the stable tables of R/tidy.R and
## the legacy list). Plots are drawn from the result by autoplot().


## 'atvar' checked: NULL, "", "NULL" mean none; "A:B" is split; each name
## must be a factor of 'modelterm', and at least one factor must be left to
## compare. NULL for a one-factor term. The user's order is kept (it orders
## the groups).
pm_check_atvar <- function(atvar, vars) {
  if (is.null(atvar) || length(vars) == 1) return(NULL)
  if (!is.character(atvar)) stop("'atvar' must be a character vector of factor names.", call. = FALSE)
  atvar <- unique(trimws(unlist(strsplit(atvar, ":", fixed = TRUE))))
  atvar <- atvar[!atvar %in% c("", "NULL")]
  if (!length(atvar)) return(NULL)
  bad <- setdiff(atvar, vars)
  if (length(bad)) {
    stop("'atvar' must name factors of modelterm '", paste(vars, collapse = ":"), "'; not: '",
         paste(bad, collapse = "', '"), "'.", call. = FALSE)
  }
  if (length(atvar) == length(vars)) {
    stop("'atvar' names every factor of modelterm: no pairs would be left to compare within its levels.",
         call. = FALSE)
  }
  atvar
}

## Arguments, options and link function: the arguments of
## predictmeans() ('a') as they
## are used by the later stages.
pm_setup <- function(a) {
  model <- a$model
  modelterm <- a$modelterm
  pm_check_modelterm(model, modelterm)
  # (class probabilities are not additive across factors that do not interact; their predictors are)
  if (!pm_cat_is(model) || identical(a$mean_type, "link")) pm_note_not_a_term(model, modelterm)
  adj <- pm_match_adj(a$adj, westfall = TRUE)

  trans <- if (is.null(a$trans)) I else a$trans
  model <- pm_means_model(model)
  # multinomial and non-parallel fits: stacked class means through a carrier (R/categorical.R)
  mean_type <- if (is.null(a$mean_type)) "conditional" else a$mean_type
  average <- a$average
  if (pm_cat_is(model)) {
    ck <- pm_cat_check(model, a)
    mean_type <- ck$mean_type
    model <- pm_cat_carrier(model, mean_type, strsplit(modelterm, ":", fixed = TRUE)[[1]], ck$parts)
    average <- "link"          # (pm_cat_est() averages the probabilities itself; settings say "response")
  } else if (mean_type %in% c("prob_diff", "link")) {
    stop("mean_type = \"", mean_type, "\" is for multinomial and non-parallel fits (nnet::multinom, VGAM::vglm).",
         call. = FALSE)
  }

  # Models with a link function other than the identity: means are
  # back-transformed with its inverse (pm_linkinv(); see ?pm_coef to support
  # other model classes), unless the user gives 'trans' (as covariatemeans(),
  # pm_cov_back()): then theirs is used, with a message, and trans = I keeps
  # the link scale. (Before 2.0.0 a user's 'trans' was silently replaced.)
  linkinv <- pm_link_trans(model)
  user_trans <- !isTRUE(a$trans_default) && !is.null(a$trans)
  link_scale <- !is.null(linkinv) && !user_trans
  trans_by <- if (identical(trans, I)) NA_character_ else "trans"
  if (link_scale) {
    trans <- linkinv
    trans_by <- "the inverse link"
  } else if (!is.null(linkinv) && !identical(trans, I)) {
    message("'trans' is used to back-transform the means instead of the model's inverse link.")
  } else if (isTRUE(a$trans_default)) {
    # a transformed response in the model formula (log(y), sqrt(y), ...,
    # R/back-transform.R): back-transformed too, as for a link; trans = I
    # keeps the model scale only
    resp_trans <- pm_response_trans(model)
    if (!is.null(resp_trans)) {
      trans <- resp_trans$fun
      trans_by <- paste0("the inverse of ", resp_trans$label, " (the response in the model formula)")
      message("The response is ", resp_trans$label, ": the means are also given back-transformed to the original ",
              "scale (trans = I for the model scale only).")
    }
  }

  # Poisson models: the back-transformed means are counts
  count <- FALSE
  if (inherits(model, c("glm", "glmerMod", "glmmTMB")) &&
      stats::family(model)$family %in% c("poisson", "quasipoisson")) {
    count <- TRUE
  }

  vars <- unlist(strsplit(modelterm, "\\:"))
  # Model frame; for model classes without one, the data recovered by
  # Kmatrix() are used instead (pm_estimate()).
  mdf <- tryCatch(model.frame(model), error = function(e) NULL)

  # predictmeans() tabulates by factors: a covariate is set with 'covariate',
  # and means along a covariate come from covariatemeans() or Kmatrix()
  if (!is.null(mdf)) {
    not_fac <- intersect(vars, names(mdf)[!vapply(mdf, function(z) is.factor(z) || is.character(z) || is.logical(z),
                                                   logical(1))])
    if (length(not_fac)) {
      stop("'", paste(not_fac, collapse = "', '"), "' ", if (length(not_fac) > 1) "are covariates" else "is a covariate",
           ": predictmeans() gives means for the levels of factors. Set its value with covariate = list(",
           not_fac[1], " = ), or use covariatemeans() or Kmatrix() for means along a covariate.", call. = FALSE)
    }
  }

  # 'compare': comparisons asked for (adj, atvar, permlist, letter_ci,
  # pairwise = TRUE): the pairs are tabulated and adjusted whatever the
  # number of means. 'pairwise': only what is printed (the p-value
  # matrices and the letters); NULL (the default) shows them when
  # comparisons are asked for, as 1.x; FALSE hides them.
  atvar <- pm_check_atvar(a$atvar, vars)
  permlist <- a$permlist
  if (pm_perm_given(permlist)) {
    # "tukey", "mvt": single-step max-|t|; "westfall": step-down (pm_perm_maxT())
    if (adj == "scheffe") {
      stop("adj = \"scheffe\" is not available for permutation p-values; use another 'adj'.", call. = FALSE)
    }
  } else if (adj == "westfall") {
    stop("adj = \"westfall\" (step-down max-|t|) needs permutation p-values ('permlist'); without them use ",
         "\"mvt\" (single step) or \"holm\".", call. = FALSE)
  }
  # letters by overlap of the comparison intervals (R/comparison.R); 'level'
  # is the significance level of the comparisons, as everywhere
  if (a$letter_ci) adj <- "none"
  if (!is.null(a$pairwise) && !(is.logical(a$pairwise) && length(a$pairwise) == 1 && !is.na(a$pairwise))) {
    stop("'pairwise' must be TRUE, FALSE or NULL.", call. = FALSE)
  }
  compare <- isTRUE(a$pairwise) || adj != "none" || !is.null(atvar) || pm_perm_given(permlist) || a$letter_ci
  pairwise <- if (is.null(a$pairwise)) compare else a$pairwise
  if (a$letter_ci) pairwise <- TRUE                   # the letters asked for are shown


  # 'level' and 'slevel' are the user's level everywhere; the printed LSDs
  # take the critical value of 'adj' at this level (pm_lsd_crit()) and the
  # mean table with a permlist (model t limits, pm_perm_ci_df()) 'ci_level' (pm_tables())
  level <- a$level

  # zero-inflated glmmTMB fits: the overall means (1 - pi) mu (R/zero-inflation.R)
  if (identical(mean_type, "overall")) {
    if (is.null(pm_zi_component(model))) {
      message("mean_type = \"overall\" applies to glmmTMB models with zero inflation; this model has none, so the ",
              "means are the usual ones.")
      mean_type <- "conditional"
    } else {
      fam <- stats::family(model)
      if (!identical(fam$link, "log") || grepl("truncated", fam$family)) {
        stop("mean_type = \"overall\" is available for zero-inflated count models with a log link (poisson, nbinom1, ",
             "nbinom2, ...), not ", fam$family, " (", fam$link, ").", call. = FALSE)
      }
      if (!is.null(a$permlist) || !is.null(a$vcov) || isTRUE(a$bias_adjust) || identical(a$bias_adjust, "taylor")) {
        stop("mean_type = \"overall\" cannot be combined with 'permlist', 'vcov' or 'bias_adjust'.", call. = FALSE)
      }
    }
  }

  # clmm class probabilities: bias_adjust = TRUE averages them over the random effects (pm_ord_pop())
  ord_pop <- inherits(model, "clmm") && isTRUE(a$bias_adjust) && mean_type %in% c("prob", "mean_class")
  if (ord_pop && !is.null(a$sigma)) {
    stop("'sigma' is not used for population-average class probabilities of a clmm: they average over the ",
         "random-effect variances of the fit.", call. = FALSE)
  }
  if (inherits(model, pm_ord_classes)) {
    # (clmm: bootstrap p-values of the latent pairs from bootmodels(), B4 G3)
    pm_ord_check(model, trans = if (isTRUE(a$trans_default)) NULL else a$trans,
                 bias_adjust = if (ord_pop) FALSE else a$bias_adjust, sigma = a$sigma, vcov = a$vcov,
                 permlist = if (!inherits(model, "clmm")) a$permlist)
  }
  if (mean_type %in% c("prob", "mean_class") && !inherits(model, c(pm_ord_classes, "pm_catfit"))) {
    stop("mean_type = \"", mean_type, "\" is for models of a class response (MASS::polr, ordinal::clm, ordinal::clmm, ",
         "VGAM::vglm, nnet::multinom).", call. = FALSE)
  }
  list(model = model, modelterm = modelterm, data = a$data, vars = vars, mdf = mdf, mean_type = mean_type,
       pairwise = pairwise, compare = compare, atvar = atvar, adj = adj, Df = a$Df, level = level, slevel = level,
       covariate = a$covariate, cov_reduce = a$cov_reduce, exclude = a$exclude, weights = a$weights,
       meandecr = if (is.logical(a$meandecr)) a$meandecr else NULL, letter_ci = a$letter_ci,
       trans = trans, trans_by = trans_by, link_scale = link_scale, count = count,
       bias_adjust = if (ord_pop) FALSE else a$bias_adjust, ord_pop = ord_pop, sigma = a$sigma, average = average,
       permlist = permlist, ncore = a$ncore, ndecimal = 4L,
       df_method = a$df_method, vcov = a$vcov, vcov_adjust = a$vcov_adjust,
       vcov_adjust_default = a$vcov_adjust_default, equivalence = a$equivalence)
}

## K matrix, predicted means, standard errors and the legacy mean/SE tables.
pm_estimate <- function(opt) {
  # (the grid of all combinations is only needed to average on the response scale)
  ctr.matrix <- Kmatrix(opt$model, opt$modelterm, opt$covariate, data=opt$data,
                        cov_reduce=opt$cov_reduce, exclude=opt$exclude, weights=opt$weights,
                        refgrid=identical(opt$average, "response") || identical(opt$mean_type, "overall") ||
                          (inherits(opt$model, "pm_catfit") && opt$mean_type != "link"))
  KK <- ctr.matrix$K
  if (anyDuplicated(ctr.matrix$fctnames[seq_along(opt$vars)])) {
    stop("predictmeans() gives one mean for each combination of the levels of '", opt$modelterm, "': give one value ",
         "for each covariate in 'covariate', or use covariatemeans() or Kmatrix() for several values.", call. = FALSE)
  }
  mdf <- if (is.null(opt$mdf)) ctr.matrix$data else opt$mdf
  inf <- pm_inference(opt$model, opt$df_method, opt$vcov_adjust, opt$vcov_adjust_default, opt$vcov, opt$permlist,
                      term_df = TRUE)
  mp <- inf$mp
  est <- pm_estimate_means(KK, ctr.matrix$fctnames, mdf, opt$vars, mp, opt$ndecimal)
  out <- list(ctr.matrix = ctr.matrix, KK = KK, mdf = mdf, response = ctr.matrix$response, inf = inf, mp = mp,
       vcov_adjust = inf$vcov_adjust, vcov_label = inf$vcov_label, vcovm = mp$vcov,
       K = est$K, label = est$label, pm = est$pm, ses = est$ses,
       ndf = est$ndf,                 # cell counts, used by the stable tables
       mt = est$mt, bkmt = est$mt,    # bkmt: for back transformed
       mean.table = est$mean.table, se.table = est$se.table)
  # the overall means of a zero-inflated model, on the log scale (R/zero-inflation.R)
  if (identical(opt$mean_type, "overall")) out <- pm_zi_overall_est(opt, out)
  if (inherits(opt$model, "pm_catfit")) out <- pm_cat_est(opt, out)       # R/categorical.R
  out
}

## Predicted means, their SEs and the legacy mean/SE tables (pure function).
## KK: K matrix of all cells; fctnames: their labels; cells with no data are
## dropped. Returns list(K, label, pm, ses, mt, ndf (cell counts), mean.table,
## se.table).
pm_estimate_means <- function(KK, fctnames, mdf, vars, mp, ndecimal) {
  label <- fctnames
  rownames(label) <- rownames(KK)
  n.table <- table(mdf[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model

  K <- pm_K_cols(KK, mp)          # columns of the coefficients (R/mymodelparm.R)
  est_flag <- attr(KK, "estimable")

  # Combinations of levels without data are left out, unless their mean can
  # be estimated (a combination of factors that do not interact in the model)
  if (any(ndf$Freq==0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]
    keep <- rownames(K) %in% rnTrt
    if (!is.null(est_flag)) keep <- keep | (est_flag[rownames(K)] %in% TRUE)
    K <- K[keep, , drop=FALSE]
    label <- label[keep, , drop=FALSE]
  }

  pm <- K %*% mp$coef
  ses <- sqrt(rowSums((K %*% mp$vcov) * K))
  # Means that cannot be estimated (they depend on aliased coefficients) are
  # NA; their differences may still be estimable and are kept.
  if (!is.null(est_flag)) {
    non_est <- !est_flag[rownames(K)]
    non_est[is.na(non_est)] <- FALSE
    if (any(non_est)) {
      pm[non_est] <- NA
      ses[non_est] <- NA
      message("Means that cannot be estimated are NA: ", paste(rownames(K)[non_est], collapse = ", "),
              " (their estimable differences are kept; see vignette(\"missing-cells\")).")
    }
  }
  mt <- data.frame(pm, ses, label)

  # xtabs() drops NA rows and then fails when levels go unused: give it 0s
  # for the non-estimable means, which are set back to NA below
  tabs <- local({
    mt$pm[is.na(mt$pm)] <- 0
    mt$ses[is.na(mt$ses)] <- 0
    list(round(xtabs(pm ~ ., mt[, c("pm", vars)], drop.unused.levels = TRUE), ndecimal),
         round(xtabs(ses ~ ., mt[, c("ses", vars)], drop.unused.levels = TRUE), ndecimal+1))
  })
  mean.table <- tabs[[1]]
  se.table <- tabs[[2]]
  no_mean <- !(xtabs(~ ., mt[, vars, drop = FALSE], drop.unused.levels = TRUE) > 0)   # combinations left out
  mean.table[no_mean] <- NA
  se.table[no_mean] <- NA
  if (anyNA(mt$pm)) {   # non-estimable means: xtabs() would show 0
    na_cells <- as.matrix(mt[is.na(mt$pm), vars, drop = FALSE])
    mean.table[na_cells] <- NA
    se.table[na_cells] <- NA
  }

  if (length(vars) > 1) {
    varsnlevel <- numeric(0)

    for (i in vars){
      varsnlevel[i] <- nlevels(mdf[, i])
    }

    tbvars <- names(sort(varsnlevel, decreasing = TRUE))
    mean.table <- ftable(mean.table, row.vars =tbvars[1], col.var=tbvars[-1])
    se.table <- ftable(se.table, row.vars =tbvars[1], col.var=tbvars[-1])
  }

  if (length(na.omit(unique(se.table)))==1) {
    se.table <- min(se.table, na.rm=TRUE)
    names(se.table) <- "All means have the same SE"
  }
  list(K = K, label = label, pm = pm, ses = ses, mt = mt, ndf = ndf,
       mean.table = mean.table, se.table = se.table)
}

## SEDs, LSDs, degrees of freedom and pairwise comparisons.
pm_compare <- function(opt, est) {
  K <- est$K
  KK <- est$KK
  label <- est$label
  vars <- opt$vars
  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)
  perm <- pm_perm_given(opt$permlist)
  # what the stage gives; 'atvar', 'mt', 'bkmt' and 'Df' start as they came
  atvar <- opt$atvar
  mt <- est$mt
  bkmt <- est$bkmt
  Df <- opt$Df
  SED.out <- LSD <- NA
  cmp_int <- NULL        # comparison intervals (R/comparison.R)
  rK <- cm <- dses <- dses.df <- dses.m <- KKvarndiff <- labKK <- pair_over <- NULL
  pm_df_raw <- Df_diff <- pairDf <- df_user <- pair_tab <- nsim <- ci_df <- NULL
  t.p.valuem <- t.p.valuep <- LSDm <- t.p.valuemGrp <- p_valueMatrix <- NULL
  resvar <- atvar.levels <- rcnplotm <- listlength <- nrow.pm <- pmlist <- pmlistTab <- pmlistLetter <- NULL
  # means over empty cells (R/empty-cells.R): their pairs are compared on the
  # cells both means have, so nothing read off the means themselves
  # (comparison intervals, letters from them) can agree with their tests
  ec <- est$ctr.matrix$empty_cells
  affected <- stats::setNames(if (is.null(ec)) rep(FALSE, nK) else unname(ec$affected[rnK] %in% TRUE), rnK)
  letter_ci <- opt$letter_ci
  result <- function() {
    list(nK = nK, rnK = rnK, perm = perm, atvar = atvar, mt = mt, bkmt = bkmt, Df = Df, SED.out = SED.out, LSD = LSD,
         cmp_int = cmp_int, rK = rK, cm = cm, dses = dses, dses.df = dses.df, dses.m = dses.m,
         KKvarndiff = KKvarndiff, labKK = labKK, pair_over = pair_over, pm_df_raw = pm_df_raw, Df_diff = Df_diff,
         ci_df = ci_df,
         pairDf = pairDf, df_user = df_user, pair_tab = pair_tab, nsim = nsim, t.p.valuem = t.p.valuem,
         t.p.valuep = t.p.valuep, LSDm = LSDm, t.p.valuemGrp = t.p.valuemGrp, p_valueMatrix = p_valueMatrix,
         resvar = resvar, atvar.levels = atvar.levels, rcnplotm = rcnplotm, listlength = listlength,
         nrow.pm = nrow.pm, pmlist = pmlist, pmlistTab = pmlistTab, pmlistLetter = pmlistLetter,
         letter_ci = letter_ci, affected = affected)
  }
  if (nK == 1) return(result())

  # The pairs are tabulated (one row per pair) when comparisons are asked
  # for (pairwise = TRUE, adj, atvar, letter_ci, permlist, equivalence), and
  # otherwise for up to pm_max_means() means. Above that (several hundred
  # lines of a breeding trial) the SEDs and LSDs are summarised from the
  # covariance of the means, without a row per pair (pm_sed_summary()).
  # Means over empty cells and non-estimable means need each pair's
  # contrast, so their SEDs take the first route too.
  tab_pairs <- isTRUE(opt$compare) || !is.null(opt$equivalence) || nK <= pm_max_means()
  by_pair <- tab_pairs || !is.null(est$ctr.matrix$empty_cells) || !is.null(attr(KK, "null_proj"))
  if (!tab_pairs) {
    message(nK, " means of '", opt$modelterm, "' (", nK*(nK - 1)/2, " pairs): the SEDs and LSDs are summarised and ",
            "the pairs are not tabulated. pairwise = TRUE (or 'adj') compares them.")
  }
  if (by_pair) {
    s <- pm_compare_sed(K, KK, label, est$ctr.matrix$fctnames, vars, est$mp$coef, est$vcovm, atvar,
                        cells = est$ctr.matrix$empty_cells)
    rK <- s[["rK"]]; cm <- s[["cm"]]; dses <- s[["dses"]]; dses.df <- s[["dses.df"]]
    KKvarndiff <- s[["KKvarndiff"]]; labKK <- s[["labKK"]]; pair_over <- s[["pair_over"]]
  } else {
    s <- pm_sed_summary(K, label, vars, est$vcovm, atvar)
  }
  SED.out <- s[["SED.out"]]
  dses.m <- s[["dses.m"]]
  if (inherits(opt$model, "pm_catfit") && opt$model$var %in% colnames(dses.m)) {
    # class probabilities: the SEDs of pairs within a class only (R/categorical.R)
    SED.out[c("Max.SED", "Min.SED", "Aveg.SED")] <- dses.m[c("Max.SED", "Min.SED", "Aveg.SED"), opt$model$var]
    # the same level of another factor means two classes: those columns are left out (also of the LSDs)
    dses.m <- dses.m[, opt$model$var, drop = FALSE]
    attr(SED.out, "For the Same Level of Factor") <- dses.m
  }
  atvar <- s[["atvar"]]              # NULL when a factor has no pairs within its levels (nested factors)

  if (!perm) {
    d <- pm_compare_df(est$inf, opt$modelterm, K, rK, mt, Df)
    Df <- d[["Df"]]; mt <- d[["mt"]]; pm_df_raw <- d[["pm_df_raw"]]; Df_diff <- d[["Df_diff"]]
    pairDf <- d[["pairDf"]]; df_user <- d[["df_user"]]
    bkmt <- mt
  } else {
    # permutation tests are of comparisons; the limits of each mean are the
    # model's t limits at conf_level (its SEs are the model's), with the
    # model's df (pm_perm_ci_df()); 'Df' of the tables stays NA
    ci_df <- pm_perm_ci_df(est$inf, opt$modelterm, K, Df)
  }
  no_atvar <- is.null(atvar) || all(atvar%in%c("NULL", ""))
  pair_cor <- NULL
  if (tab_pairs) {
    # every statistic of the pairs, once (the pairs table, the printed
    # matrices and the letters read it)
    pv <- pm_pair_pvalues(opt$model, rK, cm, dses, if (perm) NULL else Df_diff, opt$adj, perm, opt$permlist,
                          opt$ncore, est$vcovm)
    nsim <- pv[["nsim"]]
    pair_cor <- pv[["pair_cor"]]
    pair_tab <- pm_pair_table(label, if (no_atvar) NULL else atvar, rK, cm, dses, if (perm) NULL else Df_diff,
                              pv[["t.p.values"]], opt$adj, perm, opt$level, pair_cor, pv[["perm_T"]])
  }
  # The printed LSDs are the SEDs times the critical value that goes with
  # the p-values shown (pm_lsd_crit(), R/p-adjust.R): the half-width of the
  # limits that agree with them, with the same family of pairs.
  lsd <- if (perm) list(crit = NA_real_, note = NULL) else
    pm_lsd_crit(opt$adj, opt$level, Df, if (tab_pairs) attr(pair_tab, "families"), pair_cor)
  # means over empty cells whose pairs are not the differences of the shown
  # means (R/comparison.R): read their pairs, not their means
  misread <- if (tab_pairs) pm_misread_means(pair_tab, mt$pm, affected) else affected
  if (any(misread) && !perm) {
    lsd$note <- paste(c(lsd$note, paste0("The pairs of ", paste(rnK[misread], collapse = ", "), " (means over empty ",
                                         "cells) are compared on common cells: compare their pairs, not the differences ",
                                         "of their means, with the LSD.")), collapse = " ")
  }
  LSD <- pm_compare_lsd(SED.out, dses.m, vars, lsd$crit, Df, opt$slevel, opt$ndecimal, perm, lsd$note)
  # with 'atvar' only the pairs within its levels are compared: their LSDs,
  # each pair's own critical value (its df, its family), per level and over
  # all the levels (pm_family_lsd())
  if (tab_pairs && !perm && !no_atvar) {
    fl <- pm_family_lsd(pair_tab)
    if (!is.null(fl)) {
      ov <- attr(fl, "overall")
      LSD[c("Max.LSD", "Min.LSD", "Aveg.LSD")] <- round(ov[c("max", "min", "mean")], opt$ndecimal + 1)
      attr(LSD, "Within levels of atvar") <- fl
      attr(LSD, "Degree of freedom") <- round(ov[["df"]], 2)
    }
  }
  # one LSD only when all the LSDs compared are equal (to 1e-8, relative);
  # otherwise the range is the result, and a note says to read the pairs
  if (!perm && is.numeric(LSD) && all(c("Max.LSD", "Min.LSD") %in% names(LSD))) {
    rng <- unname(LSD[c("Min.LSD", "Max.LSD")])
    if (all(is.finite(rng))) {
      equal <- abs(diff(rng)) <= 1e-8 * max(abs(rng), 1) + 10^-(opt$ndecimal + 1)
      attr(LSD, "equal") <- equal
      if (!equal) {
        attr(LSD, "Note") <- paste(c(attr(LSD, "Note"), paste0(
          "The LSDs differ between pairs (", format(signif(rng[1], 4)), " to ", format(signif(rng[2], 4)),
          "): no single LSD tests every pair; use the p-values, or each pair's LSD (pairs$lsd).")), collapse = " ")
      }
    }
  }
  if (tab_pairs) cmp_int <- pm_comparison_intervals(pair_tab, mt$pm, mt$ses, if (perm) NULL else mt$Df, opt$level,
                                                    opt$adj, affected)
  if (letter_ci && any(misread)) {
    message("letter_ci: ", paste(rnK[misread], collapse = ", "), " average", if (sum(misread) == 1) "s",
            " over empty cells and ", if (sum(misread) == 1) "is" else "are", " compared on common cells, which no ",
            "interval around a mean can show; the letters are from the p-values.")
    letter_ci <- FALSE
  }

  # the printed matrices and the letters, from the pair table
  if (opt$pairwise) {
    if (no_atvar) {
      o <- pm_pairs_overall(opt$model, rnK, mt, pair_tab, nK, Df, pairDf, opt$adj, perm, nsim, lsd$crit, opt$slevel,
                            opt$meandecr)
      t.p.valuem <- o[["t.p.valuem"]]; t.p.valuep <- o[["t.p.valuep"]]; LSDm <- o[["LSDm"]]
      t.p.valuemGrp <- o[["t.p.valuemGrp"]]; p_valueMatrix <- o[["p_valueMatrix"]]
    } else {
      o <- pm_pairs_atvar(dses.df, KKvarndiff, labKK, mt, est$mdf, vars, atvar, pair_tab, opt$slevel)
      resvar <- o[["resvar"]]; atvar.levels <- o[["atvar.levels"]]; rcnplotm <- o[["rcnplotm"]]
      bkmt <- o[["bkmt"]]; listlength <- o[["listlength"]]; nrow.pm <- o[["nrow.pm"]]; pmlist <- o[["pmlist"]]
      pmlistTab <- o[["pmlistTab"]]; pmlistLetter <- o[["pmlistLetter"]]; p_valueMatrix <- o[["p_valueMatrix"]]
    }
  }
  result()
}

## The pure functions of the compare stage. Each returns a named list of the
## values it creates (only those that exist).
pm_existing <- function(nms, env) mget(intersect(nms, ls(env, all.names = TRUE)), envir = env)

## Pairwise contrasts and their SEDs (overall and within levels of each
## factor). 'atvar' becomes NULL when a factor has no pairs within its levels
## (nested factors). 'cells' (Kmatrix()$empty_cells): pairs with a mean over
## empty cells are compared on the cells both means have (R/empty-cells.R);
## pair_over says which.
pm_compare_sed <- function(K, KK, label, fctnames, vars, coef, vcovm, atvar, cells = NULL) {
  nK <- nrow(K)
  rnK <- rownames(K)
  # pairs (i, j), i > j, in column order; level labels come from the label
  # tables, not from splitting row names at ":" (levels may contain ":")
  pp <- pm_pair_index(nK)
  labK <- pm_label_matrix(label, length(vars))

  nKK <- nrow(KK)
  ppKK <- pm_pair_index(nKK)
  labKK <- pm_label_matrix(fctnames, length(vars))

  KKvarndiff <- data.frame(labKK[ppKK$i, , drop=FALSE], labKK[ppKK$j, , drop=FALSE])

  rK <- K[pp$j, , drop = FALSE] - K[pp$i, , drop = FALSE]     # mean j minus mean i
  # differences that cannot be estimated (their components along the
  # non-estimable directions do not cancel) are NA, not computed with the
  # dropped coefficients taken as zero
  non_est_pair <- rep(FALSE, nrow(rK))
  null_proj <- attr(KK, "null_proj")
  if (!is.null(null_proj) && all(rnK %in% rownames(null_proj))) {
    np <- null_proj[rnK, , drop = FALSE]
    pd <- np[pp$j, , drop = FALSE] - np[pp$i, , drop = FALSE]
    attr(pd, "scale") <- attr(null_proj, "scale")
    non_est_pair <- !pm_estimable_rows(pd, rK)
  }
  # pairs with a mean over empty cells: on the cells both means have
  pair_over <- rep(NA_character_, nrow(rK))
  if (!is.null(cells)) {
    Cf <- matrix(0, length(pp$i), nrow(KK))
    Cf[, match(rnK, rownames(KK))] <- pm_pair_matrix(nK)
    cr <- pm_common_rows(cells, Cf, rK, colnames(K))
    rK <- cr$L
    non_est_pair[cr$sub] <- cr$na[cr$sub]
    pair_over <- cr$over
  }
  cm <- rK%*%coef
  dses <- sqrt(rowSums((rK %*% vcovm) * rK))       # the diagonal of rK V rK' only
  cm[non_est_pair] <- NA
  dses[non_est_pair] <- NA

  SED.out <- if (all(is.na(dses))) c(Max.SED = NA_real_, Min.SED = NA_real_, Aveg.SED = NA_real_) else
    c(Max.SED = max(dses, na.rm = TRUE), Min.SED = min(dses, na.rm = TRUE), Aveg.SED = mean(dses, na.rm = TRUE))
  dses.df <- data.frame(labK[pp$i, , drop=FALSE], labK[pp$j, , drop=FALSE], dses, pair_id = seq_along(pp$i))
  dses.m <- NULL

  if (length(vars) > 1) {
    dses.m <- matrix(0, nrow=3, ncol=length(vars))
    colnames(dses.m) <- vars
    rownames(dses.m) <- c("Aveg.SED", "Min.SED", "Max.SED")

    for (i in seq_along(vars)) {
      varsindx <- as.character(dses.df[,i])==as.character(dses.df[, length(vars)+i])
      if(any(varsindx)){  # in case of A/B
        dv <- dses.df[varsindx, "dses"]
        dses.m[,i] <- if (all(is.na(dv))) NA else c(mean(dv, na.rm = TRUE), min(dv, na.rm = TRUE), max(dv, na.rm = TRUE))
      }else{
        dses.m[,i] <- NA
        atvar <- NULL
      }
    }
    attr(SED.out, "For the Same Level of Factor") <- dses.m
  }
  list(rK = rK, cm = cm, dses = dses, SED.out = SED.out, dses.df = dses.df,
       dses.m = dses.m, KKvarndiff = KKvarndiff, labKK = labKK, atvar = atvar, pair_over = pair_over)
}

## Above this many means, the pairs are tabulated only when comparisons are
## asked for, and no letter groups are made (pm_letter_groups()): one limit
## for both, options(predictmeans.max_means = ).
pm_max_means <- function() {
  n <- getOption("predictmeans.max_means", 100)
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 2) stop("options(predictmeans.max_means = ) must be a number of at least 2.", call. = FALSE)
  n
}

## The SEDs summarised without a row per pair (no comparisons asked for):
## the variance of mean j minus mean i is c_jj + c_ii - 2 c_ij, C = K V K'
## the covariance of the means, so all n(n-1)/2 SEDs come from one n x n
## matrix. Returns what the summaries of pm_compare_sed() return: SED.out
## (maximum, minimum, average; within the levels of each factor as an
## attribute), dses.m and 'atvar' (NULL when a factor has no pairs within its
## levels).
pm_sed_summary <- function(K, label, vars, vcovm, atvar) {
  C <- K %*% tcrossprod(as.matrix(vcovm), K)
  d <- diag(C)
  S <- sqrt(pmax(outer(d, d, "+") - 2*C, 0))
  low <- lower.tri(S)
  s3 <- function(x) c(max(x), min(x), mean(x))
  SED.out <- stats::setNames(s3(S[low]), c("Max.SED", "Min.SED", "Aveg.SED"))
  dses.m <- NULL
  if (length(vars) > 1) {
    labK <- pm_label_matrix(label, length(vars))
    dses.m <- matrix(0, nrow = 3, ncol = length(vars), dimnames = list(c("Aveg.SED", "Min.SED", "Max.SED"), vars))
    for (i in seq_along(vars)) {
      same <- low & outer(labK[, i], labK[, i], "==")
      if (any(same)) {
        dses.m[, i] <- s3(S[same])[c(3, 2, 1)]
      } else {                      # in case of A/B
        dses.m[, i] <- NA
        atvar <- NULL
      }
    }
    attr(SED.out, "For the Same Level of Factor") <- dses.m
  }
  list(SED.out = SED.out, dses.m = dses.m, atvar = atvar)
}

## The mean df of all the pairs of the rows of K, a block of pairs at a time
## (each block's contrasts are built, used and dropped).
pm_pair_df_mean <- function(inf, modelterm, K, block = 5000L) {
  pp <- pm_pair_index(nrow(K))
  n <- length(pp$i)
  tot <- 0
  n_ok <- 0L
  for (a in seq.int(1L, n, by = block)) {
    k <- a:min(n, a + block - 1L)
    d <- inf$df(K[pp$j[k], , drop = FALSE] - K[pp$i[k], , drop = FALSE], modelterm)
    tot <- tot + sum(d, na.rm = TRUE)
    n_ok <- n_ok + sum(!is.na(d))
  }
  tot/n_ok
}

## Degrees of freedom of the means and of the pairs (no permutation tests):
## the user's Df, else from the call's inference choices (R/inference.R).
## Df: one df for the LSD tables (lmer: the term's df; df varying by row:
## their mean); pairDf: the df printed with the legacy p-value matrix.
## rK = NULL (the pairs are not tabulated, pm_sed_summary()): the same Df,
## without the df of each pair.
pm_compare_df <- function(inf, modelterm, K, rK, mt, Df) {
  if (length(Df)) {
    mt$Df <- Df
    return(list(Df = Df, mt = mt, Df_diff = if (!is.null(rK)) rep_len(Df, nrow(rK)), pairDf = Df, df_user = TRUE))
  }
  pm_df_raw <- inf$df(K, modelterm)
  if (is.null(rK)) {
    term_df <- inf$per_row && !is.null(names(inf$mp$df)) && modelterm %in% names(inf$mp$df)
    Df <- if (term_df) unname(inf$mp$df[modelterm]) else if (inf$per_row) pm_pair_df_mean(inf, modelterm, K) else
      inf$df(K[1, , drop = FALSE] - K[2, , drop = FALSE], modelterm)
    if (isTRUE(Df == 0)) stop("You need provide Df for this model!")
    mt$Df <- round(pm_df_raw, 2)
    return(list(Df = Df, mt = mt, pm_df_raw = pm_df_raw))
  }
  Df_diff <- inf$df(rK, modelterm)
  # pairs that cannot be estimated (confounded factors) have NA df: leave them out
  if (inf$per_row) {
    Df <- if (!is.null(names(inf$mp$df)) && modelterm %in% names(inf$mp$df)) unname(inf$mp$df[modelterm]) else
      mean(Df_diff, na.rm = TRUE)
    pairDf <- round(mean(Df_diff, na.rm = TRUE), 2)
  } else {
    Df <- Df_diff[!is.na(Df_diff)][1]
    pairDf <- Df
  }
  if (isTRUE(Df == 0)) stop("You need provide Df for this model!")
  mt$Df <- round(pm_df_raw, 2)
  list(Df = Df, mt = mt, pm_df_raw = pm_df_raw, Df_diff = Df_diff, pairDf = pairDf)
}

## df of the confidence limits of each mean when the tests are permutation
## tests ('permlist'): the user's Df, else the model's df of each mean (the
## same df as without a permlist; Inf, normal limits, when there are none).
## Before 2.0.0 these limits were mean +/- 2 SE whatever conf_level was.
pm_perm_ci_df <- function(inf, modelterm, K, Df = NULL) {
  if (length(Df)) return(rep_len(as.numeric(Df), nrow(K)))
  d <- tryCatch(as.numeric(inf$df(K, modelterm)), error = function(e) rep(NA_real_, nrow(K)))
  d[is.na(d) | d <= 0] <- Inf
  d
}

## Critical value of the confidence limits of the means: t at 'slevel' with
## the df of each mean (ci_df with a permlist, see pm_perm_ci_df()).
pm_ci_crit <- function(slevel, df) stats::qt(1 - slevel/2, df)

## LSDs from the SEDs: SED x 'crit' (pm_lsd_crit(): the critical value that
## goes with the p-values), or 2 x SED with permutation tests.
pm_compare_lsd <- function(SED.out, dses.m, vars, crit, Df, slevel, ndecimal, perm, note = NULL) {
  if (!perm) {
    LSD <- round(crit * SED.out, ndecimal+1)
    names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
    attr(LSD, "For the Same Level of Factor") <- NULL

    if (length(vars) > 1) {
      rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
      attr(LSD, "For the Same Level of Factor") <- round(crit * dses.m, ndecimal+1)
    }

    attr(LSD, "Significant level") <- slevel
    attr(LSD, "Degree of freedom") <- round(Df, 2)
    if (!is.null(note)) attr(LSD, "Note") <- note
  }else{
    LSD <- round(2 * SED.out[1:3], ndecimal+1)
    names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")

    if (length(vars) > 1) {
      rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
      attr(LSD, "For the Same Level of Factor") <- round(2 * dses.m, ndecimal+1)
    }

    attr(LSD, "Note") <- "This is a approximate LSD (i.e. 2*SED) at 0.05 level."
  }
  LSD
}

## LSDs of the pairs compared within each level of 'atvar': each pair's SED
## times its critical value (crit_adj when 'adj' has simultaneous limits,
## else t with the pair's own df), summarised per level (by), with the number
## of pairs, the SED range and the critical values; attribute "overall": the
## same over all the levels (max, min, mean, df = mean df). NULL without
## within-level pairs.
pm_family_lsd <- function(pair_tab) {
  pt <- pair_tab[!is.na(pair_tab$by) & !is.na(pair_tab$estimate), , drop = FALSE]
  if (!nrow(pt)) return(NULL)
  crit <- if (!is.null(pt$crit_adj)) pt$crit_adj else pt$crit
  lsd <- crit * pt$std.error
  by <- factor(pt$by, levels = unique(pt$by))
  f <- function(x, g) as.vector(tapply(x, by, g))
  out <- data.frame(by = levels(by), pairs = as.vector(table(by)),
                    min.SED = f(pt$std.error, min), max.SED = f(pt$std.error, max),
                    crit = f(crit, mean), min.LSD = f(lsd, min), aveg.LSD = f(lsd, mean), max.LSD = f(lsd, max),
                    stringsAsFactors = FALSE)
  attr(out, "overall") <- c(max = max(lsd), min = min(lsd), mean = mean(lsd), df = mean(pt$df))
  out
}

## t values and raw p-values of the pairs: t with each pair's df, or by
## permutation (perm_T: the observed and permuted |t|, for max-|t|
## adjustments and permutation limits). pair_cor: the correlation of the
## pairs, for adj = "mvt".
pm_pair_pvalues <- function(model, rK, cm, dses, Df_diff, adj, perm, permlist, ncore, vcovm = NULL) {
  t.v <- cm/dses
  pair_cor <- perm_T <- nsim <- NULL
  if (!perm && adj == "mvt") pair_cor <- stats::cov2cor(as.matrix(rK %*% tcrossprod(vcovm, rK)))
  if (!perm) {
    t.p.values <- 2 * stats::pt(-abs(as.vector(t.v)), rep_len(as.numeric(Df_diff), length(t.v)))
  } else {
    nsim <- pm_perm_n(permlist)
    fam <- pm_perm_family(permlist, rK, model, abs(as.vector(t.v)), ncore, joint = TRUE)
    t.p.values <- fam$p
    perm_T <- fam[intersect(c("obs", "T", "se"), names(fam))]   # se: robust SEs (robust permlist)
  }
  pm_existing(c("t.v", "t.p.values", "nsim", "pair_cor", "perm_T"), environment())
}

## The level of 'atvar' a pair lies within, NA for a pair across its levels
## (pairs of the means in the rows of 'label', in pm_pair_index() order).
pm_pair_by <- function(label, pp, atvar) {
  by <- rep(NA_character_, length(pp$i))
  if (is.null(atvar)) return(by)
  A <- label[pp$j, atvar, drop = FALSE]
  B <- label[pp$i, atvar, drop = FALSE]
  same <- Reduce(`&`, lapply(atvar, function(v) as.character(A[[v]]) == as.character(B[[v]])))
  by[same] <- do.call(paste, c(lapply(A[same, , drop = FALSE], as.character), sep = " : "))
  by
}

## Every statistic of the pairs, computed once (pure function): the pairs
## table, the printed p-value and LSD matrices and the letter groups all
## read it. One row per pair in pm_pair_index() order (mean j minus mean i):
## by (the 'atvar' level, NA across levels), estimate, std.error, df,
## statistic, p.value, p.adj (within all the pairs, or within each level of
## 'atvar'; R/p-adjust.R), crit (critical value of each pair's limits at
## 'level': t, or the permutation quantile of |t|), crit_adj (of the
## simultaneous limits; only for methods that have them) and ci_se (the SEs
## of the limits: HC3 with a robust permlist).
pm_pair_table <- function(label, atvar, rK, cm, dses, Df_diff, p, adj, perm, level, pair_cor = NULL,
                          perm_T = NULL) {
  pp <- pm_pair_index(nrow(label))
  cm <- as.vector(cm)
  dses <- as.vector(dses)
  n <- length(cm)
  pdf <- if (perm) rep(NA_real_, n) else rep_len(as.numeric(Df_diff), n)
  pdf[is.na(cm)] <- NA
  stat <- cm/dses
  by <- pm_pair_by(label, pp, atvar)
  fams <- pm_families(if (is.null(atvar)) rep("all", n) else by, est = !is.na(cm), ij = cbind(pp$j, pp$i), L = rK)
  out <- data.frame(j = pp$j, i = pp$i, by = by, estimate = cm, std.error = dses, df = pdf, statistic = stat,
                    p.value = as.vector(p), stringsAsFactors = FALSE)
  out$p.adj <- pm_adjust_p(adj, stat, out$p.value, pdf, fams, R = pair_cor, perm_T = perm_T)
  out$crit <- if (!perm) stats::qt(1 - level/2, pdf) else if (!is.null(perm_T)) pm_perm_crit(perm_T$T, level) else NA_real_
  crit_adj <- if (!perm || !is.null(perm_T)) pm_adjust_crit(adj, level, pdf, fams, R = pair_cor, T = perm_T$T)
  if (!is.null(crit_adj)) out$crit_adj <- crit_adj
  out$ci_se <- if (!is.null(perm_T$se)) perm_T$se else dses
  attr(out, "families") <- fams
  out
}

## Overall pairwise tables, from the pair table: p-values (t value above
## the diagonal), LSDs and the matrix used for letter groups.
pm_pairs_overall <- function(model, rnK, mt, pair_tab, nK, Df, pairDf, adj, perm, nsim, lsd_crit, slevel,
                             meandecr) {
  tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
  rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
  tvm[upper.tri(tvm)] <- pair_tab$statistic
  p_valueMatrix <- NULL
  t.p.valuem[upper.tri(t.p.valuem)] <- pair_tab$p.adj

        t.p.valuep <- t.p.valuem    # for plot
        t.p.valuem <- t(t.p.valuem) + tvm
        names(t.p.valuem) <- NULL
		    diag(t.p.valuem) <- 1

		    if (!perm) {

		      if (!inherits(model, "lmerMod")) {
		        attr(t.p.valuem, "Degree of freedom") <- pairDf
		      }

		      attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, p-value (adjusted by '",
                                            adj, "' method) below the diagonal", sep="")
          LSDm.up <- lsd_crit*pair_tab$std.error
          LSDm[upper.tri(LSDm)] <- LSDm.up
          Diffm[upper.tri(Diffm)] <- pair_tab$estimate
          LSDm <- t(LSDm)+Diffm
          names(LSDm) <- NULL
          attr(LSDm,"Significant level") <- slevel
          attr(LSDm,"Degree of freedom") <- Df
          attr(LSDm,"Note") <- paste("LSDs matrix has mean differences (row-col) above the diagonal, LSDs (adjusted by '",
                                     adj, "' method) below the diagonal", sep="")
      }else{
        attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, and ", nsim, " times permutation p-value (adjusted by '", adj, "' method) below the diagonal", sep="")
      } # end of if (is.null(permlist))

      if (!is.null(meandecr) && is.logical(meandecr)){
        groupRn <- rnK[order(mt$pm, decreasing = meandecr)]
      }else{
        groupRn <- NULL
      }

		  t.p.valuemGrp <- t.p.valuem
      t.p.valuemGrp[upper.tri(t.p.valuemGrp)] <- t(t.p.valuemGrp)[upper.tri(t.p.valuemGrp)]

      if (!is.null(groupRn)){
        t.p.valuemGrp <- t.p.valuemGrp[groupRn, groupRn]
      }

      if (nrow(t.p.valuep) > 2){
        p_valueMatrix <- round(t(t.p.valuep), 4)
      }

  pm_existing(c("t.p.valuem", "t.p.valuep", "LSDm", "t.p.valuemGrp", "p_valueMatrix"), environment())
}

## Pairwise tables within each level of 'atvar', from the pair table:
## adjusted p-values, letters and the means ordered by 'atvar'.
pm_pairs_atvar <- function(dses.df, KKvarndiff, labKK, mt, mdf, vars, atvar, pair_tab, slevel) {
  p_valueMatrix <- NULL
      for (i in which(vars%in%atvar)) { # To find rows relating atvar
        KKvarndiff <- KKvarndiff[which(as.character(KKvarndiff[, i]) == as.character(KKvarndiff[, length(vars) + i])), ]
      }

      atvar.df <- f_loj_krc(KKvarndiff, dses.df, by.x=names(KKvarndiff), by.y=names(KKvarndiff))
      names(atvar.df)[seq_along(vars)] <- vars

		  for (i in vars) {       # To ensure factor vars  have the same level as original
          atvar.df[,i] <- factor(atvar.df[,i], levels=levels(mdf[, i][[1]]))
      }

      # adjusted p-values within each atvar group (pm_pair_table())
      key_atv <- do.call("paste", c(lapply(atvar.df[, atvar, drop=FALSE], as.character), sep=" : "))
      atvar.df$adj.pvalue <- pair_tab$p.adj[atvar.df$pair_id]

      rnK.df <- as.data.frame(labKK) # To find the suitable names for pm
      colnames(rnK.df) <- vars

      for (i in vars) {       # To ensure factor vars  have the same level as original
          rnK.df[,i] <- factor(rnK.df[,i], levels=levels(mdf[, i][[1]]))
      }

      rnK.df  <- rnK.df [do.call(order, rnK.df[, atvar, drop=FALSE]),]
      resvar <- vars[!vars%in%atvar]   # The rest of vars rather than at var
      rnK.df <- rnK.df[, c(atvar, resvar)]   # To ensure the right matrix name later
      atvar.levels <- unique(do.call("paste", c(rnK.df[, atvar, drop=FALSE], sep=" : ")))
      resvar.levels <- unique(do.call("paste", c(rnK.df[, resvar, drop=FALSE], sep=" : ")))
      rcnplotm <- do.call("paste", c(rnK.df[, , drop=FALSE], sep=" : ")) # row col names of image plot

      # mean table arrange by atvar
      mt.atvar <- mt[do.call(order, mt[, atvar, drop=FALSE]),]
      # no Df column with permutation p-values
      bkmt <- mt.atvar[, intersect(c(atvar, resvar, "pm", "ses", "Df"), names(mt.atvar))]
      listlength <- length(atvar.levels)
      pmlist <- pmlistTab <- pmlistLetter <- vector("list", listlength)
      nrow.pm <- length(resvar.levels)

      # the two members of each pair: first in columns 1..nv, second in nv+1..2nv
      nv <- length(vars)
      res_pos <- match(resvar, vars)
      lev1 <- do.call("paste", c(lapply(atvar.df[, res_pos, drop=FALSE], as.character), sep=" : "))
      lev2 <- do.call("paste", c(lapply(atvar.df[, nv + res_pos, drop=FALSE], as.character), sep=" : "))

      for (i in 1:listlength) {           # p-values of this atvar group, placed by level name
        idx <- which(key_atv == atvar.levels[i])
        S <- matrix(NA_real_, nrow=nrow.pm, ncol=nrow.pm, dimnames=list(resvar.levels, resvar.levels))
        if (length(idx)) {
          S[cbind(lev1[idx], lev2[idx])] <- atvar.df$adj.pvalue[idx]
          S[cbind(lev2[idx], lev1[idx])] <- atvar.df$adj.pvalue[idx]
        }
        diag(S) <- 0
        atvar.pm <- S
        atvar.pm[lower.tri(atvar.pm)] <- 0
        pmlist[[i]] <- t(atvar.pm)
        outtab <- round(as.table(t(atvar.pm)), 4)
        outtab[outtab < 0.0001] <- "0.0001"
        outtab[col(outtab)==row(outtab)] <- 1.0000
        outtab[upper.tri(outtab)] <-""
        Grpatvar.pm <- t(atvar.pm)+ atvar.pm

        # levels of this atvar group with an observed mean
        in_grp <- do.call("paste", c(lapply(mt.atvar[, atvar, drop=FALSE], as.character), sep=" : ")) == atvar.levels[i]
        obs_lev <- do.call("paste", c(lapply(mt.atvar[in_grp, resvar, drop=FALSE], as.character), sep=" : "))
        Grpatvar.letterM <- pm_letter_groups(Grpatvar.pm, resvar.levels %in% obs_lev, slevel)
        outtab <- as.table(cbind(outtab, Group=Grpatvar.letterM))
        Grpatvar.letter <- Grpatvar.letterM[resvar.levels %in% obs_lev]

	      # letters labelled by group and level, so the mean table can match them
	      if (length(Grpatvar.letter))
	        names(Grpatvar.letter) <- paste(atvar.levels[i], names(Grpatvar.letter), sep="\r")
	      pmlistLetter[[i]] <- Grpatvar.letter
        pmlistTab[[i]] <- outtab
      }

      if (nrow.pm > 2){
        p_valueMatrix <- stats::setNames(pmlist, atvar.levels)
      }

  pm_existing(c("resvar", "atvar.levels", "rcnplotm", "bkmt", "listlength", "nrow.pm", "pmlist",
                "pmlistTab", "pmlistLetter", "p_valueMatrix"), environment())
}

## The mean table: confidence limits (the comparison intervals with
## letter_ci), the response scale, and letters.
pm_tables <- function(opt, est, cmp) {
  perm <- pm_perm_given(opt$permlist)
  letter_ci <- if (is.null(cmp$letter_ci)) opt$letter_ci else cmp$letter_ci    # (pm_compare())
  cmp_shown <- if (letter_ci) cmp$cmp_int           # the limits the letters are read from
  mt_res <- pm_mean_table(cmp$mt, opt$vars, perm = perm, slevel = opt$slevel, ndecimal = opt$ndecimal,
                          cmp = cmp_shown, ci_df = cmp$ci_df)
  ci_level <- mt_res$ci_level     # level of the table's limits
  bkmt <- pm_back_transform(cmp$bkmt, opt$trans, perm = perm, slevel = opt$slevel,
                            ndecimal = opt$ndecimal, count = opt$count, cmp = cmp_shown, ci_df = cmp$ci_df)
  meanTable <- pm_final_mean_table(mt_res$table, bkmt, opt$vars, cmp$atvar, cmp$resvar, opt$pairwise, opt$meandecr,
                                   letter_ci, opt$slevel, opt$ndecimal, cmp$t.p.valuemGrp, cmp$pmlistLetter,
                                   ci_level, cmp$cmp_int)
  list(meanTable = meanTable, ci_level = ci_level, bkmt = bkmt,
       nc = ncol(bkmt))           # the plots use bkmt[, nc - 1] and bkmt[, nc] (LL, UL)
}

## The results of the stages in one list: a later stage's value replaces an
## earlier one (also by NULL).
pm_flat <- function(...) {
  all <- do.call(c, rev(list(...)))
  all[!duplicated(names(all))]
}

## Mean table with confidence limits (pure function).
## mt: data frame with columns pm, ses, Df and the factors 'vars'.
## Limits: mean -+ t SE at 'slevel', with the df of each mean; with
## permutation p-values (perm) the df are 'ci_df' (pm_perm_ci_df()) and the
## printed Df stay NA. 'cmp' (pm_comparison_intervals()): the limits are the
## comparison intervals, mean +/- cmp$mult SE, instead.
## Returns list(table, ci_level: the level of the limits).
pm_mean_table <- function(mt, vars, perm, slevel, ndecimal, cmp = NULL, ci_df = NULL) {
  meanTable <- mt
  if (!is.null(cmp)) {
    meanTable$LL <- meanTable$pm - cmp$mult * meanTable$ses
    meanTable$UL <- meanTable$pm + cmp$mult * meanTable$ses
    slevel <- round(1 - cmp$conf_level, 3)           # (local: the label of these limits)
    if (perm) meanTable$Df <- NA
  } else {
    crit <- pm_ci_crit(slevel, if (perm) ci_df else meanTable$Df)
    meanTable$LL <- meanTable$pm - crit * meanTable$ses
    meanTable$UL <- meanTable$pm + crit * meanTable$ses
    if (perm) meanTable$Df <- NA
  }
  rownames(meanTable) <- NULL
  meanTable <- meanTable[c(vars, "pm", "ses", "Df", "LL", "UL")]
  meanTable[c("pm", "LL", "UL")] <- round(meanTable[c("pm", "LL", "UL")], ndecimal)
  meanTable$ses <- round(meanTable$ses, ndecimal+1)
  colnames(meanTable) <- c(vars, "Mean", "SE", "Df", paste("LL(", (1 - slevel) * 100, "%)", sep = ""),
                           paste("UL(", (1 - slevel) * 100, "%)", sep = ""))
  list(table = meanTable, ci_level = slevel)
}

## Means and confidence limits on the response scale (pure function).
## bkmt has columns pm, ses, Df (absent with permutation p-values) and the
## factors; returns the factors with Mean, LL and UL.
pm_back_transform <- function(bkmt, trans, perm, slevel, ndecimal, count, cmp = NULL, ci_df = NULL) {
  bkmt$Mean <- trans(bkmt$pm)

  if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)){
    bkmt$Mean <- exp(bkmt$pm)
  }

  mult <- if (!is.null(cmp)) cmp$mult else pm_ci_crit(slevel, if (perm) ci_df else bkmt$Df)
  bkmt$LL <- trans(bkmt$pm - mult * bkmt$ses)
  bkmt$UL <- trans(bkmt$pm + mult * bkmt$ses)

  bkmt$pm <- bkmt$ses <- bkmt$Df <- NULL
  nc <- ncol(bkmt)
  bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], ndecimal)

  if (count) {
    bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)
    bkmt[, (nc - 2):nc][bkmt[, (nc - 2):nc] < 0] <- 0
  }

  rownames(bkmt) <- NULL
  bkmt$Trt <- NULL
  bkmt
}

## Final mean table (pure function): ordered by 'atvar', back-transformed
## columns when they differ from the means, and letter groups when pairwise.
## 'resvar', 't.p.valuemGrp' and 'pmlistLetter' are only used (and only need
## to exist) in the branches that need them. Letters at the user's level
## 'slevel'; the back-transformed limits are labelled with 'ci_level'.
pm_final_mean_table <- function(meanTable, bkmt, vars, atvar, resvar, pairwise, meandecr,
                                letter_ci, slevel, ndecimal, t.p.valuemGrp, pmlistLetter, ci_level = slevel,
                                cmp = NULL) {
  no_atvar <- is.null(atvar) || all(atvar%in%c("NULL", ""))
  if (!no_atvar && !length(resvar)) resvar <- setdiff(vars, atvar)   # (pairwise = FALSE: not set by pm_pairs_atvar())
  if (!no_atvar) {
    all_var_names <- colnames(meanTable)
    meanTable <- meanTable[c(atvar, resvar, setdiff(all_var_names, vars))]
    # atvar levels first, then the other factors; meandecr sorts the means
    # within each level of atvar
    keys <- unname(as.list(meanTable[c(atvar, resvar)]))
    if (!is.null(meandecr) && is.logical(meandecr)) {
      keys <- c(unname(as.list(meanTable[atvar])), list(if (meandecr) -meanTable$Mean else meanTable$Mean))
    }
    meanTable <- meanTable[do.call(order, keys), , drop = FALSE]
  }

  if (any(sort(meanTable$Mean) != sort(bkmt$Mean), na.rm = TRUE)) {
    # rows matched by the factor levels: bkmt can come in another order
    # (atvar first, pm_pairs_atvar())
    k <- function(d) do.call(paste, c(lapply(d[vars], as.character), sep = "\r"))
    cols <- (ncol(bkmt)-2):ncol(bkmt)
    # (bkmt can lack factor columns, e.g. a factor constant within clusters: then
    # its rows are in the order of the means)
    bk <- if (all(vars %in% names(bkmt))) bkmt[match(k(meanTable), k(bkmt)), cols] else
      bkmt[match(rownames(meanTable), rownames(bkmt)), cols]
    meanTable <- cbind(meanTable, round(bk, ndecimal))
    colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("Bk_Mean", paste("Bk_LL(", (1 - ci_level) * 100, "%)", sep = ""),
                                                                   paste("Bk_UL(", (1 - ci_level) * 100, "%)", sep = ""))
  }

  if (pairwise) {
    if (no_atvar) {
      if (!is.null(meandecr) && is.logical(meandecr)) {
        meanTable <- meanTable[order(meanTable$Mean, decreasing = meandecr), ]
      }

      if (letter_ci) {
        meanTable$LetterGrp <- pm_ci_letters(meanTable)
      } else {
        # letters only among means whose comparisons can all be estimated
        meanTable$LetterGrp <- unname(pm_letter_groups(t.p.valuemGrp, rep(TRUE, nrow(t.p.valuemGrp)), slevel))
      }

      meanTable <- list(Table=meanTable, Note=pm_letter_note(letter_ci, slevel, cmp))
    } else {
      if (letter_ci) {
        # CI overlap within each 'atvar' group
        grp <- do.call("paste", c(lapply(meanTable[atvar], as.character), sep="\r"))
        lets <- character(nrow(meanTable))
        for (g in unique(grp)) lets[grp == g] <- pm_ci_letters(meanTable[grp == g, , drop = FALSE])
        meanTable$LetterGrp <- lets
      } else {
        lets <- unlist(pmlistLetter)
        mkey <- paste(do.call("paste", c(lapply(meanTable[atvar], as.character), sep=" : ")),
                      do.call("paste", c(lapply(meanTable[resvar], as.character), sep=" : ")), sep="\r")
        meanTable$LetterGrp <- unname(lets[mkey])
      }
      meanTable <- list(Table=meanTable, Note=paste(pm_letter_note(letter_ci, slevel, cmp), " at each level of ", paste(atvar, collapse = ":"), sep=""))
    }
  }
  meanTable
}

## Letters from CI overlap (pm_ci_mcp()) for the rows of a mean table; means
## without a CI (not estimable) get no letter.
pm_ci_letters <- function(tab) {
  LL <- tab[[grep("^LL", names(tab))[1]]]
  UL <- tab[[grep("^UL", names(tab))[1]]]
  ok <- is.finite(LL) & is.finite(UL)
  if (all(ok)) return(pm_ci_mcp(LL, UL))
  out <- rep("", length(LL))
  if (any(ok)) out[ok] <- unname(pm_ci_mcp(LL[ok], UL[ok]))
  out
}

## Note under the mean table saying how the letters were made. With letter_ci
## they are read from the overlap of the comparison intervals 'cmp'
## (R/comparison.R), which the table shows as LL and UL.
pm_letter_note <- function(letter_ci, slevel, cmp = NULL) {
  if (!letter_ci) return(paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep=""))
  if (is.null(cmp)) {
    return(paste("Letter-based representation of pairwise comparisons by overlap of the confidence limits", sep=""))
  }
  paste0("Letter-based representation of pairwise comparisons by overlap of the comparison intervals (LL, UL: mean +/- ",
         signif(cmp$mult, 4), " SE, ", round(100*cmp$conf_level, 1), "% intervals): two means differ at significant level '",
         slevel, "' when their intervals do not overlap", pm_comparison_caveat(cmp))
}

## Stable tables and the returned list, from the results of the stages
## ('v', pm_flat()).
pm_assemble <- function(v) {
  # Stable result tables (means, pairs, settings).
  tidy_parts <- pm_tidy_parts(v)
  x <- stats::setNames(lapply(pm_legacy_inputs, function(n) v[[n]]), pm_legacy_inputs)
  out <- pm_add_advice(pm_finish(pm_legacy_output(x), tidy_parts), pm_means_advice(v, tidy_parts))
  # for autoplot(data = TRUE): only what the plot needs, the factors of
  # modelterm and the response on both scales (not the model: a saved result
  # would carry the whole fit)
  attr(out, "plot_data") <- pm_plot_rows(v$model, v$vars, if (identical(tidy_parts$settings$response$back, "trans")) v$trans)
  out
}

## The observations for plots: the factors of 'vars' (as character, to keep
## the result small and level order to the means), .y and .y_orig; NULL when
## the data cannot be found (autoplot() then asks for data = ).
pm_plot_rows <- function(model, vars, trans = NULL) {
  m <- if (inherits(model, "pm_aovlist")) attr(model, "aovlist") else model
  d <- tryCatch(suppressMessages(suppressWarnings(pm_data(m, trans = trans))), error = function(e) NULL)
  # (an adapter's own pm_data() method may return the data only: no plot rows)
  if (!is.data.frame(d) || !all(c(vars, ".y", ".y_orig") %in% names(d))) return(NULL)
  out <- d[c(vars, ".y", ".y_orig")]
  for (v in vars) if (is.factor(out[[v]])) out[[v]] <- as.character(out[[v]])
  if (isTRUE(all.equal(out$.y, out$.y_orig))) out$.y_orig <- NULL
  rownames(out) <- NULL
  out
}

## Values the legacy result list is built from.
pm_legacy_inputs <- c("pairwise", "atvar", "permlist", "adj", "ndecimal", "mean.table", "se.table",
                      "SED.out", "LSD", "LSDm", "t.p.valuem",
                      "meanTable", "p_valueMatrix", "listlength", "resvar", "pmlistTab", "atvar.levels", "nsim")

## The legacy result list of predictmeans() (pure function); 'x' is a named
## list of the values in pm_legacy_inputs.
pm_legacy_output <- function(x) {
  pairwise <- x$pairwise; atvar <- x$atvar; permlist <- x$permlist; adj <- x$adj
  ndecimal <- x$ndecimal; mean.table <- x$mean.table; se.table <- x$se.table
  SED.out <- x$SED.out; LSD <- x$LSD; LSDm <- x$LSDm; t.p.valuem <- x$t.p.valuem
  meanTable <- x$meanTable; p_valueMatrix <- x$p_valueMatrix; listlength <- x$listlength
  resvar <- x$resvar; pmlistTab <- x$pmlistTab; atvar.levels <- x$atvar.levels; nsim <- x$nsim

  if (pairwise) {

    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {

      # (the matrix of pairwise LSDs of 1.x is not printed: each pair's LSD is
      # in pairs$lsd, and the LSD summary goes with the p-values)
      return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                  "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4),
                  mean_table = meanTable, p_valueMatrix = p_valueMatrix))
    } else {

      outputlist <- vector("list", listlength+7)
      outputlist[[1]] <- mean.table
      outputlist[[2]] <- se.table
      outputlist[[3]] <- SED.out
      outputlist[[4]] <- LSD
      outputlist[[5]] <- paste("For variable '", paste(resvar, collapse =" and "), "' at each level of '", paste(atvar, collapse =" and "), "'", sep="")

      ## DONGWEN: I think this could be replaced with
      ## outputlist[[6:(listlength + 5)]] <- pmListTab[[1:listLength]]
      for (i in 6: (listlength+5)) {
        outputlist[[i]] <- pmlistTab[[i-5]]
      }

      outputlist[[listlength+6]] <- meanTable
      outputlist[[listlength+7]] <- p_valueMatrix

      # print(outputlist)
      if (!pm_perm_given(permlist)) {
        names(outputlist)<- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                              "LSD", paste("Pairwise comparison p-value (adjusted by '", adj, "' method)", sep=""),
                              atvar.levels, "mean_table", "p_valueMatrix")[seq_along(outputlist)]
      } else {
        names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                 "Approximated LSD"), paste("Pairwise '", nsim, "' times permuted p-value (adjusted by '", adj, "' method)", sep=""),
                               atvar.levels, "mean_table", "p_valueMatrix")[seq_along(outputlist)]
      }

            return(outputlist)
    }
  } else {
    outputlist=list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD,
                    mean_table=meanTable)
        return(outputlist)
  }
}


## Helpers for the compare stage (pure functions).

## Pairs (i, j) with i > j, ordered by i then j: (2,1), (3,1), (3,2), ...
pm_pair_index <- function(n) {
  if (n < 2) return(list(i = integer(0), j = integer(0)))
  list(i = rep(2:n, times = 1:(n - 1)), j = sequence(1:(n - 1)))
}

## Contrast matrix of the pairs: row k is mean j minus mean i.
pm_pair_matrix <- function(n) {
  pp <- pm_pair_index(n)
  CM <- matrix(0, nrow = length(pp$i), ncol = n)
  CM[cbind(seq_along(pp$i), pp$j)] <- 1
  CM[cbind(seq_along(pp$i), pp$i)] <- -1
  CM
}

## Level labels of the first 'nvar' columns of a label table, as a character
## matrix without dimnames.
pm_label_matrix <- function(lab, nvar) {
  m <- vapply(lab[seq_len(nvar)], as.character, character(nrow(lab)))
  m <- matrix(m, nrow = nrow(lab))
  dimnames(m) <- NULL
  m
}


## Letter groups from a symmetric matrix of p-values (dimnames = levels),
## with NA for comparisons that are not available (no data, or not
## estimable). Levels without data get NA. Among the others, the level with
## most unavailable comparisons is dropped, one at a time, until all the
## remaining comparisons are available; those levels get letters. A level
## left alone gets "A" when it is the only level with data (nothing to
## compare), and NA when it was left alone because its comparisons are not
## estimable.
pm_letter_groups <- function(P, has_data, slevel) {
  levs <- rownames(P)
  out <- rep(NA_character_, length(levs)); names(out) <- levs
  ok <- has_data
  # The letter algorithm (insert-absorb) can take minutes for a few hundred
  # means, and that many letters cannot be read: no letters above the limit.
  max_n <- pm_max_means()
  if (sum(ok) > max_n) {
    message("No letter groups for ", sum(ok), " means (more than ", max_n, "): they could not be read and take ",
            "long to make. The comparisons are in the 'pairs' table; options(predictmeans.max_means = ) ",
            "changes the limit.")
    return(out)
  }
  diag(P) <- 1
  while (sum(ok) > 1 && anyNA(P[ok, ok])) {
    n_na <- rowSums(is.na(P[ok, ok, drop = FALSE]))
    ok[which(ok)[which.max(n_na)]] <- FALSE
  }
  if (sum(ok) >= 2) {
    out[ok] <- multcompLetters(P[ok, ok, drop = FALSE], Letters = LETTERS, threshold = slevel)[levs[ok]]
  } else if (sum(ok) == 1 && sum(has_data) == 1) {
    out[ok] <- LETTERS[1]
  }
  out
}
