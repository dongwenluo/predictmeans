## The stages of predictmeans(), in this order. Each is an ordinary function
## that takes the results of the stages before it and returns a named list:
##   opt  pm_setup()     the arguments as used: checked and resolved options,
##                       link or response back-transformation, plot choices
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
## two readers of everything: pm_make_plots() (only when plot = TRUE; all
## statistics are computed before any plotting) and pm_assemble() (the stable
## tables of R/tidy.R and the legacy list).

## The plot settings of predictmeans(): its argument 'plot_args', a named list
## of these (the defaults are those of 1.x).
pm_plot_defaults <- list(lsd_bar = TRUE, plotord = NULL, lineplot = TRUE, plottitle = NULL,
                         plotxlab = NULL, plotylab = NULL, jitterv = 0.2, basesz = 12)

## The plot settings of a call ('plot_args'; names may be abbreviated).
## Returns the settings with the defaults for those not given.
pm_plot_args <- function(plot_args = list()) {
  nms <- names(pm_plot_defaults)
  if (!is.list(plot_args) || (length(plot_args) && (is.null(names(plot_args)) || !all(nzchar(names(plot_args)))))) {
    stop("'plot_args' must be a named list, e.g. plot_args = list(plottitle = \"Yield\", basesz = 14).", call. = FALSE)
  }
  k <- pmatch(names(plot_args), nms, duplicates.ok = TRUE)
  if (anyNA(k)) stop("'plot_args' can hold ", paste(nms, collapse = ", "), "; not ",
                     paste(names(plot_args)[is.na(k)], collapse = ", "), call. = FALSE)
  out <- pm_plot_defaults
  for (i in seq_along(k)) out[nms[k[i]]] <- list(plot_args[[i]])
  out
}

## Which plots: 'plot' as a vector of plot kinds ("mean", "bar", "pval",
## "back"); the flags the plot code reads.
pm_which_plots <- function(which.plots) {
  kinds <- c("mean", "bar", "pval", "back")
  which.plots <- match.arg(which.plots, choices = kinds, several.ok = TRUE)
  stats::setNames(as.list(kinds %in% which.plots), c("mplot", "barplot", "pplot", "bkplot"))
}

## Arguments, options, link function and plot choices: the arguments of
## predictmeans() ('a') and its plot settings ('pa', pm_plot_args()) as they
## are used by the later stages.
pm_setup <- function(a, pa) {
  model <- a$model
  modelterm <- a$modelterm
  pm_check_modelterm(model, modelterm)
  pm_note_not_a_term(model, modelterm)
  adj <- pm_match_adj(a$adj, westfall = TRUE)
  wp <- pm_which_plots(a$which.plots)

  trans <- if (is.null(a$trans)) I else a$trans
  if (inherits(model, "aovlist")){
    model <- pm_aov_carrier(model)
  }

  # Models with a link function other than the identity: means are
  # back-transformed with its inverse (pm_linkinv(); see ?pm_coef to support
  # other model classes). Otherwise the user's 'trans' is used.
  linkinv <- pm_link_trans(model)
  link_scale <- !is.null(linkinv)
  trans_by <- if (identical(trans, I)) NA_character_ else "trans"
  if (link_scale) {
    trans <- linkinv
    trans_by <- "the inverse link"
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

  # comparisons are made when asked for in any of these ways
  atvar <- if (length(vars) == 1) NULL else a$atvar
  pairwise <- a$pairwise
  permlist <- a$permlist
  if (pm_perm_given(permlist)) {
    pairwise <- TRUE
    # "tukey", "mvt": single-step max-|t|; "westfall": step-down (pm_perm_maxT())
    if (adj == "scheffe") {
      stop("adj = \"scheffe\" is not available for permutation p-values; use another 'adj'.", call. = FALSE)
    }
  } else if (adj == "westfall") {
    stop("adj = \"westfall\" (step-down max-|t|) needs permutation p-values ('permlist'); without them use ",
         "\"mvt\" (single step) or \"holm\".", call. = FALSE)
  }
  if (!is.null(atvar) && !unique(atvar%in%c("NULL", ""))){
    pairwise <- TRUE
  }
  if (adj != "none"){
    pairwise <- TRUE
  }
  # letters by overlap of the comparison intervals (R/comparison.R); 'level'
  # is the significance level of the comparisons, as everywhere
  if (a$letterCI) {
    pairwise <- TRUE
    adj <- "none"
  }

  # 'level' and 'slevel' are the user's level everywhere; the printed LSDs
  # take the critical value of 'adj' at this level (pm_lsd_crit()) and the
  # mean table with a permlist (limits +- 2 SE) 'ci_level' (pm_tables())
  level <- a$level

  # an order for the plot asks for the mean plot
  plot <- a$plot
  plotord <- pa$plotord
  if (!is.null(plotord) && !unique(plotord%in%c("NULL", ""))) {
    plot <- wp$mplot <- TRUE
    if(length(plotord) != length(vars)){
      stop(paste("plotord must be a vector of length", length(vars)))
    }
    if(is.character(plotord)){
      i <- match(plotord, vars)
      if(any(is.na(i))){
        stop(paste("The entries", paste0(plotord[is.na(i)], sep = ", "), "do not match any term in ", modelterm))
      }
      plotord <- i
    }
  }

  # zero-inflated glmmTMB fits: the overall means (1 - pi) mu (R/zero-inflation.R)
  mean_type <- if (is.null(a$mean_type)) "conditional" else a$mean_type
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

  list(model = model, modelterm = modelterm, data = a$data, vars = vars, mdf = mdf, mean_type = mean_type,
       pairwise = pairwise, atvar = atvar, adj = adj, Df = a$Df, level = level, slevel = level,
       covariate = a$covariate, cov_reduce = a$cov_reduce, exclude = a$exclude, weights = a$weights,
       meandecr = if (is.logical(a$meandecr)) a$meandecr else NULL, letterCI = a$letterCI,
       trans = trans, trans_by = trans_by, link_scale = link_scale, count = count,
       bias_adjust = a$bias_adjust, sigma = a$sigma, average = a$average,
       permlist = permlist, ncore = a$ncore, ndecimal = 4L,
       df_method = a$df_method, vcov = a$vcov, vcov_adjust = a$vcov_adjust,
       vcov_adjust_default = a$vcov_adjust_default, equivalence = a$equivalence,
       plot = plot, mplot = wp$mplot, barplot = wp$barplot, pplot = wp$pplot, bkplot = wp$bkplot,
       plotord = plotord, lsd_bar = pa$lsd_bar, lineplot = pa$lineplot,
       plottitle = pa$plottitle, plotxlab = pa$plotxlab, plotylab = pa$plotylab, jitterv = pa$jitterv,
       basesz = pa$basesz)
}

## K matrix, predicted means, standard errors and the legacy mean/SE tables.
pm_estimate <- function(opt) {
  # (the grid of all combinations is only needed to average on the response scale)
  ctr.matrix <- Kmatrix(opt$model, opt$modelterm, opt$covariate, data=opt$data,
                        cov_reduce=opt$cov_reduce, exclude=opt$exclude, weights=opt$weights,
                        refgrid=identical(opt$average, "response") || identical(opt$mean_type, "overall"))
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

  K <- KK[, mp$estimable, drop = FALSE]          # To match coef names
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
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
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
  pm_df_raw <- Df_diff <- pairDf <- df_user <- pair_tab <- nsim <- NULL
  t.p.valuem <- t.p.valuep <- LSDm <- t.p.valuemGrp <- p_valueMatrix <- NULL
  resvar <- atvar.levels <- rcnplotm <- listlength <- nrow.pm <- pmlist <- pmlistTab <- pmlistLetter <- NULL
  # means over empty cells (R/empty-cells.R): their pairs are compared on the
  # cells both means have, so nothing read off the means themselves
  # (comparison intervals, letters from them) can agree with their tests
  ec <- est$ctr.matrix$empty_cells
  affected <- stats::setNames(if (is.null(ec)) rep(FALSE, nK) else unname(ec$affected[rnK] %in% TRUE), rnK)
  letterCI <- opt$letterCI
  result <- function() {
    list(nK = nK, rnK = rnK, perm = perm, atvar = atvar, mt = mt, bkmt = bkmt, Df = Df, SED.out = SED.out, LSD = LSD,
         cmp_int = cmp_int, rK = rK, cm = cm, dses = dses, dses.df = dses.df, dses.m = dses.m,
         KKvarndiff = KKvarndiff, labKK = labKK, pair_over = pair_over, pm_df_raw = pm_df_raw, Df_diff = Df_diff,
         pairDf = pairDf, df_user = df_user, pair_tab = pair_tab, nsim = nsim, t.p.valuem = t.p.valuem,
         t.p.valuep = t.p.valuep, LSDm = LSDm, t.p.valuemGrp = t.p.valuemGrp, p_valueMatrix = p_valueMatrix,
         resvar = resvar, atvar.levels = atvar.levels, rcnplotm = rcnplotm, listlength = listlength,
         nrow.pm = nrow.pm, pmlist = pmlist, pmlistTab = pmlistTab, pmlistLetter = pmlistLetter,
         letterCI = letterCI, affected = affected)
  }
  if (nK == 1) return(result())

  # The pairs are tabulated (one row per pair) when comparisons are asked
  # for (pairwise, adj, atvar, letterCI, permlist, equivalence), and
  # otherwise for up to pm_max_means() means. Above that (several hundred
  # lines of a breeding trial) the SEDs and LSDs are summarised from the
  # covariance of the means, without a row per pair (pm_sed_summary()).
  # Means over empty cells and non-estimable means need each pair's
  # contrast, so their SEDs take the first route too.
  tab_pairs <- opt$pairwise || !is.null(opt$equivalence) || nK <= pm_max_means()
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
  atvar <- s[["atvar"]]              # NULL when a factor has no pairs within its levels (nested factors)

  if (!perm) {
    d <- pm_compare_df(est$inf, opt$modelterm, K, rK, mt, Df)
    Df <- d[["Df"]]; mt <- d[["mt"]]; pm_df_raw <- d[["pm_df_raw"]]; Df_diff <- d[["Df_diff"]]
    pairDf <- d[["pairDf"]]; df_user <- d[["df_user"]]
    bkmt <- mt
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
  if (tab_pairs) cmp_int <- pm_comparison_intervals(pair_tab, mt$pm, mt$ses, if (perm) NULL else mt$Df, opt$level,
                                                    opt$adj, affected)
  if (letterCI && any(misread)) {
    message("letterCI: ", paste(rnK[misread], collapse = ", "), " average", if (sum(misread) == 1) "s",
            " over empty cells and ", if (sum(misread) == 1) "is" else "are", " compared on common cells, which no ",
            "interval around a mean can show; the letters are from the p-values.")
    letterCI <- FALSE
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
    scale_r <- pmax(1, apply(abs(rK), 1, max))
    non_est_pair <- apply(abs(pd), 1, max) > 1e-6 * scale_r
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
        p_valueMatrix <- pmlist
      }

  pm_existing(c("resvar", "atvar.levels", "rcnplotm", "bkmt", "listlength", "nrow.pm", "pmlist",
                "pmlistTab", "pmlistLetter", "p_valueMatrix"), environment())
}

## The mean table: confidence limits (the comparison intervals with
## letterCI), the response scale, and letters.
pm_tables <- function(opt, est, cmp) {
  perm <- pm_perm_given(opt$permlist)
  letterCI <- if (is.null(cmp$letterCI)) opt$letterCI else cmp$letterCI    # (pm_compare())
  cmp_shown <- if (letterCI) cmp$cmp_int           # the limits the letters are read from
  mt_res <- pm_mean_table(cmp$mt, opt$vars, perm = perm, slevel = opt$slevel, ndecimal = opt$ndecimal,
                          cmp = cmp_shown)
  ci_level <- mt_res$ci_level     # level of the table's limits (0.05 with a permlist: +- 2 SE)
  bkmt <- pm_back_transform(cmp$bkmt, opt$trans, perm = perm, slevel = opt$slevel,
                            ndecimal = opt$ndecimal, count = opt$count, cmp = cmp_shown)
  meanTable <- pm_final_mean_table(mt_res$table, bkmt, opt$vars, cmp$atvar, cmp$resvar, opt$pairwise, opt$meandecr,
                                   letterCI, opt$slevel, opt$ndecimal, cmp$t.p.valuemGrp, cmp$pmlistLetter,
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
## With permutation p-values there are no df: limits are mean +/- 2 SE, at
## about the 0.05 level. 'cmp' (pm_comparison_intervals()): the limits are the
## comparison intervals, mean +/- cmp$mult SE, instead.
## Returns list(table, ci_level: the level of the limits).
pm_mean_table <- function(mt, vars, perm, slevel, ndecimal, cmp = NULL) {
  meanTable <- mt
  if (!is.null(cmp)) {
    meanTable$LL <- meanTable$pm - cmp$mult * meanTable$ses
    meanTable$UL <- meanTable$pm + cmp$mult * meanTable$ses
    slevel <- round(1 - cmp$conf.level, 3)           # (local: the label of these limits)
    if (perm) meanTable$Df <- NA
  } else if (!perm) {
    meanTable$LL <- meanTable$pm - qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
    meanTable$UL <- meanTable$pm + qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
  } else {
    meanTable$LL <- meanTable$pm - 2 * meanTable$ses
    meanTable$UL <- meanTable$pm + 2 * meanTable$ses
    slevel <- 0.05                     # (local: the label of these limits)
    meanTable$Df <- NA
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
pm_back_transform <- function(bkmt, trans, perm, slevel, ndecimal, count, cmp = NULL) {
  bkmt$Mean <- trans(bkmt$pm)

  if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)){
    bkmt$Mean <- exp(bkmt$pm)
  }

  if (!is.null(cmp)) {                       # comparison intervals (letterCI)
    bkmt$LL <- trans(bkmt$pm - cmp$mult * bkmt$ses)
    bkmt$UL <- trans(bkmt$pm + cmp$mult * bkmt$ses)
  } else if (!perm) {
    bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)
    bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)
  } else {
    bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)
    bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)
  }

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
                                letterCI, slevel, ndecimal, t.p.valuemGrp, pmlistLetter, ci_level = slevel,
                                cmp = NULL) {
  no_atvar <- is.null(atvar) || all(atvar%in%c("NULL", ""))
  if (!no_atvar) {
    all_var_names <- colnames(meanTable)
    meanTable <- meanTable[c(atvar, resvar, setdiff(all_var_names, vars))]
    meanTable <- with(meanTable, meanTable[order(eval(parse(text=c(atvar, resvar)[length(vars):1]))),])
  }

  if (any(meanTable$Mean!=bkmt$Mean, na.rm = TRUE)) {
    meanTable <- cbind(meanTable, round(bkmt[, (ncol(bkmt)-2):ncol(bkmt)], ndecimal))
    colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("Bk_Mean", paste("Bk_LL(", (1 - ci_level) * 100, "%)", sep = ""),
                                                                   paste("Bk_UL(", (1 - ci_level) * 100, "%)", sep = ""))
  }

  if (pairwise) {
    if (no_atvar) {
      if (!is.null(meandecr) && is.logical(meandecr)) {
        meanTable <- meanTable[order(meanTable$Mean, decreasing = meandecr), ]
      }

      if (letterCI) {
        meanTable$LetterGrp <- pm_ci_letters(meanTable)
      } else {
        # letters only among means whose comparisons can all be estimated
        meanTable$LetterGrp <- unname(pm_letter_groups(t.p.valuemGrp, rep(TRUE, nrow(t.p.valuemGrp)), slevel))
      }

      meanTable <- list(Table=meanTable, Note=pm_letter_note(letterCI, slevel, cmp))
    } else {
      if (letterCI) {
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
      meanTable <- list(Table=meanTable, Note=paste(pm_letter_note(letterCI, slevel, cmp), " at each level of ", paste(atvar, collapse = ":"), sep=""))
    }
  }
  meanTable
}

## Letters from CI overlap (ci_mcp()) for the rows of a mean table; means
## without a CI (not estimable) get no letter.
pm_ci_letters <- function(tab) {
  LL <- tab[[grep("^LL", names(tab))[1]]]
  UL <- tab[[grep("^UL", names(tab))[1]]]
  ok <- is.finite(LL) & is.finite(UL)
  if (all(ok)) return(ci_mcp(LL, UL))
  out <- rep("", length(LL))
  if (any(ok)) out[ok] <- unname(ci_mcp(LL[ok], UL[ok]))
  out
}

## Note under the mean table saying how the letters were made. With letterCI
## they are read from the overlap of the comparison intervals 'cmp'
## (R/comparison.R), which the table shows as LL and UL.
pm_letter_note <- function(letterCI, slevel, cmp = NULL) {
  if (!letterCI) return(paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep=""))
  if (is.null(cmp)) {
    return(paste("Letter-based representation of pairwise comparisons by overlap of the confidence limits", sep=""))
  }
  paste0("Letter-based representation of pairwise comparisons by overlap of the comparison intervals (LL, UL: mean +/- ",
         signif(cmp$mult, 4), " SE, ", round(100*cmp$conf.level, 1), "% intervals): two means differ at significant level '",
         slevel, "' when their intervals do not overlap", pm_comparison_caveat(cmp))
}

## All graphics of predictmeans(): p-value plots, mean and bar plots,
## back-transformed plot. Draws the requested plots and returns
## list(meanPlot, ciPlot, predictmeansBarPlot). 'x': the results of the stages
## (pm_flat()); these values of it are used.
pm_make_plots <- function(x) {
  LSD <- x[["LSD"]]; SED.out <- x[["SED.out"]]; adj <- x[["adj"]]
  atvar <- x[["atvar"]]; barplot <- x[["barplot"]]; basesz <- x[["basesz"]]
  ci_level <- x[["ci_level"]]; bkmt <- x[["bkmt"]]; bkplot <- x[["bkplot"]]
  jitterv <- x[["jitterv"]]; lineplot <- x[["lineplot"]]; lsd_bar <- x[["lsd_bar"]]
  mdf <- x[["mdf"]]; model <- x[["model"]]; modelterm <- x[["modelterm"]]
  mplot <- x[["mplot"]]; mt <- x[["mt"]]; nK <- x[["nK"]]
  nc <- x[["nc"]]; nrow.pm <- x[["nrow.pm"]]
  pairwise <- x[["pairwise"]]; plot <- x[["plot"]]; plotord <- x[["plotord"]]
  plottitle <- x[["plottitle"]]; plotxlab <- x[["plotxlab"]]; plotylab <- x[["plotylab"]]
  pm <- x[["pm"]]; pmlist <- x[["pmlist"]]; pplot <- x[["pplot"]]
  rcnplotm <- x[["rcnplotm"]]; response <- x[["response"]]
  resvar <- x[["resvar"]]; ses <- x[["ses"]]
  slevel <- x[["slevel"]]; t.p.valuep <- x[["t.p.valuep"]]; trans <- x[["trans"]]
  link_scale <- x[["link_scale"]]; vars <- x[["vars"]]
  LL <- UL <- Mean <- Trt <- NULL                       # (columns named in the plots)
  meanPlot <- ciPlot <- predictmeansBarPlot <- NULL

  # p-value level plots
  if (pairwise && nK > 1) {
    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
      if (all(nrow(t.p.valuep) > 2, pplot, plot)) {
        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Level Plot of p-value (adjusted by '", adj, "' method)\n for Pairwise Comparison", sep="")
        }

        print(PMplot(t(t.p.valuep), level=slevel, mtitle=mtitle, basesz=basesz))
      }
    } else {
      if (all(nrow.pm > 2, pplot, plot)) {
        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Adjusted p-value (by '", adj,
                          "' method)\n for Pairwise Comparison at Each Level of '",paste(atvar, collapse =" and "), "'", sep="")
        }

        print(PMplot(pmlist, level=slevel, xylabel=rcnplotm, mtitle=mtitle, basesz=basesz))
      }
    }
  }

  # means and bar plots (none when no mean can be estimated)
  if (plot && all(is.na(mt$pm))) message("No mean can be estimated: no mean plots.")
  if (plot && !all(is.na(mt$pm))) {

    if (length(vars) > 3){
      cat("\n", "There is no plot for more than three-way interaction! \n\n")
    }

    plotmt <- na.omit(mt)
    yMin <- min(plotmt[, "pm"])
    yMax <- max(plotmt[, "pm"])
    offSet <- 0.25 * (yMax - yMin)

    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
      if (lsd_bar){
        LSD_value <- LSD[3]
      }else{
        LSD_value <- SED.out[3]
      }
    }else{
      if (lsd_bar){
        LSD_value <- mean(attr(LSD,"For the Same Level of Factor")[1, atvar])
      }else{
        LSD_value <- mean(attr(SED.out, "For the Same Level of Factor")[1, atvar])
      }
    }

    if (lsd_bar){
      bar_label <- "Aveg.LSD"
    }else{
      bar_label <- "Aveg.SED"
    }

    up <- yMin + LSD_value
    lsdBar <- cbind(plotmt[, vars, drop=FALSE], up=up, yMin=yMin)
    limits <- aes(ymax = (pm + ses)*(pm > 0) + pmin(pm + ses, 0)*(pm <= 0), ymin=(pm - ses)*(pm < 0) + pmax(pm - ses, 0)*(pm >= 0))

    if (length(vars) == 1) {
      mxlab <- plotxlab

      if (is.null(plotxlab) || plotxlab%in%c("NULL", "")){
        mxlab <- paste("\n", vars, sep="")
      }

      mylab <- plotylab

      if (is.null(plotylab) || plotylab%in%c("NULL", "")){
        mylab <- paste(response, "\n", sep="")
      }

      if (mplot) {

        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for \"", vars, "\" with ", bar_label, " (", ci_level * 100, "%) Bar", sep="")
        }

        p1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm, group=1))+
            labs(title=paste(pm_wrap_title(mtitle, basesz), "\n", sep=""), x=mxlab, y=mylab)+
            lims(x= c(bar_label, levels(plotmt[, vars])), y = c(yMin - offSet, max(yMax + offSet, yMin + LSD_value + offSet))) +
            geom_point(colour="red", size=2)+
          #  geom_line(linewidth=0.5)+
            geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue") +  #data=lsdBar,
            theme_bw(basesz)

  		  if (lineplot){
  		    p1 <- p1 + geom_line(linewidth=0.5)
  		  }

        meanPlot <- p1

        print(p1)
      } # end of if mplot

      if (barplot) {


        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
        }

        dodge <- position_dodge(width=0.9)
        bp1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm))+
            labs(title=paste(pm_wrap_title(mtitle, basesz), "\n", sep=""), x=mxlab, y=mylab)+
            ylim(c(min(0, (min(pm) - max(ses)) * 1.2), max(0, (max(pm) + max(ses)) * 1.2))) +
            geom_bar(position=dodge, stat="identity", fill="papayawhip", colour="darkgreen") +
            geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())

        predictmeansBarPlot <- bp1
        print(bp1)
      }
    } else if (length(vars) == 2) {
      ##NB: I have changed this to an else statement because there shouldn't be any other way to get here. That is, length(vars) will not
      ##change from 1 to 2 between the block where length(vars)==1 is true and here, therefore the cases should be else if statements.
        if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
          plotord <- 1:2
          if (!(is.null(atvar) || all(atvar%in%c("NULL", "")))) {
            atvar_num <- which(vars %in% atvar)
            plotord <- c(atvar_num, setdiff(seq_along(vars), atvar_num))
          }
        }

        fact1 <- (vars[plotord])[1]
        fact2 <- (vars[plotord])[2]

        mxlab <- plotxlab

        if (is.null(plotxlab) || plotxlab%in%c("NULL", "")) {
          mxlab <- paste("\n", fact1, sep="")
        }

        mylab <- plotylab

        if (is.null(plotylab) || plotylab%in%c("NULL", "")) {
          mylab <- paste(response, "\n", sep="")
        }

        if (mplot) {


          mtitle <- plottitle

          if (is.null(plottitle) || plottitle%in%c("NULL", "")){
            mtitle <- paste("Predicted means for \"", fact1, "\" by \"", fact2, "\" with ",
                            paste(atvar, collapse=" "), " ", bar_label, " (", ci_level * 100, "%) Bar", sep = "")
          }

          plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))
          p2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), col=eval(parse(text = fact2))))+
            labs(title=paste(pm_wrap_title(mtitle, basesz), "\n", sep=""), x=mxlab, y=mylab)+
            lims(x= levels(plotmt[, fact1]), y = c(yMin - offSet, max(yMax + offSet, yMin + LSD_value + offSet))) +
			      geom_point(size=2)+
           # geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), linewidth=0.96)+
            geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
          #  guides(linetype = guide_legend(title = fact2))+
            guides(col = guide_legend(title = fact2))+
            theme_bw(basesz)

          if (lineplot) {
            p2 <- p2 + geom_line(aes(linetype=eval(parse(text = fact2)),
                                     col=eval(parse(text = fact2))), linewidth=0.96)+
                  guides(linetype = guide_legend(title = fact2))
          }

          meanPlot <- p2

          print(p2)
        } # end if mplot

        if (barplot) {


          mtitle <- plottitle

          if (is.null(plottitle) || plottitle%in%c("NULL", "")) {
            mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
          }

          dodge <- position_dodge(width=0.9)

          bp2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), fill=  eval(parse(text = fact2))))+
            geom_point(position=dodge) +
            geom_bar(stat = "identity", position=dodge)+
            labs( x = mxlab, y = mylab, title = pm_wrap_title(mtitle, basesz))+
            ylim(c(min(0, (min(pm) - max(ses)) * 1.1), max(0, (max(pm) + max(ses)) * 1.1))) +
            geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            scale_fill_brewer()+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
            theme(legend.position = "top")+guides(fill = guide_legend(title = fact2))
          predictmeansBarPlot <- bp2

          print(bp2)
        }
      } else if (all(length(vars) == 3, mplot)) {
        ## NB: I have changed this to else if, because neither length(vars)
        ## or mplot should have changed, and I believe this is just the third
        ## choice

        mtitle <- plottitle

        if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
          plotord <- 1:3

          if (!(is.null(atvar) || all(atvar%in%c("NULL", ""))) && length(atvar)==1){
            atvar_num <- which(vars %in% atvar)
            plotord <- c(atvar_num, setdiff(seq_along(vars), atvar_num))
          }
        }

        fact1 <- (vars[plotord])[1]
        fact2 <- (vars[plotord])[2]
        fact3 <- (vars[plotord])[3]
        plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for '", fact1, "' by '", fact2, "' for each '",
                           fact3, "'\n with ", paste(atvar, collapse=" "), " ",
                          bar_label, " (", ci_level * 100, "%) Bar\n", sep = "")
        }

        mxlab <- plotxlab

        if (is.null(plotxlab) || plotxlab%in%c("NULL", "")){
          mxlab <- paste("\n", fact1, sep="")
        }

        mylab <- plotylab

        if (is.null(plotylab) || plotylab%in%c("NULL", "")) {
          mylab <- paste(response, "\n", sep="")
        }

        p3 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=factor(eval(parse(text = fact2))), col=factor(eval(parse(text = fact2)))))+
          labs(title=paste(pm_wrap_title(mtitle, basesz), "\n", sep=""), x=mxlab, y=mylab)+
          lims(x= levels(plotmt[, fact1]), y = c(yMin-0.5*offSet, max(yMax+0.5*offSet, yMin+LSD_value+0.5*offSet))) +
          geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
		      geom_point(size=2)+
         # geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), linewidth=0.8)+
          facet_grid(eval(parse(text = paste("~",fact3, sep=""))))+
          guides(group = guide_legend(fact2))+
         # guides(linetype = guide_legend(fact2))+
          guides(col = guide_legend(fact2))+
          theme_bw(basesz)

        if (lineplot){
          p3 <- p3 + geom_line(aes(linetype=eval(parse(text = fact2)),
                                   col=eval(parse(text = fact2))), linewidth=0.8) +
            guides(linetype = guide_legend(title = fact2))
        }

        meanPlot <- p3

        print(p3)
      }
    }

  # back-transformed means with the data (works on copies of bkmt and mdf)
  if (!is.null(trans)) local({
    if (plot && bkplot && !all(is.na(mt$pm))) {

      if (!is.null(atvar) && !unique(atvar%in%c("NULL", ""))) {
        mdf <- mdf[do.call(order, mdf[, c(atvar, resvar), drop=FALSE]),]
      }

      if (response %in% names(mdf)) {    ## Transformed y before modelling

        if (inherits(mdf[, response], "factor")) {
          bky <- as.numeric(mdf[, response])-1
        } else {

          if (link_scale) {    # response scale already
            bky <- mdf[, response]

            if (!is.null(dim(mdf[, response]))){
              bky <- mdf[, response][,1]/rowSums(mdf[, response])
            }

          } else {
            bky <- trans(mdf[, response])
          }# end of if glm or glmer
        }# end of if factor
      } else {       ## Transformed y within modelling
        nresponse <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]

        if (!(nresponse %in% names(mdf))) {

          # the variable of the response in the formula (e.g. log(y + 1))
          nresponse <- intersect(all.vars(stats::formula(model)[[2]]), names(mdf))[1]
          if (is.na(nresponse)) {
            stop("The back-transformed plot needs the response variable in the data; draw the means with ",
                 "autoplot() instead.", call. = FALSE)
          }
        }

        bky <- mdf[, nresponse]
      }

      if (is.list(bky)){
        bky <- bky[[1]]
      }

      if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
        Trtn <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
        newdata2 <- data.frame(bky=bky, Trtn=Trtn)
        bkmt$Trt <- do.call("paste", c(bkmt[, vars, drop=FALSE], sep=":"))
      } else {
        Trtn <- do.call("paste", c(mdf[, c(atvar, resvar), drop=FALSE], sep=":"))
        newdata2 <- data.frame(bky=bky, Trtn=Trtn)
        bkmt$Trt <- do.call("paste", c(bkmt[, c(atvar, resvar), drop=FALSE], sep=":"))
      }

      xMax <- max(max(bkmt[, nc], na.rm=TRUE), bky, na.rm=TRUE)
      xMin <- min(min(bkmt[, nc - 1], na.rm=TRUE), bky, na.rm=TRUE)
      xoffSet <- 0.15 * (xMax - xMin)
      mtitle <- plottitle

      if (is.null(plottitle) || plottitle%in%c("NULL", "")){
        mtitle <- paste("Back Transformed Means with ", (1 - ci_level) * 100, "% CIs\n for '",
                        modelterm, "'", "\n", sep = "")
      }


      xlimv <- c(xMin - xoffSet, xMax + xoffSet)
	    bkmt[, c("Mean", "LL", "UL")] <- lapply(bkmt[, c("Mean", "LL", "UL")], as.numeric)
	    newdata2$bky <- as.numeric(newdata2$bky)
	    p <- ggplot(bkmt, aes(Mean, Trt))+
        labs(title=pm_wrap_title(mtitle, basesz), x="", y="")+
        xlim(xlimv) +
        geom_point(colour="red") + geom_errorbar(aes(xmax = UL, xmin=LL ), width=0.2, linewidth=0.8, colour="red", orientation="y") +
        scale_y_discrete(limits = rev(unique(bkmt$Trt)))+
        geom_point(aes(x=bky, y=Trtn), shape=1, position = position_jitter(width = jitterv, height = jitterv), colour="blue", alpha=0.6, data=newdata2)+
        theme_bw(basesz)
      ciPlot <<- p

      print(p)
    }

  })

  list(meanPlot = meanPlot, ciPlot = ciPlot, predictmeansBarPlot = predictmeansBarPlot)
}

## Stable tables and the returned list, from the results of the stages
## ('v', pm_flat()) and the plots.
pm_assemble <- function(v, plots) {
  # Stable result tables (means, pairs, settings).
  tidy_parts <- pm_tidy_parts(v)
  x <- stats::setNames(lapply(pm_legacy_inputs, function(n) v[[n]]), pm_legacy_inputs)
  x$predictmeansPlot <- list(meanPlot = plots$meanPlot, ciPlot = plots$ciPlot)
  x["predictmeansBarPlot"] <- list(plots$predictmeansBarPlot)
  pm_finish(pm_legacy_output(x), tidy_parts)
}

## Values the legacy result list is built from.
pm_legacy_inputs <- c("pairwise", "atvar", "permlist", "adj", "ndecimal", "mean.table", "se.table",
                      "SED.out", "LSD", "LSDm", "t.p.valuem", "predictmeansPlot", "predictmeansBarPlot",
                      "meanTable", "p_valueMatrix", "listlength", "resvar", "pmlistTab", "atvar.levels", "nsim")

## The legacy result list of predictmeans() (pure function); 'x' is a named
## list of the values in pm_legacy_inputs.
pm_legacy_output <- function(x) {
  pairwise <- x$pairwise; atvar <- x$atvar; permlist <- x$permlist; adj <- x$adj
  ndecimal <- x$ndecimal; mean.table <- x$mean.table; se.table <- x$se.table
  SED.out <- x$SED.out; LSD <- x$LSD; LSDm <- x$LSDm; t.p.valuem <- x$t.p.valuem
  predictmeansPlot <- x$predictmeansPlot; predictmeansBarPlot <- x$predictmeansBarPlot
  meanTable <- x$meanTable; p_valueMatrix <- x$p_valueMatrix; listlength <- x$listlength
  resvar <- x$resvar; pmlistTab <- x$pmlistTab; atvar.levels <- x$atvar.levels; nsim <- x$nsim

  if (pairwise) {

    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {

      if (all(!pm_perm_given(permlist), adj %in% c("none", "bonferroni"))) {
        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,ndecimal+1),
                           "Pairwise p-value" = round(t.p.valuem, 4), predictmeansPlot=predictmeansPlot,
                           predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable, p_valueMatrix=p_valueMatrix)
                return(outputlist)

      }else{

        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4),
                           predictmeansPlot=predictmeansPlot, predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable,
                           p_valueMatrix=p_valueMatrix)
                return(outputlist)
      }
    } else {

      outputlist <- vector("list", listlength+8)
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
      outputlist[[listlength+7]] <- predictmeansPlot
      outputlist[[listlength+8]] <- predictmeansBarPlot
      outputlist[[listlength+9]] <- p_valueMatrix

      # print(outputlist)
      if (!pm_perm_given(permlist)) {
        names(outputlist)<- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                              "LSD", paste("Pairwise comparison p-value (adjusted by '", adj, "' method)", sep=""),
                              atvar.levels, "mean_table", "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[seq_along(outputlist)]
      } else {
        names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                 "Approximated LSD"), paste("Pairwise '", nsim, "' times permuted p-value (adjusted by '", adj, "' method)", sep=""),
                               atvar.levels, "mean_table", "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[seq_along(outputlist)]
      }

            return(outputlist)
    }
  } else {
    outputlist=list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD, predictmeansPlot=predictmeansPlot,
                    predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable)
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
