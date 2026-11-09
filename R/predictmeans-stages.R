## The stages of predictmeans(), in this order:
##   setup          check arguments and options; link function; plot choices
##   estimate       pm_model_parms(), pm_estimate_means(): K matrix,
##                  predicted means, SEs, legacy tables
##   compare        pm_compare_sed(), pm_compare_df(), pm_compare_lsd(),
##                  pm_pair_pvalues(), pm_pairs_overall() / pm_pairs_atvar()
##   mean_table     pm_mean_table(): the mean table with confidence limits
##   backtransform  pm_back_transform(), pm_final_mean_table(): response
##                  scale; letters
##   plots          pm_make_plots(): all graphics (only when plot = TRUE)
##   assemble       pm_tidy_parts(), pm_legacy_output(), pm_finish()
## The work is done by the pure functions named above: they take their
## inputs as arguments and return their results. Each stage is a short
## quoted expression, evaluated in the state of one predictmeans() call
## ('st'), that calls them and stores what later stages use. Only setup
## works on the state directly, because it rewrites the call's arguments
## (pairwise, atvar, adj, level, the plot choices).
## All statistics are computed before any plotting, so plot = FALSE skips
## all graphics work.

## Arguments, options, link function and plot choices.
pm_stage_setup <- function(st) eval(pm_code_setup, st)
pm_code_setup <- quote({

  if (!(modelterm %in% attr(terms(model), "term.labels"))){
    stop(paste("The", modelterm, "must be exactly a term in the model (especially check the order of interaction)."))
  }

  all.null = function(...){
    all(unlist(lapply(list(...), is.null)))
  }

  if(all.null(mplot, barplot, pplot, bkplot)){
    which.plots = match.arg(which.plots, choices = c("mean", "bar", "pval", "back"), several.ok = TRUE)

    which.plots.num = c("mean", "bar", "pval", "back") %in% which.plots
    mplot = which.plots.num[1]
    barplot = which.plots.num[2]
    pplot = which.plots.num[3]
    bkplot = which.plots.num[4]
  }else{
    ## Once you remove the old arguments you remove the if statement and delete this

    old_args = c("mplot", "barplot", "pplot", "bkplot")
    is_null = vapply(list(mplot, barplot, pplot, bkplot), is.null, logical(1))
    notnull = old_args[!is_null]
    nn = length(notnull)

    if(nn == 1){
      arg_list = notnull
      verb = "is"
    }else{
      arg_list = paste(paste(notnull[-nn], collapse = ", "), "and", notnull[nn])
      verb = "are"
    }

    ## arguments not given keep their which.plots setting (the default
    ## draws all plots), so e.g. 'bar = TRUE' no longer turns the mean plot off
    which.plots = match.arg(which.plots, choices = c("mean", "bar", "pval", "back"), several.ok = TRUE)
    wp_default = c("mean", "bar", "pval", "back") %in% which.plots
    for(k in which(is_null)){
      assign(old_args[k], wp_default[k])
    }

    wp_now = c("mean", "bar", "pval", "back")[c(isTRUE(mplot), isTRUE(barplot), isTRUE(pplot), isTRUE(bkplot))]
    warnMessage = paste0(arg_list, " ", verb, " deprecated: use which.plots = c(",
                         paste0('"', wp_now, '"', collapse = ", "), ") instead.")
    warning(warnMessage, call. = FALSE)
  }

  meanPlot <- ciPlot <- predictmeansBarPlot <- NULL

  if (is.null(trans)) {
    trans <- I
  }

  # if (inherits(model, "aovlist")) stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "aovlist")){
    model <- aovlist_lmer(model)
  }

  # Models with a link function other than the identity: means are
  # back-transformed with its inverse (pm_linkinv(); see ?pm_coef to support
  # other model classes). Otherwise the user's 'trans' is used.
  linkinv <- pm_link_trans(model)
  link_scale <- !is.null(linkinv)
  if (link_scale) {
    trans <- linkinv
  }

  # Poisson models: the back-transformed means are counts
  if (inherits(model, c("glm", "glmerMod", "glmmTMB")) &&
      stats::family(model)$family %in% c("poisson", "quasipoisson")) {
    count <- TRUE
  }

  vars <- unlist(strsplit(modelterm, "\\:"))
  # Model frame; for model classes without one, the data recovered by
  # Kmatrix() are used instead (below).
  mdf <- tryCatch(model.frame(model), error = function(e) NULL)

  if(!is.null(mdf) && any(!is.element(vars, names(mdf)[sapply(mdf,is.factor)]))){
    stop(paste(vars, "must be factor(s)!"))
  }

  # option checking
  if (length(vars)==1){
    atvar <- NULL
  }

  if (!is.null(permlist) && !unique(permlist%in%c("NULL", ""))) {
    pairwise <- TRUE

    if (adj=="tukey"){
      stop("The p-value can't be adjusted by Tukey methd!")
    }
  }

  if (!is.null(atvar) && !unique(atvar%in%c("NULL", ""))){
    pairwise <- TRUE
  }

  if (adj != "none"){
    pairwise <- TRUE
  }

  if (letterCI) {
    atvar <- NULL
    pairwise <- TRUE
    adj <- "none"

    if (is.null(level)){
      slevel <- level <- 0.166
    }else{
      slevel <- level
    }
  }

  if (is.null(level)){
    slevel <- level <- 0.05
  }

  if (!is.null(plotord) && !unique(plotord%in%c("NULL", ""))) {
    plot <- mplot <- TRUE

    if(length(plotord) != length(vars)){
      stop(paste("plotord must be a vector of length", length(vars)))
    }

    if(is.character(plotord)){
      i = match(plotord, vars)

      if(any(is.na(i))){
        stop(paste("The entries", paste0(plotord[is.na(i)], sep = ", "), "do not match any term in ", modelterm))
      }else{
        plotord = i
      }
    }
  }

  if (!is.logical(meandecr)){
    meandecr <- NULL
  }

  if (!prtplt){
    newwd <- FALSE
  }

})

## K matrix, predicted means, standard errors and the legacy mean/SE tables.
pm_stage_estimate <- function(st) eval(pm_code_estimate, st)
pm_code_estimate <- quote({
  ctr.matrix <- Kmatrix(model, modelterm, covariate, data=data, prtnum=prtnum,
                        cov_reduce=cov_reduce, exclude=exclude, weights=weights)
  KK <- ctr.matrix$K
  if (is.null(mdf)) {
    mdf <- ctr.matrix$data
  }
  response <- ctr.matrix$response
  parms <- pm_model_parms(model, df_method, vcov_adjust, vcov_adjust_default)
  mp <- parms$mp
  df_type <- parms$df_type
  vcov_adj <- parms$vcov_adj
  est <- pm_estimate_means(KK, ctr.matrix$fctnames, mdf, vars, mp, ndecimal)
  K <- est$K
  label <- est$label
  pm <- est$pm
  ses <- est$ses
  ndf <- est$ndf          # cell counts, used by the stable tables
  mt <- bkmt <- est$mt    # bkmt: for back transformed
  vcovm <- mp$vcov
  mean.table <- est$mean.table
  se.table <- est$se.table
  LL <- UL <- NULL
})

## Coefficients, covariance and df method of the model (pure function).
## lmer: Kenward-Roger df for REML fits (Satterthwaite for ML fits) unless
## 'df_method' says otherwise; the adjusted covariance is computed once and
## shared by all df calculations. Returns list(mp, df_type, vcov_adj).
pm_model_parms <- function(model, df_method, vcov_adjust, vcov_adjust_default) {
  dfm <- pm_df_method(model, df_method)
  df_type <- dfm$df_type
  vcov_adj <- dfm$vcov_adj
  # containment df with the adjusted covariance only when asked for explicitly
  if (identical(df_type, "containment") && vcov_adjust && lme4::getME(model, "is_REML"))
    vcov_adj <- tryCatch(pm_vcov_kr(model), error = function(e) NULL)
  mp <- mymodelparm(model, df_type = df_type, vcov_adj = vcov_adj)
  if (vcov_adjust) {
    if (is.null(vcov_adj) || inherits(vcov_adj, "pm_satt")) {
      # the default applies to REML lmer fits only; say so only when asked for
      if (!vcov_adjust_default) warning("'vcov_adjust = TRUE' needs an lmer model fitted by REML; ignored.")
    } else {
      mp$vcov <- as.matrix(vcov_adj)
    }
  }
  list(mp = mp, df_type = df_type, vcov_adj = vcov_adj)
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

  if (any(ndf$Freq==0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rownames(K)%in%rnTrt, , drop=FALSE]
    label <- label[rownames(label)%in%rnTrt, , drop=FALSE]
  }

  pm <- K %*% mp$coef
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  # Means that cannot be estimated (they depend on aliased coefficients) are
  # NA; their differences may still be estimable and are kept.
  est_flag <- attr(KK, "estimable")
  if (!is.null(est_flag)) {
    non_est <- !est_flag[rownames(K)]
    non_est[is.na(non_est)] <- FALSE
    if (any(non_est)) {
      pm[non_est] <- NA
      ses[non_est] <- NA
      message("Means that cannot be estimated are shown as NA: ",
              paste(rownames(K)[non_est], collapse = ", "),
              ". Their pairwise differences are shown when estimable. ",
              "If factors are confounded or cells are empty, weights = \"cells\" ",
              "averages over the observed cells only.")
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
  mean.table[!(n.table)] <- NA
  se.table[!(n.table)] <- NA
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
pm_stage_compare <- function(st) eval(pm_code_compare, st)
pm_code_compare <- quote({
  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)

  if (nK == 1) {
    SED.out <- NA
    LSD <- NA
  }else {
    perm <- !(is.null(permlist) || all(permlist%in%c("NULL", "")))
    list2env(pm_compare_sed(K, KK, label, ctr.matrix$fctnames, vars, mp$coef, vcovm, adj, level, atvar),
             envir = environment())
    if (!perm) {
      list2env(pm_compare_df(model, modelterm, K, rK, mp, mt, Df, df_type, vcov_adj), envir = environment())
      bkmt <- mt  # update bkmt
    }
    LSD <- pm_compare_lsd(SED.out, dses.m, vars, level, Df, slevel, ndecimal, perm)

    if (pairwise) {
      p_valueMatrix <- NULL
      list2env(pm_pair_pvalues(model, rK, cm, dses, nK, Df, Df_diff, adj, perm, permlist, ncore),
               envir = environment())
      if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
        list2env(pm_pairs_overall(model, rnK, mt, cm, dses, t.v, t.p.values, p.tukey, nK, Df, pairDf,
                                  adj, perm, nsim, level, slevel, meandecr), envir = environment())
      } else {
        list2env(pm_pairs_atvar(dses.df, KKvarndiff, labKK, mt, mdf, vars, atvar, t.v, t.p.values,
                                Df_diff, adj, slevel), envir = environment())
      }
    }
  } # end of if (nK == 1) ... else
})

## The pure functions of the compare stage. Each returns a named list of the
## values it creates (only those that exist), which the stage stores.
pm_existing <- function(nms, env) mget(intersect(nms, ls(env, all.names = TRUE)), envir = env)

## Pairwise contrasts and their SEDs (overall and within levels of each
## factor). 'atvar' becomes NULL when a factor has no pairs within its levels
## (nested factors).
pm_compare_sed <- function(K, KK, label, fctnames, vars, coef, vcovm, adj, level, atvar) {
  nK <- nrow(K)
  rnK <- rownames(K)
  # pairs (i, j), i > j, in column order; level labels come from the label
  # tables, not from splitting row names at ":" (levels may contain ":")
  pp <- pm_pair_index(nK)
  CM <- pm_pair_matrix(nK)
  labK <- pm_label_matrix(label, length(vars))

  nKK <- nrow(KK)
  ppKK <- pm_pair_index(nKK)
  labKK <- pm_label_matrix(fctnames, length(vars))

  KKvarndiff <- data.frame(labKK[ppKK$i, , drop=FALSE], labKK[ppKK$j, , drop=FALSE])

  rK <- CM%*%K                    # calculate stats
  cm <- rK%*%coef
  vcov.contr <- rK %*% tcrossprod(vcovm, rK)
  dses <- sqrt(diag(vcov.contr))
  # differences that cannot be estimated (their components along the
  # non-estimable directions do not cancel) are NA, not computed with the
  # dropped coefficients taken as zero
  null_proj <- attr(KK, "null_proj")
  if (!is.null(null_proj) && all(rnK %in% rownames(null_proj))) {
    pd <- CM %*% null_proj[rnK, , drop = FALSE]
    scale_r <- pmax(1, apply(abs(rK), 1, max))
    non_est_pair <- apply(abs(pd), 1, max) > 1e-6 * scale_r
    cm[non_est_pair] <- NA
    dses[non_est_pair] <- NA
  }

  if (adj == "bonferroni"){
    level <- level/length(dses)
  }

  SED.out <- if (all(is.na(dses))) c(Max.SED = NA_real_, Min.SED = NA_real_, Aveg.SED = NA_real_) else
    c(Max.SED = max(dses, na.rm = TRUE), Min.SED = min(dses, na.rm = TRUE), Aveg.SED = mean(dses, na.rm = TRUE))
  dses.df <- data.frame(labK[pp$i, , drop=FALSE], labK[pp$j, , drop=FALSE], dses)
  dses.m <- NULL

  if (length(vars) > 1) {
    dses.m <- matrix(0, nrow=3, ncol=length(vars))
    colnames(dses.m) <- vars
    rownames(dses.m) <- c("Aveg.SED", "Min.SED", "Max.SED")

    for (i in 1:length(vars)) {
      varsindx <- as.character(dses.df[,i])==as.character(dses.df[, length(vars)+i])
      if(any(varsindx)){  # in case of A/B
        dses.m[,i] <- summary.default(dses.df[varsindx,"dses"])[c(4,1,6)]
      }else{
        dses.m[,i] <- NA
        atvar <- NULL
      }
    }
    attr(SED.out, "For the Same Level of Factor") <- dses.m
  }
  list(rK = rK, cm = cm, dses = dses, level = level, SED.out = SED.out, dses.df = dses.df,
       dses.m = dses.m, KKvarndiff = KKvarndiff, labKK = labKK, atvar = atvar)
}

## Degrees of freedom of the means and of the pairs (no permutation tests):
## the user's Df, else by model class.
pm_compare_df <- function(model, modelterm, K, rK, mp, mt, Df, df_type, vcov_adj) {
      if (length(Df) == 0) {

        if (inherits(model, "lme")) {
          Df <- terms(model$fixDF)[modelterm]
    		  names(Df) <- NULL
    		  mt$Df <- round(Df, 2)
    		  pairDf <- Df
    		  Df_diff <- Df
        } else if (inherits(model, "lmerMod")) {
        #  Df <- median(df_term(model, modelterm), na.rm=TRUE)
    		  Df <- mp$df[modelterm]
    		  names(Df) <- NULL
    		  pm_df_raw <- df_term(model, modelterm, ctrmatrix=K, type=df_type, vcov_adj=vcov_adj)
    		  mt$Df <- round(pm_df_raw, 2)
    		 # if (any(mt$Df <= 0)) mt$Df <- Df
    		  Df_diff <- df_term(model, modelterm, ctrmatrix=rK, type=df_type, vcov_adj=vcov_adj)
    		  pairDf <- round(mean(Df_diff), 2)
    		 # Df_diff[Df_diff <= 0] <- 1
        } else {
          # a pm_df() method can give the df of each contrast (see ?pm_coef)
          df_means <- pm_df(model, K)
          if (!is.null(df_means)) {
            pm_df_raw <- rep_len(as.numeric(df_means), nrow(K))
            mt$Df <- round(pm_df_raw, 2)
            Df_diff <- rep_len(as.numeric(pm_df(model, rK)), nrow(rK))
            Df <- mean(Df_diff)
            pairDf <- round(Df, 2)
            df_type <- "model (pm_df)"
          } else {
            Df <- mp$df
            mt$Df <- round(Df, 2)
            pairDf <- Df
            Df_diff <- Df
          }
        }

        if (Df==0){
          stop("You need provide Df for this model!")
        }
      }else {
  	    df_user <- TRUE
  	    mt$Df <- Df
  		  Df_diff <- rep_len(Df, nrow(rK))   # one df per pair, as for lmer
  		  pairDf <- Df
		  }

  pm_existing(c("Df", "mt", "pm_df_raw", "Df_diff", "pairDf", "df_type", "df_user"), environment())
}

## LSDs from the SEDs: t quantile with the df, or 2 x SED with permutation tests.
pm_compare_lsd <- function(SED.out, dses.m, vars, level, Df, slevel, ndecimal, perm) {
  if (!perm) {
    LSD <- round(qt(1 - level/2, df = Df) * SED.out, ndecimal+1)
    names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
    attr(LSD, "For the Same Level of Factor") <- NULL

    if (length(vars) > 1) {
      rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
      attr(LSD, "For the Same Level of Factor") <- round(qt(1 - level/2, df = Df) * dses.m, ndecimal+1)
    }

    attr(LSD, "Significant level") <- slevel
    attr(LSD, "Degree of freedom") <- round(Df, 2)
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

## t values and raw p-values of the pairs (t, Tukey or permutation).
pm_pair_pvalues <- function(model, rK, cm, dses, nK, Df, Df_diff, adj, perm, permlist, ncore) {
      t.v <- cm/dses

      if (!perm && adj=="tukey") {
	      if (inherits(model, "lmerMod")) {
		      p.tukey <- sapply(1:length(t.v), function(m) ptukey(sqrt(2)*abs(t.v[m]), nK, Df_diff[m], lower.tail=FALSE))
		    }else{
	        p.tukey <- ptukey(sqrt(2)*abs(t.v), nK, Df, lower.tail=FALSE)
		    }
      }

      if (!perm) {
	  	  if (inherits(model, "lmerMod")) {
		      t.p.values <- sapply(1:length(t.v), function(m) 2 * pt(-abs(t.v[m]), Df_diff[m]))
		    }else{
	        t.p.values <- 2 * pt(-abs(t.v), Df)
		    }
      }else{
        nsim <- length(permlist[[1]])

        tValue <- function(x, rK){
          cm <- rK%*%x$coef
          vcovm <- x$vcov
          vcov.contr <- rK %*% tcrossprod(vcovm, rK)
          ses <- sqrt(diag(vcov.contr))
          t.v <- cm/ses
          return(t.v)
        }

        if (.Platform$OS.type=="windows") {
          cl <- makeCluster(ncore)
          clusterExport(cl, c("tValue", "rK", "t.v"), envir = environment())
          t.TableL <- parLapplyLB(cl, permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6))
          stopCluster(cl)
        }else{
          t.TableL <- mclapply(permlist[[1]], function(x) {
            round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)
          }, mc.cores=ncore)
        }

        t.TableL <- matrix(unlist(t.TableL), ncol = length(t.TableL))
        t.p.values <-  (rowSums(t.TableL)+1)/(nsim+1)

      } # end of if (is.null(permlist))

  pm_existing(c("t.v", "p.tukey", "t.p.values", "nsim"), environment())
}

## Overall pairwise tables: p-values (t value above the diagonal), LSDs and
## the matrix used for letter groups.
pm_pairs_overall <- function(model, rnK, mt, cm, dses, t.v, t.p.values, p.tukey, nK, Df, pairDf,
                             adj, perm, nsim, level, slevel, meandecr) {
  tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
  rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
  tvm[upper.tri(tvm)] <- t.v
  p_valueMatrix <- NULL

        if (adj=="tukey"){
          t.p.valuem[upper.tri(t.p.valuem)] <- p.tukey
        }else{
          t.p.valuem[upper.tri(t.p.valuem)] <- p.adjust(t.p.values, adj)
        }

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
          LSDm.up <- qt(1-level/2, df = Df)*dses
          LSDm[upper.tri(LSDm)] <- LSDm.up
          Diffm[upper.tri(Diffm)] <- cm
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

## Pairwise tables within each level of 'atvar': adjusted p-values, letters
## and the means ordered by 'atvar'.
pm_pairs_atvar <- function(dses.df, KKvarndiff, labKK, mt, mdf, vars, atvar, t.v, t.p.values,
                           Df_diff, adj, slevel) {
  p_valueMatrix <- NULL
      dses.df$tvalue <- t.v
      dses.df$pvalue <- t.p.values
      if (adj == "tukey") {
        dses.df$pair_df <- rep_len(Df_diff, nrow(dses.df))   # df of each pair (varies for lmer)
      }

      for (i in which(vars%in%atvar)) { # To find rows relating atvar
        KKvarndiff <- KKvarndiff[which(as.character(KKvarndiff[, i]) == as.character(KKvarndiff[, length(vars) + i])), ]
      }

      atvar.df <- f_loj_krc(KKvarndiff, dses.df, by.x=names(KKvarndiff), by.y=names(KKvarndiff))
      names(atvar.df)[1:length(vars)] <- vars

		  for (i in vars) {       # To ensure factor vars  have the same level as original
          atvar.df[,i] <- factor(atvar.df[,i], levels=levels(mdf[, i][[1]]))
      }

      # adjusted p-values within each atvar group, kept in the rows they
      # belong to (the group's own rows, not split()/unlist() order)
      key_atv <- do.call("paste", c(lapply(atvar.df[, atvar, drop=FALSE], as.character), sep=" : "))
      atvar.df$adj.pvalue <- NA_real_
      for (g in unique(key_atv)) {
        idx <- which(key_atv == g)
        x <- atvar.df[idx, , drop=FALSE]
        t_valueN <- length(na.omit(x$tvalue))
        if (t_valueN < 2) {
          atvar.df$adj.pvalue[idx] <- x$pvalue
        } else if (adj=="tukey") {
          # Tukey's range distribution needs the number of means, k, in
          # this level; the group holds all k(k-1)/2 pairs among them.
          n_means <- (1 + sqrt(1 + 8*t_valueN))/2
          atvar.df$adj.pvalue[idx] <- ptukey(sqrt(2)*abs(x$tvalue), n_means, x$pair_df, lower.tail=FALSE)
        } else {
          atvar.df$adj.pvalue[idx] <- p.adjust(x$pvalue, adj)
        }
      }

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

## Mean table with confidence limits.
pm_stage_mean_table <- function(st) eval(pm_code_mean_table, st)
pm_code_mean_table <- quote({
  mt_res <- pm_mean_table(mt, vars, perm = !(is.null(permlist) || all(permlist%in%c("NULL", ""))),
                          slevel = slevel, ndecimal = ndecimal)
  meanTable <- mt_res$table
  slevel <- mt_res$slevel
})

## Mean table with confidence limits (pure function).
## mt: data frame with columns pm, ses, Df and the factors 'vars'.
## With permutation p-values there are no df: limits are mean +/- 2 SE and
## the level becomes 0.05. Returns list(table, slevel).
pm_mean_table <- function(mt, vars, perm, slevel, ndecimal) {
  meanTable <- mt
  if (!perm) {
    meanTable$LL <- meanTable$pm - qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
    meanTable$UL <- meanTable$pm + qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
  } else {
    meanTable$LL <- meanTable$pm - 2 * meanTable$ses
    meanTable$UL <- meanTable$pm + 2 * meanTable$ses
    slevel <- 0.05
    meanTable$Df <- NA
  }
  rownames(meanTable) <- NULL
  meanTable <- meanTable[c(vars, "pm", "ses", "Df", "LL", "UL")]
  meanTable[c("pm", "LL", "UL")] <- round(meanTable[c("pm", "LL", "UL")], ndecimal)
  meanTable$ses <- round(meanTable$ses, ndecimal+1)
  colnames(meanTable) <- c(vars, "Mean", "SE", "Df", paste("LL(", (1 - slevel) * 100, "%)", sep = ""),
                           paste("UL(", (1 - slevel) * 100, "%)", sep = ""))
  list(table = meanTable, slevel = slevel)
}

## Back-transformed means, and letters for the mean table.
pm_stage_backtransform <- function(st) eval(pm_code_backtransform, st)
pm_code_backtransform <- quote({
  if (!is.null(trans)) {
    Mean <- Trt <- ciPlot <- NULL
    perm <- !(is.null(permlist) || all(permlist%in%c("NULL", "")))
    bkmt <- pm_back_transform(bkmt, trans, transOff, perm = perm, slevel = slevel,
                              ndecimal = ndecimal, count = count)
    nc <- ncol(bkmt)   # the plots use bkmt[, nc - 1] and bkmt[, nc] (LL, UL)
    meanTable <- pm_final_mean_table(meanTable, bkmt, vars, atvar, resvar, pairwise, meandecr,
                                     letterCI, slevel, ndecimal, t.p.valuemGrp, pmlistLetter)
  }
})

## Means and confidence limits on the response scale (pure function).
## bkmt has columns pm, ses, Df (absent with permutation p-values) and the
## factors; returns the factors with Mean, LL and UL.
pm_back_transform <- function(bkmt, trans, transOff, perm, slevel, ndecimal, count) {
  bkmt$Mean <- trans(bkmt$pm)-transOff

  if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)){
    bkmt$Mean <- exp(bkmt$pm)-transOff
  }

  if (!perm) {
    bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)-transOff
    bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)-transOff
  } else {
    bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)-transOff
    bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)-transOff
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
## to exist) in the branches that need them.
pm_final_mean_table <- function(meanTable, bkmt, vars, atvar, resvar, pairwise, meandecr,
                                letterCI, slevel, ndecimal, t.p.valuemGrp, pmlistLetter) {
  no_atvar <- is.null(atvar) || all(atvar%in%c("NULL", ""))
  if (!no_atvar) {
    all_var_names <- colnames(meanTable)
    meanTable <- meanTable[c(atvar, resvar, setdiff(all_var_names, vars))]
    meanTable <- with(meanTable, meanTable[order(eval(parse(text=c(atvar, resvar)[length(vars):1]))),])
  }

  if (any(meanTable$Mean!=bkmt$Mean, na.rm = TRUE)) {
    meanTable <- cbind(meanTable, round(bkmt[, (ncol(bkmt)-2):ncol(bkmt)], ndecimal))
    colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("Bk_Mean", paste("Bk_LL(", (1 - slevel) * 100, "%)", sep = ""),
                                                                   paste("Bk_UL(", (1 - slevel) * 100, "%)", sep = ""))
  }

  if (pairwise) {
    if (no_atvar) {
      if (!is.null(meandecr) && is.logical(meandecr)) {
        meanTable <- meanTable[order(meanTable$Mean, decreasing = meandecr), ]
      }

      # letters only among means whose comparisons can all be estimated
      meanTable$LetterGrp <- unname(pm_letter_groups(t.p.valuemGrp, rep(TRUE, nrow(t.p.valuemGrp)), slevel))

      if (letterCI) {
        meanTable$LetterGrp <- ci_mcp(meanTable[, grepl("^LL", names(meanTable))],
                                      meanTable[, grepl("^UL", names(meanTable))])
      }

      meanTable <- list(Table=meanTable, Note=paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep=""))
    } else {
      lets <- unlist(pmlistLetter)
      mkey <- paste(do.call("paste", c(lapply(meanTable[atvar], as.character), sep=" : ")),
                    do.call("paste", c(lapply(meanTable[resvar], as.character), sep=" : ")), sep="\r")
      meanTable$LetterGrp <- unname(lets[mkey])
      meanTable <- list(Table=meanTable, Note=paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "' at each level of ", atvar, sep=""))
    }
  }
  meanTable
}

## All graphics: p-value plots, mean and bar plots, back-transformed plot.
pm_stage_plots <- function(st) eval(pm_code_plots, st)
pm_code_plots <- quote({
  plots <- pm_make_plots(mget(pm_plot_inputs, envir = environment(), ifnotfound = list(NULL)))
  meanPlot <- plots$meanPlot
  ciPlot <- plots$ciPlot
  predictmeansBarPlot <- plots$predictmeansBarPlot
})

## Values the plots are drawn from.
pm_plot_inputs <- c("LL", "LSD", "Mean", "SED.out", "Trt", "UL", "adj", "atvar", "barplot", "basesz",
                    "bkmt", "bkplot", "ciPlot", "jitterv", "lineplot", "lsd_bar", "mdf", "meanPlot",
                    "model", "modelterm", "mplot", "mt", "nK", "nc", "newwd", "nrow.pm", "pairwise",
                    "plot", "plotord", "plottitle", "plotxlab", "plotylab", "pm", "pmlist", "pplot",
                    "predictmeansBarPlot", "prtplt", "rcnplotm", "response", "responsen", "resvar",
                    "ses", "slevel", "t.p.valuep", "trans", "transOff", "link_scale", "vars")

## All graphics of predictmeans(): draws the requested plots and returns
## list(meanPlot, ciPlot, predictmeansBarPlot). 'x' is a named list of the
## values in pm_plot_inputs; the drawing code runs in a new environment that
## holds only those values.
pm_make_plots <- function(x) {
  env <- list2env(x, envir = new.env(parent = environment(pm_make_plots)))
  eval(pm_code_draw, env)
  mget(c("meanPlot", "ciPlot", "predictmeansBarPlot"), envir = env)
}

pm_code_draw <- quote({
  # p-value level plots (base graphics)
  if (pairwise && nK > 1) {
    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
      if (all(nrow(t.p.valuep) > 2, pplot, plot, prtplt)) {
        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Level Plot of p-value (adjusted by '", adj, "' method)\n for Pairwise Comparison", sep="")
        }

        PMplot(t(t.p.valuep), level=slevel, legendx=0.69, mtitle=mtitle, newwd=newwd)
      }
    } else {
      if (all(nrow.pm > 2, pplot, plot, prtplt)) {
        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Adjusted p-value (by '", adj,
                          "' method)\n for Pairwise Comparison at Each Level of '",paste(atvar, collapse =" and "), "'", sep="")
        }

        PMplot(pmlist, level=slevel, xylabel=rcnplotm, legendx=0.69, mtitle=mtitle, newwd=newwd)
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
        if (newwd){
          dev.new()
        }

        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for \"", vars, "\" with ", bar_label, " (", slevel * 100, "%) Bar", sep="")
        }

        p1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm, group=1))+
            labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
            lims(x= c(bar_label, levels(plotmt[, vars])), y = c(yMin - offSet, max(yMax + offSet, yMin + LSD_value + offSet))) +
            geom_point(colour="red", size=2)+
          #  geom_line(linewidth=0.5)+
            geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue") +  #data=lsdBar,
            theme_bw(basesz)

  		  if (lineplot){
  		    p1 <- p1 + geom_line(linewidth=0.5)
  		  }

        meanPlot <- p1

        if (prtplt){
          print(p1)
        }
      } # end of if mplot

      if (barplot) {

        if (newwd){
          dev.new()
        }

        mtitle <- plottitle

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
        }

        dodge <- position_dodge(width=0.9)
        bp1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm))+
            labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
            ylim(c(min(0, (min(pm) - max(ses)) * 1.2), max(0, (max(pm) + max(ses)) * 1.2))) +
            geom_bar(position=dodge, stat="identity", fill="papayawhip", colour="darkgreen") +
            geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())

        predictmeansBarPlot <- bp1
        if (prtplt){
          print(bp1)
        }
      }
    } else if (length(vars) == 2) {
      ##NB: I have changed this to an else statement because there shouldn't be any other way to get here. That is, length(vars) will not
      ##change from 1 to 2 between the block where length(vars)==1 is true and here, therefore the cases should be else if statements.
        if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
          plotord <- 1:2
          if (!(is.null(atvar) || all(atvar%in%c("NULL", "")))) {
            atvar_num <- which(vars %in% atvar)
            plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
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

          if (newwd){
            dev.new()
          }

          mtitle <- plottitle

          if (is.null(plottitle) || plottitle%in%c("NULL", "")){
            mtitle <- paste("Predicted means for \"", fact1, "\" by \"", fact2, "\" with ",
                            paste(atvar, collapse=" "), " ", bar_label, " (", slevel * 100, "%) Bar", sep = "")
          }

          plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))
          p2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), col=eval(parse(text = fact2))))+
            labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
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

          if (prtplt){
            print(p2)
          }
        } # end if mplot

        if (barplot) {

          if (newwd) {
            dev.new()
          }

          mtitle <- plottitle

          if (is.null(plottitle) || plottitle%in%c("NULL", "")) {
            mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
          }

          dodge <- position_dodge(width=0.9)

          bp2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), fill=  eval(parse(text = fact2))))+
            geom_point(position=dodge) +
            geom_bar(stat = "identity", position=dodge)+
            labs( x = mxlab, y = mylab, title =mtitle)+
            ylim(c(min(0, (min(pm) - max(ses)) * 1.1), max(0, (max(pm) + max(ses)) * 1.1))) +
            geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            scale_fill_brewer()+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
            theme(legend.position = "top")+guides(fill = guide_legend(title = fact2))
          predictmeansBarPlot <- bp2

          if (prtplt) {
            print(bp2)
          }
        }
      } else if (all(length(vars) == 3, mplot)) {
        ## NB: I have changed this to else if, because neither length(vars)
        ## or mplot should have changed, and I believe this is just the third
        ## choice
        if (newwd) {
          dev.new()
        }

        mtitle <- plottitle

        if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
          plotord <- 1:3

          if (!(is.null(atvar) || all(atvar%in%c("NULL", ""))) && length(atvar)==1){
            atvar_num <- which(vars %in% atvar)
            plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
          }
        }

        fact1 <- (vars[plotord])[1]
        fact2 <- (vars[plotord])[2]
        fact3 <- (vars[plotord])[3]
        plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))

        if (is.null(plottitle) || plottitle%in%c("NULL", "")){
          mtitle <- paste("Predicted means for '", fact1, "' by '", fact2, "' for each '",
                           fact3, "'\n with ", paste(atvar, collapse=" "), " ",
                          bar_label, " (", slevel * 100, "%) Bar\n", sep = "")
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
          labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
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

        if (prtplt){
          print(p3)
        }
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
            bky <- trans(mdf[, response]) - transOff
          }# end of if glm or glmer
        }# end of if factor
      } else {       ## Transformed y within modelling
        nresponse <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]

        if (!(nresponse %in% names(mdf))) {

          if (is.null(responsen) || all(responsen%in%c("NULL", ""))){
            stop(paste("Please provide suitable name for response variable using option responsen='",
                       names(mdf)[1], "'!", sep=""))
          }

          nresponse <- responsen
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
        mtitle <- paste("Back Transformed Means with ", (1 - slevel) * 100, "% CIs\n for '",
                        modelterm, "'", "\n", sep = "")
      }

      if (newwd) {
        dev.new()
      }

      xlimv <- c(xMin - xoffSet, xMax + xoffSet)
	    bkmt[, c("Mean", "LL", "UL")] <- lapply(bkmt[, c("Mean", "LL", "UL")], as.numeric)
	    newdata2$bky <- as.numeric(newdata2$bky)
	    p <- ggplot(bkmt, aes(Mean, Trt))+
        labs(title=mtitle, x="", y="")+
        xlim(xlimv) +
        geom_point(colour="red") + geom_errorbarh(aes(xmax = UL, xmin=LL ), height=0.2, linewidth=0.8, colour="red") +
        scale_y_discrete(limits = rev(unique(bkmt$Trt)))+
        geom_point(aes(x=bky, y=Trtn), shape=1, position = position_jitter(width = jitterv, height = jitterv), colour="blue", alpha=0.6, data=newdata2)+
        theme_bw(basesz)
      ciPlot <<- p

      if (prtplt) {
        print(p)
      }
    }

  })

})

## Stable tables and the returned list.
pm_stage_assemble <- function(st) eval(pm_code_assemble, st)
pm_code_assemble <- quote({
  predictmeansPlot <- list(meanPlot=meanPlot, ciPlot=ciPlot)

  # Stable result tables (means, pairs, settings) from the values above.
  tidy_parts <- pm_tidy_parts(environment())
  outputlist <- pm_legacy_output(mget(pm_legacy_inputs, envir = environment(),
                                      ifnotfound = list(NULL)))
  pm_out <- pm_finish(outputlist, tidy_parts)
})

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

      if (all(is.null(permlist) || all(permlist%in%c("NULL", "")), adj %in% c("none", "bonferroni"))) {
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
      if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {
        names(outputlist)<- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                              "LSD", paste("Pairwise comparison p-value (adjusted by '", adj, "' method)", sep=""),
                              atvar.levels, "mean_table", "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[1:length(outputlist)]
      } else {
        names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                 "Approximated LSD"), paste("Pairwise '", nsim, "' times permuted p-value (adjusted by '", adj, "' method)", sep=""),
                               atvar.levels, "mean_table", "predictmeansPlot", "predictmeansBarPlot", "p_valueMatrix")[1:length(outputlist)]
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
