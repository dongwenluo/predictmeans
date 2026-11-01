## The stages of predictmeans(). Each stage runs on the state of one
## predictmeans() call (its evaluation frame, 'st'), reading what earlier
## stages computed and adding its own results, in this order:
##   setup          check arguments and options; link function; plot choices
##   estimate       K matrix, predicted means, standard errors, legacy tables
##   compare        SEDs, LSDs, degrees of freedom, pairwise tests (overall or
##                  within atvar levels), letters within atvar levels
##   mean_table     the mean table with confidence limits
##   backtransform  back-transformed means; letters
##   plots          all graphics (only when plot = TRUE)
##   assemble       the stable tables and the returned list (pm_out)
## Each stage's code is a quoted expression evaluated in that state.
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

    ## set all NULL values to FALSE
    for(v in old_args[is_null]){
      assign(v, FALSE)
    }

    warnMessage = paste0(arg_list, " ", verb, " deprecated. Use which.plots in the future. The value will be honored for now, but may not always be")
    warning(warnMessage)
  }

  meanPlot <- ciPlot <- predictmeansBarPlot <- NULL

  if (is.null(trans)) {
    trans <- I
  }

  # if (inherits(model, "aovlist")) stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "aovlist")){
    model <- aovlist_lmer(model)
  }

  # Models with a link function: means are back-transformed with its inverse
  # (pm_linkinv(); see ?pm_coef to support other model classes).
  linkinv <- pm_linkinv(model)
  if (!is.null(linkinv)) {
    trans <- linkinv
  }

  if (inherits(model, "glm") && model$family$family %in% c("poisson", "quasipoisson")){
    count=TRUE
  }

  if (inherits(model, "glmmTMB") && model$modelInfo$family$family %in% c("poisson", "quasipoisson")){
    count=TRUE
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
  label <- ctr.matrix$fctnames
  rownames(label) <- rownames(KK)
  response <- ctr.matrix$response
  # lmer: Kenward-Roger df for REML fits (Satterthwaite for ML fits). The
  # adjusted covariance is computed once and shared by all df calculations.
  dfm <- pm_df_method(model, df_method)
  df_type <- dfm$df_type
  vcov_adj <- dfm$vcov_adj
  # containment df with the adjusted covariance only when asked for explicitly
  if (identical(df_type, "containment") && vcov_adjust && getME(model, "is_REML"))
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
  n.table <- table(mdf[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model

  K <- KK[, mp$estimable, drop = FALSE]          # To match coef names

  if (any(ndf$Freq==0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rownames(K)%in%rnTrt,]
    label <- label[rownames(label)%in%rnTrt,]
  }

  pm <- K %*% mp$coef
  vcovm <- mp$vcov
  # ses <- sqrt(diag(K %*% tcrossprod(vcovm, K)))
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  mt <- data.frame(pm, ses, label)
  LL <- UL <- NULL
  bkmt <- mt  # for back transformed

  mean.table <- round(xtabs(pm ~ ., mt[, c("pm", vars)], drop.unused.levels = TRUE), ndecimal)
  se.table <- round(xtabs(ses ~ ., mt[, c("ses", vars)], drop.unused.levels = TRUE), ndecimal+1)
  mean.table[!(n.table)] <- NA
  se.table[!(n.table)] <- NA

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

})

## SEDs, LSDs, degrees of freedom and pairwise comparisons.
pm_stage_compare <- function(st) eval(pm_code_compare, st)
pm_code_compare <- quote({
  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)
  varn1 <- varn2 <- rep(0, nK * (nK - 1)/2)

  if (nK == 1) {
    SED.out <- NA
    LSD <- NA
  }else {
    kindx <- 1:nK
    CM <-  matrix(0, nrow=nK * (nK - 1)/2, ncol=nK)
    t <- 1

    for (i in 2:nK) {                      # To construct pairwise comparison K matrix by col order
      for (j in 1:(i-1)) {
        CM[t, ] <- (kindx == j) - (kindx == i)
        varn1[t] <- rnK[i]
        varn2[t] <- rnK[j]
        t <- t+1
      }
    }

    nKK <- nrow(KK)          # To setup various matrix and row, col names
    rnKK <- rownames(KK)
    KKvarn1 <- KKvarn2 <- rep(0, nKK * (nKK - 1)/2)
    tt <- 1

    for (i in 2:nKK) {
      for (j in 1:(i-1)) {
        KKvarn1[tt] <- rnKK[i]
        KKvarn2[tt] <- rnKK[j]
        tt <- tt+1
      }
    }

    KKvarndiff <- data.frame(matrix(unlist(strsplit(KKvarn1, "\\:")), byrow=T, nrow=length(KKvarn1)),
                             matrix(unlist(strsplit(KKvarn2, "\\:")), byrow=T, nrow=length(KKvarn2)))

    rK <- CM%*%K                    # calculate stats
    cm <- rK%*%mp$coef
    vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    dses <- sqrt(diag(vcov.contr))

    if (adj == "bonferroni"){
      level <- level/length(dses)
    }

    SED.out <- c(Max.SED = max(dses), Min.SED = min(dses), Aveg.SED = mean(dses))
    dses.df <- data.frame(matrix(unlist(strsplit(varn1, "\\:")), byrow=T, nrow=length(varn1)),
                          matrix(unlist(strsplit(varn2, "\\:")), byrow=T, nrow=length(varn2)), dses)

    if (length(vars) > 1) { # all(length(vars) > 1, SED.out[1]!=SED.out[2])
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

    if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {

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

	    bkmt <- mt  # update bkmt

      LSD <- round(qt(1 - level/2, df = Df) * SED.out, ndecimal+1)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "For the Same Level of Factor") <- NULL

      if (length(vars) > 1) {
        rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
        attr(LSD, "For the Same Level of Factor") <- round(qt(1 - level/2, df = Df) * dses.m, ndecimal+1)
      } # end of if LSD

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


    if (pairwise) {
      p_valueMatrix <- NULL
      tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
      rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
      t.v <- cm/dses

      if (all(is.null(permlist) || all(permlist%in%c("NULL", "")), adj=="tukey")) {
	      if (inherits(model, "lmerMod")) {
		      p.tukey <- sapply(1:length(t.v), function(m) ptukey(sqrt(2)*abs(t.v[m]), nK, Df_diff[m], lower.tail=FALSE))
		    }else{
	        p.tukey <- ptukey(sqrt(2)*abs(t.v), nK, Df, lower.tail=FALSE)
		    }
      }

      tvm[upper.tri(tvm)] <- t.v

      if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {
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

      if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {

        if (adj=="tukey"){
          t.p.valuem[upper.tri(t.p.valuem)] <- p.tukey
        }else{
          t.p.valuem[upper.tri(t.p.valuem)] <- p.adjust(t.p.values, adj)
        }

        t.p.valuep <- t.p.valuem    # for plot
        t.p.valuem <- t(t.p.valuem) + tvm
        names(t.p.valuem) <- NULL
		    diag(t.p.valuem) <- 1

		    if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {

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

    }else{
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

      atvar.df$adj.pvalue <- unlist(lapply(split(atvar.df, atvar.df[, atvar]), function(x) {
          t_value <- x$tvalue
          t_valueN <- length(na.omit(t_value))

          if (t_valueN < 2){
            x$adj.pvalue <- x$pvalue
          } else {
            if (adj=="tukey"){
              # Tukey's range distribution needs the number of means, k, in
              # this level; the group holds all k(k-1)/2 pairs among them.
              n_means <- (1 + sqrt(1 + 8*t_valueN))/2
              x$adj.pvalue <- ptukey(sqrt(2)*abs(x$tvalue), n_means, x$pair_df, lower.tail=FALSE)
            }else{
              x$adj.pvalue <- p.adjust(x$pvalue, adj)
            }
          }
        }))

      rnK.df <- as.data.frame(matrix(unlist(strsplit(rnKK, "\\:")), byrow=T, nrow=nKK)) # To find the suitable names for pm
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
      indexlength <- nrow(atvar.df)/listlength
      nrow.pm <- length(resvar.levels)

      for (i in 1:listlength) {           # extract pvalues for each factor level
        atvar.pm <- matrix(0, nrow=nrow.pm, ncol=nrow.pm)
        atvar.pm[upper.tri(atvar.pm)] <- atvar.df$adj.pvalue[(indexlength*i-(indexlength-1)):(indexlength*i)]
        rownames(atvar.pm) <- colnames(atvar.pm) <- resvar.levels
        pmlist[[i]] <- t(atvar.pm)
        outtab <- round(as.table(t(atvar.pm)), 4)
        outtab[outtab < 0.0001] <- "0.0001"
        outtab[col(outtab)==row(outtab)] <- 1.0000
        outtab[upper.tri(outtab)] <-""
        Grpatvar.pm <- t(atvar.pm)+ atvar.pm

        if (all(!is.na(outtab))) {
    	   Grpatvar.letter <- multcompLetters(Grpatvar.pm, Letters=LETTERS, threshold=slevel)[resvar.levels]
    	   outtab <- as.table(cbind(outtab, Group=Grpatvar.letter))
    	  }else{
      	  Grpatvar.letterM <- rep(NA, nrow.pm)
      	  names(Grpatvar.letterM) <- resvar.levels
      	  atvar.rowname <- resvar.levels[!is.na(Grpatvar.pm[, 1])]
      	  Grpatvar.letterM[atvar.rowname] <- multcompLetters(Grpatvar.pm[atvar.rowname, atvar.rowname], Letters=LETTERS, threshold=slevel)[atvar.rowname]
      	  outtab <- as.table(cbind(outtab, Group=Grpatvar.letterM))
      	  Grpatvar.letter <- Grpatvar.letterM[atvar.rowname]
    	  }

	      pmlistLetter[[i]] <- Grpatvar.letter
        pmlistTab[[i]] <- outtab
      }

      if (nrow.pm > 2){
        p_valueMatrix <- pmlist
      }

    } # if (is.null(atvar))
  }# end of if(pairwise)
  } # end of if (nK == 1) ... else
})

## Mean table with confidence limits.
pm_stage_mean_table <- function(st) eval(pm_code_mean_table, st)
pm_code_mean_table <- quote({
  meanTable <- mt

  if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {
    meanTable$LL <- meanTable$pm - qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
    meanTable$UL <- meanTable$pm + qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
  }else{
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
})

## Back-transformed means, and letters for the mean table.
pm_stage_backtransform <- function(st) eval(pm_code_backtransform, st)
pm_code_backtransform <- quote({
  if (!is.null(trans)) {
    Mean <- Trt <- ciPlot <- NULL
    bkmt$Mean <- trans(bkmt$pm)-transOff

    if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)){
      bkmt$Mean <- exp(bkmt$pm)-transOff
    }

    if (is.null(permlist) || all(permlist%in%c("NULL", ""))) {
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

    if (!(is.null(atvar) || all(atvar%in%c("NULL", "")))) {
      all_var_names <- colnames(meanTable)
      meanTable <- meanTable[c(atvar, resvar, setdiff(all_var_names, vars))]
      meanTable <- with(meanTable, meanTable[order(eval(parse(text=c(atvar, resvar)[length(vars):1]))),])
    }

   # if (all(!identical(trans, function(x) x, ignore.environment=TRUE), !identical(trans, I, ignore.environment=TRUE))) {
	  if (any(meanTable$Mean!=bkmt$Mean)) {
      meanTable <- cbind(meanTable, round(bkmt[, (ncol(bkmt)-2):ncol(bkmt)], ndecimal))
      colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("Bk_Mean", paste("Bk_LL(", (1 - slevel) * 100, "%)", sep = ""),
                                                                     paste("Bk_UL(", (1 - slevel) * 100, "%)", sep = ""))
    }

	  if (pairwise) {

      if ((is.null(atvar) || all(atvar%in%c("NULL", "")))) {

        if (!is.null(meandecr) && is.logical(meandecr)) {
          meanTable <- meanTable[order(meanTable$Mean, decreasing = meandecr), ]
        }

        meanTable$LetterGrp <- multcompLetters(t.p.valuemGrp, Letters=LETTERS, threshold=slevel)

        if (letterCI) {
          meanTable$LetterGrp <- ci_mcp(meanTable[, grepl("^LL", names(meanTable))],
                                        meanTable[, grepl("^UL", names(meanTable))])
        }

        meanTable <- list(Table=meanTable, Note=paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep=""))
        # attr(meanTable, "Note") <- paste("Note: letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep="")
      } else {
	      meanTable$LetterGrp <- unlist(pmlistLetter)
        meanTable <- list(Table=meanTable, Note=paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "' at each level of ", atvar, sep=""))
	    }
	  }

  }
})

## All graphics: p-value plots, mean and bar plots, back-transformed plot.
pm_stage_plots <- function(st) eval(pm_code_plots, st)
pm_code_plots <- quote({
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

  # means and bar plots
  if (plot) {

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
    if (plot && bkplot) {

      if (!is.null(atvar) && !unique(atvar%in%c("NULL", ""))) {
        mdf <- mdf[do.call(order, mdf[, c(atvar, resvar), drop=FALSE]),]
      }

      if (response %in% names(mdf)) {    ## Transformed y before modelling

        if (inherits(mdf[, response], "factor")) {
          bky <- as.numeric(mdf[, response])-1
        } else {

          if (inherits(model, "glm") || inherits(model, "glmerMod") || inherits(model, "glmmTMB")) {
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

  if (pairwise) {

    if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {

      if (all(is.null(permlist) || all(permlist%in%c("NULL", "")), adj %in% c("none", "bonferroni"))) {
        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,ndecimal+1),
                           "Pairwise p-value" = round(t.p.valuem, 4), predictmeansPlot=predictmeansPlot,
                           predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable, p_valueMatrix=p_valueMatrix)
                pm_out <- pm_finish(outputlist, tidy_parts)

      }else{

        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4),
                           predictmeansPlot=predictmeansPlot, predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable,
                           p_valueMatrix=p_valueMatrix)
                pm_out <- pm_finish(outputlist, tidy_parts)
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

            pm_out <- pm_finish(outputlist, tidy_parts)
    }
  } else {
    outputlist=list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD, predictmeansPlot=predictmeansPlot,
                    predictmeansBarPlot=predictmeansBarPlot, mean_table=meanTable)
        pm_out <- pm_finish(outputlist, tidy_parts)
  }
  #  }

})

