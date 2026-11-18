## Runs every 2.0 line of the side-by-side code in vignette("migrating"), and checks
## that the 1.1.1 lines stop with the 2.0 spelling. Rscript dev/check-migrating.R
pkgload::load_all(".", quiet = TRUE); pdf(NULL)
Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
fm0 <- lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
fmc <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = nlme::Oats)
fs <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject), data = nlme::Dialyzer)
ck <- function(lab, e) cat(lab, ":", tryCatch({force(e); "ok"}, error = function(err) paste("ERROR", conditionMessage(err))), "\n")
ck("pm", pm <- predictmeans(fm, "nitro:Variety", pairwise = TRUE, conf.level = 0.9))
ck("autoplot pm", print(autoplot(pm, plotord = c(2, 1), plottitle = "Yield")))
ck("plot_args", predictmeans(fm, "nitro:Variety", plot = TRUE, plot_args = list(plotord = c(2, 1), plottitle = "Yield")))
ck("print digits", capture.output(print(pm, digits = 3)))
ck("bar", predictmeans(fm, "nitro", df_method = 30, plot = "bar"))
K <- rbind(c(1, -1, 0, 0)); ck("contrast", contrastmeans(fm, "nitro", ctrmatrix = K, df_method = 30))
ck("df_term", df_term(fm, "nitro", df_method = "Satterthwaite"))
ck("cm", cm <- covariatemeans(fmc, "Variety", covariate = "nitro", at = seq(0, 0.6, 0.1), conf.level = 0.9))
ck("autoplot cm", print(autoplot(cm, trellis = FALSE, mtitle = "Yield by nitrogen")))
ck("sp", sp <- semipred(fs, "QB", "pressure", at = c(0.5, 1, 1.5)))
ck("autoplot sp", print(autoplot(sp, trellis = FALSE)))
ck("K", Kmatrix(fmc, "Variety", covariate = "nitro", at = 0.3))
ck("permmodels", permmodels(fm, nperm = 49, type = 3))
ck("permlmer", permlmer(fm0, fm, nperm = 19))
ck("residplot", residplot(fm, re_level = 2, id = 3))
ck("CookD", CookD(fm, id = 5))
ck("PMplot", PMplot(pm$p_valueMatrix, conf.level = 0.99))
ck("varcomp", varcomp(fm, conf.level = 0.9))
ck("R2", R2_model(fm))
ck("letters", print(predictmeans(fm, "nitro", pairwise = TRUE, letterCI = TRUE)$mean_table))
# 1.1.1 side must stop with the right message
for (e in list(quote(predictmeans(fm, "nitro", mplot = FALSE, barplot = TRUE)),
               quote(covariatemeans(fmc, "Variety", covariate = "nitro", trellis = FALSE)),
               quote(semipred(semireg = fs, modelterm = "QB", covariate = "pressure")),
               quote(residplot(fm, level = 2)), quote(PMplot(pm$p_valueMatrix, level = 0.01)),
               quote(permanova.lmer(fm)), quote(R2_glmm(fm)), quote(Kmatrix(fmc, "Variety", covariate = "nitro", prtnum = TRUE))))
  cat(deparse(e)[1], "->", tryCatch({eval(e); "NO ERROR"}, error = function(err) conditionMessage(err)), "\n")
