test_that("letters: two means share a letter exactly when they don't differ (40 means)", {
  set.seed(5)
  k <- 40
  m <- sort(rnorm(k, sd = 3))
  P <- 2*pnorm(-abs(outer(m, m, "-"))/sqrt(2))
  dimnames(P) <- list(paste0("t", 1:k), paste0("t", 1:k))
  lets <- multcompLetters(P, Letters = c(LETTERS, letters), threshold = 0.05)
  sets <- lapply(strsplit(lets, ""), function(x) setdiff(x, " "))
  share <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) length(intersect(sets[[i]], sets[[j]])) > 0))
  expect_equal(share[upper.tri(share)], (P >= 0.05)[upper.tri(P)])
})

test_that("letterCI with atvar: CI-overlap letters within each atvar group", {
  data(Oats, package = "nlme", envir = environment())
  Oats$nitro <- factor(Oats$nitro)
  fm <- nlme::lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Oats)
  r <- predictmeans(fm, "nitro:Variety", atvar = "Variety", letterCI = TRUE, plot = FALSE)
  tab <- r$mean_table$Table
  for (v in levels(Oats$Variety)) {
    g <- tab[tab$Variety == v, ]
    expect_equal(unname(g$LetterGrp), unname(ci_mcp(g[["LL(83.4%)"]], g[["UL(83.4%)"]])))
  }
  expect_match(r$mean_table$Note, "overlap of the 83.4% CIs", fixed = TRUE)
  expect_match(r$mean_table$Note, "at each level of Variety$")
  expect_equal(r$means$group[match(paste(tab$Variety, tab$nitro), paste(r$means$Variety, r$means$nitro))],
               as.character(tab$LetterGrp))
})

test_that("letterCI: a mean without a CI gets no letter instead of an error", {
  tab <- data.frame(Mean = c(1, 2, NA, 9), "LL(83.4%)" = c(0, 1, NA, 8), "UL(83.4%)" = c(2, 3, NA, 10),
                    check.names = FALSE)
  lets <- pm_ci_letters(tab)
  expect_equal(trimws(lets[3]), "")
  expect_equal(trimws(lets[-3]), trimws(unname(ci_mcp(c(0, 1, 8), c(2, 3, 10)))))
})
