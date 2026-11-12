## permindex() and the 'perms' argument of permmodels() / permlmer()

test_that("permindex(): same permutations as 1.x for sorted data, right ones for unsorted data", {
  d <- data.frame(block = rep(1:3, each = 6), group = rep(rep(1:2, each = 3), 3), y = 1:18)
  # the 1.x algorithm (data sorted by block)
  Bsplit <- split(seq_len(18), as.character(d$block))
  set.seed(7 + 1)
  old <- replicate(20, unlist(lapply(Bsplit, function(x) x[sample.int(length(x))]), use.names = FALSE))
  old <- old[, !duplicated(apply(old, 2, paste, collapse = ",")), drop = FALSE]
  expect_identical(permindex(d, block = "block", nsim = 20, seed = 7), old)
  # unsorted data: rows take values only from their own block (1.x indices
  # referred to the sorted data, so here they would mix blocks)
  o <- sample(18)
  du <- d[o, ]
  iu <- permindex(du, block = "block", nsim = 20, seed = 7)
  expect_true(all(apply(iu, 2, function(p) all(du$block[p] == du$block))))
  # whole groups within blocks: each group's rows come from one group of the same block
  ig <- permindex(du, block = "block", group = "group", nsim = 10, seed = 3)
  expect_true(all(apply(ig, 2, function(p) {
    src <- interaction(du$block, du$group)[p]
    all(du$block[p] == du$block) && all(tapply(src, interaction(du$block, du$group), function(s) length(unique(s))) == 1)
  })))
  # groups of unequal size cannot change places
  expect_error(permindex(d[-1, ], block = "block", group = "group", nsim = 5), "same number of rows")
  # duplicates dropped: a block of two rows has two permutations
  expect_identical(ncol(permindex(data.frame(b = c(1, 1)), block = "b", nsim = 50, seed = 1)), 2L)
  # the seed leaves the user's stream alone
  set.seed(9); a <- runif(1)
  set.seed(9); invisible(permindex(d, nsim = 5, seed = 1)); b <- runif(1)
  expect_identical(a, b)
})

test_that("permmodels() and permlmer() use the permutations given in 'perms'", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro, data = Oats)
  idx <- permindex(Oats, block = "Block", nsim = 39, seed = 2)
  r <- permmodels(fm, perms = idx, prt = FALSE, ncore = 1)
  Fobs <- anova(fm)["nitro", "F value"]
  Fs <- apply(idx, 2, function(p) {
    dd <- Oats; dd$yield <- Oats$yield[p]
    anova(lm(yield ~ Block + nitro, data = dd))["nitro", "F value"]
  })
  expect_equal(as.numeric(r$ANOVA["nitro", "Perm_p_value"]),
               round((sum(round(Fs, 6) >= round(Fobs, 6)) + 1) / (ncol(idx) + 1), 4))
  # Freedman-Lane with the same permutations, and its pairs
  rf <- permmodels(fm, perms = idx, prt = FALSE, method = "freedman-lane")
  expect_identical(rf$perm$nperm, ncol(idx))
  pm <- predictmeans(fm, "nitro", permlist = rf, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_equal(pm$pairs$p.value[1], unname(rf$COEFFICENTS["nitro0.2", "Perm_p_value"]), tolerance = 1e-4)
  expect_error(permmodels(fm, perms = idx[-1, ], prt = FALSE, ncore = 1), "one row per observation")
  bad <- idx; bad[1, 1] <- bad[2, 1]
  expect_error(permmodels(fm, perms = bad, prt = FALSE, ncore = 1), "permutation of 1 to")

  # lme: the rows are reordered inside; results do not depend on the order of the data
  skip_on_cran()
  o <- sample(nrow(Oats))
  Os <- Oats[o, ]
  m1 <- lme(yield ~ nitro, random = ~ 1 | Block, data = Oats)
  m2 <- lme(yield ~ nitro, random = ~ 1 | Block, data = Os)
  i1 <- permindex(Oats, block = "Block", nsim = 19, seed = 4)
  i2 <- i1[o, ]; i2[] <- match(i1[o, ], o)        # the same permutations, for the shuffled rows
  r1 <- permmodels(m1, perms = i1, prt = FALSE, ncore = 1)
  r2 <- permmodels(m2, perms = i2, prt = FALSE, ncore = 1)
  expect_identical(r1$ANOVA$Perm_p_value, r2$ANOVA$Perm_p_value)
  # permlmer
  l0 <- lme4::lmer(yield ~ nitro + (1 | Block), data = Oats)
  l1 <- lme4::lmer(yield ~ nitro + (1 | Block/Variety), data = Oats)
  p <- permlmer(l0, l1, perms = i1, ncore = 1)
  expect_true(p$`Perm-p`[2] > 0 && p$`Perm-p`[2] <= 1)
})
