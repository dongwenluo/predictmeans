## smZ() was removed in 2.0: smooths are written in the formula of semireg()
## and semireg_tmb() (vignette("migrating") shows each smZ() form beside its
## formula form; dev/check-smz-migration.R checked that they give the same
## fits).
smZ <- function(...) pm_defunct_fun("smZ")
