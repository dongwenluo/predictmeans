# Which denominator df are reliable? (2026-11-13)

The maintainer's rule: each default is the most reliable method on the
evidence. This study measures the type I error at the 5% level of every df
method of predictmeans: two-sided t tests of single contrasts (a 5% test is a
95% CI covering 95%) and F tests of terms, all true null hypotheses. Each
cell has 1000 data sets; a true 5% lies in 3.2–6.8% (99% interval).

Scripts (each writes one `.rds` per cell; `dev/sim-df-summary.R <dir>`
writes `results.md`):

| script | models | results |
|---|---|---|
| `dev/sim-df.R` | split plots (`lmer`, `lme`), repeated measures (`gls` AR(1) and unstructured, `lme`) | `dev/sim-df/results.md` |
| `dev/sim-df-kr.R` | follow-up: Kenward–Roger with the expected or the observed information when values are missing (`lmer` random slopes, `gls`) | `dev/sim-df-kr/results.md` |
| `dev/sim-df-glmm.R` | GLMMs (`glmer`): Poisson split plots, binomial clusters, binary RCBD | `dev/sim-df-glmm/results.md` |

Methods: `KR` Kenward–Roger with the expected information (pbkrtest; at the
time the `lmer` default, now `df_method = "KR-expected"`); `KRSAS` with the
observed information (SAS `DDFM=KR`; then `"KR-SAS"`, now `"Kenward-Roger"`); `Satt` Satterthwaite; `cont`
containment (SAS `DDFM=CONTAIN`); `BW` between-within; `nlme` nlme's df of the
term (the `lme` default); `nlmeC` the smallest nlme df of the coefficients a
contrast uses; `resid` the residual df (the `gls` and GLMM default); `z` no
df (normal, chi-square).

## Findings

**Split plots (`lmer`, `lme`), balanced or with 15–30% missing or a whole
plot missing.** Kenward–Roger (either information), Satterthwaite,
containment and nlme's df all hold 5% (2.9–6.6%). Between-within is
conservative for sub-plot effects (tB 3.1–4.0%, interaction F 2.6–2.9%). The
normal approximation is liberal (up to 13.3%). The current `lmer` default is
right.

**Repeated measures, unstructured covariance (`gls`), 12 or 24 subjects.** The
current `gls` default, the residual df (30), is liberal: t tests 6.2–7.3%,
F tests 15.5–15.9%; with 10% of the values missing 7.6–10.0% and
24.1–26.5%. Satterthwaite is right for t tests but its F tests (lmerTest's
rule) give 8.9–9.9% (15.4–15.7% with missing values). Kenward–Roger holds
4.2–5.7%.

**Kenward–Roger with missing values: expected or observed information.** With
unstructured repeated measures and 10% missing, the expected information
(pbkrtest) gave 8.5–8.7% for the F tests and 7.3% for a within-subject
contrast; the observed information (SAS) 5.5–6.6%. Kenward and Molenberghs
(1998) show that with missing data the expected information is not valid
(it treats the pattern of missing values as fixed); the observed information
is. The follow-up (`dev/sim-df-kr.R`) confirms it: unstructured, 24 subjects,
20% missing: F tests 6.2–6.9% (expected) against 5.6–5.7% (observed);
AR(1) with unequal variances, 20% missing: 4.5–5.7% against 4.4–5.3%. In
`lmer` models the two differ little (split plot with 30% missing: 3.2–5.3%
against 3.4–5.3%; random slopes with 20% missing: F tests 5.4–5.6% against
4.5–5.0%; random intercept with 20% missing: the same). Nowhere was the
observed information worse. Satterthwaite F tests are liberal whenever the
covariance has several parameters (AR(1) with unequal variances 6.6–8.4%,
unstructured 9.7–11.9%); between-within is liberal there too (up to 12.8%).

**GLMMs (`glmer`).** The current df, the residual df (n − parameters), are
liberal for effects that vary between clusters: whole-plot factor of a
Poisson split plot 9.0–10.7% (t) and 12.1–15.2% (F); a treatment applied to
whole clusters 8.4–12.8%. Containment (SAS GLIMMIX's default) helps for the
split plot (7.3–8.1%) but not for the cluster treatment (12.7%): SAS's rule
looks at the names in the model, and `(1 | clus)` does not name the
treatment. nlme's df by coefficient (`nlmeC`, from the design: the level of
the random effects at which each coefficient varies) is the best overall:
whole plot 6.4–8.1%, cluster treatment 6.4–7.7%, treatment within clusters
4.1–5.0%, within-cluster effects 3.6–4.5%. Between-within is as good between
clusters but conservative within (1.5–2.9%). What remains liberal (6–8% with
4–6 blocks or 4 clusters per group) is the Wald SE of a GLMM with few
clusters; `bootmodels()` or `permmodels()` are the checks there.

## Recommendations

1. `lme` and `gls` fitted by REML: Kenward–Roger by default (as `lmer`), with
   adjusted SEs; nlme's df and the residual df stay available as
   `df_method = "containment"`. The nlme Kenward–Roger code now works subject
   by subject (sparse block-diagonal covariance): an unstructured `gls` with
   2000 observations takes 3 s (was about 70 s).
2. One Kenward–Roger: the observed information (as SAS) for every model,
   because it stays valid with missing values; the expected information
   (pbkrtest, emmeans) as an option.
3. GLMMs: the df of the level at which each contrast varies (nlme's rule by
   coefficient), computed from the design, instead of the residual df —
   **done** (`df_method = "inner-outer"`, `R/df-innerouter.R`). The rule is
   computed from ranks, which equals nlme for nested random effects and
   extends it to crossed ones; checked by simulation
   (`dev/sim-df-glmm/results.md`, `gl_pois_crossed`, `gl_binom_crossed`): a
   treatment on the columns of a rows x columns design 7.0-7.2% against
   9.8-13.0% with the residual df or containment; effects within cells
   3.7-5.4% (both).

Decisions 1 and 2 were also made (2026-11-13): Kenward–Roger with the
observed information is the default for REML `lmer`, `lme` and `gls`.
4. Satterthwaite F tests of several df are liberal with unstructured
   covariances; say so in the help and prefer Kenward–Roger for F tests.
