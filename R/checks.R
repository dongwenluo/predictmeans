## Argument checks with messages that say what to type instead.

## The p-value adjustment method, matched in any letter case ("Tukey" ->
## "tukey", "bh" -> "BH"). 'tukey = FALSE' where only p.adjust() methods work.
pm_match_adj <- function(adj, tukey = TRUE, westfall = FALSE) {
  # "westfall" (step-down permutation max-|t|) where a permlist can be given
  choices <- c("none", if (tukey) "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
               "sidak", "scheffe", "mvt", if (westfall) "westfall")
  hit <- if (is.character(adj) && length(adj) == 1) choices[tolower(choices) == tolower(adj)] else character(0)
  if (length(hit) == 0) {
    stop("'adj' must be one of ", paste0('"', choices, '"', collapse = ", "), ", not ",
         paste(deparse(adj), collapse = ""), ".", call. = FALSE)
  }
  hit
}

## 'modelterm' is one string: fixed-effect variables of the model joined by
## ":" ("A", "A:B"). It need not be a term of the model: any combination of
## its factors gives means (Kmatrix() says which variables exist).
pm_check_modelterm <- function(model, modelterm) {
  if (!is.character(modelterm) || length(modelterm) != 1 || is.na(modelterm) || !nzchar(modelterm)) {
    stop("'modelterm' must be one character string such as \"A\" or \"A:B\".", call. = FALSE)
  }
  invisible(modelterm)
}

## A combination of factors that is not a term of the model: its means are
## what the model implies for it, the sum of the effects of its parts.
pm_note_not_a_term <- function(model, modelterm) {
  labs <- tryCatch(attr(stats::terms(model), "term.labels"), error = function(e) NULL)
  if (is.null(labs)) return(invisible(NULL))
  key <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(v) paste(sort(pm_term_vars(v)), collapse = ":"), "")
  if (key(modelterm) %in% key(labs)) return(invisible(NULL))
  message("'", modelterm, "' is not a term of the model: its factors do not interact there, so these means are the ",
          "sums of their effects (the differences between the levels of one factor are the same at every level ",
          "of the others).")
  invisible(NULL)
}

## The df of a term from a set of df by term ('one_df', named by term label):
## the term's own df; for a combination of variables that is not a term of
## the model (A:B in y ~ A + B), the smallest df of the model's terms made of
## its variables. NULL when no term of the model is made of them.
pm_term_df <- function(one_df, term) {
  vars_of <- function(x) strsplit(x, ":", fixed = TRUE)
  d <- one_df[intersect(term, names(one_df))]
  if (length(d)) return(min(unname(d)))
  tv <- unique(pm_term_vars(unlist(vars_of(term))))
  labs <- setdiff(names(one_df), "(Intercept)")
  inside <- vapply(vars_of(labs), function(v) all(pm_term_vars(v) %in% tv), logical(1))
  if (!any(inside)) return(NULL)
  min(unname(one_df[labs[inside]]))
}
