## Argument checks with messages that say what to type instead.

## The p-value adjustment method, matched in any letter case ("Tukey" ->
## "tukey", "bh" -> "BH"). 'tukey = FALSE' where only p.adjust() methods work.
pm_match_adj <- function(adj, tukey = TRUE) {
  choices <- c("none", if (tukey) "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
               "sidak", "scheffe", "mvt")
  hit <- if (is.character(adj) && length(adj) == 1) choices[tolower(choices) == tolower(adj)] else character(0)
  if (length(hit) == 0) {
    stop("'adj' must be one of ", paste0('"', choices, '"', collapse = ", "), ", not ",
         paste(deparse(adj), collapse = ""), ".", call. = FALSE)
  }
  hit
}

## 'modelterm' must be a fixed-effect term of the model, written as in the
## model. For the same factors in another order or letter case, say which.
pm_check_modelterm <- function(model, modelterm) {
  labs <- attr(stats::terms(model), "term.labels")
  if (is.character(modelterm) && length(modelterm) == 1 && modelterm %in% labs) return(invisible(modelterm))
  term <- paste(modelterm, collapse = " ")
  key <- function(x) vapply(strsplit(tolower(x), ":", fixed = TRUE),
                            function(v) paste(sort(trimws(v)), collapse = ":"), "")
  hit <- labs[key(labs) == key(term)]
  hint <- if (length(hit)) {
    paste0(" Did you mean '", hit[1], "'?")
  } else {
    paste0(" The terms of the model are: ", paste0("'", labs, "'", collapse = ", "), ".")
  }
  stop("'", term, "' is not a term of the model.", hint, call. = FALSE)
}
