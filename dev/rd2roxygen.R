# Convert man/*.Rd to roxygen blocks (one R/doc-<name>.R per page).
# Used once for the roxygen2 migration; the blocks have since been moved next
# to their functions (data sets in R/data.R), so don't run it again.
pkg <- commandArgs(TRUE)[1]
bs <- "\\"
rdtext <- function(x) {
  tag <- attr(x, "Rd_tag")
  if (is.null(tag) || tag == "LIST") return(paste0("{", kids(x), "}"))
  if (tag %in% c("TEXT", "RCODE", "VERB", "UNKNOWN", "USERMACRO")) return(gsub("%", "\\%", paste(unlist(x), collapse = ""), fixed = TRUE))
  if (tag == "COMMENT") return("")
  if (substr(tag, 1, 1) == "#") return(paste0(tag, " ", kids(x), "#endif\n"))
  if (substr(tag, 1, 1) == bs) {
    opt <- attr(x, "Rd_option")
    if (!is.null(opt)) tag <- paste0(tag, "[", kids(list(opt)), "]")
    if (!length(x)) return(tag)
    islist <- vapply(x, function(z) is.list(z) && (is.null(attr(z, "Rd_tag")) || attr(z, "Rd_tag") == "LIST"), TRUE)
    if (all(islist)) return(paste0(tag, paste0("{", vapply(x, kids, ""), "}", collapse = "")))
    return(paste0(tag, "{", kids(x), "}"))
  }
  kids(x)
}
kids <- function(x) paste(vapply(x, rdtext, ""), collapse = "")
trim <- function(s) gsub("^\\s+|\\s+$", "", s)
rox <- function(s) { l <- strsplit(s, "\n", fixed = TRUE)[[1]]; paste0("#'", ifelse(nzchar(l), " ", ""), l) }
map <- c(value = "return", note = "note", references = "references", author = "author", seealso = "seealso",
         format = "format", source = "source")
for (f in list.files(file.path(pkg, "man"), "\\.Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(f)
  tags <- sub("^\\\\", "", vapply(rd, function(x) attr(x, "Rd_tag"), ""))
  get <- function(tag) rd[tags == tag]
  name <- trim(kids(get("name")[[1]]))
  out <- character(0); add <- function(...) out <<- c(out, ...)
  add(rox(trim(kids(get("title")[[1]]))), "#'")
  for (d in get("description")) add("#' @description", rox(trim(kids(d))))
  for (d in get("details")) add("#' @details", rox(trim(kids(d))))
  u <- get("usage"); if (length(u)) add("#' @usage", rox(trim(kids(u[[1]]))))
  for (a in get("arguments")) for (it in a) if (identical(attr(it, "Rd_tag"), "\\item"))
    add(rox(paste0("@param ", trim(kids(it[[1]])), " ", trim(kids(it[[2]])))))
  for (tg in names(map)) for (v in get(tg)) add(paste0("#' @", map[[tg]]), rox(trim(kids(v))))
  for (s in get("section")) add(paste0("#' @section ", trim(kids(s[[1]])), ":"), rox(trim(kids(s[[2]]))))
  ex <- get("examples"); if (length(ex)) add("#' @examples", rox(sub("^\n", "", sub("\\s+$", "", kids(ex[[1]])))))
  al <- vapply(get("alias"), function(a) trim(kids(a)), "")
  add(paste0("#' @name ", name))
  if (length(setdiff(al, name))) add(paste("#' @aliases", paste(setdiff(al, name), collapse = " ")))
  is_pkg <- any(vapply(get("docType"), function(d) trim(kids(d)) == "package", TRUE))
  for (d in get("docType")) if (trim(kids(d)) != "package") add(paste("#' @docType", trim(kids(d))))
  for (k in get("keyword")) add(paste("#' @keywords", trim(kids(k))))
  for (k in get("concept")) add(paste("#' @concept", trim(kids(k))))
  other <- setdiff(unique(tags), c("name","title","description","details","usage","arguments",names(map),"section","examples","alias","docType","keyword","concept","TEXT","COMMENT","encoding"))
  if (length(other)) message(name, ": unhandled ", paste(other, collapse = " "))
  add(if (is_pkg) '"_PACKAGE"' else "NULL")
  writeLines(c(paste0("# Help page ", name, " (roxygen2 source; man/", basename(f), " is generated)."), out),
             file.path(pkg, "R", paste0("doc-", name, ".R")))
}
