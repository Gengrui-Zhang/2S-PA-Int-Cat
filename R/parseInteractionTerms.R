# Interaction-branch

#' Unconstrained product indicator method
#'
#' @param model A string variable describing the structural path model with interaction term(s),
#'              in \code{lavaan} syntax.
#' @return A list of interaction pair
#' @export

parseInteractionTerms <- function(model) {
  lines <- unlist(strsplit(model, "\n"))
  lines <- lines[!grepl(":=", lines)]
  inter_lines <- grep(":", lines, value = TRUE)
  inter_lines <- inter_lines[grepl("~", inter_lines)]
  inter_terms <- c()
  for (line in inter_lines) {
    clean_line <- gsub(" ", "", line)
    rhs <- unlist(strsplit(clean_line, "~"))[2]
    rhs_terms <- unlist(strsplit(rhs, "\\+"))
    inters <- rhs_terms[grepl(":", rhs_terms)]
    inters <- gsub("^[^*]*\\*", "", inters)
    inter_terms <- c(inter_terms, inters)
  }
  inter_terms <- unique(inter_terms)
  interpairs <- function (inter_terms) {
    inter_vars <- list()
    for (i in seq_along(inter_terms)) {
      vars <- unlist(strsplit(inter_terms[i], ":"))
      inter_vars[[i]] <- vars
      names(inter_vars)[i] <- paste0("inter_pair_", i)
    }
    return(inter_vars)
  }
  pairs <- interpairs(inter_terms)
  return(pairs)
}
