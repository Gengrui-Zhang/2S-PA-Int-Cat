# Interaction-branch

#' Unconstrained product indicator method
#'
#' @param model A string variable describing the structural path model with interaction term(s),
#'              in \code{lavaan} syntax.
#' @return A list of indicators of interaction pair
#' @export

parseIndicators <- function(model) {
  # Step 1: Remove comments and trim whitespace
  lines <- unlist(strsplit(model, "\n"))
  lines <- gsub("#.*", "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]  # remove empty lines
  
  # Step 2: Accumulate =~ lines, supporting multiline
  measurement_lines <- c()
  current <- ""
  for (line in lines) {
    if (grepl("=~", line)) {
      if (nzchar(current)) measurement_lines <- c(measurement_lines, current)
      current <- line
    } else if (grepl("~|~~|:=", line)) {
      if (nzchar(current)) {
        measurement_lines <- c(measurement_lines, current)
        current <- ""
      }
    } else {
      current <- paste(current, line)
    }
  }
  if (nzchar(current)) measurement_lines <- c(measurement_lines, current)
  
  # Step 3: Extract indicators
  indicators <- list()
  for (line in measurement_lines) {
    parts <- strsplit(line, "=~")[[1]]
    if (length(parts) != 2) next
    latent <- trimws(parts[1])
    rhs <- trimws(parts[2])
    
    rhs_items <- trimws(unlist(strsplit(rhs, "\\+")))
    
    # Strip loadings on either side: 1*x1 or x1*1 → x1
    rhs_items <- gsub("^.*\\*", "", rhs_items)  # remove prefix if "1*x1"
    rhs_items <- gsub("\\*.*$", "", rhs_items)  # remove suffix if "x1*1"
    
    indicators[[latent]] <- rhs_items
  }
  
  return(indicators)
}

