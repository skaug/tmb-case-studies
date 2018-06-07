library(TMB)
#' This function writes out files such that they can be illustrated in rmarkdown.
#' @param x The file name
#' @param linesToInclud A vector giving the line numbers which are to be written
include_source <- function(x,linesToInclud=seq(1:9999)) {
    ans <- readLines(x)
    linesToInclud = linesToInclud[1:min(length(linesToInclud), length(ans))]
    ans = ans[linesToInclud]
    ans = ans[!is.na(ans)]
    return(cat(paste(ans, collapse="\n")))
}

