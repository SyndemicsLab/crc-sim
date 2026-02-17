#' Function to extract the ground truth from data created through the \code{create.data} method
#'
#' @param DT data.table from \code{create.data}
#' @param capture list: strings of captures
#' @param group character of strata column to extract on
#'
#' @import data.table
#' @returns list the same length of \code{group}
#'
#' @keywords internal

extract.groundTruth <- function(DT,
                                capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP"),
                                group){
  tmp <- NULL
  DT <- DT[, tmp := rowSums(.SD), .SDcols = capture
           ][tmp == 0,
             ][, tmp := NULL]


  out <- as.list(DT[["N_ID"]])
  if(!missing(group)) {
    names <- as.list(DT[[group]])
    names(out) <- names
  } else names(out) <- "base"

  return(out)
}
