#' Attempts recovery of the ground truth - fundamentally a wrapper for Syndemics::crc
#'
#' @importFrom Syndemics crc
#' @param DT data.table from the \code{create.data} step
#' @param capture list: strings of captures
#' @param group stratification string to group by
#' @param suppress numeric: maximum value to suppress at
#' @param method string: selection for spatial capture recapture model - either "poisson" or "negbin"
#' @param formula.selection string: selection for formula decision - either "aic", "corr", or "stepwise"
#' @param opts.stepwise list: list containing 'direction' of 'forward' 'backward' or 'both', and 'threshold': p value threshold for stepwise selection
#'
#' @importFrom Syndemics crc
#' @keywords internal
simulate <- function(DT, capture = c("APCD", "BSAS", "Casemix", "Death", "Matris", "PMP"), group, suppress,
                     method, formula.selection, opts.stepwise){
  N_ID <- NULL
  tmp <- NULL

  DT <- DT[, tmp := rowSums(.SD), .SDcols = capture
           ][tmp != 0,
             ][, tmp := NULL]

  if(suppress) DT <- DT[!N_ID %in% 1:suppress, ]

  if(!missing(group)) {
    DT.list <- c()
    n_groups <- sort(unique(DT[[group]]))

    for(i in seq_along(n_groups)){
      DT.list[[i]] <- DT[get(group) == n_groups[i],
                         ][, paste0(group) := NULL]
    }
    out.list <- lapply(DT.list, function(x) Syndemics::crc(x, "N_ID", binary.variables = capture, method, formula.selection, opts.stepwise = opts.stepwise))
    return(out.list)
  } else {
    out <- Syndemics::crc(DT, "N_ID", binary.variables = capture, method, formula.selection, opts.stepwise = opts.stepwise)
    return(out)
  }
}
