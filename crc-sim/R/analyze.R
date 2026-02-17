#' Function to run the analysis
#'
#' @param DT data.table: data.table created from \code{create.data}
#' @param suppress int: supression cap, defaults to 10
#' @param method string: selection for spatial capture recapture model - either "poisson" or "negbin"
#' @param formula.selection string: selection for formula decision - either "aic", "corr", or "stepwise"
#' @param opts.stepwise list: list containing 'direction' of 'forward' 'backward' or 'both', and 'threshold': p value threshold for stepwise selection
#'
#' @import data.table
#' @keywords internal

analyze <- function(DT, suppress,
                    method = "poisson", formula.selection = "stepwise", opts.stepwise = list(direction = "forward",
                                                                                             threshold = 0.5,
                                                                                             verbose = FALSE)){
  captures <- names(DT)[!names(DT) %in% c("N_ID", "strata", "tmp")]

  groundTruth <- extract.groundTruth(DT, capture = captures, group = "strata")
  model <- simulate(DT, capture = captures, group = "strata", suppress = suppress,
                    method, formula.selection, opts.stepwise)

  out <- compare(model, groundTruth)
  return(out)
}
