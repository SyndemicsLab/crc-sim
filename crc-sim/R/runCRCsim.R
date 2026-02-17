#' Function to run CRC sim - returns a list of data.tables resulting from \code{analyze}
#'
#' @param nboot int: number of bootstraps
#' @param ncores int: number of cores to use
#' @param p_captures list: list of capture probabilities
#' @param p_strata list: list of strata probabilities - for nonstratified just use '1'
#' @param seed int: starting seed
#'
#' @importFrom doFuture registerDoFuture
#' @importFrom doRNG registerDoRNG
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @import data.table
#'
#' @export

runCRC <- function(nboot, ncores, p_captures, p_strata, seed = 2024){
  Model <- NULL
  Method <- NULL

  ncores <<- ncores
  evalq(
    {
      doFuture::registerDoFuture()
      doRNG::registerDoRNG()
      future::plan(future::multisession(workers = ncores), gc = TRUE)
    },
    envir = .GlobalEnv
  )
  rm(ncores, envir = .GlobalEnv)

  seeds <- future.apply::future_lapply(1:nboot, function(x) .Random.seed, future.seed = seed)
  output <- future.apply::future_lapply(1:nboot, function(x) {
    n <- 3e5
    suppression <- 10

    config <- list(
      f0.05 = list(direction = "forward", threshold = 0.05),
      f0.1 = list(direction = "forward", threshold = 0.1),
      b0.05 = list(direction = "backward", threshold = 0.05),
      b0.1 = list(direction = "backward", threshold = 0.1),
      fb0.05 = list(direction = "both", threshold = 0.05),
      fb0.1 = list(direction = "both", threshold = 0.1)
    )

    DT <- create.data(n, p_captures, p_strata)
    gc(); gc()

    pois <- lapply(config, function(x){
      analyze(DT, suppress = suppression,
              method = "poisson", formula.selection = "stepwise", opts.stepwise = c(x, verbose = FALSE))
    })
    pois_data <- rbindlist(pois, idcol = c("Model", names(pois)))
    pois_data <- pois_data[, Model := paste0(gsub("b", "Backward-", gsub("f", "Forward-", gsub("fb", "Both-", Model))))]

    nb <- lapply(config, function(x){
      analyze(DT, suppress = suppression,
              method = "negbin", formula.selection = "stepwise", opts.stepwise = c(x, verbose = FALSE))
    })
    nb_data <- rbindlist(nb, idcol = c("Model", names(nb)))
    nb_data <- nb_data[, Model := paste0(gsub("b", "Backward-", gsub("f", "Forward-", gsub("fb", "Both-", Model))))]

    data <- rbind(nb_data[, Method := "NB"], pois_data[, Method := "Poisson"])
    gc()
    return(data)
  })
  out <- rbindlist(output, idcol = "Run")
  return(out)
}
