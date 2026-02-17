#' Function to compare between simulation results and ground truth
#'
#' @param sim list of simulations done by \code{simulate}
#' @param groundTruth list of ground truths by \code{extract.groundTruth}
#'
#' @import data.table
#' @keywords internal

compare <- function(sim, groundTruth) {
    compare.list <- c()
    if (length(groundTruth) > 1) {
        for (i in seq_along(sim)) {
            aic <- sim[[i]]$AIC
            est <- sim[[i]]$estimate
            est.lci <- sim[[i]]$lower_ci
            est.uci <- sim[[i]]$upper_ci

            group <- sort(names(groundTruth))[i]
            gt <- as.integer(groundTruth[[group]])

            out.list <- list(
                aic = aic,
                estimate = est,
                lower_ci = est.lci,
                upper_ci = est.uci,
                ground = gt,
                group = group
            )
            compare.list[[i]] <- out.list
        }
        out <- rbindlist(compare.list)
    } else {
        aic <- sim$AIC
        est <- sim$estimate
        est.lci <- sim$lower_ci
        est.uci <- sim$upper_ci
        group <- "base"
        gt <- as.integer(groundTruth[1])

        out <- data.table(
            aic = aic,
            estimate = est,
            lower_ci = est.lci,
            upper_ci = est.uci,
            ground = gt,
            group = group
        )
    }
    return(out)
}
