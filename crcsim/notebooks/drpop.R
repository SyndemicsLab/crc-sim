################################################################################
# File: drpop.r                                                                #
# Project: crcsim                                                              #
# Created Date: 2026-09-02                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-09-02                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

popsize <- function(
    data,
    K = 2,
    j,
    k,
    margin = 0.005,
    filterrows = FALSE,
    nfolds = 5,
    funcname = c("rangerlogit"),
    sl.lib = c(
        "SL.gam",
        "SL.glm",
        "SL.glm.interaction",
        "SL.ranger",
        "SL.glmnet"
    ),
    getnuis,
    q1mat,
    q2mat,
    q12mat,
    idfold,
    TMLE = TRUE,
    PLUGIN = TRUE,
    Nmin = 100,
    ...
) {
    if (!missing(j) & !missing(k)) {
        if (j == k) {
            k <- j %% K + 1
            warning(paste0("Selected lists are identical. Using k = ", k, "."))
        }
        if (j > k) {
            j <- j + k
            k <- j - k
            j <- j - k
            warning("Switching j and k to ensure j < k.")
        }
    }
    if (missing(getnuis) & missing(q1mat) & missing(q2mat) & missing(q12mat)) {
        if (!missing(j) & missing(k)) {
            if (j == K) {
                return(popsize_base(
                    data,
                    K = K,
                    k0 = j,
                    filterrows = filterrows,
                    funcname = funcname,
                    nfolds = nfolds,
                    margin = margin,
                    sl.lib = sl.lib,
                    Nmin = Nmin,
                    TMLE = TMLE,
                    PLUGIN = PLUGIN,
                    ...
                ))
            } else {
                return(popsize_base(
                    data,
                    K = K,
                    j0 = j,
                    filterrows = filterrows,
                    funcname = funcname,
                    nfolds = nfolds,
                    margin = margin,
                    sl.lib = sl.lib,
                    Nmin = Nmin,
                    TMLE = TMLE,
                    PLUGIN = PLUGIN,
                    ...
                ))
            }
        } else if (missing(j) & !missing(k)) {
            if (k < K) {
                return(popsize_base(
                    data,
                    K = K,
                    j0 = k,
                    filterrows = filterrows,
                    funcname = funcname,
                    nfolds = nfolds,
                    margin = margin,
                    sl.lib = sl.lib,
                    Nmin = Nmin,
                    TMLE = TMLE,
                    PLUGIN = PLUGIN,
                    ...
                ))
            } else {
                return(popsize_base(
                    data,
                    K = K,
                    k0 = k,
                    filterrows = filterrows,
                    funcname = funcname,
                    nfolds = nfolds,
                    margin = margin,
                    sl.lib = sl.lib,
                    Nmin = Nmin,
                    TMLE = TMLE,
                    PLUGIN = PLUGIN,
                    ...
                ))
            }
        } else {
            return(popsize_base(
                data,
                K = K,
                j0 = j,
                k0 = k,
                filterrows = filterrows,
                funcname = funcname,
                nfolds = nfolds,
                margin = margin,
                sl.lib = sl.lib,
                Nmin = Nmin,
                TMLE = TMLE,
                PLUGIN = PLUGIN,
                ...
            ))
        }
    }
    K <- 2
    n <- nrow(data)
    if (missing(j)) {
        j <- 1
    }
    if (missing(k)) {
        k <- 2
    }
    if (!missing(getnuis)) {
        if (class(getnuis) == "data.frame") {
            q1mat <- subset(
                getnuis,
                select = grep(colnames(getnuis), pattern = "q1$")
            )
            colnames(q1mat) <- stringr::str_remove_all(colnames(q1mat), "\\.q1")
            q2mat <- subset(
                getnuis,
                select = grep(colnames(getnuis), pattern = "q2$")
            )
            colnames(q2mat) <- stringr::str_remove_all(colnames(q1mat), "\\.q2")
            q12mat <- subset(
                getnuis,
                select = grep(colnames(getnuis), pattern = "q12$")
            )
            colnames(q12mat) <- stringr::str_remove_all(
                colnames(q1mat),
                "\\.q12"
            )
        } else {
            q1mat <- getnuis$q1mat
            q2mat <- getnuis$q2mat
            q12mat <- getnuis$q12mat
            idfold <- getnuis$idfold
        }
    }
    stopifnot(!is.null(q1mat) & !is.null(q2mat) & !is.null(q12mat))
    funcname <- colnames(q12mat)
    if (missing(idfold) | is.null(idfold)) {
        idfold <- rep(1, n)
    }
    nfolds <- max(idfold)
    stopifnot(!is.null(dim(data)))
    if (!informat(data = data, K = K)) {
        data <- reformat(data = data, capturelists = 1:K)
    }
    data <- as.data.frame(data)
    N <- nrow(data)
    stopifnot(N > 1)
    colnames(data) <- c(
        paste("L", 1:K, sep = ""),
        paste("x", 1:(ncol(data) - K), sep = "")
    )
    psiinv_summary <- matrix(
        0,
        nrow = K * (K - 1) / 2,
        ncol = 3 *
            length(funcname)
    )
    rownames(psiinv_summary) <- paste0(j, ",", k)
    colnames(psiinv_summary) <- paste(
        rep(funcname, each = 3),
        c("PI", "DR", "TMLE"),
        sep = "."
    )
    var_summary <- psiinv_summary
    ifvals <- matrix(
        NA,
        nrow = N * K * (K - 1) / 2,
        ncol = length(funcname) +
            1
    )
    colnames(ifvals) <- c("listpair", funcname)
    ifvals[, "listpair"] <- rep(rownames(psiinv_summary), each = N)
    nuis <- matrix(
        NA,
        nrow = N * K * (K - 1) / 2,
        ncol = 3 * length(funcname) + 1
    )
    colnames(nuis) <- c(
        "listpair",
        paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = ".")
    )
    nuis <- as.data.frame(nuis)
    sapply(nuis, "class")
    nuis[, "listpair"] <- ifvals[, "listpair"]
    if (TMLE) {
        nuistmle <- nuis
    }
    psiinvmat <- matrix(NA, nrow = nfolds, ncol = 3 * length(funcname))
    colnames(psiinvmat) <- paste(
        rep(funcname, each = 3),
        c("PI", "DR", "TMLE"),
        sep = "."
    )
    varmat <- psiinvmat
    for (folds in 1:nfolds) {
        List2 <- data[idfold == folds, ]
        yj <- List2[, paste("L", j, sep = "")]
        yk <- List2[, paste("L", k, sep = "")]
        for (func in funcname) {
            q12 <- q12mat[idfold == folds, func]
            q1 <- pmin(pmax(q12, q1mat[idfold == folds, func]), 1)
            q2 <- pmax(
                q12 / q1,
                pmin(q2mat[idfold == folds, func], 1 + q12 - q1, 1)
            )
            nuis[
                idfold == folds,
                paste(func, c("q12", "q1", "q2"), sep = ".")
            ] <- cbind(q12, q1, q2)
            gammainvhat <- q1 * q2 / q12
            psiinvhat <- mean(gammainvhat, na.rm = TRUE)
            phihat <- gammainvhat *
                (yk / q2 + yj / q1 - yj * yk / q12) -
                psiinvhat
            ifvals[idfold == folds, func] <- phihat
            Qnphihat <- mean(phihat, na.rm = TRUE)
            psiinvhat.dr <- max(psiinvhat + Qnphihat, 1)
            psiinvmat[folds, paste(func, c("PI", "DR"), sep = ".")] <- c(
                psiinvhat,
                psiinvhat.dr
            )
            sigmasq <- var(phihat, na.rm = TRUE)
            varmat[folds, paste(func, c("PI", "DR"), sep = ".")] <- sigmasq / N
            datmat <- as.data.frame(cbind(
                yj,
                yk,
                yj * yk,
                q1 -
                    q12,
                q2 - q12,
                q12
            ))
            datmat[, 4:6] <- cbind(apply(datmat[, 4:6], 2, function(u) {
                return(pmin(pmax(u, margin), 1 - margin))
            }))
            colnames(datmat) <- c("yj", "yk", "yjk", "q10", "q02", "q12")
            if (TMLE) {
                tmle <- tmle(datmat = datmat, margin = margin, K = 2, ...)
            } else {
                tmle <- list(error = TRUE)
            }
            if (tmle$error) {
                warning("TMLE did not run or converge.")
                psiinvmat[folds, paste(func, "TMLE", sep = ".")] <- NA
                varmat[folds, paste(func, "TMLE", sep = ".")] <- NA
            } else {
                datmat <- tmle$datmat
                q12 <- pmax(datmat$q12, margin)
                q1 <- pmin(datmat$q12 + datmat$q10, 1)
                q2 <- pmax(
                    pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1),
                    q12 / q1
                )
                nuistmle[
                    idfold == folds,
                    paste(func, c("q12", "q1", "q2"), sep = ".")
                ] <- cbind(q12, q1, q2)
                gammainvhat <- q1 * q2 / q12
                psiinvhat.tmle <- mean(gammainvhat, na.rm = TRUE)
                phihat <- gammainvhat *
                    (yj /
                        q1 +
                        yk / q2 -
                        yj *
                            yk /
                            q12) -
                    psiinvhat.tmle
                Qnphihat <- mean(phihat, na.rm = TRUE)
                psiinvmat[
                    folds,
                    paste(func, "TMLE", sep = ".")
                ] <- psiinvhat.tmle
                sigmasq <- var(phihat, na.rm = TRUE)
                varmat[folds, paste(func, "TMLE", sep = ".")] <- sigmasq / N
            }
        }
    }
    psiinv_summary[paste0(j, ",", k), ] <- colMeans(psiinvmat, na.rm = TRUE)
    var_summary[paste0(j, ",", k), ] <- colMeans(varmat, na.rm = TRUE)
    result <- list(
        psi = 1 / psiinv_summary,
        sigma = sqrt(N * var_summary),
        n = round(N * psiinv_summary),
        sigman = sqrt(
            N^2 * var_summary + N * psiinv_summary * (psiinv_summary - 1)
        ),
        cin.l = round(pmax(
            N *
                psiinv_summary -
                1.96 *
                    sqrt(
                        N^2 *
                            var_summary +
                            N * psiinv_summary * (psiinv_summary - 1)
                    ),
            N
        )),
        cin.u = round(
            N *
                psiinv_summary +
                1.96 *
                    sqrt(
                        N^2 *
                            var_summary +
                            N *
                                psiinv_summary *
                                (psiinv_summary -
                                    1)
                    )
        )
    )
    result <- Reduce(
        function(...) merge(..., by = c("listpair", "Var2")),
        lapply(1:length(result), function(i) {
            reshape2::melt(
                result[[i]],
                value.name = names(result)[i],
                varnames = c("listpair", "Var2")
            )
        })
    )
    result <- tidyr::separate(
        data = result,
        col = "Var2",
        into = c("model", "method"),
        sep = "\\."
    )
    if (!TMLE) {
        result <- result[result$method != "TMLE", ]
    }
    if (!PLUGIN) {
        result <- result[result$method != "PI", ]
    } else {
        warning(
            "Plug-in variance is not well-defined. Returning variance evaluated using DR estimator formula"
        )
    }
    ifvals <- as.data.frame(ifvals)
    ifvals$listpair <- paste0(j, ",", k)
    nuis <- as.data.frame(nuis)
    nuis$listpair <- paste0(j, ",", k)
    object <- list(
        result = result,
        N = N,
        ifvals = as.data.frame(ifvals),
        nuis = as.data.frame(nuis),
        idfold = idfold
    )
    if (TMLE) {
        nuistmle <- as.data.frame(nuistmle)
        nuistmle$listpair <- paste0(j, ",", k)
        object$nuistmle <- as.data.frame(nuistmle)
    }
    class(object) <- "popsize"
    return(invisible(object))
}

popsize_base <- function(
    data,
    K = 2,
    j0,
    k0,
    filterrows = FALSE,
    funcname = c("rangerlogit"),
    nfolds = 5,
    margin = 0.005,
    sl.lib = c(
        "SL.gam",
        "SL.glm",
        "SL.glm.interaction",
        "SL.ranger",
        "SL.glmnet"
    ),
    Nmin = 500,
    TMLE = TRUE,
    PLUGIN = TRUE,
    ...
) {
    requireNamespace("dplyr", quietly = TRUE, warn.conflicts = FALSE)
    requireNamespace("tidyr")
    l <- ncol(data) - K
    n <- nrow(data)
    stopifnot(!is.null(dim(data)))
    if (!informat(data = data, K = K)) {
        data <- reformat(data = data, capturelists = 1:K)
    }
    data <- na.omit(data)
    if (filterrows) {
        data <- data[which(rowSums(data[, 1:K]) > 0), ]
    }
    data <- as.data.frame(data)
    N <- nrow(data)
    stopifnot(N > 1)
    if (l >= 0 & N < Nmin) {
        l <- 0
        warning(cat(
            "Insufficient number of observations for doubly-robust estimation."
        ))
    }
    conforminglists <- apply(data[, 1:K], 2, function(col) {
        return(setequal(col, c(0, 1)))
    })
    if (sum(conforminglists) < 2) {
        stop("Data is not in the required format or lists are degenerate.")
        return(NULL)
    }
    if (sum(conforminglists) < K) {
        message(cat(
            "Lists ",
            which(conforminglists == FALSE),
            " are not in the required format."
        ))
    }
    if (!missing(j0)) {
        list1_vec <- j0
    } else {
        list1_vec <- c(1:(K - 1))
    }
    if (!missing(k0)) {
        list2_vec <- k0
    } else {
        list2_vec <- c(1:K)
    }
    if (l == 0) {
        colnames(data) <- c(paste("L", 1:K, sep = ""))
        listpair <- unlist(sapply(list1_vec, function(j1) {
            sapply(
                setdiff(list2_vec, list1_vec[list1_vec <= j1]),
                function(k1) {
                    return(paste(min(j1, k1), ",", max(j1, k1), sep = ""))
                }
            )
        }))
        psiinv <- data.frame(listpair = listpair)
        psiinv$psiin <- NA
        psiinv$sigma <- NA
        for (j in list1_vec) {
            j0 <- j
            if (!setequal(data[, j], c(0, 1))) {
                next
            }
            for (k in setdiff(list2_vec, list1_vec[list1_vec <= j0])) {
                if (!setequal(data[, k], c(0, 1))) {
                    next
                }
                if (j0 > k) {
                    j <- k
                    k <- j0
                } else {
                    j <- j0
                }
                q1 <- mean(data[, j])
                q2 <- mean(data[, k])
                q12 <- mean(data[, j] * data[, k])
                psiinv[psiinv$listpair == paste0(j, ",", k), ]$psiin <- pmax(
                    q1 * q2 / q12,
                    1
                )
                psiinv[psiinv$listpair == paste0(j, ",", k), ]$sigma <- sqrt(
                    q1 * q2 * pmax(q1 * q2 - q12, 0) * (1 - q12) / q12^3 / N
                )
            }
        }
        result <- psiinv %>%
            mutate(
                psi = 1 / psiin,
                sigma = sqrt(N) *
                    sigma,
                n = round(N * psiin),
                sigman = sqrt(
                    N^2 *
                        sigma^2 +
                        N * psiin * (psiin - 1)
                ),
                cin.l = round(pmax(
                    N *
                        psiin -
                        1.96 *
                            sqrt(
                                N^2 *
                                    sigma^2 +
                                    N *
                                        psiin *
                                        (psiin -
                                            1)
                            ),
                    N
                )),
                cin.u = round(
                    N *
                        psiin +
                        1.96 *
                            sqrt(
                                N^2 *
                                    sigma^2 +
                                    N * psiin * (psiin - 1)
                            )
                )
            ) %>%
            as.data.frame()
        result <- subset(result, select = -c(psiin))
        object <- list(result = result, N = N)
        class(object) <- "popsize"
        return(object)
    } else {
        colnames(data) <- c(
            paste("L", 1:K, sep = ""),
            paste("x", 1:(ncol(data) - K), sep = "")
        )
        if (nfolds > 1 & nfolds > N / 50) {
            nfolds <- pmax(floor(N / 50), 1)
            cat(
                "nfolds is reduced to ",
                nfolds,
                " to have sufficient test data.\n"
            )
        }
        listpair <- unlist(sapply(list1_vec, function(j1) {
            sapply(
                setdiff(list2_vec, list1_vec[list1_vec <= j1]),
                function(k1) {
                    return(paste(min(j1, k1), ",", max(j1, k1), sep = ""))
                }
            )
        }))
        psiinv_summary <- matrix(
            0,
            nrow = length(listpair),
            ncol = 3 *
                length(funcname)
        )
        rownames(psiinv_summary) <- listpair
        colnames(psiinv_summary) <- paste(
            rep(funcname, each = 3),
            c("PI", "DR", "TMLE"),
            sep = "."
        )
        var_summary <- psiinv_summary
        ifvals <- matrix(
            NA,
            nrow = N * length(listpair),
            ncol = length(funcname) +
                1
        )
        colnames(ifvals) <- c("listpair", funcname)
        ifvals[, "listpair"] <- rep(rownames(psiinv_summary), each = N)
        nuis <- matrix(
            NA,
            nrow = N * length(listpair),
            ncol = 3 *
                length(funcname) +
                1
        )
        colnames(nuis) <- c(
            "listpair",
            paste(rep(funcname, each = 3), c("q12", "q1", "q2"), sep = ".")
        )
        nuis <- as.data.frame(nuis)
        sapply(nuis, "class")
        nuis[, "listpair"] <- ifvals[, "listpair"]
        if (TMLE) {
            nuistmle <- nuis
        }
        permutset <- sample(1:N, N, replace = FALSE)
        for (j in list1_vec) {
            j0 <- j
            if (!setequal(data[, j], c(0, 1))) {
                next
            }
            for (k in setdiff(list2_vec, list1_vec[list1_vec <= j0])) {
                if (!setequal(data[, k], c(0, 1))) {
                    next
                }
                if (j0 > k) {
                    j <- k
                    k <- j0
                } else {
                    j <- j0
                }
                psiinvmat <- matrix(
                    numeric(0),
                    nrow = nfolds,
                    ncol = 3 * length(funcname)
                )
                colnames(psiinvmat) <- paste(
                    rep(funcname, each = 3),
                    c("PI", "DR", "TMLE"),
                    sep = "."
                )
                varmat <- psiinvmat
                ifvalsfold <- matrix(
                    numeric(0),
                    nrow = N,
                    ncol = length(funcname)
                )
                colnames(ifvalsfold) <- funcname
                nuisfold <- matrix(
                    numeric(0),
                    nrow = N,
                    ncol = 3 *
                        length(funcname)
                )
                colnames(nuisfold) <- paste(
                    rep(funcname, each = 3),
                    c("q12", "q1", "q2"),
                    sep = "."
                )
                nuistmlefold <- nuisfold
                idfold <- rep(1, N)
                for (folds in 1:nfolds) {
                    if (nfolds == 1) {
                        List1 <- data
                        List2 <- List1
                        sbset <- 1:N
                    } else {
                        sbset <- ((folds - 1) *
                            ceiling(N / nfolds) +
                            1):(folds * ceiling(N / nfolds))
                        sbset <- sbset[sbset <= N]
                        List1 <- data[permutset[-sbset], ]
                        List2 <- data[permutset[sbset], ]
                        idfold[permutset[sbset]] <- folds
                    }
                    yj <- List2[, paste("L", j, sep = "")]
                    yk <- List2[, paste("L", k, sep = "")]
                    overlapjk <- mean(List1[, j] * List1[, k])
                    if (overlapjk < margin) {
                        warning(cat(
                            "Overlap between the lists ",
                            j,
                            " and ",
                            k,
                            " is less than ",
                            margin,
                            ".\n",
                            sep = ""
                        ))
                    }
                    for (func in funcname) {
                        qhat <- try(
                            get(paste0("qhat_", func))(
                                List.train = List1,
                                List.test = List2,
                                K,
                                j,
                                k,
                                margin = margin,
                                ...
                            ),
                            silent = TRUE
                        )
                        if ("try-error" %in% class(qhat)) {
                            next
                        }
                        q12 <- qhat$q12
                        q1 <- pmin(pmax(q12, qhat$q1), 1)
                        q2 <- pmax(q12 / q1, pmin(qhat$q2, 1 + q12 - q1, 1))
                        nuisfold[
                            permutset[sbset],
                            paste(func, c("q12", "q1", "q2"), sep = ".")
                        ] <- cbind(q12, q1, q2)
                        gammainvhat <- q1 * q2 / q12
                        psiinvhat <- mean(gammainvhat, na.rm = TRUE)
                        phihat <- gammainvhat *
                            (yk /
                                q2 +
                                yj / q1 -
                                yj *
                                    yk /
                                    q12) -
                            psiinvhat
                        ifvalsfold[permutset[sbset], func] <- phihat
                        Qnphihat <- mean(phihat, na.rm = TRUE)
                        psiinvhat.dr <- max(psiinvhat + Qnphihat, 1)
                        psiinvmat[
                            folds,
                            paste(func, c("PI", "DR"), sep = ".")
                        ] <- c(psiinvhat, psiinvhat.dr)
                        sigmasq <- var(phihat, na.rm = TRUE)
                        varmat[
                            folds,
                            paste(func, c("PI", "DR"), sep = ".")
                        ] <- sigmasq / N
                        datmat <- as.data.frame(cbind(
                            yj,
                            yk,
                            yj *
                                yk,
                            q1 - q12,
                            q2 - q12,
                            q12
                        ))
                        colnames(datmat) <- c(
                            "yj",
                            "yk",
                            "yjk",
                            "q10",
                            "q02",
                            "q12"
                        )
                        if (TMLE) {
                            tmle <- tmle(
                                datmat = datmat,
                                margin = margin,
                                K = K,
                                ...
                            )
                        } else {
                            next
                        }
                        if (tmle$error) {
                            warning("TMLE did not run or converge.")
                            psiinvmat[
                                folds,
                                paste(func, "TMLE", sep = ".")
                            ] <- NA
                            varmat[folds, paste(func, "TMLE", sep = ".")] <- NA
                        } else {
                            datmat <- tmle$datmat
                            q12 <- pmax(datmat$q12, margin)
                            q1 <- pmin(datmat$q12 + datmat$q10, 1)
                            q2 <- pmax(
                                pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1),
                                q12 / q1
                            )
                            nuistmlefold[
                                permutset[sbset],
                                paste(func, c("q12", "q1", "q2"), sep = ".")
                            ] <- cbind(q12, q1, q2)
                            gammainvhat <- q1 * q2 / q12
                            psiinvhat.tmle <- mean(gammainvhat, na.rm = TRUE)
                            phihat <- gammainvhat *
                                (yj / q1 + yk / q2 - yj * yk / q12) -
                                psiinvhat.tmle
                            Qnphihat <- mean(phihat, na.rm = TRUE)
                            psiinvmat[
                                folds,
                                paste(func, "TMLE", sep = ".")
                            ] <- psiinvhat.tmle
                            sigmasq <- var(phihat, na.rm = TRUE)
                            varmat[
                                folds,
                                paste(func, "TMLE", sep = ".")
                            ] <- sigmasq / N
                        }
                    }
                }
                psiinv_summary[paste0(j, ",", k), ] <- colMeans(
                    psiinvmat,
                    na.rm = TRUE
                )
                var_summary[paste0(j, ",", k), ] <- colMeans(
                    varmat,
                    na.rm = TRUE
                )
                ifvals[
                    ifvals[, "listpair"] == paste0(j, ",", k),
                    colnames(ifvals) != "listpair"
                ] <- ifvalsfold
                nuis[
                    nuis[, "listpair"] == paste0(j, ",", k),
                    colnames(nuis) != "listpair"
                ] <- nuisfold
                if (TMLE) {
                    nuistmle[
                        nuistmle[, "listpair"] == paste0(j, ",", k),
                        colnames(nuistmle) != "listpair"
                    ] <- nuistmlefold
                }
            }
        }
        result <- list(
            psi = 1 / psiinv_summary,
            sigma = sqrt(
                N *
                    var_summary
            ),
            n = round(N * psiinv_summary),
            sigman = sqrt(
                N^2 *
                    var_summary +
                    N *
                        psiinv_summary *
                        (psiinv_summary -
                            1)
            ),
            cin.l = round(pmax(
                N *
                    psiinv_summary -
                    1.96 *
                        sqrt(
                            N^2 *
                                var_summary +
                                N *
                                    psiinv_summary *
                                    (psiinv_summary -
                                        1)
                        ),
                N
            )),
            cin.u = round(
                N *
                    psiinv_summary +
                    1.96 *
                        sqrt(
                            N^2 *
                                var_summary +
                                N * psiinv_summary * (psiinv_summary - 1)
                        )
            )
        )
        result <- Reduce(
            function(...) merge(..., by = c("listpair", "Var2")),
            lapply(1:length(result), function(i) {
                reshape2::melt(
                    result[[i]],
                    value.name = names(result)[i],
                    varnames = c("listpair", "Var2")
                )
            })
        )
        result <- separate(
            data = result,
            col = "Var2",
            into = c("model", "method"),
            sep = "\\."
        )
        if (!TMLE) {
            result <- result[result$method != "TMLE", ]
        }
        if (!PLUGIN) {
            result <- result[result$method != "PI", ]
        } else {
            warning(
                "Plug-in variance is not well-defined. Returning variance evaluated using DR estimator formula"
            )
        }
        object <- list(
            result = result,
            N = N,
            ifvals = as.data.frame(ifvals),
            nuis = as.data.frame(nuis),
            idfold = idfold
        )
        if (TMLE) {
            object$nuistmle <- as.data.frame(nuistmle)
        }
        class(object) <- "popsize"
        return(object)
    }
}

qhat_logit <- function(
    List.train,
    List.test,
    K = 2,
    j = 1,
    k = 2,
    margin = 0.005,
    ...
) {
    stopifnot(ncol(List.train) > K)
    if (missing(List.test)) {
        List.test <- List.train
    }
    colnames(List.train) <- c(
        paste("L", 1:K, sep = ""),
        paste("x", 1:(ncol(List.train) - K), sep = "")
    )
    colnames(List.test) <- c(
        paste("L", 1:K, sep = ""),
        paste("x", 1:(ncol(List.train) - K), sep = "")
    )
    fitj0 <- try(glm(
        formula(paste("L", j, "*(1 - L", k, ") ~.", sep = "")),
        family = binomial(link = "logit"),
        data = List.train[,
            c(j, k, (K + 1):ncol(List.train))
        ]
    ))
    fit0k <- try(glm(
        formula(paste("L", k, "*(1 - L", j, ") ~.", sep = "")),
        family = binomial(link = "logit"),
        data = List.train[,
            c(j, k, (K + 1):ncol(List.train))
        ]
    ))
    fitjk <- try(glm(
        formula(paste("L", j, "*L", k, " ~.", sep = "")),
        family = binomial(link = "logit"),
        data = List.train[,
            c(j, k, (K + 1):ncol(List.train))
        ]
    ))
    if ("try-error" %in% c(class(fitj0), class(fit0k), class(fitjk))) {
        warning("One or more fits with logistic regression failed.")
        return(NULL)
    } else {
        q12 <- pmax(
            predict(fitjk, newdata = List.test, type = "response"),
            margin
        )
        q1 <- pmax(
            q12 + predict(fitj0, newdata = List.test, type = "response"),
            margin
        )
        q2 <- pmax(
            q12 + predict(fit0k, newdata = List.test, type = "response"),
            margin
        )
        return(list(q1 = q1, q2 = q2, q12 = q12))
    }
}

qhat_rangerlogit <- function(
    List.train,
    List.test,
    K = 2,
    j = 1,
    k = 2,
    margin = 0.005,
    ...
) {
    requireNamespace("ranger", quietly = TRUE)
    requireNamespace("nnls")
    l <- ncol(List.train) - K
    if (missing(List.test)) {
        List.test <- List.train
    }
    stopifnot(l > 0)
    colnames(List.train) <- c(
        paste("L", 1:K, sep = ""),
        paste("x", 1:(ncol(List.train) - K), sep = "")
    )
    colnames(List.test) <- c(
        paste("L", 1:K, sep = ""),
        paste("x", 1:(ncol(List.train) - K), sep = "")
    )
    fitj <- ranger(
        formula(paste("factor(L", j, ") ~.", sep = "")),
        data = List.train[, -c(1:K)[-j]],
        probability = TRUE,
        classification = TRUE
    )
    fitk <- ranger(
        formula(paste("factor(L", k, ") ~.", sep = "")),
        data = List.train[, -c(1:K)[-k]],
        probability = TRUE,
        classification = TRUE
    )
    fitjk <- ranger(
        formula(paste("factor(L", j, "*L", k, ") ~.", sep = "")),
        data = List.train[, c(j, k, K + 1:l)],
        probability = TRUE,
        classification = TRUE
    )
    if ("try-error" %in% c(class(fitj), class(fitk), class(fitjk))) {
        warning("One or more fits with GAM regression failed.")
        return(NULL)
    } else {
        q12.r <- pmax(
            predict(fitjk, data = List.test, type = "response")$predictions[,
                "1"
            ],
            margin
        )
        q1.r <- pmax(
            predict(fitj, data = List.test, type = "response")$predictions[,
                "1"
            ],
            margin
        )
        q2.r <- pmax(
            predict(fitk, data = List.test, type = "response")$predictions[,
                "1"
            ],
            margin
        )
        ql <- qhat_logit(
            List.train = List.train,
            List.test = List.test,
            K = K,
            j = j,
            k = k,
            margin = margin
        )
        q12.l <- ql$q12
        q1.l <- ql$q1
        q2.l <- ql$q2
        requireNamespace("nnls")
        A <- cbind(1, q12.r, q12.l)
        coef <- coef(nnls(
            b = List.test[, paste0("L", j)] *
                List.test[,
                    paste0("L", k)
                ],
            A = A
        ))
        coef <- coef / sum(coef)
        q12 <- A %*% coef
        A <- cbind(1, q1.r, q1.l)
        coef <- coef(nnls(b = List.test[, paste0("L", j)], A = A))
        coef <- coef / sum(coef)
        q1 <- pmax(A %*% coef, q12)
        A <- cbind(1, q2.r, q2.l)
        coef <- coef(nnls(b = List.test[, paste0("L", k)], A = A))
        coef <- coef / sum(coef)
        q2 <- pmin(pmax(A %*% coef, q12), 1 - q1 + q12)
        return(list(q1 = q1, q2 = q2, q12 = q12))
    }
}
