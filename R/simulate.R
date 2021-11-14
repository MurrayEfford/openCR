################################################################################
## openCR 1.3.6
## simulate.R 
## 2018-05-28 intervals argument for runsim.spatial
## 2018-11-02 revised sumsims
## 2018-11-27 revised sumsims to allow method = "none"
## 2019-04-07 runsim.nonspatial, runsim.spatial return NULL for failed replicates
## 2021-05-13 sumsims failed with 'non-numeric argument to binary operator'
################################################################################
sim.nonspatial <- function(N, turnover = list(), p, nsessions, noccasions = 1, 
                           intervals = NULL, recapfactor = 1, seed = NULL, 
                           savepopn = FALSE, ...)  {

    if ((length(N)==2) & (length(p)==2)) {
        ## 2-class mixture
        if (!is.null(seed)) set.seed(seed)
        CH1 <- sim.nonspatial(N[1], turnover, p[1], nsessions, noccasions, intervals,
                              recapfactor, seed = NULL, savepopn = FALSE, ...)
        CH2 <- sim.nonspatial(N[2], turnover, p[2], nsessions, noccasions, intervals,
                              recapfactor, seed = NULL, savepopn = FALSE, ...)
        CH <- abind(CH1, CH2, along = 1)
        rownames(CH) <- 1:nrow(CH)
        class(CH) <- "capthist"
        sessionlabels(CH) <- sessionlabels(CH1)
        intervals(CH) <- intervals(CH1)
        CH
    }
    else {
        if (is.null(intervals)) {
            intervals <- rep(c(rep(0, noccasions-1),1), nsessions)[-nsessions*noccasions]
            if (nsessions<=1) stop ("sim.nonspatial requires nsessions > 1")
        }
        primary <- primarysessions(intervals)
        nsessions <- length(unique(primary))   ## number primary
        S <- length(primary)                   ## number secondary

        ## generate nsessions populations with required turnover
        primaryintervals <- intervals[match(2:nsessions, primary) - 1]
        turnover$phi <- turnover$phi^primaryintervals
        turnover$lambda <- turnover$lambda^primaryintervals
        core <- make.grid()
        if (!is.null(seed)) set.seed(seed)
        pop <- sim.popn (D = NA, Nbuffer = N[1], Ndist = "fixed", core = core,
                         nsessions = nsessions, details = turnover, ...)

        ## convert to logical array
        superN <- as.numeric(tail(row.names(pop[[nsessions]]),1))
        alive <- matrix(FALSE, nrow = superN, ncol = S)
        rownames(alive) <- 1:superN
        for (j in 1:nsessions) {
            alive[rownames(pop[[j]]), primary == j] <- TRUE
        }

        ## sample for each secondary session
        p <- rep(p, length.out = S)
        p <- matrix(p, nrow = superN, ncol = S, byrow = TRUE)
        p <- p*alive
        one <- function (px) {
            ch <- runif(S) < px
            if (any(ch) & (recapfactor != 1)) {
                first <- match(TRUE, ch>0)
                if (first<S) {
                    later <- (first+1):S
                    px1 <- px[later] * recapfactor
                    ch[later] <- runif(length(later)) < px1
                }
            }
            ch
        }
        CH <- t(apply(p, 1, one))
        CH <-  array(CH, dim = c(superN, S, 1))
        rownames(CH) <- 1:superN

        ## drop null capture histories
        caught <- apply(CH,1,sum) > 0
        CH <- CH[caught,,,drop = FALSE]

        ## return capthist object
        class(CH) <- "capthist"
        sessionlabels(CH) <- 1:nsessions
        intervals(CH) <- intervals
        if (savepopn) attr(CH, 'popn') <- pop
        CH
    }
}
################################################################################

runsim.nonspatial <- function (nrepl = 100, seed = NULL, ncores = NULL, 
                               fitargs = list(), extractfn = predict, ...) {
    onesim <- function(r) {
        CH <- sim.nonspatial(...)  # does not pass seed
        fitargs$capthist <- CH
        fit <- try(do.call(openCR.fit, fitargs))
        if (!inherits(fit, 'try-error')) {   ## 2019-04-07
            out <- extractfn(fit)
            attr(out, 'eigH') <- fit$eigH
            attr(out, 'fit') <- fit$fit
            message("Completed replicate ", r, "  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
        }
        else {
            out <- NULL
            message("Failed replicate ", r, "  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
        }
        out
    }
    list(...)
    if (!is.null(ncores)) fitargs$ncores <- 1       # use multiple cores only at level of replicates
    message("Start  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
    ## note use of || for first item evaluation 2018-04-26
    if (is.null(ncores) || (ncores == 1)) {
        if (!is.null(seed)) set.seed(seed)
        out <- lapply(1:nrepl, onesim)
    }
    else {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        clusterSetRNGStream(clust, seed)
        clusterExport(clust, c('fitargs','extractfn'), environment())
        clusterEvalQ(clust, library(openCR))
        out <- parSapply (clust, 1:nrepl, onesim, simplify = FALSE)
        stopCluster(clust)
    }
    message("Finish ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
    out
}
################################################################################

runsim.spatial <- function(nrepl = 100, seed = NULL, ncores = NULL,
                           popargs = list(), detargs = list(), fitargs = list(),
                           extractfn = predict, intervals = NULL)  {

    onesim <- function(r) {
        detargs$popn <- do.call(sim.popn, popargs)
        detargs$renumber <- FALSE
        if (is.null(detargs$traps)) detargs$traps <- popargs$core
        fitargs$capthist <- do.call(sim.capthist, detargs)
        intervals(fitargs$capthist) <- intervals
        fit <- try(do.call(openCR.fit, fitargs))
        if (!inherits(fit, 'try-error')) {   ## 2019-04-07
            out <- extractfn(fit)
            if (is.list(fit)) {              ## 2019-06-10 
                attr(out, 'eigH') <- fit$eigH
                attr(out, 'fit') <- fit$fit
            }
            message("Completed replicate ", r, "  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
        }
        else {
            out <- NULL
            message("Failed replicate ", r, "  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
        }
        out
    }
    popargs$seed <- NULL
    detargs$seed <- NULL
    if (!is.null(ncores)) fitargs$ncores <- 1     # use multiple cores only at level of replicates
    message("Start  ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
    ## note use of || for first item evaluation 2018-04-26
    if (is.null(ncores) || (ncores == 1)) {
        if (!is.null(seed)) set.seed(seed)
        out <- lapply(1:nrepl, onesim)
    }
    else {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        clusterSetRNGStream(clust, seed)
        
        clusterExport(clust, c('popargs','detargs','fitargs','extractfn', 'intervals'), environment())
        out <- parSapply (clust, 1:nrepl, onesim, simplify = FALSE)
        stopCluster(clust)
    }
    message("Finish ", format(Sys.time(), "%H:%M:%S %d %b %Y"))
    out
}
################################################################################

sumsims <- function (sims, parm = 'phi', session = 1, dropifnoSE = TRUE, 
                     svtol = NULL, maxcode = 3, true = NULL) {
    if (length(session)>1) {
        results <- lapply(session, sumsims, sims=sims, parm = parm, dropifnoSE = dropifnoSE,
                          svtol = svtol, maxcode = maxcode, true = true)
        names(results) <- paste('session',session)
        results
    }
    else {
        if (is.null(sims)) {
            data.frame(cbind(median = NA, mean = NA, sd = NA, n = NA))
        }
        else {
            ## adapted 2018-10-28 to use output from either predict or summary
            if (is.null(sims[[1]]$predicted))
                predicted <- sims ## assume extractfn = predict
            else
                predicted <- lapply(sims, '[[', 'predicted')
            if (!parm %in% names(predicted[[1]])) {
                stop ("Parameter ", parm, " not reported (try ", paste(names(predicted[[1]]), collapse=', '), ")")
            }
            mat <- t(sapply(predicted, function(x) unlist(x[[parm]][session,])))
            class(mat) <- 'numeric'   ## 2021-05-13
            if (dropifnoSE)
                OK <- !is.na(mat[,3]) & mat[,3]>0   ## 0 condition added 2019-06-15
            else
                OK <- rep(TRUE, nrow(mat))

            if (attr(sims[[1]], "fit")$method != "none") {
                if (is.na(maxcode)) maxcode <- Inf
                codes <- sapply(lapply(sims, attr, 'fit'), '[[', 'code')
                if (is.null(codes[[1]]))
                    codes <- sapply(lapply(sims, attr, 'fit'), '[[', 'convergence')
                # 2018-11-27
                codes <- unlist(codes)
                if (any(is.null(codes)))
                    stop ("problem finding convergence codes")
                OK <- OK & (codes<= maxcode)    
            }
            if (!is.null(svtol)) {
                nbeta <- function(x) length(attr(x, 'fit')$par)
                NP <- sapply(sims, nbeta)  # nrepl vector, all same
                eigH <- sapply(sims, attr, 'eigH')          # matrix np x nrepl
                rank <- apply(eigH>svtol,2,sum)
                OK <- OK & (rank == NP)
            }
            if (ncol(mat)>4) {    # drop 'session' column if present (not superD)
                mat <- mat[OK, c('estimate', 'SE.estimate', 'lcl','ucl'), drop = FALSE]  ## 2018-05-05, 2018-05-28, 2018-11-02
            }
            mat <- as.data.frame(mat)
            mat$RSE <- mat$SE.estimate / mat$estimate
            mat$CI.length <- mat$ucl - mat$lcl
            if (!is.null(true)) {
                mat$Bias <- mat$estimate - true
                mat$RB <- mat$Bias / true
            }
            out <- data.frame(median = apply(mat, 2, median, na.rm = TRUE),
                              mean = apply(mat, 2, mean, na.rm = TRUE),
                              sd = apply(mat, 2, sd, na.rm = TRUE),
                              n = apply(mat, 2, function(x) sum(!is.na(x))))
            if (!is.null(true)) {
                nr <- nrow(out)
                out$rRMSE <- round(c(mean((mat$estimate - true)^2, na.rm = TRUE), rep(NA,nr-1))^0.5 / true,4)
                out$COV <- c(mean((true>mat$lcl) & (true< mat$ucl), na.rm = TRUE), rep(NA,nr-1))
            }
            return(out)
        }
    }
}
################################################################################

runsim.RMark <- function (nrepl = 100, model = 'CJS', model.parameters = NULL, 
                          extractfn, seed = NULL, ...) {
    onesim <- function(r) {
        CH <- sim.nonspatial(...)
        ch <- RMarkInput(CH)
        fit <- RMark::mark(ch, model = model, model.parameters = model.parameters,
                           invisible = TRUE, output = FALSE)
        extractfn(fit)
    }
    if (!requireNamespace('RMark'))
        stop ("Please install package RMark")
    if (!all (nchar(Sys.which(c('mark.exe', 'mark64.exe', 'mark32.exe'))) < 2)) {
        if (!is.null(seed)) set.seed(seed)
        lapply(1:nrepl, onesim)
    } else {
        message ("MARK executable not found; set e.g. MarkPath = 'c:/Mark/'")
    }
}
################################################################################

