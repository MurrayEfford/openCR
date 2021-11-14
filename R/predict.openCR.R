###############################################################################
## package 'openCR'
## predict.openCR.R
## moved from methods 2017-12-25
## 2019-04-11 ignores contrasts arg of lpredictor?
## 2019-04-13 failed with single parameter when others fixed
## 2021-04-21 stratified
## 2021-07-30 makeNewData() renamed from openCR.make.newdata
###############################################################################

predict.openCRlist <- function (object, newdata = NULL, se.fit = TRUE,
  alpha = 0.05, savenew = FALSE, ...) {
  lapply(object, predict, newdata, se.fit, alpha, savenew, ...)
}
predict.openCR <- function (object, newdata = NULL, se.fit = TRUE, alpha = 0.05,
  savenew = FALSE, ...) {
  if (is.null(object$fit)) {
    warning ("empty (NULL) object")
    return(NULL)
  }
  if (is.null(newdata)) {
    # newdata <- openCR.make.newdata (object, ...)
    newdata <- makeNewData (object, ...)
  }
  nstrata <- length(strata(object$capthist))
  nsess <- sapply(primaryintervals(object), function(x) sum(x>0) + 1)
  
  if (! 'stratum' %in% names(newdata)) {
    newdata <- data.frame(stratum = factor(rep(1,nrow(newdata)), 
      levels = 1:nstrata), newdata)
  }
  newdata <- stringsAsFactors(newdata)   ## 2021-05-13
  
  if (! 'session' %in% names(newdata)) {
    # not guaranteed to work for varying session number among strata
    newdata <- data.frame(session = factor(rep(1,nrow(newdata)),
      levels = 1: max(nsess)), newdata)
  }
  if ('h2' %in% names(newdata)) newdata1 <- newdata[newdata$h2==1,] 
  else if ('h3' %in% names(newdata)) newdata1 <- newdata[newdata$h3==1,] 
  else newdata1 <- newdata  

  # used only for mlogit
  invm <- function(x,j) {
    tempmat <- matrix(x, nrow = nsess[j])
    as.numeric(apply(tempmat,2,invmlogit))
  }
  beta <- complete.beta(object)
  beta.vcv <- complete.beta.vcv(object)
  if (is.null(beta.vcv)) se.fit <- FALSE
  getfield <- function (x) {
    lpredictor (newdata = newdata, model = object$model[[x]],
      indx = object$parindx[[x]], beta = beta, field = x,
      beta.vcv = beta.vcv, validlevels = object$design$validlevels)
  }
  predict <- sapply (object$realnames, getfield, simplify = FALSE)
  
  if(se.fit) {
    realvcv <- vcov(object, realnames = names(predict), newdata = newdata, byrow = TRUE)
    nreal <- length(predict)
    realvcv <- array(unlist(realvcv), dim=c(nreal, nreal, nrow(newdata)))
    realSE <- apply(realvcv,3, function(x) suppressWarnings(sqrt(diag(x))))
    # bug fix 2019-04-13
    # if single parameter, apply() returns vector instead of matrix
    if (is.null(dim(realSE))) realSE <- matrix(realSE, nrow = 1)
    ####################
    rownames(realSE) <- names(predict)
  }
  if (!is.null(predict$pmix)) {
    nmix <- object$details$nmix
    temp <- matrix(predict$pmix$estimate, ncol = nmix)
    temp2 <- apply(temp, 1, function(est) logit(invmlogit(est)))
    predict$pmix$estimate <- as.numeric(t(temp2))
    predict$pmix$se <- NA    ## uncertain
  }
  z <- abs(qnorm(1-alpha/2))
  for (i in names(predict)) {
    nc <- ncol(predict[[i]])
    if (i %in% c('superN','superD')) {
      ## only one row per stratum for the 'super' parameters
      strata <- match(unique(predict[[i]]$stratum), predict[[i]]$stratum)
      predict[[i]] <- predict[[i]][strata,]
      predict[[i]]$session <- NULL
      predict[[i]]$t <- NULL
      predict[[i]]$h2 <- NULL
      predict[[i]]$h3 <- NULL
    }
    else if (i %in% c('tau')) {
      M <- object$details$M
      predict[[i]] <- predict[[i]][c(1:M, M),]
      predict[[i]][M+1,] <- rep(NA,3)
      predict[[i]]$session <- 1:(object$details$M+1)
      rownames(predict[[i]]) <- 1:(object$details$M+1)
    }
    else if (i %in% c('pmix')) {
      predict[[i]] <- predict[[i]][predict[[i]]$session==1,,drop=FALSE]
      predict[[i]]$session <- NULL
    }
    else {
      # replace factor session with integer codes
      if (!is.null(predict[[i]]$session)) predict[[i]]$session <- as.numeric(predict[[i]]$session) 
      predict[[i]]$t <- NULL
      # drop spurious mixtures
      if (('h2' %in% names(predict[[i]])) & (!('h2' %in% all.vars(object$model[[i]])))) {
        predict[[i]] <- predict[[i]][predict[[i]]$h2 == 1,]
        predict[[i]]$h2 <- NULL
      }
      if (('h3' %in% names(predict[[i]])) & (!('h3' %in% all.vars(object$model[[i]])))) {
        predict[[i]] <- predict[[i]][predict[[i]]$h3 == 1,]
        predict[[i]]$h3 <- NULL
      }
    }
    lpred  <- predict[[i]][,'estimate']
    if (i == 'b') {
      templist <- split(lpred, newdata1$stratum)
      tmp <- mapply(invm, templist, 1:length(templist), SIMPLIFY=FALSE)
      predict[[i]]$estimate <- as.numeric(unlist(tmp))
    }
    else
      if (i == 'tau') {
        tempmat <- matrix(lpred, nrow=object$details$M+1)
        predict[[i]]$estimate <- as.numeric(apply(tempmat,2,invmlogit))
      }
    else
      if (i == 'pmix') {
        predict[[i]]$estimate <- invmlogit(lpred)
      }
    else
      predict[[i]]$estimate <- untransform(lpred, object$link[[i]])
    if (se.fit) {
      selpred <- predict[[i]][,'se']
      if (i %in% c('b', 'tau')) {
        predict[[i]]$SE.estimate <- rep(NA,length(selpred))
        ## doubtful
        ## warning("doubtful confidence limits for 'b' or 'tau' model in predict() - to be tested")
        
        templist <- split(lpred-z*selpred, newdata1$stratum)
        tmp <- mapply(invm, templist, 1:length(templist), SIMPLIFY=FALSE)
        predict[[i]]$lcl <- as.numeric(unlist(tmp))
        templist <- split(lpred+z*selpred, newdata1$stratum)
        tmp <- mapply(invm, templist, 1:length(templist), SIMPLIFY=FALSE)
        predict[[i]]$ucl <- as.numeric(unlist(tmp))
      }
      else {
        se <- realSE[i,]
        se[is.na(predict[[i]]$estimate)] <- NA
        predict[[i]]$SE.estimate <- se[1:nrow(predict[[i]])]   ## only 1 for superN, superD
        predict[[i]]$lcl <- untransform(lpred-z*selpred, object$link[[i]])
        predict[[i]]$ucl <- untransform(lpred+z*selpred, object$link[[i]])
      }
      predict[[i]][is.na(predict[[i]])] <- NA
    }
    else {
      predict[[i]]$SE.estimate <- rep(NA, nrow(predict[[i]])) 
      predict[[i]]$lcl <- rep(NA, nrow(predict[[i]])) 
      predict[[i]]$ucl <- rep(NA, nrow(predict[[i]])) 
    }
    predict[[i]]$se <- NULL
    if (i == 'superN') {
      if (object$distribution == 'binomial') {
        if (nstrata>1) {
          sumfreq <- sapply(object$capthist, function(x) sum(covariates(x)$freq))
        }
        else {
          sumfreq <- sum(covariates(object$capthist)$freq)
        }
        ncf <- sumfreq[unique(newdata$stratum)]
        if (nrow(predict[[i]]) != nstrata) {
          warning ("predict.openCR struggles with superN ", 
            " when more than one row per stratum")
        }
        predict[[i]][,c('estimate','lcl','ucl')] <- 
          sweep(predict[[i]][,c('estimate','lcl','ucl')], MARGIN = 2, 
            FUN = '+', STATS = ncf)
      }
    }
    # primary session labels by stratum
    # sessionlabels is always stored by openCR.fit as a list
    sessnames <- lapply(object$sessionlabels, unlist)
    if (!is.null(sessnames) && !is.null(predict[[i]]$session) ) {
      nsess <- max(sapply(sessnames, length))
      sessnames <- lapply(sessnames, function(x) c(x, rep('', nsess-length(x))))
      namematrix <- matrix(unlist(sessnames), nrow = nstrata, byrow = TRUE)
      rn <- namematrix[cbind(as.numeric(predict[[i]]$stratum),
        as.numeric(predict[[i]]$session))]
      predict[[i]]$session <- rn
    }
  }
  if (savenew) attr(predict, 'newdata') <- newdata
  # drop 'stratum' if not a stratified analysis
  if (is.null(object$stratified) || !object$stratified) {
    predict <- lapply(predict, function(x) x[,names(x)!='stratum'])
  }
  predict
}
############################################################################################

