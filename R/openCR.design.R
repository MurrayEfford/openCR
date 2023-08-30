###############################################################################
## package 'openCR'
## openCR.design.R
## 2011 12 29
## 2012-06-08 robust design (J from intervals, not ncol(capthist))
## 2012-12-20 JSSARET
## 2017-11-19 robust design intervals
## 2017-12-10 expanded for secondary-session effects
## 2018-01-23 timecov added; does sessioncov work?
## 2018-02-23 Bsession added, bsession corrected
## 2018-10-29 CJSp1 argument
## 2018-11-22 bk, Bk added; bsession etc. redefined
## 2018-11-23 tidy up
## 2018-12-21 trap covariates
## 2020-10-19 agecov for recoding age effects
## 2021-04-18 stratified
## 2021-05-12 test for time-varying trap covariates allows for stratification
## 2021-09-27 usage bug fixed (used should have been binary)

## openCR.design()

################################################################################

openCR.design <- function (capthist, models, type, naive = FALSE, 
  stratumcov = NULL, sessioncov = NULL, timecov = NULL, agecov = NULL, 
  dframe = NULL, contrasts = NULL, initialage = 0, minimumage = 0, 
  maximumage = 1, agebreaks = NULL, CJSp1 = FALSE, ...) {
  
  ## Generate design matrix, reduced parameter array, and parameter index array (PIA)
  ## for each parameter
  ## 'capthist' must be of class 'capthist' or 'list'
  
  findvars.MS <- function (cov, vars, dimcov, scov = FALSE) {
    ## function to add covariates to a design data frame 'dframe'
    ## based closely on function of same name in 'secr'
    ## cov may be a dataframe or list of dataframes, one per stratum (R > 1),
    ## if list, then require predictors to appear in all strata
    ## uses pad1 from utility.R and insertdim from secr 
    
    if (is.null(cov) | (length(cov)==0) | (length(vars)==0)) return()
    else {
      found <- ''
      if (!is.data.frame(cov)) {
        ## therefore assume cov is a list, one per stratum
        ## care required as strata may differ in n, S, K
        ## solution is to pad to max [n,S,K] over strata
        if (!is.list(cov) | (R==1))
          stop ("irregular covariates; check stratum structure")
        cov <- lapply(cov, stringsAsFactors)
        covnames <- lapply(cov, names)
        varincov <- sapply(covnames, function(nam) vars %in% nam)
        if (length(vars)>1) found <- vars[apply(varincov,1,all)]
        else found <- vars[all(varincov)]
        for (variable in found) {
          ## if factor, check all levels the same 
          if (any(sapply(cov, function(x) is.factor(x[,variable])))) {
            baselevels <- levels(cov[[1]][,variable])
            lev <- lapply(cov, function(x) levels(x[,variable]))
            if (!all(sapply(lev[-1], function(x) identical(x,baselevels)))) {
              print(lev)
              stop("covariate factor levels differ between strata")
            }
          }
          ## Pad on first dimension. Is this sufficient for e.g. b?
          onestratum <- function(x,i) {
            xv <- x[,variable]
            if (scov) xv <- rep(xv, secondarysessions[[i]])
            pad1(xv, dims[dimcov[1]])
          }
          vals <- mapply(onestratum, cov, 1:length(cov))
          vals <- unlist(vals)
          dframe[,variable] <<- insertdim (vals, dimcov, dims)
        }
      }
      else
      {
        cov <- stringsAsFactors(cov)
        found <- names(cov) %in% vars
        if (is.data.frame(cov) & any(found)) {
          found <- names(cov)[found]
          values <- as.data.frame(cov[,found])
          names(values) <- found
          if (length(values)>0) {
            for (variable in found) {
              vals <- values[,variable]
              if (scov) vals <- rep(vals, secondarysessions)
              dframe[,variable] <<- insertdim (vals, dimcov, dims)
            }
          }
        }
      }
      vars <<- vars[!(vars %in% found)]
    }
  }
  
  # not tested for openCR 2.0
  findvars.covtime <- function (covindices, vars) {
    ## function to add time-specific trap covariates to a 
    ## design data frame 'dframe'
    ## covindices should be a list of numeric or character index vectors, one component per session
    if (R>1) stop ("time-varying detector covariates not implemented for stratified designs")
    dimcov <- c(2,3)   ## animal, secondarysession
    if (length(covindices[[1]]) != J)
      stop ("require one index per primary session")
    covnames <- names(covindices)
    found <- covnames[covnames %in% vars]
    vars <<- vars[!(vars %in% found)]

    for (variable in found) {
      firstcol <- zcov[,covindices[[1]][1]]
      factorlevels <- NULL
      if (is.factor(firstcol)) {
        ## all must have same levels!!
        factorlevels <- levels(firstcol)
      }
      
      getvals <- function (indices, zcov) {
        notOK <- is.na(zcov[,indices])
        if (any(notOK)) {
          warning ("covariate missing values set to -1")
          zcov[,indices][notOK] <- -1
        }
        mat <- as.matrix(zcov[,indices]) ## detectors x occasions
      }
      vals <- getvals(covindices[[variable]], zcov)
      vals <- vals[,primarysessions(intervals)]
      vals <- unlist(vals)  
      if (!is.null(factorlevels)) {
        vals <- factor(vals, factorlevels)
      }
      dframe[,variable] <<- secr::insertdim (vals, dimcov, dims)
    }
  }
  #--------------------------------------------------------------------------------
  
  models$settle <- NULL                     # handled separately
  
  npar     <- length(models)                # real parameters
  parnames <- names(models)                 # c('p','phi') for CJS
  vars     <- unique (unlist(sapply (models, all.vars)))
  MS       <- ms(capthist)    # multistratum
  trps     <- traps(capthist)
  trapcov  <- covariates(trps)
  if (MS) {
    R <- length(capthist)                                    # nstrata
    n <- max(sapply(capthist, nrow))                         # max over strata
    S <- max(sapply(capthist, ncol))                         # max over strata
    K <- if (!grepl('secr', type)) 1
    else max(sapply(trps, nrow))                             # max over strata
    stratumlevels <- session(capthist)
  }
  else {
    R <- 1
    n <- nrow(capthist)
    S <- ncol(capthist)
    K <- if (grepl('secr', type)) nrow(trps) else 1
    stratumlevels <- '1'
  }
  
  nmix  <- get.nmix(models)
  zcov  <- covariates(capthist)            # stratum-specific individual covariates
  agefactor <- NULL
  
  if (MS) {
    intervals <- lapply(capthist, intervals)
    if (is.null(intervals)) intervals <- lapply(sapply(capthist, ncol)-1, rep, x=1)
    J                 <- lapply(intervals, function(x) sum(x>0)+1)
    primarysession    <- lapply(intervals, primarysessions)
    secondarysessions <- lapply(primarysession, tabulate)
    getj              <- lapply(J, seq, from = 1)
    firstofsession    <- mapply(match, getj, primarysession, SIMPLIFY = FALSE )
    validlevels       <- lapply(J, getvalidlevels, type=type, parnames=parnames, CJSp1=CJSp1)
  }
  else {
    intervals         <- attr(capthist, 'intervals')
    if (is.null(intervals)) intervals <- rep(1, ncol(capthist)-1)
    J                 <- sum(intervals>0) + 1 
    primarysession    <- primarysessions(intervals)
    secondarysessions <- tabulate(primarysession)
    firstofsession    <- match(1:J, primarysession)
    validlevels       <- getvalidlevels (type, parnames, J, CJSp1)
  }
  #--------------------------------------------------------------------------
  
  if (!grepl('secr', type) & any(.openCRstuff$traplearnedresponses %in% vars))
    stop ("cannot use detector-specific predictor with non-spatial model")
  
  if (sum(.openCRstuff$learnedresponses %in% vars) > 1)
    stop ("model should not use more than one type of behavioural response")
  
  if (any(.openCRstuff$learnedresponses %in% vars) &
      packageDescription("openCR")$Version<"1.4.0")  # sunset < 1.4.0 for this msg
    warning ("learned response models were re-defined in version 1.3 - check vignette")
  
  #--------------------------------------------------------------------------
  
  sessioncov <- stdcovlist(sessioncov, 'scov', R, unlist(J))
  timecov    <- stdcovlist(timecov, 'tcov', R, if(MS) sapply(capthist, ncol) else ncol(capthist))
  agecov     <- stdcovlist(agecov, 'acov', R, (maximumage - minimumage + 1))
  stratumcov <- stdcovlist(stratumcov, 'stratumcov', 1, NULL)
  #--------------------------------------------------------------------------
  dims <- c(R,n,S,K,nmix)       # virtual dimensions
  dframenrow <- prod(dims)    # number of rows
  autovars <- c(.openCRstuff$learnedresponses, 'stratum', 'session', 't', 
    'tt', 'Session', 'h2', 'h3', 'age', 'Age', 'Age2')
  #--------------------------------------------------------------------------
  # user-specified dframe
  if (is.null(dframe)) {
    dframe <- data.frame(matrix(nrow=dframenrow, ncol=0))
    dframevars <- ""
  }
  else {
    if (nrow(dframe) !=  dframenrow )
      stop ("dframe should have ", R*n*S*K*nmix, " rows ( R*n*S*K*nmix )")
    dframevars <- names(dframe)
  }
  #--------------------------------------------------------------------------
  # stratum
  dframe$stratum <- factor( insertdim (stratumlevels, 1, dims),
    levels = stratumlevels)
  #--------------------------------------------------------------------------
  if (MS) {
    sess <- mapply(rep, x = getj, times = secondarysessions, SIMPLIFY = FALSE)
    sess <- unlist(lapply(sess, pad1, S))
    dframe$session <- factor(secr::insertdim (sess,c(3,1),dims))
    dframe$Session <- secr::insertdim (sess-1,c(3,1),dims)
    dframe$tt <- factor(secr::insertdim (1:S, c(3,1), dims))
    if ('t' %in% vars) {
      dframe$t <- dframe$session
    }
  }
  else {
    dframe$session <- factor(secr::insertdim (rep(1:J, secondarysessions), 3, dims))
    dframe$Session <- secr::insertdim (rep(0:(J-1), secondarysessions), 3, dims)
    dframe$tt <- factor(secr::insertdim (1:S, 3, dims))
    ## t as synonym of session
    if ('t' %in% vars) {
      dframe$t <- factor(secr::insertdim (rep(1:J, secondarysessions), 3, dims))
    }
  }
  
  if (any(c('age','Age','Age2', names(agecov)) %in% vars)) {
    
    if (MS) {
      temp <- lapply(capthist, age.matrix, initialage, minimumage, maximumage)
      temp <- lapply(temp, padarray, c(n,S))
      age <- unlist(temp)
    }
    else {
      age <- age.matrix(capthist, initialage, minimumage, maximumage) # one stratum
    }
    
    if ('age' %in% vars) {
      # labels not used
      agefactor <- cut(age, breaks = agebreaks, right = FALSE)
      dframe$age <- factor(secr::insertdim (agefactor, c(2,3,1), dims))
    } 
    if ('Age' %in% vars) {
      dframe$Age <- secr::insertdim (age, c(2,3,1), dims)
    } 
    if ('Age2' %in% vars) {
      dframe$Age2 <- secr::insertdim (age^2, c(2,3,1), dims)
    } 
    for (i in names(agecov)) {
      if (i %in% vars) {
        agecovi <- agecov[,i][age-minimumage+1]
        dframe[,i] <- secr::insertdim (agecovi, c(2,3,1), dims)
      }
    }
  }
  #--------------------------------------------------------------------------
  
  ## behavioural response fields
  
  makeb <- function (caphist) {      ## global response
    temp0 <- apply(abs(caphist), 1:2, sum)
    ## convert to n x (S+1) CH
    t(apply(temp0, 1, prevcapt))
  }
  ## by primary session
  makebJ <- function (caphist) {      ## global response
    temp0 <- apply(abs(caphist), 1:2, sum)
    temp0 <- t(apply(temp0, 1, tapply, primarysession, sum))
    temp0 <- t(apply(temp0, 1, prevcapt))
    ## convert to n x (S+1) CH
    t(apply(temp0, 1, rep, secondarysessions))
  }
  #--------------------------------------------------------------------------
  if ('bsession' %in% vars) {
    if (naive) dframe$bsession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
        ## apply within each primary session
        unlist(lapply( split(x, primarysession), prevcapt1))
      }
      if (MS) {
        temp <- lapply(capthist, makeb)
        temp <- lapply(temp, padarray, c(n,S))
        temp <- unlist(temp)
      }
      else temp <- makeb(capthist)  # one stratum
      dframe$bsession <- secr::insertdim (as.vector(temp), c(2,3,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('Bsession' %in% vars) {
    if (naive) dframe$Bsession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, x[-length(x)]>0)
        ## apply within each primary session
        unlist(lapply( split(x, primarysession), prevcapt1))
      }
      if (MS) {
        temp <- lapply(capthist, makeb)
        temp <- lapply(temp, padarray, c(n,S))
        temp <- unlist(temp)
      }
      else temp <- makeb(capthist)  # one stratum
      dframe$Bsession <- secr::insertdim (as.vector(unlist(temp)), c(2,3,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('b' %in% vars) {
    if (naive) dframe$b <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
      if (MS) {
        temp <- lapply(capthist, makeb)
        temp <- lapply(temp, padarray, c(n,S))
        temp <- unlist(temp)
      }
      else temp <- makeb(capthist)  # one stratum
      dframe$b <- secr::insertdim (as.vector(temp), c(2,3,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('B' %in% vars) {
    if (naive) dframe$B <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
        ## apply within each primary session
        capt1 <- unlist(lapply( split(x, primarysession), prevcapt1))
        pcapt <- unique(primarysession[x>0])
        ## EITHER same primary as ct and later OR next primary
        (capt1 & (primarysession %in% pcapt)) | (primarysession %in% (pcapt+1))
      }
      if (MS) {
        temp <- lapply(capthist, makeb)
        temp <- lapply(temp, padarray, c(n,S))
      }
      else temp <- makeb(capthist)  # one stratum
      dframe$B <- secr::insertdim (as.vector(temp), c(2,3,1), dims)
    }
  }
  #------------------------------------------------
  ## individual trap-specific responses
  
  makebk <- function (caphist) {     
    if (nrow(caphist)==0) 
      array(dim = c(0,S,K))
    else {
      temp <- apply(abs(caphist), c(1,3), prevcapt)
      aperm(temp, c(2,1,3))
    }
  }
  
  #------------------------------------------------
  if ('bksession' %in% vars) {
    if (naive) dframe$bksession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
        ## apply within each primary session
        unlist(lapply( split(x, primarysession), prevcapt1))
      }
      if (MS) {
        temp <- lapply(capthist, makebk)
        temp <- lapply(temp, padarray, c(n,S,K))
      }
      else temp <- makebk(capthist)  # one stratum
      dframe$bksession <- secr::insertdim(temp, c(2,3,4,1), dims)  
    }
  }
  
  #------------------------------------------------
  if ('Bksession' %in% vars) {
    if (naive) dframe$Bksession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) c(FALSE, x[-S]>0)
      if (MS) {
        temp <- lapply(capthist, makebk)
        temp <- lapply(temp, padarray, c(n,S,K))
      }
      else temp <- makebk(capthist)  # one stratum
      dframe$Bksession <- secr::insertdim(temp, c(2,3,4,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('bk' %in% vars) {
    if (naive) dframe$bk <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
      if (MS) {
        temp <- lapply(capthist, makebk)
        temp <- lapply(temp, padarray, c(n,S,K))
      }
      else temp <- makebk(capthist)  # one stratum
      dframe$bk <- secr::insertdim(temp, c(2,3,4,1), dims)  
    }
  }
  
  #------------------------------------------------
  
  ## 2018-11-22
  ## Detector response (k) is an approximation because
  ## "naive" state refers to animals not detectors --
  ## undocumented for now (k matches secr)
  
  makek <- function (caphist) {      ## trap responds to capture of any animal
    temp <- apply(abs(caphist), c(2,3), sum) # occasion x trap
    apply(temp, 2, prevcapt)
  }
  
  #------------------------------------------------
  if ('k' %in% vars) {
    if (naive) dframe$k <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
      if (MS) {
        temp <- lapply(capthist, makek)
        temp <- lapply(temp, padarray, c(S,K))
      }
      else temp <- makek(capthist)  # one stratum
      dframe$k <- secr::insertdim(temp, c(3,4,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('ksession' %in% vars) {
    if (naive) dframe$ksession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, cumsum(x[-length(x)])>0)
        ## apply within each primary session
        unlist(lapply( split(x, primarysession), prevcapt1))
      }
      if (MS) {
        temp <- lapply(capthist, makek)
        temp <- lapply(temp, padarray, c(S,K))
      }
      else temp <- makek(capthist)  # one stratum
      dframe$ksession <- secr::insertdim(temp, c(3,4,1), dims)
    }
  }
  
  #------------------------------------------------
  if ('Ksession' %in% vars) {
    if (naive) dframe$Ksession <- rep(FALSE, dframenrow)
    else {
      prevcapt <- function(x) {
        prevcapt1 <- function(x) c(FALSE, x[-length(x)]>0)
        ## apply within each primary session
        unlist(lapply( split(x, primarysession), prevcapt1))
      }
      if (MS) {
        temp <- lapply(capthist, makek)
        temp <- lapply(temp, padarray, c(S,K))
      }
      else temp <- makek(capthist)  # one stratum
      dframe$Ksession <- secr::insertdim(temp, c(3,4,1), dims)
    }
  }
  #------------------------------------------------
  ## h2 or h3
  if (nmix > 1) {
    mixture <- paste('h',nmix,sep='')
    dframe[,mixture] <- secr::insertdim(factor(1:nmix), 5, dims)
  }
  
  #--------------------------------------------------------------------------
  ## all autovars should have now been dealt with
  vars <- vars[!(vars %in% c(autovars, names(agecov), dframevars))]
  
  #--------------------------------------------------------------------------
  
  # add miscellaneous covariates
  
  findvars.MS (stratumcov, vars, 1)
  findvars.MS (zcov, vars, 2)
  findvars.MS (timecov, vars, 3)
  findvars.MS (trapcov, vars, 4) 
  findvars.MS (sessioncov, vars, 3, scov = TRUE)  ## expands correctly 2018-05-07
  findvars.MS (agecov, vars, 2, scov = TRUE)      ## new 2020-10-19
  
  # time-varying trap covariates
  tvc <- timevaryingcov(capthist)
  if (!is.null(unlist(tvc)) && (length(vars)>0)) {
    findvars.covtime (tvc, vars)
  }
  
  if (length(vars)>0) {
    if (!is.null(zcov)) {
      if (is.data.frame(zcov))
        znames <- names(zcov)
      else
        znames <- unlist(lapply(zcov, names))
    }
    stop ("covariate(s) ", paste(vars,collapse=","), " not found")
  }
  make.designmatrix <- function (formula, prefix, ...) {
    # combine formula and dframe to generate design matrix
    if (is.null(formula)) {
      list (model = NULL, index = rep(1,dframenrow))
    }
    else {
      # # adjust for unidentifiable parameters
      # dframe <- adjustlevels(prefix, dframe, validlevels)
      # localcontrasts <- contrasts[names(contrasts) %in% all.vars(formula)]
      # if (length(localcontrasts)==0) localcontrasts <- NULL
      # tempmat <- model.matrix(formula, data = dframe, contrasts.arg = localcontrasts, ...)
      
      tempmat <- get.model.matrix(formula, prefix, dframe, validlevels, contrasts, ...)
      
      ## drop pmix beta0 column from design matrix
      if (prefix=='pmix') tempmat <- tempmat[,-1,drop=FALSE]
      ## temp <- secr::make.lookup (tempmat)   # retain unique rows
      temp <- makelookupcpp (tempmat)   # retain unique rows   ## 2018-11-06
      list (model=temp$lookup, index=temp$index)
    }
  }
  dframe[is.na(dframe)] <- 0
  # list with one component per real parameter
  # each of these is a list with components 'model' and 'index'
  designMatrices <- sapply (1:length(models), simplify=FALSE,
    function (x) make.designmatrix(models[[x]], names(models[x])))
  names(designMatrices) <- names(models)
  
  ## dim(indices) = c(n*S*K*nmix, npar)
  indices <- sapply (designMatrices, function(x) x$index)
  indices <- matrix(unlist(indices), ncol = npar)
  
  # retain just the 'model' components of 'designMatrices'
  designMatrices <- lapply (designMatrices, function(x)x$model )
  
  # prefix column names in 'designMatrices' with parameter name
  for (i in 1:npar) {
    colnames(designMatrices[[i]]) <- paste (parnames[i], '.',
      colnames(designMatrices[[i]]), sep='')
  }
  
  # repackage indices to define unique combinations of parameters
  ##indices2 <- secr::make.lookup(indices)
  indices2 <- makelookupcpp(indices)
  
  #--------------------------------------------------------------------
  # PIA = Parameter Index Array
  #       index to row of parameterTable for a given R,n,s,K,nmix
  # dim(parameterTable) = c(uniqueparcomb, npar)
  #       index to row of designMatrix for each real parameter
  #--------------------------------------------------------------------
  PIA <- array(indices2$index, dim = dims)
  PIAJ <- njx(PIA, primarysession)
  
  parameterTable <- indices2$lookup
  colnames(parameterTable) <- parnames
  
  #--------------------------------------------------------------------
  # Zero the index of trap+time pairs that were 'not set'
  # the external C code checks for this and sets p(detection) to zero
  #--------------------------------------------------------------------
  if (grepl('secr', type)) {
    for (stratum in 1:R) {
      if (MS) {
        # used <- usage(trps)[[stratum]]
        # binary from 2021-09-27
        used <- usage(trps)[[stratum]] > 0
      }
      else {
        # used <- usage(trps) 
        # binary from 2021-09-27
        used <- usage(trps) > 0
      }
      if ((!is.null(unlist(used))) & (length(used)>0)) {
        allused <- unlist(used)
        if (!is.null(allused)) {
          if (any(!allused)) {
            used <- padarray(t(used), c(S,K))
            PIA[stratum,,,,] <- PIA[stratum,,,,] * rep(rep(used,rep(n,S*K)),nmix)
          }
        }
      }
    }
  }
  individual <- individualcovariates(PIA)
  #--------------------------------------------------------------------
  
  list(
    designMatrices = designMatrices, 
    parameterTable = parameterTable, 
    PIA = PIA,
    PIAJ = PIAJ, 
    validlevels = validlevels, 
    individual = individual, 
    agelevels = levels(agefactor))
}
############################################################################################

## PIA has dim = c(R, nc, ss, kk, xx)
## primarysession gives primary session (j) for each secondary session (s)
# returns list with components 
#    lookup - matrix of unique rows, each row ss * kk long, containing indices as in PIA to 
#    the rows of realparval0
#    index  vector of length nc * jj * xx; values are indices to rows of lookup

# this could also be used in secr to speed up esa

# BUT NOT WORKING 2018-02-13
njxlookup <- function (PIA, primarysession) {
  dims <- lapply(dim(PIA), seq, from=1)
  names(dims) <- c('n','s','k','x')
  df <- data.frame(pia = as.numeric(PIA), do.call(expand.grid, dims))
  ss <- dim(PIA)[2]
  kk <- dim(PIA)[3]
  names(df) <- c('pia', 'n','s','k','x')
  df$j <- formatC(primarysession[df$s], width=3, flag="0")
  df$n <- formatC(df$n, width=4, flag="0")
  splitter <- apply(df[,c('x','j','n')], 1, paste, collapse='.')
  splitdf <- split(df[,c('pia','s')], splitter)
  fixpiask <- function (x) {
    pia <- matrix(0,ss,kk)
    pia[x$s,] <- x$pia
    as.numeric(pia)
  }
  piask <- lapply(splitdf, fixpiask)
  njxIA <- do.call(rbind, piask)
  ## lookup <- secr::make.lookup(njxIA)
  lookup <- makelookupcpp(njxIA)
  lookup
}

# direct index to session PIA in njx array 2018-02-10
# slow - must be a better way 2018-04-11
njx <- function (PIA, primarysession) {
  R <- dim(PIA)[1]
  n <- dim(PIA)[2]
  S <- dim(PIA)[3]
  xx <- dim(PIA)[5]
  maxJ <- max(unlist(primarysession))
  out <- array(1, dim = c(R,n,maxJ,xx))
  if (R>1) {
    for (i in 1:R) {
      J <- max(primarysession[[i]])
      s1 <- pad1(match(1:J, primarysession[[i]]), maxJ)
      out[i,,,] <- PIA[i,,s1,1,]
    }
  }
  else {
    J <- max(primarysession)
    s1 <- match(1:J, primarysession)
    out[1,,,] <- PIA[1,,s1,1,]
  }
  out
}
# 
# n <- dim(PIA)[1]
# J <- max(primarysession)
# xx <- dim(PIA)[4]
# # instead 2018-04-11   
# s1 <- match(1:J, primarysession)
# piaj <- PIA[,s1,1,]
# array(piaj, dim = c(n,J,xx))
# }