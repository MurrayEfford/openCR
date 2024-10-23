\name{Internal}

\alias{prwi}
\alias{prwisecr}
\alias{PCH1}
\alias{PCH1secr}
\alias{pradelloglik}
\alias{cyclic.fit}

\title{Internal Functions}

\description{ 

Functions called by \code{openCR.fit} when \code{details$R == TRUE}, and some others 

}

\usage{

prwi (type, n, x, jj, cumss, nmix, w, fi, li, openval, PIA, PIAJ, intervals, CJSp1)

prwisecr (type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li, gk, openval, 
    PIA, PIAJ, binomN, Tsk, intervals, h, hindex, CJSp1, moveargsi, 
    movementcode, sparsekernel, edgecode, usermodel, kernel = NULL, 
    mqarray = NULL, cellsize = NULL, r0)

PCH1 (type, x, nc, cumss, nmix, openval0, PIA0, PIAJ, intervals)

PCH1secr (type, individual, x, nc, jj, cumss, kk, mm, openval0, PIA0, PIAJ, gk0,
    binomN, Tsk, intervals,  moveargsi, movementcode, sparsekernel, edgecode, 
    usermodel, kernel, mqarray, cellsize, r0) 

pradelloglik (type, w, openval, PIAJ, intervals)

cyclic.fit (..., maxcycle = 10, tol = 1e-5, trace = FALSE) 

}

\arguments{
  \item{type}{character}
  \item{n}{integer index of capture history}
  \item{x}{integer index of latent class}
  \item{jj}{integer number of primary sessions}
  \item{cumss}{integer vector cumulative number of secondary sessions at start of 
      each primary session}
  \item{nmix}{integer number of latent classes}
  \item{w}{array of capture histories}
  \item{fi}{integer first primary session}
  \item{li}{integer last primary session}
  \item{openval}{dataframe of real parameter values (one unique combination per row)}
  \item{PIA}{parameter index array (secondary sessions)}
  \item{PIAJ}{parameter index array (primary sessions)}
  \item{intervals}{integer vector}
  \item{h}{numeric 3-D array of hazard (mixture, mask position, hindex)}
  \item{hindex}{integer n x s matrix indexing h for each individual, secondary session}
  \item{CJSp1}{logical; should CJS likelihood include first primary session?}
  
  \item{moveargsi}{integer 2-vector for index of move.a, move.b (negative if unused)}
  \item{movementcode}{integer 0 static, 1 uncorrelated etc.}
  \item{sparsekernel}{logical; if TRUE then only cardinal and intercardinal axes are included}
  \item{edgecode}{integer 0 none, 1 wrap, 2 truncate}
  \item{usermodel}{function to fill kernel}
  \item{kernel}{dataframe with columns x,y relative coordinates of kernel cell centres}
  \item{mqarray}{integer matrix}
  \item{cellsize}{numeric length of side of kernel cell}
  \item{r0}{numeric; effective radius of zero cell for movement models (usually 0.5)}

  \item{gk}{real array}
  \item{Tsk}{array detector usage}

  \item{openval0}{openval for naive animals}
  \item{PIA0}{PIA for naive animals}
  \item{individual}{logical; TRUE if model uses individual covariates}
  \item{gk0}{gk for naive animals}
  
  \item{nc}{number of capture histories}
  \item{kk}{number of detectors}
  \item{mm}{number of points on habitat mask}
  \item{binomN}{code for distribution of counts (see \code{\link[secr]{secr.fit}})}
  
  \item{\dots}{named arguments passed to \code{\link{openCR.fit}} or \code{predict} (see extractFocal)}
  \item{maxcycle}{integer maximum number of cycles (maximizations of a given parameter)}
  \item{tol}{absolute tolerance for improvement in log likelihood}
  \item{trace}{logical; if TRUE a status message is given at each maximization}
  
}

\details{

  \code{cyclic.fit} implements cyclic fixing more or less as described by
  Schwarz and Arnason (1996) and used by Pledger et al. (2010). The
  intention is to speed up maximization when there are many (beta)
  parameters. However, fitting is slower than with a single call to
  \code{\link{openCR.fit}}, and the function is here only as a curiosity
  (it is not exported in 1.2.0).

}

\value{

\code{cyclic.fit} returns a fitted model object of class `openCR'.

Other functions return numeric components of the log likelihood.

}

\seealso{

\code{\link{openCR.fit}}

}

\references{

Pledger, S., Pollock, K. H. and Norris, J. L. (2010) Open
  capture--recapture models with heterogeneity: II. Jolly-Seber
  model. \emph{Biometrics} \bold{66}, 883--890.

Schwarz, C. J. and Arnason, A. N. (1996) A general methodology for the
  analysis of capture-recapture experiments in open
  populations. \emph{Biometrics} \bold{52}, 860--873.

}

\examples{

\dontrun{

openCR:::cyclic.fit(capthist = dipperCH, model = list(p~t, phi~t), tol = 1e-5, trace = TRUE)

}

}