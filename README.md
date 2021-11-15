# openCR
Open Population Capture-Recapture

Functions for non-spatial open population analysis by Cormack-Jolly-Seber (CJS) and Jolly-Seber-
Schwarz-Arnason (JSSA) methods, and by spatially explicit extensions of these methods. The
methods build on Schwarz and Arnason (1996), Borchers and Efford (2008) and Pledger et al.
(2010) (see vignette for more comprehensive references and likelihood). The parameterisation of
JSSA recruitment is flexible (options include population growth rate λ, per capita recruitment f
and seniority γ). Spatially explicit analyses may assume home-range centres are static or allow
dispersal between primary sessions according to various probability kernels, including bivariate
normal (BVN) and bivariate t (BVT).

**openCR** 2.2.1 is available on [CRAN](https://CRAN.R-project.org/package=openCR).

See also --

www.otago.ac.nz/density/openCR.html

[openCR-vignette.pdf](https://www.otago.ac.nz/density/pdfs/openCR-vignette.pdf)

[openCR-manual.pdf](https://www.otago.ac.nz/density/pdfs/openCR-manual.pdf)

The code here is under development. It may be installed using
```
devtools::install_github("MurrayEfford/openCR")
```

Installation is rather slow, and produces voluminous messages because of the compilation of C++ code.


## References

  Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
  likelihood methods for capture--recapture studies. *Biometrics*
  **64**, 377--385.

  Efford, M. G. and Schofield, M. R. (2020) A spatial open-population capture--recapture model.
  *Biometrics* **76**, 392--402.

  Pledger, S., Pollock, K. H. and Norris, J. L. (2010) Open
  capture--recapture models with heterogeneity: II. Jolly-Seber
  model. *Biometrics* **66**, 883--890.

  Schwarz, C. J. and Arnason, A. N. (1996) A general methodology for the
  analysis of capture-recapture experiments in open
  populations. *Biometrics* **52**, 860--873.
