## Function which computes the density of a mixture of Student-t densities
## __input__
## theta : [Nxk matrix] of values
## mit   : [list] of mixture information (default: univariate Cauchy)
## log    : [logical] log output (default: TRUE)
## __output__
## [Nx1 vector] of density evaluated at theta
## __20080329__
'dMit' <- function(theta, mit=list(), log=TRUE)
  {
    if (missing(theta))
      stop ("'theta' is missing in 'dMit'")
    H <- length(mit$p)
    if (H==0) 
      { ## default for the mixture
        warning ("'mit' not well defined in 'dMit'; set to default")
        mit <- list(p=1, mu=as.matrix(0), Sigma=as.matrix(1), df=1)
        H <- 1
      }
    if (is.vector(theta))
      { ## if vector is supplied instead of a matrix
        if (ncol(mit$mu)==1)
          { ## univariate density
            theta <- as.matrix(theta)
          }
        else
          { ## multivariate density but evaluated at a single point
            theta <- matrix(theta, nrow=1)
          }
      }
    r <- tmp <- 0
    for (h in 1:H)
      { 
        tmp <- exp( log(mit$p[h])+ fn.dmvt(theta, mit$mu[h,], mit$Sigma[h,], mit$df, log=TRUE) )
        r <- r+tmp
      }
    if (log)
      r <- log(r)
    as.vector(r)
  }
