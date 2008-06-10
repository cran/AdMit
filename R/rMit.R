## Function which generates draws from a mixture of Student-t densities
## __input__
## N        : [integer>0] number of draws (default: 1)
## mit      : [list] of mixture information (default: univariate Cauchy)
##  $p      : [Hx1 vector] of probabilities
##  $mu     : [Hxk matrix] of mean vectors (in row)
##  $Sigma  : [Hxk^2 matrix] of scale matrices (in row)
##  $df     : [integer>0] degrees of freedom parameter
## __output__
## [Nxk matrix] of draws
## __20080429__
'rMit' <- function(N=1, mit=list())
  {
    H <- length(mit$p)
    if (H==0)
      { ## default
        warning ("'mit' not well defined in 'rMit'; set to default")
        mit <- list(p=1, mu=as.matrix(0), Sigma=as.matrix(1), df=1)
        H <- 1
      }
    comp <- sample(1:H, N, prob=mit$p, replace=TRUE)
    r <- NULL
    for (h in 1:H)
      { 
        nh <- length(comp[comp==h]) ## determine the nb of each component
        tmp <- NULL
        if (nh>0)
          tmp <- fn.rmvt(nh, mit$mu[h,], mit$Sigma[h,], mit$df)
        r <- rbind(r, tmp)
      }
    ## randomize the output (iid sample)
    r <- as.matrix(r[sample(1:N),])
    if (ncol(r)==1) 
      r <- as.vector(r)
    return(r)
  }



