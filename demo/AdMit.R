'wait' <- function()
  {
    t <- readline("\nPlease 'q' to quit the demo or any other key to continue...\n")
    if (t == "q") stop ("end of the demo")
  }

## Initialization
options(digits=4, max.print=40)

## Gelman and Meng (2001) kernel function
'GelmanMeng' <- function(x, A=1, B=0, C1=3, C2=3, log=TRUE)
  {
    if (is.vector(x))
      x <- matrix(x, nrow=1)
    r <- -.5 * (A*x[,1]^2*x[,2]^2 + x[,1]^2 + x[,2]^2
                - 2*B*x[,1]*x[,2] - 2*C1*x[,1] - 2*C2*x[,2])
    if (!log)
      r <- exp(r)
    as.vector(r)
  }
wait()

## Contour plot of the Gelman and Meng (2001) kernel function.
'PlotGelmanMeng' <- function(x1, x2)
  {
    GelmanMeng(cbind(x1,x2), log=FALSE)
  }
x1 <- x2 <- seq(from=-1, to=6, by=0.1)
z <- outer(x1, x2, FUN=PlotGelmanMeng)
contour(x1, x2, z, nlevel=20, las=1, lwd=2, col=rainbow(20),
        xlab=expression(X[1]), ylab=expression(X[2]))
abline(a=0, b=1, lty='dotted')
wait()

## Use AdMit to the Gelman and Meng (2001) kernel function.
set.seed(1234)
outAdMit <- AdMit(GelmanMeng, mu0=c(0,0.1))
print(outAdMit)
wait()

## Contour plot of the mixture approximation obtained with AdMit.
'PlotMit' <- function(x1, x2, mit)
  {
    dMit(cbind(x1, x2), mit=mit, log=FALSE)
  }
z <- outer(x1, x2, FUN=PlotMit, mit=outAdMit$mit)
contour(x1, x2, z, nlevel=20, las=1, lwd=2, col=rainbow(20),
        xlab=expression(X[1]), ylab=expression(X[2]))
abline(a=0, b=1, lty='dotted')
wait()

## Contour plot of the components of the mixture approximation
## obtained with AdMit.
par(mfrow=c(2,2))
for (h in 1:4)
  {
    mith <- list(p=1,
                 mu=outAdMit$mit$mu[h,,drop=FALSE],
                 Sigma=outAdMit$mit$Sigma[h,,drop=FALSE],
                 df=outAdMit$mit$df)
    z <- outer(x1, x2, FUN=PlotMit, mit=mith)
    contour(x1, x2, z, las=1, nlevel=20, lwd=2, col=rainbow(20),
            xlab=expression(X[1]), ylab=expression(X[2]))
    abline(a=0, b=1, lty='dotted')
    title(main=paste("component nr.", h))
  }
wait()

## Use importance sampling with the mixture approximation
## as the importance density.
outAdMitIS <- AdMitIS(KERNEL=GelmanMeng, mit=outAdMit$mit)
print(outAdMitIS)
wait()

## Use an alternative 'G' function in importance sampling
## for computing the covariance matrix.
'G.cov' <- function(theta, mu)
  {
    'G.cov_sub' <- function(x)
      (x-mu) %*% t(x-mu)
    
    theta <- as.matrix(theta)
    tmp <- apply(theta, 1, G.cov_sub)
    if (length(mu)>1)
      t(tmp)
    else
      as.matrix(tmp)
  }
outAdMitIS <- AdMitIS(KERNEL=GelmanMeng, G=G.cov, mit=outAdMit$mit, mu=c(1.459,1.459))
print(outAdMitIS)
V <- matrix(outAdMitIS$ghat,2,2)
print(V)
cov2cor(V)
wait()

## Use independence Metropolis-Hastings algorithm with
## the mixture approximation as the candidate density.
outAdMitMH <- AdMitMH(KERNEL=GelmanMeng, mit=outAdMit$mit)
print(outAdMitMH)
wait()

## Use some functions of the package 'coda' to obtain summaries from the MCMC output.
draws <- as.mcmc(outAdMitMH$draws[1001:1e5,])
colnames(draws) <- c("X1","X2")
summary(draws)$stat
summary(draws)$stat[,3]^2 / summary(draws)$stat[,4]^2
plot(draws)
wait()
    

  
    
  
