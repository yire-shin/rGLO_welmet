
library(ismev)
library(lmomco)
library(EnvStats) 

#------------------------------------------------------------------------------
test.glo = welmet.rglo(xdat= bevern[,2:4], numr=3, alpha.qm=0.5, w.Hlme1=0.5, 
                       w.Hqm= 0.5, numB= 200, qqt= c(.98,.99,.995))

# --- main program -------------------------------------------               

welmet.rglo = function(xdat, numr=NULL, alpha.qm=0.5, w.Hlme1=0.5, w.Hqm= 0.5,
                       numB= 500, qqt= c(.95,.99,.995), Sinv, com.Sinv=T)
{

  z=list(); mle.rglo=list()
  numq=length(qqt)
  if( is.null(numr) ) numr = dim(xdat)[2]
  
  nsample=nrow(xdat)  # nrow(xdat) = sample size n
                      # ncol(xdat) = r
                         
# ++++++ lme1 +++++++++++++
  
  lme1 = lme1.glo(datr=xdat, qqt=qqt)
  
  z$lme1.rl = lme1$lme1.rl
  z$lme1.theta = lme1$lme1.theta
  theta.BM = lme1$lme1.theta
  
# --- Parametric bootstrap for S inverse ------
  
    if(com.Sinv==T){
  
    theta =matrix(NA, numB,3); lme.rl =matrix(NA, numB,numq)
#    Bid=seq(1,nsample)

    for (ib in 1:numB){
      
#      sam.id = sample(Bid, size=nrow(xdat), replace=TRUE)

      Bsam = gen.rglo.park(par= z$lme1.theta, sim_r=numr, sim_n=nsample)
      
#      Bsam = xdat[sam.id,]

      theta[ib,1:3] = lme1.glo(datr=Bsam, qqt=qqt)$lme1.theta
      lme.rl[ib,1:numq] = quaglo(qqt, vec2par(theta[ib,1:3],'glo'))
    }
  
    Hth = cov(theta)
    if( det(Hth) <= 0 ){
      cat("trouble in cov of theta.lme1","\n")
    }
    Sinv = solve(Hth)
    z$PBse.lme1.theta = sqrt( c(Hth[1,1], Hth[2,2], Hth[3,3]) )
    
    Hrl = cov(lme.rl)
    for (kq in 1:numq){
      z$PBse.lme1.rl[kq] = sqrt(Hrl[kq,kq])
    }
    
    }  # if com.Sinv

# +++++ welmet_QM ++++++++++
    
  wel.qm = QM.trsf2.glo(xdat, numr=numr, alpha.qm=alpha.qm, 
                        quant=qqt, theta.BM=theta.BM, Sinv=Sinv)

  z$welmet = wel.qm

# +++++++  rmle +++++++++  

    mle.rglo = rglo.fit.park(xdat, r=numr, num_inits = 30, show=F)
    
    if( (mle.rglo$conv != 0) | abs(mle.rglo$mle[3]) >= 1.0 )  {
      
      mle.rglo = rglo.fit.park(xdat, r=numr, num_inits = 40, 
                               show=F, method="L-BFGS-B")
    }

    z$rmle.rl = quaglo(qqt, vec2par(mle.rglo$mle,'glo'))
    z$rmle.theta = mle.rglo$mle

# ++++  mle1 +++++++++++++  

    mle.1glo = rglo.fit.park(xdat, r=1, num_inits = 30, show=F)
    
    if( (mle.rglo$conv != 0) | abs(mle.rglo$mle[3]) >= 1.0 )  {
      
      mle.1glo = rglo.fit.park(xdat, r=1, num_inits = 40, 
                               show=F, method="L-BFGS-B")
    }
    
    # if( mle.1glo$mle[3] <= -1.0)  mle.1glo$mle[3] = -.999
    # if( mle.1glo$mle[3] >= 1.0)  mle.1glo$mle[3] = .999

     z$mle1.rl = quaglo(qqt, vec2par(mle.1glo$mle,'glo'))
     z$mle1.theta = mle.1glo$mle
  
# +++++ Hybrid of rmle and lme +++++++

    z$rl.Hlme1 = w.Hlme1*z$rmle.rl + (1-w.Hlme1)*z$lme1.rl
    z$theta.Hlme1 = w.Hlme1*z$rmle.theta + (1-w.Hlme1)*z$lme1.th
  
    z$rl.Hqm = w.Hqm*z$rmle.rl + (1-w.Hqm)*wel.qm$welmet.rl
    z$theta.Hqm = w.Hqm*z$rmle.theta + (1-w.Hqm)*wel.qm$welmet.th

  return(z)
}
#--------------------------------------------------------------
 com.rl.glo = function(xdat, td.hap, td.cbd, numr2=NULL, quant=NULL,
                  theta.BM=NULL, Sinv=NULL){

  z=list()

  numq=length(quant)
  cbd.rl = matrix(NA, numq, numr2);
  ma.rl = matrix(NA, numq, numr2)
  rk.glo = matrix(NA, nrow=numr2, ncol=3)
  ma.theta = rep(NA, 3)
  ld=matrix(NA, numr2, 3); gld= rep(NA, numr2)

#  td.hap.new = c(xdat[,1], td.hap)
  td.cbd.new = cbind(xdat[,1], td.cbd)
  
  for (kw in 1:(numr2) ) {
         kung = lmoms(td.cbd.new[,kw], nmom=3)
         
         if( are.lmom.valid(kung, checkt3t4=TRUE) == F){
           cbd.rl[,kw] = NA
           cat(" Invalid L-moms", "\n")
         }else{

           rk.glo[kw,1:3] = parglo(kung, checklmom=F)$para     # Lme

         } # end if mom.valid
         
       cbd.rl[,kw] = quaglo(quant, vec2par(rk.glo[kw,1:3],type="glo")  )

       ld[kw,1:3] = theta.BM[1:3]-rk.glo[kw,1:3]
       gld[kw] = exp( -( t(ld[kw,1:3]) %*% Sinv %*% ld[kw,1:3] )/2 )
  } # end for kw
  
  gld[1] =1.0

  id =seq(1,numr2)
  numid=length(id)

  if(numid==0){
    cat("check LME for BM data= ",  rk.glo[1,1:3],"\n")
    stop
    
  }else{
    
    wlme=rep(0, numr2)
    wlme[id] = gld[id]/sum(gld[id])
    
    z$gld = -2*log(gld)
    z$di = gld
    z$wlme = wlme
    z$welmet.rl = cbd.rl[,id] %*% wlme[id]     # rl with welmet
    z$welmet.th = wlme[id] %*% rk.glo[id,]     # theta with welmet

    wt2= rep(1/numid, numr2)
    z$ma.rl = cbd.rl[,id] %*% wt2[id]     # rl with simple average
    z$ma.theta = wt2[id] %*% rk.glo[id,]  # theta with simple average
    
    z$each.rl = cbd.rl      # rl for each component
    z$each.theta = rk.glo

  }
  return(z)
 }
#------------------------------------------------------------
QM.trsf2.glo = function(xdat, numr=NULL, alpha.qm=0.5, quant=NULL,
                  theta.BM=NULL, Sinv=NULL) { 
  
  z=list()
  alpha = alpha.qm
  nsample = nrow(xdat)
  ndim = dim(xdat)[2]

  td.fin=matrix(NA, nrow=nsample, ncol=numr-1)
  numq=length(quant)

  td.hap=NULL
  td.cbd=NULL
  r1= xdat[,1]
  
  r1.glo = parglo(lmoms(r1, nmom=3),checklmom=F)

  for (kq in 1:(numr-1)) {

    r2= xdat[,kq]
    r3= xdat[,kq+1]

# using empirical cdf 'pemp' function from EnvStat package ------------------   
    
     r2c=  pemp(r2, r2)*alpha + pemp(r2, r3)*(1-alpha)

    r2c[ which(r2c >= 0.999999) ] = 0.999
    r2c[ which(r2c <= 0.000001) ] = 0.001
    
    td.fin[,kq] = quaglo(r2c, r1.glo)
    
    td.hap= c(td.hap, td.fin[,kq])
    td.cbd= cbind(td.cbd, td.fin[,kq] )
    
  } # end for kq

    z= com.rl.glo(xdat, td.hap, td.cbd, numr2=numr, quant=quant,
              theta.BM=theta.BM, Sinv=Sinv)

  return(z)
}
#-------------------------------------------------------------
gen.rglo.park = function(par, sim_r, sim_n){
  
  umat  <- matrix(runif(sim_n *sim_r),c(sim_n, sim_r))
  umat2 =  matrix(NA, sim_n, sim_r)

  umat2[,1]=umat[,1]
  for(ir in 2:sim_r){ 
      umat2[,ir] = umat2[,ir-1]* umat[,ir]^(1/ir)
  }

  sample=  quaglo(umat2[,1:sim_r],vec2par(par,'glo') )
  
  return(sample)
}
# #-------------------------------------------------------------
#--------------------------------------------------------------
  ginit.max.glo <-function(data,ntry){
    
    n=ntry
    init <-matrix(rep(0,n*3),ncol=3)
    
    lmom_init = lmoms(data,nmom=5)
    lmom_est <- parglo(lmom_init)
    
    init[1,1]    <-lmom_est$para[1]
    init[1,2]    <-lmom_est$para[2]
    init[1,3]    <-lmom_est$para[3]
    
    maxm1=ntry; maxm2=maxm1-1
    init[2:maxm1,1] <- init[1,1]+ rnorm(n=maxm2,mean=0,sd = 5)
    init[2:maxm1,2] <- abs( init[1,2]+ rnorm(n=maxm2,mean=5,sd = 5)) +1
    init[2:maxm1,3] <- runif(n=maxm2,min= -0.5,max=0.5)
#    init[2:maxm1,2] = max(0.1, init[2:maxm1,2])
    
    return(init)
  }
#----------------------------------------------------------
lme1.glo = function (datr=NULL, qqt=NULL){
  
  z=list()
  kung= lmoms(datr[,1], nmom=3)
  
  if( are.lmom.valid(kung, checkt3t4=TRUE) == F){
    cat("----trouble in lme1.gev----","\n")
    z$lme1.rl = NA
    
  }else{
    hos= parglo(kung, checklmom=F)
    
    z$lme1.rl= quaglo(qqt, hos)
    z$lme1.theta = hos$para
    
  }
  return(z)
}
#----------------------------------------------------------------
rglo.fit.park <- function(xdat, r = dim(xdat)[2], ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                          mulink = identity, siglink = identity, shlink = identity, num_inits = 30,
                          muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE,
                          maxit = 1000, lower=c(-Inf, 0, -1.0), upper=c(Inf, Inf, 1.0),
                          method = "Nelder-Mead", ...){
  
  # method = "L-BFGS-B",
  # park changed the method to "L-BFGS-B", and new lower and upper
  
  options(digits=8)
  z <- list()
  
  # Determine the number of parameters for each component (mu, sigma, xi, h)
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  
  z$trans <- FALSE
  
  # Generate parameter names based on the length of each list
  mu_names <- if (is.null(mul)) "mu" else c("mu", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma", paste0("sigma", seq_len(npsc - 1)))
  xi_names <- if (is.null(shl)) "xi" else c("xi", paste0("xi", seq_len(npsh - 1)))
  
  # Set initial values based on L-moments of the data and user-specified predictors
  glopar <- lmomco::parglo(lmomco::lmoms(xdat[, 1]))$para
  
  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(glopar)
  
  
  # Set up the mu matrix and initial values if each component (mu, sigma, xi, h) is provided
  
  if(is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( muinit)) muinit <- glopar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if( is.null( muinit)) muinit <- c(glopar[1], rep(0, length(mul)))
  }
  if(is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( siginit)) siginit <- glopar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if( is.null( siginit)) siginit <- c(glopar[2], rep(0, length(sigl)))
  }
  if(is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( shinit)) shinit <- glopar[3]
  }  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
    if( is.null( shinit)) shinit <- c(glopar[3], rep(0, length(shl)))
  }
  
  
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  
  z1 <- as.matrix(xdat[, 1],ncol=1)
  zr <- as.matrix(xdat[, r],ncol=1)
  
  init <- c(muinit, siginit, shinit)
  names(init) <- c(mu_names, sigma_names, xi_names)
  
  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(init)
  
  for (i in 2:num_inits) {
    # Apply random noise to each component of init to create diversity in starting points
    new_init <- init + c(
      stats::rnorm(npmu, mean = 0, sd = 1),  # Random noise for mu parameters
      abs(stats::rnorm(npsc, mean = 0, sd = 1)),  # Random noise for sigma parameters
      stats::rnorm(npsh, mean = 0, sd = 0.5) # Random noise for xi parameters
    )
    init_list[[i]] <- new_init
  }
  
  # Define the log-likelihood function for the case when r = 1
  
  #------------------------------------------------------------------------------------
  glo.lik.park <- function(a) {
    
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    nsam =nrow(xdat)
    
    y <- (xdat[,1] - mu)/sc  #park changed xdat
    y <- 1 - xi * y
    
    f <- 1 + y^(1/xi)
    
    F_ <- 1-(1-1/f)
    
    if (any(y <= 0,na.rm=T) || any(sc <= 0,na.rm=T) || any(f <=0,na.rm=T) || any(F_ >1,na.rm=T))
      return(10^6)
    
    nllh = nsam*(log(sc)) + sum(log(f) * 2) + sum(log(y) * (1 - 1/xi))
    nllh[1] 
    
  }
  #----------------------------------------------------------------------------------------
  # Define the log-likelihood function for the case when r > 1
  rglo.lik.park <- function(a) {
    
    mu <- mulink(drop(mumat %*% (a[1:npmu])))
    sc <- siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    xi <- shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
    
    ri  <- (r-seq(1:(r))) # r-i
    cr  <- (1+ri)    # c_r
    nsam =nrow(xdat)
    
    # constraints 1 #
    
    if (any(sc <= 0,na.rm=T) | any(cr < 0,na.rm=T) ) return(10^6)
    
    y <- 1 - xi * (xdat[,1:r] - mu)/sc   #park changed xdat
    f <- 1 + (1 - xi * (zr - mu)/sc)^(1/xi)
    
    # constraints 2,3,4 #
    
    if (any(y<= 0, na.rm = TRUE)  || any(f<= 0, na.rm = TRUE) ) return(10^6)
    
    f2 = (1 - 1/xi) * log(y)
    
    w <-  f2 - log(cr)
    
    if(r>1){ wrs <- rowSums(w, na.rm = TRUE) 
    
    }else{
      wrs=w
    }
    
    nllh = sum((r + 1) * log(f), na.rm=T) + sum(wrs) + nsam*(r*log(sc))
    nllh[1]
    
  }
  
  # Apply optimization on each set of initial values and retain results
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      stats::optim(init, glo.lik.park, hessian=F, method=method, 
                   control=list(maxit=maxit)) # Park changed to Hessian= F
    } else {
      stats::optim(init, rglo.lik.park, hessian=F, method=method, 
                   control=list(maxit=maxit)) # Park changed to Hessian= F
    }
  })
  
  # Collect optimization results and filter out invalid results
  optim_value <- data.frame(
    num = 1:length(optim_results),
    nllh = sapply(optim_results, function(res) res$value)
    #,
    # grad = sapply(optim_results, function(res) {
    #   sum(abs(if (r == 1) numDeriv::grad(glo.lik, res$par) else numDeriv::grad(rglo.lik, res$par)))
    # })
  )
  
  
  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  #  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  best_result <- optim_results[[optim_value$num[1]]]
  
  # Extract and store the best-fit parameter values
  mu <- drop(mumat %*% (best_result$par[1:npmu]))
  sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  xi <- drop(shmat %*% (best_result$par[seq(npmu + npsc + 1, length = npsh)]))
  
  # Store results in the output list with dynamic parameter names
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle  <- best_result$par
  #  z$cov  <- solve(best_result$hessian)  # Park changed to #
  #  z$se   <- sqrt(diag(z$cov))    # Park changed to #
  z$vals <- cbind(mu, sc, xi)
  z$r    <- r
  if(show) {
    if(z$trans)
      print(z[c(2, 3, 4)])
    #else print(z[4])
    if(!z$conv)
      print(z[c(4, 5, 7, 9)])
  }
  
  class(z) <- "rglo.fit"
  invisible(z)
}
