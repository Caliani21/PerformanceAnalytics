###############################################################################
# Functions to perform component risk calculations on portfolios of assets.
#
# Copyright (c) 2004-2020 Kris Boudt and Brian G. Peterson
# This R package is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
###############################################################################
# $Id$
###############################################################################


.setalphaprob = function (p)
{
  if ( p >= 0.51 ) {
    # looks like p was a percent like .99
    alpha <- 1-p
  } else {
    alpha = p
  }
  
  return (alpha)
}

pvalJB = function(R)
{
  m2 = centeredmoment(R,2)
  m3 = centeredmoment(R,3)
  m4 = centeredmoment(R,4)
  skew = .skewEst(m3, m2)
  exkur = .kurtEst(m4, m2)
  JB = ( length(R)/6 )*( skew^2 + (1/4)*(exkur^2) )
  out = 1-pchisq(JB,df=2)
}

VaR.Gaussian =  function(R,p)
{
  alpha = .setalphaprob(p)
  
  columns = ncol(R)
  for(column in 1:columns) {
    r = as.vector(na.omit(R[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric") 
    # location = apply(R,2,mean);
    m2 = centeredmoment(r,2)
    VaR = - mean(r) - qnorm(alpha)*sqrt(m2)
    VaR=array(VaR)
    if (column==1) {
      #create data.frame
      result=data.frame(VaR=VaR)
    } else {
      VaR=data.frame(VaR=VaR)
      result=cbind(result,VaR)
    }
  }
  colnames(result)<-colnames(R)
  return(result)
}

ES.Gaussian =  function(R,p)
{
  alpha = .setalphaprob(p)
  columns = ncol(R)
  for(column in 1:columns) {
    r = as.vector(na.omit(R[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric") 
    # location = apply(R,2,mean);
    m2 = centeredmoment(r,2)
    GES = - mean(r) + dnorm(qnorm(alpha))*sqrt(m2)/alpha
    GES=array(GES)
    if (column==1) {
      #create data.frame
      result=data.frame(GES=GES)
    } else {
      GES=data.frame(GES=GES)
      result=cbind(result,GES)
    }
  }
  colnames(result)<-colnames(R)
  return(result)
}

VaR.CornishFisher =  function(R,p)
{
  alpha = .setalphaprob(p)
  z = qnorm(alpha)
  columns = ncol(R)
  for(column in 1:columns) {
    r = as.vector(na.omit(R[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric") 
    
    # location = apply(r,2,mean);
    m2 = centeredmoment(r,2)
    m3 = centeredmoment(r,3)
    m4 = centeredmoment(r,4)
    skew = .skewEst(m3, m2)
    exkurt = .kurtEst(m4, m2)
    
    h = z + (1/6)*(z^2 -1)*skew + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2
    
    
    VaR = - mean(r) - h*sqrt(m2)
    VaR=array(VaR)
    if (column==1) {
      #create data.frame
      result=data.frame(VaR=VaR)
    } else {
      VaR=data.frame(VaR=VaR)
      result=cbind(result,VaR)
    }
  }
  colnames(result)<-colnames(R)
  return(result)
}

ES.CornishFisher =  function(R,p,c=2)
{
  alpha = .setalphaprob(p)
  p = alpha
  z = qnorm(alpha)
  
  columns = ncol(R)
  for(column in 1:columns) {
    r = as.vector(na.omit(R[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric") 
    # location = apply(R,2,mean);
    m2 = centeredmoment(r,2)
    m3 = centeredmoment(r,3)
    m4 = centeredmoment(r,4)
    skew = .skewEst(m3, m2)
    exkurt = .kurtEst(m4, m2)
    
    h = z + (1/6)*(z^2 -1)*skew
    if(c==2){ h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2};
    
    # MES = dnorm(h)
    # MES = MES + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
    # MES = MES +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
    # MES = MES + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
    MES <- dnorm(h) * (1 + h^3 * skew / 6.0 +
                         (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew^2 / 72 +
                         (h^4 - 2 * h^2 - 1) * exkurt / 24)
    MES = - mean(r) + (sqrt(m2)/p)*MES
    MES=array(MES)
    if (column==1) {
      #create data.frame
      result=data.frame(MES=MES)
    } else {
      MES=data.frame(MES=MES)
      result=cbind(result,MES)
    }
  }
  colnames(result)<-colnames(R)
  return(result)
}

operES.CornishFisher =  function(R,p,c=2)
{
  alpha = .setalphaprob(p)
  p = alpha
  z = qnorm(alpha)
  columns = ncol(R)
  for(column in 1:columns) {
    r = as.vector(na.omit(R[,column]))
    if (!is.numeric(r)) stop("The selected column is not numeric") 
    #location = apply(R,2,mean);
    m2 = centeredmoment(r,2)
    m3 = centeredmoment(r,3)
    m4 = centeredmoment(r,4)
    skew = .skewEst(m3, m2)
    exkurt = .kurtEst(m4, m2)
    
    h = z + (1/6)*(z^2 -1)*skew
    if(c==2){ h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2};
    
    # MES = dnorm(h)
    # MES = MES + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
    # MES = MES +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
    # MES = MES + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
    MES <- dnorm(h) * (1 + h^3 * skew / 6.0 +
                         (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew^2 / 72 +
                         (h^4 - 2 * h^2 - 1) * exkurt / 24)
    MES = - mean(r) - (sqrt(m2))*min( -MES/alpha , h )
    MES=array(MES)
    if (column==1) {
      #create data.frame
      result=data.frame(MES=MES)
    } else {
      MES=data.frame(MES=MES)
      result=cbind(result,MES)
    }
  }	
  colnames(result)<-colnames(R)
  return(result)
}

# Definition of statistics needed to compute Gaussian and modified VaR and ES for the return R of portfolios
# and to compute the contributions to portfolio downside risk, made by the different positions in the portfolio.
#----------------------------------------------------------------------------------------------------------------


#' Portfolio moments
#' 
#' computes the portfolio second, third and fourth central moments given the 
#' multivariate comoments and the portfolio weights. The gradient functions
#' compute the gradient of the portfolio central moment with respect to the
#' portfolio weights, evaluated in the portfolio weights.
#' 
#' For documentation on the coskewness and cokurtosis matrices, we refer to
#' ?CoMoments. Both the full matrices and reduced form can be the used as
#' input for the function related to the portfolio third and fourth central
#' moments.
#' 
#' @name portfolio-moments
#' @concept co-moments
#' @concept moments
#' @param w vector of length p with portfolio weights
#' @param sigma portfolio covariance matrix of dimension p x p
#' @param M3 matrix of dimension p x p^2, or a vector with 
#' (p * (p + 1) * (p + 2) / 6) unique coskewness elements
#' @param M4 matrix of dimension p x p^3, or a vector with 
#' (p * (p + 1) * (p + 2) * (p + 3) / 12) unique coskewness elements
#' @author Kris Boudt, Peter Carl, Dries Cornilly, Brian Peterson
#' @seealso \code{\link{CoMoments}} \cr \code{\link{ShrinkageMoments}} \cr \code{\link{EWMAMoments}}
#' \cr \code{\link{StructuredMoments}} \cr \code{\link{MCA}}
###keywords ts multivariate distribution models
#' @examples
#' 
#' data(managers)
#' 
#' # equal weighted portfolio of managers
#' p <- ncol(edhec)
#' w <- rep(1 / p, p)
#' 
#' # portfolio variance and its gradient with respect to the portfolio weights
#' sigma <- cov(edhec)
#' pm2 <- portm2(w, sigma)
#' dpm2 <- derportm2(w, sigma)
#' 
#' # portfolio third central moment and its gradient with respect to the portfolio weights
#' m3 <- M3.MM(edhec)
#' pm3 <- portm3(w, m3)
#' dpm3 <- derportm3(w, m3)
#' 
#' # portfolio fourth central moment and its gradient with respect to the portfolio weights
#' m4 <- M4.MM(edhec)
#' pm4 <- portm4(w, m4)
#' dpm4 <- derportm4(w, m4)
#' 
#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
portm2 <- function(w, sigma)
{
  return(as.numeric(t(w)%*%sigma%*%w)) #t(w) for first item?
}

#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
derportm2 <- function(w, sigma)
{
  return(2*sigma%*%w)
}

#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
portm3 <- function(w, M3)
{
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- M3.mat2vec(M3)
  .Call('M3port', w, M3, length(w), PACKAGE="PerformanceAnalytics")
}

#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
derportm3 <- function(w, M3)
{
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- M3.mat2vec(M3)
  as.matrix(.Call('M3port_grad', w, M3, length(w), PACKAGE="PerformanceAnalytics"), ncol = 1)
}

#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
portm4 <- function(w, M4)
{
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- M4.mat2vec(M4)
  .Call('M4port', w, M4, length(w), PACKAGE="PerformanceAnalytics")
}

#'@useDynLib PerformanceAnalytics
#'@export
#'@rdname portfolio-moments
derportm4 <- function(w, M4)
{
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- M4.mat2vec(M4)
  as.matrix(.Call('M4port_grad', w, M4, length(w), PACKAGE="PerformanceAnalytics"), ncol = 1)
}

Portmean = function(w,mu)
{
  return( list( t(w)%*%mu , as.vector(w)*as.vector(mu) , as.vector(w)*as.vector(mu)/t(w)%*%mu) )
}

Portsd =  function(w,sigma)
{
  pm2 = portm2(w,sigma)
  dpm2 = derportm2(w,sigma)
  dersd = (0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = dersd*as.vector(w)
  names(contrib) = names(w)
  pct_contrib = contrib/sqrt(pm2)
  names(pct_contrib) = names(w)
  # check
  if( abs( sum(contrib)-sqrt(pm2))>0.01*sqrt(pm2)) { print("error") 
  } else {
    ret<-list(  sqrt(pm2) , contrib , pct_contrib )
    names(ret) <- c("StdDev","contribution","pct_contrib_StdDev")
  }
  return(ret)
}

VaR.Gaussian.portfolio =  function(p,w,mu,sigma)
{
  alpha = .setalphaprob(p)
  p=alpha
  location = sum(w * mu)
  pm2 = portm2(w,sigma)
  dpm2 = derportm2(w,sigma)
  VaR = - location - qnorm(alpha)*sqrt(pm2)
  derVaR = - as.vector(mu)- qnorm(alpha)*(0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = derVaR*as.vector(w)
  names(contrib) = names(w)
  pct_contrib = contrib/VaR
  names(pct_contrib) = names(w)
  if( abs( sum(contrib)-VaR)>0.01*abs(VaR)) { stop("contribution does not add up") } else {
    ret<-list( VaR  ,  contrib  , pct_contrib ) 
    names(ret) = c("VaR","contribution","pct_contrib_VaR")
  }
  return(ret)
  
}

kernel = function( x , h )
{
  return( apply( cbind( rep(0,length(x)) , 1-abs(x/h) ) , 1 , 'max' ) );
}

VaR.kernel.portfolio =  function( R, p, w )
{
  alpha = .setalphaprob(p)
  T = dim(R)[1]; N = dim(R)[2];
  portfolioreturn = c();
  for( t in 1:T ){ portfolioreturn = c( portfolioreturn , sum(w*R[t,]) ) }
  bandwith = 2.575*sd(portfolioreturn)/(T^(1/5)) ;
  CVaR = c();
  VaR = -quantile( portfolioreturn , probs = alpha );
  weights = kernel(x= (-VaR-portfolioreturn) , h=bandwith);
  correction  = VaR/sum(weights*portfolioreturn)
  for( i in 1:N ){ CVaR = c( CVaR , sum( weights*R[,i] ) ) }
  # check
  CVaR = w*CVaR
  #print( sum(CVaR) ) ; print( sum( weights*portfolioreturn)  )
  CVaR = CVaR/sum(CVaR)*VaR
  pct_contrib = CVaR/VaR
  names(CVaR)<-colnames(R)
  names(pct_contrib)<-colnames(R)
  ret= list( VaR  ,  CVaR  , pct_contrib  )
  names(ret) = c("VaR","contribution","pct_contrib_VaR")
  return(ret)
}

ES.kernel.portfolio= function( R, p, w )
{#WARNING incomplete
  VAR<-VaR.kernel.portfolio( R, p, w )
  
  #I'm sure that using Return.portfolio probably makes more sense here...
  T = dim(R)[1]; N = dim(R)[2];
  portfolioreturn = c();
  for( t in 1:T ){ portfolioreturn = c( portfolioreturn , sum(w*R[t,]) ) }
  
  PES<-mean(portfolioreturn>VAR$VaR)
  
  
}

ES.Gaussian.portfolio =  function(p,w,mu,sigma)
{
  alpha = .setalphaprob(p)
  p=alpha
  location = sum(w * mu)
  pm2 = portm2(w,sigma)
  dpm2 = derportm2(w,sigma)
  ES = - location + dnorm(qnorm(alpha))*sqrt(pm2)/alpha
  derES = - as.vector(mu) + (1/p)*dnorm(qnorm(alpha))*(0.5*as.vector(dpm2))/sqrt(pm2);
  contrib = as.vector(w)*derES;
  names(contrib) = names(w)
  pct_contrib = contrib/ES
  names(pct_contrib) = names(w)
  if( abs( sum(contrib)-ES)>0.01*abs(ES)) { stop("contribution does not add up") } 
  else {
    ret = list(  ES  ,  contrib ,  pct_contrib  )
    names(ret) = c("ES","contribution","pct_contrib_ES")
    return(ret)
  }
}

VaR.CornishFisher.portfolio = function(p,w,mu,sigma,M3,M4)
{
  alpha = .setalphaprob(p)
  p=alpha
  z = qnorm(alpha)
  location = sum(w * mu)
  pm2 = portm2(w,sigma)
  dpm2 = as.vector( derportm2(w,sigma) )
  pm3 = portm3(w,M3)
  dpm3 = as.vector( derportm3(w,M3) )
  pm4 = portm4(w,M4)
  dpm4 = as.vector( derportm4(w,M4) )
  
  skew = .skewEst(pm3, pm2)
  exkurt = .kurtEst(pm4, pm2)
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  
  MVaR = - location - h*sqrt(pm2)
  
  # derGausVaR = - as.vector(mu)- qnorm(alpha)*(0.5*as.vector(dpm2))/sqrt(pm2);
  # derMVaR = derGausVaR + (0.5*dpm2/sqrt(pm2))*( -(1/6)*(z^2 -1)*skew  - (1/24)*(z^3 - 3*z)*exkurt + (1/36)*(2*z^3 - 5*z)*skew^2 )
  # derMVaR = derMVaR + sqrt(pm2)*( -(1/6)*(z^2 -1)*derskew  - (1/24)*(z^3 - 3*z)*derexkurt + (1/36)*(2*z^3 - 5*z)*2*skew*derskew  )
  derh = (z^2 - 1.0) * derskew / 6 + (z^3 - 3 * z) * derexkurt / 24 - (2 * z^3 - 5 * z) * skew * derskew / 18;
  derMVaR = -as.vector(mu) - h * as.vector(dpm2) / (2.0 * sqrt(pm2)) - sqrt(pm2) * derh;
  
  contrib = as.vector(w)*as.vector(derMVaR)
  pct_contrib = contrib/MVaR
  names(contrib) <- names(w)
  names(pct_contrib) <- names(w)
  if( abs( sum(contrib)-MVaR)>0.01*abs(MVaR)) { stop("contribution does not add up") } 
  else {
    ret=(list(   MVaR  ,  contrib, pct_contrib  ) )
    names(ret) = c("MVaR","contribution","pct_contrib_MVaR")
    return(ret)
  }
}


derIpower = function(power,h)
{
  # probably redundant now
  fullprod = 1;
  
  if( (power%%2)==0 ) #even number: number mod is zero
  {
    pstar = power/2;
    for(j in c(1:pstar))
    {
      fullprod = fullprod*(2*j)
    }
    I = -fullprod*h*dnorm(h);
    
    for(i in c(1:pstar) )
    {
      prod = 1;
      for(j in c(1:i) )
      {
        prod = prod*(2*j)
      }
      I = I + (fullprod/prod)*(h^(2*i-1))*(2*i-h^2)*dnorm(h)
    }
  }else{
    pstar = (power-1)/2
    for(j in c(0:pstar) )
    {
      fullprod = fullprod*( (2*j)+1 )
    }
    I = -fullprod*dnorm(h);
    
    for(i in c(0:pstar) )
    {
      prod = 1;
      for(j in c(0:i) )
      {
        prod = prod*( (2*j) + 1 )
      }
      I = I + (fullprod/prod)*(h^(2*i)*(2*i+1-h^2) )*dnorm(h)
    }
  }
  return(I)
}

CF_moments_RSM <- function(S_raw, K_raw) {
  S <- abs(S_raw)
  K <- K_raw
  
  if (S < 1e-6) {
    return(list(Sc = 0, Kc = K, Case = "S ~ 0 (No Correction)"))
  } 
  
  if (K >= 5 && K <= 40) {
    if (S >= 0.5 && S <= 2.2) {
      Kc <- -5.962 + 21.53*S^0.5 - 1.548*K^0.5 - 26.52*S + 1.820*K + 11.08*S^1.5 - 0.0443*K^1.5 - 2.564*S^0.5*K + 5.739*S*K^0.5 + 0.342*S^2 + 0.00162*K^2 - 3.773*S^1.5*K^0.5 + 0.880*S*K + 0.0328*S^0.5*K^1.5 + 0.000901*S*K^2 + 0.0717*S^2*K - 0.0216*S^1.5*K^1.5 - 0.721*log(S)*log(K) + 0.349*log(S)*K + 0.0928*S*log(K) + 0.366*S^-1 - 0.555*K^-1
      Sc <- -1.816 + 6.812*S^0.5 - 0.577*K^0.5 - 8.636*S + 0.508*K + 4.235*S^1.5 - 0.00685*K^1.5 - 0.848*S^0.5*K + 2.671*S*K^0.5 - 0.0969*S^2 - 0.000304*K^2 - 1.259*S^1.5*K^0.5 + 0.226*S*K + 0.0191*S^0.5*K^1.5 + 0.000196*S*K^2 + 0.0249*S^2*K - 0.00666*S^1.5*K^1.5 - 0.105*log(S)*log(K) + 0.0987*log(S)*K - 0.845*S*log(K) + 0.135*S^-1 - 0.416*K^-1
      case <- "Case 1"
    } else {
      Kc <- 0.0832 + 0.0451*S^0.5 + 0.732*K^0.5 - 0.601*S + 0.124*K + 0.396*S^0.5*K^0.5 + 1.261*S^1.5 - 0.0195*K^1.5 - 0.0704*S^0.5*K - 0.528*S*K^0.5 - 0.198*S^2 + 0.00181*K^2 - 0.122*S^1.5*K^0.5 + 0.0836*S*K + 0.000231*S^0.5*K^1.5 + 9.56e-05*S*K^2 + 0.0133*S^2*K - 0.00373*S^1.5*K^1.5 - 0.0305*log(S)*log(K) + 0.0029*log(S)*K + 0.24*S*log(K) - 0.000296*S^-1 - 0.444*K^-1
      Sc <- -0.0189 + 0.161*S^0.5 + 0.0215*K^0.5 + 0.453*S + 0.00139*K - 0.0862*S^0.5*K^0.5 + 0.326*S^1.5 - 8.51e-06*K^1.5 - 0.00168*S^0.5*K + 0.23*S*K^0.5 - 0.0136*S^2 - 2.32e-06*K^2 - 0.129*S^1.5*K^0.5 - 0.000326*S*K - 0.000151*S^0.5*K^1.5 + 4.93e-05*S*K^2 + 0.00662*S^2*K - 0.000649*S^1.5*K^1.5 + 0.00396*log(S)*log(K) + 0.000457*log(S)*K - 0.221*S*log(K) + 0.000228*S^-1 - 0.025*K^-1
      case <- "Case 2"
    }
  } else { 
    if (S >= 0.5) { 
      Kc <- 1.749 - 6.604*K^0.5 + 3.425*S + 1.313*K + 7.491*S^0.5*K^0.5 - 11.83*S^1.5 - 0.858*K^1.5 + 9.011*S^2 + 0.141*K^2 - 3.346*S^1.5*K^0.5 + 0.638*S*K + 0.11*S^0.5*K^1.5 - 0.124*S*K^2 - 0.642*S^2*K + 0.499*S^1.5*K^1.5 - 0.517*log(S)*log(K) - 0.65*log(S)*K + 0.834*S*log(K) + 0.136*S^-1 + 0.0989*K^-1
      Sc <- 2.111 - 3.498*K^0.5 - 2.87*S - 0.123*K + 3.836*S^0.5*K^0.5 + 2.956*S^1.5 - 0.162*K^1.5 + 2.008*S^2 + 0.037*K^2 - 4.884*S^1.5*K^0.5 + 1.72*S*K - 0.153*S^0.5*K^1.5 - 0.00138*S*K^2 + 0.239*S^2*K - 0.0883*S^1.5*K^1.5 - 0.227*log(S)*log(K) - 0.436*log(S)*K + 0.7*S*log(K) - 0.0739*S^-1 + 0.0414*K^-1
      case <- "Case 3"
    } else if (S >= 0.25) { 
      Kc <- -1.612 + 1.894*S^0.5 + 1.938*K^0.5 + 0.273*K - 1.018*S^0.5*K^0.5 - 4.22*S^1.5 - 0.141*K^1.5 + 2.164*S^2 + 0.0247*K^2 + 2.786*S^1.5*K^0.5 - 0.454*S*K + 0.0381*S^0.5*K^1.5 - 0.0392*S*K^2 - 0.862*S^2*K + 0.307*S^1.5*K^1.5 + 0.103*log(S)*log(K) + 0.0341*log(S)*K - 0.481*S*log(K) + 0.0164*S^-1 - 0.00817*K^-1
      Sc <- 0.172 + 0.132*S^0.5 - 0.296*K^0.5 - 0.0415*K + 0.346*S^0.5*K^0.5 + 1.491*S^1.5 - 0.0327*K^1.5 + 0.134*S^2 + 0.00278*K^2 - 1.33*S^1.5*K^0.5 + 0.249*S*K + 0.0333*S^0.5*K^1.5 - 0.00129*S*K^2 + 0.205*S^2*K - 0.0597*S^1.5*K^1.5 - 0.0109*log(S)*log(K) - 0.0507*log(S)*K + 0.114*S*log(K) - 0.00419*S^-1 + 0.00152*K^-1
      case <- "Case 4"
    } else { 
      Kc <- -0.304 + 0.743*S^0.5 + 0.597*K^0.5 - 1.662*S + 0.676*K - 1.073*S^0.5*K^0.5 + 0.226*S^1.5 - 0.299*K^1.5 + 0.49*S^0.5*K + 2.314*S*K^0.5 + 0.463*S^2 + 0.0432*K^2 - 0.234*S^1.5*K^0.5 - 0.891*S*K - 0.0254*S^0.5*K^1.5 - 0.00616*S*K^2 - 0.272*S^2*K + 0.205*S^1.5*K^1.5 + 0.00942*log(S)*log(K) - 0.00642*log(S)*K - 0.164*S*log(K) - 0.0000209*S^-1 + 0.00151*K^-1
      Sc <- 0.00512 - 0.024*S^0.5 - 0.00778*K^0.5 + 1.277*S + 0.00499*K + 0.0386*S^0.5*K^0.5 - 0.114*S^1.5 - 0.000479*K^1.5 - 0.0336*S^0.5*K - 0.483*S*K^0.5 + 0.265*S^2 - 0.000052*K^2 - 0.0857*S^1.5*K^0.5 + 0.109*S*K + 0.00708*S^0.5*K^1.5 - 0.00487*S*K^2 - 0.0332*S^2*K + 0.0161*S^1.5*K^1.5 - 0.00027*log(S)*log(K) + 0.000262*log(S)*K + 0.0513*S*log(K) + 0.000000429*S^-1 + 0.00011*K^-1
      case <- "Case 5"
    }
  }
  return(list(Sc = sign(S_raw) * Sc, Kc = Kc, Case = case))
}

ES.CornishFisher.portfolio =  function(p,w,mu,sigma,M3,M4)
{
  alpha = .setalphaprob(p)
  p=alpha
  z = qnorm(alpha)
  location = sum(w * mu)
  pm2 = portm2(w,sigma)
  dpm2 = as.vector( derportm2(w,sigma) )
  pm3 = portm3(w,M3)
  dpm3 = as.vector( derportm3(w,M3) )
  pm4 = portm4(w,M4)
  dpm4 = as.vector( derportm4(w,M4) )

  new_moments <- CF_moments_RSM(S_raw = .skewEst(pm3, pm2), K_raw = .kurtEst(pm4, pm2))
  skew = new_moments$Sc # .skewEst(pm3, pm2)
  exkurt = new_moments$Kc # .kurtEst(pm4, pm2)
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  # E = dnorm(h)
  # E = E + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  # E = E +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  # E = E + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  # E = E/alpha
  E <- dnorm(h) * (1 + h^3 * skew / 6.0 +
                     (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew^2 / 72 +
                     (h^4 - 2 * h^2 - 1) * exkurt / 24) / alpha
  MES = - location + sqrt(pm2)*E
  
  # derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
  # derE = (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*derexkurt
  # derE = derE +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*derskew
  # derE = derE + (1/36)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*skew*derskew
  # X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt
  # X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew
  # X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
  # derE = derE+derh*X  # X is a scalar
  # derE = derE/alpha
  # derMES = derMES + sqrt(pm2)*derE
  derMES = -mu + E * dpm2 / (2 * sqrt(pm2)) - sqrt(pm2) * E * h * derh + dnorm(h) * sqrt(pm2) / alpha *
    (derh * (h^2 * skew / 2 + (6.0 * h^5 - 36 * h^3 + 18 * h) * skew * skew / 72 +
               (4 * h^3 - 4 * h) * exkurt / 24) + h^3 * derskew / 6 +
       (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew * derskew / 36 + (h^4 - 2 * h^2 - 1) * derexkurt / 24)
  contrib = as.vector(w)*as.vector(derMES)
  names(contrib) = names(w)
  pct_contrib = contrib/MES
  names(pct_contrib) = names(w)
  if( abs( sum(contrib)-MES)>0.01*abs(MES)) { stop("contribution does not add up") } 
  else {
    ret= list(   MES , contrib , pct_contrib) 
    names(ret) = c("MES","contribution","pct_contrib_MES")
    return(ret)
  }
}


operES.CornishFisher.portfolio = function(p,w,mu,sigma,M3,M4)
{
  alpha = .setalphaprob(p)
  p=alpha
  z = qnorm(alpha)
  location = sum(w * mu)
  pm2 = portm2(w,sigma)
  dpm2 = as.vector( derportm2(w,sigma) )
  pm3 = portm3(w,M3)
  dpm3 = as.vector( derportm3(w,M3) )
  pm4 = portm4(w,M4)
  dpm4 = as.vector( derportm4(w,M4) )
  
  skew = .skewEst(pm3, pm2)
  exkurt = .kurtEst(pm4, pm2)
  
  derskew = ( 2*(pm2^(3/2))*dpm3 - 3*pm3*sqrt(pm2)*dpm2 )/(2*pm2^3)
  derexkurt = ( (pm2)*dpm4 - 2*pm4*dpm2    )/(pm2^3)
  
  h = z + (1/6)*(z^2 -1)*skew
  h = h + (1/24)*(z^3 - 3*z)*exkurt - (1/36)*(2*z^3 - 5*z)*skew^2;
  
  derh = (1/6)*(z^2 -1)*derskew + (1/24)*(z^3 - 3*z)*derexkurt - (1/18)*(2*z^3 - 5*z)*skew*derskew
  
  # E = dnorm(h)
  # E = E + (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*exkurt
  # E = E +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*skew;
  # E = E + (1/72)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*(skew^2)
  # E = E/alpha
  E <- dnorm(h) * (1 + h^3 * skew / 6.0 +
                     (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew^2 / 72 +
                     (h^4 - 2 * h^2 - 1) * exkurt / 24) / alpha
  MES = - location - sqrt(pm2)*min(-E,h)
  
  if(-E<=h){
    # derMES = -mu + 0.5*(dpm2/sqrt(pm2))*E
    # derE = (1/24)*(   Ipower(4,h) - 6*Ipower(2,h) + 3*dnorm(h)   )*derexkurt
    # derE = derE +  (1/6)*(   Ipower(3,h) - 3*Ipower(1,h)   )*derskew
    # derE = derE + (1/36)*(  Ipower(6,h) -15*Ipower(4,h)+ 45*Ipower(2,h) - 15*dnorm(h) )*skew*derskew
    # X = -h*dnorm(h) + (1/24)*(  derIpower(4,h) - 6*derIpower(2,h) -3*h*dnorm(h)  )*exkurt
    # X = X + (1/6)*( derIpower(3,h) - 3*derIpower(1,h) )*skew
    # X = X + (1/72)*( derIpower(6,h) - 15*derIpower(4,h) + 45*derIpower(2,h) + 15*h*dnorm(h)  )*skew^2
    # derE = derE+derh*X  # X is a scalar
    # derE = derE/alpha
    # derMES = derMES + sqrt(pm2)*derE
    derMES = -mu + E * dpm2 / (2 * sqrt(pm2)) - sqrt(pm2) * E * h * derh + dnorm(h) * sqrt(pm2) / alpha *
      (derh * (h^2 * skew / 2 + (6.0 * h^5 - 36 * h^3 + 18 * h) * skew * skew / 72 +
                 (4 * h^3 - 4 * h) * exkurt / 24) + h^3 * derskew / 6 +
         (h^6 - 9 * h^4 + 9 * h^2 + 3) * skew * derskew / 36 + (h^4 - 2 * h^2 - 1) * derexkurt / 24)
  } else {
    derMES = -mu - 0.5*(dpm2/sqrt(pm2))*h - sqrt(pm2)*derh
  }
  contrib = as.vector(w)*as.vector(derMES)
  names(contrib) = names(w)
  pct_contrib = contrib/MES
  names(pct_contrib) = names(w)
  if( abs( sum(contrib)-MES)>0.01*abs(MES)) { stop("contribution does not add up") } 
  else {
    ret= list(   MES , contrib , pct_contrib) 
    names(ret) = c("MES","contribution","pct_contrib_MES")
    return(ret)
  }
}

ES.historical = function(R,p) {
  alpha = .setalphaprob(p)
  for(column in 1:ncol(R)) {
    r = na.omit(as.vector(R[,column]))
    q = quantile(r,probs=alpha)
    exceedr = r[r<q]
    hES = (-mean(exceedr))
    if(is.nan(hES)){
      warning(paste(colnames(R[,column]),"No values less than VaR observed.  Setting ES equal to VaR."))
      hES=-q
    }
    
    hES=array(hES)
    if (column==1) {
      #create data.frame
      result=data.frame(hES=hES)
    } else {
      hES=data.frame(hES=hES)
      result=cbind(result,hES)
    }
  }	
  colnames(result)<-colnames(R)
  return(result)
}

ES.historical.portfolio = function(R,p,w)
{
  hvar = as.numeric(VaR.historical.portfolio(R,p,w)$hVaR)
  T = dim(R)[1]
  N = dim(R)[2]
  c_exceed = 0;
  r_exceed = 0;
  realizedcontrib = rep(0,N);
  for(t in c(1:T) )
  {
    rt = as.vector(R[t,])
    rp = sum(w*rt)
    if(rp<= -hvar){
      c_exceed = c_exceed + 1;
      r_exceed = r_exceed + rp;
      for( i in c(1:N) ){
        realizedcontrib[i] =realizedcontrib[i] + w[i]*rt[i] }
    }
  }
  pct.contrib=as.numeric(realizedcontrib)/r_exceed ;
  names(pct.contrib)<-names(w)
  # TODO construct realized contribution
  
  ret <- list(-r_exceed/c_exceed,c_exceed,pct.contrib) 
  names(ret) <- c("-r_exceed/c_exceed","c_exceed","pct_contrib_hES")
  return(ret)
}

VaR.historical = function(R,p)
{
  alpha = .setalphaprob(p)
  for(column in 1:ncol(R)) {
    r = na.omit(as.vector(R[,column]))
    rq = -quantile(r,probs=alpha)
    if (column==1) {
      result=data.frame(rq=rq)
    } else {
      rq=data.frame(rq=rq)
      result=cbind(result,rq)
    }
  }
  colnames(result)<-colnames(R)
  return(result)
} 

VaR.historical.portfolio = function(R,p,w=NULL)
{
  alpha = .setalphaprob(p)
  rp = Return.portfolio(R,w, contribution=TRUE)
  hvar = -quantile(zerofill(rp$portfolio.returns),probs=alpha)
  names(hvar) = paste0('hVaR ',100*(1-alpha),"%")
  
  # extract negative periods, 
  zl<-rp[rp$portfolio.returns<0,]
  
  # and construct weighted contribution
  zl.contrib <- colMeans(zl)[-1]
  ratio <- -hvar/sum(colMeans(zl))
  zl.contrib <- ratio * zl.contrib
  
  # and construct percent contribution
  zl.pct.contrib <- (1/sum(zl.contrib))*zl.contrib
  
  ret=list(hVaR = hvar,
           contribution = zl.contrib,
           pct_contrib_hVaR = zl.pct.contrib)
  
  return(ret)
}

.skewEst <- function(m3, m2) {
  if (isTRUE(all.equal(m2, 0.0))) {
    return(0.0)
  } else {
    return(m3 / m2^(3/2))
  }
}

.kurtEst <- function(m4, m2) {
  if (isTRUE(all.equal(m2, 0.0))) {
    return(0.0)
  } else {
    return(m4 / m2^2 - 3)
  }
}

###############################################################################
# R (https://r-project.org/) Econometrics for Performance and Risk Analysis
#
# Copyright (c) 2004-2020 Peter Carl and Brian G. Peterson and Kris Boudt
#
# This R package is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id$
#
###############################################################################
