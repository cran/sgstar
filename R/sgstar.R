#' @title Fit Seasonal Generalized Space Time Autoregressive Model
#'
#' @description sgstar function return the parameter estimation of Seaonal Generalized Space Time Autoregressive Model by using Generalized Least Square (GLS)
#'
#' @param data A dataframe that contain timeseries data with k column as space and n rows as time.
#' @param w a spatial weight, matrix ncol(data) * ncol(data) with diagonal = 0.
#' @param p an autoregressive order, value must be greater than 0.
#' @param ps an autoregressive order for seasonal, value must be greater than 0.
#' @param s an order of the seasonal period
#'
#' @return sgstar returns output with detail are shown in the following list :
#' \item{Coefficiens}{coefficiens parameter model for each location}
#' \item{Fitted.Values}{ a dataframe with fit value for each location based on model}
#' \item{Residual}{a dataframe that contain residual,that is response minus fitted values based on model}
#' \item{Performance}{a dataframe containing the following objects:}
#' \itemize{
#'    \item MSE      : Mean Squared Error (MSE) for all the data combined.
#'    \item RMSE     : Root Mean Squared Error (RMSE) for all the data combined.
#'    \item AIC      : a Version of Akaike's Information Criterion (AIC)
#'    \item Rsquared : R^2, the ‘fraction of variance explained by the model’.
#'   }
#' \item{p}{an autoregressive order}
#' \item{ps}{an autoregressive order for seasonal}
#' \item{s}{an order of the seasonal period}
#' \item{weight}{a spatial weight}
#' \item{data}{a dataset that used in modeling}
#'
#' @references Setiawan, Suhartono, and Prastuti M.(2016).S GSTAR-SUR for Seasonal Spatio Temporal Data Forecasting. Malaysian Journal Of Mathematical Sciences.10.<Corpus ID :189955959>.
#'
#' @export
#' @importFrom utils stack
#' @importFrom stats lm
#' @importFrom stats cov
#'
#' @examples
#' library(sgstar)
#' data("coords")
#' data("simulatedata")
#'
#' #create weight matrix using distance inverse matrix
#'
#' z<-dist(coords,method = "euclidean")
#' z <- as.matrix(z)
#'
#' matriksd <- 1/z
#' matriksd[is.infinite(matriksd)] <- 0
#'
#' matriksd_w <- matriksd / rowSums(as.data.frame(matriksd))
#'
#' fit <- sgstar(data = simulatedata, w = matriksd_w, p = 2,ps = 1, s =4)
#' fit
#'
#'
#'
#'
#'
#'
sgstar <- function(data,w,p,ps,s){

  y <- data

  zt<-(stack(as.data.frame(t(y)))[,1])
  k <- ncol(y)
  n <- nrow(y)

  MD<-function(zt){
    M<-matrix(0,(n*k),k)
    z=0
    for( i in 1:(n)){
      for(j in 1:k){
        z<-z+1
        M[z,j]<-zt[z]
      }
    }
    M   }

  M1 <- MD(zt)

  W1<-kronecker(diag(n), w)
  ztw<-W1%*%zt
  M2<-MD(ztw)
  zt<-as.matrix(zt)


  MS <- matrix(0,(n*k-k*ps*s),k*ps)
  MSW <- matrix(0,(n*k-k*ps*s),k*ps)
  for (h in 1:ps){
    MS[(1:(n*k-k*ps*s)),(k*(h-1)+1):(k*((h-1)+1))] <- M1[(k*(ps-h)*s+1):(n*k-(k*h*s)),]
    MSW[(1:(n*k-k*ps*s)),(k*(h-1)+1):(k*((h-1)+1))] <- M2[(k*(ps-h)*s+1):(n*k-(k*h*s)),]
  }

  MA<-matrix(0,(n*k-k*ps*s),k*p)
  MAW<-matrix(0,(n*k-k*ps*s),k*p)
  for (i in 1:p){
    MA[(1:(n*k-ps*s*k)),(k*(i-1)+1):(k*((i-1)+1))]<-M1[((ps*s-i)*k+1):(n*k-(k*i)),]
    MAW[(1:(n*k-ps*s*k)),(k*(i-1)+1):(k*((i-1)+1))]<-M2[((ps*s-i)*k+1):(n*k-(k*i)),]
  }
  ZT<-zt[-(1:(k*ps*s)),]

  space <- colnames(y)

  MA <- data.frame(MA)
  colnamesMA <- c(rep("",p*k))
  for ( i in 1:p ){
    for(j in 1:k ){
      colnamesMA[(i-1)*k+j] <- paste("psi",i,"0",".",space[j],sep ="")
    }
  }
  colnames(MA) <- colnamesMA

  MAW <- data.frame(MAW)
  colnamesMAW <- c(rep("",p*k))
  for ( i in 1:p ){
    for(j in 1:k ){
      colnamesMAW[(i-1)*k+j] <- paste("psi",i,"1",".",space[j],sep ="")
    }
  }
  colnames(MAW) <- colnamesMAW

  MS <- data.frame(MS)
  colnamesMS <- c(rep("",ps*k))
  for ( i in 1:ps ){
    for(j in 1:k ){
      colnamesMS[(i-1)*k+j] <- paste("psi",i,"0",".",s,".",space[j],sep ="")
    }
  }
  colnames(MS) <- colnamesMS

  MSW <- data.frame(MSW)
  colnamesMSW <- c(rep("",ps*k))
  for ( i in 1:ps ){
    for(j in 1:k ){
      colnamesMSW[(i-1)*k+j] <- paste("psi",i,"1",".",s,".",space[j],sep ="")
    }
  }
  colnames(MSW) <- colnamesMSW


  XT<-MA
  WXT<-MAW
  XST <- MS
  WXST <- MSW
  GSTAR<-data.frame(ZT,XT,WXT,XST,WXST)

  fitlm <- lm(ZT~.-1,data=GSTAR)
  residlm <-matrix(fitlm$residuals,ncol = k,byrow=T)

  varcov <- cov(residlm)
  #varcov <- sqrt(varcov)

  mat.I <- diag(nrow(residlm))
  sigma <- kronecker(mat.I,varcov)

  Y <- as.matrix(GSTAR$ZT)
  X <- as.matrix(GSTAR[,-1])
  sigma_inv <- solve(sigma)
  XTX_inv <- solve(t(X) %*% sigma_inv %*% X)
  betahat <- XTX_inv %*% t(X) %*% sigma_inv %*% Y

  sd_betahat <- sqrt(diag(XTX_inv))

  Yhat <- X%*%betahat
  ehat <- Y - Yhat
  param <- length(betahat)/k

  ehat_ <- as.data.frame(matrix(ehat,ncol = k, byrow = T))
  Yhat_ <- as.data.frame(matrix(Yhat,ncol = k, byrow = T))
  Y_    <- as.data.frame(matrix(Y,   ncol = k, byrow = T))

  MSE <- vector(length = ncol(ehat_))
  RMSE <- vector(length = ncol(ehat_))
  AIC <- vector(length = ncol(ehat_))
  R2 <- vector(length = ncol(ehat_))
  temp_SSE <- vector(length = ncol(ehat_))
  temp_SST <- vector(length = ncol(ehat_))
  for(i in 1:ncol(ehat_)){
    temp_et <- ehat_[,i]
    temp_SSE[i] <- sum(temp_et^2)
    temp_SST[i] <- sum(  (Y_[,i] - mean(Y_[,i]))^2 )
    MSE[i] <- temp_SSE[i]/length(temp_et)
    RMSE[i] <- sqrt(MSE[i])
    R2[i] <- 1-(temp_SSE[i]/temp_SST[i])
    obs <- length(temp_et)
    AIC[i] <- exp(param/obs)*temp_SSE[i]/obs
  }
  SSE_total <- sum(temp_SSE)
  MSE_total <- SSE_total/(length(ehat))
  RMSE_total <- sqrt(MSE_total)
  SST_total <- sum(temp_SST)
  R2_total <- 1-(SSE_total/SST_total)


  AIC_total <- (exp(param/obs))*(SSE_total/length(ehat))

  coef <- matrix(c(betahat,sd_betahat),ncol = 2)
  row.names(coef) <- row.names(betahat)
  colnames(coef) <- c("Coefficients","sd")


  Performance_space<-data.frame(MSE=MSE, RMSE=RMSE,AIC=AIC, Rsquare=R2)
  Performance_total <- data.frame(MSE=MSE_total,RMSE = RMSE_total,AIC=AIC_total,Rsquare=R2_total)

  Performance <- rbind(Performance_space,Performance_total)
  rownames(Performance) <- c(space,"total")
  Zt<-as.data.frame(matrix(Y,k,byrow=T))
  colnames(Zt) <- space
  Zhat<-Yhat_
  colnames(Zhat) <- space
  Residual<-ehat_
  colnames(Residual) <- space


  hasil <- list(Coefficients = coef, Fitted.values= Zhat, Residual = Residual, Performance = Performance, p = p , ps = ps , s = s, data=y, weight = w )
  return(hasil)
}
