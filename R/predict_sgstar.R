#' @title Predict for Seasonal GSTAR model.
#'
#' @description Predicting value based on Sgstar object
#'
#' @param formula an object from the output from sgstar() function.
#' @param n_time number of steps ahead for which prediction is required.
#'
#' @return returns output a dataframe that shown predict value based on model, with rows as time and column that shown for each location.
#'
#' @references Setiawan, Suhartono, and Prastuti M.(2016).S GSTAR-SUR for Seasonal Spatio Temporal Data Forecasting. Malaysian Journal Of Mathematical Sciences.10.<Corpus ID :189955959>.
#'
#' @export
#'
#' @examples
#' library(sgstar)
#' data("coords")
#' data("simulatedata")
#'
#' #create weight matrix using distance inverse matrix
#' z<-dist(coords,method = "euclidean")
#' z <- as.matrix(z)
#'
#' matriksd <- 1/z
#' matriksd[is.infinite(matriksd)] <- 0
#'
#' matriksd_w <- matriksd / rowSums(as.data.frame(matriksd))
#'
#'
#' fit <- sgstar(data = simulatedata, w = matriksd_w, p = 2,ps = 1, s =4)
#'
#' #predicting for 12 time ahead
#' predict.fit <-predict_sgstar(fit,12)


predict_sgstar<- function(formula, n_time){
  c<- n_time
  n <- nrow(formula$data)
  k <- ncol(formula$data)
  ps <- formula$ps
  p <- formula$p
  s <- formula$s
  w <- formula$w
  zt<-(stack(as.data.frame(t(formula$data)))[,1])

  MDpred<-function(zt){
    M<-matrix(0,(n*k)+c*k,k)
    z=0
    for( i in 1:(n+c)){
      for(j in 1:k){
        z<-z+1
        M[z,j]<-zt[z]
      }
    }
    M   }

  zt[(n*k+1):(n*k+c*k)] <- 0
  MSpred <- matrix(0,c*k,k*ps)
  MSWpred <- matrix(0,c*k,k*ps)
  MApred <- matrix(0,c*k,k*p)
  MAWpred <- matrix(0,c*k,k*p)

  #colnames
  space <- colnames(formula$data)

  MApred <- data.frame(MApred)
  colnamesMApred <- c(rep("",p*k))
  for ( i in 1:p ){
    for(j in 1:k ){
      colnamesMApred[(i-1)*k+j] <- paste("psi",i,"0","(",space[j],")",sep ="")
    }
  }
  colnames(MApred) <- colnamesMApred

  MAWpred <- data.frame(MAWpred)
  colnamesMAWpred <- c(rep("",p*k))
  for ( i in 1:p ){
    for(j in 1:k ){
      colnamesMAWpred[(i-1)*k+j] <- paste("psi",i,"1","(",space[j],")",sep ="")
    }
  }
  colnames(MAWpred) <- colnamesMAWpred

  MSpred <- data.frame(MSpred)
  colnamesMSpred <- c(rep("",ps*k))
  for ( i in 1:ps ){
    for(j in 1:k ){
      colnamesMSpred[(i-1)*k+j] <- paste("psi",i,"0","(",s,")","(",space[j],")",sep ="")
    }
  }
  colnames(MSpred) <- colnamesMSpred

  MSWpred <- data.frame(MSWpred)
  colnamesMSWpred <- c(rep("",ps*k))
  for ( i in 1:ps ){
    for(j in 1:k ){
      colnamesMSWpred[(i-1)*k+j] <- paste("psi",i,"1","(",s,")","(",space[j],")",sep ="")
    }
  }
  colnames(MSWpred) <- colnamesMSWpred

  for(j in 1:c){

    M1pred <- MDpred(zt)
    W1pred<-kronecker(diag(n+c), w)
    ztwpred<-W1pred%*%zt
    M2pred<-MDpred(ztwpred)

    for (h in 1:ps){

      MSpred[(k*(j-1)+1):(j*k),(k*(h-1)+1):(k*h)] <- M1pred[(n*k-(k*h*s)+(j-1)*k+1):(n*k-(k*h*s)+j*k),]
      MSWpred[(k*(j-1)+1):(j*k),(k*(h-1)+1):(k*h)] <- M2pred[(n*k-(k*h*s)+(j-1)*k+1):(n*k-(k*h*s)+j*k),]
    }

    for (i in 1:p){
      MApred[(k*(j-1)+1):(j*k),(k*(i-1)+1):(k*i)]<-
        M1pred[(n*k-(k*i)+(j-1)*k+1):(n*k-(k*i)+(j*k)),]
      MAWpred[(k*(j-1)+1):(j*k),(k*(i-1)+1):(k*i)]<-
        M2pred[(n*k-(k*i)+(j-1)*k+1):(n*k-(k*i)+(j*k)),]
    }

    Xpred <- data.frame(MApred,MAWpred,MSpred,MSWpred)
    Xpred <- as.matrix(Xpred)
    #sum(is.na(Xpred))
    coef <- as.matrix(formula$Coefficient)
    ypred <- Xpred%*%coef
    zt[(n*k+1):(n*k+c*k)] <- as.vector(ypred)

  }

  ypred <- as.data.frame(matrix(ypred,ncol = k,byrow = T))
  colnames(ypred) <- space
  return(ypred)
}
