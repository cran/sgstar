#' @title Timeseries Plot for Model
#'
#' @description Plotting line chart dataset and fit.values of the Seasonal GSTAR model.
#' @param formula an object from the output from sgstar() function.
#'
#' @return  returns output a list that shown line chart for each location.
#'
#' @export
#' @importFrom  tidyr gather
#' @import dplyr
#' @import ggplot2
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
#' plot1 <- plot_sgstar(fit)
#'
#'
#'

plot_sgstar <- function(formula){

  z_pred <- matrix(NA, nrow =formula$s, ncol=ncol(formula$data))
  colnames(z_pred) <- colnames(formula$data)

  time <- 1:(nrow(formula$data)+formula$s)
  dat <-  rbind(formula$data,z_pred)


  z <- matrix(NA, nrow = (nrow(formula$data)-nrow(formula$Fitted.values)), ncol=ncol(formula$data))
  colnames(z) <- colnames(formula$data)



  fit <- rbind(z,formula$Fitted.values,z_pred)
  colnames(fit) <- paste(colnames(fit),"fit",sep = ".")

  temp_pred <- predict_sgstar(formula,formula$s)
  z1 <- matrix(NA, nrow = nrow(formula$data)-1, ncol=ncol(formula$data))
  pred0<-formula$Fitted.values[nrow(formula$Fitted.values),]
  colnames(z1) <- colnames(formula$data)
  pred <- rbind(z1,pred0,temp_pred)
  colnames(pred) <- paste(colnames(pred),"predict",sep = ".")
  plott <- list()
  temp <- list()
  categories <- c()
  value <- c()

  data_full <- data.frame(time,dat,fit,pred)
  for (q in 1:ncol(formula$data)){
    temp[[q]] <- data_full %>%
      select(c(1,1+q,ncol(formula$data)+q+1,2*ncol(formula$data)+q+1))%>%
      gather(categories,value,c(2,3,4))

    plott[[q]]<- ggplot(temp[[q]], aes(x = time, y = value, group=categories))+
      geom_line(aes(color=categories),na.rm = TRUE) + geom_point(aes(color=categories),na.rm = TRUE)

  }
  return(plot <- plott)
}
