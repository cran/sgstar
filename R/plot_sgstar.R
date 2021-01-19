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
plot_sgstar <- function(formula){

  time <- 1:nrow(formula$data)
  dataset <- cbind(time,formula$data)


  z <- matrix(NA, nrow = (nrow(formula$data)-nrow(formula$Fitted.values)), ncol=ncol(formula$data))
  colnames(z) <- colnames(formula$data)
  fit <- rbind(z,formula$Fitted.values)

  plott <- list()
  temp <- list()
  categories <- c()
  value <- c()
  for (q in 1:ncol(formula$data)){
    temp[[q]] <- data.frame(dataset,fit) %>%
      select(c(1,1+q,ncol(formula$data)+q+1))%>%
      gather(categories,value)

    plott[[q]]<- ggplot(temp[[q]], aes(x = time, y = value, group=categories))+
      geom_line(aes(color=categories)) + geom_point(aes(color=categories))

  }
  return(plot <- plott)
}
