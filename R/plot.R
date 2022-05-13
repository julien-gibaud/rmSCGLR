if(getRversion()>="2.15.1") {
  # remove warnings due to ggplot2 syntax strangeness
  utils::globalVariables(c("x","y","label","angle","hjust"))
}

#' @title rmSCGLR generic plot
#' @export
#' @importFrom stats cor
#' @importFrom grid circleGrob gpar
#' @import ggplot2
#' @param x an object from rmSCGLR
#' @param thresold correlations with the two components of the plane lower than this threshold will be ignored
#' @param group a integer indicating the group
#' @param plan a size-2 vector indicating which components are plotted
#' @param lin.pred a logical. Should we draw linear predictor
#'
#' @return an object of class \code{\link{ggplot2}}
#'
#' @examples \dontrun{
#'
#'
#' # load sample data
#' data <- genus
#' data <- as.data.frame(apply(data, 2, as.numeric ))
#'
#' # get variable names from dataset
#' n <- names(data)
#' ny <- n[grep("^gen",n)]
#' nx <- n[-grep("^gen",n)]
#' na <- c("geology")
#' nx <- nx[!nx%in%c("geology", "surface")]
#'
#' # build multivariate formula
#' form <- multivariateFormula(Y = ny, X = nx, A = na)
#'
#' # define family
#' fam <- rep("poisson", length(ny))
#'
#' # run function
#' H <- c(2,2)
#' met <- methodSR.RMSCGLR(l=4, s=0.1, t=0.5)
#' res <- ResponseMixtureSCGLR(formula=form, data=data, H=H,
#'                            family=fam, method=met, offset = data$surface)
#'
#' # plot the results
#' plot_RMSCGLR(x=res, thresold=0.5, group=1, plan=c(1,2))
#' }
#'
plot_RMSCGLR <- function(x, thresold=0, group=1, plan=c(1,2), lin.pred=FALSE){
  if (class(x) != "rmSCGLR")
    stop("This plot function need an rmSCGLR result")

  labels.offset <- 0.01
  names <- colnames(x$comp)
  names.group <- paste("^", "G", group, sep = "")
  names <- names[grep(names.group, names)]
  inertia <- cor(x$X, x$comp[,names])
  inertia <- inertia^2
  inertia <- colMeans(inertia)

  p <- qplot()+
    coord_fixed()+ theme(axis.title = element_text(size = 30),
                         plot.title = element_text(size = 30)) +
    labs(title = paste("Component plane (", plan[1], ",", plan[2], ") for the group ", group, sep = "")) +
    # thicker x unit arrow
    xlab(paste("SC", plan[1], " (", round(100*inertia[plan[1]],2), "%", ")", sep = "")) +
    geom_hline(yintercept=0)+
    geom_segment(aes(x=-1.1,xend=1.1,y=0,yend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))+
    # thicker y unit arrow
    ylab(paste("SC", plan[2], " (", round(100*inertia[plan[2]],2), "%", ")", sep = "")) + #theme(axis.title.y = element_text(size = 20)) +
    geom_vline(xintercept=0)+
    geom_segment(aes(y=-1.1,yend=1.1,x=0,xend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))
  p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(fill=NA)),-1,1,-1,1)
  p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(lty=2, fill=NA)),
                             -thresold,thresold,-thresold,thresold)

  co1 <- as.data.frame(cor(x$X, x$comp[,names][,c(plan[1], plan[2])]))
  names(co1) <- c("x", "y")
  co1$norm <- sqrt(co1$x^2+co1$y^2)
  co1$label <- names(as.data.frame(x$X))
  co1$arrows.color <- 'black'
  co1$labels.color <- 'black'
  co1$labels.size <- 6
  co1$angle <- atan2(co1$y,co1$x)*180/pi
  co1$hjust <- ifelse(abs(co1$angle)>90,1,0)
  co1$angle <- ifelse(abs(co1$angle)>90,co1$angle+180,co1$angle)

  co1 <- co1[co1$norm>thresold,]

  p <- p + geom_segment(
    aes(x=0,y=0,xend=x,yend=y),
    data=co1,
    color=co1$arrows.color,
    arrow=arrow(length=unit(0.02,"npc"))
  )

  p <- p + geom_text(
    aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
    data=co1,
    color=co1$labels.color,
    size=co1$labels.size
  )
  # browser()
  if(lin.pred==TRUE){
    co2 <- as.data.frame(cor(x$eta[[group]][,x$cluster==group], x$comp[,names][,c(plan[1], plan[2])]))
    # browser()
    names(co2) <- c("x", "y")
    co2$norm <- sqrt(co2$x^2+co2$y^2)
    co2$label <- names(as.data.frame(x$Y))[x$cluster==group]
    co2$arrows.color <- 'red'
    co2$labels.color <- 'red'
    co2$labels.size <- 6
    co2$angle <- atan2(co2$y,co2$x)*180/pi
    co2$hjust <- ifelse(abs(co2$angle)>90,1,0)
    co2$angle <- ifelse(abs(co2$angle)>90,co2$angle+180,co2$angle)

    co2 <- co2[co2$norm>thresold,]

    p <- p + geom_segment(
      aes(x=0,y=0,xend=x,yend=y),
      data=co2,
      color=co2$arrows.color,
      arrow=arrow(length=unit(0.02,"npc"))
    )

    p <- p + geom_text(
      aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
      data=co2,
      color=co2$labels.color,
      size=co2$labels.size
    )

  }

  return(list(p=p))
}
