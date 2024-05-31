

load.posterior.varsamples <- function(file, burnin=1000, thin=1){
  x <- load.to.object(file)
  nall <- dim(x$alpha)[3] # get number of timesteps
  s <- seq(from=burnin+1, to=nall, by=thin)
  # y ~ litter.size + study  + z       + (1|cohort)                + (1|kinship)            + (z|kinship)
  #     alpha1        alpha2   alpha3    u ~ N(, sigma_sq_lambda)    N(, sigma_sq_beta[1])    N(, sigma_sq_beta[2])
  data.frame(
    g      = c(x$sigma_sq_beta[1,,s]),
    gxt    = c(x$sigma_sq_beta[2,,s]),
    cohort = c(x$sigma_sq_lambda[,s]),
    noise  = c(x$sigma_sq[,s]))
}

## calc.h2 <- function(mc, type=NULL){  # mc = data.frame of mcmc variance samples: g, gxt, cohort, noise
##   if (4!=ncol(mc)){
##     stop("Bad input format\n")
##   }
##   if (is.null(type)){
##     stop("Specify type")
##   }
##   if ("gcor"==type){
##     x <- (mc$g - mc$gxt) / (mc$g + mc$gxt)
##     attr(x, "range") <- c(-1,1)
##     return (x)
##   }
##   total <- mc$g + mc$gxt + mc$cohort + mc$noise
##   x <- switch(type,
##               g      = mc$g/total,
##               gxt    = mc$gxt/total,
##               cohort = mc$cohort/total,
##               noise  = mc$noise/total,
##               allg   = (mc$g + mc$gxt) / total,
##               gcoh   = (mc$g + mc$cohort) / total) 
##   attr(x, "range") <- c(0,1)
##   if (is.null(x)){ stop("Unknown type:", type) }
##   x
## }

## get.funct.range <- function(type){
##   if (type %in% c("g", "gxt", "cohort", "noise", "allg")){
##     return (c(0,1))
##   }
##   if ("gcor"==type){
##     return (c(-1,1))
##   }
##   stop("Unknown type:", type)
## }

plot.post <- function(d, prior.lsp, main="", xlim=c(0,1), ylim=NULL, digits=3,
                      xlab="intraclass correlation coefficient (heritability)"){
  n <- 100
  xx <- (0:(n-1))/(n-1)*diff(xlim) + xlim[1]
  yy <- dlogspline(xx, prior.lsp)
  for (vc in names(d)){
    yy <- c(yy, dlogspline(xx, d[[vc]]$post$lsp))
    bf <- d[[vc]]$post$bf
    if (!is.null(bf)){
      if (1==length(bf)){
        d[[vc]]$string <- paste0(d[[vc]]$string, "  (BF=", round(bf, digits=digits), ")")
      } else if (2==length(bf)){
        d[[vc]]$string <- paste0(d[[vc]]$string, 
                                 "\nBF.0=", round(bf[1], digits=digits), 
                                 "\nBF.1=", round(bf[2], digits=digits))
      }
    }
  }
  if (is.null(ylim)){
    ylim <- c(0, max(yy))
  }
  plot(prior.lsp, xlim=xlim, ylim=ylim, lty=2, las=1, 
       xlab=xlab, ylab="posterior density", main=main)
  for (vc in names(d)){
    plot(d[[vc]]$post$lsp, xlim=xlim, lty=1, add=TRUE, col=d[[vc]]$col)
  }
  legend("topright", bty="n",
          legend=c("prior", sapply(d, function(x){x$string})),
          col=c("black", sapply(d, function(x){x$col})),
          lty=c(2, rep(1, length(d)))
          )
}

plot.logBFs <- function(res, ...){
  oldpar <- par()
  par(mar=c(5, 20, 5, 2))
  plot(res$logBF, nrow(res):1, axes=FALSE, ylab="", type="o", pch=c(1,1,19,19), ...)
  axis(2, at=nrow(res):1, labels=res$phenotype, las=1)
  abline(v=c(-1,0,1), lty=2, col=c("gray", "black", "gray"))
  axis(1)
  axis(3)
}

plot.logBFs.paired <- function(phenos, y1, y2, col=c("red", "blue"), xlab="logBF against zero", ...){
  oldpar <- par()
  par(mar=c(5, 20, 5, 2))
  n <- length(phenos)
  plot(y1, n:1, axes=FALSE, xlim=range(c(y1,y2)), ylab="", type="o", pch=1, col=col[1], xlab=xlab, ...)
  points(y2, n:1, ylab="", type="o", pch=1, col=col[2])
  axis(2, at=n:1, labels=phenos, las=1)
  abline(v=c(-1,0,1), lty=2, col=c("gray", "black", "gray"))
  axis(1)
  axis(3)
  box()
  par(oldpar)
}


plot.logBFs.paired.multi <- function(phenos, data,
                                col=c("red", "blue"),
                                xlim=NULL,
                                gap=0.1, # gap between plots; assumes each plot is 1 unit wide
                                ...){
  # expects dat=list(name=, mat=, xlim=)
  oldpar <- par()
  par(mar=c(5, 20, 5, 2))
  n <- length(phenos)
  K <- length(data)
  xstart <- (1:K - 1)*(1+gap)
  xend   <- xstart + 1
  plot(c(xstart[1], xend[K]), c(1,n), xlim=c(0,xend[K]), axes=FALSE, type="n", xlab="", ylab="")
  for (k in 1:K){
    dat <- data[[k]]
    mtext(dat$name, side=3, at=xstart[k]+0.5, line=2, font=2)
    mtext(dat$xlab, side=1, at=xstart[k]+0.5, line=2, font=1)
    segments(x0=rep(xstart[k], n), y0=1:n, x1=rep(xend[k], n), y1=1:n, col="gray90")
    d.range <- range(dat$mat)
    if (!is.null(xlim)){
      d.range <- xlim
    }
    if (!is.null(dat$xlim)){
      d.range <- dat$xlim
    }
    d2x <- function(d){ (d - d.range[1]) / diff(d.range) + xstart[k] }
    abline(v=d2x(c(-1,0,1)), lty=2, col=c("gray", "black", "gray"))
    x.mat <- d2x(dat$mat) 
    for (j in 1:ncol(x.mat)){
      points(x.mat[,j], n:1, type="o", pch=1, col=col[j])
    }
    d.ticks <- pretty(d.range)
    d.ticks <- d.ticks[d.ticks >= d.range[1] & d.ticks <= d.range[2]]
    axis(1, at=d2x(d.ticks), labels=d.ticks, las=1)
    axis(3, at=d2x(d.ticks), labels=d.ticks, las=1)
  }
  axis(2, at=n:1, labels=phenos, las=1)
  par(oldpar)
}

# 
# data is a K-length list of subplots. Each of the K subplots is a list with fields
#   including:
#      name : the name of the x-axis label
#      mat  : a 3d array with dimensions n x S x J, for n observations measured on
#             S confidence interval statistics (mean, median, lower50, upper50, lower95,
#             upper95), for J conditions.
plot.ci.multi <- function(phenos,
                          data,
                          col=c("red", "blue"),
                          xlim=NULL,
                          gap=0.1, # gap between plots; assumes each plot is 1 unit wide
                          vgap=0.1,
                          ...){
  oldpar <- par()
  par(mar=c(5, 20, 5, 2))
  n <- length(phenos)
  K <- length(data)
  xstart <- (1:K - 1)*(1+gap)
  xend   <- xstart + 1
  plot(c(xstart[1], xend[K]), c(1,n), xlim=c(0,xend[K]), axes=FALSE, type="n", xlab="", ylab="")
  for (k in 1:K){
    dat <- data[[k]]
    mtext(dat$name, side=3, at=xstart[k]+0.5, line=2, font=2)
    mtext(dat$xlab, side=1, at=xstart[k]+0.5, line=2, font=1)
    segments(x0=rep(xstart[k], n), y0=1:n, x1=rep(xend[k], n), y1=1:n, col="gray90")
    d.range <- range(dat$mat)
    if (!is.null(xlim)){
      d.range <- xlim
    }
    if (!is.null(dat$xlim)){
      d.range <- dat$xlim
    }
    d2x <- function(d){ (d - d.range[1]) / diff(d.range) + xstart[k] }
    if (!is.null(dat$dashed)){
      abline(v=d2x(dat$dashed), lty=2, col=c("gray"))
    }
    if (!is.null(dat$full)){
      abline(v=d2x(dat$full), col="black")
    }
    x.mat <- d2x(dat$mat)
    y.offset <- scale(1:dim(dat$mat)[3])*vgap
    for (j in 1:dim(dat$mat)[3]){
      y <- n:1 + y.offset[j]
      ci95 <- x.mat[, 5:6 ,j]
      segments(x0=ci95[, 1], x1=ci95[, 2], y0=y, y1=y, lwd=1, col=col[j])
      ci50 <- x.mat[, 3:4, j]
      segments(x0=ci50[, 1], x1=ci50[, 2], y0=y, y1=y, lwd=3, col=col[j])
      points(x.mat[,1,j], y, type="p", pch="|", col=col[j])
      points(x.mat[,2,j], y, type="p", pch="|", col="white")
    }
    d.ticks <- pretty(d.range)
    d.ticks <- d.ticks[d.ticks >= d.range[1] & d.ticks <= d.range[2]]
    axis(1, at=d2x(d.ticks), labels=d.ticks, las=1)
    axis(3, at=d2x(d.ticks), labels=d.ticks, las=1)
  }
  axis(2, at=n:1, labels=phenos, las=1)
  par(oldpar)
}




## calc.BF <- function(prior.lsp, post.lsp, xvalue=0, plotit=TRUE,
##                     pheno.string="phenotype", 
##                     param.string="", ylim=NULL, xlim=0:1, digits=3){
##   # estimate density at x
##   prior.at.x <- dlogspline(xvalue, prior.lsp)
##   post.at.x  <- dlogspline(xvalue, post.lsp)
##   # calc Bayes factor
##   BF <- prior.at.x / post.at.x
##   # plotting
##   if (plotit){
##     n <- 100
##     xx <- (0:(n-1))/(n-1) * diff(xlim) - xlim[1]
##     yy.prior <- dlogspline(xx, prior.lsp)
##     yy.post  <- dlogspline(xx, post.lsp)
##     if (is.null(ylim)){
##       ylim <- c(0, max(c(yy.prior, yy.post)))
##     }
##     plot(post.lsp, xlim=xlim, ylim=ylim, lty=1, las=1, xlab=param.string, ylab="density",
##         main=paste0(pheno.string, "\nBF for ", param.string, " is ", round(BF, digits=digits)))
##     plot(prior.lsp, xlim=xlim, lty=2, add=TRUE)
##   }
##   data.frame(
##     phenotype=pheno.string, 
##     parameter=param.string,
##     xvalue=xvalue, 
##     prior.at.x=prior.at.x,
##     post.at.x=post.at.x,
##     BF=BF,
##     logBF=log10(BF))
## }
