calc.BF <- function(prior.lsp, post.lsp, xvalue=0, plotit=FALSE,
                    pheno.string="phenotype", 
                    param.string="", ylim=NULL, xlim=0:1, digits=3){
  # estimate density at x
  prior.at.x <- dlogspline(xvalue, prior.lsp)
  post.at.x  <- dlogspline(xvalue, post.lsp)
  # calc Bayes factor
  BF <- prior.at.x / post.at.x
  # plotting
  if (plotit){
    n <- 100
    xx <- (0:(n-1))/(n-1) * diff(xlim) - xlim[1]
    yy.prior <- dlogspline(xx, prior.lsp)
    yy.post  <- dlogspline(xx, post.lsp)
    if (is.null(ylim)){
      ylim <- c(0, max(c(yy.prior, yy.post)))
    }
    plot(post.lsp, xlim=xlim, ylim=ylim, lty=1, las=1, xlab=param.string,
         ylab="density",
         main=paste0(pheno.string, "\nBF for ", param.string, " is ",
                     round(BF, digits=digits)))
    plot(prior.lsp, xlim=xlim, lty=2, add=TRUE)
  }
  data.frame(
    phenotype=pheno.string, 
    parameter=param.string,
    xvalue=xvalue, 
    prior.at.x=prior.at.x,
    post.at.x=post.at.x,
    BF=BF,
    logBF=log10(BF))
}
