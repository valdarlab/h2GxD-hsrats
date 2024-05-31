calc.h2 <- function(mc, type=NULL){
  # mc = data.frame of mcmc variance samples: g, gxt, cohort, noise
  if (4!=ncol(mc)){
    stop("Bad input format\n")
  }
  if (is.null(type)){
    stop("Specify type")
  }
  if ("gcor"==type){
    x <- (mc$g - mc$gxt) / (mc$g + mc$gxt)
    attr(x, "range") <- c(-1,1)
    return (x)
  }
  total <- mc$g + mc$gxt + mc$cohort + mc$noise
  x <- switch(type,
              g      = mc$g/total,
              gxt    = mc$gxt/total,
              cohort = mc$cohort/total,
              noise  = mc$noise/total,
              allg   = (mc$g + mc$gxt) / total,
              gcoh   = (mc$g + mc$cohort) / total) 
  attr(x, "range") <- c(0,1)
  if (is.null(x)){ stop("Unknown type:", type) }
  x
}
