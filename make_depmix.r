# NOT RUN {
require(gamlss)
require(gamlss.dist)
require(depmixS4)

library("datasets")
library("gamlss")
library("gamlss.dist")
library("depmixS4")
data(speed)
rt <- speed$rt

# define a response class which only contains the standard slots, no additional slots
setClass("exgaus", contains="response")

# define a generic for the method defining the response class

setGeneric("exgaus", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("exgaus"))

# define the method that creates the response class

setMethod("exgaus", 
    signature(y="ANY"), 
    function(y, pstart=NULL,fixed=NULL, ...) {
        y <- matrix(y,length(y))
		x <- matrix(1)
		parameters <- list()
		npar <- 3
		if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
		if(!is.null(pstart)) {
		if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
		  parameters$mu <- pstart[1]
		  parameters$sigma <- log(pstart[2])
		  parameters$nu <- log(pstart[3])
        }
        mod <- new("exgaus",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
        mod
	}
)

setMethod("show","exgaus",
    function(object) {
        cat("Model of type exgaus (see ?gamlss for details) \n")
        cat("Parameters: \n")
        cat("mu: ", object@parameters$mu, "\n")
        cat("sigma: ", object@parameters$sigma, "\n")
        cat("nu: ", object@parameters$nu, "\n")
    }
)

setMethod("dens","exgaus",
 function(object,log=FALSE) {
   lambda <- object@parameters$lambda +
            object@parameters$gamma1 * cos(2 * pi * object@x / period) +
            object@parameters$gamma2 * sin(2 * pi * object@x / period)
  
  # Utiliser dpois pour calculer la densité
  dpois(object@y, lambda = lambda, log = log) + object@parameters$alpha * object@a
    
  }
)

setMethod("getpars","response",
    function(object,which="pars",...) {
        switch(which,
            "pars" = {
                parameters <- numeric()
                parameters <- unlist(object@parameters)
                pars <- parameters
            },
            "fixed" = {
                pars <- object@fixed
            }
        )
        return(pars)
    }
)

setMethod("setpars","exgaus",
    function(object, values, which="pars", ...) {
        npar <- npar(object)
        if(length(values)!=npar) stop("length of 'values' must be",npar)
        # determine whether parameters or fixed constraints are being set
		nms <- names(object@parameters)
		switch(which,
		  "pars"= {
		      object@parameters$mu <- values[1]
		      object@parameters$sigma <- values[2]
		      object@parameters$nu <- values[3]
		      },
		  "fixed" = {
		      object@fixed <- as.logical(values)
		  }
		  )
        names(object@parameters) <- nms
        return(object)
    }
)

setMethod("fit","exgaus",
    function(object,w) {
        if(missing(w)) w <- NULL
        y <- object@y
        fit <- gamlss(y~1,weights=w,family=exGAUS(),
			control=gamlss.control(n.cyc=100,trace=FALSE),
			mu.start=object@parameters$mu,
			sigma.start=exp(object@parameters$sigma),
			nu.start=exp(object@parameters$nu))
		pars <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients)
		object <- setpars(object,pars)
		object
	}
)

setMethod("predict","exgaus", 
    function(object) {
        ret <- object@parameters$mu
        return(ret)
    }
)

rModels <- list(
  list(
	  exgaus(rt,pstart=c(5,.1,.1)),
		GLMresponse(formula=corr~1, data=speed, 
		family=multinomial("identity"), pstart=c(0.5,0.5))
	),
	list(
		exgaus(rt,pstart=c(6,.1,.1)),
		GLMresponse(formula=corr~1, data=speed, 
		family=multinomial("identity"), pstart=c(.1,.9))
	)
)

trstart=c(0.9,0.1,0.1,0.9)
transition <- list()
transition[[1]] <- transInit(~Pacc,nstates=2,data=speed,pstart=c(trstart[1:2],0,0))
transition[[2]] <- transInit(~Pacc,nstates=2,data=speed,pstart=c(trstart[3:4],0,0))

instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,family=multinomial("identity"), data=data.frame(rep(1,3)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=c(168,134,137), 
homogeneous=FALSE)

fm3 <- fit(mod,emc=em.control(rand=FALSE))
summary(fm3)
# }
# NOT RUN {
# }