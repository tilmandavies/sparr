bivariate.density <- function(pp,h0,hp=NULL,adapt=FALSE,resolution=128,gamma.scale="geometric",edge=c("uniform","diggle","none"),intensity=FALSE,trim=5,xy=NULL,pilot.density=NULL,leaveoneout=FALSE,parallelise=NULL,davies.baddeley=NULL,verbose=TRUE){
	if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")

	W <- Window(pp)
	if(!is.null(xy)){
	  xy <- checkxy(xy)
	  dimyx <- NULL
	  resolution <- length(xy$x)
	} else {
	  resolution <- checkit(resolution,"'resolution'")
	  dimyx <- rep(resolution,2)
	}

	edge <- checkedge(edge)
	h0 <- checkit(h0,"'h0'")
	if(!is.null(trim)&&!is.na(trim)) trim <- checkit(trim,"'trim'")
	
	n <- npoints(pp)

	if(adapt){
		if(is.null(hp)) hp <- h0
		else hp <- checkit(hp,"'hp'")
		
		if(leaveoneout) return(bivden.LOO(pp,h0,hp,gamma.scale,trim,resolution,parallelise))
    
		pd <- pilot.density
		pilot.data <- pp
		if(!is.null(pd)){
			if(is.im(pd)){
				if(is.null(xy)){
			    		if(!all(dim(pd)==resolution)) stop("'pilot.density' image resolution must strictly have 'resolution' x 'resolution' pixels")
				} else {
			    		if((!all(pd$xcol==xy$x))||(!all(pd$yrow==xy$y))) stop("'pilot.density' xcol and yrow must strictly match coords in 'xy'")
				}
		  	pilot.density[pd<=0] <- min(pd[pd>0])
			} else if(is.ppp(pd)){
				pilot.data <- pd
				if(!identical_windows(Window(pp),Window(pilot.data))) stop("'pilot.density' window must be identical to 'pp' window")
				pilot.density <- density(pilot.data,sigma=hp,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),dimyx=dimyx,xy=xy,positive=TRUE)
			} else {
				stop("'pilot.density' must be an object of class \"im\" or \"ppp\"")
			}
		} else {
			pilot.density <- density(pp,sigma=hp,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),dimyx=dimyx,xy=xy,positive=TRUE)
		}
		
		pilot.density.spec <- safelookup(pilot.density,pp,warn=FALSE)
		pi.int <- integral(pilot.density)
		pilot.density <- pilot.density/pi.int
		pilot.density.spec <- pilot.density.spec/pi.int
		pspec <- pilot.density.spec^(-0.5)
		gamma <- processgamma(gamma.scale,safelookup(pilot.density,pilot.data,warn=FALSE)) #'was: processgamma(gamma.scale,pilot.density.spec)
		gs <- gspd <- exp(mean(log(pspec))) 
		if(!is.null(pd)) gspd <- exp(mean(log(safelookup(pilot.density,pilot.data,warn=FALSE)^(-0.5))))
		
		# PREVIOUS TRIMMING REGIMEN #
		# h.spec <- h0*pilot.density.spec^(-0.5)/gamma
		# h.hypo <- h0*pilot.density^(-0.5)/gamma
		# if(is.null(trim)) beta.h <- 5*median(h.spec,na.rm=TRUE)
		# else if(is.na(trim)) beta.h <- max(h.hypo,na.rm=TRUE)
		# else beta.h <- trim
		# h.spec[h.spec>beta.h] <- beta.h
		# h.hypo[h.hypo>beta.h] <- beta.h
		# h.hypo.mat <- as.matrix(h.hypo)
		
		# NEW TRIMMING REGIMEN #
		# h.spec <- h0*pmin(pspec/gamma,trim)  ### Generalised below for numeric gamma argument vals. Trimming is universally determined by the geometric mean 'gs', regardless of 'gamma.scale' ###
		# h.hypo <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5)))/gamma,trim),resolution,resolution),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		
		h.spec <- h0*pmin(pspec,trim*gspd)/gamma
		h.hypo <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5))),trim*gspd),resolution,resolution)/gamma,xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		h.hypo.mat <- as.matrix(h.hypo)

		if(!is.null(davies.baddeley)){
			db <- checkdb(davies.baddeley)
			db.result <- adens(x=pp,bwim=h.hypo,bwpts=h.spec,resolution=resolution,intensity=intensity,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),weights=NULL,hstep=db[1],qstep=db[2],qres=ifelse(is.na(db[3]),resolution,db[3]),verbose)
			result <- list(z=db.result$result,h0=h0,hp=hp,h=h.spec,him=h.hypo,q=db.result$edg,gamma=gamma,geometric=gs,pp=pp)
			class(result) <- "bivden"		
			return(result)
		}

		evalxy <- as.matrix(expand.grid(pilot.density$xcol,pilot.density$yrow))
		notin <- !inside.owin(x=evalxy[,1],y=evalxy[,2],w=W)
		surf <- rep(NA,nrow(evalxy))
		ef <- NULL
		
		if(edge=="uniform"){
			qhz <- rep(NA,resolution^2)
			if(verbose) pb <- txtProgressBar(0,nrow(evalxy))
			for(i in 1:nrow(evalxy)){
				ht <- h.hypo.mat[which(pilot.density$yrow==evalxy[i,2]),which(pilot.density$xcol==evalxy[i,1])]
				if(is.na(ht)) next
		      
		    	gxy <- ht^(-2)*(exp(-0.5*rowSums((cbind(evalxy[,1]-evalxy[i,1],evalxy[,2]-evalxy[i,2])/ht)^2))/(2*pi))
		    	gxy[notin] <- NA
		    	qhz[i] <- integral(im(matrix(gxy,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow))
		      
		    	uxy <- cbind(evalxy[i,1]-pp$x,evalxy[i,2]-pp$y)/h.spec
		    	ivals <- h.spec^(-2)*(exp(-0.5*rowSums(uxy^2))/(2*pi))
		    
		    	if(!intensity) surf[i] <- mean(ivals)/qhz[i]
		    	else surf[i] <- sum(ivals)/qhz[i]
		    	if(verbose) setTxtProgressBar(pb,i)
			}
			if(verbose) close(pb)
			ef <- im(matrix(qhz,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		}
		
		if(edge=="diggle"){
		  qx <- rep(1,n)
		  if(verbose) pb <- txtProgressBar(0,n+nrow(evalxy))
		  for(i in 1:n){
		    pxy <- h.spec[i]^(-2)*(exp(-0.5*rowSums((cbind(evalxy[,1]-pp$x[i],evalxy[,2]-pp$y[i])/h.spec[i])^2))/(2*pi))
		    pxy[notin] <- NA
		    qx[i] <- integral(im(matrix(pxy,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow))
		    if(verbose) setTxtProgressBar(pb,i)
		  }
		  
		  for(i in 1:nrow(evalxy)){
		    uxy <- cbind(evalxy[i,1]-pp$x,evalxy[i,2]-pp$y)/h.spec
		    ivals <- h.spec^(-2)*(exp(-0.5*rowSums(uxy^2))/(2*pi))
		      
		    if(!intensity) surf[i] <- mean(ivals/qx)
		    else surf[i] <- sum(ivals/qx)
		    if(verbose) setTxtProgressBar(pb,n+i)
		  }
		  if(verbose) close(pb)
		  ef <- qx
		}
		  
		if(edge=="none"){
		  if(verbose) pb <- txtProgressBar(0,nrow(evalxy))
		  for(i in 1:nrow(evalxy)){
		    uxy <- cbind(evalxy[i,1]-pp$x,evalxy[i,2]-pp$y)/h.spec
		    ivals <- h.spec^(-2)*(exp(-0.5*rowSums(uxy^2))/(2*pi))
		    if(!intensity) surf[i] <- mean(ivals)
		    else surf[i] <- sum(ivals)
		    if(verbose) setTxtProgressBar(pb,i)
		  }
		  if(verbose) close(pb)
		}
		
		surf[notin] <- NA
		surf <- im(matrix(surf,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		
	} else {
    h.spec <- rep(h0,n)
	  h.hypo <- ef <- NULL
	  gs <- gamma <- NA
    
	  dens <- density.ppp(pp,sigma=h0,dimyx=dimyx,xy=xy,edge=(edge=="diggle"||edge=="uniform"),diggle=(edge=="diggle"),spill=1)
		surf <- dens$raw[W,drop=FALSE]
		ef <- dens$edg[W,drop=FALSE]
		
		if(edge=="diggle"){
		  ef <- safelookup(ef,pp,warn=FALSE)
		} else if(edge=="uniform"){
		  surf <- surf/ef
		  surf[surf<0] <- 0
		}
		ef[ef>1] <- 1
		if(!intensity) surf <- surf/integral(surf)
	}
	
	result <- list(z=surf,h0=h0,hp=hp,h=h.spec,him=h.hypo,q=ef,gamma=gamma,geometric=gs,pp=pp,fromms=FALSE)
	class(result) <- "bivden"		
	
	return(result)
}