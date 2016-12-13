#~ akima spline for log birth rate

require(akima)
require(Rcpp)
require(BH)
#~ sourceCpp('semiparDemoProcessModel0.cpp')

.transform.splinePoints <- function(sp, r)
{
	nsp <- sp / tail(sp,1)
	dsp <- diff(nsp) + exp(r)
	c(sp[1], sp[1] + (tail(sp,1) - sp[1]) * cumsum(dsp) / sum(dsp))
}

.spline.points <- function(theta, t0, t1)
{
	pnames <- names(theta)
	sps_names <- pnames[ grepl('splinePoint', pnames) ]
	np <- sum(grepl('akima', names(theta)) )
	if (np==1){
		return(t0)
	} else if(np==2){
		return(c(t0, t1))
	}
	if ((length(sps_names) + 2)!=np) stop('number spline point parameters must equal number akima parameters - 2') #NOTE number parameters here is np - 2
	
	# first element fixed at 1/2, sets scale; other values restricted (0,1)
	x <-  cumsum( c(1/2 , exp(theta[sps_names] ) / ( 1 + exp(theta[sps_names])) ) ) 
	c( 0, x / tail(x, 1) )*(t1-t0)  + t0
}

gen.demo.model3 <- function( t0, bdt, npoints = 5, model_gamma = NA)
{
	#print(splinePoints)
	splineNames <- paste( 'akima', 1:npoints, sep='')
	function( theta, x0, t0, t1, res = 1000, tfin = NULL)
	{
		# reqd parameters akima1...akima_npoints, gamma
		tfin <- ifelse(is.null(tfin), t1, tfin )
		times <- seq(t0, tfin, length.out = res)
		
		sps <- .spline.points( theta, t0, t1 )
		if (length(sps)==1){
			betas <- rep( exp(theta[splineNames]), length(times))
		} else {
			betas <- exp(aspline( sps, theta[splineNames], xout = times  )$y)
			betas[ times > t1 ] <- exp( tail( theta[splineNames] ,1 ) )
		}
		if (!is.na(model_gamma)) {
			y <- as.vector( solve_semiPar0( t0, max(t1,tfin), res, x0[1], betas, model_gamma ))
		} else if (!is.na(theta['gamma'])) {
			y <- as.vector( solve_semiPar0( t0, max(t1,tfin), res, x0[1], betas, theta['gamma'] ))
		} else{
			stop('model_gamma or theta[gamma] must be defined')
		}
		list( times = rev(times)
		, f = rev( betas * y  )
		, y = rev( y) )
	}
}

gen.demo.model2 <- function( t0, bdt, npoints = 5, model_gamma = NA)
{
	intimes<- bdt$maxSampleTime - bdt$heights[(bdt$n+1):(bdt$n+bdt$Nnode)]
	if (npoints >= 3){
		splinePoints <- c( t0,  quantile( intimes , probs = 1:(npoints-2)/(npoints-1) ), bdt$maxSampleTime )
	} else if ( npoints == 2){
		splinePoints <- c( t0, bdt$maxSampleTime )
	} else if (npoints == 1){
		splinePoints <- t0
	} else{
		stop( 'npoints must be integer >= 1')
	}
	#print(splinePoints)
	splineNames <- paste( 'akima', 1:npoints, sep='')
	function( theta, x0, t0, t1, res = 1000)
	{
		# reqd parameters akima1...akima_npoints, gamma
		times <- seq(t0, t1, length.out = res)
		
		sps <- splinePoints
		if ('transformSplinePoints' %in% names(theta)){
			sps <- .transform.splinePoints(splinePoints, theta['transformSplinePoints'])
		} 
		
		if (length(splinePoints)==1){
			betas <- rep( exp(theta[splineNames]), length(times))
		} else {
			betas <- exp(aspline( sps, theta[splineNames], xout = times  )$y)
		}
		if (!is.na(model_gamma)) {
			y <- as.vector( solve_semiPar0( t0, t1, res, x0[1], betas, model_gamma ))
		} else if (!is.na(theta['gamma'])) {
			y <- as.vector( solve_semiPar0( t0, t1, res, x0[1], betas, theta['gamma'] ))
		} else{
			stop('model_gamma or theta[gamma] must be defined')
		}
		list( times = rev(times)
		, f = rev( betas * y  )
		, y = rev( y) )
	}
}



