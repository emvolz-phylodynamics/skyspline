# this version uses gen.demo.model3 , flexible spline points

.lrt.accept.h1 <- function(fit0, fit1, np0=1,np1=2){
	# assume fit1 has one extra parameter, test if fit1 is sign. better
	if (is.null(fit0) | is.null(fit1)){
		return( TRUE )
	}
	if ( 2*((-fit1$value) - (-fit0$value) ) > qchisq(.95,np1-np0) ){
		return(TRUE)
	}
	FALSE
}

.gen.start.splinePoint.parms <- function(np){
	if (np<3) return(NULL)
	setNames( rep(0, np-2)
	 , paste( sep='', 'splinePoint', 1:(np-2))
	)
}

#ML version
fit.skyspline.mle3 = fit.skyspline.ml <- function(bdts
  , death_rate_guess 
  , t0 = 0
  , R0guess = 2
  , cotype=c('com12', 'mscom')
  , Ne_guess = NA
  , y0_guess = NA
  , np_range = NULL
  , priors = list() # named list of functions ln.d(x) 
  , trace =FALSE
  , ... # passed to likelihood
){
	if (class(bdts)[1]=='DatedTree'){
		bdt = bdts
	} else if(class(bdts)[1]=='phylo'){
		n <- length(bdts$tip.label)
		x <- setNames( dist.nodes( bdts)[ (n+1), 1:n ], bdts$tip.label )
		bdt = bdts <- DatedTree( bdts, x )
	}else if (class(bdts)[1]=='list'){
		for (k in 1:length(bdts)){
			bdt <- bdts[[k]]
			if(class(bdt)=='phylo'){
				n <- length(bdt$tip.label)
				x <- setNames( dist.nodes( bdt)[ (n+1), 1:n ], bdt$tip.label)
				bdts[[k]] <- DatedTree( bdt, x )
			}
		}
		bdt <- bdts[[1]]
	} else{
		stop('bdts must be DatedTree or list of DatedTree')
	}
	if (is.null(np_range)){
		np_range <- c(1,Inf)
	}
	done <- FALSE
	np <- np_range[1]
	est_death_rate <- ( 'lngamma' %in% names(priors) | 'death_rate' %in% names(priors) )
	if (!est_death_rate){
		dm <- gen.demo.model3( t0, bdt, npoints = np, model_gamma = death_rate_guess)
	} else{
		dm <- gen.demo.model3( t0, bdt, npoints = np)
	}
	sp_start <- .gen.start.splinePoint.parms(np)
	
	if (!('R0'%in% names(priors))) priors[[length(priors)+1]] <- function(x) dexp( x, R0guess/2, log=T)
	
	if (is.na(y0_guess)) y0_guess <- 1
	if (is.na(Ne_guess)) Ne_guess <- 1 / death_rate_guess
	
	if (cotype=='mscom'){
		a_start <- log( R0guess * death_rate_guess )
		anames <- paste( sep='', 'akima', 1:np)
		theta_start <- c( setNames( rep(a_start , np ), anames) 
		 , sp_start
		 , lnNe = log(Ne_guess)
		 , lny0 = log(1)
		)
	} else if( cotype == 'com12'){
		a_start <- log( R0guess * death_rate_guess )
		anames <- paste( sep='', 'akima', 1:np)
		theta_start <- c( setNames( rep(a_start , np ), anames) 
		 , sp_start
		 , lny0 = log(1)
		)
	}
	if (est_death_rate){
		theta_start <- c( theta_start, lngamma = unname(log(death_rate_guess)) )
	}
	
	ln.dprior <- function(x){
		lnd <- 0
		if (length(priors)>0){
			for( xn in names(x)){
				if (xn %in% names(priors)){
					lnd <- lnd + priors[[xn]]( x[xn] )
				}
			}
			if ('R0' %in% names(priors)){
				avals <- exp( x[grepl('akima', names(x))] )
				g <- ifelse( 'lngamma' %in% names(x), exp(x['lngamma']), death_rate_guess)
				R0s <- avals / g
				lnd <- lnd + sum(sapply(R0s, function(R0) priors$R0(R0) )) 
			}
			if ('death_rate' %in% names(priors)){
				lnd <- lnd + priors$death_rate(exp(x['lngamma']) )
			}
		}
		lnd
	}
	
	.of.mscom <- function( x, ...)
	{
		Ne <- exp(x['lnNe'])
		y0 <- exp(x['lny0'])
		if (est_death_rate) x['gamma'] <- unname(exp(x['lngamma']))
		rv <- -colik.mscom(bdts
			  , x
			  , dm
			  , y0
			  , t0 = t0
			  , Ne = Ne
			  , maxk=10
			  , ...
			)
		if (trace) print(c(exp=exp(x), -rv))
		rv - ln.dprior(x)
	}
	
	.of.com12 <- function( x, ...)
	{
		y0 <- exp(x['lny0'])
		if (est_death_rate) x['gamma'] <- unname(exp(x['lngamma']))
		rv <- -colik(bdts
			  , x
			  , dm
			  , y0
			  , t0 = t0
			  , ...
			)
		if (trace) print(c(exp=exp(x), -rv))
		rv - ln.dprior(x)
	}
		
	if (cotype=='mscom'){
		of <- .of.mscom
	} else if(cotype=='com12'){
		of <- .of.com12
	} else{
		stop('Coalescent model (cotype) must be a supported type')
	}
	
	fit0 <-  optim( par = theta_start
	 , fn = of
	 , ...
	 , method='Nelder-Mead'
	 , control = list( reltol = 1e-6, trace=trace, maxit = 1e3)
	)
	fit1 <- NA
	.dm <- dm
	if  ((np+1) > tail(np_range,1) ){ 
		done <- TRUE
		np <- np + 1 # for output
	}
	while(!done){
		np <- np + 1
		.dm <- dm
		#dm <- gen.demo.model3( t0, bdt, npoints = np, model_gamma = death_rate)
		if (!est_death_rate){
			dm <- gen.demo.model3( t0, bdt, npoints = np, model_gamma = death_rate_guess)
		} else{
			dm <- gen.demo.model3( t0, bdt, npoints = np)
		}
		sp_start <- .gen.start.splinePoint.parms(np)
			if (cotype=='mscom'){
				a_start <- log( R0guess * death_rate_guess )
				anames <- paste( sep='', 'akima', 1:np)
				theta_start <- c( setNames( rep(a_start , np ), anames) 
				 , sp_start
				 , lnNe = log(Ne_guess)
				 , lny0 = log(1)
				)
			} else if( cotype == 'com12'){
				a_start <- log( R0guess * death_rate_guess )
				anames <- paste( sep='', 'akima', 1:np)
				theta_start <- c( setNames( rep(a_start , np ), anames) 
				 , sp_start
				 , lny0 = log(1)
				)
			}
			if (est_death_rate){
				theta_start <- c( theta_start, lngamma = unname(log(death_rate_guess)) )
			}
		
		fit1 <-  optim( par = theta_start
			 , fn = of
			 , ...
			 , method='Nelder-Mead'
			 , control = list( reltol = 1e-6, trace=0, maxit=1e3)
			)
		if ( !.lrt.accept.h1(fit0, fit1) ) {
			done <- TRUE
		} else{
			fit0 <- fit1
			cat (paste('number splines points', np-1, 'rejected... incrementing...\n'))
			if  ((np+1) > tail(np_range,1) ){ done <- TRUE}
		}
	}
	cat('Number of spline points: \n')
	print(np-1)
	y0 <- exp( fit0$par['lny0'])
	x <- fit0$par
	if (est_death_rate){
		x <- c( x, gamma = unname(exp(x['lngamma'])))
		death_rate <-  x['gamma']
	} else{
		death_rate <- death_rate_guess
	}
	tfgy <- tryCatch( {.dm(x, y0, t0, bdt$maxSampleTime, res = 1e2) }, error = function(e) browser() )
	demo.history <- tryCatch( data.frame(times = tfgy$times, pop.size=tfgy$y, reproduction.number= tfgy$f / death_rate/ tfgy$y)
	 , error = function(e) NA )
	list( fit = fit0 
	 , par = x
	 , demo.model = .dm 
	 , demo.history = demo.history
	 , Ne = ifelse( 'lnNe' %in% names(fit0$par) , unname(exp(fit0$par['lnNe'] )), NA)  
	 , numberSplinePoints = np - 1
	 , bdts = bdts
	 , bdt = bdt
	 , t0 = t0
	 , death_rate_guess = death_rate_guess
	 , cotype = cotype
	 , R0guess = R0guess
	 , Ne_guess = Ne_guess
	 , y0_guess = y0_guess
	 , priors = priors
	 , est_death_rate = est_death_rate
	 #, additional_parameters = ... 
	)
}

parboot.skyspline.mle <- function(fit, nreps = 2e2, tfin = NULL, ...)
{
	x0 <- unname(exp(fit$fit$par['lny0']))
	if (fit$est_death_rate){
		g <- unname( exp( fit$fit$par['lngamma'] ))
		fit$fit$par['gamma'] <- g
	} else{
		g <- fit$death_rate_guess
		fit$fit$par['gamma'] <- g
	}
	
	if ('lnNe' %in% names(fit$fit$par)){
		Ne <- unname( exp(fit$fit$par['lnNe']) ) 
	} else{
		Ne <- (1/g) / 1e2 #small
	}
	t0 <- fit$t0
	theta <- c()
	for (k in 1:nreps)
	{
		sbdt <- tryCatch({
			sim.Mmodel0(fit$bdt$sampleTimes, fit$fit$par, fit$demo.model
			  , x0 = x0
			  , t0 = t0
			  , Ne = Ne #wh, scalar
			  , maxk = 10 # approx max lines in hosts
			  , res = 1e3
			)
		}, error = function(e) NA)
		if (!any(is.na(sbdt))){
			.fit <- tryCatch({ fit.skyspline.mle3(sbdt
			  , t0
			  , death_rate_guess = fit$death_rate_guess
			  , R0guess = fit$R0guess
			  , cotype=fit$cotype
			  , Ne_guess = fit$y0_guess
			  , y0_guess = fit$y0_guess
			  , np_range = fit$numberSplinePoints
			  , priors = fit$priors
			  , ... # passed to likelihood
			)}, error = function(e) list(fit=list(par=NULL)) #list(fit=list(par=rep(NA, ncol(theta))))
			)
			theta <- rbind( theta, .fit$fit$par )
		}
		print( paste( 'rep', k, 'complete', date()))
	}
	
	vcv <- cov( theta )
	CIs <- setNames( lapply( colnames(theta), function(pn){
		dx <- sqrt( vcv[pn,pn] ) * 1.96
		c( fit$fit$par[pn] - dx, fit$fit$par[pn] + dx )
	}), colnames( theta ))
	
	## traj
	Ys <- c()
	Rs <- c()
	cumF <- c()
	require(mvtnorm)
	if ('gamma' %in% names( fit$par )){
		mpar <- fit$par[-which(names(fit$par) == "gamma")]
	} else{
		mpar <- fit$par 
	}
	rmvnorm_theta <- rmvnorm(nreps, mean = mpar, sigma =vcv) #NOTE fit$par may include gamma
	if (any(is.na(rmvnorm_theta))){
		rmvnorm_theta <- rmvnorm(nreps, mean = mpar, sigma =diag(diag(vcv))) #NOTE fit$par may include gamma
	}
	if (any(is.na(rmvnorm_theta))){
		stop('Parboot failed; NA values in fit.')
	}
	for (k in 1:nreps){
		.theta <- setNames( rmvnorm_theta[k, ], colnames(theta))
		if (fit$est_death_rate){
			.theta <- c( .theta, gamma = unname(exp( .theta['lngamma'])))
		} 
		.x0 <- unname( exp( .theta['lny0']))
		tfgy <- fit$demo.model( .theta, .x0, t0, fit$bdt$maxSampleTime , res = 1e2, tfin = tfin)
		Ys <- cbind( Ys, tfgy$y )
		if (fit$est_death_rate){
			Rs <- cbind(Rs, tfgy$f / tfgy$y / .theta['gamma']  ) 
		} else{
			Rs <- cbind(Rs, tfgy$f / tfgy$y / fit$death_rate_guess ) 
		}
		times <- tfgy$t
		tstep <- abs(times[1] - times[2] )
		cumF <- cbind( cumF 
		 , rev( cumsum(tstep * rev(tfgy$f)))
		)
	}
	list( CIs = CIs
	 , population_size = sapply( 1:nrow(Ys), function(k) quantile( Ys[k,], probs = c(.5, .025, .975) ) )
	 , R.t = sapply( 1:nrow(Rs), function(k) quantile( Rs[k,], probs = c(.5, .025, .975) ) )
	 , cumulative_births = sapply( 1:nrow(cumF), function(k) quantile( cumF[k,], probs = c(.5, .025, .975) ) )
	 , times = times
	 , VCV = vcv
	 , theta = rmvnorm_theta
	)
}
#######################################
skyspline.metrop.hastings3 =  skyspline.metrop.hastings<- function(
	  bdts
	  , nsteps = 1e4
	  , sd_prop_lnR = .25 # scale for mvnorm proposal of R
	  , sd_prop_logDeathRate = NA # user must define
	  , sd_prop_Ne = 1/10
	  , sd_prop_lny0 = 1/4
	  , sd_prop_splinePoint = 1/2
	  , start_lnR = log(2)
	  , start_log_deathRate = NA # user must define
	  , start_log_Ne = NA # user must define 
	  , start_log_y0 = log(1)
	  , death_rate_logprior = function(x) stop('Informative prior must be defined')
	  , R0_logprior = function(x) dlnorm(x, log(2), 1, log=T)
	  , Ne_logprior = function(x) 0# log(1/x)
	  , y0_logprior = function(x) dexp(1, rate=1, log = T)
	  , splinePoint_logprior = function(x) dbeta(x, 2,2,log=TRUE)
	  , cotype=c('com12', 'mscom')
	  , numberSplinePoints = 4
	  , t0 = 0
	  , reportFreq = 100
	  , thin = 100
	  , ...
){
	if (class(bdts)[1]=='DatedTree'){
		bdt = bdts
	} else if (class(bdts)[1]=='list'){
		bdt <- bdts[[1]]
	} else{
		stop('bdts must be DatedTree or list of DatedTree')
	}
	done <- FALSE
	np <- numberSplinePoints
	dm <- gen.demo.model3( t0, bdt, npoints = np)
		
	akimaNames <- paste( sep='', 'akima', 1:np)
	sp_start <- .gen.start.splinePoint.parms(np)
	theta_start <- c( setNames(log(rep( exp(start_lnR) * exp(start_log_deathRate) , np)), akimaNames)
	  , sp_start
	  , lny0 = unname( start_log_y0)
	  , lngamma = unname( start_log_deathRate )
	) 
	if (cotype == 'mscom'){
		theta_start <- c( theta_start, lnNe = unname(start_log_Ne))
	}
	
	log.prior.theta <- function(theta)
	{
		lpt <- death_rate_logprior(exp(theta['lngamma'])) + 
		y0_logprior(exp(theta['lny0'])) + 
		sum(sapply( akimaNames, function(an) {
			R <- exp(theta[an]) / exp(theta['lngamma'])
			R0_logprior( R )
		} )) + 
		ifelse( cotype == 'mscom', Ne_logprior(exp(theta['lnNe'])) , 0)
		
		if (length(sp_start) > 0){ 
			lpt <- lpt + sum(sapply(names(sp_start), function(spname){
				x <- exp(theta[spname]) / ( 1 + exp(theta[spname]))
				splinePoint_logprior( x )
			})) 
		}
		lpt
	}
	#TODO
	q.theta <- function(theta, istep){
		i <- istep %% length(theta) + 1
		pname <- names(theta)[i]
		if ( grepl('akima', pname)){
			theta[pname] <- theta[pname] + rnorm(1, 0, sd_prop_lnR)
		} else if (pname=='lny0') {
			theta[pname] <- theta[pname] + rnorm(1, 0, sd_prop_lny0 )
		} else if (pname=='lngamma') {
			theta[pname] <- theta[pname] + rnorm(1, 0, sd_prop_logDeathRate)
		} else if( pname == 'lnNe'){
			theta[pname] <- theta[pname] + rnorm(1, 0, sd_prop_Ne )
		} else if (grepl('splinePoint', pname)){
			theta[pname] <- theta[pname] + rnorm( 1, 0, sd_prop_splinePoint)
		}
		theta
	}
	
	.of.mscom <- function( x, ...)
	{
		Ne <- exp(x['lnNe'])
		y0 <- exp(x['lny0'])
		x['gamma'] <-  unname( exp(x['lngamma']))
		rv <- colik.mscom(bdts
			  , x
			  , dm
			  , y0
			  , t0 = t0
			  , Ne = Ne
			  , maxk=10
			  , ...
			)
		#print(c(exp=exp(x), -rv))
		unname(rv)
	}
	
	.of.com12 <- function( x, ...)
	{
		y0 <- exp(x['lny0'])
		x['gamma'] <-  unname( exp(x['lngamma']))
		rv <- colik(bdts
			  , x
			  , dm
			  , y0
			  , t0 = t0
			  , ...
			)
		#print(c(exp=exp(x), -rv))
		unname(rv)
	}
	
	if (cotype=='mscom'){
		of <- .of.mscom
	} else if(cotype=='com12'){
		of <- .of.com12
	} else{
		stop('Coalescent model (cotype) must be a supported type')
	}
	
	theta <- theta_start
	theta_trace <- c()# matrix(NA, nrow = length(theta), ncol = nsteps)
	of_trace <- c()
	o_trace <- c()
#~ browser()
	
	oprior <- log.prior.theta(theta)
	ofval <- of(theta, ...) 
	o <- ofval + oprior
	
	nAccept <- 0
	for (istep in 1:nsteps){		
		.theta <- q.theta(theta, istep)
		
		.oprior <- log.prior.theta(.theta)
		.ofval <- of(.theta, ...) 
		.o <- .ofval + .oprior
		
		.o <- ifelse( is.na(.o), -Inf, .o)
		
		if ( (istep %% reportFreq) == 0){
			cat('istep \n', istep, '\n')
			print( theta)
			print(.theta)
			cat( 'likelihood\n')
			print(unname(c(o, .o)))
			cat('prior\n')
			print(unname(c(oprior, .oprior)))
			cat('objective\n')
			print(unname(c(ofval, .ofval)))
			cat('acceptance probability\n')
			print(nAccept / istep)
			cat('################################################\n\n')
		}
		
#~ x <- tryCatch({exp(.o - o) > runif(1)}, error = function(e) browser() )
		if ( exp(.o - o) > runif(1) )
		{
			theta <- .theta
			o <- .o
			ofval <- .ofval
			oprior <- .oprior
			nAccept <- nAccept + 1
		}
		if ( (istep %% thin) ==0){
			#theta_trace[, istep] <- theta
			theta_trace <- cbind( theta_trace, theta)
			of_trace <- c( of_trace, ofval )
			o_trace <- c( o_trace, o)
		}
	}
	theta_trace <- t(theta_trace)
	colnames(theta_trace) <- names(theta)
	list( theta = theta_trace, objFun = of_trace, logLik=o_trace, acceptProb = nAccept / nsteps, demo.model =dm )
}
