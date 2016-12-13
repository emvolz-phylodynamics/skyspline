
################################################################
# DatedTree
DatedTree <- function( phylo
  , sampleTimes
  , tol = 1e-6
  , minEdgeLength = 0
  , roundEdgeLengthDown = 0)
{
	if (is.null(names(sampleTimes))) stop('sampleTimes vector must have names of tip labels')
	
	# resolve any multifurcations 
	phylo <- tryCatch( { multi2di( phylo ) }, error = function(e) phylo )
	
	phylo$sampleTimes <- sampleTimes[phylo$tip.label]
	
	if (roundEdgeLengthDown > 0 ){
		phylo$edge.length[ phylo$edge.length < roundEdgeLengthDown ] <- 0
	}
	
	phylo$n = n <- length(sampleTimes)
	Nnode <- phylo$Nnode
	# compute heights, ensure consistency of sample times and branch lengths
	phylo$maxSampleTime   <- max(phylo$sampleTimes)
	heights <- rep(NA, (phylo$Nnode + length(phylo$tip.label)) )
	heights[1:length(phylo$sampleTimes)] <- phylo$maxSampleTime - phylo$sampleTimes
	curgen <- 1:length(phylo$sampleTimes)
	edgeLengthChange <- TRUE 
	while (edgeLengthChange)
	{
		edgeLengthChange <- FALSE
		while( length(curgen) > 0) { 
			nextgen <- c()
			icurgenedges <- which(  phylo$edge[,2] %in% curgen  )
			for (i in icurgenedges){
				u<- phylo$edge[i,1]
				v<- phylo$edge[i,2]
				if (!is.na(heights[u])){ # ensure branch lengths consistent
					if ( heights[u] > 0 & abs(heights[u] - (phylo$edge.length[i] + heights[v])) > tol )
					{ #
					  stop( 'Tree is poorly formed. Branch lengths incompatible with sample times.')
					} else if ( 0!=(heights[u] - (phylo$edge.length[i] + heights[v]) ) ){
						edgeLengthChange <- TRUE 
					}
					phylo$edge.length[i] <- max(0, max(minEdgeLength, heights[u] - heights[v] ) )
					if ( phylo$edge.length[i] < roundEdgeLengthDown ){
						phylo$edge.length[i] <- 0
						edgeLengthChange <- TRUE 
					}
					heights[u] <- heights[v]  + phylo$edge.length[i]
				} else{
					heights[u] <- phylo$edge.length[i] + heights[v]
				}
				
				nextgen <- c(nextgen, u)
			}
			curgen <- unique(nextgen)
		}
	}
	phylo$heights <- heights
	phylo$maxHeight <- max(phylo$heights)
	#phylo$heights <- signif( phylo$heights, digits = floor( 1 / phylo$maxHeight /10 )  +  6 ) #
	phylo$parentheights = phylo$parentheight <- sapply( 1:(n+Nnode), function(u){
		i <- which( phylo$edge[,2]== u)
		if (length(i)!=1) return( NA )
		phylo$heights[ phylo$edge[i,1] ]
	})
	
	phylo$root <- which.max( phylo$heights)
	
	ix <- sort( phylo$sampleTimes, decreasing = TRUE, index.return=TRUE)$ix
	phylo$sortedSampleHeights <- phylo$maxSampleTime - phylo$sampleTimes[ix]
	# parents and daughters 
	phylo$parent = phylo$parents <- sapply(1:(phylo$n+phylo$Nnode), function(u) {
		i <-  which(phylo$edge[,2]==u)
		if (length(i) == 0) return (NA)
		a <- phylo$edge[i ,1]
		if (length(a) > 1) stop("Tree poorly formed; node has more than one parent.")
		a
	}) 
	phylo$daughters <- t(sapply( 1:(phylo$n+phylo$Nnode), function(a){
		uv <- phylo$edge[which(phylo$edge[,1]== a),2]
		if (length(uv)==0) uv <- c(NA, NA)
		#if (length( uv)!=2) print( uv )
		uv
	}))
	
	class(phylo) <- c("DatedTree", "phylo")
	phylo
}



################################################################


.lambda <- function(A, Ne, g_x){
	#~ 	 lambda = (A/(2*Ne)) * g'' / g'
	mk <- length(g_x)
	K <- 1:mk
	gpp <- sum( K[2:mk] * K[1:(mk-1)]  * g_x[2:mk] ) #k(k-1)pk
	gp <- sum( K * g_x) # k pk
	(A/(2*Ne)) * gpp / gp
}

#~ .update.g_x.node <- function( g_x, A )
#~ {
#~ 	mk <- length(g_x)
#~ 	K <- 1:mk
#~ 	gpp <- sum( K[2:mk] * K[1:(mk-1)]  * g_x[2:mk] ) #k(k-1)pk
#~ 	gp <- sum( K * g_x) # k pk
#~ 	w <- gpp*A - gpp ...
#~ }

# note this version does not have sampled ancestors
.update.g_x.sample <- function( g_x, A)
{
	maxk <- length(g_x)
	m1 <- sum( 1:maxk * g_x )
	B <- A / m1
	g_x[1] <- (1/(B+1)) * (B*g_x[1] + 1 )
	g_x[2:maxk] <- g_x[2:maxk] * B / (B + 1 )
	g_x / sum( g_x )
}

################################################################################
.mean.log.x <- function(x)
{
	mx <- max(x )
	if (is.infinite( mx)) return(mx)
	log( sum( exp( x - mx)) ) + mx
}

################################################################################
colik.mscom <- function(trees
 , theta
 , demographic.process.model
 , x0
 , t0
  , Ne #wh, scalar
  , maxk # approx max lines in hosts
  , res = 1e3
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , logd.sim.prior = NULL #function(tfgy, theta)
) {
	if (class(trees)[1]=='DatedTree'){
		tree <- trees
	} else if(class(trees)[1]=='list'){
		tree <- trees[[1]]
	} else{
		stop('trees must be DatedTree or list of DatedTree')
	}
	if ( tree$maxHeight >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res) 
	if (class(trees)[1]=='DatedTree'){
		return(colik.mscom.tfgy( trees, tfgy, theta, Ne, maxk, res, timeOfOriginBoundaryCondition, maxHeight, forgiveAgtY, AgtY_penalty, returnTree, logd.sim.prior=logd.sim.prior))
	}
	else if (class(trees)[1]=='list'){
		lls <- sapply( trees, function(tre){
			colik.mscom.tfgy( tre, tfgy, theta, Ne, maxk, res, timeOfOriginBoundaryCondition, maxHeight, forgiveAgtY, AgtY_penalty, returnTree, logd.sim.prior=logd.sim.prior)
		})
		return(.mean.log.x( lls ) )
	}
}
colik.mscom.tfgy <- function(tree
  , tfgy # list with compoents times, f, y; time in decreasing order
  , theta
  , Ne #wh, scalar
  , maxk # approx max lines in hosts
  , res = 1e3
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = .2 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 10 # penalises likelihood if A > Y
  , returnTree = FALSE
  , logd.sim.prior = NULL #function(tfgy, theta)
) 
{
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	#tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res) 
	# list (heights, f, y)
	if (tfgy[[1]][1] < tfgy[[1]][2]) stop('tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample')
	t0 <- tail( tfgy[[1]], 1 )
	
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		ih <- min( length(tfgy$times), 1+floor( length(tfgy$times) * h / (tree$maxSampleTime- t0)  ) )
		list( f = unname(tfgy$f[ih])
		 , y = unname(tfgy$y[ih])
		)
	}
	
	if (is.null( tree$n ) ) tree$n <- length( tree$sampleTimes)
		
	# 
	eventTimes <- unique( sort(tree$heights) )
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]
	}
	L <- 0
	
	extantAtEvent_nodesAtHeight <- eventTimes2extant( eventTimes, tree$heights, tree$parentheight ) #
	extantAtEvent_list <- extantAtEvent_nodesAtHeight[[1]]
	nodesAtHeight <- extantAtEvent_nodesAtHeight[[2]]
	
	#variables to track progress at each node
	tree$A <- c() #h->A
	tree$lnS <- c()
	tree$lnr <- c()
	tree$ih <- c()
	tree$g_x0 <- c()
	tree$g_x <- c()
	tree$tfgy <- tfgy 
	
	loglik <- 0
	if (!is.null(logd.sim.prior)){
		loglik <- logd.sim.prior(tfgy, theta)
		#cat('sim prior value:\n')
		#print(loglik)
	}
	
	g_x <- rep(0, maxk)
	g_x[1] <- 1 # initially all hosts have single lineage
	for (ih in 1:(length(eventTimes)-1))
	{
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
		
		#get A0, process new samples, calculate state of new lines
		extantLines <- extantAtEvent_list[[ih]]
		A0 = A <- length(extantLines)
		Lambda_g <- solve_dgdt0( tfgy$times, tfgy$f, tfgy$y
		  , g_x
		  , L
		  , h0
		  , h1
		  , A
		  , Ne
		  , tree$maxSampleTime - t0 #tree$maxHeight # should match time axis on tfgy
		)
		L <- Lambda_g$Lambda
		g_x <- pmax(as.vector(Lambda_g$g), 0)
		g_x <- g_x / sum(g_x)
		
		# clean output
		if (is.nan(L)) {L <- Inf}
		if (sum(is.nan(A)) > 0) A <- A0
		
		#if applicable: update ustate & calculate lstate of new line
		newNodesAndSamples <- nodesAtHeight[[ih+1]]
		newNodes <- newNodesAndSamples[newNodesAndSamples > tree$n] 
		newSamples <- setdiff( newNodesAndSamples, newNodes )
		
		m1 <- sum( 1:maxk * g_x )
		B <- A / m1
		
		# incorporate samples, update g
		.A <- A + 1 - 1
		for (u in newSamples){
			g_x <- .update.g_x.sample(g_x, .A)
			.A <- .A + 1
		}
		
		{
			YmA <- (sum(fgy$y) - B)
			if (any(is.na(YmA))) {
				if (returnTree){
					return(list( loglik = -Inf
					 , tree = tree ))
				} else{
					return(-Inf)
				}
			}
			if (YmA < 0){
				if (length(extantLines)/length(tree$tip.label)  > forgiveAgtY){
					L <- Inf
				} else{
					L <- L + L * abs(YmA) * AgtY_penalty
				}
			}
			
			for (alpha in newNodes)
			{
				tree$g_x0 <- rbind( tree$g_x0, g_x)
				lambda_gx_list <- update_lambda_gx_node2( A,  Ne,  g_x, fgy$y) #update aug 26
				#tree$coalescentRates[alpha] <- .lambda(A, Ne, g_x)
				tree$coalescentRates[alpha] <- lambda_gx_list$lambda
				
				if (is.nan(L))
				{
					warning('is.nan(L)')
					L <- (h1 - h0) * tree$coalescentRates[alpha] 
				}
				
				# trace
				tree$A <- rbind( tree$A, A)
				tree$lnS <- c( tree$lnS, -L)
				tree$lnr <- c( tree$lnr, log(tree$coalescentRates[alpha]) )
				tree$ih <- c( tree$ih, h1)
				tree$g_x <- rbind( tree$g_x, g_x )
				
				# update lik
				loglik <- loglik + log( tree$coalescentRates[alpha] ) - L ;
				if (is.infinite( loglik)){
					if (returnTree){
						return(list( loglik = loglik
						 , tree = tree ))
					} else{
						return(loglik)
					}
				}
			} 
			
			if (length(newNodes) > 0){
				g_x <- as.vector( lambda_gx_list$g_x ) # a resolution of problem of concurrent nodes
				# not a great resolution: g_x should be modified more than once (for each node), but doing that makes small likelihoods
				# since survival terms are also shared, I am going with this approx
				L <- 0
			}
		}
	}
	
	if (returnTree){
		return(list( loglik = loglik
		 , tree = tree ))
	}
	return( loglik )
}


sim.Mmodel0 <- function(sampleTimes, theta, demographic.process.model, x0, t0
  , Ne #wh, scalar
  , maxk # approx max lines in hosts
  , res = 1e3
  , integrationMethod='lsoda'
) {
	if (is.null(names(sampleTimes))){
		names(sampleTimes) <- paste('t', 1:length(sampleTimes), sep='')
	}
	ddt <- (max(sampleTimes) - t0)/res # default step
	
	maxSampleTime <- max(sampleTimes)
	
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	tfgy <- demographic.process.model( theta, x0, t0, maxSampleTime, res = res) 
		
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		ih <- min( length(tfgy$times), 1+floor( length(tfgy$times) * h / (maxSampleTime- tail(tfgy$times,1))  ) )
		list( f = unname( tfgy$f[ih] )
		 , y = unname( tfgy$y[ih] )
		)
	}
	
	mst <- max(sampleTimes)
	ssh <- mst - sort( sampleTimes , decreasing = TRUE)
	n <- length(ssh)
	samplesAdded <- 0
	nodesAdded <- 0
	
	# 
	edge <- matrix( NA, nrow = (n + n -2 ), ncol =2) # no root edge
	edge.length <- rep(NA, n + n - 2)
	edgesAdded <- 0
	heights <- rep(NA, n + n - 1)
	
	#
	L <- 0
	
	extant <- rep(FALSE, n + n -1 )
	A <- 0
		
	g_x <- rep(0, maxk)
	g_x[1] <- 1 # initially all hosts have single lineage
	
	h0 <- 0
	h1 <- 0
	nextSampleHeight <- 0
	while( samplesAdded < n | nodesAdded < (n-1) )
	{
		h0 <- h1
		h1 <- h0 + ddt
		
		# incorp samples at base of interval
		if ( h0 >= nextSampleHeight){
			samplesAdded <- samplesAdded + 1
			extant[samplesAdded] <- TRUE
			heights[samplesAdded] <- h0
			A <- A + 1
			
			# next sample height
			if ( samplesAdded < n ){
				nextSampleHeight <- ssh[ samplesAdded + 1]
			} else{
				nextSampleHeight <- Inf
			}
			
			#update g_x
			g_x <- .update.g_x.sample(g_x, A)
			
			# update end of interval
			h1 <- min( h0 + ddt, nextSampleHeight )
		}
		
		#if (h1 > (mst - t0)) break; 
		if (h1 > h0 ){
			L <- 0
			Lambda_g <- solve_dgdt0( tfgy$times, tfgy$f, tfgy$y
			  , g_x
			  , L
			  , h0
			  , h1
			  , A
			  , Ne
			  , mst - t0  # maxheight
			)
			L <- Lambda_g$Lambda
			g_x <- pmax(as.vector(Lambda_g$g), 0)
			g_x <- g_x / sum(g_x)
			
			# add any new nodes at h1
			numNewNodes <- rpois(1,  L )
			if (numNewNodes > 0){
				fgy <- get.fgy(h1)
				for (k in 1:numNewNodes ){
					#update g_x
					lambda_gx_list <- update_lambda_gx_node( A,  Ne,  g_x, fgy$y)
					g_x <- as.vector( lambda_gx_list$g_x )
					
					# pick two at random
					uv <- sample( which(extant), size = 2, replace = FALSE)
					nodesAdded <- nodesAdded + 1
					a <- nodesAdded + n
					
					#update extant and heights
					heights[a] <- h1
					extant[a] <- TRUE
					extant[uv[1]] <- FALSE
					extant[uv[2]] <- FALSE
					A <- A - 1
					
					# add edges
					edgesAdded <- edgesAdded + 1
					edge[ edgesAdded ,1] <- a
					edge[ edgesAdded ,2] <- uv[1]
					edge.length[edgesAdded] <- heights[a] - heights[uv[1]]
					
					edgesAdded <- edgesAdded + 1
					edge[ edgesAdded ,1] <- a
					edge[ edgesAdded ,2] <- uv[2]
					edge.length[edgesAdded] <- heights[a] - heights[uv[2]]
				}
			}
		}
	}
	
	if (nodesAdded <(n-1)){
		#browser()
		# make polytomy?
		NA
	}
	
	tre <- list( edge = edge, edge.length = edge.length
	 , Nnode = n - 1
	 , tip.label = names(ssh) )
	class(tre) <- 'phylo'
	tre <- read.tree( text = write.tree(tre ))
	DatedTree( tre, sampleTimes , tol = ddt)
}


########################################################################
.solve_CoM12L_deSolve <- function(times, fs, ys, L, h0, h1, A, maxSampleTime) 
{
	.dLdt <- function(t, y, parms, ... ){
		ih <- min( length(times), 1+floor( length(times) * t / (maxSampleTime- tail(times,1))  ) )
		ff <- fs[ih]
		yy <- ys[ih]
		list( (A * (A-1)/2) * 2 * ff / yy^2 )
	}
	o <- ode( y = L, times = c(h0, h1) , func = .dLdt, method = 'lsoda')
	o[2,2] 
}

########################################################################
# CoM12 model (unstructured)
colik <- function(trees
  , theta
  , demographic.process.model, x0, t0, res = 1e3
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = .2 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 10 # penalises likelihood if A > Y
) {
	if (class(trees)[1]=='DatedTree'){
		tree <- trees
	} else if(class(trees)[1]=='list'){
		tree <- trees[[1]]
	} else{
		stop('trees must be DatedTree or list of DatedTree')
	}
	if ( tree$maxHeight >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res) 
	if (class(trees)[1]=='DatedTree'){
		return(colik.tfgy( trees, tfgy
		 , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
		 , maxHeight = maxHeight 
			  , forgiveAgtY = forgiveAgtY
			  , AgtY_penalty = AgtY_penalty))
	}
	else if (class(trees)[1]=='list'){
		lls <- sapply( trees, function(tre){
			colik.tfgy( tre
			  , tfgy
			  , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
			  , maxHeight = maxHeight 
			  , forgiveAgtY = forgiveAgtY
			  , AgtY_penalty = AgtY_penalty
			)
		})
		return(.mean.log.x( lls ) )
	}
}
colik.tfgy  <- function(tree
  , tfgy
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
) 
{
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	if (tfgy[[1]][1] < tfgy[[1]][2]) stop('tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample')
	t0 <- tail( tfgy[[1]], 1 )
	
	g_hres <- length(tfgy$times)
	g_treeT <- tfgy$times[1] - tail(tfgy$times,1)
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		#ih <- min( length(tfgy$times), 1+floor( length(tfgy$times) * h / (tree$maxSampleTime- tail(tfgy$times,1))  ) )
		ih <- 1+ max(0, min( g_hres-1, floor(g_hres*h/g_treeT ) ))
		list( f = unname( tfgy$f[ih] )
		 , y = unname( tfgy$y[ih] )
		)
	}
	
	if (is.null( tree$n ) ) tree$n <- length( tree$sampleTimes)
	eventTimes <- unique( sort(tree$heights) )
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]
	}
	S <- 1
	L <- 0
	
	extantAtEvent_nodesAtHeight <- eventTimes2extant( eventTimes, tree$heights, tree$parentheight ) #1-2 millisec
	extantAtEvent_list <- extantAtEvent_nodesAtHeight[[1]]
	nodesAtHeight <- extantAtEvent_nodesAtHeight[[2]]
	#variables to track progress at each node
	tree$A <- c() #h->A
	tree$lnS <- c()
	tree$lnr <- c()
	tree$ih <- c()
	loglik <- 0
	for (ih in 1:(length(eventTimes)-1))
	{
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
		#fgy <- get.fgy(h0)
		
		#get A0, process new samples, calculate state of new lines
		extantLines <- extantAtEvent_list[[ih]]
		A0 = A <- length(extantLines)
		
		# quick boost solver: 
		L <- solve_CoM12L(tfgy$t, tfgy$f, tfgy$y, L, h0, h1, A) 
		#L <- .solve_CoM12L_deSolve(tfgy$t, tfgy$f, tfgy$y, 0, h0, h1, A, tree$maxSampleTime) 
		
		# clean output
		if (is.na(L)) {return(-Inf)}
		
		#if applicable: update ustate & calculate lstate of new line
		newNodes <- nodesAtHeight[[ih+1]]
		newNodes <- newNodes[newNodes > tree$n] 
		
		{
			YmA <- (sum(fgy$y) - length(extantLines))
			if (is.na(YmA)) return(-Inf)
			if (YmA < 0){
				if (length(extantLines)/length(tree$tip.label)  > forgiveAgtY){
					L <- Inf
				} else{
					L <- L + L * abs(YmA) * AgtY_penalty
				}
			}
			
			for (alpha in newNodes )
			{
				# NOTE if nodes are concurrent, the likelihood terms share L
				tree$coalescentRates[alpha] = corate <- ((A * (A - 1)) / 2) * 2 * fgy$f / max(1, fgy$y^2 - fgy$y)
				
				# trace
				tree$A <- rbind( tree$A, A)
				tree$lnS <- c( tree$lnS, -L )
				tree$lnr <- c( tree$lnr, log(corate) )
				tree$ih <- c( tree$ih, h1)
				
				# update lik
				loglik <- loglik + log( corate ) - L ;
				if (is.infinite( loglik)){
						return(loglik)
				}
			}  
		}
		
		# if coalescent occurred, reset cumulative hazard function
		if (length(newNodes) > 0) {
			L <- 0 
		}
	}
	return( loglik)
} 
