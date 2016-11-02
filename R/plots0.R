require(ggplot2)

plot.mcmc.population.size <- function(omh, t0, t1, burnin_percent = 20, sample_size = 2e2, res = 1e2
 , quantiles = c(.025, .975)
 , log = TRUE
 , ...
)
{
	# plot size, cumulative infections, and R(t)
	pnames <- colnames(omh$theta)
	i <- sample( floor(nrow(omh$theta)*(burnin_percent/100)):nrow(omh$theta), replace=F, size = sample_size)
	thetas <- omh$theta[i, ]
	
	tfgys <- lapply( 1:nrow(thetas), function(k){
		.theta <- thetas[k,]
		.theta['gamma'] <- unname(exp(.theta['lngamma']))
		y0 <- unname(exp(.theta['lny0']))
		omh$demo.model(.theta, y0, t0, t1, res =res )
	})
	
	Y <- (sapply( tfgys, function(tfgy) tfgy$y))
	Yq <- t(sapply(1:nrow(Y), function(k) quantile( Y[k,] , prob = c(.5, quantiles))))
	t <- tfgys[[1]]$t
	
	yplot <- 
	qplot( t, Yq[, 1], xlab = 'Time', ylab = 'Population size' , geom='path', ...) + 
	  geom_ribbon( aes( x = t, ymin = Yq[,2], ymax = Yq[,3]), fill='blue', alpha = .15  ) 
	
	if (log){
		yplot <- yplot + scale_y_log10()
	}
	yplot
}

plot.mcmc.R.t <- function(omh, t0, t1, burnin_percent = 20, sample_size = 2e2, res = 1e2
 , quantiles = c(.025, .975)
 , log = FALSE
 , ...
)
{
	# plot size, cumulative infections, and R(t)
	pnames <- colnames(omh$theta)
	i <- sample( floor(nrow(omh$theta)*(burnin_percent/100)):nrow(omh$theta), replace=F, size = sample_size)
	thetas <- omh$theta[i, ]
	
	Rs <- sapply( 1:nrow(thetas), function(k){
		.theta <- thetas[k,]
		.theta['gamma'] <- unname(exp(.theta['lngamma']))
		y0 <- unname(exp(.theta['lny0']))
		tfgy <- omh$demo.model(.theta, y0, t0, t1, res =res )
		unname(rev( tfgy$f / tfgy$y / .theta['gamma']) )
	})
	
	xq <- t(sapply(1:nrow(Rs), function(k) quantile( Rs[k,] , prob = c(.5, quantiles))))
	t <- seq(t0, t1, length.out = res )
	
	qp <- 
	qplot( t, xq[, 1], xlab = 'Time', ylab = 'Reproduction number' , geom='path', ...) + 
	  geom_ribbon( aes( x = t, ymin = xq[,2], ymax = xq[,3]), fill='blue', alpha = .15  ) 
	
	if (log){
		qp <- qp + scale_y_log10()
	}
	qp
}



plot.mcmc.cumulative.births <- function(omh, t0, t1, burnin_percent = 20, sample_size = 2e2, res = 1e2
 , quantiles = c(.025, .975)
 , log = TRUE
 , ...
)
{
	# plot size, cumulative infections, and R(t)
	pnames <- colnames(omh$theta)
	i <- sample( floor(nrow(omh$theta)*(burnin_percent/100)):nrow(omh$theta), replace=F, size = sample_size)
	thetas <- omh$theta[i, ]
	
	t <- seq(t0, t1, length.out = res )
	X <- sapply( 1:nrow(thetas), function(k){
		.theta <- thetas[k,]
		.theta['gamma'] <- unname(exp(.theta['lngamma']))
		y0 <- unname(exp(.theta['lny0']))
		tfgy <- omh$demo.model(.theta, y0, t0, t1, res =res )
		with(tfgy,{
			f <- rev(f)
			dx <- t[2] - t[1]
			cumsum(f) * dx
		})
	})
	
	xq <- t(sapply(1:nrow(X), function(k) quantile( X[k,] , prob = c(.5, quantiles))))
	
	qp <- qplot( t, xq[, 1], xlab = 'Time', geom='path', ...) + 
	  geom_ribbon( aes( x = t, ymin = xq[,2], ymax = xq[,3]), fill='blue', alpha = .15  ) 
	
	if (log){
		qp <- qp + scale_y_log10()
	}
	qp
}



plot.pbfit <- function (pbfit,log = TRUE, ...) 
{
    Yq <- t(pbfit$population_size)
    t <- pbfit$times
    yplot <- qplot(t, Yq[, 1], xlab = "Time", ylab = "Population size", 
        geom = "path", ...) + geom_ribbon(aes(x = t, ymin = Yq[, 
        2], ymax = Yq[, 3]), fill = "blue", alpha = 0.15)
    if (log) {
        yplot <- yplot + scale_y_log10()
    }
    yplot
}
