	

# Likelihood functions

Log.Lik = function(family, data, params, it = NULL) {
	
	if (is.null(it)) { it = ' ' }

	source('~/work/Infl_VN/MCMC_0501/Code/MCMC/Running_cpp_MCMC.R', local = TRUE)
	cpp.out = cpp.out[(nrow(cpp.out)-3649):nrow(cpp.out), ]
	cpp.agg = aggregate(cpp.out[ , c('inc1', 'inc2', 'inc3')], 
		by = list('week' = ceiling(c(1:(nrow(cpp.out))/7))), FUN = sum)
	cpp.agg = cpp.agg[1:520, ]
	
	scale.fac = params['report']
	
	
	Log.lik = 0
	for (cc in 1:ncol(data)) {
		if (family == 'poisson') {
			Lik = log(dpois(round(data[ , cc]), 
				lambda = cpp.agg[ , cc+1] * scale.fac))
		} else if (family == 'normal') {
			Lik = dnorm(data[ , cc], 
				mean = cpp.agg[ , cc+1] * scale.fac, 
				sd = sqrt(params['norm_var']),  
				log = TRUE)
		} else if (family == 'SSE') {
			Lik = -1 * ((data[ , cc] - 
				(cpp.agg[ , cc+1] * scale.fac))^2)
		} else if (family == 'gamma') {
			sc.mod = (cpp.agg[ , cc+1] + 0.01) * scale.fac
			Lik = dgamma(data[ , cc] + 0.01, 
				shape = sc.mod * sc.mod / params['norm_var'], 
				rate = sc.mod / params['norm_var'], log = TRUE
			)
		}
		
		
		Lik[which(Lik == Inf)] = 0
		Lik[which(Lik == -Inf)] = -2 * max(25, abs(Lik[which(is.finite(Lik))]), na.rm = TRUE)
		Lik = Lik[which(is.finite(Lik))]
		
		
		Log.lik = Log.lik + sum(Lik, na.rm = TRUE)
		
	}



	if (length(unique(cpp.agg[which(!is.na(data[ , 1])), 2])) == 1 ||
		length(unique(cpp.agg[which(!is.na(data[ , 2])), 3])) == 1 || 
		length(unique(cpp.agg[which(!is.na(data[ , 3])), 4])) == 1) {

		Log.lik = Log.lik - 100000

	} else if ((sum(!is.na(data[ , 1])) > 0 && max(cpp.agg[ , 2] * scale.fac, na.rm = TRUE) < 0.2) || 
		(sum(!is.na(data[ , 2])) > 0 && max(cpp.agg[ , 3] * scale.fac, na.rm = TRUE) < 0.2) ||
		(sum(!is.na(data[ , 3])) > 0 && max(cpp.agg[ , 4] * scale.fac, na.rm = TRUE) < 0.2)) {

		Log.lik = Log.lik - 100000

	}
	

	return(Log.lik) 
	# return(list('Log.lik' = Log.lik, 'cpp.out' = cpp.out))
	
}



Log.Lik.tau = function(taus, mean.tau, sd.tau) {
	
	Log.lik = 0
	
	Lik = log(dnorm(taus, mean = mean.tau, sd = sd.tau))
	Lik[which(!is.finite(Lik))] = -1000
	Log.lik = Log.lik + sum(Lik, na.rm = TRUE)
	
	return(Log.lik)
	
}



########
# Priors
########


Log.prior = function(family, value, params) {
	
	if (family == 'unif') {
		pr = dunif(value, min = params['min'], max = params['max'])
	} else if (family == 'normal') {
		pr = dnorm(value, mean = params['mean'], sd = params['sd'])
	} else if (family == 'poisson') {
		pr = dpois(round(value), lambda = params['lambda'])	
	} else if (family == 'exponential') {
		pr = dexp(value, rate = params['rate'])
	} else if (family == 'gamma') {
		pr = dgamma(value, shape = params['shape'], rate = params['rate'])
	} else {
		pr = NA
	}
	
	return(log(pr))
	
}


