

###
### 5. Set items for MCMC
###



num.chains = 4

lik.fams = as.list(rep('normal', num.chains))

chain.length = 750000
if (file.num == 0) { chain.length = 150000 }
burn_in = 0

if (file.num < 3) { 
	recheck.cov = rep(chain.length * 2, num.chains) 
	use.cur = FALSE
} else if (file.num > 14) {
	recheck.cov = rep(chain.length * 2, num.chains)
	use.cur = TRUE
} else if (file.num > 8) {
	recheck.cov = rep(100000, num.chains)
	use.cur = TRUE
} else {
	recheck.cov = rep(50000, num.chains)
	use.cur = TRUE
}

# cat('\n\n\n', date(), '\n\n\n', file = '~/work/Infl_VN/MCMC_0501/oohgurl.txt')

###
### 6. Begin Metropolis-Hastings algorithm
###

MH_chain = function(ch.num, chain.length, lik.fam, burn.in, recheck.cov, 
	search_params_names, search_params_index, 
	cur.params, all_params_sigma, all_params_bounds, 
	priors_fams, priors_pars, 
	mv.pars, use.cur, type,  
	dat) {
		
	s2.shaby = (2.4^2) / length(search_params_index)
	
	for(it in 1:chain.length) {

		start.time = proc.time()
		
		# Establish first parameter/posterior value
		
		if (it == 1) {
		
			if (paste('cur_cov', type, ch.num, '_', file.num-1,  '.csv', sep = '') %in% 
				list.files(path = paste('~/work/Infl_VN/MCMC_0501/MCMC_', type, '/cov', sep = '')) && 
				use.cur == TRUE) {
				cov = read.csv(paste('~/work/Infl_VN/MCMC_0501/MCMC_', type, '/cov/cur_cov', type,  
					ch.num, '_', file.num-1, '.csv', sep = ''))
				rownames(cov) = cov$X
				cov$X = NULL; cov$X.1 = NULL
				cov = as.matrix(cov)
			} else {
				cov = matrix(0, nrow = length(search_params_index), ncol = length(search_params_index))
				diag(cov) = (all_params_sigma[search_params_index])^2
				rownames(cov) = colnames(cov) = search_params_names
			}
			

			R.denom = Log.Lik(family = lik.fam, data = dat, params = cur.params, it = it)
			if (prod(c('mu.tau', 'sigma.tau') %in% search_params_names) == 1) {
				R.denom = R.denom + Log.Lik.tau(taus = cur.params[grep('tau.', names(cur.params))], 
					mean.tau = cur.params['mu.tau'], sd.tau = cur.params['sigma.tau'])
			}	
			
			for (j in 1:length(search_params_index)) {
				R.denom = R.denom + 
					Log.prior(family = priors_fams[[search_params_index[j]]], 
					params = priors_pars[[search_params_index[j]]], 
					value = cur.params[search_params_index[j]])
			}
			
			post = rbind(c(0, cur.params, R.denom, NA), c(0, cur.params, R.denom, NA))
			colnames(post) = c('it', names(cur.params), 'post', 'time')
			
		}
		
		# Propose new parameters, make sure they are within bounds
		
		new.params = cur.params
		new.params[search_params_index] = round(rmvnorm(1, mean = cur.params[search_params_index], 
			sigma = cov), 8)
			
		while(prod(new.params[search_params_index] > all_params_bounds[[1]][search_params_index]) != 1 ||
			prod(new.params[search_params_index] < all_params_bounds[[2]][search_params_index]) != 1 || 
			sum(new.params[grep('tau.', names(new.params))]) > 3650 ||
			sum(new.params[grep('tau.', names(new.params))] <= new.params['epidur']) > 0) {
		
			# if (ch.num == 5) {
			# cat(it, '\t', new.params[search_params_index], '\n', 
			# 	'\t', new.params[search_params_index] > all_params_bounds[[1]][search_params_index], '\n', 
			# 	'\t', new.params[search_params_index] < all_params_bounds[[2]][search_params_index], '\n', 
			# 	'\t', sum(new.params[grep('tau.', names(new.params))]) < 3650, '\n',  
			# 	'\t', new.params[grep('tau.', names(new.params))] > new.params['epidur'], '\n\n', 
			# 	file = '~/work/Infl_VN/MCMC_0501/oohgurl.txt', append = TRUE) 
			# }		

			new.params[search_params_index] = round(rmvnorm(1, mean = cur.params[search_params_index], 
				sigma = cov), 8)

		}	
		
		
		# Evaluate posterior for new parameters
		
		R.num = Log.Lik(family = lik.fam, data = dat, params = new.params, it = it)
		
		if (prod(c('mu.tau', 'sigma.tau') %in% search_params_names) == 1) {
			R.num = R.num + Log.Lik.tau(taus = new.params[grep('tau.', names(new.params))], 
				mean.tau = new.params['mu.tau'], sd.tau = new.params['sigma.tau'])
		}
		for (j in 1:length(search_params_index)) {
			R.num = R.num + 
				Log.prior(family = priors_fams[[search_params_index[j]]], 
				params = priors_pars[[search_params_index[j]]], 
				value = new.params[search_params_index[j]])
		}
		
		
		
		# Compare posteriors, accept/reject
		
		R = R.num - R.denom
		
		prob = runif(1, 0, 1)
		if (exp(R) > prob) {
			post = rbind(post, c(it, new.params, R.num, (proc.time()-start.time)[3]))
			cur.params = new.params
			R.denom = R.num
		}
		
		# Update covariance matrix if needed
		
		if (it %% recheck.cov == 0 && length(mv.pars) > 0) {
			
			c0 = 1; c1 = 0.8 # c0 > 0 and c1 in (0, 1] per Shaby and Wells
			
			zz = post[which(post[ , 1] > (it - recheck.cov) & post[ , 1] <= it), ]
			if (!is.null(nrow(zz)) && nrow(zz) > 0.01*recheck.cov && nrow(zz) < 0.9*recheck.cov) {
				
				for (ss in 1:length(mv.pars)) {
					
					mv.index = which(names(cur.params) %in% mv.pars[[ss]])
					
					r.hat = nrow(zz) / recheck.cov
					Sigma.hat = var(zz[ , mv.index+1], na.rm = TRUE)
					
					g1 = (recheck.cov/it)^c1
					g2 = c0 * g1
					
					s2.shaby = exp(log(s2.shaby) + (g2 * (r.hat - 0.234)))

					cov.new = cov[which(rownames(cov) %in% mv.pars[[ss]]), 
						which(colnames(cov) %in% mv.pars[[ss]])] + 
						(g1 * (Sigma.hat - cov[which(rownames(cov) %in% mv.pars[[ss]]), 
						which(colnames(cov) %in% mv.pars[[ss]])]))
					
					cov[which(rownames(cov) %in% mv.pars[[ss]]), 
						which(colnames(cov) %in% mv.pars[[ss]])] = cov.new
					
				}
				
				diag(cov)[which(diag(cov) == 0)] = 
					all_params_sigma[names(diag(cov))[which(diag(cov) == 0)]]
				
			} else if (is.null(nrow(zz)) || (!is.null(nrow(zz)) && nrow(zz) < 0.01*recheck.cov)) { 
				cov = cov * 0.9
			} else if (!is.null(nrow(zz)) && nrow(zz) > 0.9 * recheck.cov) {
				cov = cov * 1.1
			}
		}

		if (it %in% round(seq(chain.length/10, chain.length*9/10, len = 9))) {
			save(post, file = paste('~/work/Infl_VN/MCMC_0501/MCMC_', type, '/midchain', type, ch.num, '.Rdata', sep = ''))
		}

	}	# goes with for (it in 1:chain.length)
	
	write.csv(cov, paste('~/work/Infl_VN/MCMC_0501/MCMC_', type, '/cov/cur_cov', type, ch.num, '_', file.num, '.csv', sep = ''))
	
	
	# return(rbind(post0, post))
	return(post) # return(list('post' = post, 'cov' = cov)) # return(post)
		

} # closes function




