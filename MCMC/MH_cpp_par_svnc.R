
date()

rm(list = ls())


library(tictoc)
library(mvtnorm)
library(parallel)


file.num = 26
	type = 'svnc'

source('~/work/Infl_VN/MCMC_0501/Code/MCMC/Likelihoods_priors_cpp.R')


# make.need = TRUE
if (exists('make.need')) {
	setwd('~/work/Infl_VN/MCMC_0501/Code/cpp_v4')
	system('make')
}


tic()

###
### 1. Read in data
###

load('~/work/Infl_VN/Data/SVN_Detr_Ag.Rdata')
dat = dat.south

dat$Time_cts = dat$Year + ((dat$Week-1)/52)
dat = dat[order(dat$Time_cts), ]
dat$ILI_p = dat$ILI_p_H1 + dat$ILI_p_B + dat$ILI_p_H3
dat$ILI_p_A = dat$ILI_p_H1 + dat$ILI_p_H3





###
### 	2. Define all of the parameters you are going to need in your model; 
###    some will be fixed; some will be fit by the MCMC
###			Serves as initial values
###

all_params = c(
	'N' = 1e6, 
	'Na' = 0.8 * 1e6, 
	'Nb' = 0.2 * 1e6, 
	'beta_H1' = 0.1950, # runif(1, 0.19, 0.23), # 0.26, 
	'beta_B' = 0, # 0.23, 
	'beta_H3' = 0, 
	'nu_denom' = 5, 
	'amp' = 0.35, 
	'rho_denom' = 1616, #sample(c(400:1200), 1), 
	'sigma12' = 0.75, 
	'sigma13' = 0.50, 
	'sigma23' = 0.75, 
	'etaab' = 0.01, 
	'etaba' = 0.05, 
	'phi.1' = 115, #sample(c(1:400), 1), # 250, 
	'tau.1' = 165, #sample(c(100:650), 1), # 500, 
	'tau.2' = 220, #sample(c(100:650), 1), # 500, 
	'tau.3' = 404, #sample(c(100:650), 1), # 300, 
	'tau.4' = 286, #sample(c(100:650), 1), # 400, 
	'tau.5' = 120, #sample(c(100:650), 1), # 200, 
	'tau.6' = 330, #sample(c(100:650), 1), # 400, 
	'tau.7' = 416, #sample(c(100:650), 1), # 400, 
	'tau.8' = 250, #sample(c(100:650), 1), # 350, 
	'tau.9' = 132, #sample(c(100:650), 1), # 350, 
	'tau.10' = 222, #sample(c(100:650), 1), # 350, 
	'tau.11' = 168, #sample(c(100:650), 1), # 350, 
	'tau.12' = 166, #sample(c(100:650), 1), # 350, 
	'tau.13' = 340, #sample(c(100:650), 1), # 350, 
	'mu.tau' = 465, 
	'sigma.tau' = 100, 
	'report' = 0.096, 
	'norm_var' = 0.1^2, 
	'epidur' = 25, #sample(c(30:150), 1), 
	'zeroprevlevel' = 6, #sample(c(0:500), 1),
	'immig' = 302 #sample(seq(10, 500, by = 10), 1)
	)




if (file.num > 1) {
	
	load(paste('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_', file.num-1, '.Rdata', sep = ''))
	
	start.pars = list('p1' = pps[[1]][which(pps[[1]][ , 'it'] == max(pps[[1]][ , 'it']))[1], names(all_params)], 
		'p2' = pps[[2]][which(pps[[2]][ , 'it'] == max(pps[[2]][ , 'it']))[1], names(all_params)], 
		'p3' = pps[[3]][which(pps[[3]][ , 'it'] == max(pps[[3]][ , 'it']))[1], names(all_params)], 
		'p4' = pps[[4]][which(pps[[4]][ , 'it'] == max(pps[[4]][ , 'it']))[1], names(all_params)]) 
	
	rm(pps)

} else if (file.num == 1) {

	load(paste('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_', file.num-1, '.Rdata', sep = ''))
	
	start.pars = list('p1' = pps[[1]][which(pps[[1]][ , 'post'] == max(pps[[1]][ , 'post']))[1], names(all_params)], 
		'p2' = pps[[2]][which(pps[[2]][ , 'post'] == max(pps[[2]][ , 'post']))[1], names(all_params)], 
		'p3' = pps[[3]][which(pps[[3]][ , 'post'] == max(pps[[3]][ , 'post']))[1], names(all_params)], 
		'p4' = pps[[4]][which(pps[[4]][ , 'post'] == max(pps[[4]][ , 'post']))[1], names(all_params)]) 
	
	rm(pps)	

} else if (file.num == 0) {

	# load('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0509.Rdata')
	load('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_0.Rdata')

	start.pars = list()
	for (i in 1:4) {
		ss = sample(c(1:4), 1); ss = i
		start.pars[[i]] = pps[[ss]][which(pps[[ss]][ , 'post'] == max(pps[[ss]][ , 'post']))[1], names(all_params)]

	}

	rm(pps)

}




###
### 	3. Initialize the left/right bounds and sd values
###			for the prior distributions on all the parameters
###		Sd values govern MCMC movement sizes
###

all_params_bounds = c(
	1e6 - 1, 1e6 + 1, 1, # N  
	(0.8 * 1e6) - 1, (0.8 * 1e6) + 1, 1, #Na 
	(0.2 * 1e6) - 1, (0.2 * 1e6) + 1, 1, #Nb
	0.14, 1.0, 0.01, # beta_H1 
	0.14, 1.0, 0.01, # beta_B 
	0.14, 1.0, 0.01, # beta_H3 
	3, 8, 0.5, # nu_denom 
	0.01, 1.5, 0.01, # amp 
	365, 365*6, 10, # rho_denom 
	0.05, 0.95, 0.1, # sigma12 
	0.05, 0.95, 0.1, # sigma13 
	0.05, 0.95, 0.1, # sigma23 
	0, 0.8, 0.05, # etaab 
	0, 0.8, 0.05, # etaba 
	1, 450, 5, # phi.1
	20, 700, 5, # tau.1
	20, 700, 5, # tau.2
	20, 700, 5, # tau.3
	20, 700, 5, # tau.4
	20, 700, 5, # tau.5
	20, 700, 5, # tau.6
	20, 700, 5, # tau.7
	20, 700, 5, # tau.8
	20, 700, 5, # tau.9
	20, 700, 5, # tau.10
	20, 700, 5, # tau.11
	20, 700, 5, # tau.12
	20, 700, 5, # tau.13
	20, 700, 5, # mu.tau
	1, 200, 2,  # sigma.tau
	0.0001, 0.1, 0.0002,  # reporting param
	0.001, 1, 0.1,  # norm_var
	20, 180, 2, # epidur
	1, 1000, 5, # zeroprevlevel
	30, 1000, 5 # imig
	)


if (length(all_params_bounds) != 3*length(all_params)) {
	cat('check your parameters and bounds!!!', '\n\n')
}


all_params_leftbound = all_params_rightbound = all_params_sigma = list()

all_params_leftbound = all_params_bounds[seq(1, length(all_params_bounds), by = 3)]
all_params_rightbound = all_params_bounds[seq(2, length(all_params_bounds), by = 3)]
all_params_sigma = all_params_bounds[seq(3, length(all_params_bounds), by = 3)]




###
### 	4. Define the search parameters' priors
###



prior_family = c(
	'N' = NA, 
	'Na' = NA, 
	'Nb' = NA, 
	'beta_H1' = 'unif', 
	'beta_B' = 'unif', 
	'beta_H3' = 'unif', 
	'nu_denom' = 'unif', 
	'amp' = 'exponential', 
	'rho_denom' = 'unif', 
	'sigma12' = 'unif', 
	'sigma13' = 'unif', 
	'sigma23' = 'unif', 
	'etaab' = 'unif', 
	'etaba' = 'unif', 
	'phi.1' = 'unif', 
	'tau.1' = 'unif', 
	'tau.2' = 'unif', 
	'tau.3' = 'unif', 
	'tau.4' = 'unif', 
	'tau.5' = 'unif', 
	'tau.6' = 'unif', 
	'tau.7' = 'unif', 
	'tau.8' = 'unif', 
	'tau.9' = 'unif', 
	'tau.10' = 'unif', 
	'tau.11' = 'unif', 
	'tau.12' = 'unif', 
	'tau.13' = 'unif', 
	'mu.tau' = 'unif', 
	'sigma.tau' = 'unif', 
	'report' = 'unif', 
	'norm_var' = 'unif', 
	'epidur' = 'unif', 
	'zeroprevlevel' = 'unif', 
	'immig' = 'unif'
	)
	
	
	
prior_params = list(
	'N' = NA, 
	'Na' = NA, 
	'Nb' = NA, 
	'beta_H1' = c('min' = 0.14, 'max' = 1), 
	'beta_B' = c('min' = 0.14, 'max' = 1), 
	'beta_H3' = c('min' = 0.14, 'max' = 1), 
	'nu_denom' = c('min' = 3, 'max' = 8), 
	'amp' = c('rate' = 1), 
	'rho_denom' = c('min' = 180, 'max' = 365*6), 
	'sigma12' = c('min' = 0, 'max' = 1), 
	'sigma13' = c('min' = 0, 'max' = 1), 
	'sigma23' = c('min' = 0, 'max' = 1), 
	'etaab' = c('min' = 0, 'max' = 1),  
	'etaba' = c('min' = 0, 'max' = 1), 
	'phi.1' = c('min' = 0, 'max' = 450),  
	'tau.1' = c('min' = 20, 'max' = 700),  
	'tau.2' = c('min' = 20, 'max' = 700),  
	'tau.3' = c('min' = 20, 'max' = 700),  
	'tau.4' = c('min' = 20, 'max' = 700),  
	'tau.5' = c('min' = 20, 'max' = 700),  
	'tau.6' = c('min' = 20, 'max' = 700),  
	'tau.7' = c('min' = 20, 'max' = 700),  
	'tau.8' = c('min' = 20, 'max' = 700),  
	'tau.9' = c('min' = 20, 'max' = 700),  
	'tau.10' = c('min' = 20, 'max' = 700),  
	'tau.11' = c('min' = 20, 'max' = 700),  
	'tau.12' = c('min' = 20, 'max' = 700),  
	'tau.13' = c('min' = 20, 'max' = 700),  
	'mu.tau' = c('min' = 20, 'max' = 700), 
	'sigma.tau' = c('min' = 1, 'max' = 200), 
	'report' = c('min' = 0.0001, 'max' = 0.1), 
	'norm_var' = c('min' = 0.1, 'max' = 2), 
	'epidur' = c('min' = 20, 'max' = 180), 
	'zeroprevlevel' = c('min' = 0, 'max' = 1000), 
	'immig' = c('min' = 30, 'max' = 1000)
	)



if ((prod(names(all_params) %in% names(prior_family)) != 1) || 
	(prod(names(all_params) %in% names(prior_params)) != 1)) {
		
	cat('\n check your priors and parameters \n\n')
}





###
### 5. Set items for MCMC
###


file.name = paste('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_', file.num, '.Rdata', sep = '')

search_params_names = c('beta_H1', 'amp', 'rho_denom', 
   	'phi.1', 'tau.1', 'tau.2', 'tau.3', 'tau.4', 'tau.5', 'tau.6', 'tau.7', 'tau.8', 'tau.9',
	'tau.10', 'tau.11', 'tau.12', 'tau.13',  
 	'mu.tau', 'sigma.tau', 'report', 'epidur', 'zeroprevlevel', 'immig')
search_params_index = which(names(all_params) %in% search_params_names)


source('~/work/Infl_VN/MCMC_0501/Code/MCMC/MCMC_Code.R')

mv.pars = list(c('beta_H1'), c('amp'), c('rho_denom'), c('phi.1'), 
	c('tau.1', 'tau.2', 'tau.3', 'tau.4', 'tau.5', 'tau.6', 'tau.7','tau.8', 'tau.9', 
	'tau.10', 'tau.11', 'tau.12', 'tau.13'), 
	c('mu.tau'), c('sigma.tau'), c('report'), c('epidur'), c('zeroprevlevel'), c('immig'))




###
### 6. Begin Metropolis-Hastings algorithm
###


# if (3 == 4) {

pp = MH_chain(ch.num = 0, chain.length = 200, 
 	burn.in = 0,
	lik.fam = 'normal',  
	recheck.cov = 500, 
	search_params_names = search_params_names, 
	search_params_index = search_params_index, 
	cur.params = start.pars[[1]], 
	all_params_sigma = all_params_sigma, 
	all_params_bounds = list(all_params_leftbound, all_params_rightbound), 
	priors_fams = prior_family, 
	priors_pars = prior_params, 
	mv.pars = mv.pars, 
	use.cur = use.cur,
	type = 'svnc', 
	dat = as.matrix(cbind(dat$ILI_p, rep(NA, nrow(dat)), rep(NA, nrow(dat)))))

# }
	
	
toc()




# if (3 == 4) {
tic()
pps = list()

pps = mclapply(1:num.chains, function(x) {
	
	pp = MH_chain(ch.num = x, chain.length = chain.length, 
		burn.in = 0, 
		lik.fam = lik.fams[[x]], 
		recheck.cov = recheck.cov[x], 
		search_params_names = search_params_names, 
		search_params_index = search_params_index, 
		cur.params = start.pars[[x]], 
		all_params_sigma = all_params_sigma, 
		all_params_bounds = list(all_params_leftbound, all_params_rightbound), 
		priors_fams = prior_family, 
		priors_pars = prior_params, 
		mv.pars = mv.pars, 
		use.cur = use.cur,
		type = 'svnc',  
		dat = as.matrix(cbind(dat$ILI_p, rep(NA, nrow(dat)), rep(NA, nrow(dat)))))
	
	}, mc.cores = num.chains, mc.preschedule = TRUE)

toc()

save(pps, file = file.name)



if (file.num >= 3) { 

	rm(pps)

	all.data = list()
	load('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_1.Rdata')
	for (i in 1:length(pps)) { all.data[[i]] = pps[[i]] } 

	rm(pps)

	for (i in 2:file.num) {

		load(paste('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501_', i, '.Rdata', sep = ''))
		
		for (j in 1:length(all.data)) { all.data[[j]] = rbind(all.data[[j]], pps[[j]]) }

		rm(pps)

	}

	for (i in 1:length(all.data)) {
		all.data[[i]] = all.data[[i]][which(all.data[[i]][ , 'it'] > 0), ]
	}

	save(all.data, file = paste('~/work/Infl_VN/MCMC_0501/MCMC_svnc/MCMC_svnc_0501.Rdata', sep = ''))

	load(file.name)

}



for (i in 1:length(pps)) { cat(i, '\t', dim(pps[[i]]), '\n') }
for (i in 1:length(pps)) { cat(type, '	', file.num, '	', i, '	', dim(pps[[i]]), '
', file = 'outputdim.txt', append = TRUE) }


date()

times = NULL
for (i in 1:length(pps)) { times = c(times, pps[[i]][ , 'time']) }
summary(times)

type = 'svnc'; source('~/work/Infl_VN/MCMC_0501/Code/plotting.R')


# } # goes with 3 == 4

date()


