

setwd('~/work/Infl_VN/MCMC_0501/Code/cpp_v4')

	
	if (exists('params')) { 
		all_params = params 
	} else {
		cat('\n\n The drawn parameters arent working! \n\n')
	}
		
		
	if ('nu_denom' %in% names(all_params)) {
		all_params['nu'] = 1 / all_params['nu_denom']
	}
	if ('rho_denom' %in% names(all_params)) {
		all_params['rho'] = 1 / all_params['rho_denom']
	}
		
	num.phis = length(grep('phi.', names(all_params)))
	for (i in grep('tau.', names(all_params))) {
		new.phi = all_params[paste('phi.', num.phis, sep = '')] + all_params[i]
		all_params[paste('phi.', num.phis+1, sep = '')] = new.phi
		num.phis = num.phis + 1
	}
	
	



sys.command = './odesim'
if ('beta_H1' %in% names(all_params)) {
	sys.command = paste(sys.command, '-beta1', all_params['beta_H1'])
}
if ('beta_B' %in% names(all_params)) {
	sys.command = paste(sys.command, '-beta2', all_params['beta_B'])
}
if ('beta_H3' %in% names(all_params)) {
	sys.command = paste(sys.command, '-beta3', all_params['beta_H3'])
}
if ('sigma12' %in% names(all_params)) {
	sys.command = paste(sys.command, '-sigma12', all_params['sigma12'])
}
if ('sigma13' %in% names(all_params)) {
	sys.command = paste(sys.command, '-sigma13', all_params['sigma13'])
}
if ('sigma23' %in% names(all_params)) {
	sys.command = paste(sys.command, '-sigma23', all_params['sigma23'])
}
if ('amp' %in% names(all_params)) {
	sys.command = paste(sys.command, '-amp', all_params['amp'])
}
if ('nu_denom' %in% names(all_params)) {
	sys.command = paste(sys.command, '-nu_denom', all_params['nu_denom'])
}
if ('rho_denom' %in% names(all_params)) {
	sys.command = paste(sys.command, '-rho_denom', all_params['rho_denom'])
}
if ('epidur' %in% names(all_params)) {
	sys.command = paste(sys.command, '-epidur', all_params['epidur'])
}
if (length(grep('phi.', names(all_params))) > 0) {
	sys.command = paste(sys.command, '-phi')
	for (i in grep('phi.', names(all_params))) { 
		sys.command = paste(sys.command, all_params[i])
	}
}
if ('zeroprevlevel' %in% names(all_params)) {
	sys.command = paste(sys.command, '-zeroprevlevel', all_params['zeroprevlevel'])
}
if ('immig' %in% names(all_params)) {
	sys.command = paste(sys.command, '-immig', all_params['immig'])
}

# cpp.run = system(command = paste(sys.command, '> throwaway.txt'))
# cpp.out = tryCatch(
# 	read.table('throwaway.txt', header = FALSE), 
# 	error = function(e) {
# 		cp.run.2 = system(command = paste(sys.command, '> throwaway.txt'))
# 		cpp.out = read.table('throwaway.txt', header = FALSE)
# 	}
# )

cpp.run = pipe(sys.command, 'rt')
cpp.out = read.delim(cpp.run, sep = '\t', header = FALSE)
close.connection(cpp.run)



# names(cpp.out) = c('t', 'trend', 'R1', 'R2', 'R3', 'I1', 'I2', 'I3', 'J1', 'J2', 'J3', 'S', 'N')
names(cpp.out)[ncol(cpp.out)-4] = 'J1'
names(cpp.out)[ncol(cpp.out)-3] = 'J2'
names(cpp.out)[ncol(cpp.out)-2] = 'J3'

cpp.out$inc1 = c(0, cpp.out$J1[-1] - cpp.out$J1[-nrow(cpp.out)])
cpp.out$inc2 = c(0, cpp.out$J2[-1] - cpp.out$J2[-nrow(cpp.out)])
cpp.out$inc3 = c(0, cpp.out$J3[-1] - cpp.out$J3[-nrow(cpp.out)])

# cpp.out = cpp.out[(nrow(cpp.out)-3649):nrow(cpp.out), ]

setwd('~/work/Infl_VN/MCMC_0501')



