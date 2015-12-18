# Simulation study to check the precision of estimates of microbial distributions
# Output is stored in RDS file
# Output is the total variation distances

library( snowfall )
library( optparse )
source( "../MCMC_fns.R" )

option_list = list(
	make_option( c("-n", "--n.pop"), action = "store", default = 22, type = "integer", help = "No. of populations" ),
	make_option( c("-p", "--n.species"), action = "store", default = 68, type = "integer", help = "No. of species" ),
	make_option( c("-f", "--n.tru.factor"), action = "store", default = 6, type = "integer", help = "No. of true factors" ),
	make_option( c("-b", "--n.block"), action = "store", default = 2, type = "integer", help = "No. of population blocks" ),
	make_option( c("-s", "--a.er"), action = "store", default = 3, type = "double", help = "First param of Gamma prior for er" ),
	make_option( c("-r", "--b.er"), action = "store", default = 1, type = "double", help = "Second param of Gamma prior for er" ),
	make_option( c("-d", "--Y.sd"), action = "store", default = 0.5, type = "double", help = "Similarity between populations" ),
	make_option( c("-c", "--alpha"), action = "store", default = 10, type = "double", help = "Concentration param of DP" ),
	make_option( c("-q", "--Q.power"), action = "store", default = 2, type = "double", help = "Power of Q component in data generating process (Q\neq 2->misspecied)" ),
	make_option( c("-a", "--rep.tot"), action = "store", default = 1, type = "integer", help = "number of replication" ),
	make_option( c("-l", "--read.max"), action = "store", default = 100, type = "integer", help = "max number of read depth" ),
	make_option( c("-e", "--read.inc"), action = "store", default = 10, type = "integer", help = "increment in read depth" ),
	make_option( "--step", action = "store", default = 30000, type = "integer", help = "MCMC steps" ),
	make_option( "--thin", action = "store", default = 15, type = "integer", help = "MCMC thinning" ),
	make_option( "--save.path", action = "store", default = NA, type = "character", help = "Cache path" ),
	make_option( "--res.path", action = "store", default = NA, type = "character", help = "Result path" )
)
opt = parse_args(OptionParser(option_list=option_list))

#sigma follows beta(alpha/p,1/2)
sigma.value = seq(0.001,0.999,0.001)
tmp = c(0,pbeta( sigma.value, opt$alpha/opt$n.species, 1/2-opt$alpha/opt$n.species ))
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

# new data simulation module to simulate biological samples with compound symmetric covariance sturcture
sim.CDP.shrink.test = function( n, p, m, K, a.er, b.er, sigma.value, sigma.weight, Y.sd = 0.5 ){
  require( mvtnorm )
  #we consider K block simulation stucture
  sigma = sample( sigma.value, p, replace = T, sigma.weight )
  X = matrix( rnorm( p*m ), nrow = m )
  #split the m factors into K blocks
  factor.indx.group = split( 1:m, cut( 1:m, breaks = K ) )
  pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) )
  Y = matrix( 0, nrow = m, ncol = n )
  
  for( x in 1:K ){
	cor.mt = diag( 1, length(pop.indx.group[[x]]) )
    cor.mt[upper.tri(cor.mt)] = Y.sd
    cor.mt[lower.tri(cor.mt)] = Y.sd 
    raw = rmvnorm( length(factor.indx.group[[x]]), mean=rep(1,length(pop.indx.group[[x]])), sigma=cor.mt )
    Y[factor.indx.group[[x]],pop.indx.group[[x]]] = raw
  }
  
  er = 1/rgamma( 1, a.er, b.er)
  Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
  
  return( list( sigma = sigma, Q = Q,
                X.tru=X, Y.tru = Y, er = er ) )
}

prediction.pw.ind = function( sim.res, read.depth, save.path, Q.power ){
	res.cache = c()
	final.weights = sim.res$sigma*(sim.res$Q*(sim.res$Q>0))^Q.power
	
	for( i in 1:length( read.depth ) ){
	cat( sprintf( "Read depth is %d\n", read.depth[i] ) )
    #simulate data
    dat.ind = apply( final.weights, 2, function(x) rmultinom( 1, read.depth[i], prob=x ) )
    #Bayesian method
    hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = 10 )
    start = list( sigma = sample( sigma.value, size = 68, replace = T, prob = sigma.prior ),
                  T.aug = colSums( dat.ind ),
                  Q = matrix( 0.5, nrow = 68, ncol = 22 ),
                  X = matrix( rnorm( 68*hyper$m ), nrow = hyper$m ),
                  Y = matrix( rnorm( 22*hyper$m ), nrow = hyper$m ),
                  er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                  delta = c( rgamma( 1, shape = 3, rate = 1 ), rgamma( hyper$m - 1, shape = 4, rate = 1 ) ),
                  phi = matrix( rgamma( 22*hyper$m, shape = 3/2, rate = 3/2 ), nrow = 22 ) )
    main.mcmc.shrink( dat.ind, start, hyper, sigma.value, sigma.prior, paste( save.path, "res", sep = "/"), save.obj = c("sigma","Q"), step = opt$step, thin=opt$thin )
    res = lapply( list.files( save.path, pattern = "res", full.names=T ), readRDS )
    bayes.est.ls = lapply( res, function(x){
      weight = x$sigma*(x$Q)^2*(x$Q>0)
      t(t(weight)/colSums(weight) )
    })
	bayes.est.mt = array( unlist( bayes.est.ls ), dim = c( dim( bayes.est.ls[[1]] ), length( bayes.est.ls ) ) )
	bayes.est.mean = apply( bayes.est.mt, 1:2, mean )
	bayes.diff = bayes.est.mean - t(t(final.weights)/colSums(final.weights))
	bayes.error = apply( bayes.diff, 2, function(x) sum(x[x>0]) )
    
    #direct normalization method
    dat.diff = t(t(dat.ind)/colSums(dat.ind)) - t(t(final.weights)/colSums(final.weights))
    dat.error = apply( dat.diff, 2, function(x) sum(x[x>0]) )
    
    #gain in each population for bayesian method
    res.cache = rbind( res.cache, dat.error - bayes.error )
  }
  res.cache
}

read.depth = seq( opt$read.inc, opt$read.max, opt$read.inc )

sfInit( parallel=TRUE, cpus=opt$rep.tot )
sfExportAll()

res.all = sfLapply( 1:opt$rep.tot, function(rep){
	cat( sprintf( "Replication %d\n", rep ) )
	save.path = paste( opt$save.path, rep, sep = "_" )
	
	sim.res = sim.CDP.shrink.test( opt$n.pop, opt$n.species, opt$n.tru.factor, opt$n.block, opt$a.er, opt$b.er, sigma.value, sigma.prior, opt$Y.sd )
	
	prediction.pw.ind( sim.res, read.depth, save.path, opt$Q.power )
} )
sfStop()

saveRDS( res.all, opt$res.path )
