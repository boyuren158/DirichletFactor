# script used to generate MCMC results for Ravel's vaginal microbiome data

source("../MCMC_fns.R")
# "ravel.rds" is a data matrix with species in row and biological samples in column
dat.use.filter = readRDS( 'ravel.rds' )

#discretize prior of sigma
alpha = 10
sigma.value = seq(0.001,0.999,0.001)
tmp = c( 0,pbeta( sigma.value, alpha/p, 1/2-alpha/p ) )
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
plot( sigma.value, sigma.prior )

hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = 10 )
get.start = function( data, sigma.value, sigma.prior, hyper ){
  list( sigma = sample( sigma.value, nrow( data ), replace = T),
        T.aug = rep( 1000, ncol( data ) ),
        Q = matrix( 1, nrow = nrow( data ), ncol = ncol( data ) ),
        X = matrix( rnorm( hyper$m*nrow( data ) ), nrow = hyper$m ),
        Y = matrix( rnorm( hyper$m*ncol( data ) ), nrow = hyper$m ),
        er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
        delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ),
                   rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
        phi = matrix( rgamma( hyper$m*ncol( data ),
                              shape = hyper$nv/2, rate = hyper$nv/2 ),
                      nrow = ncol( data ) ) )
}

start.ravel = list( sigma = sample( sigma.value, nrow(  data.use.filter ), replace = T),
                   T.aug = rep( 1000, ncol(  data.use.filter) ),
                   Q = matrix( 1, nrow = nrow(  data.use.filter), ncol = ncol(  data.use.filter) ),
                   X = matrix( rnorm( hyper$m*nrow(  data.use.filter ) ), nrow = hyper$m ),
                   Y = matrix( rnorm( hyper$m*ncol(  data.use.filter ) ), nrow = hyper$m ),
                   er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                   delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ),
                              rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
                   phi = matrix( rgamma( hyper$m*ncol(  data.use.filter ),
                                         shape = hyper$nv/2, rate = hyper$nv/2 ),
                                 nrow = ncol(  data.use.filter ) ) )

source( "myMCMC_fn_vec.R" )
main.mcmc.shrink( as.matrix( dat.use.filter ), start.ravel, hyper, sigma.value, sigma.prior, save_path = "ravel_new/res", step = 50000, thin = 10 )
