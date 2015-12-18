# Simulation study when prior is the data generating process
# 68 species and 22 samples
# This script is designed to be called from command-line
# Instead of loading into R directly

source( "../MCMC_fns.R" )
library( snowfall )
library( getopt )
opt = getopt( c(
	"start.idx", "a", 1, "integer",
	"end.idx", "b", 1, "integer",
	"step", "s", 1, "integer",
	"cores", "c", 1, "integer",
	"truth.file", "t", 1, "character",
	"res.file", "r", 1, "character",
	"tru.factor", "n", 1, "integer",
	"tru.block", "m", 1, "integer"
	) )
	
sfInit( parallel=TRUE, cpus=opt$cores )

# Discretize prior distribution of sigma
# Based on beta distribution beta(alpha/p, 0.5-alpha/p)
# Asymptotically gives you Poisson process

sigma.value = seq(0.001,0.999,0.001)
tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

sfExportAll()

sfLapply( seq( opt$start.idx, opt$end.idx ), function(x){
  sim.res.shrink.test = sim.CDP.shrink.test( c(1000,10000,100000), 22, 68, opt$tru.factor, opt$tru.block, 3, 1, sigma.value, sigma.prior )
  saveRDS( sim.res.shrink.test, paste( opt$truth.file, x, sep = "" ) )
  # hyper$m is the maximal number of factors we consider
  hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = 10 )
  for( id in 1:length( sim.res.shrink.test$data ) ){
	dat.ind = sim.res.shrink.test$data[[id]]
	start = list( sigma = sample( sigma.value, size = length( sim.res.shrink.test$sigma ), replace = T, prob = sigma.prior ), 
                T.aug = colSums( dat.ind ),
                Q = matrix( 0.5, nrow = 68, ncol = 22 ),
                X = matrix( rnorm( hyper$m*68 ), nrow = hyper$m ),
                Y = matrix( rnorm( hyper$m*22 ), nrow = hyper$m ),
                er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                delta = c( rgamma( 1, shape = 3, rate = 1 ), rgamma( hyper$m-1, shape = 4, rate = 1 ) ),
                phi = matrix( rgamma( hyper$m*22, shape = 3/2, rate = 3/2 ), nrow = 22 ) )
	main.mcmc.shrink( dat.ind, start, 
                    hyper, sigma.value, sigma.prior, paste( opt$res.file, x, "/sim", "_", id, sep = "" ), save.obj = c("Y","er"), step = opt$step )
  }
})

sfStop()
