# script used to test the self-consistent algorihtm and generate Figure S2.
# the truth file should contain the values of parameters used to generate the data

source("fast_fns.R")
source("../utilities.R")
library( "MCMCpack" )
library( "tmvtnorm" )

all.est = list(50)
all.tru = list(50)
for( indx in 1:50 ){
  print( indx )
  sim = readRDS(paste("path_to_truth_file_",indx,sep=""))
  data.use = sim$data[rowSums(sim$data)>0,]
  data.use.norm = t(t(data.use)/colSums(data.use))
  sigma.est = rowMeans( apply( data.use.norm, 2, function(x) x/(1-x)*(nrow(data.use.norm)-1)  ) )
  Q.est = sqrt(data.use.norm/sigma.est)
  corr.tru = cov2cor( diag( rep( sim$er, ncol(sim$Y.tru) ) ) + t(sim$Y.tru)%*%sim$Y.tru )
  all.tru[[indx]] = corr.tru
  
  # estimate the correlation matrix in a pairwise manner to reduce computation burden
  mt.sim = diag( rep( 1, ncol( sim$Q ) ) )
  for( i in 1:(ncol( sim$Q )-1) ){
    for( j in (i+1):ncol( sim$Q ) ){
      print( c(i,j) )
      X.full.zero = Q.est[,c(i,j)]
      tryCatch({
        corr.est = cov2cor( get_corr_est(X.full.zero,max.iter = 100) )
        mt.sim[i,j] = corr.est[1,2]
        mt.sim[j,i] = corr.est[1,2]},
        error = function(e){cat(c(i,j,indx))})
    }
  }
  all.est[[indx]] = mt.sim
}

saveRDS( all.est, "fast_all_est.rds")
saveRDS( all.tru, "fast_all_truth.rds")

all.est = readRDS( "fast_all_est.rds" )
all.tru = readRDS( "fast_all_truth.rds" )

a = sapply( 1:50, function(x) rv.coef( all.est[[x]], all.tru[[x]] ) )

mt.sim[lower.tri(mt.sim)] = corr.tru[lower.tri(corr.tru)]
gg.heatmap( mt.sim, "Truth", "Estimated" )
