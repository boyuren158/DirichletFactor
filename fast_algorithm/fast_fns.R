# modules used to carry out self consistent estimate of S in Sec 3.1

library( "MCMCpack" )
library( "tmvtnorm" )
library( "Matrix" )

# functions used to generate covariance matrix of flattened Q
get.expand.Sigma = function( Sigma, data.ratio ){
  n.species = nrow(data.ratio)
  value = rep( Sigma, each = n.species )
  ij.ind = expand.grid(  1:nrow(Sigma), 1:ncol(Sigma) )
  ij.ind.expand = Reduce( rbind, lapply( 1:nrow(ij.ind), function(x){
    cbind( (ij.ind[x,1]-1)*n.species + (1:n.species), 
           (ij.ind[x,2]-1)*n.species + (1:n.species) )
  }) )
  sparseMatrix( i = ij.ind.expand[,1], j = ij.ind.expand[,2], x = value ) 
}

# functions used to carry out the EM-algorithm in Sec 3.1
get_corr_est = function(X.full.zero, max.iter = 50, error = 1e-3){
  X.neg.ind = (1:length(X.full.zero))[c(X.full.zero)==0]
  X.pos.ind = (1:length(X.full.zero))[c(X.full.zero)>0]
  
  Sigma.start = diag( apply( X.full.zero, 2, function(x) sd(x[x>0]) ) )
  
  for( i in 1:max.iter ){
    Sigma.expand = get.expand.Sigma( Sigma.start, X.full.zero )
    
    S11 = Sigma.expand[X.neg.ind,X.neg.ind]
    S12 = Sigma.expand[X.neg.ind,X.pos.ind]
    S22 = Sigma.expand[X.pos.ind,X.pos.ind]
    
    Sigma.cond = S11 - S12%*%solve(S22)%*%t(S12)
    mu.cond = S12%*%solve(S22)%*%X.full.zero[X.pos.ind]
    
    X.neg.all = rtmvnorm( 1000, mean = as.vector(mu.cond), 
                          sigma = Sigma.cond, 
                          upper = rep(0,length(mu.cond)), 
                          algorithm = "gibbs" )
    X.full.ls = lapply( 1:nrow(X.neg.all), 
                        function(x){
                          X.neg = X.neg.all[x,]
                          X.full = X.full.zero
                          X.full[X.neg.ind] = X.neg
                          X.full
                        } )
	# replace t(x)%*%x by cov(x) is in principal the same, but has better finite sample performance
    Sigma.est = Reduce( "+", lapply( X.full.ls, function(x) cov(x) ) )/length( X.full.ls )
    
    #print( max(abs(corr.est-corr.start)) )
    if( max(abs(cov2cor(Sigma.est)-cov2cor(Sigma.start)))<=1e-3 )
      break
    Sigma.start = Sigma.est
  }
  Sigma.est
}
