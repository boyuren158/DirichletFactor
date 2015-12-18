# script used to generate Figure 2(f)

source("../MCMC_fns.R")
source("../utilities.R")

sim.for.contour = function( lscounts, n, p, m, K, a.er, b.er, sigma.value, sigma.weight, strength ){
  #we consider K block simulation stucture
  #browser()
  sigma = sample( sigma.value, p, replace = T, sigma.weight )
  X = matrix( rnorm( p*m ), nrow = m )
  #split the m factors into K blocks
  pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) )
  Y = matrix( 0, nrow = m, ncol = n )
  
  for( x in 1:K){
    Y[,pop.indx.group[[x]]] = (2*x - (K+1))*strength + 
      matrix( rnorm( m*length(pop.indx.group[[x]]) ), 
              ncol = length(pop.indx.group[[x]]) )
  }
  
  er = 1/rgamma( 1, a.er, b.er)
  Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
  final.weights = sigma*(Q*(Q>=0))
  data = lapply( lscounts, function(counts) apply( final.weights, 2, function(x) rmultinom( 1, counts, prob = x ) ) )
  
  return( list( sigma = sigma, Q = Q, data = data,
                X.tru=X, Y.tru = Y, er = er ) )
}

get_RV_inner = function( all.cov ){
  RV_mt = matrix( 1, nrow = length( all.cov ), ncol = length( all.cov ) )
  for( i in 1:(length( all.cov )-1) ){
    print( i )
    for( j in (i+1):length(all.cov) ){
      rv_ij = rv.coef( all.cov[[i]], all.cov[[j]] )
      RV_mt[i,j] = rv_ij
      RV_mt[j,i] = rv_ij
    }
  }
  return( RV_mt )
}

statis = function( all.cov, n.vis=2 ){
  RV_mt = get_RV_inner( all.cov )
  RV_eigen = eigen( RV_mt )
  cat("1st PC of RV explained\n")
  cat( RV_eigen$values[1]/sum(RV_eigen$values) )
  
  wt_comp = abs( RV_eigen$vector[,1] )/sum( abs( RV_eigen$vector[,1] ) )
  mt_comp = Reduce( '+', lapply( 1:length(wt_comp), function(x) all.cov[[x]]*wt_comp[x] ) )
  
  
  #get first two axes of comp mt
  eigen_comp = eigen( mt_comp )
  print( eigen_comp$values[1:n.vis] )
  list( coord = Reduce( rbind, lapply( all.cov, function(x) 
    x%*%eigen_comp$vectors[,1:n.vis]%*%diag(1/sqrt(eigen_comp$values[1:n.vis])) ) ),
    ratio = eigen_comp$values[1]/eigen_comp$values[2] )
}

plot_statis = function( statis_co, n, n.class, label=F ){
  require( ggplot2 )
  n.rep = nrow(statis_co)/n
  plot.data = data.frame( x = statis_co[,1], y = statis_co[,2], 
                          sub = rep( 1:n, n.rep ),
                          class = rep( cut( 1:n, n.class, labels = 1:n.class ), n.rep ) )
  contr = ggplot( data = plot.data, aes( x = x, y = y, color = class, group = sub ) ) + 
    geom_density2d()
  if( label ){
    centroids1=tapply(plot.data$x,plot.data$sub,mean)
    centroids2=tapply(plot.data$y,plot.data$sub,mean)
    centroids.df=data.frame( Comp1=centroids1, Comp2=centroids2, sample= 1:n )
    contr=contr+
      with(centroids.df, annotate(geom="text",label = sample , x = Comp1, y = Comp2),size=2.5)
  }
  contr + theme_bw() + 
    theme( axis.title = element_text( size = 15 ), 
                              axis.text = element_text( size = 12 ),
                              legend.title = element_text( size = 15),
                              legend.text = element_text( size = 12 ) ) + 
    xlab("Compromise axis 1") + ylab("Compromise axis 2")
}

#set prior for sigma
alpha = 10
n = 22
p = 68
sigma.value = seq(0.001,0.999,0.001)
tmp = c( 0,pbeta( sigma.value, alpha/p, 1/2-alpha/p ) )
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
plot( sigma.value, sigma.prior )

#simulate data
dat.sim = sim.for.contour( 1000, n, p, 
                           m = 3, K = 2, a.er = 1, 
                           b.er = 0.3, sigma.value,
                           sigma.prior, strength = 3 )

#do MCMC
hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = 22 )
for( i in 1:length(dat.sim$data) ){
  print(i)
  i = 1
  start = list( sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ),
                T.aug = colSums( dat.sim$data[[i]] ),
                Q = matrix( 0.5, nrow = p, ncol = n ),
                X = matrix( rnorm( p*hyper$m ), nrow = hyper$m ),
                Y = matrix( rnorm( n*hyper$m ), nrow = hyper$m ),
                er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                delta = c( rgamma( 1, shape = 3, rate = 1 ), 
                           rgamma( hyper$m-1, shape = 4, rate = 1 ) ),
                phi = matrix( rgamma( n*hyper$m, shape = 3/2, rate = 3/2 ), nrow = n ) )
  main.mcmc.shrink( dat.sim$data[[i]], start, hyper, sigma.value, sigma.prior, 
                    paste( "C:/Users/Renboyu/Desktop/sim_contour/", "res", sep = "/"), 
                    save.obj = c("sigma","Q", "Y","er"), step = 50000, thin=5 )
  all.res = lapply( list.files( "C:/Users/Renboyu/Desktop/sim_contour/", 
                                pattern = "res", full.names = T), readRDS )
  all.corr = lapply( all.res, function(x) cov2cor( t(x$Y)%*%x$Y + diag( rep( x$er, ncol(x$Y) ) ) ) )
  all.corr.use = all.corr[sample(1:length(all.corr),size = 1000, replace = F)]
  
  statis.res = statis( all.corr.use, 2 )
  
  pic = plot_statis( statis.res$coord, n, 2, label = T )
  ggsave( paste( "C:/Users/Renboyu/Desktop/factor_sim/contour/plot_", i, ".pdf", sep = "" ), pic )
}
