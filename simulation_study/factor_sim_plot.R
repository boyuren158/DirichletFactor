# script to plot figure 2(a-d)
# Assume the output of factor_sim.R/factor_sim_mis.R are collected in three RDS object
# for all three true factor scenarios (true factor = 3, 6, 9), with names:
# factor_3.rds/factor_3_mis.rds
# factor_6.rds/factor_6_mis.rds
# factor_9.rds/factor_9_mis.rds
# factor_3.rds, for example, stores a list of lists. Each entry includes posterior samples of saved model parameters at one iteration

source( "../utilites.R" )

library( ggplot2 )

#calculate the posterior mean correlation
get_corr = function( all.res ){
  lapply( all.res, function(y){ lapply( y, function(x){
  post.corr = lapply( x, function(mt) cov2cor( t(mt$Y)%*%mt$Y + diag( rep( mt$er, ncol( mt$Y ) ) ) ) )
  post.corr.mt = array( unlist( post.corr ), dim = c( dim( post.corr[[1]] ), length( post.corr ) ) )
  apply( post.corr.mt, c(1,2), mean )} ) } )
}

# get rv-coefficients for estimates of correlation matrix
compare_correlation = function( filename ){
  all.res = readRDS( filename )
  cat( "Done readin\n" )
  corr.est = get_corr( lapply( all.res, function(x) x$res ) )
  cat( "Done corr est\n" )
  corr.tru = lapply( lapply( all.res, function(x) x$truth ), function(x) cov2cor( t(x$Y.tru)%*%x$Y.tru + diag( rep( x$er, ncol(x$Y.tru) ) ) ) )
  
  sapply( 1:length(corr.est), function(x){
    sapply( corr.est[[x]], function( mt ) rv.coef( mt, corr.tru[[x]] ) )
  })
}

corr.mt.est.3 = compare_correlation( "factor_3" )
corr.mt.est.6 = compare_correlation( "factor_6" )
corr.mt.est.9 = compare_correlation( "factor_9" )
corr.mt.est.3.mis = compare_correlation( "factor_3_mis" )
corr.mt.est.6.mis = compare_correlation( "factor_6_mis" )
corr.mt.est.9.mis = compare_correlation( "factor_9_mis" )

#boxplots to summarize all the RV coeffients
plot.corr.est = data.frame( rv = ( c( t(corr.mt.est.3), t(corr.mt.est.6), t(corr.mt.est.9) ) ),
                            total = as.factor( rep( rep( c(1000,10000,100000), each = ncol( corr.mt.est.3 ) ), 3 ) ),
                            tru.factor = as.factor( rep( c(3,6,9), each = 3*ncol( corr.mt.est.3 ) ) ),
                            group = as.factor( rep( 1:9, each = ncol( corr.mt.est.3 ) ) ) )
plot.corr.est.mis = data.frame( rv = ( c( t(corr.mt.est.3.mis), t(corr.mt.est.6.mis), t(corr.mt.est.9.mis) ) ),
                                total = as.factor( rep( rep( c(1000,10000,100000), each = ncol( corr.mt.est.3 ) ), 3 ) ),
                                tru.factor = as.factor( rep( c(3,6,9), each = 3*ncol( corr.mt.est.3 ) ) ),
                                group = as.factor( rep( 1:9, each = ncol( corr.mt.est.3 ) ) ) )

saveRDS(list(plot.corr.est, plot.corr.est.mis), "plot_corr.rds")
tmp = readRDS( "plot_corr.rds" )
plot.corr.est = tmp[[1]]
plot.corr.est.mis = tmp[[2]]

#add sqrt of rv to make better visualization
p.corr = ggplot( data = plot.corr.est, aes( x = total, y = rv, group = group, dodge = tru.factor, fill = tru.factor ), guide = F ) + 
  geom_boxplot() + xlab("Total counts") + ylab("Estimate accuracy") + scale_y_continuous(limits=c(0.5, 1)) + theme_bw() +
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ),
         legend.title = element_text( size = 15 ), legend.text = element_text( size = 15 ),
         plot.title = element_text( size = 18, face = "bold") ) +
  labs( fill = "No. factors") + ggtitle( "Correctly specified model" )

p.corr.mis = ggplot( data = plot.corr.est.mis, aes( x = total, y = rv, group = group, dodge = tru.factor, fill = tru.factor ) ) + 
  geom_boxplot() + xlab("Total counts") + ylab("Estimate accuracy") + scale_y_continuous(limits=c(0.5, 1)) + theme_bw() +
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ),
         legend.title = element_text( size = 15 ), legend.text = element_text( size = 15 ),
         plot.title = element_text( size = 18, face = "bold")) +
  labs( fill = "No. factors") + ggtitle( "Misspecified model" )


#all pca
get_pca_plot = function( filename ){
  all.res = readRDS( filename )
  cat( "Done reading.\n" )
  lapply( all.res, function(x){
    lapply( x$res, function(y) 
      sapply( y, function(z){tmp = princomp(t(z$Y));tmp$sdev^2/sum(tmp$sdev^2)}) )
  })
}

pca.3 = get_pca_plot( "factor_3" )
pca.6 = get_pca_plot( "factor_6" )
pca.9 = get_pca_plot( "factor_9" )
pca.3.mis = get_pca_plot( "factor_3_mis" )
pca.6.mis = get_pca_plot( "factor_6_mis" )
pca.9.mis = get_pca_plot( "factor_9_mis" )

pca.3.sub = array( unlist( lapply( pca.3, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]
rm( pca.3 )

pca.6.sub = array( unlist( lapply( pca.6, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]
pca.9.sub = array( unlist( lapply( pca.9, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]

pca.3.sub.mis = array( unlist( lapply( pca.3.mis, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]
pca.6.sub.mis = array( unlist( lapply( pca.6.mis, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]
pca.9.sub.mis = array( unlist( lapply( pca.9.mis, function(x) 
  apply( x[[1]], 1, function(vec) c( mean(vec), quantile( vec, probs = c(0.025,0.975 ) ) ) ) ) ),
  dim = c( 3,10,50 ) )[1,,]

pca.plot.data = data.frame( ve = sqrt( c( pca.3.sub, pca.6.sub, pca.9.sub ) ), 
                            comp = as.factor( rep( rep( 1:10, 50 ), 3) ),
                            tru.factor = as.factor( rep( c(3,6,9), each = 500 ) ) )
pca.plot.data.mis = data.frame( ve = sqrt( c( pca.3.sub.mis, pca.6.sub.mis, pca.9.sub.mis ) ), 
                            comp = as.factor( rep( rep( 1:10, 50 ), 3) ),
                            tru.factor = as.factor( rep( c(3,6,9), each = 500 ) ) )

x = rep( c( 0.75, 1, 1.25 ), 10 ) + rep( 0:9, each = 3 ) - 0.125
xend = rep( c( 0.75, 1, 1.25 ), 10 ) + rep( 0:9, each = 3 ) + 0.125
y = rep( 0.6, 30 )
yend = rep( 0.6, 30 )

saveRDS( list(pca.plot.data, pca.plot.data.mis), "pca_res.rds" )
tmp = readRDS( "pca_res.rds" )
pca.plot.data = tmp[[1]]
pca.plot.data.mis = tmp[[2]]

p.pca = ggplot( data = pca.plot.data, aes( x = comp, y = ve, dodge = tru.factor, fill = tru.factor ) ) + 
  geom_boxplot(outlier.size = 0) + xlab("Principal component") + 
  ylab("Sqrt variance explained") + theme_bw() +
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ),
         legend.title = element_text( size = 15 ), legend.text = element_text( size = 15 ),
         plot.title = element_text( size = 18, face = "bold" ) ) +
  labs( fill = "No. factors") + ggtitle( "Correctly specified model")

p.pca.mis = ggplot( data = pca.plot.data.mis, aes( x = comp, y = ve, dodge = tru.factor, fill = tru.factor ) ) + 
  geom_boxplot(outlier.size = 0) + xlab("Principal component") + 
  ylab("Sqrt variance explained") + theme_bw() +
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ),
         legend.title = element_text( size = 15 ), legend.text = element_text( size = 15 ),
         plot.title = element_text( size = 18, face = "bold" )) +
  labs( fill = "No. factors") + ggtitle( "Misspecified model" )

get_pca_tru_plot = function( n, all.tru ){
  raw = sapply( all.tru, function(x){
    res = princomp( t( x$Y ) )
    res$sdev^2/sum(res$sdev^2)
  })
  plot.raw = apply( raw, 1, function(x) c(mean(x), quantile(x, probs = c(0.025,0.975) ) ) )
  ret.mt = matrix( 0, nrow = 3, ncol = n )
  ret.mt[1:nrow(plot.raw),1:ncol(plot.raw)] = plot.raw
  data.frame( mean.tru = ret.mt[1,], lower.tru = ret.mt[2,], upper.tru = ret.mt[3,] )
}