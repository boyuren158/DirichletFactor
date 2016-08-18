# script used to perform MCMC on GlobalPatterns data (collapsed to genus level)
# also used to generate figure 4

source("MCMC_fns.R")
library( vegan )
library( phyloseq )
library( DistatisR )
library( fpc )

data( GlobalPatterns )
gp_genus = tax_glom( GlobalPatterns, "Genus")
gp_metadata = sample_data( GlobalPatterns )
data = otu_table( gp_genus )

alpha = 10
n = ncol(data)
p = nrow(data)
sigma.value = seq(0.001,0.999,0.001)
tmp = c( 0,pbeta( sigma.value, alpha/p, 1/2-alpha/p ) )
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
plot( sigma.value, sigma.prior )

hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = n )
start = list( sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ),
              T.aug = colSums( data ),
              Q = matrix( 0.5, nrow = p, ncol = n ),
              X = matrix( rnorm( p*hyper$m ), nrow = hyper$m ),
              Y = matrix( rnorm( n*hyper$m ), nrow = hyper$m ),
              er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
              delta = c( rgamma( 1, shape = 3, rate = 1 ), 
                         rgamma( hyper$m-1, shape = 4, rate = 1 ) ),
              phi = matrix( rgamma( n*hyper$m, shape = 3/2, rate = 3/2 ), nrow = n ) )

main.mcmc.shrink( data, start, hyper, sigma.value, sigma.prior, 
                  paste( "C:/Users/Renboyu/Desktop/gp_MCMC/", "res", sep = "/"), 
                  save.obj = c("sigma","Q", "Y","er"), step = 50000, thin=10 )

# plot module				  
# get result (all MCMC runs)
all.res = lapply( list.files( "gp_MCMC/", pattern = "res_", full.names = T ),
                  readRDS )

#get bray-curtis distance matrix
all.bc = lapply( all.res, function(x){
  weights = x$Q^2*(x$Q>0)*x$sigma
  w.norm = t(weights)/colSums(weights)
  vegdist( w.norm, method = "bray" )
})

use.idx = sample( 1:length(all.res), 1000, replace = F )

#distatis
all.bc.ls = lapply( all.bc[use.idx], as.matrix )
all.bc.mt = array( unlist( all.bc.ls ), dim = c( dim( all.bc.ls[[1]] ), length( all.bc.ls ) ) )
bc.distatis = distatis( all.bc.mt, nfact2keep = 3 )
bc.ev = eigen(bc.distatis$res4Splus$Splus)$values
bc.prop = bc.ev[1:3]/sum(bc.ev)
bc.coord = apply( bc.distatis$res4Splus$PartialF, 2, rbind )

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

pcoa_credible = function( all.coord, labels, size = 2.5, annotation = F ){
  labels.color = gg_color_hue( length(levels(labels) ) )

  n.rep = nrow( all.coord )/length( labels )
  plot.data = data.frame( x = all.coord[,1], y = all.coord[,2], 
                          group = rep( 1:length(labels), n.rep),
                          type = rep( labels, n.rep ) )
  contr = ggplot() + geom_density2d (data = plot.data, aes( x=x, y=y, group = group, color =  type) ) +
    xlim(-0.4,0.4) + ylim(-0.4,0.4) +
    theme_bw() + scale_color_discrete( guide = F )
  if( annotation ){
    x.annote = tapply( plot.data$x, plot.data$group, mean )
    y.annote = tapply( plot.data$y, plot.data$group, mean )
    annote.data = data.frame( x=x.annote, y=y.annote, 
                              color = labels.color[as.numeric(labels)],
                              labels = as.character( labels ) )
    contr = contr + geom_point( data = annote.data, aes( x=x, y=y, color = labels )  ) +
      #with(annote.data, annotate(geom="segment", xend = x, yend = y, x = x+0.05, y = y+0.05, arrow=arrow() ) ) + 
      with(annote.data, annotate(geom="text", x = x+0.01 , y = y, color = color, label = labels, size = 8) )
  }
  contr
}

p.bc12 = pcoa_credible( bc.coord[,1:2], gp_metadata$SampleType, annotation = T, size = 0.5 ) +
  xlab(paste("Compromise axis 1 (",round(bc.prop[1]*100,1),"%)",sep="")) + 
  ylab(paste("Compromise axis 2 (",round(bc.prop[2]*100,1),"%)",sep="")) + 
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ) ) +
  coord_fixed()
p.bc13 = pcoa_credible( bc.coord[,c(1,3)], gp_metadata$SampleType, annotation = T, size = 0.5 ) +
  xlab(paste("Compromise axis 1 (",round(bc.prop[1]*100,1),"%)",sep="")) + 
  ylab(paste("Compromise axis 3 (",round(bc.prop[3]*100,1),"%)",sep="")) +
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ) ) +
  coord_fixed()
p.bc23 = pcoa_credible( bc.coord[,c(2:3)], gp_metadata$SampleType, annotation = T, size = 0.5 ) +
  xlab(paste("Compromise axis 2 (",round(bc.prop[2]*100,1),"%)",sep="")) + 
  ylab(paste("Compromise axis 3 (",round(bc.prop[3]*100,1),"%)",sep="")) + 
  theme( axis.text = element_text( size = 15 ), axis.title = element_text( size = 18 ) ) +
  coord_fixed()

ggsave("../writeup/gp_figure/raw/b12.pdf", p.bc12, width = 7, height = 7)
ggsave("../writeup/gp_figure/raw/b13.pdf", p.bc13, width = 7, height = 7)
ggsave("../writeup/gp_figure/raw/b23.pdf", p.bc23, width = 7, height = 7)

#clustering plot
clustering_posterior = function( all.dist ){
  res.cluster = lapply( all.dist, function(x){
    x = as.matrix(x)
    cluster = pamk( x, krange = seq(2,ncol(x)-1), diss = T )
    cluster.id = cluster$pamobject$clustering
    outer( cluster.id, cluster.id, "==" )
  })
  res.cluster.mt = array( unlist( res.cluster ), dim = c( dim( res.cluster[[1]] ), length( res.cluster ) ) )
  apply( res.cluster.mt, 1:2, mean )
}

#clustering using bc
bc.clustering = clustering_posterior( all.bc )

#plot the clustering scheme
clustering_plot = function( res.mt, labels ){
  require( reshape2 )
  require( stringr )
  
  res.mt.order = res.mt[order( labels ),order( labels )]
  plot.res = melt( res.mt.order )
  axis.labels = levels( labels )
  axis.pos = cumsum( table( labels ) ) + 0.5
  axis.text.pos = as.vector( axis.pos - table( labels )/2 )
  x.delim = data.frame( x = axis.pos, xend = axis.pos, yend = rep( 0.5, length(axis.pos) ), y = rep( 0, length(axis.pos) ) )
  y.delim = data.frame( y = axis.pos, yend = axis.pos, xend = rep( 0.5, length(axis.pos) ), x = rep( 0, length(axis.pos) ) )
  
  ggplot() + geom_tile( data = plot.res, aes( x=Var1, y=Var2, fill = value ), color = "black" ) + 
    scale_fill_gradient( low = "white", high = "red" ) + guides(fill=guide_legend(title="Posterior\nprobability")) +
    annotate( geom="text", x = axis.text.pos, y = 0, label = axis.labels, size = 6, angle = 30, hjust = 1 ) + 
    annotate( geom="text", y = axis.text.pos, x = 0, label = axis.labels, size = 6, angle = 30, hjust = 1 ) + 
    theme_bw() + coord_fixed() + 
    theme( axis.text = element_blank( ),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.title = element_text(size = 18,face = "bold"),
           legend.text = element_text(size = 15) ) + 
    geom_segment( data = x.delim, aes(x = x, xend = xend, y = y, yend = yend ) ) + 
    geom_segment( data = y.delim, aes(x = x, xend = xend, y = y, yend = yend ) )
}

p.bc.cluster = clustering_plot( bc.clustering, gp_metadata$SampleType )
ggsave( "cluster.pdf", p.bc.cluster, width = 7, height = 7 )
