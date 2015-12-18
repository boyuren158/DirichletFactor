# script used to generate figure 2(e)
# results from prediction_pw.R are assumed to be in the same folder
# named after their parameter specification
# e.g. "sim05_mis_0_5.rds" means correlation between biological samples is 0.5, misspecified data generation model
# a = 0.5

library( ggplot2 )

read.depth = seq(10,100,10)
all.files = list.files( "./", pattern = "rds" )
all.res = lapply( all.files, readRDS )
names( all.res ) = sapply( all.files, function(x) strsplit(x, split = ".", fixed = T )[[1]][1] )
all.gains = lapply( all.res, function(x) rowMeans( sapply( x, rowMeans ) ) )

plot.data = data.frame( read.depth = rep( read.depth,length(all.gains) ), 
                        gain = unlist( all.gains ),
                        similarity = rep( as.factor(c(0.5,0.75,0.95)), each = 4*length(read.depth) ),
                        Q.power = rep( rep( as.factor(c(2,1,0.5,3)), each = length(read.depth) ), 3 ) )

cor.pw = ggplot( data = plot.data, aes(x=read.depth, y=gain, color=similarity, shape = Q.power) ) + 
  geom_line() + geom_point(size=5) + theme_bw() +
  xlab("Read depth") + ylab("Gain in estimate accuracy" ) + scale_shape_manual( values = c(15:17,4)) +
  scale_x_continuous( breaks = read.depth ) + 
  theme( axis.title = element_text( size = 18 ), axis.text= element_text( size = 15 ),
         legend.title = element_text( size = 15, face = "bold" ), plot.title = element_text( size = 15, face = "bold") )
cor.pw
ggsave( "prediction_pw_cor.pdf", cor.pw, width = 8, height = 6.5 )
