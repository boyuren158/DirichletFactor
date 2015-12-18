mt.trace = function( A ){
  sum(diag(A))
}

rv.coef = function( A, B ){
  mt.trace(A%*%B)/sqrt(mt.trace(A%*%A)*mt.trace(B%*%B))
}

gg.heatmap <- function( cormt, x.lab="", y.lab="", title="" ){
  require( ggplot2 )
  require( reshape2 )
  plot.dat = melt( cormt )
  ggplot( dat = plot.dat, aes( x=Var1, y=Var2, fill=value ) ) + geom_tile(color = "black") + 
    scale_fill_gradient2( low="steelblue", mid = "white", high = "red", limits=c(-1, 1) ) + 
    theme_bw() + coord_fixed() + 
    theme( axis.title = element_text(face="bold"), plot.title = element_text(face="bold"),
           axis.text = element_blank(), axis.ticks = element_blank(),
           panel.border = element_blank(), panel.grid = element_blank() ) + 
    xlab(x.lab) + ylab(y.lab) + ggtitle( title )
}

multiplot <- function(..., plotlist=NULL, file, cols=1, byrow = F, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),byrow = byrow,
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}