# Library of R functions to use duringin the course
# basic function to compute FST

library(ggplot2)
library(vegan)
library (reshape2)
fst <- function (p) {
  varID = rownames(p)
  fst = apply( p, 1, function (x) var(x)/(mean(x)*(1-mean(x))) )
  fst[ fst > 1 ] = 1
  varID = varID [complete.cases( fst )]
  varID = varID [!is.na(fst)]
  fst[complete.cases( fst )]
  fst = na.omit( fst)
  fst = as.vector( fst )
  names(fst) = varID
  return ( fst =  fst )
}

# function to compute pairwie-FST
pwFst <- function (p){
  #check if the input is a matrix
  if ( is.matrix (p) == 1){
    return (paste (substitute(x),"is a matrix"))
  }
  
  #create an empty matrix filled with 1
  mfst = matrix( rep( 0, ncol(p)^2), nrow = ncol(p), ncol = ncol(p) )
  popName = colnames(p)
  colnames(mfst) = popName
  rownames(mfst) = popName
  max = ncol(p)
  
  #compute pairwise Fst
  for (i in 1:max ){
    for ( j in i+1:max ){
      if ( j>max || i == j ){
        break
      }
      else{
        fst =  fst( data.frame( p[,i], p[,j] ) )
        mfst[i,j] = median(fst)
        mfst[j,i] = median(fst)
      }
    }
  }
  return (as.data.frame(mfst))
}


get_lower_tri<-function(m){
  m[upper.tri(m)] <- NA
  return(m)
}

plotFST<-function(pairewiseFST){
  fstTable = data.frame( t(combn(names(get_lower_tri(pairwiseFST)),2)), fst=t(pairwiseFST)[lower.tri(pairwiseFST)] )
  colnames(fstTable) = c("Population1", "Population2", "Fst")
  Fst<-c();c<-NULL
  p = ggplot(data = fstTable, aes_(x=~Population1, y=~Population2, fill = ~Fst))+
    geom_tile(color = "white")+
    geom_text(aes(label = round(Fst, 2)))+
    scale_fill_gradient2(low = "white", high = "darkred", midpoint = 0, limit = c(0,1)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    coord_fixed()
  plot(p)
}

LK <- function (p){
  n = ncol(p)
  fst = fst(p)
  lk = (n-1)*fst/mean(fst)
  pv = pchisq (q = lk, df = n-1, lower.tail = F)
  qv = p.adjust (pv, method = "BH")
  d = cbind( names(fst), fst, lk, pv,  qv )
  colnames(d) = c( "varID", "Fst", "LK", "p_value","q_value" )
  d = as.data.frame(d)
  d = as.data.frame(apply(d, 2, function(x) as.numeric(as.character(x))))
  return ( d )
}

# Do and plot Mantel test 
MantelPlot <- function(EnvMatrix,FstMatrix)  {
  
  Genetic_Distance = melt(get_lower_tri(FstMatrix))
  Genetic_Distance[is.na(Genetic_Distance)] = 0
  Genetic_Distance = Genetic_Distance[Genetic_Distance[,3] !=0,]
  
  Env_Distance = melt(get_lower_tri(EnvMatrix))
  Env_Distance[is.na(Env_Distance)] = 0
  Env_Distance = Env_Distance[Env_Distance[,3] !=0,]
  
  
  res.mantel = mantel(vegdist(FstMatrix),vegdist(EnvMatrix),
                      permutations = 10000,na.rm = TRUE)
  
  
  toplot = data.frame(paste(Genetic_Distance[,1],Genetic_Distance[,2],sep="-"),
                      Genetic_Distance[,3],Env_Distance[,3])
  colnames(toplot)=c("Population","Fst","Env")
  
  print( ggplot(toplot,aes(x=toplot[,3],y=Fst,label=Population)) + 
           geom_smooth(method="lm",size=2,level = 0.95) + 
           geom_point(shape=18, color="black", size=4) +
           geom_label_repel(force=10) +
           theme(
             panel.background = element_rect(fill = "transparent",colour = NA)# bg of the panel
             , plot.background = element_rect(fill = "transparent", colour = NA) # bg of the plot
             , panel.grid.major = element_blank() # get rid of major grid
             , panel.grid.minor = element_blank() # get rid of minor grid
             , axis.line = element_line(colour = "black",size = 2)
             , axis.text = element_text(size = 15,face="bold")
             , axis.title=element_text(size=14,face="bold")
             , plot.title = element_text(hjust = 0.5)) +
           labs(title = paste("Bivariate plot of Median-Pairwise Fst vs ",env,sep=""),
                subtitle = paste(paste("Mantel r statistic = ",round(res.mantel$statistic,digits=2),sep=""),
                                 paste("pvalue = ",formatC(res.mantel$signif,format = "e",digits=2),sep=""),sep="\n")) +
           stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE,xpos = max(toplot[,3])/2,ypos = 0.12) )
}