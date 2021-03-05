# Library of R functions to use duringin the course
# basic function to compute FST

library(ggplot2)

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

