##Mantel Osimilis
library(ggplot2)
library(vegan)

MantelSimilis <- function(EnvMatrix,FstMatrix)  {
  
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
MantelSimilis(Env_parameters_filtered,Median_fst_Matrix)
