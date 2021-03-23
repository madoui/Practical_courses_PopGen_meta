# Library of R functions to use duringin the course
# basic function to compute FST

library(ggplot2)
library(vegan)
library(reshape2)
library(metaVaR)
library(MM4LMM)
library(viridis)
library(pcadapt)
library(ggrepel)

# filtering dataset

is.nan.data.frame <- function(x) {
  
  do.call(cbind, lapply(x, is.nan))
}
CreateFreqMatrix <- function(d,poplist) {
  matrixf = matrix(0,nrow = nrow(d),ncol = length(poplist))
  n = 0
  for (i in poplist) {
    n = n + 1 
    x = d[,grep(as.character(i),colnames(d),value = T)]
    matrixf[,n] = x[,2]/(x[,1]+x[,2])
  }
  colnames(matrixf) = poplist
  rownames(matrixf) = rownames(d)
  return(as.matrix(matrixf))
}
CreateExpressionMatrix <- function(d,poplist) {
  matrixf = matrix(0,nrow = nrow(d),ncol = length(poplist))
  n = 0
  for (i in poplist) {
    n = n + 1 
    x = d[,grep(as.character(i),colnames(d),value = T)]
    matrixf[,n] = x[,4]/(x[,3]+x[,4])
  }
  colnames(matrixf) = poplist
  rownames(matrixf) = d[,1]
  return(as.matrix(matrixf))
}
FilteringCov <- function(MainData,poplist,dev=2, minCov=5, maxCov = 150) {
  CreateCovMatrix <- function(d,poplist) {
    matrixcov = matrix(0,nrow = nrow(d),ncol = length(poplist))
    n = 0
    for (i in poplist) {
      n = n + 1  
      matrixcov[,n] = d[,grep(paste("GB",as.character(i),sep=""),colnames(d),value = T)] +
        d[,grep(paste("GA",as.character(i),sep=""),colnames(d),value = T)]
      
    }
    colnames(matrixcov) = paste("p",poplist,sep="")
    rownames(matrixcov) = rownames(d)
    return(as.matrix(matrixcov))
  }
  subsetCov2 <- function  (cov , dev = 2, minCov = 4, maxCov = 150){
    N = nrow(cov)
    M = ncol(cov)
    sd = apply (cov , 2, sd)
    med = apply (cov , 2, median)
    covOk = matrix (rep ( 0, N*M), nrow = N, ncol = M )
    rownames (covOk) = rownames (cov)
    for (i in 1:N){
      for (j in 1:M){
        if (i > N || j > M){break}
        if ( cov[i,j] > med[j]-dev*sd[j] && cov[i,j] < med[j]+dev*sd[j] && cov[i,j] > minCov && cov[i,j] < maxCov ){
          covOk[i,j] = 1
        }
      }
    }
    ok = apply ( covOk , 1, function (x) all ( x == 1 ) )
    cov = cov[ ok == 1, ]
    return (cov)
  }
  MatrixCov = CreateCovMatrix(MainData,poplist)
  New = subsetCov2(MatrixCov,dev,minCov,maxCov)
  data = MainData[rownames(MainData) %in% rownames(New),]
  return(data)
}
FilteringFreqExpression <- function(d, poplist,GAF=0.1, TAF=0.05) {
  
  #Filtre >GAF/TAF dans au moins une station
  MinAF <- function(matrix,df,minAF=0.1) {
    test = matrix(0, nrow=length(rownames(df)),ncol=length(colnames(matrix)))
    rownames(test)= rownames(df)
    for (i in 1:length(colnames(matrix))) {
      for (j in 1:nrow(matrix)) {
        
        if (matrix[j,i] >= minAF ) {
          test[j,i] = 1
        }
        else {
          test[j,i] = 0
        }
      }}
    
    c = data.frame(rownames(df),apply(test,1,sum))
    c_ok = c[c[,2] !=0,]
    
    data = df[rownames(df) %in% c_ok[,1],]
    return(data)
  }
  fmatrix= CreateFreqMatrix(d,poplist)
  w = MinAF(fmatrix,d,GAF)
  fmatrix= CreateFreqMatrix(w,poplist)
  w = MinAF(1-fmatrix,w,GAF)
  
  ematrix <- CreateExpressionMatrix(w,poplist)
  w = MinAF(ematrix,w,TAF)
  ematrix= CreateExpressionMatrix(w,poplist)
  w = MinAF(1-ematrix,w,TAF)
  
  print(nrow(w))
  return(w)
}


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
  #if ( is.matrix (p) == 1){
  #  return (paste (substitute(x),"is a matrix"))
  #}
  
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

plotFST<-function(pairwiseFST){
  fstTable = data.frame( t(combn(names(get_lower_tri(pairwiseFST)),2)), 
                         fst=t(pairwiseFST)[lower.tri(pairwiseFST)] )
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
  
  Genetic_Distance = melt(get_lower_tri(as.matrix(FstMatrix)))
  Genetic_Distance[is.na(Genetic_Distance)] = 0
  Genetic_Distance = Genetic_Distance[Genetic_Distance[,3] !=0,]
  
  diag(EnvMatrix) = NA
  Env_Distance = na.omit(melt(get_lower_tri(as.matrix(EnvMatrix))))
  
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
           labs(title = "Bivariate plot of Median-Pairwise Fst vs distance",
                subtitle = paste(paste("Mantel r statistic = ",round(res.mantel$statistic,digits=2),sep=""),
                                 paste("pvalue = ",formatC(res.mantel$signif,format = "e",digits=2),sep=""),sep="\n")))
}

# MVS functions ; variation partitioning and plotting

LinearMixedModelMVS = function(MVS,LagMatrix,EnvList) {
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  Results = list()
  
  Fst_matrix = MVS@pwFst
  subpoplist = MVS@pop
  MVSname = MVS@name
  print(MVSname)
  
  diag.boolean<-FALSE
  EnvList2 = list()
  for (env in 1:length(EnvList)) {
    e = EnvList[[env]]
    t = as.data.frame(e[e[,1] %in% subpoplist,-1])
    row.names(t) = e[e[,1] %in% subpoplist,1]
    t = as.matrix(t[order(rownames(t)),])
    row.names(t) = subpoplist
    t_matrix = as.matrix(dist(t,"euclidean")) 
    EnvList2[[env]] = t_matrix[upper.tri(t_matrix,diag=diag.boolean)]
    names(EnvList2)[env] = names(EnvList)[env]
  }
  
  
  Y<-Fst_matrix[upper.tri(Fst_matrix,diag=diag.boolean)]
  X.lagrangian<-c(LagMatrix[subpoplist,subpoplist][upper.tri(LagMatrix[subpoplist,
                                                                       subpoplist],diag=diag.boolean)])
  
  X = as.data.frame(matrix(ncol=length(EnvList2)+2,nrow=length(Y),0))
  X[,1] = Y
  X[,2] = scale(X.lagrangian)
  colnames(X)[1:2] = c("Fst","Lagrangian")
  for (w in 1:length(EnvList2)) {
    X[,w+2] = scale(EnvList2[[w]])
    colnames(X)[w+2] = names(EnvList2)[w]
  }
  
  
  # table1 = stargazer(lm(Fst~Temperature,data=X),
  #                    lm(Fst~Lagrangian,data=X),
  #                    lm(Fst~Silicium, data=X),
  #                    lm(Fst~Salinity, data=X),
  #                    lm(Fst~.,data=X), type = 'text')
  ZL = lapply(EnvList2,function(X) as.matrix(data.frame(scale(X))))
  ZL[[length(ZL)+1]] = as.matrix(data.frame(scale(X.lagrangian)))
  ZL[[length(ZL)+1]] = diag(length(Y))
  names(ZL)[c(length(EnvList2)+1,length(EnvList2)+2)] = c("Lagrangian","Error")
  
  VL  = lapply(ZL,function(X) X = matrix(1))
  VL$Error = diag(length(Y))
  Model<-MMEst(Y=Y,X=NULL,Cofactor = NULL,
               ZList=ZL,   
               VarList=VL,CritVar = 0.001,
               CritLogLik = 0.01,MaxIter = 100000)
  #AnovaTest(Model , TestedCombination=NULL , Type = "TypeI" ,
  #                Cofactor = NULL , X = NULL , formula = NULL , VarList = VL)
  
  ##Trace Fst
  TotalVariance=var(Y)*(length(Y)-1)/length(Y)
  variances<-c(sum(Model$NullModel$VarBeta,Model$NullModel$Sigma2),
               Model$NullModel$VarBeta,Model$NullModel$Sigma2) 
  MM.table<-data.frame(variances=variances,percentage=variances/sum(variances[-1]))
  row.names(MM.table)<-c("Fst","FixedEffect",names(ZL)[-length(ZL)],"Error")
  table2 = MM.table
  #knitr::kable(MM.table)
  
  Results = list(MVSname = paste(MVSname,sep="_"), Model = Model, Data = X, 
                 MME=table2,TotalVariance = TotalVariance)
  print(table2)
  return(Results) 
} 

PlotLMM = function(Results) {
  print(ggplot(data= Results$MME[-c(1,2),], 
               aes(x=rownames(Results$MME[-c(1,2),]),
                   y=percentage)) + 
          geom_bar(stat = "identity",fill="goldenrod",col="black") + 
          theme_bw() + ggtitle(Results$MVSname) + xlab("Parameters") + 
          ylab("Proportion of explained variation") + ylim(values=c(0,1))) 
}






