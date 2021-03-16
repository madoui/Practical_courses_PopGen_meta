
######MVS and variation partitioning#####
library("metaVaR")
library("ggplot2")
library("MM4LMM")
library("reshape2")
library("viridis")

###Read MVS
MVS.List = readMvsList("data/MVS")

###Fst analysis

for(i in 1:4){
plotFST(pwFst(MVS.List[[i]]@freq))
}

###Read Environmental table

EnvTable = read.table("data/WorldOceanAtlas_data.txt",header = T)
envlist = list(Temperature = EnvTable[,c(1,7)],
             Salinity = EnvTable[,c(1,8)],
             Silicate = EnvTable[,c(1,4)],
             Po4 = EnvTable[,c(1,6)], 
             No3 = EnvTable[,c(1,5)])

LagrangianMatrix = read.table("data/Lagrangian_Matrix.txt",header = T)
ggplot(data = melt(as.matrix(LagrangianMatrix)), aes(Var1,Var2, fill = value)) + 
  geom_tile(color = "white",names ="Travel Time (days)") + scale_fill_viridis() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=90,vjust = 0.8 ), 
        axis.text.y = element_text(vjust = 0.8 )) +
  xlab("") + ylab("")

###Linear Mix Model & plot
LinearMixedModelMVS = function(MVSList,LagMatrix,EnvList) {
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  Results = list()
  for (i in 1:length(MVSList)) { 
    Fst_matrix = MVSList[[i]]@pwFst
    subpoplist = MVSList[[i]]@pop
    MVSname = MVSList[[i]]@name
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
    
    Results[[i]] = list(MVSname = paste(MVSname,sep="_"), Model = Model, Data = X, 
                        MME=table2,TotalVariance = TotalVariance)
    print(table2)
    
  } 
  return(Results) }
PlotLMM = function(Results) {
    print(ggplot(data= Results$MME[-c(1,2),], 
    aes(x=rownames(Results$MME[-c(1,2),]),
        y=percentage)) + 
    geom_bar(stat = "identity",fill="goldenrod",col="black") + 
    theme_bw() + ggtitle(Results$MVSname) + xlab("Parameters") + 
      ylab("Proportion of explained variation") + ylim(values=c(0,1))) 
}


VariancePartition.Results = LinearMixedModelMVS(MVS.List,LagrangianMatrix,envlist)

for (i in 1:4) {
PlotLMM(VariancePartition.Results[[i]])
}



