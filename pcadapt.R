library(pcadapt)
library(ggplot2)
library(ggrepel)

###Produce 3 plots : PCA, screeplot, qqplot
###Results : list with 1°) pcadapt raw output 2°) Qvalues for all variants 3°) Top variants
PcadaptAnalysis <- function(FreqMatrix, ColVector,poplist,threshold=0.1) {
  x <- t(as.matrix(FreqMatrix))
  filename <- read.pcadapt(x,type="pool")
  y <- pcadapt(filename,min.maf= 0.05)
  scores = data.frame(y$scores)
  rownames(scores)= poplist
  scores$col = ColVector[rownames(scores)]
  PEV = paste(round((y$singular.values^2)*100,digits=1),"%",sep="")
  
  pdf("Structure.pdf")
  for (j in 2:(length(colnames(FreqMatrix))-1)) {
    print(ggplot(data = scores, aes(x=scores[,j-1],y=scores[,j])) + theme_bw() + 
            theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
            ggtitle("Principal Component Analysis") + geom_label_repel(color= scores$col,label=rownames(scores)) +
            xlab(paste(paste("PC",j-1,sep=""),"(",PEV[j-1],")",sep="")) + ylab(paste(paste("PC",j,sep="")," (",PEV[j],")",sep="")) +
            theme(legend.position="none") + geom_vline(xintercept = 0, lty = 3 ) + geom_hline(yintercept = 0,lty = 3) )
  }
  dev.off()
  pdf("Pcadapt_Screeplot.pdf")
  plot(y,option="screeplot")
  dev.off()
  pdf("Pcadapt_QQplot.pdf")
  print(plot(y,option="qqplot"))
  dev.off()
  
  PCAdapt_pvalue <- y$pvalue
  PCAdapt_qvalue <- p.adjust(PCAdapt_pvalue,method="BH")
  Selection = data.frame(rownames(FreqMatrix),PCAdapt_qvalue)
  Selection = Selection[Selection[,2] < threshold,]
  Selection= na.omit(Selection)
  colnames(Selection)= c("varID","q-value")
  
  Res = list(y,PCAdapt_qvalue,Selection)
  return(Res)
  write.table(Selection,"Pcadapt_Results.txt",quote=F,row.names = F)
}


ColVector = c("155" = "dodgerblue","158" = "dodgerblue","178" =" forestgreen",
              "206" = "orange","208" = "orange","209" = "orange","210" = "black")
  
poplist= c("155","158","178","206","208","209","210")  

PcadaptAnalysis(Freqmatrix,ColVector,poplist)

