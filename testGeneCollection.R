if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")

library(GSEABase)
library(celldex)
library(singscore)
library(ggplot2)


data(sample.ExpressionSet) # from Biobase
egs <- GeneSet(sample.ExpressionSet[201:250,], setName="Sample")
egs

head(geneIds(egs))
details(egs)


ref <- BlueprintEncodeData()

probs = c(.1,.25,.33333333,.5,.6666666,.75,.9)
diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2,3,4,5)
types=unique(colData(ref)$label.main)
#TODO: add dependencies
message('get quantiles...')
q = lapply(types,function(y){
  q = apply(assays(ref)$logcounts,1,function(z) quantile(z,probs,na.rm=TRUE))
})

message('create all signatures...')
ntypes= length(types)
#signature= data.frame(ref@colData)
listIn= list()
rankData <- rankGenes(ref)
for (i in 1:ntypes) {
  
  for (diff in diff_vals) {
    for (j in 1:round(length(probs)/2+0.5)) {
      diff_genes = lapply(q,function(x) q[[i]][j,]>x[length(probs)-j,]+diff)
      output <- matrix(unlist(diff_genes), nrow = ntypes, byrow = TRUE)
      for (n in (ntypes-1):(ntypes-3)) {
        g = colSums(output)>=n
        if (sum(g)>7 & sum(g)<201){
          gs=GeneSet(rownames(ref)[g], 
                     setName=sprintf("%%%s%%j%g%%d%g%%n%g",toupper(colnames(ref)[i]),round(probs[j]*100),diff,ntypes-n,
        sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE))
          listIn<-c(listIn, gs)
          
          score <-  simpleScore(rankData, geneIds(gs), centerScore = TRUE)$TotalScore
          score.sing = matrix(NA,ncol(rankData),length(gs))
          for (i in 1:length(gs))
            score.sing[,i] <- simpleScore(rankData, geneIds(gs), centerScore = TRUE)$TotalScore
          score = t(score)
          rownames(score) = names(gs)
          colnames(score) = colnames(rankData)
          score[is.na(score.sing)] = 0
          
        }
        # write.table(rownames(ref)[g],file=file.path(working.dir,'signatures',
        #                                              sprintf("%%%s%%j%g%%d%g%%n%g.txt",toupper(colnames(ref)[i]),round(probs[j]*100),diff,ntypes-n)),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
      }
    }
  }
}
#gsc= GeneSetCollection(listIn)
# 
# score <-  simpleScore(rankData, geneIds(gsc[[i]]), centerScore = TRUE)$TotalScore
# score = t(score)

df = data.frame(CellType=ref$label.main, Score = t(score.sing[1,]))
ggplot(df,aes(x=CellType,y=score,fill=CellType))+geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
       


