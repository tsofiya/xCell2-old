library(GSVA)
library(GSEABase)
library(pracma)
library(RColorBrewer)
library(pheatmap)

read.types.dependencies = function(file.name) {
  con  <- file(file.path(file.name), open = "r")
  out <- list()
  i=1
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(oneLine, "\t"))
    out$types[i] = vec[1]
    n = max(which(!(vec=="")))
    if (n<2)
      out$dep[i] = list(types="")
    else
      out$dep[i] = list(types=vec[2:n])
    i=i+1
  }
  close(con)
  out
}

create.signatures = function(ref,dependencies){
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
  signature= data.frame(ref@colData)
  listIn= list()
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

          }
          # write.table(rownames(ref)[g],file=file.path(working.dir,'signatures',
          #                                              sprintf("%%%s%%j%g%%d%g%%n%g.txt",toupper(colnames(ref)[i]),round(probs[j]*100),diff,ntypes-n)),sep="\t",row.names=FALSE,quote =FALSE,col.names = FALSE)
        }
      }
    }
  }
  gsc= GeneSetCollection(listIn)
  return gsc


}

score.ref.signatures = function(ref,working.dir) {
  dir.create(file.path(working.dir, 'scores'), showWarnings = FALSE)
  egc = read.signatures.dir(paste0(file.path(working.dir,'signatures')))
  scores = gsva(as.matrix(assays(ref)$logcounts),egc,method='ssgsea',ssgsea.norm=FALSE)
  write.table(scores,file=file.path(paste0(working.dir,"scores",paste(ref@colData$label.main,"_ssgsea.txt"))),sep="\t",row.names=TRUE,quote =FALSE,col.names = NA)
}

# TODO: delete this when replacing files with objects.
read.signatures.dir = function(path) {
  sig = list()
  files = list.files(path=path)
  for (i in 1:length(files)) {
    sig[i] = GeneSet(scan(paste0(path,'/',files[i]), what="", sep="\n"), setName=files[i])
  }
  signatures <- GeneSetCollection(sig)
  #toGmt(signatures,'~/Documents/signatures/Mouse/signatures.txt')
  signatures
}
