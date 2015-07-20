setwd('~/Desktop/komalrclust/ngs/MAGnet_RNAseq/PhasedGenotypes/')
samples = read.csv('~/Desktop/komalrclust/ngs/MAGnet_RNAseq/phenoData_HumanHeart.csv')

count.fun <- function(x)
{
  x.tot = nrow(x)
  x.sig = nrow(x[which(x$p.adjust<0.05),])
  dat = data.frame(total=x.tot,sig=x.sig,perc.sig=x.sig/x.tot*100)
  return(dat)
}

tmp = data.frame()
for(i in 1:22)
{
  # get counts
  chr = paste('chr',i,'_ase',sep = '')
  setwd(paste('~/Desktop/komalrclust/ngs/MAGnet_RNAseq/PhasedGenotypes/',chr,'/ase_counts/',sep=''))
  dat.list = list.files(pattern="*.counts")
  dat = do.call("cbind", lapply(dat.list, function(fn) data.frame(Filename=fn, read.table(fn,stringsAsFactors=FALSE))))
  row.names(dat) = dat$V1
  dat = dat[,grep("V2", colnames(dat))]
  names(dat) = paste(sub("^([^_.]*).*", "\\1", dat.list))
  dat = dat[grep('^ENSG',rownames(dat)),]
  
  # manipulate 
  dat.t = dat
  colnames(dat.t) = sub('.counts','',dat.list)
  dat.t$Gene = rownames(dat)
  dat.t = melt(dat.t,id.vars = 'Gene')
  dat.t$variable = as.character(dat.t$variable)
  dat.t = dat.t[order(dat.t$Gene,dat.t$variable),]
  dat.t$hap = sub('.*_','',dat.t$variable)
  dat.t$variable = sub('_.*','',dat.t$variable)
  dat.tt = dcast(data = dat.t, formula = Gene+variable~hap,value.var = 'value')
  
  # chisq
  dat.ttt = dat.tt[which(rowSums(dat.tt[3:4])>=20),]
  dat.ttt = cbind(dat.ttt,t(apply(dat.ttt[,c(3:4)],1,function(x) with(chisq.test(x[1:2]),c(statistic,p.value=p.value)))))
  colnames(dat.ttt)[c(2,5)] = c('sample_name','X.squared')
  dat.ttt$p.adjust <- p.adjust(dat.ttt$p.value,method="fdr")
  
  # count significant ind.
  count.res = ddply(.data=dat.ttt,.variables="Gene",.fun=count.fun)
  count.res = merge(count.res,gtf,by.x="Gene",by.y='V5')
  count.res = count.res[order(count.res$perc.sig,count.res$sig,decreasing = T),]
  
  ch = paste('chr',i,sep='')
  count.res = count.res[which(count.res$V1==ch),]
  tmp = rbind(tmp,count.res)
}
colnames(tmp) = c('Ensembl','Total.Samples','Significant','Perc.Significant','Chr','Start','End','Strand','Gene','Biotype')

tmp = tmp[order(tmp$Perc.Significant,tmp$Significant,decreasing = T),]

