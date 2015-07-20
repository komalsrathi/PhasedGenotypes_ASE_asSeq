library(reshape2)
library(ggplot2)

####################################### get allele specific counts ####################################### 
######################### counts #########################
setwd("~/Desktop/komalrclust/ngs/MAGnet_RNAseq/PhasedGenotypes/chr19_ase/ase_counts/")
samples = read.csv('../../../bin/sample_phenoData.csv')

dat.list = list.files(pattern="*.counts")
dat = do.call("cbind", lapply(dat.list, function(fn) data.frame(Filename=fn, read.table(fn,stringsAsFactors=FALSE))))
row.names(dat) = dat$V1
dat = dat[,grep("V2", colnames(dat))]
names(dat) = paste(sub("^([^_.]*).*", "\\1", dat.list))
dat = dat[grep('^ENSG',rownames(dat)),]
######################### counts #########################

######################################### get annotation ###################################
ann <- read.delim('~/Desktop/komalrclust/ngs/UNC_HBE/gtf_files/Homo_sapiens.GRCh37.75.protein_linc_nonCodingNovels_geneid_genename.txt',header=F)
ann <- ann[grep('^ENSG',ann$V1),]
######################################### get annotation ###################################

# dat.t = dat[grep('ENSG00000170873|ENSG00000249816',rownames(dat)),] # MTSS1 & LINC00964
# dat.t = dat[grep('ENSG00000167874|ENSG00000173917|ENSG00000120093|ENSG00000188716',rownames(dat)),] # TMEM88, HOXB2, HOXB3, UPD
dat.t = dat[grep('ENSG00000130816|ENSG00000022556|ENSG00000269699|ENSG00000198300|ENSG00000268654',rownames(dat)),]
dat.tt = as.data.frame(t(dat.t))
dat.tt$hap = rep(x=c('hap1','hap2'),64)
dat.tt$sample = colnames(dat)


dat.mtss1 = dcast(data = dat.tt,formula = sample~hap,value.var = 'ENSG00000170873')
dat.mtss1$ratio = dat.mtss1$hap1/(dat.mtss1$hap1+dat.mtss1$hap2)
dat.mtss1 = merge(dat.mtss1,samples[,c(1,2,4)],by.x='sample',by.y='sample_name')
dat.mtss1$ratio[is.na(dat.mtss1$ratio)] <- 0

dat.lnc00964 = dcast(data = dat.tt,formula = sample~hap,value.var = 'ENSG00000249816')
dat.lnc00964$ratio = dat.lnc00964$hap1/(dat.lnc00964$hap1+dat.lnc00964$hap2)
dat.lnc00964 = merge(dat.lnc00964,samples[,c(1,2,4)],by.x='sample',by.y='sample_name')
dat.lnc00964$ratio[is.na(dat.lnc00964$ratio)] <- 0


write.csv(dat.lnc00964,'../../LINC00964_ASE.csv',quote=F,row.names=F)
write.csv(dat.mtss1,'../../MTSS1_ASE.csv',quote=F,row.names=F)
####################################### get allele specific counts ####################################### 

# position 5831 snp of interest
dt <- (getsnpList.R)
############################## make dot plot ################################################

snpfiles <- list.files(path = '~/Desktop/komalrclust/ngs/MAGnet_RNAseq/PhasedGenotypes/hetSNP/snpList_files/',pattern="*.txt",full.names = T)

vector.of.homozygotes = c('')
for(i in 1:64)
{
  snpfiles.dt = read.delim(snpfiles[[i]],header=F)
  if(length(grep('125857359',snpfiles.dt$V2))==0)
  {
    names = sub('.*_','',snpfiles[[i]])
    names = sub('.txt','',names)
    vector.of.homozygotes = c(vector.of.homozygotes,names)
  }  
}

vector.of.homozygotes = vector.of.homozygotes[2:length(vector.of.homozygotes)]

dat.mtss1$heterozygote = ifelse(dat.mtss1$sample %in% vector.of.homozygotes,'no','yes')
dat.lnc00964$heterozygote = ifelse(dat.lnc00964$sample %in% vector.of.homozygotes,'no','yes')

dat.mtss1 = merge(dat.mtss1,dt,by='sample')
dat.mtss1$allele = paste(dat.mtss1$hap1.allele,dat.mtss1$hap2.allele,sep='')
dat.mtss1$allele[dat.mtss1$allele=='TG'] <- 'GT'

dat.lnc00964 = merge(dat.lnc00964,dt,by='sample')
dat.lnc00964$allele = paste(dat.lnc00964$hap1.allele,dat.lnc00964$hap2.allele,sep='')
dat.lnc00964$allele[dat.lnc00964$allele=='TG'] <- 'GT'

ggplot(dat.mtss1, aes(y=ratio,x=allele,color=status)) + geom_point(cex=3) + ggtitle('ASE MTSS1 (rs12541595)\n') +
  ylab('hap1/hap1+hap2\n') + xlab('') + theme(axis.text.x =element_text(colour = 'black',size=14),
                                              axis.text.y=element_text(colour="black",size=14),
                                              title=element_text(size=16),
                                              legend.text = element_text(color='black',size=14),
                                              legend.title = element_text(color='black',size=14))

ggplot(dat.lnc00964,aes(y=ratio,x=allele,color=status)) + geom_point(cex=3) + ggtitle('ASE LINC00964 (rs12541595)\n') +
  ylab('hap1/hap1+hap2\n') + xlab('') + theme(axis.text.x =element_text(colour = 'black',size=14),
                                              axis.text.y=element_text(colour="black",size=14),
                                              title=element_text(size=16),
                                              title=element_text(size=16),
                                              legend.text = element_text(color='black',size=14),
                                              legend.title = element_text(color='black',size=14)) 


ggplot(dat.lnc00964,aes(y=ratio,x=status,fill=status,color=status)) + 
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5, position='dodge',show_guide=F) + facet_grid(~allele) + ggtitle('ASE LINC00964 (rs12541595)\n') +
  ylab('hap1/hap1+hap2\n') + xlab('') + theme(axis.text.x =element_text(colour = 'black',size=14),
                                              axis.text.y=element_text(colour="black",size=14),
                                              title=element_text(size=16),
                                              title=element_text(size=16),
                                              legend.text = element_text(color='black',size=14),
                                              legend.title = element_text(color='black',size=14),
                                              strip.text = element_text(color='black',size=14)) 

ggplot(dat.mtss1,aes(y=ratio,x=status,fill=status,color=status)) + 
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5, position='dodge',show_guide=F) + facet_grid(~allele) + ggtitle('ASE MTSS1 (rs12541595)\n') +
  ylab('hap1/hap1+hap2\n') + xlab('') + theme(axis.text.x =element_text(colour = 'black',size=14),
                                              axis.text.y=element_text(colour="black",size=14),
                                              title=element_text(size=16),
                                              title=element_text(size=16),
                                              legend.text = element_text(color='black',size=14),
                                              legend.title = element_text(color='black',size=14),
                                              strip.text = element_text(color='black',size=14)) 


############################## make dot plot ################################################



######################################### Chi Square Test ###############################################
library(plyr)
dat.t = dat
colnames(dat.t) = sub('.counts','',dat.list)
dat.t$Gene = rownames(dat)
dat.t = melt(dat.t,id.vars = 'Gene')
dat.t$variable = as.character(dat.t$variable)
dat.t = dat.t[order(dat.t$Gene,dat.t$variable),]
dat.t$hap = sub('.*_','',dat.t$variable)
dat.t$variable = sub('_.*','',dat.t$variable)
dat.tt = dcast(data = dat.t, formula = Gene+variable~hap,value.var = 'value')

# chisq test
dat.ttt = dat.tt[which(rowSums(dat.tt[3:4])>=20),]
dat.ttt = cbind(dat.ttt,t(apply(dat.ttt[,c(3:4)],1,function(x) with(chisq.test(x[1:2]),c(statistic,p.value=p.value)))))
colnames(dat.ttt)[c(2,5)] = c('sample_name','X.squared')
dat.ttt$p.adjust <- p.adjust(dat.ttt$p.value,method="fdr")

# summation of chisq pvalues
chisq.sum <- function(arg1)
{
  X.squared.sum <- sum(arg1$X.squared)
  n  <- length(unique(arg1$sample_name))
  if(n>1)
  {
    p.value <- pchisq(q=X.squared.sum,df=n-1,lower.tail=FALSE)
    dt <- rbind(data.frame(n,X.squared.sum,p.value))
    return(dt)
  }
}

results <- ddply(.data=dat.ttt,.variables="Gene",.fun=chisq.sum)
results <- results[order(results$X.squared.sum,decreasing=T),]
results.sig <- results[which(results$p.value<0.05),]
write.csv(results,'~/Desktop/komalrclust/ngs/MAGnet_RNAseq/PhasedGenotypes/chr19_ase/chr19_ASE_chisq.csv',quote=F,row.names=F)
######################################### Chi Square Test ###############################################



######################################### get significant individuals ###################################
count.fun <- function(x)
{
  x.tot = nrow(x)
  x.sig = nrow(x[which(x$p.adjust<0.05),])
  dat = data.frame(total=x.tot,sig=x.sig,perc.sig=x.sig/x.tot*100)
  return(dat)
}
count.res = ddply(.data=dat.ttt,.variables="Gene",.fun=count.fun)
count.res = merge(count.res,ann,by.x='Gene',by.y='V1')
count.res = count.res[order(count.res$perc.sig,count.res$sig,decreasing = T),]
######################################### get significant individuals ###################################


