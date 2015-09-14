############################################### get SNP lists for asSeq ################################################
library(data.table)
library(reshape2)

# for this script you only need one input that is a snplist.txt file
# one chr:position per line

# read in list of snps you want to look up haplotypes for
snplist <- read.table('snplist.txt',header=F,stringsAsFactors=F)

# make a list of chromosomes that are in the snplist
chrs <- colsplit(snplist$V1,':',c('chr','pos'))
chrs <- unique(paste('chr',chrs$chr,sep=''))

tmp <- data.frame()
# loop through all chromosomes in chrs
for(i in 1:length(chrs))
{
  setwd('/fujfs/d2/MAGnet_RNAseq/ASE/asSeq_analysis/')
  
  # make a directory
  chr <- chrs[i]
  dirname <- paste(chrs[i],'_ase',sep = '')
  
  # set path to that directory
  setwd(dir = dirname)
  
  # get the snp positions & names (here names=positions)
  anno1 <- fread(sprintf('/fujfs/d1/Genotype/Illumina/imputed/Illumina_OEE_FWD_%s.info',chr),header=T)
  anno1$coordinate = sub('.*:','',anno1$SNP)
  snp1 <- as.character(anno1$SNP)

  # get the number of samples
  ff1 <- sprintf("%s.1.hapLabel", chr)
  cmd = sprintf("wc -l %s", ff1)
  wc1 = system(cmd, intern=TRUE)
  wc1 = unlist(strsplit(wc1, split="\\s+"))[1]
  nn1 = as.numeric(wc1)/2
  
  # loop through all samples
  for(k in 1:nn1){
    message("  ", k, " ", date())
    ff2  = sprintf("%s.%d.%d.tmp", ff1, 2*k-1, 2*k)
    cmd1 = sprintf("sed -n '%d,%d p' %s > %s", 2*k-1, 2*k, ff1, ff2)
    system(cmd1)
    
    # read in haplotypes
    dat1 = scan(ff2, what=character())
    cmd1 = sprintf("rm %s", ff2)
    system(cmd1)
    
    # do some checking
    if(dat1[1] != dat1[4]){ stop("sample names do not match\n") }
    if(dat1[2] != "HAPLO1"){ 
      stop(sprintf("expect HAPLO1 but see %s\n", dat1[2])) 
    }
    if(dat1[5] != "HAPLO2"){ 
      stop(sprintf("expect HAPLO2 but see %s\n", dat1[5])) 
    }
    
    # get all haplotype 1 and 2 data 
    hap1 = unlist(strsplit(dat1[3], split=""))
    hap2 = unlist(strsplit(dat1[6], split="")) 
    
    # check the length of snps and haplotypes match or not
    if(length(hap1) != length(snp1)){
      stop("length of haplotype 1 does not match annotation\n")
    }
    if(length(hap2) != length(snp1)){
      stop("length of haplotype 2 does not match annotation\n")
    }
    
    # get positions
    pos <- which(anno1$SNP %in% snplist$V1)
    
    # get snps at that positions
    snps <- anno1[pos,SNP]
    
    # get haplotypes for all snps in snplist
    tmp <- rbind(tmp, data.frame(sample = dat1[4], hap1.allele = hap1[pos], hap2.allele = hap2[pos], SNP=snps))
  }
    tmp$sample = sub('.*_','',tmp$sample)
}
