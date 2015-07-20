############################################### get SNP lists for asSeq ################################################

# make a list of chromosomes 1 to 22
chrs <- paste('chr',1:22,sep='')

# loop through chr1 to 22
for(i in 1:length(chrs))
{
  setwd('/fujfs/d3/home/komalr/ngs/MAGnet_RNAseq/PhasedGenotypes')
  
  # make a directory
  chr <- chrs[i]
  cmd <- sprintf('mkdir %s_ase',chr)
  system(cmd)
  dirname <- unlist(strsplit(x = cmd,split=' '))[2]
  
  # set path to that directory
  dirpath <- paste('/fujfs/d3/home/komalr/ngs/MAGnet_RNAseq/PhasedGenotypes/',dirname,sep='')
  setwd(dir = dirpath)
  
  # get the snp positions & names (here names=positions)
  anno1 <- read.table(sprintf('/fujfs/d1/Genotype/Illumina/imputed/Illumina_OEE_FWD_%s.info',chr),header=T)
  anno1$coordinate = sub('.*:','',anno1$SNP)
  snp1 <- as.character(anno1$SNP)
  
  # sort the haplotype file by sample name
  cmd = sprintf('zcat /fujfs/d1/Genotype/Illumina/imputed/Illumina_OEE_FWD_%s.hapLabel.gz | sort -k1,1 > %s.hapLabel',chr,chr)
  system(cmd)
  
  # get the number of samples
  # split by sample name
  cmd = sprintf('grep -f /fujfs/d3/home/komalr/ngs/MAGnet_RNAseq/bin/sample_list.txt %s.hapLabel > %s.1.hapLabel',chr,chr)
  system(cmd)
  ff1 <- sprintf("%s.1.hapLabel", chr)
  cmd = sprintf("wc -l %s", ff1)
  wc1 = system(cmd, intern=TRUE)
  wc1 = unlist(strsplit(wc1, split="\\s+"))[1]
  nn1 = as.numeric(wc1)/2
  
  # make a snplist file folder
  system('mkdir snpList_files')
  
  # loop through all samples
  # tmp <- data.frame()
  for(k in 1:nn1){
    message("  ", k, " ", date())
    ff2  = sprintf("%s.%d.%d.tmp", ff1, 2*k-1, 2*k)
    cmd1 = sprintf("sed -n '%d,%d p' %s > %s", 2*k-1, 2*k, ff1, ff2)
    system(cmd1)
    
    # --------------------------------------------------------------------
    # read in haplotypes and do some checking
    # --------------------------------------------------------------------
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
    
    hap1 = unlist(strsplit(dat1[3], split=""))
    hap2 = unlist(strsplit(dat1[6], split="")) 
    
    if(length(hap1) != length(snp1)){
      stop("length of haplotype 1 does not match annotation\n")
    }
    if(length(hap2) != length(snp1)){
      stop("length of haplotype 2 does not match annotation\n")
    }
    
    # get haplotypes for all snps
    # tmp <- rbind(tmp, data.frame(sample = dat1[4], hap1.allele = hap1, hap2.allele = hap2))
    
    # define alleles
    alleles = c("A", "C", "G", "T")
    
    # --------------------------------------------------------------------
    # extract SNPs with heterozygous genotypes
    # --------------------------------------------------------------------
    wdiff  = which(hap1 != hap2 & hap1 %in% alleles & hap2 %in% alleles)
    
    hap1A  = hap1[wdiff]
    hap1B  = hap2[wdiff]
    
    pos1   = anno1$coordinate[wdiff]
    
    hetDat = data.frame(chr=rep(chr, length(pos1)), pos=pos1, hap1A, hap1B)
    hetDat = hetDat[order(pos1),]
    
    sam1 = dat1[1]
    sam1 = sub('.*_','',sam1)
    ff2 = sprintf("snpList_files/hetSNP_%s.txt", sam1)
    
    write.table(hetDat, file = ff2, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
#   tmp$SNP = rep(anno1$coordinate,64,each = T)
#   tmp$sample = sub('.*_','',tmp$sample)
#   tmp.file <- paste(chr,'_haplotypes.csv',sep='')
#   write.csv(tmp, file = tmp.file, row.names=F, quote=F)
#   cmd = sprintf('bzip2 %s_haplotypes.csv',chr)
#   system(cmd)
}
