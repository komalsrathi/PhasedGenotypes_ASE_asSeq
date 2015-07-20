library(plyr)
library(asSeq)
library(doMC)
registerDoMC(4)

chrs <- paste('chr',1:22,sep='')
chrs <- chrs[c(2,10,14,16)]

extractreads <- function(x)
{
  # reset directory
  setwd('/fujfs/d3/home/komalr/ngs/MAGnet_RNAseq/PhasedGenotypes/')
  
  # set directory to chromosome
  chr <- x
  setwd(paste(chr,'_ase',sep=''))
  
  # make output directory
  system('mkdir extractAsReads_output')
  
  # set path to that directory
  setwd(dir = 'extractAsReads_output')
  
  # get input files
  getsnpList <- list.files(path='../snpList_files',pattern='.txt',full.names=T)
  getbamList <- list.files(path = '../extractAsReads_input',pattern='.sorted.bam',full.names=T)
  for(i in 1:length(getsnpList))
  {
    input <- getbamList[[i]]
    snpList <- getsnpList[[i]]
    outputTag <- sub('.txt','',sub('../snpList_files/hetSNP_','',snpList))
    extractAsReads(input, snpList, outputTag = outputTag, flag2kp = 0, flag2rm = 1796,prop.cut = 0.5,
                   min.avgQ = 10, min.snpQ = 10, phred = 33, skip = 0)
    print(paste(outputTag,'completed',sep=' '))
  }
  system('rm *_hapN.bam')
}

# run function
l_ply(.data = chrs, .fun = extractreads,.parallel = TRUE)
