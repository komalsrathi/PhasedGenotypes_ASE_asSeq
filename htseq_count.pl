#!/usr/bin/perl -w
use Parallel::ForkManager;
use Data::Dumper;

#parameter file name as commandline argument
my $paramfile = $ARGV[0];

#read in param file, PARAM is the file handle, and create a hash of all the parameters.
#split by => and make the first column, the key and the second column, the value. 
open PARAM, $paramfile or die print $!;
my %param;
while(<PARAM>)
{
        chomp;
        @r = split('=>');
        #print "@r";
        $param{$r[0]}=$r[1];
}

#open file with filename which is value of the key FASTALIST
open FILE, $param{'FASTALIST'} or die

#create an array samples
my @samples;

#splitting each line based on ',' and storing in an array @r
#pushing the reference of this array in another array @samples
while(<FILE>){
        chomp;
        my @r = split(',');
        push(@samples,\@r);
}

#run parallel jobs
my $pm=new Parallel::ForkManager($param{'JOBS'});

# for loop to go through each chromosome
foreach my $c (1..16,18,20..22)
{
    $pm->start and next;

    # assume PROJECTNAME directory is PhasedGenotype dir
    my $chr = "chr$c";
    $chr = $chr.'_ase';
    $param{'INPUTDIR'} = $param{'PROJECTNAME'}.'/'.$chr.'/extractAsReads_output';
    
    # output directory for counts
    $param{'OUTPUTDIR'} = $param{'PROJECTNAME'}.'/'.$chr.'/ase_counts';
    system("mkdir $param{'OUTPUTDIR'}") unless (-d $param{'OUTPUTDIR'});
    
    # new fork manager
    my $pfm=new Parallel::ForkManager($param{'THREADS'});
    
    # for loop to go through each sample
    foreach (@samples)
    {	
        $pfm->start and next;
        
            my $htseqcounthap1 = "samtools view -f 0x0002 $param{'INPUTDIR'}/".$_->[2]."_hap1.bam | awk '!/\\t\\*\\t/' - | htseq-count -q -s reverse - $param{'GTF'} > $param{'OUTPUTDIR'}/".$_->[2]."_hap1.counts";
            print $htseqcounthap1,"\n";
            system($htseqcounthap1);
	
            my $htseqcounthap2 = "samtools view -f 0x0002 $param{'INPUTDIR'}/".$_->[2]."_hap2.bam | awk '!/\\t\\*\\t/' - | htseq-count -q -s reverse - $param{'GTF'} > $param{'OUTPUTDIR'}/".$_->[2]."_hap2.counts";
            print $htseqcounthap2,"\n";
            system($htseqcounthap2);
            
            print "$_->[2] completed\n";

        $pfm->finish;
    }
    
    print "chr$c completed\n";
    
    $pm->finish;
}
