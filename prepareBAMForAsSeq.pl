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
foreach my $c (2,10,14,16)
{
    $pm->start and next;

    # assume PROJECTNAME directory is PhasedGenotype dir
    # create the chromosome directory if not already made
    my $chr = "chr$c";
    my $filename = $chr.'_ase';
    my $chrdir = $param{'PROJECTNAME'}.'/'.$filename;
    system("mkdir $chrdir") unless (-d $chrdir);
    
    # create extractAsReads directory in the chromosome directory
    my $outdir = $chrdir."/extractAsReads_input";
    system("mkdir $outdir") unless (-d $outdir);
    
    # new fork manager
    my $pfm=new Parallel::ForkManager($param{'THREADS'});
    
    # for loop to go through each sample
    foreach (@samples)
    {
        $pfm->start and next;

                # step 1:
                my $samids="samtools view $param{'INDIR'}/".$_->[2]."_filter_sort.bam $chr | cut -f1 | sort -k1,1 | uniq > $outdir/".$_->[2]."_IDs.txt";
                print $samids,"\n";
                system($samids);

                # step 2:
                my $grepids="samtools view $param{'INDIR'}/".$_->[2]."_filter_sort.bam | LC_ALL=C grep -w -F -f $outdir/".$_->[2]."_IDs.txt - > $outdir/".$_->[2]."_sample.sam";
                print $grepids,"\n";
                system($grepids);

                # step 3:
                my $getheader="samtools view -H $param{'INDIR'}/".$_->[2]."_filter_sort.bam > $outdir/".$_->[2]."_new.sam";
                print $getheader,"\n";
                system($getheader);

                # step 4:
                my $catsams="cat $outdir/".$_->[2]."_new.sam $outdir/".$_->[2]."_sample.sam | awk '!/\\t\\*\\t/' - > $outdir/".$_->[2].".sam";
                print $catsams,"\n";
                system($catsams);

                # step 5
                my $samtobam="samtools view -bS $outdir/".$_->[2].".sam > $outdir/".$_->[2].".bam";
                print $samtobam,"\n";
                system($samtobam);

                # step 6
                my $bamsort="samtools sort -n $outdir/".$_->[2].".bam $outdir/".$_->[2].".sorted";
                print $bamsort,"\n";
                system($bamsort);
		
                # step 7 remove unwanted files
                my $removefiles="rm $outdir/".$_->[2].".bam $outdir/".$_->[2].".sam $outdir/".$_->[2]."_new.sam $outdir/".$_->[2]."_sample.sam $outdir/".$_->[2]."_IDs.txt";
                print $removefiles,"\n";
                system($removefiles);

                print "$_->[2] completed\n";

        $pfm->finish;
    }
    
    $pm->finish;
}
