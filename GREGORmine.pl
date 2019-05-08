#!/usr/bin/perl
# GREGORmine.pl by _CST_ started on 150508 May 8, 2015
use strict; use warnings;

###!!Note that bedtools IntersectBed can do this much more quickly and efficiently than my script!! I wish I knew about bedtools before I made this...##

#--- Purpose: ----------------------------------------------------------------------------------------------------------------------
#We would like to use statistical modeling to identify chromatin marks/genomic features that contribute to a given snp or genomic locus being important for platelet development.
#We know that several snps met genome wide significance in Gieger et al - the question is what features may have contributed to this.
#This program will mine through GREGOR output data to put together a big tab-delimited spreadsheet that has the following columns:
# snp(chr:pos)	MinorAlleleFrequency	DistanceToNearestGene	#LDbuddies(with r2>0.8 as specified in GREGOR)	mark	mark	mark   etc...
#Note that the "marks" will be either 0 or 1 to say whether that snp was found to overlap with a given mark
#This data can then be used easily for statistical applications and modeling in R (or anywhere else)

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

##########################################
##### Specifics: #####

# This program takes the results of a GREGOR run (http://genome.sph.umich.edu/wiki/GREGOR) and mines it for data on Index SNPs and the proxy SNPs identified by GREGOR. It prints out those information into a file called Mined_output.txt that will reside within your GREGOR_results directory.
# The output from this program can be used for regression modeling or other purposes for which having data on which snps reside in which genomic features is useful. 

#The actual output is a table - each row is a SNP and the columns are printed at the top of the table. The columns are MAF, Distance to nearest gene, # of LD buddies, whether or not it was an Index SNP or Neighbor SNP (Index SNP values are '1', Neighbors are '0'), and all of the chromatin marks that you fed into GREGOR. '1'  means the snp overlaps a given signal or chromatin mark, '0' means the snp is not found overlapping that mark. 

#####################################################################################################
#####################################################################################################
#####################################################################################################

#Usage statement
die "usage: GREGORmine.pl by _CST_ 150508. Navigate to the GREGOR_results directory and Insert the following on the command line:\nperl <Path to GREGORmine.pl> <GREGOR results directory> <Output file> <Example _short.bed file directory>\n Make sure results directory has a '/' at the end.\n Make sure all the files you want to analyze end in _short.bed! And change the names of those _short.bed files to include all useful information like Mark, Feature, Cell Type, Genome Build, etc. prior to running GREGORmine.\n" unless @ARGV==3;


# Identify the GREGOR Results Directory
my $resultsdir;
if ($ARGV[0] =~ /\/$/) { #$ARGV[0] = GREGOR results directory, this if/else makes sure there's a backslash at end of directory for later
    $resultsdir = "$ARGV[0]";
}
else {$resultsdir = "$ARGV[0]/"}

opendir(DIR, $resultsdir) or die "cannot open $resultsdir!\n";

#Open a file for output
my $outputfile = $ARGV[1];
open(OUT, ">", "$resultsdir$outputfile") or die "cannot open your Output file\n"; #This will be opened within your GREGOR results directory, which is why it's important to navigate there prior to run!


#Mark time that program started running
my $ds1 = localtime();
print "GREGORmine.pl run started at $ds1\n";

###################################################################
###Beginning with the GREGOR_results directory, need to compile a list of directories that include results > [Mark]_short.bed > index.snp.on.chr##.in.bed.txt / neighbor.snp.in.bed.txt

#Directory is ALL of the short.bed files that describe whether a given snp was found overlapping with a given mark
#note that there is a *_short.bed file for each track (e.g. histone mark)

my %bighash = (); #where I will store all data points

my $goodtxtfiles; #will be values for $bedhash{$markbed} later (to ONLY include txt files with info about index or neighbor snps within bed)
my @idxtxtfilearray; #an array of the goodtxtfiles for index snps
my @nbrtxtfilearray; #an array of text files for neighbor snps

my $dircont = `ls -lh | awk -F ' ' ' {print(\$9)} '`; #calls Unix and commands to open file list (ls -lh) and print 9th column, which is the file name
my @dirlist = split '\n', $dircont;
my $markbed; #bed file for a given mark

my @bedlist; #list of the "good" bed files that I want to look at later

foreach $markbed (@dirlist) {

    unless ($markbed =~ /^index.SNP.and.LD.for.top/) { #there is always a file called index.SNP.and.LD.for.top.2.bed that should NOT be included
	if ($markbed =~ /bed$|Peak$/) { #name of each file that I want to look at should end in _short.bed
	    push @bedlist, $markbed; # my $scalar = scalar(@bedlist) #this works 
	}
    }
}

#Make a list of all the txt files within a each mark's _short.bed file (there is a txt file for each chromosome)
my $examplebed = $ARGV[2]; #this is used to find all of the text files within the bed directories
my $markbedcont = `ls -lh $resultsdir$examplebed/ | awk -F ' ' ' {print(\$9)} '`; #Need to input an example bed file here!! e.g. CTCF_K562_HG19_short.bed
my @txtfilelist = split '\n', $markbedcont;

foreach my $txtfile (@txtfilelist) {
    if ($txtfile =~ /index.snp.LD.on.chr\d+.in.bed.txt$/) {
	push (@idxtxtfilearray, $txtfile);
    }
    if ($txtfile =~ /neighbor.on.chr\d+.LDbuddy.in.bed.txt$/) {
	push (@nbrtxtfilearray, $txtfile);
    }
}

####At this point I have organized useful directories for use later to figure out whether index/neighbor snps were found within a given Mark's bed file 
##############################################################################

my $ds1a = localtime();
print "Directories organized at $ds1a\n";

##################################################################################
#### The other important Directory to parse is for basic annotation (snp, maf, dist, #ofLDbuddies)

# Identify and open the file within the GREGOR results directory found at Results/index_SNP/annotated.index.snp.txt
my $annotation = "$resultsdir/index_SNP/annotated.index.snp.txt"; #this should always be the path to the annotated.index.snp.txt file
open(my $annotationfile, "<$annotation") or die "cannot open annotation file\n";

# Want to extract the first 4 columns from the annotation file, which are snp (in chr:pos format) maf dist LDbuddies
	#Will want to specify this as a hash with snp as key and others as values

#hash variables, note that for some reason using plain text (i.e., "MAF") works but "$MAF" doesn't!
my $maf; #minor allele frequency
my $dist; #distance to nearest gene
my $LDbuddies; #the number of snps in LD with any given snp (query or proxy)
my $snp= keys %bighash; #everything will be referenced to the $snp, and the $snp is in chr:pos format
my @snparray; #need this for removing duplicates later
#my $indexpos; #position

while (my $line = <$annotationfile>) { #this file is located at GREGOR_results file/index_SNP/annotated.index.snp.txt

	if ($line =~ /^[0-9]/) {
	my @array = split("\t", $line);

	$snp = $array[0]; #this works

	$bighash{$snp}{'MAF'} = $array[1];
	$bighash{$snp}{'Dist'} = $array[2];
	$bighash{$snp}{'LDBuddies'} = $array[3];
	$bighash{$snp}{'IndexSNP'} = 1; #all of the snps in this file are Index SNPs
	push @snparray, $snp;

	my @indexsnpsplit = split("\:", $snp);
	$bighash{$snp}{'IndexPos'} = $indexsnpsplit[1]; #need to know index snp position later (do not need chromosome # but could get that here too)
	}
}

close($annotationfile);
closedir(DIR);

my $ds1c = localtime();
print "Annotation data parsed at $ds1c\n";

########### Basic annotation data is now parsed and organized into hash structure
####################################################################################

#Initialize all the values as 0 for Index SNPs
 
foreach $markbed (@bedlist) { #these seem to be in alphabetical order when done this way, but regardless should be in same order as above
	foreach $snp (sort keys %bighash) {
	    $bighash{$snp}{$markbed} = '0';
	}
}

#Now go through each bedfile and textfile and change 0s to 1s if the snp does indeed overlap with a chromatin mark

foreach $markbed (@bedlist) { #these seem to be in alphabetical order when done this way, but regardless should be in same order as above
    
    foreach my $txtfile (@idxtxtfilearray) { #first is CTCF_HG19_short.bed
	
	open(IN, "<", "$resultsdir$markbed/$txtfile") or die "cannot open $resultsdir$markbed/$txtfile\n"; #Now this is only getting opened once!
	
	while (my $line = <IN>) {		
		
	    my @array = split("\t", $line);	       

	    foreach $snp (sort keys %bighash) {

		if ($array[0] =~ /$snp/ and $array[2] =~ /$bighash{$snp}{'IndexPos'}/) { #Column 1 are snps in chr:pos and Column 3 is position only!		
		    $bighash{$snp}{$markbed} = '1';

#	    if (defined $bighash{$array[0]}) {

		#The line below is for stringent inclusion...The snp itself has to be within mark bed to count!
	#	if ($array[2] =~ /$bighash{$array[0]}{'IndexPos'}/) { #Column 1 are snps in chr:pos and Column 3 is position only!
	#	    $bighash{$snp}{$markbed} = '1';
		}		
		
	    }
	}
	close(IN);	       		
    }
}

my $dsinitialized = localtime();
print "Index SNPs initialized at $dsinitialized!\nWorking to print now...\n";


#foreach $snp (sort keys %bighash) {
#
#		#Use this first option on the next line below for less stringent inclusion - note that a snp will get a '1' even if it's just another snp in LD that overlaps with the mark
#		if ($array[0] =~ /$snp/) { #Column 1 are snps in chr:pos and Column 3 is position only!
#
#		#Use this second option on the line below for more stringent inclusion (i.e. snp itself has to be within mark bed to count!
#    if ($array[0] =~ /$snp/ and $array[2] =~ /$bighash{$snp}{'IndexPos'}/) { #Column 1 are snps in chr:pos and Column 3 is position only!		
#	$bighash{$snp}{$markbed} = '1';
#    }
#    
#}
#	}
#close(IN);	       	
#    }

#}

################# Print Index SNP part of the table #################

#Print column headers
print OUT "SNP\tDist\tIndexSNP\tLDBuddies\tMAF" or die "trying to print on closed filehandle Index SNP output\n";
foreach $markbed (@bedlist) { #these seem to be in alphabetical order when done this way, but should be in same order as done above for analysis
    print OUT "\t$markbed";
}
print OUT "\n";

foreach $snp (sort keys %bighash) {
    print OUT "$snp\t$bighash{$snp}{Dist}\t$bighash{$snp}{IndexSNP}\t$bighash{$snp}{LDBuddies}\t$bighash{$snp}{MAF}";
    
    foreach $markbed (@bedlist) {
	print OUT "\t$bighash{$snp}{$markbed}";
    }
    
    print OUT "\n";
}


#################### Index SNP table end ##################
#############################################################
my $ds2 = localtime();
print "Index SNP analysis complete at $ds2!\nWorking on neighbor_snp data now...\n";

#####################################################################################################
#####################################################################################################
########## The rest of the program analyzes NEIGHBOR snps ###################







##################################################################################
##### Now I will parse out basic information for NEIGHBOR snps

my $nbrsnp= keys %bighash; #everything will be referenced to the $nbrsnp, and the $nbrsnp is in chr:pos format
my $hashref; #will use when initializing the nbr hash refs

opendir(DIR, $resultsdir) or die "cannot open $resultsdir!\n";

my $nbranno = "$resultsdir/neighbor_SNP/index.snp.neighbors.txt"; #this should always be the path to the annotated.index.snp.txt file
open(my $nbrannotationfile, "<$nbranno") or die "cannot open nbr annotation file\n";

my @bignbrarray; #These are just the nbrsnps associated with Index SNPs. The cube contains MANY more nbr snps

while (my $line = <$nbrannotationfile>) { 

    if ($line =~ /^[0-9]/) {
	chomp($line);
	my @array = split("\t", $line);
	
	my $somenbrsnps = $array[5];
	
	$somenbrsnps =~ s/\|/\t/g;
	my @somenbrarray = split("\t", $somenbrsnps);

	foreach $nbrsnp (@somenbrarray) {
	    push (@bignbrarray, $nbrsnp) unless exists $bighash{$nbrsnp}; #won't include index snps as nbrs!
	}	
    }
}

# Remove duplicates from @bignbrarray.

#my @nbrsnparray = uniq(@bignbrarray); #There were ~4% nbr snps duplicated in a given run before I added this! ##Can't use this unless install List::More Utils, so instead...

my %seen;
@seen{@bignbrarray} = ();
my @nbrsnparray = sort keys %seen;



my $val = scalar(@nbrsnparray); #There are 180158 neighbor snps in the 150427K562 run, which is consistent with ~2500 per index even though I only wanted min 500 proxies per index. 811082 nbr snps if 73 Index w 10000 nbrs per Index (150507All..)
print "The total number of unique nbr snps in your run is: $val\n";

close($nbrannotationfile);

#my $nbrsnpfile = `mkdir Nbr.SNP.list`;
#open(OUTA, ">", "Nbr.SNP.list");
#print OUTA "keys %seen\n"; #the goal here is to have a list to put into bedtoolsIntersect, which might work way faster than this...
#exit();



my $cubecont = `ls -lh $resultsdir/cube | awk -F ' ' ' {print(\$9)} '`; #cube contains all of the basic maf/dist/ldnum data for the neighbor snps (and seemingly a ton more snps too)
my @cubelist = split '\n', $cubecont;

my $cubepos; #position on chromosome for each nbr snp

foreach my $cubefile (@cubelist) {
    if ($cubefile =~ /chr\d+.maf.dist.ldnum.txt$/) { #these are the files that contain the basic info for the nbr snps
    open (INFO, "<", "$resultsdir/cube/$cubefile") or die "cannot open $cubefile\n" ;

    while (my $line = <INFO>) {
	chomp($line);
	if ($line =~ /^[0-9]/) {
	    my @array = split("\t", $line);

	    $cubefile =~ /chr(\d+).maf.dist.ldnum.txt$/;
	    my $cubechr = $1;
	    my $cubepos = $array[0];
	    $nbrsnp = "$cubechr:$cubepos";
	    
	    $bighash{$nbrsnp}{'NbrPos'} = $cubepos; #need this later to ensure nbrsnp ITSELF is in bed, not just something in LD
	    $bighash{$nbrsnp}{'MAF'} = $array[1]; #$maf; #Note that this only works for me when the second key is "MAF", cannot be "$MAF"
	    $bighash{$nbrsnp}{'Dist'} = $array[2]; #dist;
	    $bighash{$nbrsnp}{'LDBuddies'} = $array[3]; #LDbuddies;
	    $bighash{$nbrsnp}{'IndexSNP'} = 0;

	}
    }
    close(INFO);    

    my $cubetime = localtime();
    print "finished $cubefile at $cubetime\n"; #this takes <5min to get through all the cube files for the 150427K562 GREGOR run

    }
}


########### Basic annotation data is now parsed and organized into hash structure
####################################################################################
my $ds3 = localtime();
print "nbr annotation organization and basic information compiled at $ds3.\nStarting nbrsnp initialization now.\n";

####################################################################################
#### Search and print markbed data for neighbor snps just the same as Index SNPs, with the slight change of @nbrtxtfilearray instead of @idxtxtfilearray
################################################                                                                                                                         
###print out the contents of the table for neighbor snps                                                                          

###########################################################                                                                                                              
############ Neighbor SNP table ##############################                                                                                             


#Initialize all the values as 0 for Neighbor SNPs

#################
###this is from Lorenz method                 

#initialize hash of hashes                                                                                                                                     
###note that I think this will work, but am foregoing for now (as of 180823) - AstlePLT used this on 8.23.18
#foreach my $markbed (@bedlist) {
#    print "found markbed\n";
#    foreach my $nbrsnp (@nbrsnparray) {
#	print "found nbrnsnp\n";
#        $hashref -> {$nbrsnp} -> {$markbed} =0;
#    }
#}

#foreach my $a (@bedfiles) {
 #   foreach my $n (@nbrsnps) {
  #      $hashref -> {$a} -> {$n} =0;
   # }
#}

my $ds4x = localtime();
print "maybe faster?\nNbrsnps initialized at $ds4x.\n"; #will this work faster? 

########################
#########################
 
foreach $markbed (@bedlist) { #these seem to be in alphabetical order when done this way, but regardless should be in same order as above
	foreach $nbrsnp (@nbrsnparray) {
	    $bighash{$nbrsnp}{$markbed} = '0';
	}
}

my $ds4 = localtime();
print "Nbrsnps initialized at $ds4.\n"; #initialization takes 2 min with 180K nbrs, 1 bedfiles (73 index, 150506Tijssen)

#Now go through each bedfile and textfile and change 0s to 1s if the snp does indeed overlap with a chromatin mark

foreach $markbed (@bedlist) { #these seem to be in alphabetical order when done this way, but regardless should be in same order as above
    
    foreach my $txtfile (@nbrtxtfilearray) { #first is CTCF_HG19_short.bed #will need a separate loop for neighbor snps using @nbrtxtfilearray
	
	open(IN, "<", "$resultsdir$markbed/$txtfile") or die "cannot open $txtfile\n";	     #Now this is only getting opened once!
	
	while (my $line = <IN>) {		
		
	    my @array = split("\t", $line);	       

	    if (defined $bighash{$array[0]}) {
		
		#The following line and if statement ensures strictly nbr snps that themselves overlap with the mark are given '1'
		if ($array[2] =~ /$bighash{$array[0]}{'NbrPos'}/) { # This ensures that the nbr snp ITSELF overlaps the mark, not just another snp that's in LD with the nbr snp!
		$bighash{$array[0]}{$markbed} = '1';
		}
	    }
	}
	close(IN);	       		
    }
}

my $ds5 = localtime();
print "Now actually starting to print nbrsnp table at $ds5.\n"; #from ds4 until ds5 takes 1 second with 180158 nbr snps and 1bedfile (73index, 150506 Tijssen)

foreach $nbrsnp (@nbrsnparray) {
    print OUT "$nbrsnp\t$bighash{$nbrsnp}{'Dist'}\t$bighash{$nbrsnp}{'IndexSNP'}\t$bighash{$nbrsnp}{'LDBuddies'}\t$bighash{$nbrsnp}{'MAF'}\t"; #Get error message for uninitialized values- bc I'm looking up @nbrsnparray rather than keys bighash?
   
    foreach $markbed (@bedlist) {
	print OUT "$bighash{$nbrsnp}{$markbed}\t";
    }
    
    print OUT "\n";
}

#################### Neighbor SNP table end ##################                                                                                       
#############################################################

closedir(DIR);
close(OUT);

my $ds6 = localtime();
print "GREGORmine.pl complete at $ds6!\nCheck your GREGOR_results directory for $outputfile to see your results.\n"; #printing the nbrsnp portion of the table takes 3 min with 73 index, 180158 nbrs and 1bedfile (Tijssen)

exit();


########## On 150604 CST removed the code from below, as all duplicated snps should already have been removed in the code above this ########
#################################################
#################################################


############### Old Last Step ##############
##### this should produce an empty file that would otherwise contain snps that are duplicated in your output ###
#### ie any Index SNPs that ended up as Neighbor snps, or any Neighbor snps that got repeated ###

# Need to ensure that no Index SNPs ended up as Neighbor SNPs (which could happen for a different Index SNP, since GREGOR doesn't care about stuff like that. GREGOR finds Neighbor SNPs for each Index SNP independently

my $ds7 = localtime();
print "On to the last step!\nEnsuring no Index SNPs ended up listed as Neighbor SNPs at $ds7\n";  

#Open Mined_output.txt file to read                                                   
open(IN, "<", "$outputfile");

#Open Output file 
open(OUT2, ">", "Duplicates_in_$outputfile");

# Identify any snps that were used more than once

my %count = (); #individual snp identifiers
#my $total = 0; #total number of times a snp showed up in outputfile

#count the number of times each snp shows up

while (my $line = <IN>) {

    if ($line =~ /^[0-9]/) {
	my @array = split("\t", $line);
	
	my $wtf = $array[0];
	
	if (exists $count{$wtf}) {$count{$wtf}++}
	else                     {$count{$wtf} = 1}
#	$total++; #running tally of times a snp shows up in outputfile
    }
}

#report any snps that were used more than once

print OUT2 "Duplicated SNP ID\t#times used\n";

foreach my $wtf (sort keys %count) {
    if ($count{$wtf} > 1) {
    printf OUT2 "%s\t%d\n", $wtf, $count{$wtf};
    }
}

exit();
