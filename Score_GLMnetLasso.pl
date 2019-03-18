# !/usr/bin/perl
# Score_GLMnetLasso.pl by _CST_ 150719
use strict; #use warnings;

################################
# Purpose: This program uses results from GREGOR/GREGORmine.pl as well as GLMnet_Lasso and creates a score for each snp. This is useful for identifying snps that might be interesting to look into biologically (ie with genome editing) but which might not have been previously identified in GWAS as statistically significantly associated with a given trait. 
################################

##############################################
##### Specifics ##############################
#    The whole path to getting results is to run GREGOR using some set of tracks/marks with some set of snps (ie, whole HG19 snp set vs Gieger et al statistically significant snps) ... GREGOR will output its file which is convoluted but useable. 
#    The next step is to run GREGORmine.pl on those GREGOR results to mine out information and output a large spreadsheet txt file with a single snp in a row and then a 0 or 1 depending on whether that snp is found within a given track (Tracks/Marks are in columns)
#    Then in a separate workflow, use the GREGOR results as input into GLMnet_Lasso, an R program that gives weights (beta values) to each mark. 
#   The final step is carried out here. Using the GREGORmine.pl output text file, this program creates a new file with a snp\tScore. The Score is a multiplier based on the Lasso beta value assigned to each mark. If a snp is found within a mark, the score for that mark is 1*WEight. If it's not, the score is 0*Weight.
# You will need to have a file that specifies feature name and beta multiplier (aka coefficient, for example from Lasso)  
##############################################
##############################################

die "\nUsage: Score_GLMnetLasso.pl by _CST_ 150719\nThis program takes input from GREGOR, GREGORmine.pl and GLMnet_Lasso (w coefficients) and outputs a new tab delimited file of SNP..Index..Score\n\nInput requires a GREGOR_Mined file with rows of snps and columns of features in a 0/1 pattern describing nonoverlapping/overlapping features.\nInput also requires a file with Feature name and Coefficient/Multiplier in tab delimited format\nOutput is a file that contains a column of snps next to a column of Scores for those snps.\n\nInput perl <path to Score_GLMnetLasso.pl> <GREGORmine.pl results file path> <Input File for Coefficients> <Column # to find Feature Name> <Col# to find Coef value> <Column number delineating Index/Proxy in Mined file> <Output txt file name for 2 column snp..Index..Score output>\n\n\n" unless @ARGV==6;

# Input files
my $m = $ARGV[0];
open(my $mined, "<", "$m") or die "cannot open $ARGV[0]\n"; #GREGORMine.pl results file

my $C = $ARGV[1];
open(my $coeffile, "<", "$C") or die "cannot open $ARGV[1]\n"; #File of coefficients

my $Featcolnum = $ARGV[2]-1;
my $Coefcolnum = $ARGV[3]-1;
my $Indexcolnum = $ARGV[4]-1;

#print "Within file $ARGV[1]\nyour features are in column $ARGV[2]\n and coefficients are in column $ARGV[3]\n\nwithin file $ARGV[0]\n your Index SNPs 0/1 is in column $ARGV[4]\n\n\n";

print "Within file $ARGV[1]\nyour features are in column $ARGV[2]\n and coefficients are in column $ARGV[3]\n\n\n";

# Output files 
my $o = $ARGV[5];
open(OUT, ">", "$o");

# Go through Coef file and link feature names with coefficient values

my $feat;
my $coef;
my %feathash;
my @feat = keys %feathash;

print "You are scoring on the following Features (w Coefficients):\n";
print OUT "Scoring based on:\nFeature\tCoefficient\n";

while (my $line = <$coeffile>) {    
    my @array = split('\s+', $line);

    if ($array[$Coefcolnum] != 0) {
    $feat = $array[$Featcolnum];
    $feathash{$feat}{'coef'} = $array[$Coefcolnum]; #now useful features set as feature>coef in hash

    print "$feat\t$array[$Coefcolnum]\n"; #this is a feature that will be scored
    print OUT "$feat\t$array[$Coefcolnum]\n"; #this is a feature that will be scored
    }
}
close($coeffile);

print OUT "\n\n";

###################
# Associate feature names with column position in Mined file

my $score;

print OUT "SNP\tIndex\tScore\n";
#print OUT "SNP\tScore\n"; #if don't want Index/nonIndex column

while (my $line = <$mined>) { #only want to do this for first line (column name header)

    my @array = split('\s+', $line);

    if ($line =~ /^[SNP|snp]/) {
	chomp($line);

	for (my $i=0; $i<scalar(@array); $i++) { 
	    
	    if (defined $feathash{$array[$i]}) {

		$feathash{$array[$i]}{'Col Number'} = $i; #column names linked with column number so I can find them in each row
		print "Column $i is being scored\t ($feathash{$array[$i]})\n"; #if you want, remove # before to give the column numbers you are scoring on
	    }
	}
    }

    if ($line =~ /^[0-9|chr|rs]/) {
	chomp($line);
#	print "found one\n\n";
	print OUT "$array[0]\t"; #this is snp name
	print OUT "$array[$Indexcolnum]\t"; #this is 1 for Index SNPs, and 0 for all others --- turn this OFF if you don't want this column and change your Header printOUT parameters above

	$score = 0; #will be total score for snp when coefficients multiplied

	foreach $feat (sort keys %feathash) {

#	    print "Feature $feat is in column $feathash{$feat}{'Col Number'}\n"; #this works
#	    print "$feat\t$array[$feathash{$feat}{'Col Number'}]*$feathash{$feat}{'coef'}\n"; #works

	    $score = $score + $array[$feathash{$feat}{'Col Number'}]*$feathash{$feat}{'coef'};

#	    print "$feat\t$array[$feathash{$feat}{'Col Number'}] x $feathash{$feat}{'coef'}\n"; #this works
	}

	print OUT "$score\n";
#	print "$score\n";# score calculated as above using sum of: $array[feature column number] * feature coefficient. All steps above worked on test runs
    }
}
close(OUT);
print "Your results can be found in folder $ARGV[5]\n";
exit();



###############################
################################
################################

# Old notes from prior trials of this. Probably all junk!









#6Individual marks from glmnet 1.9-8
#	print "You are scoring on the following marks:\n$array[15]\n$array[17]\n$array[19]\n$array[44]\n$array[45]\n$array[46]\n";

#6 Pairwise marks from glmnet 2.0-2 on 150826
#	print "You are scoring on the following marks:\n$array[83]\n$array[84]\n$array[85]\n$array[550]\n$array[552]\n$array[557]\n$array[557]\n$array[605]\n$array[617]\n$array[621]\n$array[632]\n";
#    }

#    if ($line =~ /^[0-9]/) {

#### !! NEED TO EDIT THE FOLLOWING LINE FOR EACH NEW DATASET !! #################


                      #6 INDIVIDUAL 150826- h3k36, h3k27ac, h3k4me2, fli1, gata1, gata2 
#    my $score = $array[17]*7.052358e-01 + $array[44]*3.773137e-01 + $array[45]*7.118857e-01 ;# + $array[44]*6.901514e-01 + $array[45]*9.100699e-01 + $array[46]*3.633894e-01 ; #I am not including scoring for Dist, LDBuddies or MAF here but I could

                      #6 Pairwise scores 150827
#    my $score = $array[15]*7.676002e-01 + $array[51]*2.842274e+00 + $array[53]*2.319708e-01 + $array[54]*3.615937e-01 + $array[57]*4.659684e-01 + $array[58]*-3.817886e-01 ; 

## 160204 Coefficients from MasterRun1-9

#	my $score = $array[83]*0.261 + $array[84]*0.629302338 + $array[85]*0.319805273 + $array[550]*0.435105213 + $array[552]*0.122608963 + $array[557]*0.472358863 + $array[605]*0.182827048 + $array[617]*0.548542538 + $array[621]*0.117955013 + $array[632]*0.546909938 ;

#    print OUT "$array[0]\t$score\n"; #This creates a 2 column snp\tscore txt file as it goes along. To add an IndexSNP column, put in $array[2] to be printed out somewhere

#    }
#}

###################
#150619 Analysis- these marks and coefficients were NOT based on Strict analysis!!##

#[thomc@scisub GREGOR-results_150507Alltracks_Pruned]$ perl /project/voight_datasets/blood/thomc/thomc_perl_scripts/Score_GLMnetLasso.pl /project/voight_datasets/blood/thomc/GREGOR-results/GREGOR-results_150507Alltracks_Pruned/150611_Minedoutput_Strict_150507Alltracks_Pruned_10Knbr.txt 150507AllTracks_10Knbrs_Scored.txt
#You are scoring on the following marks:
#H3k27ac_K562_HG19_short.bed
#H3k36me3_K562_HG19_short.bed
#H3k9ac_K562_HG19_short.bed
#Hmm_4_Strong_Enhancer_K562_HG19_short.bed
#Hmm_7_Weak_Enhancer_Hepg2_HG19_short.bed
#MACS_FLI1_vs_rIgG_PrimaryMK_hg19_short.bed
#MACS_GATA2_vs_rIgG_PrimaryMK_hg19_short.bed
#MACS_SCL_vs_rIgG_PrimaryMK_hg19_short.bed


#150820 Analysis - these marks and coefficients are based on Strict analysis and can be used in further analyses! ####
#[thomc@red01 GREGORmine_Output]$ perl /project/voight_datasets/blood/thomc/thomc_perl_scripts/Score_GLMnetLasso.pl /project/voight_datasets/blood/thomc/thomc_results/GREGORmine_Output/150611_Minedoutput_Strict_150507Alltracks_Pruned_10Knbr.txt 150820_Scored_150507MinedOutputStrict.txt
#You are scoring on the following marks:
#H3k27ac_K562_HG19_short.bed
#H3k36me3_K562_HG19_short.bed
#H3k4me2_K562_HG19_short.bed
#MACS_FLI1_vs_rIgG_PrimaryMK_hg19_short.bed
#MACS_GATA1_vs_rIgG_PrimaryMK_hg19_short.bed
#MACS_GATA2_vs_rIgG_PrimaryMK_hg19_short.bed

## There are 182 snps with perfect scores now! Note that the snps used here are simply the Index snps and proxy snps that GREGOR identified in the 150507 AllTracks run!


# **Need to score the HaemGenMPV GREGOR results to get a Genome wide idea of scores...but first need to Mine those results


#150827 pairwise Analysis




############ The following statistics are for ################### 
#### Dropbox>150610_cvfit_auc_Coef
#### and
#### /project/voight_datasets/blood/thomc/GREGOR-results/GREGOR-results_150507Alltracks_Pruned/150611_Minedoutput_Strict_150507Alltracks_Pruned_10Knbr.txt

#    my $score = $array[15]*2.551812e-01 + $array[17]*9.360492e-01 + $array[23]*1.139106e-02 + $array[33]*3.370834e-01 + $array[38]*3.796785e-03 + $array[44]*5.403801e-01 + $array[46]*3.042230e-01 + $array[48]*4.082557e-01 ; #I am not including scoring for Dist, LDBuddies or MAF here but I could


