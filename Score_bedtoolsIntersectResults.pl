# !/usr/bin/perl
# Score_GLMnetLasso.pl by _CST_ 150719 >>> Score.bedtoolsIntersectResults.pl by CST 170618
use strict; use warnings;

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

die "\nUsage: Score_GLMnetLasso.pl by _CST_ 150719\nThis program takes input from GREGOR, GREGORmine.pl and GLMnet_Lasso (w coefficients) and outputs a new tab delimited file of SNP..Index..Score\n\nInput requires a GREGOR_Mined file with rows of snps and columns of features in a 0/1 pattern describing nonoverlapping/overlapping features.\nInput also requires a file with Feature name and Coefficient/Multiplier in tab delimited format\nOutput is a file that contains a column of snps next to a column of Scores for those snps.\n\nInput perl <path to Score_GLMnetLasso.pl> <GREGORmine.pl results file path> <Input File for Coefficients> <Column # to find Feature Name> <Col# to find Coef value> <Output txt file name for 2 column snp..Index..Score output>\n\nNote that there is no delineation between Index/Proxy here, as opposed to Score_GLMnetLasso.pl\n\n\n" unless @ARGV==5;

# Input files
my $m = $ARGV[0];
open(my $mined, "<", "$m") or die "cannot open $ARGV[0]\n"; #Mined results file - a table with rows as snps and columns 0/1 for overlap w features

my $C = $ARGV[1];
open(my $coeffile, "<", "$C") or die "cannot open $ARGV[1]\n"; #File of coefficients, for example 160222.LassoCoefs.txt

my $Featcolnum = $ARGV[2]-1;
my $Coefcolnum = $ARGV[3]-1;
#my $Indexcolnum = $ARGV[4]-1; #not applicable for this

#print "Within file $ARGV[1]\nyour features are in column $ARGV[2]\n and coefficients are in column $ARGV[3]\n\nwithin file $ARGV[0]\n your Index SNPs 0/1 is in column $ARGV[4]\n\n\n";

print "Within file $ARGV[1]\nyour features are in column $ARGV[2]\n and coefficients are in column $ARGV[3]\n\n\n";

# Output files 
my $o = $ARGV[4];
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

#print OUT "SNP\tIndex?\tScore\n";
print OUT "SNP\tScore\n"; #if don't want Index/nonIndex column

while (my $line = <$mined>) { #only want to do this for first line (column name header)

    my @array = split('\s+', $line);
#    print "\n\n@array\n\n";

    if ($line =~ /^[SNP|snp]/) {
	chomp($line);

	for (my $i=1; $i<scalar(@array); $i++) {
#	    print "Column $i is $feathash{$array[$i]}\n";

	    if (defined $feathash{$array[$i]}) {

		$feathash{$array[$i]}{'Col Number'} = $i; #column names linked with column number so I can find them in each row
#		print "Column $i $array[$i] and is being scored as _($feathash{$array[$i]})_\n"; #if you want, remove # before to give the column numbers you are scoring on
	    }
	}
    }

    if ($line =~ /^[0-9|chr|rs]/) {
	chomp($line);
#	print "found one\n/n";
	print OUT "$array[0]\t"; #this is snp name

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
