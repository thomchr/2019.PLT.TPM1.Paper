#!usr/bin/perl
# Pval.vs.ScorePlot.pl by _cst_ started 160925
use strict;use warnings;

# Purpose: The intent is to create an R plot comparing Pval and Score of a large series of snps. The utility will be to see if there's a correlation between pval of a snp in a GWAS and Score assigned to the given snp 

#Specifics:
# This uses a file that has SNP (as chr# pos) and Score: for example, GREGOR_..160218/160403_.. now moved to /project/voight_datasets/blood/thomc/thomc_results/Pval.vs.Score_Plot/160403.GenomeMPV.Scored.160222LassoCoefs.txt.Padded.sorted.cp
#This also uses a file that has SNP and pval: for example, from /project/voight_datasets/GWAS.. moved to /project/voight_datasets/blood/thomc/thomc_results/Pval.vs.Score_Plot/MPV.txt
# and file that links SNPs rsID and chr/pos: for example, from cthom@coruscant and ben's snp_lookup.py moved to /project/voight_datasets/blood/thomc/thomc_results/Pval.vs.Score_Plot/MPV_locations_all.txt

die "/nUsage: Pval.vs.ScorePlot.pl by _cst_ works to link SNPs with pval and Scores to create a plot comparing these two values.\n\nPlease enter:\nperl Pval.vs.ScorePlot.pl <SNP_Score file> <SNP_pval file> <Linker from snp.lookup.py> <output file> <other output file for missing data>\n\n" unless @ARGV==5;

my $scorefile1 = $ARGV[0];
my $pvalfile1 = $ARGV[1];
my $linkfile1 = $ARGV[2];
my $outfile = $ARGV[3];
my $outfile_2 = $ARGV[4];

#variables
my %snpscore = ();
my $chrpos;
my $score;
my @chrpos = keys %snpscore;

my %snppval = ();
my $rsid;
my $pval;

my %linkhash = ();
my @rsid = keys %linkhash;

#open and hash rsID with pval
open(my $scorefile, "<", "$scorefile1");

while (my $line = <$scorefile>) {
    chomp($line);
    my @array = split /\t/, $line;

    $chrpos = $array[0];
    push(@chrpos, $array[0]);

    $snpscore{$chrpos} = "$array[1]";
  #  print "$chrpos, $snpscore{$chrpos}\n";
}
close($scorefile);

print "Finished with Score hash\n\n\n";

#open and hash chr.pos with Score
open(my $pvalfile, "<", "$pvalfile1");

while (my $line = <$pvalfile>) {
    chomp($line);
    my @array = split("\t+", $line);

    $rsid = $array[0];

    $snppval{$rsid} = "$array[5]";
 #   print "$rsid, $snppval{$rsid}\n";
}
close($pvalfile); 

print "Finished w Pval hash\n\n";

#Open and print each linker line with associated pval and score in same line
open(my $linkfile, "<", "$linkfile1");

while (my $line = <$linkfile>) {
    chomp($line);
    my @array = split("\t+", $line);
    $rsid = $array[4];
    push(@rsid, $rsid);

    $linkhash{$rsid} = $array[0];
#    print "$rsid, $linkhash{$rsid}\n";
}
close($linkfile);

open(OUT, ">", "$outfile");
open(OUT2, ">", "$outfile_2");

print OUT "rsID\tChrPos\tPval\tScore\n";
print OUT2 "rsID\tChrPos\tPval\tScore\n";

foreach my $rsid (@rsid) {
    if (defined $snppval{$rsid} and $snpscore{$linkhash{$rsid}}) {
    print OUT "$rsid\t$linkhash{$rsid}\t$snppval{$rsid}\t$snpscore{$linkhash{$rsid}}\n"; #I think the problem has to do with genome version for Scored snps
    }

    else {print OUT2 "$rsid\t$linkhash{$rsid}\t$snppval{$rsid}\t$snpscore{$linkhash{$rsid}}\n";}
}
close(OUT);
close(OUT2);

exit();