#!/usr/bin/perl
# GLMnet_Lasso_Rcode_creation.pl by _CST_ 160118
use strict; use warnings;

# to make Rcode but swap values for outputs to avoid manual input each time

# makes Rcode templates to be run using R CMD BATCH [...].R


die "Usage: GLMnet_Lasso_Rcode_creation.pl by  _CST_ created 160118. Please enter:\nperl GLMnet_Lasso_Rcode_creation.pl> <date> <run #> <#features> <exact feature column numbers (i.e., '2,4,5,6:18') this does NOT need parentheses but does need to quotes to work>\n and output will be a file named Date.GLMnet_Lasso_Rcode_#col#.columns_run#run#.R" unless @ARGV==4;

my $date = $ARGV[0];
my $run = $ARGV[1];
my $colnum = $ARGV[2];
my $cols = $ARGV[3];
my $totalcols = 3 + $colnum;

my $output = "$date.GLMnet_Lasso_Rcode_$colnum.columns_run$run.R";

open(my $out, ">", "$output");

print $out "setwd(\"/project/voight_ML/thomc/thomc_results/DEMO/GREGOR-results.DEMO\")
Sys.sleep(2)
z <- read.table(file=\"DEMO.Mined.txt\", header=T)
Sys.sleep(15)
library(glmnet)
Sys.sleep(2)
sink(\"$date.GLMnet_Lasso_Rcode_$colnum.columns_run$run.txt\")
cat(\"\\nnrow: \")
nrow(z)
cat(\"ncol: \")
ncol(z)
x <- data.matrix(z[,c($cols)])
y <- data.matrix(z[,3])
pfactor = matrix(nrow = 1, ncol = $totalcols)
pfactor[1,] = 1
pfactor[1,c(1:3)] = 0
cvfit = cv.glmnet(x, y, penalty.factor = pfactor, family = \"binomial\", type.measure = \"auc\")
Sys.sleep(15)
cat(\"print(cvfit)\\n\")
print(cvfit)
cat(\"coef(cvfit, s=lambda.1se)\\n\")
coef(cvfit, s=\"lambda.1se\")
cat(\"coef(cvfit, s=lambda.min)\\n\")
coef(cvfit, s=\"lambda.min\")
cat(\"AUC saved as $date.$colnum.factors.run$run.AUCplot.pdf\\n\")
pdf(\"$date.$colnum.factors.run$run.AUCplot.pdf\")
plot(cvfit)
dev.off()
savehistory(file=\"Routput.Rhistory\")
save.image(\"script-output.RData\")
cat(\"Script completed\\n\\n\")
sink()
detach(\"package:glmnet\")
q()";

close($out);
exit();
