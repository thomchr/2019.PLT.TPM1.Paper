# 2019.PLT.TPM1.Paper

Enclosed are DataFiles and Scripts associated with Thom et al, Machine learning-based identification and cellular validation of Tropomyosin 1 as a genetic inhibitor of hematopoiesis. 

Please direct questions to thomc@email.chop.edu


# Overview
This workflow is designed to take sentinel GWAS SNPs (or any 'positive control' loci) and a set of chromatin features, create a table of overlap information and ultimately output a machine learning prediction model including features that accurately discriminate your original GWAS SNPs from a set of controls. 

There are too many chromatin features to replicate this entire workflow from the paper here, so a DEMO is provided that can walk you through the computational pipeline. 

# You will need the following software packages in addition to the small perl scripts created in this Github repo
- R (https://www.r-project.org/)
- Perl (https://www.perl.org/get.html)
- GREGOR (https://genome.sph.umich.edu/wiki/GREGOR)
- Glmnet (https://cran.r-project.org/web/packages/glmnet/index.html)

#############
##DEMO
Here are the steps to work through the DEMO, which will create a LASSO-based model for 85 sentinel SNPs (vs controls) choosing from 7 chromatin features:

#Identify SNP/feature overlaps and controls for each sentinel SNP
1. Download Giegersnpsonly.txt #SNPs associated with platelet trait variation in Gieger et al 2011
2. Download the chromatin feature tracks and file 'DEMO.FeatureList'. You will need to specify correct directories in the DEMO.FeatureList where the chromatin features get downloaded 
3. Download the template .conf file (you will need to change directories to make it work)
4. Run GREGOR
perl /your/path/to/GREGOR/script/GREGOR.pl --conf thomc.190601.DEMO.conf #this small demo should take <15min to run

The output from this GREGOR file needs to be 'Mined' to create a table documenting overlaps between SNPs and chromatin features
The script GREGORmine.pl will output this table for you...

#Generate a table compiling SNP/feature overlaps
1. Download GREGORmine.pl from Scripts
2. From within the GREGOR-results output from before, run the following (should take less than 10min for DEMO)
perl /project/voight_ML/thomc/thomc_scripts/thomc_perl_scripts/GREGORmine.pl /project/voight_ML/thomc/thomc_results/DEMO/GREGOR-results.DEMO/ DEMO.Mined.txt wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak

Output should be a table with feature names and SNP characteristics as headers, and 0/1 values describing overlaps. MAF, Dist to Nearest Gene and Number of LD SNPs are in there too

#Run a LASSO model
1. Download GLMnet_Lasso_Rcode_creation.DEMO.pl from Scripts (this will help create the LASSO script to run)
2. Create an R code to run GLMnet
perl /project/voight_ML/thomc/thomc_scripts/thomc_perl_scripts/GLMnet_Lasso_Rcode_creation.DEMO.pl 190601 1 7 '2,4:12'

    Usage: GLMnet_Lasso_Rcode_creation.DEMO.pl by  _CST_ created 160118. Please enter:
     perl GLMnet_Lasso_Rcode_creation.pl> <date> <run #> <#features> <exact feature column numbers (i.e., '2,4,5,6:18') this         does NOT need parentheses but does need to quotes to work>
     and output will be a file named Date.GLMnet_Lasso_Rcode_#col#.columns_run#run#.R
  
 ###################
################
3. Run the output R script to make a LASSO model (#look at input/output directories, you might want to change them!)
R CMD BATCH 190601.GLMnet_Lasso_Rcode_7.columns_run1.R

   Your output here will show the following (should take <5 min)
   a. A cvfit file that has lots of info in it
   b. An AUC curve showing model fits associated with number of features included (at each level of lambda)
   c. 

#You can then Score SNPs genome-wide based on the LASSO features and coefficients
- Use Score..pl to do this

#Then make ROC curves to estimate accuracy of your model. 
- If you download the relevant files and run with script here, you should be able to reproduce data from our paper. 






