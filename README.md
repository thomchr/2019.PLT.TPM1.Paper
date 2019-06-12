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

# Hardware requirements
This package should run on a 'normal' desktop computer

# OS requirements
The package development version is tested on Linux operating systems.

# Installation Guide
The Data Files branch of this repo has relevant files. Time to install should be <10min. All of the files used to generate data for this paper are publicly available from Geo Datasets (ENCODE, BROAD, Tijssen et al, Paul et al)

The Scripts branch of this repo contains fairly basic perl scripts to run the software above.

The DEMO below walks through the workflow to generate a LASSO model based on SNPs and chromatin features.

# DEMO
This will create a LASSO-based model for 8 sentinel platelet trait GWAS SNPs (Gieger 2011) vs matched controls, choosing from 13 chromatin features:

#Identify SNP/feature overlaps and controls for each sentinel SNP
1. Download Giegersnpsonly.txt #SNPs associated with platelet trait variation in Gieger et al 2011, GREGOR finds 79 of them
2. Download the chromatin feature tracks and file 'DEMO.FeatureList'. You will need to specify correct directories in the DEMO.FeatureList where the chromatin features get downloaded 
3. Download the template .conf file (you will need to change directories to make it work)
4. Run GREGOR
perl /project/voight_ML/lorenzk/GREGOR/script/GREGOR.pl --conf thomc.190601.DEMO.conf 
(this should take <10min to run)

The output from this GREGOR file needs to be 'Mined' to create a table documenting overlaps between SNPs and chromatin features
The script GREGORmine.pl will output this table for you...

#Generate a table compiling SNP/feature overlaps
1. Download GREGORmine.pl from Scripts
2. From within the GREGOR-results output from before, run the following (should take less than 10min for DEMO)
perl /project/voight_ML/thomc/thomc_scripts/thomc_perl_scripts/GREGORmine.pl /project/voight_ML/thomc/thomc_results/DEMO/GREGOR-results.DEMO/ DEMO.Mined.txt wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak

Output will be a table DEMO.Mined.txt
   This is a table with feature names and SNP characteristics as headers; MAF, Dist to Nearest Gene, and # LD SNPs, as well as 0/1 values describing overlaps comprise the rest of the table.

#Run a LASSO model
1. Download GLMnet_Lasso_Rcode_creation.DEMO.pl from Scripts (this will help create the LASSO script to run)
2. Create an R code to run GLMnet LASSO
perl /project/voight_ML/thomc/thomc_results/DEMO/GLMnet_Lasso_Rcode_creation.DEMO.pl 190601 1 7 '2,4:12'

    Usage: GLMnet_Lasso_Rcode_creation.DEMO.pl by  _CST_ created 160118. Please enter:
     perl GLMnet_Lasso_Rcode_creation.pl> <date> <run #> <#features> <exact feature column numbers (i.e., '2,4:18') this         does NOT need parentheses but does need to quotes to work>
     and output will be a file named Date.GLMnet_Lasso_Rcode_#col#.columns_run#run#.R
  
3. Run the output R script to make a LASSO model (#look at input/output directories, you might want to change them!)
R CMD BATCH 190601.GLMnet_Lasso_Rcode_7.columns_run1.R

   Your output here will show the following (should take <5 min)
   a. A cvfit file that has lots of info in it
   b. An AUC curve showing model fits associated with number of features included (at each level of lambda)
   c. 

You've now created a LASSO model and identified relevant features with associated coefficients for Gieger et al 2011 platelet trait variation GWAS! We have found that there will be some variation in model-selected features and coefficients between runs. 

-- end of demo --

#You can then Score SNPs genome-wide based on the LASSO features and coefficients
- Use Score.GLMnet.pl or Score_bedtoolsIntersectResults.pl to do this
- this is useful for accuracy estimates, etc. 

#Once you score genome-wide (or any subset you like), you can make ROC curves to estimate accuracy of your model. 
- If you download the relevant SNP score files from different models (LASSO, GWAVA, DeepSEA, CADD) and run with ROCR..R script, you will be able to reproduce data from our paper. 






