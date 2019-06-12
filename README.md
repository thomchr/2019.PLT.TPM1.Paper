# 2019.PLT.TPM1.Paper

Enclosed are DataFiles and Scripts associated with Thom et al, Machine learning-based identification and cellular validation of Tropomyosin 1 as a genetic inhibitor of hematopoiesis. 

Please direct questions to thomc@email.chop.edu


# Overview
This workflow is designed to take sentinel GWAS SNPs (or any 'positive control' loci) and a set of chromatin features, create a table of overlap information and ultimately output a machine learning prediction model including features that accurately discriminate your original GWAS SNPs from a set of controls. 

There are too many chromatin features to replicate this entire workflow from the paper here, so a DEMO is provided that can walk you through the computational pipeline. 

# You will need the following software packages to run the small perl scripts created in this Github repo
- R (https://www.r-project.org/)
- Perl (https://www.perl.org/get.html)
- GREGOR (https://genome.sph.umich.edu/wiki/GREGOR)
- Glmnet (https://cran.r-project.org/web/packages/glmnet/index.html)

# Hardware requirements
This package should run on a 'normal' desktop computer

# OS requirements
These scripts were tested on Linux operating systems.

# Installation Guide
The Data Files branch of this repo has relevant files. Time to install should be <10min. All of the files used to generate data for this paper are publicly available from Geo Datasets (ENCODE, BROAD, Tijssen et al, Paul et al)

The Scripts branch of this repo contains fairly basic perl scripts to run the software above.

The DEMO below walks through the workflow to generate a LASSO model based on SNPs and chromatin features.

# DEMO
This will take you through the steps we used to create our LASSO-based prediction model. This DEMO uses 8 sentinel platelet trait GWAS SNPs (Gieger 2011) vs matched controls, choosing from 7 chromatin features. This is not a big enough data set to get a very good model, but will at least walk you through the practical steps we used:

#Identify SNP/feature overlaps and controls for each sentinel SNP
   1. Download Giegersnpsonly.txt #SNPs associated with platelet trait variation in Gieger et al 2011, GREGOR finds 79 of them
   2. Download the chromatin feature tracks and file 'DEMO.FeatureList'. You will need to specify correct directories in the DEMO.FeatureList where the chromatin features get downloaded 
   3. Download the template .conf file (you will need to change directories to make it work)
   4. Run GREGOR
      perl /project/voight_ML/lorenzk/GREGOR/script/GREGOR.pl --conf thomc.190601.DEMO.conf 
      (this should take <10min to run)

   The output from this GREGOR file needs to be 'Mined' to create a table documenting overlaps between SNPs and chromatin     features
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
  
-- this will probably give you an 'empty model' because of a convergence issue, or a bad model (see AUC pdf). This is because there are not enough features in this small DEMO, but this at least gives you a feel for practical aspects of how we created our model -- 

# To recreate data presented in the paper

#You can then Score SNPs genome-wide based on the LASSO features and coefficients
- Use Score.GLMnet.pl or Score_bedtoolsIntersectResults.pl to do this
- this is useful for accuracy estimates, etc. 

#Once you score genome-wide (or any subset you like), you can make ROC curves to estimate accuracy of your model. 
- If you download the relevant SNP score files from different models (LASSO, GWAVA, DeepSEA, CADD) and run with ROCR..R script, you will be able to reproduce data from our paper. 






