# 2019.PLT.TPM1.Paper

Enclosed are DataFiles and Scripts associated with Thom et al, Machine learning-based identification and cellular validation of Tropomyosin 1 as a genetic inhibitor of hematopoiesis. 

Please direct questions to thomc@email.chop.edu


# Overview
This workflow is designed to take sentinel GWAS SNPs (or any 'positive control' loci) and a set of chromatin features, create a table of overlap information and ultimately output a machine learning prediction model including features that accurately discriminate your original GWAS SNPs from a set of controls. 

There are too many chromatin features to replicate this entire workflow from the paper here, so a DEMO is provided that can walk you through initial stages of our computational pipeline. You can also run LASSO on the actual data we used to generate our model with datasets included here. 

# You will need the following software packages to run the small perl scripts we created
- R (https://www.r-project.org/)
- Perl (https://www.perl.org/get.html)
- GREGOR (https://genome.sph.umich.edu/wiki/GREGOR)
- Glmnet (https://cran.r-project.org/web/packages/glmnet/index.html)

# Hardware requirements
This package should run on any 'normal' desktop computer

# OS requirements
These scripts were tested on Linux operating systems.

# Installation Guide
The Data Files branch of this repo has relevant files. Time to install should be <10min. All of the chromatin features tracks and files used to generate data for this paper are publicly available from Geo Datasets (ENCODE, BROAD, Tijssen et al, Paul et al)

The Scripts branch of this repo contains fairly basic perl scripts that utilize the software above.

# DEMO
There is a small DEMO that walks through the workflow to generate a LASSO model input based on 8 sentinel platelet trait GWAS SNPs (Gieger 2011) vs matched controls, choosing from 7 chromatin features. 
This is not a big enough data set to get a very good model! But it represents the practical steps we used to make ours.

# To recreate data presented in the paper

#To make a LASSO model using actual data sets and chromatin features overlaps (this takes hours to run)
	Download Giegersnps_Mined.txt
	Run LASSO.DEMO.R
		R CMD BATCH LASSO.DEMO.R

   Output will be analogous to: 
	Lasso_output.txt
	LASSO.output.AUCplot.pdf
 

#To score SNPs genome-wide
We ran Score_GLMnetLasso.pl to analyze the snp147 genome based on the features and coefficients in our model.
The output is too big for Github, but we are happy to provide you with this file if you’d like.
Relevant subsets of this scored genome that were used for subsequent analyses are included (e.g., as controls for Gieger73.LASSO.forROCR.txt)

#To create ROCR plots to estimate accuracy of those score predictions, 

1. For Training set, download the following Data Files: 
	Gieger73.LASSO.forROCR.txt
	Gieger73.GWAVA.forROCR.txt
	Gieger73.DeepSea.forROCR.txt
	Gieger73.CADD.forROCR.txt

	Download and run from Scripts: ROCR.submission.Gieger73.R

2. For Holdout set, download the following Data Files: 
	Gieger8.LASSO.forROCR.txt
	Gieger8.GWAVA.forROCR.txt
	Gieger8.DeepSea.forROCR.txt
	Gieger8.CADD.forROCR.txt

	Download and run from Scripts: ROCR.submission.Gieger8holdouts.R

3. For Validation set, download the following Data Files: 
	AstlePLT.LASSO.forROCR.txt
	AstlePLT.GWAVA.forROCR.txt
	AstlePLT.DeepSea.forROCR.txt
	AstlePLT.CADD.forROCR.txt

	Download and run from Scripts: ROCR.submission.AstlePLT_Giegerpaper.R

#To create plots that help visualize overlay gene locations and chromatin feature tracks/peaks, we used Gviz. Sample script that generated data found in the paper can be seen at: 
	PIK3CG.Gviz.Rcode.R
	TPM1.Gviz.Rcode.R
	
   We are happy to help but recommend you play around with this and your own chromatin feature tracks. There are also very nice tutorials for Gviz online. 

#To identify associations between model scores and GWAS pvalues, we used the following Script. Summary statistics used for these analyses are not publicly available (Gieger 2011 study) but you could certainly adapt this script for other related GWAS: 
	Pval.vs.Score.pl







