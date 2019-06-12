# 2019.PLT.TPM1.Paper

This is a small DEMO that walks you through the steps we used to create our penalized regression model
It does not make a very good model, because there aren't many chromatin features!

Please direct questions to thomc@email.chop.edu

#Identify SNP/feature overlaps and controls for each sentinel SNP
   1. Download 8rsids.DEMO #8 SNPs associated with platelet trait variation in Gieger et al 2011
   2. Download the chromatin feature tracks and file 'DEMO.FeatureList'. You will need to specify correct directories in the DEMO.FeatureList where the chromatin features get downloaded 
   3. Download the template .conf file (you will need to change directories to make it work)
   4. Run GREGOR
      perl /your/path/to/GREGOR/script/GREGOR.pl --conf DEMO.conf 
      (this should take <10min to run)

   The output from this GREGOR file needs to be 'Mined' to create a table documenting overlaps between SNPs and chromatin     features
   The script GREGORmine.pl will output this table for you...

#Generate a table compiling SNP/feature overlaps
   1. Download GREGORmine.pl from Scripts
   2. From within the GREGOR-results output from before, run the following (should take less than 10min for DEMO)
perl /your/path/to/GREGORmine.pl /your/path/to/GREGOR-results.DEMO/ DEMO.Mined.txt wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak

   Output will be a table DEMO.Mined.txt
      This is a table with feature names and SNP characteristics as headers; MAF, Dist to Nearest Gene, and # LD SNPs, as well as 0/1 values describing overlaps comprise the rest of the table.

#Run a LASSO model
   1. Download GLMnet_Lasso_Rcode_creation.DEMO.pl from Scripts (this will help create the LASSO script to run)
   2. Create an R code to run GLMnet LASSO
      perl /your/path/go/GLMnet_Lasso_Rcode_creation.DEMO.pl 190601 1 7 '2,4:12'

   Usage: GLMnet_Lasso_Rcode_creation.DEMO.pl by  _CST_ created 160118. Please enter:
     perl GLMnet_Lasso_Rcode_creation.pl> <date> <run #> <#features> <exact feature column numbers (i.e., '2,4:18') this         does NOT need parentheses but does need to quotes to work>
     and output will be a file named Date.GLMnet_Lasso_Rcode_#col#.columns_run#run#.R
  
-- this will probably give you an 'empty model' because of a convergence issue, or perhaps a bad model (see DEMO.7.factors.run1.AUCplot.pdf). This is because there are not enough features in this small DEMO, but this at least gives you a feel for practical aspects of how we created our model -- 

