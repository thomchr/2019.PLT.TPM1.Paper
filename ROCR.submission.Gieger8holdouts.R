#### Do this stuff first ####
#ssh -X thomc and xbash
#Put this in myself into Unix command line and then the rest runs on its own

#bsub -o out%J -e error%J -m Interdictor R CMD BATCH ROCR.submission.R

setwd("/project/voight_ML/thomc/thomc_results/ROCR/")
LASSO <- read.table(file="Gieger8.LASSO.forROCR", header=T)
dim(LASSO)
str(LASSO)

GWAVA <- read.table(file="Gieger8.GWAVA.forROCR", header=T)
dim(GWAVA)
str(GWAVA)

CADD <- read.table(file="Gieger8.CADD.forROCR", header=T)
dim(CADD)
str(CADD)

DeepSea <- read.table(file="Gieger8.DeepSea.forROCR", header=T)
dim(DeepSea)
str(DeepSea)

###########
#Plot an ROC curve

library(ROCR)

#Lasso
pred <- prediction( LASSO$LASSO, LASSO$Index )
perf <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

#FATHMMKL
#pred <- prediction( FATHMMMKL$fathmmMKL_coding_score, FATHMMMKL$Index )
#perf2 <- performance( pred, "tpr", "fpr" )
#auc.perf = performance(pred, measure = "auc")
#auc.perf@y.values

#CADD
pred <- prediction( CADD$CADD_phred, CADD$Index )
perf3 <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

#GWAVA
pred <- prediction( GWAVA$GWAVA_Score, GWAVA$Index)
perf4 <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

#DeepSea
pred <- prediction( DeepSea$neglogFSS, DeepSea$Index)
perf5 <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

pdf("Gieger8holdouts.ROCR.LASSO.DeepSea.CADD.GWAVA.performance.pdf")
plot( perf, main="Gieger8 Holdouts Index vs 500 random snp147 ctrls", col='black') #LASSO
plot( perf3, add = TRUE, col='red') #CADD
plot( perf4, add = TRUE, col='purple') #GWAVA purple
plot( perf5, add = TRUE, col='green') #DeepSea green
#plot( perf6, add = TRUE, col='yellow') # eigen yellow
legend("bottomright", legend=c("LASSO - AUC 0.82", "CADD - AUC 0.68", "GWAVA - AUC 0.76", "DeepSea - AUC 0.63"), col=c("black", "red", "purple", "green"), lty=1:2, cex=0.8)
dev.off()
