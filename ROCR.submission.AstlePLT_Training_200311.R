#### Do this stuff first ####
#ssh -X thomc and xbash
#Put this in myself into Unix command line and then the rest runs on its own

#bsub -o out%J -e error%J -m Interdictor R CMD BATCH ROCR.submission.R

###############
# what is new in this version vs 191108 is the addition of more controls (1328000, same as in training), ~2000 per sample, to appeal to reviewer 4
##############

setwd("/project/voight_ML/thomc/thomc_results/ROCR/")
LASSO <- read.table(file="AstlePLT.LASSO.cat.Training.forROCR", header=T) 
#dim(LASSO)
#lapply(LASSO, class)
#str(LASSO)

GWAVA <- read.table(file="AstlePLT.GWAVA.cat.Training.forROCR", header=T)
#dim(GWAVA)
#str(GWAVA)

CADD <- read.table(file="AstlePLT.CADD.cat.Training.forROCR", header=T)
#dim(CADD)
#str(CADD)

DeepSea <- read.table(file="AstlePLT.Training.DeepSEA.forROCR", header=T)
#dim(DeepSea)
#str(DeepSea)

###########
#Plot an ROC curve

library(ROCR)

#Lasso
#data(LASSO)
pred1 <- prediction( LASSO$LASSO_Score, LASSO$Index )
perf1 <- performance( pred1, "tpr", "fpr" )
auc.perf1 = performance(pred1, measure = "auc")
auc.perf1@y.values
#LASSOAUC <- bquote(.(auc.perf1@y.values)) 

#Optimal cut point (https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/)
opt.cut = function(perf1, pred1){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf1@x.values, perf1@y.values, pred1@cutoffs)
}
print(opt.cut(perf1, pred1))

#partial AUC
pauc.perf1 = performance(pred1, measure = "auc", fpr.stop=0.2)
pauc.perf1@y.values

#proc.perf1 = pROC(pred, fpr.stop=0.1) #can't find this program pROC, don't think I need it anyway

#CADD
#data(CADD)
pred <- prediction( CADD$PHRED, CADD$Index )
perf3 <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

#GWAVA
pred4 <- prediction( GWAVA$GWAVA_Score, GWAVA$Index)
perf4 <- performance( pred4, "tpr", "fpr" )
auc.perf4 = performance(pred4, measure = "auc")
auc.perf4@y.values
#GWAVAAUC <- bquote(auc.perf4@y.values)

#partial AUC
pauc.perf4 = performance(pred4, measure = "auc", fpr.stop=0.2)
pauc.perf4@y.values

#proc.perf4 = pROC(pred4, fpr.stop=0.1) #can't find this program	pROC, don't think I need it anyway


#DeepSea
pred <- prediction( DeepSea$neglogFSS, DeepSea$Index)
perf5 <- performance( pred, "tpr", "fpr" )
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

pdf("200311.AstlePLT.Training.ROCR.LASSO.GWAVA.performance.pdf")
plot( perf1, main="Astle PLT Training vs 100 ctrls each", col='black') #LASSO
#plot( perf2, add = TRUE, col='blue') #FATHMMKL
#plot( perfa, add = TRUE, col='orange') #LASSOvs1328000 ctrls
plot( perf3, add = TRUE, col='purple') #CADD purple
plot( perf4, add = TRUE, col='red') #GWAVA red
#plot( perfb, add = TRUE, col='black') #GWAVA1328000ctrls
plot( perf5, add = TRUE, col='green') #DeepSea green
#plot( perf6, add = TRUE, col='yellow') # eigen yellow
#plot( perf7, add = TRUE, col='orange') # eigen_pc oranage

legend("bottomright", legend=c("LASSO - AUC 0.799", "GWAVA - AUC 0.745", "CADD - AUC 0.581", "DeepSEA - AUC 0.629"), col=c("black", "red", "purple", "green"), lty=1:2, cex=0.8)
#legend("bottomright", legend=c("LASSO", "CADD", "GWAVA", "DeepSea", "FATHMMKL", "EIGEN", "EIGEN_PC"), col=c("black", "red", "purple", "green", "blue", "yellow", "orange"), lty=1:2, cex=0.8)
dev.off()

pdf("200311.AstlePLT.Training.ROCR.LASSO.GWAVA_fprstop0.2_performance.pdf")
plot(perf1, main="Astle PLT Training vs 54,000 random snp147 ctrls", xlim=c(0,0.2), col='black') #LASSO
#plot( perfa, add = TRUE, col='orange') #LASSOvs1328000 ctrls
plot(perf4, add = TRUE, col='red') #GWAVA red
#plot( perfb, add = TRUE, col='black') #GWAVA1328000ctrls
legend("bottomright", legend=c("LASSO - pAUC 0.0831", "GWAVA - pAUC 0.0723"), col=c("black", "red"), lty=1:2, cex=0.8)
dev.off()

#plot(perf, main="ROCR fingerpainting toolkit", colorize=TRUE, xlab="Mary's axis", ylab="", box.lty=7, box.lwd=5, box.col="gold", lwd=17, colorkey.relwidth=0.5, xaxis.cex.axis=2, xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=2, yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1, xaxis.lwd=2, yaxis.lwd=3, yaxis.col.axis="orange", cex.lab=2, cex.main=2)

