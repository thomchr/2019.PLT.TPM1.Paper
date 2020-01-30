#load package including smaller subitems, as specified in Gviz.pdf
library(Gviz) 

chr <- 'chr13' #this is how everything downstream can be generic. Can just change the 'chr#' for each
gen <- 'hg19'
loc <- 33164418

library(GenomicRanges)
gtrack <- GenomeAxisTrack()
#itrack <- IdeogramTrack(bands=read.delim("ideogramtrack", col.names=c("chrom", "chromStart", "chromEnd", "name", "gieStain")), genome="hg19")

#The Gviz.pdf way to plot GeneRegion:
grtrack <- GeneRegionTrack(range="Homo_sapiens.GRCh37.55.protein_coding.No-AC.gtf", genome = "hg19", name = "Genes")

#class(grtrack)
#nrow(grtrack)
#head(grtrack)

#load the relevant tracks into Gviz
#To get bedGraph files ready:
#note that I had to sed 1d this file to get rid of header and make .bedGraph

setwd("/project/voight_ML/thomc/thomc_results/ScorePlottingR/Gviz/Gviz.AstleModels")

#GeneRegion Track - had to download my own because auto-lookup to UCSC kept timing out (firewall issue)

###now go to Gviz directory for score tracks

WGplt <- DataTrack(range = "AstlePLT.WG.bedGraph", genome = gen, type = "l", chromosome = chr, name = "Whole Genome Scores", fill.histogram = "gray", col.histogram = "gray", ylim=c(-1,5))

WGred <- DataTrack(range = "AstleRedCell.WG.bedGraph", genome = gen, type = "l", chromosome = chr, name = "Whole Genome Scores", fill.histogram = "gray", col.histogram = "gray", ylim=c(-1,5))

ProxyPlt <- DataTrack(range = "AstlePLT.r7.bedGraph", genome = gen, type = "l", chromosome = chr, name = "LDlink EURr2=0.7 Scores", fill.histogram = "cyan", col.histogram = "cyan", ylim=c(-1,5))

ProxyRed <- DataTrack(range = "AstleRedCell.r7.bedGraph", genome = gen, type = "l", chromosome = chr, name = "LDlink EURr2=0.7 Scores", fill.histogram = "cyan", col.histogram = "cyan", ylim=c(-1,5))

IndexPlt <- DataTrack(range = "AstlePLT.Index.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GWAS SNP Scores", fill.histogram = "black", col.histogram = "black", ylim=c(-1,5)) #make sure to specify this ylim to keep consistent so that when overlayed the scores match on y axis

IndexRed <- DataTrack(range = "AstleRedCell.Index.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GWAS SNP Scores", fill.histogram = "black", col.histogram = "black", ylim=c(-1,5)) #make sure to specify this ylim to keep consistent so that when overlayed the scores match on y axis

HSC <- DataTrack(range = "HSC.bedGraph", genome = gen, type = "l", chromosome = chr, name = "HSC", fill.histogram = "blue", col.histogram = "blue") #, ylim=c(-1,5))
CD34 <- DataTrack(range = "CD34.bedGraph", genome = gen, type = "l", chromosome = chr, name = "HSC", fill.histogram = "orange", col.histogram = "orange") #, ylim=c(-1,5))
MEP <- DataTrack(range = "MEP.bedGraph", genome = gen, type = "l", chromosome = chr, name = "HSC", fill.histogram = "red", col.histogram = "red") #, ylim=c(-1,5))
Ery <- DataTrack(range = "Erythro.bedGraph", genome = gen, type = "l", chromosome = chr, name = "HSC", fill.histogram = "dark red", col.histogram = "dark red") #, ylim=c(-1,5))

ylims <- extendrange(range(c(values(WGplt), values(WGred), values(ProxyPlt), values(ProxyRed), values(IndexPlt), values(IndexRed)))) #to get them both on y axis
ot <- OverlayTrack(trackList = list(WGplt, ProxyPlt, IndexPlt), ylim = c(-1,5))
ot2 <- OverlayTrack(trackList = list(WGred, ProxyRed, IndexRed), ylim = c(-1,5))


pdf("PDS5B.100kb.wWGS.locus.pdf")
ub <- loc + 100000
lb <- loc - 100000
plotTracks(list(ot2, ot, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.10kb.wWGS.locus.pdf")
ub <- loc + 10000
lb <- loc - 10000
plotTracks(list(ot2, ot, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 



#############
# no whole genome scores, just proxies
ot <- OverlayTrack(trackList = list(ProxyPlt, IndexPlt), ylim = c(-1,5))
ot2 <- OverlayTrack(trackList = list(ProxyRed, IndexRed), ylim = c(-1,5))


pdf("PDS5B.1Mb.locus.pdf")
ub <- loc + 1000000
lb <- loc - 1000000
plotTracks(list(WGplt, WGred, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.100kb.WGS.locus.pdf")
ub <- loc + 100000
lb <- loc - 100000
plotTracks(list(WGplt, WGred, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = TRUE, shape = "arrow", transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.100kb.locus.pdf")
ub <- loc + 100000
lb <- loc - 100000
plotTracks(list(ot, grtrack, HSC, CD34, MEP, Ery, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.10kb.locus.pdf")
ub <- loc + 10000
lb <- loc - 10000
plotTracks(list(ot, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.1kb.locus.pdf")
ub <- loc + 1000
lb <- loc - 1000
plotTracks(list(ot2, ot, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 



quit()




#####################
#to get ylim max 2 need to recharacterize each of the Score tracks
######################

#WGscores <- DataTrack(range = "160403.GenomeMPV.Scored.160222LassoCoefs.txt.Padded.sorted.cp.bedGraph", genome = gen, type = "l", chromosome = chr, name = "Whole Genome Scores", fill.histogram = "gray", col.histogram = "gray", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial')
#class(WGscores)

#GWASproxyScores <- DataTrack(range = "plinkproxies.EURr2_0.7_assocSNPs.bedGraph", genome = gen, type = "l", chromosome = chr, name = "EUR r2=0.7 Scores", fill.histogram = "cyan", col.histogram = "cyan", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial')
#class(GWASproxyScores)
#used to be 190317.LDLink.rs3809566.EURr2_0.7_assocSNPs.bedGraph

#GWASscores <- DataTrack(range = "Gieger73Index.Location.Score.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GWAS SNP Scores", fill.histogram = "black", col.histogram = "black", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial') #make sure to specify this ylim to keep consistent so that when overlayed the scores match on y axis

#ot <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores, WGscores), ylim = c(0,2))
#ot2 <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores), ylim = c(0,2))

#pdf("PDS5B.70kb.SNPsGRMarks.locus.pdf")
#plotTracks(list(ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63320000, to = 63390000, type = "histogram", scale = 0.25, labelPos = 'below') 
#dev.off()

#pdf("PDS5B.100kb.wWGS.locus.pdf")
#plotTracks(list(ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
#dev.off() 

pdf("PDS5B.100kb.locus.pdf")
ub <- loc + 100000
lb <- loc - 100000
plotTracks(list(ot2, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.100kb.SNPs.GeneRegiononly.locus.pdf")
plotTracks(list(ot2, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.10kb.locus.pdf")
lb <- loc - 10000
ub <- loc + 10000
plotTracks(list(ot2, grtrack, gtrack), chromosome = chr, from = lb, to = ub, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

#pdf("PDS5B.1kb.locus.pdf")
#plotTracks(list(ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63341496, to = 63342495, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
#dev.off() 

pdf("PDS5B.1Mb.GRonly.locus.pdf")
plotTracks(list(ot2, grtrack, gtrack), chromosome = chr, from = 62841996, to = 63841996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = TRUE, shape = "arrow", transcriptAnnotation = 'gene_name')
dev.off() 

pdf("PDS5B.1Mb.GRonly.noannotation.locus.pdf")
plotTracks(list(ot2, grtrack, gtrack), chromosome = chr, from = 62841996, to = 63841996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = TRUE, shape = "arrow")
dev.off() 

quit()
