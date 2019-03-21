#load package including smaller subitems, as specified in Gviz.pdf
library(Gviz) 

chr <- 'chr15' #this is how everything downstream can be generic. Can just change the 'chr#' for each
gen <- 'hg19'

library(GenomicRanges)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(bands=read.delim("ideogramtrack", col.names=c("chrom", "chromStart", "chromEnd", "name", "gieStain")), genome="hg19")

#The Gviz.pdf way to plot GeneRegion:
grtrack <- GeneRegionTrack(range="Homo_sapiens.GRCh37.55.protein_coding.No-AC.gtf", genome = "hg19", name = "Genes")

class(grtrack)
nrow(grtrack)
head(grtrack)

#load the relevant tracks into Gviz
#To get bedGraph files ready:
#note that I had to sed 1d this file to get rid of header and make .bedGraph

setwd("/project/voight_ML/thomc/thomc_results/ScorePlottingR/Gviz")

#GeneRegion Track - had to download my own because auto-lookup to UCSC kept timing out (firewall issue)

GATA1 <- DataTrack(range="GSM607949_MK_GATA1_hg19.cp.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GATA1", fill.histogram = "purple", col.histogram = "purple") 
class(GATA1)

FLI1 <- DataTrack(range = "GSM607952_MK_FLI1_hg19.cp.bedGraph", genome = gen, type = "l", chromosome = chr, name = "FLI1", fill.histogram = "purple", col.histogram = "purple")
class(FLI1)

###now go to Gviz directory for score tracks

WGscores <- DataTrack(range = "160403.GenomeMPV.Scored.160222LassoCoefs.txt.Padded.sorted.cp.bedGraph", genome = gen, type = "l", chromosome = chr, name = "Whole Genome Scores", fill.histogram = "gray", col.histogram = "gray", ylim=c(0,2.5))
class(WGscores)

GWASproxyScores <- DataTrack(range = "190317.LDLink.rs3809566.EURr2_0.7_assocSNPs.bedGraph", genome = gen, type = "l", chromosome = chr, name = "LDlink EURr2=0.7 Scores", fill.histogram = "cyan", col.histogram = "cyan", ylim=c(0,2.5))
class(GWASproxyScores)

GWASscores <- DataTrack(range = "Gieger73Index.Location.Score.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GWAS SNP Scores", fill.histogram = "black", col.histogram = "black", ylim=c(0,2.5)) #make sure to specify this ylim to keep consistent so that when overlayed the scores match on y axis

ylims <- extendrange(range(c(values(GWASproxyScores), values(GWASscores), values (WGscores)))) #to get them both on y axis
ot <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores, WGscores), ylim = c(0,2.5))
ot2 <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores), ylim = c(0,2.5))

GmH3K4m1 <- DataTrack(range = "wgEncodeBroadHistoneGm12878H3k04me1StdSigV2.bigWig", genome = gen, type = "l", chromosome = chr, name = "Gm12878 H3K4me1", fill.histogram = "red", col.histogram = "red")

GmH3K4m2 <- DataTrack(range = "wgEncodeBroadHistoneGm12878H3k4me2StdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "Gm12878 H3K4me2", fill.histogram = "red", col.histogram = "red")

GmH3K27ac <- DataTrack(range = "wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "Gm12878 H3K27ac", fill.histogram = "red", col.histogram = "red")

HuvH3K4m2 <- DataTrack(range = "wgEncodeBroadHistoneHuvecH3k4me2StdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "HUVEC H3K4me2", fill.histogram = "orange", col.histogram = "orange")

K562H3K36m3 <- DataTrack(range = "wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "K562 H3K36me3", fill.histogram = "darkblue", col.histogram = "darkblue")

K562H3K79m2 <- DataTrack(range = "wgEncodeBroadHistoneK562H3k79me2StdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "K562 H3K79me2", fill.histogram = "darkblue", col.histogram = "darkblue")

K562RbBP5 <- DataTrack(range = "wgEncodeBroadHistoneK562Rbbp5a300109aStdSig.bigWig", genome = gen, type = "l", chromosome = chr, name = "K562 RbBP5", fill.histogram = "darkblue", col.histogram = "darkblue")

pdf("TPM1.100kb.wWGS.locus.pdf")
plotTracks(list(itrack, ot2, ot, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.ylim2.5.locus.pdf")
plotTracks(list(ot, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off()


#####################
#to get ylim max 2 need to recharacterize each of the Score tracks
######################

WGscores <- DataTrack(range = "160403.GenomeMPV.Scored.160222LassoCoefs.txt.Padded.sorted.cp.bedGraph", genome = gen, type = "l", chromosome = chr, name = "Whole Genome Scores", fill.histogram = "gray", col.histogram = "gray", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial')
class(WGscores)

GWASproxyScores <- DataTrack(range = "190317.LDLink.rs3809566.EURr2_0.7_assocSNPs.bedGraph", genome = gen, type = "l", chromosome = chr, name = "EUR r2=0.7 Scores", fill.histogram = "cyan", col.histogram = "cyan", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial')
class(GWASproxyScores)

GWASscores <- DataTrack(range = "Gieger73Index.Location.Score.bedGraph", genome = gen, type = "l", chromosome = chr, name = "GWAS SNP Scores", fill.histogram = "black", col.histogram = "black", ylim=c(0,2), col.sampleNames="black",background.title="white",col.axis="black",col.title="black", fontfamily='Arial') #make sure to specify this ylim to keep consistent so that when overlayed the scores match on y axis

ot <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores, WGscores), ylim = c(0,2))
ot2 <- OverlayTrack(trackList = list(GWASproxyScores, GWASscores), ylim = c(0,2))

pdf("TPM1.70kb.SNPsGRMarks.locus.pdf")
plotTracks(list(ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63320000, to = 63390000, type = "histogram", scale = 0.25, labelPos = 'below') 
dev.off()

pdf("TPM1.100kb.wWGS.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.100kb.locus.pdf")
plotTracks(list(itrack, ot2, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.100kb.SNPs.GeneRegiononly.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, gtrack), chromosome = chr, from = 63291996, to = 63391996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.10kb.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63336996, to = 63346996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.1kb.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, GATA1, FLI1, K562RbBP5, K562H3K36m3, K562H3K79m2, GmH3K4m1, GmH3K4m2, GmH3K27ac, HuvH3K4m2, gtrack), chromosome = chr, from = 63341496, to = 63342495, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = 'meta', transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.1Mb.GRonly.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, gtrack), chromosome = chr, from = 62841996, to = 63841996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = TRUE, shape = "arrow", transcriptAnnotation = 'gene_name')
dev.off() 

pdf("TPM1.1Mb.GRonly.noannotation.locus.pdf")
plotTracks(list(itrack, ot2, grtrack, gtrack), chromosome = chr, from = 62841996, to = 63841996, type = "histogram", scale = 0.25, labelPos = 'below', collapseTranscripts = TRUE, shape = "arrow")
dev.off() 

quit()