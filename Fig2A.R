#producing input files
#multiBigwigSummary BED-file --BED c_elegans.WS220.gene.annotations.sorted.classes.new_TSS_point_1kb.bed -b MACS141_e5_DPY27_N2_MxEmb_SE30_CJ39_CJ132_SE172_average.bw -o DPY27_N2_MxEmb_SE30_CJ39_CJ132_SE172_average_InputSubt_1KbWBTSS.npz --outRawCounts DPY27_N2_MxEmb_SE30_CJ39_CJ132_SE172_average_InputSubt_1KbWBTSS_raw_counts.txt

####packages######
library(ggplot2)
#library(ggrepel)

####functions#####
`%!in%` = Negate(`%in%`)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x)  | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

####run###########
#read in input files
N2WBTSS<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/DPY27_vs_GeneExpression_plots/DPY27_N2_MxEmb_SE30_CJ39_CJ132_SE172_average_InputSubt_1KbWBTSS_raw_counts.txt',header=T)
colnames(N2WBTSS)<-c("chr","start","end","DPY27embAvg")
rnaSeq<-read.delim(file='~/Desktop/Ercan_Lab/Annotations/ChromatinPaper/N2_MxEmb_RNAseq_max.txt',header=T)
rnaSeq<-rnaSeq[,-c(8:12)]
groTSS<-read.delim(file='~/Desktop/Ercan_Lab/Annotations/ChromatinPaper/Meyer_GroSeq_c_elegans.WS220.WS230.gene.annotation.TSS.all.sorted_points_1kb_chr.bed',header=F)
groTSS<-groTSS[,c(1,2,3,4,6)]
colnames(groTSS)<-c("chr","start","end",'gene','strand')
WBTSS<-read.delim(file='~/Desktop/Ercan_Lab/Annotations/ChromatinPaper/c_elegans.WS220.gene.annotations.sorted.classes.new_TSS_point.txt',header=F)
colnames(WBTSS)<-c("chr","start","end",'gene','strand')
ref<-read.table(file='~/Desktop/Ercan_Lab/Chromatin_Paper/Plot_annotations/c.elegans.WS220.gene.wbid.chr.start.stop.txt',header=T,stringsAsFactors=F)

#make sure all files have gene and WBID values so they can be matched against each other
groTSS$WBID<-ref$WBID[match(groTSS$gene,ref$gene)]
WBTSS$WBID<-ref$WBID[match(WBTSS$gene,ref$gene)]
N2WBTSS$gene<-WBTSS$gene
N2WBTSS$WBID<-WBTSS$WBID

#find which genes are expressed and which are not
#expressed: FPKM > 1 and in the GroSeq dataset
#not expressed: FPKM = 0 and not in the GroSeq dataset
#maternally loaded: FPKM > 1 and not in the GroSeq dataset
#4816 genes are excluded from this analysis because they have FPKM values between 0 and 1
rnaSilent<-rnaSeq[which(rnaSeq$mixed.embs.N2.mean.FPKM==0),]
rnaExpressed<-rnaSeq[which(rnaSeq$mixed.embs.N2.mean.FPKM > 1),]

expressedGenes<-WBTSS[which(WBTSS$WBID %in% groTSS$WBID & WBTSS$WBID %in% rnaExpressed$Gene_WB_ID),]
nrow(expressedGenes[which(expressedGenes$chr=='chrX'),])
expressedGenes$cat<-'\nExpressed \n(n=678)'
silentGenes<-WBTSS[which(WBTSS$WBID %!in% groTSS$WBID & WBTSS$WBID %in% rnaSilent$Gene_WB_ID),]
nrow(silentGenes[which(silentGenes$chr=='chrX'),])
silentGenes$cat<-'\nNot Expressed \n(n=367)'
maternalGenes<-WBTSS[which(WBTSS$WBID %!in% groTSS$WBID & WBTSS$WBID %in% rnaExpressed$Gene_WB_ID),]
nrow(maternalGenes[which(maternalGenes$chr=='chrX'),])
maternalGenes$cat<-'\nMaternally Loaded \n(n=1056)'

genes<-rbind(expressedGenes,silentGenes,maternalGenes)

#assign expression categories to the TSS in the N2 DPY27 binding input
N2WBTSS$cat<-as.factor(genes$cat[match(N2WBTSS$WBID,genes$WBID)])

#get the data for the plot: on the X chr and either expressed or not expressed
dataToPlot<-N2WBTSS[which(N2WBTSS$chr=='chrX' & !is.na(N2WBTSS$cat)),]

#transform dpy27 binding by log10 to make it easier to see all data on the plot
dataToPlot$DPY27embAvglog10<-log10(dataToPlot$DPY27embAvg)

#reorder the boxes so they are in the correct order in the plot
dataToPlot$cat <- factor(dataToPlot$cat, levels = c("\nExpressed \n(n=678)","\nNot Expressed \n(n=367)","\nMaternally Loaded \n(n=1056)"))

#get the TSS that are outliers
#most interesting outliers: not expressed but high DCC binding or low DCC binding but expressed
expToPlot<-dataToPlot[which(dataToPlot$cat=='\nExpressed \n(n=678)'),]
silToPlot<-dataToPlot[which(dataToPlot$cat=='\nNot Expressed \n(n=367)'),]
matToPlot<-dataToPlot[which(dataToPlot$cat=='\nMaternally Loaded \n(n=1056)'),]

expOutlier<-expToPlot[which(is_outlier(expToPlot$DPY27embAvg)==T),]
silOutlier<-silToPlot[which(is_outlier(silToPlot$DPY27embAvg)==T),]
matOutlier<-matToPlot[which(is_outlier(matToPlot$DPY27embAvg)==T),]

outliers<-rbind(expOutlier,silOutlier,matOutlier)
outliers$gene_name<-rnaSeq$Gene_Public_Name[match(outliers$WBID,rnaSeq$Gene_WB_ID)]



#plot the boxplot
bxplot<-ggplot(dataToPlot, aes(x=cat,y=DPY27embAvglog10,fill=cat)) + geom_boxplot(notch=TRUE) + theme_classic() + 
  scale_fill_manual(values=c("chartreuse4","dodgerblue3","red")) + 
  labs(fill="",y= "log10(DPY27 ChIP Enrichment)", title= "DPY27 Binding at X Gene TSS in WT MxEmb",x='ChrX Genes') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.title = element_text(size=10))+
  #geom_text(data = outliers, aes(label = gene_name),nudge_x = 0.05,hjust = 0)+
  stat_boxplot(geom ='errorbar')


print(bxplot)
ggsave("~/Desktop/GroTSS_and_RNAseq_Exp_vs_dpy27N2MxEmb_bxplot_fixedNnumbers.pdf", width = 8, height = 7)


outFile<-outliers
outFile$cat <- gsub('\n', '', outFile$cat)
write.table(outFile,file='~/Desktop/DPY27enrichOutliers_expByGroSeq_and_RNAseq.txt',sep='\t',quote=F,col.names=T,row.names=F)




##plot a histogram of the data
hist <- ggplot(dataToPlot, aes(DPY27embAvg,col=cat)) +
  geom_line(aes(y=..count..),stat="density") + theme_classic() + scale_color_manual(values=c("dodgerblue3","chartreuse4")) +
  labs(col="ChrX Genes",x= "DPY27 ChIP Enrichment", y="Number of Genes", title= "DPY27 Binding at X Gene TSS in WT MxEmb") +  
  theme(plot.title = element_text(size=10))
print(hist)

ggsave("~/Desktop/GroTSS_and_RNAseq_Exp_vs_dpy27N2MxEmb_histogram.pdf", width = 7, height = 6)


##
#try other gene annotations

geneType<-read.delim(file='~/Desktop/Rechtsteiner.gene.set.AWG.040210.mes-4-1_condensindataJan2013.txt',header=T)
geneType<-geneType[,c(1:4,13,17:25)]

ubiqExp<-geneType[which(geneType$ubiquitous.expressed==1),]
ubiqExp$cat<-'Ubiquitous Expressed \n(n=2580)'
silent<-geneType[which(geneType$silent==1),]
silent$cat<-'Silent \n(n=415)'
germlineEnrich<-geneType[which(geneType$germline.enriched==1),]
germlineEnrich$cat<-'Germline Enriched \n(n=2243)'
spermatogenesis<-geneType[which(geneType$spermatogenesis==1),]
spermatogenesis$cat<-'Spermatogenesis \n(n=861)'
germlineExp<-geneType[which(geneType$germline.expressed==1),]
germlineExp$cat<-'Germline Expressed \n(n=4693)'
germlineSpecific<-geneType[which(geneType$germline.specific==1),]
germlineSpecific$cat<-'Germline Specific \n(n=2243)'
somaSpecificAND<-geneType[which(geneType$soma.specific..AND.==1),]
somaSpecificAND$cat<-'Soma Specific AND \n(n=323)'
somaSpecificOR<-geneType[which(geneType$soma.specific..OR.==1),]
somaSpecificOR$cat<-'Soma Specific OR \n(n=1182)'
embExp<-geneType[which(geneType$embryo.expressed==1),]
embExp$cat<-'Embryo Expressed \n(n=797)'

genesToPlot<-rbind(ubiqExp,silent,germlineEnrich,spermatogenesis,germlineExp,germlineSpecific,somaSpecificAND,somaSpecificOR,embExp)

genesToPlot$DPY27embAvg<-N2WBTSS$DPY27embAvg[match(genesToPlot$Wormbase.ID,N2WBTSS$WBID)]
genesToPlot$DPY27embAvglog10<-log10(genesToPlot$DPY27embAvg)

#plot the boxplot
bxplot<-ggplot(genesToPlot, aes(x=cat,y=DPY27embAvglog10,fill=cat)) + stat_boxplot(geom ='errorbar', width = 0.4) + geom_boxplot(notch=F) + 
  theme_classic() + 
  #scale_fill_manual(values=c("chartreuse4","dodgerblue3","red")) + 
  labs(fill="",y= "log10(DPY27 ChIP Enrichment)", title= "DPY27 Binding at X Gene TSS in WT MxEmb",x='ChrX Genes') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.title = element_text(size=10)) #+
  #geom_text(data = outliers, aes(label = gene_name),nudge_x = 0.05,hjust = 0)


print(bxplot)
ggsave("~/Desktop/RechtsteinerGenes_vs_dpy27N2MxEmb_log10_bxplot.pdf", width = 8, height = 7)




