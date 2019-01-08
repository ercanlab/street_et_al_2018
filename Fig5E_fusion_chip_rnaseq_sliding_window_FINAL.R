##########################################
###plot DPY27 binding, RNA-seq, and H3K27ac in sliding windows

########input##############################################################################
#module load deeptools/3.0.2
#multiBigwigSummary BED-file --BED /scratch/las821/Meyer_GroSeq_c_elegans.WS220.WS230.gene.annotation.TSS.all.sorted_points_1kb_chr.bed -b /scratch/las821/ChromatinPaperAnalysis/Fig6B/15eh1_L3_H3K27ac_SE351_SE352_ext140_ext141_average.bw -o test.npz --outRawCounts 15eh1_L3_H3K27ac_SE351_SE352_ext140_ext141_average_1KbGroTSS_raw_counts.txt
#multiBigwigSummary BED-file --BED /scratch/las821/Meyer_GroSeq_c_elegans.WS220.WS230.gene.annotation.TSS.all.sorted_points_1kb_chr.bed -b /scratch/las821/ChromatinPaperAnalysis/Fig6B/N2_L3_H3K27ac_FE1_CTTGTA_L001_ME5_CAGATC_L005_5054_average.bw -o test.npz --outRawCounts N2_L3_H3K27ac_FE1_CTTGTA_L001_ME5_CAGATC_L005_5054_average_1KbGroTSS_raw_counts.txt

########functions##############################################################################

get_relative_diff<- function(df,mut,wt) {
  rel_diff<-c()
  for (i in 1:nrow(df)){
    rel_diff[i]<-((df[i,mut]-df[i,wt])/min(abs(df[i,wt]),abs(df[i,mut])))
  }
  return(rel_diff)
}

getSlidingWindow<-function(window,step){
  chrLengths<-data.frame(chr=c("chrI","chrII","chrIII","chrIV","chrV","chrX"),length=c(15072423,15279345,13783700,17493793,20924149,17718866),stringsAsFactors=TRUE)
  output<-data.frame(chr=character(),start=integer(),end=integer(),stringsAsFactors=TRUE)
  for (i in 1:nrow(chrLengths)){
    len<-chrLengths[i,2]
    Chr<-chrLengths[i,1]
    start<-as.integer(seq(0,(len-(window-(step+1))),step))
    chr<-rep(Chr,length(start))
    end<-as.integer((start)+window)
    temp<-data.frame(chr,start,end)
    temp[length(start),3]=len
    output<-rbind(output,temp)
  }
  return(output)
}

getAverages<-function(data,dataCol,slidingWindows){
  output<-slidingWindows
  output$avg<-c()
  
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  
  for (i in 1:nrow(slidingWindows)){
    chr<-slidingWindows[i,1]
    start<-slidingWindows[i,2]
    end<-slidingWindows[i,3]
    avg<-mean(data[which(data[,1]==chr & data[,2]>start & data[,3]<end),dataCol],na.rm=TRUE)
    output[i,"avg"]<-avg
  }
  out<-output[which(!is.nan(output$avg)),]
  return(out)
}

tTest<-function(data,windows){
  a<-windows
  b<-data
  a$ttest<-c()
  b$group<-c(rep("other",nrow(b)))
  for (i in 1:nrow(a)){
    b[which(b$chr== a[i,'chr'] & b$start > a[i,"start"] & b$stop < a[i,"end"]),"group"]<- "test"
    z<-b[which(b$group=='test'),]
    if (nrow(z) > 2){
      a[i,"ttest"]<-t.test(log2ratioZscore~group,data=b,alternative = "two.sided")$p.value
      b[which(b$start > a[i,"start"] & b$stop < a[i,"end"]),"group"]<- "other"
    }
    else{
      a[i,"ttest"]<-'NA'
    }
  }
  return(a)
}

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

#######run################################################################################################
fusH3K27ac<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/fusion_plots/ChIPseq_plot_input_data/15eh1_L3_H3K27ac_SE351_SE352_ext140_ext141_InputSubt_1KbGroTSS_raw_counts.txt',header=T)
colnames(fusH3K27ac)<-c("chr","start","stop","fusAvg")
wtH3K27ac<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/fusion_plots/ChIPseq_plot_input_data/N2_L3_H3K27ac_fe1_me5_5054_InputSubt_1KbGroTSS_raw_counts.txt',header=T)
colnames(wtH3K27ac)<-c("chr","start","stop","wtAvg")
fusH3K4me3<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/fusion_plots/ChIPseq_plot_input_data/H3K4me3_N2_15eh1_L3_InputSubt_1kbGroTSS_raw_counts.txt',header=T)
colnames(fusH3K4me3)<-c("chr","start","stop",'wtAvg',"fusAvg")

H3K4me3Data<-fusH3K4me3[which(fusH3K4me3$chr=="chrV" & fusH3K4me3$start > 5000000),]
H3K4me3Data$log2ratio<-log2(H3K4me3Data$fusAvg/H3K4me3Data$wtAvg)
H3K4me3DataToPlot<-H3K4me3Data[which(H3K4me3Data$log2ratio!="NaN"),] #remove TSS with a NaN log2ratio
H3K4me3DataToPlot$log2ratioZscore<-(H3K4me3DataToPlot$log2ratio-mean(H3K4me3DataToPlot$log2ratio))/sd(H3K4me3DataToPlot$log2ratio) #get the Zscore of the log2ratio

H3K27acData<-fusH3K27ac[which(fusH3K27ac$chr=="chrV" & fusH3K27ac$start > 5000000),]
wtH3K27ac<-wtH3K27ac[which(wtH3K27ac$chr=="chrV" & wtH3K27ac$start > 5000000),]
H3K27acData$wtAvg<-wtH3K27ac$wtAvg
H3K27acData$log2ratio<-log2(H3K27acData$fusAvg/H3K27acData$wtAvg)
H3K27acDataToPlot<-H3K27acData[which(H3K27acData$log2ratio!="NaN"),] #remove TSS with a NaN log2ratio
H3K27acDataToPlot$log2ratioZscore<-(H3K27acDataToPlot$log2ratio-mean(H3K27acDataToPlot$log2ratio))/sd(H3K27acDataToPlot$log2ratio) #get the Zscore of the log2ratio

windows<-getSlidingWindow(100000,10000) #get sliding windows with 100000bp window and 10000bp step
spreadPlotWind<-windows[which(windows$chr=="chrV"),]
spreadPlotWind<-windows[which(windows$chr=="chrV" & windows$start > 5000000),]

avgSlidWindH3K27ac<-getAverages(H3K27acDataToPlot,"log2ratioZscore",spreadPlotWind) #get the averages of the log2ratioZscores within each window
avgSlidWindH3K27ac$mid<-avgSlidWindH3K27ac$start+((avgSlidWindH3K27ac$end-avgSlidWindH3K27ac$start)/2)

avgSlidWindH3K4me3<-getAverages(H3K4me3DataToPlot,"log2ratioZscore",spreadPlotWind) #get the averages of the log2ratioZscores within each window
avgSlidWindH3K4me3$mid<-avgSlidWindH3K4me3$start+((avgSlidWindH3K4me3$end-avgSlidWindH3K4me3$start)/2)

H3K27acdataTtest<-tTest(H3K27acDataToPlot,spreadPlotWind)
H3K27acdataTtest$avgZscore<-getAverages(H3K27acDataToPlot,"log2ratioZscore",H3K27acdataTtest)$avg
H3K27acdataTtest<-H3K27acdataTtest[which(!is.nan(H3K27acdataTtest$avgZscore)),]
H3K27acsigUpWindows<-H3K27acdataTtest[which(H3K27acdataTtest$avgZscore > mean(H3K27acdataTtest$avgZscore) & H3K27acdataTtest$ttest < 0.05),]
H3K27acsigUpWindows$mid<-H3K27acsigUpWindows$start+((H3K27acsigUpWindows$end-H3K27acsigUpWindows$start)/2)
H3K27acsigDownWindows<-H3K27acdataTtest[which(H3K27acdataTtest$avgZscore < mean(H3K27acdataTtest$avgZscore) & H3K27acdataTtest$ttest < 0.05),]
H3K27acsigDownWindows$mid<-H3K27acsigDownWindows$start+((H3K27acsigDownWindows$end-H3K27acsigDownWindows$start)/2)

H3K4me3dataTtest<-tTest(H3K4me3DataToPlot,spreadPlotWind)
H3K4me3dataTtest$avgZscore<-getAverages(H3K4me3DataToPlot,"log2ratioZscore",H3K4me3dataTtest)$avg
H3K4me3dataTtest<-H3K4me3dataTtest[which(!is.nan(H3K4me3dataTtest$avgZscore)),]
H3K4me3sigUpWindows<-H3K4me3dataTtest[which(H3K4me3dataTtest$avgZscore > mean(H3K4me3dataTtest$avgZscore) & H3K4me3dataTtest$ttest < 0.05),]
H3K4me3sigUpWindows$mid<-H3K4me3sigUpWindows$start+((H3K4me3sigUpWindows$end-H3K4me3sigUpWindows$start)/2)
H3K4me3sigDownWindows<-H3K4me3dataTtest[which(H3K4me3dataTtest$avgZscore < mean(H3K4me3dataTtest$avgZscore) & H3K4me3dataTtest$ttest < 0.05),]
H3K4me3sigDownWindows$mid<-H3K4me3sigDownWindows$start+((H3K4me3sigDownWindows$end-H3K4me3sigDownWindows$start)/2)

inputFusion<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/fusion_plots/ChIPseq_plot_input_data/15eh1_L3_DPY27_SEA149_SEA151_ext140_ext141_average_rawCounts_1kbGroTSS.txt',header=T)
colnames(inputFusion)<-c("chr","start","stop","avg")
inputFusion<-inputFusion[which(inputFusion$chr=="chrV" & inputFusion$start > 5000000),]
inputWT<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/fusion_plots/ChIPseq_plot_input_data/DPY27_JL_N2_L3_SE262_SEA153_average_rawCounts_1kbGroTSS.txt',header=T)
colnames(inputWT)<-c("chr","start","stop","avg")
inputWT<-inputWT[which(inputWT$chr=="chrV" & inputWT$start > 5000000),]

avgSlidWindFus<-getAverages(inputFusion,"avg",spreadPlotWind) #get the averages of the fusion chip enrichment within each window, remove windows with NaN avg
avgSlidWindWT<-getAverages(inputWT,"avg",spreadPlotWind) #get the averages of the wt chip enrichment within each window, remove windows with NaN avg
avgSlidWindWT$mid<-avgSlidWindWT$start+((avgSlidWindWT$end-avgSlidWindWT$start)/2)
avgSlidWindFus$mid<-avgSlidWindFus$start+((avgSlidWindFus$end-avgSlidWindFus$start)/2)
avgSlidWindFus<-avgSlidWindFus[which(!is.nan(avgSlidWindFus$avg)),]

inDat<-read.delim(file='~/Desktop/Ercan_Lab/DESeq/15eh1_vs_N2/results/deseq.summaryoverview.txt',header=T,stringsAsFactors = F)
a<-rownames(inDat)
inDat$gene<-a
rnaSeqData<-inDat[,c(5,7)]
colnames(rnaSeqData)<-c("log2","gene")

ref<-read.table(file='~/Desktop/Ercan_Lab/Chromatin_Paper/Plot_annotations/c.elegans.WS220.gene.wbid.chr.start.stop.txt',header=T,stringsAsFactors=F)

rnaSeqData$chr<-ref$chr[match(rownames(rnaSeqData),ref$gene)]
rnaSeqData$start<-ref$start[match(rownames(rnaSeqData),ref$gene)]
rnaSeqData$end<-ref$stop[match(rownames(rnaSeqData),ref$gene)]
rnaSeqData<-rnaSeqData[c(3,4,5,1,2)]
rnaSeqData<-rnaSeqData[which(rnaSeqData$chr=="chrV" & rnaSeqData$start > 5000000),]
avgSlidWindRNA<-getAverages(rnaSeqData,"log2",spreadPlotWind)
avgSlidWindRNA$mid<-avgSlidWindRNA$start+((avgSlidWindRNA$end-avgSlidWindRNA$start)/2)

colnames(rnaSeqData)<-c('chr','start','stop','log2ratioZscore','gene')
RNAdataTtest<-tTest(rnaSeqData,spreadPlotWind)
RNAdataTtest$avgZscore<-getAverages(rnaSeqData,"log2ratioZscore",RNAdataTtest)$avg
RNAdataTtest<-RNAdataTtest[which(!is.nan(RNAdataTtest$avgZscore)),]
RNAsigUpWindows<-RNAdataTtest[which(RNAdataTtest$avgZscore > mean(RNAdataTtest$avgZscore) & RNAdataTtest$ttest < 0.05),]
RNAsigUpWindows$mid<-RNAsigUpWindows$start+((RNAsigUpWindows$end-RNAsigUpWindows$start)/2)
RNAsigDownWindows<-RNAdataTtest[which(RNAdataTtest$avgZscore < mean(RNAdataTtest$avgZscore) & RNAdataTtest$ttest < 0.05),]
RNAsigDownWindows$mid<-RNAsigDownWindows$start+((RNAsigDownWindows$end-RNAsigDownWindows$start)/2)

c<-H3K27acData[which(H3K27acData$log2ratio!="NaN" & H3K27acData$chr=='chrV' & H3K27acData$start>5000000),]

pdf("~/Desktop/DPY27_binding_and_Zscorelog2_INsubt_H3K4me3_15eh1_over_N2_1KbGroTSS_and_RNAseq_log2_15eh1_overN2_100kb10kbSlidWind_log2ratioZscore_pval0_05.pdf", width=13, height=5) #plot the fusion and wt averages across the spreading region

par(mar=c(5,9,5,16))
options(scipen=999)
#plot DPY27 binding
plot(avgSlidWindFus$mid,avgSlidWindFus$avg,type='n',xlab='ChrV coordinate',ylab='DPY27 ChIP enrichment',main='DPY27 and H3K4me3 ChIP enrichment and Expression levels in X;V and WT',ylim = c(-3,6))
#legend("topright",inset=c(-0.375,-.25),c("X;V DPY27 ChIP \nenrichment","Zscore log2(X;V/WT) \nH3K27ac ChIP enrichment","log2(X;V/WT) \nExpression levels","H3K27ac Sig. Up \n(pval < 0.05)","H3K27acSig. Down \n(pval < 0.05)"),lty=c(1,1,1,1,1),lwd=c(1,1,1,4,4),pch=c(NA,NA,NA,NA,NA),col=c('#3333CC',"chartreuse4","red",'lightblue','deeppink'),xpd=TRUE,bty="n",pt.cex = 1,cex=0.8,y.intersp=2) #make a legend
#legend("topright",inset=c(-0.375,-.25),c("X;V DPY27 ChIP \nenrichment","Zscore log2(X;V/WT) \nH3K27ac ChIP enrichment","log2(X;V/WT) \nExpression levels","H3K27ac Sig. Up \n(pval < 0.05)","H3K27acSig. Down \n(pval < 0.05)","NaN log2 TSS","number log2 TSS"),lty=c(1,1,1,1,1,0,0),lwd=c(1,1,1,4,4,0,0),pch=c(NA,NA,NA,NA,NA,19,19),col=c('#3333CC',"chartreuse4","red",'lightblue','deeppink',"purple","orange"),xpd=TRUE,bty="n",pt.cex = 1,cex=0.8,y.intersp=2) #make a legend
legend("topright",inset=c(-0.375,-.25),c("X;V DPY27 ChIP \nenrichment","Zscore log2(X;V/WT) \nH3K4me3 ChIP enrichment","log2(X;V/WT) \nExpression levels","H3K4me3 Sig. Up \n(pval < 0.01)","H3K4me3 Sig. Down \n(pval < 0.01)","TSS"),lty=c(1,1,1,1,1,0),lwd=c(1,1,1,4,4,0),pch=c(NA,NA,NA,NA,NA,19),col=c('#3333CC',"chartreuse4","red",'lightblue','deeppink',"orange"),xpd=TRUE,bty="n",pt.cex = 1,cex=0.8,y.intersp=2) #make a legend
lines(smooth.spline(avgSlidWindFus$mid,avgSlidWindFus$avg,nknots=150),col="#3333CC",lwd=2,type='l') 
abline(h=0,lty='dashed',col='black') #abline at 0
#overlay plot of H3K27ac enrichment
par(new=T)
plot(avgSlidWind$mid,avgSlidWind$avg,type='n',xlab='',ylab='',main='',ylim = c(-3,6),axes = F)
lines(smooth.spline(avgSlidWind$mid,avgSlidWind$avg,nknots=150),col="chartreuse4",lwd=2,type='l')
#axis(side = 2)
mtext("Zscore of log2(X;V/WT) \nH3K4me3 ChIP enrichment \nor\n", side = 2, line = 3)
segments(x0=(sigUpWindows$start+200000),x1=(sigUpWindows$end-200000),y0=-2,y1=-2,col='lightblue',lwd=4) #significant up windows
segments(x0=(sigDownWindows$start+200000),x1=(sigDownWindows$end-200000),y0=-2,y1=-2,col='deeppink',lwd=4) #significant down windows
#overlay plot of RNAseq ratio
par(new=T)
plot(avgSlidWindRNA$mid,avgSlidWindRNA$avg,type='n',xlab='',ylab='',main='',ylim = c(-2,4),axes = F)
lines(smooth.spline(avgSlidWindRNA$mid,avgSlidWindRNA$avg,nknots=150),col="red",lwd=2,type='l') 
axis(side = 4)
mtext("log2(X;V/WT) Expression levels", side = 4, line = 2)
#segments(x0=(b$start),x1=(b$stop),y0=-1.9,y1=-1.9,col='purple',lwd=4) #TSSs with NaN log2 ratios
segments(x0=(c$start),x1=(c$stop),y0=-1.9,y1=-1.9,col='orange',lwd=4) #TSSs with number log2 ratios


dev.off()

