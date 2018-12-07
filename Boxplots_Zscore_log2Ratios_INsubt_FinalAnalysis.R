############################################################################################################
###script to get boxplots of zscore normalized log2 ratio of ChIP signal at peaks or features such as TSS###
############################################################################################################
#written by Lena Street
#las821@nyu.edu
#5/11/18

#Boxplot: boxplot of the Zscored log2 ratio between mutant and wt chip-seq signal at locations of interest (ie. TSS, peaks) for each chromosome
#         bottom of box is 25 percentile, top is 75 percentile
#         center line is median
#         upper whisker extends from the hinge to the largest value no further than 1.5 * IQR from the hinge (where IQR is the inter-quartile range, or distance between the first and third quartiles), the same happens in the other direction for the bottom whisker 
#         outliers are not plotted
#         notch shows the 95% confidence interval of the median (median +/- 1.58*IQR/sqrt(n))
#         p-values from a two-tailed unpaired t-test are plotted above the chromosomes testing if the median ratio for that chromosome is significantly different from the median ratios of the other autosomes

##########input###############################################################################

#the input raw counts file is the output of the deeptools tool multiBigWigSummary
#multiBigwigSummary BED-file --BED /scratch/las821/Meyer_GroSeq_c_elegans.WS220.WS230.gene.annotation.TSS.all.sorted_points_1kb_chr.bed -b H3K4me3_CB428_average.bw -o temp.npz --outRawCounts H3K4me3_CB428_N2_emb_raw_counts.txt
#edit the output txt file to convert CHROMOSOME_ to chr if necessary

########packages###############################################################################

library(ggplot2)

########functions##############################################################################

prep_fusion_XV<- function(x) {
  x$chrFusion <- as.character(x$chr)
  x$chrFusion[which(x$chr=="chrV" & x$start > 17000000)] <- "chrV_Mb17"
  x$chrFusion[which(x$chr=="chrV" & x$start > 18000000)] <- "chrV_Mb18"
  x$chrFusion[which(x$chr=="chrV" & x$start > 19000000)] <- "chrV_Mb19"
  x$chrFusion[which(x$chr=="chrV" & x$start > 20000000)] <- "chrV_Mb20"
  x$chrFusion<-as.factor(x$chrFusion)
  
  return(x)
}

prep_pval_fusion_XV<- function(x) {
  x$bpChrI <- x$chr
  levels(x$bpChrI)[levels(x$bpChrI)=="chrII"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrIII"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrIV"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrV"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrX"] <- NA
  
  x$bpChrII <- x$chr
  levels(x$bpChrII)[levels(x$bpChrII)=="chrI"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrIII"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrIV"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrV"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrX"] <- NA
  
  x$bpChrIII <- x$chr
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrI"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrII"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrIV"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrV"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrX"] <- NA
  
  x$bpChrIV <- x$chr
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrI"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrII"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrIII"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrV"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrX"] <- NA
  
  x$bpChrV <- x$chrFusion
  levels(x$bpChrV)[levels(x$bpChrV)=="chrI"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrII"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrIII"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrIV"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrV_Mb17"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrV_Mb18"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrV_Mb19"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrV_Mb20"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrX"] <- NA
  
  x$bpChrV_17 <- x$chrFusion
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrI"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrII"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrIII"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrIV"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrV"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrV_Mb18"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrV_Mb19"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrV_Mb20"] <- "Aut"
  levels(x$bpChrV_17)[levels(x$bpChrV_17)=="chrX"] <- NA
  
  x$bpChrV_18 <- x$chrFusion
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrI"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrII"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrIII"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrIV"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrV"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrV_Mb17"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrV_Mb19"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrV_Mb20"] <- "Aut"
  levels(x$bpChrV_18)[levels(x$bpChrV_18)=="chrX"] <- NA
  
  x$bpChrV_19 <- x$chrFusion
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrI"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrII"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrIII"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrIV"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrV"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrV_Mb17"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrV_Mb18"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrV_Mb20"] <- "Aut"
  levels(x$bpChrV_19)[levels(x$bpChrV_19)=="chrX"] <- NA
  
  x$bpChrV_20 <- x$chrFusion
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrI"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrII"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrIII"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrIV"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrV"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrV_Mb17"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrV_Mb18"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrV_Mb19"] <- "Aut"
  levels(x$bpChrV_20)[levels(x$bpChrV_20)=="chrX"] <- NA
  
  x$bpChrX <- x$chr
  levels(x$bpChrX)[levels(x$bpChrX)=="chrI"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrII"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrIII"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrIV"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrV"] <- "Aut"
  
  return(x)
}

t_test_fusion_XV <- function(x, cut_off_floor, cut_off_ceiling) {
  a<-x[which(x$log2ratioZscore < cut_off_ceiling & x$log2ratioZscore > cut_off_floor),]
  chrIt<-t.test(log2ratioZscore~bpChrI,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIIt<-t.test(log2ratioZscore~bpChrII,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIIIt<-t.test(log2ratioZscore~bpChrIII,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIVt<-t.test(log2ratioZscore~bpChrIV,data=a,alternative = "two.sided",na.action="na.exclude")
  chrVt<-t.test(log2ratioZscore~bpChrV,data=a,alternative = "two.sided",na.action="na.exclude")
  chrV_17t<-t.test(log2ratioZscore~bpChrV_17,data=a,alternative = "two.sided",na.action="na.exclude")
  chrV_18t<-t.test(log2ratioZscore~bpChrV_18,data=a,alternative = "two.sided",na.action="na.exclude")
  chrV_19t<-t.test(log2ratioZscore~bpChrV_19,data=a,alternative = "two.sided",na.action="na.exclude")
  chrV_20t<-t.test(log2ratioZscore~bpChrV_20,data=a,alternative = "two.sided",na.action="na.exclude")
  chrXt<-t.test(log2ratioZscore~bpChrX,data=a,alternative = "two.sided",na.action="na.exclude")
  temp<-c(chrIt$p.value,chrIIt$p.value,chrIIIt$p.value,chrIVt$p.value,chrVt$p.value,chrV_17t$p.value,chrV_18t$p.value,chrV_19t$p.value,chrV_20t$p.value,chrXt$p.value)
  temp1<-c()
  for(i in 1:length(temp)) {
    if(round(temp[i], digits=3) < 0.001) {temp1[i]<-formatC(temp[i], format = "e", digits = 2)}
    else {temp1[i]<-round(temp[i],digits=3)}
  }
  chrFusion<-c('chrI','chrII','chrIII','chrIV','chrV','chrV_Mb17','chrV_Mb18','chrV_Mb19','chrV_Mb20','chrX')
  ttest<-data.frame(chrFusion=as.factor(chrFusion),pval=temp1)
  return(ttest)
}

box_plot_fusion_XV <- function(y_title,
                               main_title,
                               floor_limit,
                               ceiling_limit,
                               y_scale_tick){
  
  bxplot<-ggplot(rawCounts, aes(x=chrFusion,y=log2ratioZscore,fill=chrFusion)) + geom_boxplot(outlier.shape = NA,notch=FALSE) + 
    theme_classic() + scale_fill_manual(values=c("#3399FF","#6699FF","#3366CC","#3366FF","#3333CC",
                                                 "lightsteelblue1","slategray1","lightsteelblue2","slategray2","#FF3300")) + 
    labs(fill="chr",y=as.character(y_title), title= as.character(main_title)) +
    scale_y_continuous(breaks=seq(floor_limit, ceiling_limit, y_scale_tick), limits=c(floor_limit, ceiling_limit)) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5)+ theme(plot.title = element_text(size=10)) + 
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    geom_text(data = ttest, aes(x = chrFusion, y = ceiling_limit, label = pval)) 
  
  return(bxplot)
}

prep_pval<- function(x) {
  x$bpChrI <- x$chr
  levels(x$bpChrI)[levels(x$bpChrI)=="chrII"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrIII"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrIV"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrV"] <- "Aut"
  levels(x$bpChrI)[levels(x$bpChrI)=="chrX"] <- NA
  
  x$bpChrII <- x$chr
  levels(x$bpChrII)[levels(x$bpChrII)=="chrI"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrIII"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrIV"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrV"] <- "Aut"
  levels(x$bpChrII)[levels(x$bpChrII)=="chrX"] <- NA
  
  x$bpChrIII <- x$chr
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrI"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrII"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrIV"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrV"] <- "Aut"
  levels(x$bpChrIII)[levels(x$bpChrIII)=="chrX"] <- NA
  
  x$bpChrIV <- x$chr
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrI"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrII"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrIII"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrV"] <- "Aut"
  levels(x$bpChrIV)[levels(x$bpChrIV)=="chrX"] <- NA
  
  x$bpChrV <- x$chr
  levels(x$bpChrV)[levels(x$bpChrV)=="chrI"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrII"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrIII"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrIV"] <- "Aut"
  levels(x$bpChrV)[levels(x$bpChrV)=="chrX"] <- NA
  
  x$bpChrX <- x$chr
  levels(x$bpChrX)[levels(x$bpChrX)=="chrI"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrII"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrIII"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrIV"] <- "Aut"
  levels(x$bpChrX)[levels(x$bpChrX)=="chrV"] <- "Aut"
  
  return(x)
}

t_test <- function(x, cut_off_floor, cut_off_ceiling) {
  a<-x[which(x$log2ratioZscore < cut_off_ceiling & x$log2ratioZscore > cut_off_floor),]
  chrIt<-t.test(log2ratioZscore~bpChrI,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIIt<-t.test(log2ratioZscore~bpChrII,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIIIt<-t.test(log2ratioZscore~bpChrIII,data=a,alternative = "two.sided",na.action="na.exclude")
  chrIVt<-t.test(log2ratioZscore~bpChrIV,data=a,alternative = "two.sided",na.action="na.exclude")
  chrVt<-t.test(log2ratioZscore~bpChrV,data=a,alternative = "two.sided",na.action="na.exclude")
  chrXt<-t.test(log2ratioZscore~bpChrX,data=a,alternative = "two.sided",na.action="na.exclude")
  temp<-c(chrIt$p.value,chrIIt$p.value,chrIIIt$p.value,chrIVt$p.value,chrVt$p.value,chrXt$p.value)
  temp1<-c()
  for(i in 1:length(temp)) {
    if(round(temp[i], digits=3) < 0.001) {temp1[i]<-formatC(temp[i], format = "e", digits = 2)}
    else {temp1[i]<-round(temp[i],digits=3)}
  }
  chr<-c('chrI','chrII','chrIII','chrIV','chrV','chrX')
  ttest<-data.frame(chr=as.factor(chr),pval=temp1)
  return(ttest)
}

box_plot <- function(y_title,
                     main_title,
                     floor_limit,
                     ceiling_limit,
                     y_scale_tick){

  bxplot<-ggplot(rawCounts, aes(x=chr,y=log2ratioZscore,fill=chr)) + stat_boxplot(geom ='errorbar', width = 0.4) + geom_boxplot(outlier.shape = NA,notch=TRUE) + 
          theme_classic() + scale_fill_manual(values=c("#3399FF","#6699FF","#3366CC","#3366FF","#3333CC","#FF3300")) + 
          labs(fill="",y=as.character(y_title), title= as.character(main_title)) +
          scale_y_continuous(breaks=seq(floor_limit, ceiling_limit, y_scale_tick), limits=c(floor_limit, ceiling_limit)) +
          geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5)+ theme(plot.title = element_text(size=10)) + 
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
          geom_text(data = ttest, aes(x = chr, y = ceiling_limit, label = pval),size=4)

return(bxplot)
}


#######run################################################################################################

inputFile<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/Histogram_and_boxplots/Plot_input_data/RNAPolII_GW638_N2_emb_raw_counts.txt',header=T)

samples<-colnames(inputFile)
samples

inputData<-inputFile
# inputData<-inputFile[,c("X..chr.","X.start.","X.end.")]
# cols<-samples[grep("H3K27me1",samples)]
# inputData<-cbind(inputData,inputFile[,cols])
# cols
colnames(inputData)<-c("chr","start","end",'N2','GW638')

# #if you need to get the top 1% of data
# inputData<-inputData[which(is.nan(inputData$N2)==FALSE),]
# n<-1
# inputData<-inputData[inputData$N2 > quantile(inputData$N2,prob=1-n/100,na.rm=T),]

#get the log2 ratio
inputData$log2ratio<-log2(inputData$GW638/inputData$N2)

NaNs<-nrow(inputData[which(is.na(inputData$log2ratio)),])
XNaNs<-nrow(inputData[which(is.na(inputData$log2ratio) & inputData$chr=='chrX'),])
a<-inputData[which(is.na(inputData$log2ratio)),]
inputData_1<-inputData[which(!is.na(inputData$log2ratio)),]


# ###
# #if the data includes the XV fusion strain 15eh1
# fusIntermed<-prep_fusion_XV(inputData)
# rawCounts<-prep_pval_fusion_XV(fusIntermed)
# 
# rawCounts$log2ratioZscore<-(rawCounts$log2ratio-mean(rawCounts$log2ratio))/sd(rawCounts$log2ratio)
# 
# ttest<-t_test_fusion_XV(rawCounts, -Inf, Inf)
# 
# bxplot<-box_plot_fusion_XV("Zscore of log2 ratio X;V/WT",
#                  "Zscore of log2 ratio of DPY27 ChIP-seq signal at 1kb around Gro-Seq TSS",-3,8,1)
# print(bxplot)
# ggsave("~/Desktop/DPY27_15eh1_N2_emb_at1kb_around_Gro-seq_TSS_boxplot_log2RatioZscore.pdf", width = 9, height = 7)


###
#if the data does not include fusion strains
rawCounts<-prep_pval(inputData_1)

rawCounts$log2ratioZscore<-(rawCounts$log2ratio-mean(rawCounts$log2ratio))/sd(rawCounts$log2ratio)

ttest<-t_test(rawCounts, -Inf, Inf)

#plot the boxplot
bxplot<-box_plot("Zscore of log2 (GW638 / N2) \nRNA PolII ChIP-seq",
                 "Zscore of log2 ratio of RNA PolII signal at 1kb Gro-seq TSS",-4,4,1)
print(bxplot)
ggsave("~/Desktop/RNApolII_GW638_N2_emb_at_1kbGroTSS_boxplot_log2RatioZscore_INsubt.pdf", width = 5, height = 7)


#200bp around Vector LIN-15 summits
#1kb around Gro-Seq TSS
#options(scipen=-100, digits=4)

# pval<-levels(ttest$pval)[ttest$pval]
# if(is.null(pval)==T){pval<-ttest$pval}
# 
# #summaryFile<-read.delim(file='~/Desktop/Summary_Zscore_log2ratio_INsubt_boxplotData_1kb_around_Gro-seq_TSS.txt',header=T,stringsAsFactors = F)
# #summaryFile<-setNames(data.frame(matrix(ncol = 12, nrow = 0)), c("ChIP","NumStrain","DenStrain","stage","NaNs","XNaNs","chrI","chrII","chrIII","chrIV","chrV","chrX"))
# newRow<-c("H3K9me3","N2 Vector RNAi","N2","emb",NaNs,XNaNs,pval[1],pval[2],pval[3],pval[4],pval[5],pval[6])
# summaryFile[nrow(summaryFile)+1,] <- newRow
# 
# write.table(summaryFile,file="~/Desktop/Summary_Zscore_log2ratio_INsubt_boxplotData_WBTSS.txt",sep="\t",quote=F,col.names=T,row.names=F)

#options(scipen=999)

#get output txt file with the zscores for each location of interest (ie. TSS)
#finalFile<-read.delim(file="~/Desktop/log2ratioZscore_boxplotData_1kb_around_Gro-seq_TSS.txt",header=T)
# outTable<-(inputFile[,1:3])
# outTable<-(inputData[,1:3]) #use when getting top 1%
# colnames(outTable)<-c("chr","start","stop")
a<-prep_pval(a)
a$log2ratioZscore<-NA
b<-rbind(rawCounts,a)
b$index <- as.numeric(row.names(b))
b<-b[order(b$index), ]
outTable$PQN85_dpy27rnai_over_N2_emb<-round(b$log2ratioZscore,digits=3)
write.table(outTable,file="~/Desktop/PQN85_200bpN2embSplitSummit_Zscore_log2ratio_boxplotData.txt",sep="\t",quote=F,col.names=T,row.names=F)

