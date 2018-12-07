inDat<-read.delim(file='~/Desktop/Ercan_Lab/DESeq/VC199_vs_N2/results/deseq.summaryoverview.txt',header=T,stringsAsFactors = F)
a<-rownames(inDat)
inDat$gene<-a
#inDat$log2FC<-log2((inDat[,3]-inDat[,2])/inDat[,2])
colnames(inDat)
inputData<-inDat[,c(18,24)]
colnames(inputData)<-c("log2","gene")

ref<-read.table(file='~/Desktop/Ercan_Lab/Chromatin_Paper/Plot_annotations/c.elegans.WS220.gene.wbid.chr.start.stop.txt',header=T,stringsAsFactors=F)

inputData$chr<-ref$chr[match(rownames(inputData),ref$gene)]
inputData$start<-ref$start[match(rownames(inputData),ref$gene)]
inputData$end<-ref$stop[match(rownames(inputData),ref$gene)]

chrLengths<-c(15072423,15279345,13783700,17493793,20924149,17718866)
rnaSeqWind<-data.frame(chr=c("chrI","chrI","chrI","chrII","chrII","chrII","chrIII","chrIII","chrIII","chrIV","chrIV","chrIV","chrV","chrV","chrV","chrX","chrX","chrX"))
rnaSeqWind$start<-rep(NA,18)
rnaSeqWind$stop<-rep(NA,18)

#0.25Mb windows
for(i in seq(1, nrow(rnaSeqWind), 3)){
  a<-rep(chrLengths, each = 3)
  rnaSeqWind[i,"start"]<-1
  rnaSeqWind[i,"stop"]<-250000
  rnaSeqWind[i+1,"start"]<-((a[i]/2)-125000)
  rnaSeqWind[i+1,"stop"]<-((a[i]/2)+124999)
  rnaSeqWind[i+2,"start"]<-(a[i]-250000)
  rnaSeqWind[i+2,"stop"]<-a[i]
}


#0.5Mb windows
for(i in seq(1, nrow(rnaSeqWind), 3)){
  a<-rep(chrLengths, each = 3)
  rnaSeqWind[i,"start"]<-1
  rnaSeqWind[i,"stop"]<-500000
  rnaSeqWind[i+1,"start"]<-((a[i]/2)-250000)
  rnaSeqWind[i+1,"stop"]<-((a[i]/2)+249999)
  rnaSeqWind[i+2,"start"]<-(a[i]-500000)
  rnaSeqWind[i+2,"stop"]<-a[i]
}

#1Mb windows
for(i in seq(1, nrow(rnaSeqWind), 3)){
  a<-rep(chrLengths, each = 3)
  rnaSeqWind[i,"start"]<-1
  rnaSeqWind[i,"stop"]<-1000000
  rnaSeqWind[i+1,"start"]<-((a[i]/2)-500000)
  rnaSeqWind[i+1,"stop"]<-((a[i]/2)+499999)
  rnaSeqWind[i+2,"start"]<-(a[i]-1000000)
  rnaSeqWind[i+2,"stop"]<-a[i]
}

rnaSeqWind$Chr<-c("chrI_L","chrI_M","chrI_R","chrII_L","chrII_M","chrII_R","chrIII_L","chrIII_M","chrIII_R","chrIV_L","chrIV_M","chrIV_R","chrV_L","chrV_M","chrV_R","chrX_L","chrX_M","chrX_R")

#when using YPT41- 20Kb are removed from the left end of chrII in the fusion
rnaSeqWind[4,"start"]<-20000
rnaSeqWind[4,"stop"]<-270000

#add 0.5 Mb spreading windows in 15eh1
for(i in 19){
  rnaSeqWind[i,"chr"]<-"chrV"
  rnaSeqWind[i,"start"]<-18924149
  rnaSeqWind[i,"stop"]<-19424149
  rnaSeqWind[i,"Chr"]<-"chrV_17"
  rnaSeqWind[i+1,"chr"]<-"chrV"
  rnaSeqWind[i+1,"start"]<-19424149
  rnaSeqWind[i+1,"stop"]<-19924149
  rnaSeqWind[i+1,"Chr"]<-"chrV_18"
  rnaSeqWind[i+2,"chr"]<-"chrV"
  rnaSeqWind[i+2,"start"]<-19924149
  rnaSeqWind[i+2,"stop"]<-20424149
  rnaSeqWind[i+2,"Chr"]<-"chrV_19"
}

#add Mb spreading windows in 15eh1
for(i in 19){
  rnaSeqWind[i,"chr"]<-"chrV"
  rnaSeqWind[i,"start"]<-16924149
  rnaSeqWind[i,"stop"]<-17924149
  rnaSeqWind[i,"Chr"]<-"chrV_17"
  rnaSeqWind[i+1,"chr"]<-"chrV"
  rnaSeqWind[i+1,"start"]<-17924149
  rnaSeqWind[i+1,"stop"]<-18924149
  rnaSeqWind[i+1,"Chr"]<-"chrV_18"
  rnaSeqWind[i+2,"chr"]<-"chrV"
  rnaSeqWind[i+2,"start"]<-18924149
  rnaSeqWind[i+2,"stop"]<-19924149
  rnaSeqWind[i+2,"Chr"]<-"chrV_19"
}


#add 0.5 Mb spreading windows YPT41
for(i in 19){
  rnaSeqWind[i,"chr"]<-"chrII"
  rnaSeqWind[i,"start"]<-520000
  rnaSeqWind[i,"stop"]<-1020000
  rnaSeqWind[i,"Chr"]<-"chrII_2"
  rnaSeqWind[i+1,"chr"]<-"chrII"
  rnaSeqWind[i+1,"start"]<-1020000
  rnaSeqWind[i+1,"stop"]<-1520000
  rnaSeqWind[i+1,"Chr"]<-"chrII_3"
  rnaSeqWind[i+2,"chr"]<-"chrII"
  rnaSeqWind[i+2,"start"]<-1520000
  rnaSeqWind[i+2,"stop"]<-2020000
  rnaSeqWind[i+2,"Chr"]<-"chrII_4"
}


#add Mb spreading windows YPT41
for(i in 19){
  rnaSeqWind[i,"chr"]<-"chrII"
  rnaSeqWind[i,"start"]<-1020000
  rnaSeqWind[i,"stop"]<-2020000
  rnaSeqWind[i,"Chr"]<-"chrII_2"
  rnaSeqWind[i+1,"chr"]<-"chrII"
  rnaSeqWind[i+1,"start"]<-2020000
  rnaSeqWind[i+1,"stop"]<-3020000
  rnaSeqWind[i+1,"Chr"]<-"chrII_3"
  rnaSeqWind[i+2,"chr"]<-"chrII"
  rnaSeqWind[i+2,"start"]<-3020000
  rnaSeqWind[i+2,"stop"]<-4020000
  rnaSeqWind[i+2,"Chr"]<-"chrII_4"
}

rnaSeqWind<-rnaSeqWind[order(rnaSeqWind[,1], rnaSeqWind[,2]),]

dataToPlot<-c()
for(i in 1:nrow(rnaSeqWind)){
    chrom<-rnaSeqWind$chr[i]
    start<-rnaSeqWind$start[i]
    end<-rnaSeqWind$stop[i]
    temp<-inputData[which((inputData$chr == chrom & inputData$start >=start & inputData$start <=end) | (inputData$chr == chrom & inputData$end >=start & inputData$end <=end)),]
    chr<-rep(rnaSeqWind[i,"Chr"],nrow(temp))
    temp$Chr<-chr
    out<-cbind(temp$Chr,temp$log2,temp$chr)
    dataToPlot<-rbind(dataToPlot,out)
  }
dataToPlot<-as.data.frame(dataToPlot)
colnames(dataToPlot)<-c("Chr","log2","chr")
dataToPlot<-dataToPlot[which(!is.na(dataToPlot$log2)),]
#dataToPlot$Chr<-as.factor(dataToPlot$Chr)
dataToPlot$log2<-as.numeric(as.character(dataToPlot$log2))


#get p vaues
up<-c()
down<-c()
pval<-c()
all.up<-length(which(inputData$log2 > 0))
all.down<-length(which(inputData$log2 < 0))
for(i in 1:nrow(rnaSeqWind)){
  chrom<-rnaSeqWind$chr[i]
  start<-rnaSeqWind$start[i]
  end<-rnaSeqWind$stop[i]
  temp<-inputData[which((inputData$chr == chrom & inputData$start >=start & inputData$start <=end) | (inputData$chr == chrom & inputData$end >=start & inputData$end <=end)),]
  up[i]<-length(which(temp$log2 > 0))
  down[i]<-length(which(temp$log2 < 0))
  mat<-matrix(c(up[i],down[i], all.up,all.down),ncol=2,byrow=T)
  pval[i]<-as.numeric(fisher.test(mat,alternative='two.sided')$p.value)
}

for(i in 1:length(pval)) {
  numb<-as.numeric(pval[i])
  if(round(numb, digits=3) < 0.001) {pval[i]<-formatC(numb, format = "e", digits = 2)}
  else {pval[i]<-round(numb,digits=3)}
}

rnaSeqWind$up<-up
rnaSeqWind$down<-down
rnaSeqWind$pval<-pval


median(dataToPlot$log2)

#set factor levels for non-fusion sample
dataToPlot$Chr <- factor(dataToPlot$Chr, levels = c("chrI_L","chrI_M","chrI_R","chrII_L","chrII_M","chrII_R","chrIII_L","chrIII_M","chrIII_R","chrIV_L","chrIV_M","chrIV_R","chrV_L","chrV_M","chrV_R","chrX_L","chrX_M","chrX_R"))

#set factor levels for 15eh1 spreading
dataToPlot$Chr <- factor(dataToPlot$Chr, levels = c("chrI_L","chrI_M","chrI_R","chrII_L","chrII_M","chrII_R","chrIII_L","chrIII_M","chrIII_R","chrIV_L","chrIV_M","chrIV_R","chrV_L","chrV_M","chrV_17","chrV_18","chrV_19","chrV_R","chrX_L","chrX_M","chrX_R"))

#set factor levels for YPT41 spreading
dataToPlot$Chr <- factor(dataToPlot$Chr, levels = c("chrI_L","chrI_M","chrI_R","chrII_L","chrII_2","chrII_3","chrII_4","chrII_M","chrII_R","chrIII_L","chrIII_M","chrIII_R","chrIV_L","chrIV_M","chrIV_R","chrV_L","chrV_M","chrV_R","chrX_L","chrX_M","chrX_R"))

#plot the boxplot
library("ggplot2")
box_plot <- function(y_title,
                     main_title,
                     floor_limit,
                     ceiling_limit,
                     y_scale_tick){
  
  #xTicks<-c("left 250 Kb","middle 250 Kb","right 250 Kb")
  #xTicks<-c("left Mb","middle Mb","right Mb")
  xTicks<-c("left 500Kb","middle 500Kb","right 500Kb")
  #xTicks<-c("left 500Kb","middle 500Kb","right 500Kb","1 500Kb","2 500Kb","3 500Kb","4 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb")
  #xTicks<-c("left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","right 500Kb","left 500Kb","middle 500Kb","39 500Kb","40 500Kb","41 500Kb","42 500Kb","left 500Kb","middle 500Kb","right 500Kb")
  #xTicks<-c("left Mb","middle Mb","right Mb","1 Mb","2 Mb","3 Mb","4 Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb")
  #xTicks<-c("left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","right Mb","left Mb","middle Mb","17 Mb","18 Mb","19 Mb", "20 Mb","left Mb","middle Mb","right Mb")
  xTicks<-rep(xTicks,6)
  bxplot<-ggplot(dataToPlot, aes(x=Chr,y=log2,fill=chr)) + stat_boxplot(geom ='errorbar', width = 0.5,lwd=0.3) + geom_boxplot(lwd=0.3,outlier.shape = NA) + 
    theme_classic() + scale_fill_manual(values=c("#3399FF","#6699FF","#3366CC","#3366FF","#3333CC","#FF3300")) + 
    labs(fill="chr",y=as.character(y_title), title= as.character(main_title)) +
    scale_y_continuous(breaks=seq(floor_limit, ceiling_limit, y_scale_tick), limits=c(floor_limit, ceiling_limit)) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5)+ theme(plot.title = element_text(size=10)) + 
    scale_x_discrete(labels=xTicks)+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
    geom_text(data = rnaSeqWind, aes(x = Chr, y = (ceiling_limit-0.05), label = pval),size=3,angle = 90) 
  
  return(bxplot)
}

bxPlot<-box_plot("log2(VC199 L3 counts/N2 L3 counts)","",-2.5,3,0.5)
print(bxPlot)

ggsave("~/Desktop/RNAseq_VC199_vs_N2_L3_chr500Kb_boxplot.pdf", width = 10, height = 8)



