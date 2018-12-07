library(reshape2)
library(ggplot2)
library(scales)
library(corr)

inputData<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/Histogram_and_boxplots/Zscore_log2Ratio_INsubt/Zscore_values/CB428_Comparison_1kbGroTSS_Zscore_log2ratio.txt',header=T, stringsAsFactors = F)
inputData<-inputData[which(inputData$chr=='chrX'),]
dpy27testDataAll<-inputData[,c(4:16)]
dpy27testDataN2<-inputData[,c(4,6,7,8,9,11,13,15,16)]
dpy27testDataVect<-inputData[,c(5,6,7,8,9,10,12,14,16)]
CB428testData<-inputData[,c(4:16)]
CB428testData_1<-inputData[,c(4,5,7,8,10,12,13,16)]

write.table(dpy27testDataAll,file='~/Desktop/Ercan_Lab/Chromatin_Paper/Correlation_heatmaps/dpy27testDataAll_wholeGenome.txt',quote = F,row.names =F,sep='\t',col.names = T)
write.table(dpy27testDataN2,file='~/Desktop/Ercan_Lab/Chromatin_Paper/Correlation_heatmaps/dpy27testDataN2_wholeGenome.txt',quote = F,row.names =F,sep='\t',col.names = T)
write.table(dpy27testDataVect,file='~/Desktop/Ercan_Lab/Chromatin_Paper/Correlation_heatmaps/dpy27testDataVect_wholeGenome.txt',quote = F,row.names =F,sep='\t',col.names = T)
write.table(CB428testData_1,file='~/Desktop/Ercan_Lab/Chromatin_Paper/Correlation_heatmaps/CB428testData_X_1.txt',quote = F,row.names =F,sep='\t',col.names = T)



Corrdpy27testDataAll<-round(cor(CB428testData_1,use = "everything",method = 'spearman'),2)

#write.matrix(Corrdpy27testDataAll,file='~/Desktop/test.mat',sep=' ')

melted_Corrdpy27testDataAll<-melt(Corrdpy27testDataAll)

dpy27testDataAllHeatmap<-ggplot(data = melted_Corrdpy27testDataAll, aes(x=Var1, y=Var2,fill=value)) +
  geom_tile(color = "black") + #scale_fill_brewer(palette=13) +
  #scale_fill_gradientn(colours=c("firebrick", "red2", "orange","yellow","lightblue", "blue", "darkblue"),values=rescale(c(-1,-0.6,-0.3,0,0.3, 0.6, 1)), guide="colorbar")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "lightgoldenrod1", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
  theme_minimal() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
print(dpy27testDataAllHeatmap)



dpy27testDataAll_1<-as.matrix(inputData[,c(5,9,10,18,19,26)])

