########packages###############################################################################

library(ggplot2)
library(ggpmisc)

########functions##############################################################################

# scatter_plot <- function(df, x, y, x_title, y_title, main_title, x_min, x_max, y_min, y_max, axis_scale_tick){
#   
#   my.formula <- df$y ~ df$x
#   
#   splot<-ggplot(df, aes_string(x=x, y=y)) + geom_point(col="darkgrey") + theme_classic() + 
#     geom_smooth(method=lm,aes(fill="chr"),fullrange=TRUE) + 
#     stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
#     labs(x=as.character(x_title),y=as.character(y_title),
#          title=as.character(main_title)) + theme(plot.title = element_text(size=10)) +
#     theme(legend.position="none") + scale_fill_manual(values=c("#3399FF")) +
#     scale_y_continuous(breaks=seq(y_min, y_max, axis_scale_tick), limits=c(y_min, y_max)) +
#     scale_x_continuous(breaks=seq(x_min, x_max, axis_scale_tick), limits=c(x_min, x_max)) #+
#     #geom_hline(yintercept=0, linetype="dashed", color = "black") +
#     #geom_vline(xintercept=0, linetype="dashed", color = "black") +
#     
#   return(splot)
# }

#######run################################################################################################
options(scipen=999)

##Fig 2B
inputFile<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/DPY27_vs_GeneExpression_plots/N2_Gro_Expression_Data_ChrX_DCCbinding.txt',header=T)

forPlot<-inputFile[,c(7,8,9,19)]
colnames(forPlot)<-c("exp",'explog10','VectRNAiexp','dpy27')

forPlot$VectRNAiexpLog10<-log10(forPlot$VectRNAiexp)

splot<-scatter_plot(forPlot,explog10,dpy27,
                    "Expression",
                    "DPY27 Binding",
                    "",
                    -10,10,-10,10,2)
print(splot)
ggsave("~/Desktop/scatter_log2ratio_H3K27ac_DPY27_DPY27RNAi_CB428_N2_emb_at1kb_around_Gro-seq_TSS.pdf", width = 5, height = 5)

my.formula <- forPlot$dpy27 ~ forPlot$explog10

splot<-ggplot(forPlot, aes(explog10, dpy27)) + geom_point(col="darkgrey") + theme_classic() + 
  xlim(0,2.5) + ylim(0,30)+
  geom_smooth(method=lm,aes(fill="chr"),fullrange=TRUE) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..rr.label..)), parse = TRUE) +
  #stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
  labs(x="Expression \nLog10 (GRO-seq gene body count Kruesi et al.)",y="DPY-27 binding \n(Average ChIP-seq -200 bp-TSS)",title="") + theme(plot.title = element_text(size=10)) +
  theme(legend.position="none") + scale_fill_manual(values=c("#3399FF")) #+
  #scale_y_continuous(breaks=seq(0, 30, 10), limits=c(0, 30)) +
  #scale_x_continuous(breaks=seq(1, 100, 10), limits=c(1, 100))
print(splot)

ggsave("~/Desktop/dpy27Binding_vs_LOG10groGeneBodyCountExp_scatter.pdf", width = 7, height = 5)



################################################################################
################################################################################

##Fig 2D

inputFile<-read.delim(file='~/Desktop/Ercan_Lab/Chromatin_Paper/DPY27_vs_GeneExpression_plots/N2_Gro_Expression_Data_ChrX_DCCbinding.txt',header=T)

forPlot2<-inputFile[,c(32,33)]
colnames(forPlot2)<-c("L3_divBy_emb_GroFPKM",'L3_divBy_emb_dpy27_N2')
forPlot2$log2FPKM<-log2(forPlot2$L3_divBy_emb_GroFPKM)
forPlot2$log2DPY27<-log2(forPlot2$L3_divBy_emb_dpy27_N2)

my.formula <- forPlot2$log2DPY27 ~ forPlot2$log2FPKM

splot<-ggplot(forPlot2, aes(log2FPKM, log2DPY27)) + geom_point(col="darkgrey") + theme_classic() + 
  xlim(-4,4) + ylim(-4,4)+
  geom_smooth(method=lm,aes(fill="chr"),fullrange=TRUE,se=T) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..rr.label..)), parse = TRUE) +
  #stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
  labs(x="Change in expression \nLog2 L3/Emb (GRO-seq gene body count Kruesi et al.)",y="Change in DPY-27 binding \nLog2 L3/Emb (Avg. ChIP-seq -200 bp-TSS)",title="") + 
  theme(axis.title=element_text(size=20),axis.text.y = element_text(size=14),axis.text.x = element_text(size=14)) +
  theme(legend.position="none") + scale_fill_manual(values=c("#3399FF")) #+
  #geom_hline(yintercept=0, linetype="dashed", color = "black") +
  #geom_vline(xintercept=0, linetype="dashed", color = "black")

print(splot)

ggsave("~/Desktop/Change_dpy27Binding_vs_LOG10groGeneBodyCountExp_scatter.pdf", width = 9, height = 7)

