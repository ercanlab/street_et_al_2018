import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
import pandas as pd

#To plot the correlation of the Zscore log2(mut/wt) values

#input data are the Zscore values calculated by the boxplots script
inputData = pd.read_table('/Users/lenastreet/Desktop/Ercan_Lab/Chromatin_Paper/Histogram_and_boxplots/Zscore_log2Ratio_INsubt/Zscore_values/CB428_Comparison_1kbGroTSS_Zscore_log2ratio.txt')

#to plot X chromosome data
forPlot=inputData.loc[inputData.chr == 'chrX']

#to plot autosomal data
forPlot=inputData.loc[inputData.chr != 'chrX']

#select only the data you want to plot
forPlot = forPlot.drop(forPlot.columns[[0,1,2,5,8,10,13,14]], axis=1)

#calculate the spearman rank correlations
corr=round(forPlot.corr(method='spearman'),2)

#plot the spearman rank correlation values as a heatmap with a scale from -1 to 1
ax = sns.clustermap(corr, linewidth=0.5,cmap="RdYlBu",annot=True, annot_kws={"size": 7},vmin=-1, vmax=1,linecolor='black')
plt.savefig("/Users/lenastreet/Desktop/CB428_over_N2_Zscorelog2ratio_AutData_ChIPenrich_at_1kbGroTSS_SpearmanCorr_heatmap.pdf",bbox_inches='tight')
plt.show()