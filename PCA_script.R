
##コマンドラインからバッチ処理として起動した場合に引数が渡せるか？
#  $ Rscript script.r arg1 arg2 arg3 ...

##svg出力をなんとかできないか？
#　→svg出力はillustratorでは開くとおかしいが、webブラウザでは普通に開ける
# 　→Post Scriptだと

##出力ファイルの出力先は？
# →デフォルトではgetwd() で得られるディレクトリに保存される。
##load the tsv file

#path <- commandArgs(TRUE)
path <- '/Users/cellbiology/Desktop/PESI_DATA/Non_tumor_tissues/funalyzer_140530_non_tumor_liver/non_tumor_liver(59)mz(10.0-1999.8).tsv'
raw_data <- read.csv(path,sep='\t',skip=1)


##pure intensity, except for  'm/z' column
samples <- raw_data[3:19894, colnames(raw_data) != "m.z"]

##string vector of mz column
mz <- as.character(raw_data$m.z)

#compress data; bin width 0.1 to 1.0
n <- length(samples)
start_index <- 4  #start from m/z 10.4
end_index <- 19894 #end at m/z 1999.5

for (i in 1:n){
  for (j in 1:1989){
    
  }
  
  
}


##execute pca using "prcomp" function with the data scaled
is_scaled <- TRUE
pca <- prcomp(t(samples), scale = is_scaled)
pc_summary <- summary(pca)

##set variables
sdev = pca$sdev
loading <- as.data.frame(pca$rotation) * pca$sdev
#loading["PC1"] <- loading$PC1 / sd(samples$X1) 
PC <- as.data.frame(pca$x) #dataframe made of all principal components
PC1 <- PC$PC1
PC2 <- PC$PC2
PC3 <- PC$PC3

PC1_2 <- PC[c("PC1","PC2")] #PC1 vs PC2
PC1_3 <- PC[c("PC1","PC3")] #PC1 vs PC3
PC2_3 <- PC[c("PC2","PC3")] #PC2 vs PC3


##plotting PC plot; the first 20 samples are arbitrarily chosen from 50 samples
num_start <- 1
num_end <- 20


cairo_ps('PC1_2.ps')
plot(PC1_2[num_start:num_end,], xlab=sprintf("PC1 (%.1f %s)", pc_summary$importance[2,1]*100,"%" ),ylab=sprintf("PC2 (%.1f %s) ", pc_summary$importance[2,2]*100,"%"))
title(main='PCA plot: PC1 vs PC2')
text(PC1_2[num_start:num_end,1],PC1_2[num_start:num_end,2],pos=3,cex=0.5)
dev.off()


cairo_ps('PC1_3.ps')
plot(PC1_3[num_start:num_end,], xlab=sprintf("PC1 (%.1f %s)", pc_summary$importance[2,1]*100,"%" ),ylab=sprintf("PC3 (%.1f %s) ", pc_summary$importance[2,3]*100,"%"))
title(main='PCA plot: PC1 vs PC3')
text(PC1_3[num_start:num_end,1],PC1_3[num_start:num_end,2],pos=3,cex=0.5)
dev.off()


cairo_ps('PC2_3.ps')
plot(PC2_3[num_start:num_end,],xlab=sprintf("PC2 (%.1f %s) ", pc_summary$importance[2,2]*100,"%"),ylab=sprintf("PC3 (%.1f %s) ", pc_summary$importance[2,3]*100,"%"))
title(main='PCA: plot: PC2 vs PC3')
text(PC2_3[num_start:num_end,1],PC2_3[num_start:num_end,2],pos=3,cex=0.5)
dev.off()

##plotting loading plot of PC1 ;
#svg(filename='loading.svg')
png(filename='loading.png',width=1500, height=500)
#cairo_ps('loading.ps',width=1500, height=500)
plot(rowSums (samples, na.rm = FALSE, dims = 1),type='l',axes=FALSE, xlab = "",ylab="",col="brown")
axis(4)
par(new=T)
plot(loading$PC1,type='l',main='Loading plot of PC1' ,xlab="",ylab="Factor loading",col="blue")
legend("topright",legend=c("Factor Loading","Intensity"),lty=c(1,0),fill=c(NA,"purple"),col=c("light blue","purple"), bg="gray90")
#mtext("m/z", side = 1,line = 3,family="Times-Italic", cex=1.5)
mtext("m/z", side = 1,line = 3, cex=1.5)
mtext("Intensity", side=4, line=3)
dev.off()



##scree plotに書き直す
##plotting scree plot
#svg(filename='scree_plot.svg')
#png(filename='proportion.png')
cairo_ps('scree_plot.ps')
plot(pc_summary$sdev**2,type='b', main='Scree plot', xlab="Component number",ylab="Eigenvalue",col="dark green")
dev.off()

##plotting proportion of variance
svg(filename='proportion.svg')
cairo_ps('proportion.ps')
plot(pc_summary$importance[2,]*100,type='b',ylim=range(0,100), main='Proportion', xlab="Component number",ylab="Propotion(%)",col="dark red")
par(new=T)
plot(pc_summary$importance[3,]*100,ylim=range(0,100),type='b',col="purple",xlab='',ylab='')
dev.off()

##plotting cummulative proportion
#svg('cumulative_proportion.svg')
#png('cumulative_proportion.png')
cairo_ps('cumulative_proportion.ps')
plot(pc_summary$importance[3,]*100,type='b',ylim=range(0,100),main='Cumulutive Proportion',xlab="i th principal component",ylab="Propotion(%)",col="purple")
dev.off()


dev.off()
dev.off()
dev.off()


