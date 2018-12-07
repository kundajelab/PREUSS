rm(list=ls())
library(ggplot2)
library(onehot)
library(Hmisc)
library(corrplot)
library(gplots)

library(optparse)
options(warn=-1)

## Read in user arguments ##
option_list=list(
  make_option(c("--feature_mat",type="character",default=NULL,help="Path to feature matrix")),
  make_option(c("--output_prefix",type="character",default=NULL,help="Directory to store output files"))
  )
opt_parser=OptionParser(option_list=option_list)
opt=parse_args(opt_parser)


if (is.null(opt$feature_mat)){
  print_help(opt_parser)
  stop("Please supply feature matrix with --feature_mat flag", call.=FALSE)
}
if (is.null(opt$output_prefix)){
  print_help(opt_parser)
  stop("Please supply output prefix with --output_prefix flag", call.=FALSE)
}

df=read.table(opt$feature_mat,header=TRUE,sep='\t')
df$cur_id=NULL
df$source=NULL
df$editing_level=NULL

encoder=onehot(df)
onehot_df <- predict(encoder,df)
print("one-hot-encoded factors of data matrix:")
print(nrow(onehot_df))
print(ncol(onehot_df))


#compute correlation matrix (spearman)
spearman_cor <- rcorr(as.matrix(onehot_df),type=c("spearman"))
write.table(spearman_cor$r,file=paste(opt$output_prefix,"SpearmanCorrelation.tsv",sep=""),col.names=TRUE,row.names=TRUE,sep='\t')
#compute correlation matrix (pearson)
pearson_cor <- rcorr(as.matrix(onehot_df),type=c("pearson"))
write.table(pearson_cor$r,file=paste(opt$output_prefix,"PearsonCorrelation.tsv",sep=""),col.names=TRUE,row.names=TRUE,sep='\t')

#Plot Correlation Heatmap 
# Insignificant correlation are crossed
svg(paste(opt$output_prefix,"SpearmanCorrNEIL1.svg",sep=""),height=15,width=15)
print(corrplot(spearman_cor$r, type="upper", order="hclust", 
         p.mat = spearman_cor$P, sig.level = 0.01, insig = "blank", tl.col = "black", tl.srt = 45,tl.cex = 0.5,title = "Spearman Correlation of Features for NEIL1"))
dev.off() 

svg(paste(opt$output_prefix,"PearsonCorrNEIL1.svg",sep=""),height=15,width=15)
print(corrplot(pearson_cor$r, type="upper", order="hclust", 
         p.mat = pearson_cor$P, sig.level = 0.01, insig = "blank", tl.col = "black", tl.srt = 45,tl.cex = 0.5,title = "Pearson Correlation of Features for NEIL1"))
dev.off() 

#Generate heatmap 
col<- colorRampPalette(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))(100)
svg(paste(opt$output_prefix,"SpearmanHeatmap.svg",sep=""),height=15,width=15)
print(heatmap.2(as.matrix(spearman_cor$r),
          col=col,
          Rowv=TRUE,
          Colv=TRUE,
          scale="none",
          trace="none",
          margins = c(15,15)
          ))
dev.off() 

svg(paste(opt$output_prefix,"PearsonHeatmap.svg",sep=""),height=15,width=15)
print(heatmap.2(as.matrix(pearson_cor$r),
                col=col,
                Rowv=TRUE,
                Colv=TRUE,
                scale="none",
                trace="none",
                margins=c(15,15)
))
dev.off() 
