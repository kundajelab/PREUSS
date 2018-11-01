rm(list=ls())
library(gplots)
library(dendextend)
library(circlize)


#####INPUTS########
adar_editing_levels="../data_from_xin/Neil1_Editing_Level_For_Anna_2018-01-23.csv"

#FOR STRUCTURE IDENTITY
#distance_matrix="web_beagle_alignment_clustering/str_id/beagle.str_id.txt"
#output_image="web_beagle_dendrogram_str_id.png"

#FOR SEQUENCE IDENTITY
distance_matrix="web_beagle_alignment_clustering/seq_id/beagle.seq_id.txt"
output_image="web_beagle_dendrogram_seq_id.png"

#FOR STRUCTURE SIMILARITY
#distance_matrix="web_beagle_alignment_clustering/str_sim/beagle.str_sim.txt"
#output_image="web_beagle_dendrogram_str_sim.png"
###################


data=read.table(distance_matrix,header=TRUE,sep='\t',row.names=1)
names(data)=rownames(data)
data=na.omit(data)
data=as.matrix(data)

# create a dendrogram
hc <- hclust(dist(t(data)))
dend <- as.dendrogram(hc)
hc_labels=as.numeric(sapply(strsplit(labels(dend),"-"), `[`, 1))

#find the corresponding adar editing levels 
#order to match hc_labels 
adar_edit=read.table(adar_editing_levels,sep='\t',na.strings = "NA",header=TRUE)
adar_edit=adar_edit[hc_labels,]


png(output_image,height=10,width=10,units='in',res=300)
col_fun = colorRamp2(c(
                       min(na.omit(adar_edit$Ave_Editing_Level)), 
                       mean(na.omit(adar_edit$Ave_Editing_Level)), 
                       max(na.omit(adar_edit$Ave_Editing_Level))), 
                     c("blue", "yellow", "red"))
circos.par(cell.padding = c(0, 0, 0, 0),track.margin=c(0.01,0))
circos.initialize('a', xlim = c(0,292))

circos.track(ylim=c(0,20),track.height=0.3,panel.fun=function(x,y){
  circos.text(1:292-1,rep(0,292),labels(dend),facing="reverse.clockwise",adj=c(1,0),niceFacing=TRUE,cex=0.4)
},bg.border=NA)

adar_edit[is.na(adar_edit)]=-10
circos.track(ylim=c(0,1),track.height=0.2, bg.border = NA, panel.fun = function(x, y) {
  sector.index = CELL_META$sector.index
  m = t(adar_edit$Ave_Editing_Level)
  names(m)="Ave_Editing_Level"
  col_mat = col_fun(as.vector(m))
  col_mat[m==-10]="#000000"
  nr = nrow(m)
  nc = ncol(m)
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                                 1:nc, rep(nr - i + 1, nc), 
                                 border = col_mat, col = col_mat)
  }
})

max_height = attr(dend, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend, max_height = max_height)
},  bg.border = NA,track.height=0.3)
vals=seq(0,1,0.1)
legend(-.1, .3, legend=c('NA',as.character(vals)),
       col=c("#000000",col_fun(vals)), lty=1,cex=0.6)
circos.clear()
dev.off()
