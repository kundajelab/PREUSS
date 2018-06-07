rm(list=ls())
data=read.table("beagle.str_sim.txt",header=TRUE,sep='\t',row.names=1)
names(data)=rownames(data)
data=na.omit(data)
library(gplots)
data=as.matrix(data)
png("web.beagle.png",width=15,height=6,units='in',res=600)
heatmap.2(data,
            col=redblue(30),
            symm=FALSE,
            trace='none',
            Colv=TRUE,
            Rowv=TRUE,
            cexCol=0.3,
            cexRow=0.3,
            density.info="none",
            na.color = '#000000')
dev.off()
# install.packages("dendextend")
# install.packages("circlize")
library(dendextend)
library(circlize)

# create a dendrogram
hc <- hclust(dist(t(data)))
dend <- as.dendrogram(hc)

png("web_beagle_dendrogram.png",height=14,width=14,units='in',res=400)
col_fun = colorRamp2(c(min(data), mean(data), max(data)), c("blue", "yellow", "red"))
circos.par(cell.padding = c(0, 0, 0, 0),track.margin=c(0.01,0))
circos.initialize('a', xlim = c(0,292))

circos.track(ylim=c(0,20),track.height=0.01,panel.fun=function(x,y){
  circos.text(1:292-0.5,rep(0,292),labels(dend),facing="reverse.clockwise",adj=c(1,0),niceFacing=TRUE,cex=0.4)
},bg.border=NA)

circos.track(ylim=c(0,10),track.height=0.2, bg.border = NA, panel.fun = function(x, y) {
  sector.index = CELL_META$sector.index
  m = data
  dend = dend
  m2 = m[, order.dendrogram(dend)]
  col_mat = col_fun(m2)
  nr = nrow(m2)
  nc = ncol(m2)
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                1:nc, rep(nr - i + 1, nc), 
                border = col_mat[i, ], col = col_mat[i, ])
  }
})

max_height = attr(dend, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend, max_height = max_height)
},  bg.border = NA,track.heigh=0.5)

circos.clear()
dev.off()
