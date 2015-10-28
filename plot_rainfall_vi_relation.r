
#need to run after residual analysis
# dataset used as follows:
#calculate correlation and p value
tcw.rf.cor2 = apply(tcw.1.df, 1, COR.test)
tcb.rf.cor2 = apply(tcb.1.df, 1, COR.test)
tcg.rf.cor2 = apply(tcg.1.df, 1, COR.test)
tca.rf.cor2 = apply(tca.1.df, 1, COR.test)
evi.rf.cor2 = apply(evi.1.df, 1, COR.test)
ndwi.rf.cor2 = apply(ndwi.1.df, 1, COR.test)

#set class and color
rel_arg1 <- list(at=seq(1,8,1), labels=c("< 0","0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5","0.5 - 0.6","> 0.6")) #these are the class names
rel_labels=c("< 0","0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5","0.5 - 0.6","> 0.6")
library(RColorBrewer)
rel_color1=brewer.pal(8,"YlGnBu")

#set threshold
threshold =  c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

#for TCW
tcw.rf.cor.grd = Point2raster(tcw.rf.cor2[1,], TCW2.dry[[1]])
tcw.rf.cor.grd.sig = Point2raster(tcw.rf.cor2[2,], TCW2.dry[[1]])
tcw.rf.cor.grd2 = recls2(tcw.rf.cor.grd, threshold = threshold)

#for TCB
tcb.rf.cor.grd = Point2raster(tcb.rf.cor2[1,], TCW2.dry[[1]])
tcb.rf.cor.grd.sig = Point2raster(tcb.rf.cor2[2,], TCW2.dry[[1]])
tcb.rf.cor.grd2 = recls2(tcb.rf.cor.grd, threshold = threshold)

#for TCG
tcg.rf.cor.grd = Point2raster(tcg.rf.cor2[1,], TCW2.dry[[1]])
tcg.rf.cor.grd.sig = Point2raster(tcg.rf.cor2[2,], TCW2.dry[[1]])
tcg.rf.cor.grd2 = recls2(tcg.rf.cor.grd, threshold = threshold)

#for TCA
tca.rf.cor.grd = Point2raster(tca.rf.cor2[1,], TCW2.dry[[1]])
tca.rf.cor.grd.sig = Point2raster(tca.rf.cor2[2,], TCW2.dry[[1]])
tca.rf.cor.grd2 = recls2(tca.rf.cor.grd, threshold = threshold)

#for EVI
evi.rf.cor.grd = Point2raster(evi.rf.cor2[1,], TCW2.dry[[1]])
evi.rf.cor.grd.sig = Point2raster(evi.rf.cor2[2,], TCW2.dry[[1]])
evi.rf.cor.grd2 = recls2(evi.rf.cor.grd, threshold = threshold)

#for NDWI
ndwi.rf.cor.grd = Point2raster(ndwi.rf.cor2[1,], TCW2.dry[[1]])
ndwi.rf.cor.grd.sig = Point2raster(ndwi.rf.cor2[2,], TCW2.dry[[1]])
ndwi.rf.cor.grd2 = recls2(ndwi.rf.cor.grd, threshold = threshold)

corr = stack(tcb.rf.cor.grd2,tcg.rf.cor.grd2,tcw.rf.cor.grd2,tca.rf.cor.grd2,evi.rf.cor.grd2,ndwi.rf.cor.grd2)
names(corr) <- c("TCB", "TCG","TCW","TCA","EVI","NDWI")

#raster plot methods
png(file = ".\\NBAR_results3\\rainfall_vi_relationship.png", width = 4000, height = 3000, units = "px", res = 300)

par(mfrow=c(3,2),mar=c(0,0,0,0))

for (i in 1:6){
plot(corr[[i]],col = rel_color1, 
     #axis.arg=arg
     #xlim=c(-2500000, -500000), 
     #ylim=c(ext@ymin, ext@ymax),
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)

text(x=-5.5, y=9.75, names(trd)[i], cex = 2)
plot(county_b_ghana, add = TRUE)
plot(regions, border = "red", add = TRUE)

}

legend(x = -13.5, y = 7.5, legend = rel_labels, fill = rel_color1,  cex = 2,  box.lwd = 0, box.col = "white",bg = "transparent")

dev.off()

#rastervis plot
breaks2_change <- 0:8        
legendbrks2_change <- 1:8 - 0.5

png(file = ".\\NBAR_results3\\rainfall_vi_relationship-2.png", width = 4000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(corr,
          maxpixels = nrow(corr)*ncol(corr),
          at= breaks2_change, margin=FALSE,
          col.regions= rel_color1,
          colorkey= list(labels= list(labels= rel_labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions, col = "red", lwd = 1.5)) 

dev.off()



