
clasnames <- c("Water",
               "RainForest",
               "Woody Savannas",
               "Savannas",
               "Grasslands",
               "Cropland/Natural Mosaic",
               "Others")

#after run calculating all indices and plot trend

lc_rc1 = lc_rc*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-29,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 4)

lc_rc2 = lc_rc*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-29,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 4)

lc_rc3 = lc_rc*100 + NDWI.trd2.grd 
lc_rc3.df = data.frame(freq(lc_rc3))[-29,]
lc_rc3.df$lc = floor(lc_rc3.df$value/100)
lc_rc3.df$trend = lc_rc3.df$value - lc_rc3.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc3.df, FUN = "sum")
lc_rc3.df$prop = 100*lc_rc3.df$count/rep(lc_sum$count, each = 4)

lc_rc4 = lc_rc*100 + TCG.trd2.grd 
lc_rc4.df = data.frame(freq(lc_rc4))[-29,]
lc_rc4.df$lc = floor(lc_rc4.df$value/100)
lc_rc4.df$trend = lc_rc4.df$value - lc_rc4.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc4.df, FUN = "sum")
lc_rc4.df$prop = 100*lc_rc4.df$count/rep(lc_sum$count, each = 4)

lc_rc5 = lc_rc*100 + TCA.trd2.grd 
lc_rc5.df = data.frame(freq(lc_rc5))[-29,]
lc_rc5.df$lc = floor(lc_rc5.df$value/100)
lc_rc5.df$trend = lc_rc5.df$value - lc_rc5.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc5.df, FUN = "sum")
lc_rc5.df$prop = 100*lc_rc5.df$count/rep(lc_sum$count, each = 4)


lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                  data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))),
                  data.frame(lc_rc3.df, VI = rep("NDWI", nrow(lc_rc3.df))),
                 data.frame(lc_rc4.df, VI = rep("TCG", nrow(lc_rc4.df))),
                 data.frame(lc_rc5.df, VI = rep("TCA", nrow(lc_rc5.df))))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 1 | lc_rc.df$lc == 7| lc_rc.df$lc == 5), ] #Remove class 1: water; and 7: other class
lc_rc.df = lc_rc.df[-which(lc_rc.df$trend == 4), ] #remove Not Calculated classes

library(reshape2)
library(ggplot2)

lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative","Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- clasnames[as.numeric(levels(lc_rc1.df.long$lc))]

ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover") + 
  theme(legend.position="none")
  
ggsave(".\\NBAR_results3\\trend_biomes.png", width = 10, height = 7.5, units = "in")


#plot number of obs used to calculate trends
TCW2.dry.trd = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.grd")
Numb.recls = recls2(TCW2.dry.trd[[3]], threshold = c(4,6,8,10,12))

threshold = c(4,6,8,10,12)
color = "Spectral"

clscolor_change = brewer.pal(length(threshold)+2,color)[c(2:(length(threshold)+2))]
clscolor_change = c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837")

breaks2_change <- 0:(length(threshold)+1)        
legendbrks2_change <- 1:(length(threshold)+1) - 0.5

clasnames_change = c(paste("<=", threshold[1], sep = " "))
for(j in 1:(length(threshold)-1)){
  
  clasnames_change = c(clasnames_change, paste(threshold[j], "-", threshold[j+1], sep = " "))
  
}
clasnames_change[length(threshold)+1] <- c(paste(">", threshold[length(threshold)], sep = " "))

png(file = ".\\NBAR_results3\\trend.numb.wa.png", width = 2000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.3, lines=2)
levelplot(Numb.recls,
          maxpixels = nrow(Numb.recls)*ncol(Numb.recls),
          at= breaks2_change, margin=FALSE,
          col.regions= clscolor_change,
          colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()


