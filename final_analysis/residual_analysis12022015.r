# go to github https://github.com/liuzh811/
# ForestDegradationWestAfrica/blob/master/residual_analysis_all_indices.r
# change the 
county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#1. read into TRMM data
fn = list.files(path = "D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2", pattern = "*.nc$")
r = stack(paste("D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2\\", fn, sep = ""))
#r = crop(r, county_b.geo) #crop the dataset to western Africa boundary
r = crop(r, county_b_ghana) #crop the dataset to western Africa boundary
names(r) <- fn

#2. calulate rainfall 1, 2,3,4,5,6 month before dry season, and during the dry season
rf.dry = list()  # no lag
rf.1mon = list() #1 month before
rf.2mon = list() #2 month before
rf.3mon = list() #3 month before
rf.4mon = list()
rf.5mon = list()
rf.6mon = list()
rf.7mon = list()
rf.8mon = list() #rainy season rainfall
rf.9mon = list()
rf.10mon = list()


for(i in 2000:2014){
  idx.1 = which(as.numeric(substr(names(r), 13,16))==i & as.numeric(substr(names(r), 18,19))==11 & as.numeric(substr(names(r), 21,22))==15)
  idx.1mon = (idx.1-30):(idx.1+60)
  idx.2mon = (idx.1-60):(idx.1+30)
  idx.3mon = (idx.1-90):idx.1
  idx.8mon = (idx.1-240):idx.1
  
  rf.1mon[[i-1999]] <- calc(r[[idx.1mon]], sum) 
  rf.2mon[[i-1999]] <- calc(r[[idx.2mon]], sum) 
  rf.3mon[[i-1999]] <- calc(r[[idx.3mon]], sum) 
  rf.8mon[[i-1999]] <- calc(r[[idx.8mon]], sum) 
  
  idx.4mon = (idx.1-120):idx.1
  idx.5mon = (idx.1-150):idx.1
  idx.6mon = (idx.1-180):idx.1
  idx.7mon = (idx.1-210):idx.1
  idx.9mon = (idx.1-270):idx.1
  idx.10mon = (idx.1-300):idx.1
  
  if(i < 2014){idx.dry = idx.1:(idx.1+135)} else (idx.dry = idx.1:(idx.1+70))
  rf.4mon[[i-1999]] <- calc(r[[idx.4mon]], sum)
  rf.5mon[[i-1999]] <- calc(r[[idx.5mon]], sum)
  rf.6mon[[i-1999]] <- calc(r[[idx.6mon]], sum)
  rf.7mon[[i-1999]] <- calc(r[[idx.7mon]], sum)
  rf.9mon[[i-1999]] <- calc(r[[idx.9mon]], sum)
  rf.10mon[[i-1999]] <- calc(r[[idx.10mon]], sum)
  rf.dry[[i-1999]] <- calc(r[[idx.dry]], sum)
  
}

rf.1mon = stack(rf.1mon)
rf.2mon = stack(rf.2mon)
rf.3mon = stack(rf.3mon)
rf.8mon = stack(rf.8mon)

rf.4mon = stack(rf.4mon)
rf.5mon = stack(rf.5mon)
rf.6mon = stack(rf.6mon)
rf.7mon = stack(rf.7mon)
rf.9mon = stack(rf.9mon)
rf.10mon = stack(rf.10mon)
rf.dry = stack(rf.dry)

# 3. residual analysis
# 3. find the strongest relationship between rainfall and VI
# 3.1. get points
TCW2.dry = stack(".\\NBAR_results4\\TCW2.dry.wa.grd")
EVI2.dry = stack(".\\NBAR_results4\\EVI2.dry.wa.grd")

pts.sp = Ex.pts.all(TCW2.dry[[1]])

#extract rainfall
rf.1mon.df = raster::extract(rf.1mon, pts.sp)
rf.2mon.df = raster::extract(rf.2mon, pts.sp)
rf.3mon.df = raster::extract(rf.3mon, pts.sp)
rf.8mon.df = raster::extract(rf.8mon, pts.sp)

rf.4mon.df = raster::extract(rf.4mon, pts.sp)
rf.5mon.df = raster::extract(rf.5mon, pts.sp)
rf.6mon.df = raster::extract(rf.6mon, pts.sp)
rf.7mon.df = raster::extract(rf.7mon, pts.sp)
rf.9mon.df = raster::extract(rf.9mon, pts.sp)
rf.10mon.df = raster::extract(rf.10mon, pts.sp)
rf.dry.df = raster::extract(rf.dry, pts.sp)

#extract vi
tcw.df = raster::extract(TCW2.dry, pts.sp)
evi.df = raster::extract(EVI2.dry, pts.sp)

#combine rainfall and vi
rf.df = data.frame(rf.dry = apply(rf.dry.df, 2,median, na.rm = T),
                   rf.1mon = apply(rf.1mon.df, 2,median, na.rm = T),
                   rf.2mon = apply(rf.2mon.df, 2,median, na.rm = T),
                   rf.3mon = apply(rf.3mon.df, 2,median, na.rm = T),
                   rf.4mon = apply(rf.4mon.df, 2,median, na.rm = T),
                   rf.5mon = apply(rf.5mon.df, 2,median, na.rm = T),
                   rf.6mon = apply(rf.6mon.df, 2,median, na.rm = T),
                   rf.7mon = apply(rf.7mon.df, 2,median, na.rm = T),
                   rf.8mon = apply(rf.8mon.df, 2,median, na.rm = T),
                   rf.9mon = apply(rf.9mon.df, 2,median, na.rm = T),
                   rf.10mon = apply(rf.10mon.df, 2,median, na.rm = T))


colnames <- c("No.lag","1.Month.lag","2.Month.lag","3.Month.lag",
              "4.Month.lag","5.Month.lag","6.Month.lag",
              "7.Month.lag","8.Month.lag","9.Month.lag","10.Month.lag","VI")

rf.tcw = data.frame(rf.df, tcw = apply(tcw.df, 2,median, na.rm = T)) 
colnames(rf.tcw) <- colnames

rf.evi = data.frame(rf.df, evi = apply(evi.df, 2,median, na.rm = T)) 
colnames(rf.evi) <- colnames


write.csv(rf.tcw,".\\NBAR_results4\\rainfall.TCW.relation.csv")
write.csv(rf.evi,".\\NBAR_results4\\rainfall.EVI.relation.csv")

library(psych)
png(file = ".\\NBAR_results4\\rainfall.TCW.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.tcw)

dev.off()


png(file = ".\\NBAR_results4\\rainfall.EVI.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.evi)

dev.off()


##############################################################
##residual analysis, use 1 month before season rainfall
tcw.1.df = data.frame(tcw.df, rf.1mon.df)
tcw.rf.cor = apply(tcw.1.df, 1, fun.cor2)                                #based on strongest relationship
tcw.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tcw.rf.cor[1,], FUN="*") #predict from rainfall data
tcw.df.pred = sweep(tcw.df.pred, MARGIN = 1, STATS = tcw.rf.cor[2,], FUN="+") #add the intercept
tcw.df.residual = tcw.df - tcw.df.pred                                        #calculate residual
tcw.df.residual.trd <- apply(tcw.df.residual, 1, mk1)                             #calculate residual trend
#calculate residual trend

evi.1.df = data.frame(evi.df, rf.1mon.df)
evi.rf.cor = apply(evi.1.df, 1, fun.cor2)                                #based on strongest relationship
evi.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = evi.rf.cor[1,], FUN="*") #predict from rainfall data
evi.df.pred = sweep(evi.df.pred, MARGIN = 1, STATS = evi.rf.cor[2,], FUN="+") #add the intercept
evi.df.residual = evi.df - evi.df.pred                                        #calculate residual
evi.df.residual.trd <- apply(evi.df.residual, 1, mk1)                             #calculate residual trend


# 3.4 convert to raster, and classify, save
#for tcw
tcw.df.residual.trd.value = as.numeric(tcw.df.residual.trd[1,])                                    
tcw.df.residual.trd.p = as.numeric(tcw.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tcw.rf.cor[1,], raster = TCW2.dry[[1]])

tcw.df.residual.trd.value = Point2raster(tcw.df.residual.trd.value, raster = TCW2.dry[[1]])
tcw.df.residual.trd.p = Point2raster(tcw.df.residual.trd.p, raster = TCW2.dry[[1]])

tcw.df.residual.trd.value2 = tcw.df.residual.trd.value
tcw.df.residual.trd.value2[tcw.df.residual.trd.value < 0 & tcw.df.residual.trd.p <= 0.1] = 1
tcw.df.residual.trd.value2[tcw.df.residual.trd.value >= 0 & tcw.df.residual.trd.p <= 0.1] = 2
tcw.df.residual.trd.value2[tcw.df.residual.trd.p > 0.1] = 3
tcw.df.residual.trd.value2[is.na(tcw.df.residual.trd.p)] = 4
#save results
writeRaster(tcw.df.residual.trd.value2,".\\NBAR_results4\\tcw.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results4\\tcw.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for evi
evi.df.residual.trd.value = as.numeric(evi.df.residual.trd[1,])                                    
evi.df.residual.trd.p = as.numeric(evi.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(evi.rf.cor[1,], raster = TCW2.dry[[1]])

evi.df.residual.trd.value = Point2raster(evi.df.residual.trd.value, raster = TCW2.dry[[1]])
evi.df.residual.trd.p = Point2raster(evi.df.residual.trd.p, raster = TCW2.dry[[1]])

evi.df.residual.trd.value2 = evi.df.residual.trd.value
evi.df.residual.trd.value2[evi.df.residual.trd.value < 0 & evi.df.residual.trd.p <= 0.1] = 1
evi.df.residual.trd.value2[evi.df.residual.trd.value >= 0 & evi.df.residual.trd.p <= 0.1] = 2
evi.df.residual.trd.value2[evi.df.residual.trd.p > 0.1] = 3
evi.df.residual.trd.value2[is.na(evi.df.residual.trd.p)] = 4
#save results
writeRaster(evi.df.residual.trd.value2,".\\NBAR_results4\\evi.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results4\\evi.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#################################
## plot residual
TCW.trd2.grd.res = raster(".\\NBAR_results4\\tcw.df.residual.trd.wa2.tif")
EVI.trd2.grd.res = raster(".\\NBAR_results4\\evi.df.residual.trd.wa2.tif")

county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")

TCW.trd2.grd.res = TCW.trd2.grd.res*county_b_ghana.r
EVI.trd2.grd.res = EVI.trd2.grd.res*county_b_ghana.r

trd.res = stack(TCW.trd2.grd.res, EVI.trd2.grd.res)
names(trd.res) <- c("TCW","EVI")

breaks2_change <- 0:4        
legendbrks2_change <- 1:4 - 0.5

arg1 <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels = c("Negative","Positive","No Trend","Not Calculated")
color1=c("#e66101", "#1a9641","#ffffbf","#2b83ba")


png(file = ".\\NBAR_results4\\residual.trend.wa.png", width = 4000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(trd.res,
          maxpixels = nrow(trd.res)*ncol(trd.res),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5)) 

dev.off()

#calulate the correlation between strongly rainfall and VIs
# based on correlations, the 2 month lag had the strongest correlations
tcw.2.df = data.frame(tcw.df, rf.2mon.df)
evi.2.df = data.frame(evi.df, rf.2mon.df)

tcw.rf.cor.test = apply(tcw.2.df, 1, COR.test)    
tcw.rf.corr = Point2raster(tcw.rf.cor.test[1,], county_b_ghana.r)
tcw.rf.corr.sig = Point2raster(tcw.rf.cor.test[2,], county_b_ghana.r)
tcw.rf.corr = tcw.rf.corr*county_b_ghana.r
tcw.rf.corr.sig = tcw.rf.corr.sig*county_b_ghana.r

evi.rf.cor.test = apply(evi.2.df, 1, COR.test)    
evi.rf.corr = Point2raster(evi.rf.cor.test[1,], county_b_ghana.r)
evi.rf.corr.sig = Point2raster(evi.rf.cor.test[2,], county_b_ghana.r)
evi.rf.corr = evi.rf.corr*county_b_ghana.r
evi.rf.corr.sig = evi.rf.corr.sig*county_b_ghana.r

threshold =  c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
tcw.rf.corr2 = recls2(tcw.rf.corr, threshold = threshold)
evi.rf.corr2 = recls2(evi.rf.corr, threshold = threshold)

s = stack(tcw.rf.corr2, evi.rf.corr2)
names(s) <- c("TCW", "EVI")

clscolor_change = brewer.pal(9,"Greens")[c(2:9)]
breaks2_change <- 0:8        
legendbrks2_change <- 1:8 - 0.5

clasnames_change <- c("< 0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6",">0.6")

png("D:/users/Zhihua/MODIS/NBAR_results4/rainfall_vi_relationship.png",height = 3000, width = 5000, res = 300, units = "px")

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(s,
          maxpixels = nrow(s)*ncol(s),
          at= breaks2_change, margin=FALSE,
          col.regions= clscolor_change,
          colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) 

dev.off()

