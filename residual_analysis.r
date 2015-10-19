#10/19/2015
#residual trend analsysi for all the indices in the whole west africa

#####################################################################################################
#####   based on one month lag: from de Wasseige et al (2003) and Zhou et al (2014)    ##############
#####################################################################################################

## Section 2: residual trend analysis
## 1. select the suitable temporal scale for rainfall summarization
# 1.1 de Wasseige et al (2003 Agricultural and Forest Meteorology) described that LAI lag one month after dry season starts
#     in moist deciduous Congo–Guinean forest. [cited by Zhou et al 2014 nature congo forest]
# 1.2 Zhou et al. (2014) find April-May-June EVI correlated strongly with April-May-June and March-Apr-May rainfall.
#     [no lag & one month before]
#1.3 Enquist, B. J. & Enquist, C. A. F. Long-term change within a neotropical forest: assessing differential 
#    functional and floristic responses to disturbance and drought. Glob. Change Biol. 17, 1408–1424 (2011) 
#   Enquist and Enquist note that the tropical forest had different sensistivity to rainfall change, therefore we 
#  examined up to 3 month before the dry season


setwd("D:\\users\\Zhihua\\MODIS")
library("MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)
library(Kendall)

#read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia"),]

county_b_ghana = extent(-13.5, 1.2, 4.35, 9.5)

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
TCW2.dry = stack(".\\NBAR_results3\\TCW2.dry.wa.grd")
TCG2.dry = stack(".\\NBAR_results3\\TCG2.dry.wa.grd")
TCB2.dry = stack(".\\NBAR_results3\\TCB2.dry.wa.grd")
TCA2.dry = stack(".\\NBAR_results3\\TCA2.dry.wa.grd")
EVI2.dry = stack(".\\NBAR_results3\\EVI2.dry.wa.grd")
NDWI2.dry = stack(".\\NBAR_results3\\NDWI2.dry.wa.grd")

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
tcg.df = raster::extract(TCG2.dry, pts.sp)
tcb.df = raster::extract(TCB2.dry, pts.sp)
tca.df = raster::extract(TCA2.dry, pts.sp)
evi.df = raster::extract(EVI2.dry, pts.sp)
ndwi.df = raster::extract(NDWI2.dry, pts.sp)

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
rf.tcg = data.frame(rf.df, tcg = apply(tcg.df, 2,median, na.rm = T)) 
colnames(rf.tcg) <- colnames
rf.tcb = data.frame(rf.df, tcb = apply(tcb.df, 2,median, na.rm = T)) 
colnames(rf.tcb) <- colnames
rf.tca = data.frame(rf.df, tca = apply(tca.df, 2,median, na.rm = T)) 
colnames(rf.tca) <- colnames
rf.evi = data.frame(rf.df, evi = apply(evi.df, 2,median, na.rm = T)) 
colnames(rf.evi) <- colnames
rf.ndwi = data.frame(rf.df, ndwi = apply(ndwi.df, 2,median, na.rm = T)) 
colnames(rf.ndwi) <- colnames

write.csv(rf.tcw,".\\NBAR_results3\\rainfall.TCW.relation.csv")
write.csv(rf.tcg,".\\NBAR_results3\\rainfall.TCG.relation.csv")
write.csv(rf.tcb,".\\NBAR_results3\\rainfall.TCB.relation.csv")
write.csv(rf.tca,".\\NBAR_results3\\rainfall.TCA.relation.csv")
write.csv(rf.evi,".\\NBAR_results3\\rainfall.EVI.relation.csv")
write.csv(rf.ndwi,".\\NBAR_results3\\rainfall.NWDI.relation.csv")

library(psych)
png(file = ".\\NBAR_results3\\rainfall.TCW.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.tcw)

dev.off()

png(file = ".\\NBAR_results3\\rainfall.TCG.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.tcg)

dev.off()

png(file = ".\\NBAR_results3\\rainfall.TCB.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.tcb)

dev.off()

png(file = ".\\NBAR_results3\\rainfall.TCA.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.tca)

dev.off()

png(file = ".\\NBAR_results3\\rainfall.EVI.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.evi)

dev.off()

png(file = ".\\NBAR_results3\\rainfall.NDWI.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.ndwi)

dev.off()

##############################################################
##residual analysis, use 1 month before season rainfall
tcw.1.df = data.frame(tcw.df, rf.1mon.df)
tcw.rf.cor = apply(tcw.1.df, 1, fun.cor2)                                #based on strongest relationship
tcw.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tcw.rf.cor[1,], FUN="*") #predict from rainfall data
tcw.df.pred = sweep(tcw.df.pred, MARGIN = 1, STATS = tcw.rf.cor[2,], FUN="+") #add the intercept
tcw.df.residual = tcw.df - tcw.df.pred                                        #calculate residual
tcw.df.residual.trd <- apply(tcw.df.residual, 1, mk1)                             #calculate residual trend

tcb.1.df = data.frame(tcb.df, rf.1mon.df)
tcb.rf.cor = apply(tcb.1.df, 1, fun.cor2)                                #based on strongest relationship
tcb.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tcb.rf.cor[1,], FUN="*") #predict from rainfall data
tcb.df.pred = sweep(tcb.df.pred, MARGIN = 1, STATS = tcb.rf.cor[2,], FUN="+") #add the intercept
tcb.df.residual = tcb.df - tcb.df.pred                                        #calculate residual
tcb.df.residual.trd <- apply(tcb.df.residual, 1, mk1)                             #calculate residual trend

tcg.1.df = data.frame(tcg.df, rf.1mon.df)
tcg.rf.cor = apply(tcg.1.df, 1, fun.cor2)                                #based on strongest relationship
tcg.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tcg.rf.cor[1,], FUN="*") #predict from rainfall data
tcg.df.pred = sweep(tcg.df.pred, MARGIN = 1, STATS = tcg.rf.cor[2,], FUN="+") #add the intercept
tcg.df.residual = tcg.df - tcg.df.pred                                        #calculate residual
tcg.df.residual.trd <- apply(tcg.df.residual, 1, mk1)                             #calculate residual trend

tca.1.df = data.frame(tca.df, rf.1mon.df)
tca.rf.cor = apply(tca.1.df, 1, fun.cor2)                                #based on strongest relationship
tca.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tca.rf.cor[1,], FUN="*") #predict from rainfall data
tca.df.pred = sweep(tca.df.pred, MARGIN = 1, STATS = tca.rf.cor[2,], FUN="+") #add the intercept
tca.df.residual = tca.df - tca.df.pred                                        #calculate residual
tca.df.residual.trd <- apply(tca.df.residual, 1, mk1)                             #calculate residual trend

evi.1.df = data.frame(evi.df, rf.1mon.df)
evi.rf.cor = apply(evi.1.df, 1, fun.cor2)                                #based on strongest relationship
evi.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = evi.rf.cor[1,], FUN="*") #predict from rainfall data
evi.df.pred = sweep(evi.df.pred, MARGIN = 1, STATS = evi.rf.cor[2,], FUN="+") #add the intercept
evi.df.residual = evi.df - evi.df.pred                                        #calculate residual
evi.df.residual.trd <- apply(evi.df.residual, 1, mk1)                             #calculate residual trend

ndwi.1.df = data.frame(ndwi.df, rf.1mon.df)
ndwi.rf.cor = apply(ndwi.1.df, 1, fun.cor2)                                #based on strongest relationship
ndwi.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = ndwi.rf.cor[1,], FUN="*") #predict from rainfall data
ndwi.df.pred = sweep(ndwi.df.pred, MARGIN = 1, STATS = ndwi.rf.cor[2,], FUN="+") #add the intercept
ndwi.df.residual = ndwi.df - ndwi.df.pred                                        #calculate residual
ndwi.df.residual.trd <- apply(ndwi.df.residual, 1, mk1)                             #calculate residual trend


# 3.4 convert to raster, and classify, save
#for tcw
tcw.df.residual.trd.value = as.numeric(tcw.df.residual.trd[1,])                                    
tcw.df.residual.trd.p = as.numeric(tcw.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tcw.rf.cor[1,], raster = TCW2.dry[[1]])
corr[lcc.grd != 1] = NA

tcw.df.residual.trd.value = Point2raster(tcw.df.residual.trd.value, raster = TCW2.dry[[1]])
tcw.df.residual.trd.p = Point2raster(tcw.df.residual.trd.p, raster = TCW2.dry[[1]])

tcw.df.residual.trd.value2 = tcw.df.residual.trd.value
tcw.df.residual.trd.value2[tcw.df.residual.trd.value < 0 & tcw.df.residual.trd.p <= 0.1] = 1
tcw.df.residual.trd.value2[tcw.df.residual.trd.value >= 0 & tcw.df.residual.trd.p <= 0.1] = 2
tcw.df.residual.trd.value2[tcw.df.residual.trd.p > 0.1] = 3
tcw.df.residual.trd.value2[lcc.grd == 1&is.na(tcw.df.residual.trd.p)] = 4
#save results
writeRaster(tcw.df.residual.trd.value2,".\\NBAR_results3\\tcw.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\tcw.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for tcb
tcb.df.residual.trd.value = as.numeric(tcb.df.residual.trd[1,])                                    
tcb.df.residual.trd.p = as.numeric(tcb.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tcb.rf.cor[1,], raster = tcb2.dry[[1]])
corr[lcc.grd != 1] = NA

tcb.df.residual.trd.value = Point2raster(tcb.df.residual.trd.value, raster = tcb2.dry[[1]])
tcb.df.residual.trd.p = Point2raster(tcb.df.residual.trd.p, raster = tcb2.dry[[1]])

tcb.df.residual.trd.value2 = tcb.df.residual.trd.value
tcb.df.residual.trd.value2[tcb.df.residual.trd.value < 0 & tcb.df.residual.trd.p <= 0.1] = 1
tcb.df.residual.trd.value2[tcb.df.residual.trd.value >= 0 & tcb.df.residual.trd.p <= 0.1] = 2
tcb.df.residual.trd.value2[tcb.df.residual.trd.p > 0.1] = 3
tcb.df.residual.trd.value2[lcc.grd == 1&is.na(tcb.df.residual.trd.p)] = 4
#save results
writeRaster(tcb.df.residual.trd.value2,".\\NBAR_results3\\tcb.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\tcb.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for tcg
tcg.df.residual.trd.value = as.numeric(tcg.df.residual.trd[1,])                                    
tcg.df.residual.trd.p = as.numeric(tcg.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tcg.rf.cor[1,], raster = tcg2.dry[[1]])
corr[lcc.grd != 1] = NA

tcg.df.residual.trd.value = Point2raster(tcg.df.residual.trd.value, raster = tcg2.dry[[1]])
tcg.df.residual.trd.p = Point2raster(tcg.df.residual.trd.p, raster = tcg2.dry[[1]])

tcg.df.residual.trd.value2 = tcg.df.residual.trd.value
tcg.df.residual.trd.value2[tcg.df.residual.trd.value < 0 & tcg.df.residual.trd.p <= 0.1] = 1
tcg.df.residual.trd.value2[tcg.df.residual.trd.value >= 0 & tcg.df.residual.trd.p <= 0.1] = 2
tcg.df.residual.trd.value2[tcg.df.residual.trd.p > 0.1] = 3
tcg.df.residual.trd.value2[lcc.grd == 1&is.na(tcg.df.residual.trd.p)] = 4
#save results
writeRaster(tcg.df.residual.trd.value2,".\\NBAR_results3\\tcg.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\tcg.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for tca
tca.df.residual.trd.value = as.numeric(tca.df.residual.trd[1,])                                    
tca.df.residual.trd.p = as.numeric(tca.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tca.rf.cor[1,], raster = tca2.dry[[1]])
corr[lcc.grd != 1] = NA

tca.df.residual.trd.value = Point2raster(tca.df.residual.trd.value, raster = tca2.dry[[1]])
tca.df.residual.trd.p = Point2raster(tca.df.residual.trd.p, raster = tca2.dry[[1]])

tca.df.residual.trd.value2 = tca.df.residual.trd.value
tca.df.residual.trd.value2[tca.df.residual.trd.value < 0 & tca.df.residual.trd.p <= 0.1] = 1
tca.df.residual.trd.value2[tca.df.residual.trd.value >= 0 & tca.df.residual.trd.p <= 0.1] = 2
tca.df.residual.trd.value2[tca.df.residual.trd.p > 0.1] = 3
tca.df.residual.trd.value2[lcc.grd == 1&is.na(tca.df.residual.trd.p)] = 4
#save results
writeRaster(tca.df.residual.trd.value2,".\\NBAR_results3\\tca.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\tca.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for evi
evi.df.residual.trd.value = as.numeric(evi.df.residual.trd[1,])                                    
evi.df.residual.trd.p = as.numeric(evi.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(evi.rf.cor[1,], raster = evi2.dry[[1]])
corr[lcc.grd != 1] = NA

evi.df.residual.trd.value = Point2raster(evi.df.residual.trd.value, raster = evi2.dry[[1]])
evi.df.residual.trd.p = Point2raster(evi.df.residual.trd.p, raster = evi2.dry[[1]])

evi.df.residual.trd.value2 = evi.df.residual.trd.value
evi.df.residual.trd.value2[evi.df.residual.trd.value < 0 & evi.df.residual.trd.p <= 0.1] = 1
evi.df.residual.trd.value2[evi.df.residual.trd.value >= 0 & evi.df.residual.trd.p <= 0.1] = 2
evi.df.residual.trd.value2[evi.df.residual.trd.p > 0.1] = 3
evi.df.residual.trd.value2[lcc.grd == 1&is.na(evi.df.residual.trd.p)] = 4
#save results
writeRaster(evi.df.residual.trd.value2,".\\NBAR_results3\\evi.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\evi.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

#for ndwi
ndwi.df.residual.trd.value = as.numeric(ndwi.df.residual.trd[1,])                                    
ndwi.df.residual.trd.p = as.numeric(ndwi.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(ndwi.rf.cor[1,], raster = ndwi2.dry[[1]])
corr[lcc.grd != 1] = NA

ndwi.df.residual.trd.value = Point2raster(ndwi.df.residual.trd.value, raster = ndwi2.dry[[1]])
ndwi.df.residual.trd.p = Point2raster(ndwi.df.residual.trd.p, raster = ndwi2.dry[[1]])

ndwi.df.residual.trd.value2 = ndwi.df.residual.trd.value
ndwi.df.residual.trd.value2[ndwi.df.residual.trd.value < 0 & ndwi.df.residual.trd.p <= 0.1] = 1
ndwi.df.residual.trd.value2[ndwi.df.residual.trd.value >= 0 & ndwi.df.residual.trd.p <= 0.1] = 2
ndwi.df.residual.trd.value2[ndwi.df.residual.trd.p > 0.1] = 3
ndwi.df.residual.trd.value2[lcc.grd == 1&is.na(ndwi.df.residual.trd.p)] = 4
#save results
writeRaster(ndwi.df.residual.trd.value2,".\\NBAR_results3\\ndwi.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results3\\ndwi.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

##compare residual trend and raw trend
trd.grd2 = raster(".\\NBAR_results2\\trend.wa.2.tif")
trd.comparsion = trd.grd2 == tcw.df.residual.trd.value2
freq(trd.comparsion)
