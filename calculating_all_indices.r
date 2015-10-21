#10/19/2015
#calculate trends for all the indices

setwd("D:\\users\\Zhihua\\MODIS")
library("MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)
library(Kendall)

#3. read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia"),]

county_b_ghana = extent(-13.5, 1.2, 4.35, 9.5)

###################################################################################
## Section 1: calculate the trend of vegetation trend
#1 read bands info
band1 <- preStack(path=".\\NBAR", pattern="*Band1.tif$")
band2 <- preStack(path=".\\NBAR", pattern="*Band2.tif$")
band3 <- preStack(path=".\\NBAR", pattern="*Band3.tif$")
band4 <- preStack(path=".\\NBAR", pattern="*Band4.tif$")
band5 <- preStack(path=".\\NBAR", pattern="*Band5.tif$")
band6 <- preStack(path=".\\NBAR", pattern="*Band6.tif$")
band7 <- preStack(path=".\\NBAR", pattern="*Band7.tif$")
YearDOY.band = substr(band1, 17, 23)

YearDOYModis = data.frame(Year = as.numeric(substr(YearDOY.band,1,4)), DOY = as.numeric(substr(YearDOY.band,5,7)))

#' find the index for dry season [Nov 15 – Apr 1, ]
EarlyDry.idx = c(which(YearDOYModis$Year == 2000 & YearDOYModis$DOY>318), which(YearDOYModis$Year >= 2001 & (YearDOYModis$DOY>318 |YearDOYModis$DOY<90)))
EarlyDry.idx.M = YearDOYModis[EarlyDry.idx,]
EarlyDry.idx.M$Year[which(EarlyDry.idx.M$DOY>318)] = EarlyDry.idx.M$Year[which(EarlyDry.idx.M$DOY>318)]+1 #reassign hydro-year

#read only dry season data
band1 = stack(band1[EarlyDry.idx])
band2 = stack(band2[EarlyDry.idx])
band3 = stack(band3[EarlyDry.idx])
band4 = stack(band4[EarlyDry.idx])
band5 = stack(band5[EarlyDry.idx])
band6 = stack(band6[EarlyDry.idx])
band7 = stack(band7[EarlyDry.idx])

# band albedo quality
band.qc = preStack(path=".\\NBAR_QC", pattern="*Albedo_Band_Quality.tif$")
YearDOY.band.qc = substr(band.qc , 20, 26)
band.qc = stack(band.qc[EarlyDry.idx])

#crop the dataset to study area
band11 = list()
band21 = list()
band31 = list()
band41 = list()
band51 = list()
band61 = list()
band71 = list()
band.qc1 = list()

for(i in 1:nlayers(band1)){
  
  band1.1 = crop(band1[[i]], county_b_ghana)
  band2.1 = crop(band2[[i]], county_b_ghana)
  band3.1 = crop(band3[[i]], county_b_ghana)
  band4.1 = crop(band4[[i]], county_b_ghana)
  band5.1 = crop(band5[[i]], county_b_ghana)
  band6.1 = crop(band6[[i]], county_b_ghana)
  band7.1 = crop(band7[[i]], county_b_ghana)
  band.qc.1 = crop(band.qc[[i]], county_b_ghana)
  
  band11[[i]] = band1.1
  band21[[i]] = band2.1
  band31[[i]] = band3.1
  band41[[i]] = band4.1
  band51[[i]] = band5.1
  band61[[i]] = band6.1
  band71[[i]] = band7.1
  band.qc1[[i]] = band.qc.1
  
print(paste("Finish Cropping data for ", i, "of", nlayers(band1), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

band11 = stack(band11)
band21 = stack(band21)
band31 = stack(band31)
band41 = stack(band41)
band51 = stack(band51)
band61 = stack(band61)
band71 = stack(band71)
band.qc1 = stack(band.qc1)

rm(list = c("band1","band2","band3","band4","band5","band6","band7","band.qc"))

# calculate wetness indices, using different level of QC values

TCW1 = list() # select 0,1
TCW2 = list() # select 0,1,2
TCW3 = list() # select 0,1,2,3

TCG1 = list() #TC greenness 
TCG2 = list() #TC greenness
TCG3 = list() #TC greenness

TCB1 = list() #TC brightness
TCB2 = list() #TC brightness
TCB3 = list() #TC brightness

TCA1 = list() #TC angle 
TCA2 = list() #TC angle 
TCA3 = list() #TC angle 

NDWI1 = list() #normalized difference moisture index
NDWI2 = list() #normalized difference moisture index
NDWI3 = list() #normalized difference moisture index

EVI1 = list() #normalized difference moisture index
EVI2 = list() #normalized difference moisture index
EVI3 = list() #normalized difference moisture index

for (i in 1:nlayers(band11) ){
  
  #calculate indices
  brightness = 0.4395*band11[[i]]+0.5945*band21[[i]]+0.246*band31[[i]]+0.3918*band41[[i]]+0.3506*band51[[i]]+0.2136*band61[[i]]+0.2678*band71[[i]]
  greenness = -0.4064*band11[[i]]+0.5129*band21[[i]]-0.2744*band31[[i]]-0.2893*band41[[i]]+0.4882*band51[[i]]-0.0036*band61[[i]]-0.4169*band71[[i]]
  wetness = 0.1147*band11[[i]]+0.2489*band21[[i]]+0.2408*band31[[i]]+0.3132*band41[[i]]-0.3122*band51[[i]]-0.6416*band61[[i]]-0.5087*band71[[i]]
  EVI.1 = 2.5*(band21[[i]]*0.0001-band11[[i]]*0.0001)/(1+band21[[i]]*0.0001-band11[[i]]*0.0001)
  TCA.1 = atan(greenness/brightness)
  NDWI.1 = (band21[[i]]-band61[[i]])/(band21[[i]]+band61[[i]])
  
  qc.r = QAfind.mt2(band.qc1[[i]])
  
  #select 0,1
  qc.r1 = qc.r[[2]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  
  brightness1 = brightness*qc.r1
  greenness1 = greenness*qc.r1
  wetness1 = wetness*qc.r1
  EVI.11 = EVI.1*qc.r1
  TCA11 = TCA.1*qc.r1
  NDWI11 = NDWI.1*qc.r1
  
  TCB1[[i]] = brightness1
  TCG1[[i]] = greenness1
  TCW1[[i]] = wetness1
  TCA1[[i]] = TCA11
  EVI1[[i]] = EVI.11
  NDWI1[[i]] = NDWI11
  
  #select 0,1,2
  qc.r1 = qc.r[[3]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  
  brightness1 = brightness*qc.r1
  greenness1 = greenness*qc.r1
  wetness1 = wetness*qc.r1
  EVI.11 = EVI.1*qc.r1
  TCA11 = TCA.1*qc.r1
  NDWI11 = NDWI.1*qc.r1
  
  TCB2[[i]] = brightness1
  TCG2[[i]] = greenness1
  TCW2[[i]] = wetness1
  TCA2[[i]] = TCA11
  EVI2[[i]] = EVI.11
  NDWI2[[i]] = NDWI11

  #select 0,1,2,3
  qc.r1 = qc.r[[4]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  
  brightness1 = brightness*qc.r1
  greenness1 = greenness*qc.r1
  wetness1 = wetness*qc.r1
  EVI.11 = EVI.1*qc.r1
  TCA11 = TCA.1*qc.r1
  NDWI11 = NDWI.1*qc.r1
  
  TCB3[[i]] = brightness1
  TCG3[[i]] = greenness1
  TCW3[[i]] = wetness1
  TCA3[[i]] = TCA11
  EVI3[[i]] = EVI.11
  NDWI3[[i]] = NDWI11
  
  print(paste("Finish calculating TC for ", i, "of", nlayers(band11), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

rm(list = c("brightness","greenness","wetness","EVI.1","TCA.1","NDWI.1","qc.r","qc.r1"))
rm(list = c("brightness1","greenness1","wetness1","EVI.11","TCA11","NDWI11"))

TCW1 = stack(TCW1)
TCW2 = stack(TCW2)
TCW3 = stack(TCW3)

TCG1 = stack(TCG1)
TCG2 = stack(TCG2)
TCG3 = stack(TCG3)

TCB1 = stack(TCB1)
TCB2 = stack(TCB2)
TCB3 = stack(TCB3)

TCA1 = stack(TCA1)
TCA2 = stack(TCA2)
TCA3 = stack(TCA3)

EVI1 = stack(EVI1)
EVI2 = stack(EVI2)
EVI3 = stack(EVI3)

NDWI1 = stack(NDWI1)
NDWI2 = stack(NDWI2)
NDWI3 = stack(NDWI3)

# calculate annual median TCW
TCW1.dry = list() 
TCW2.dry = list()
TCW3.dry = list()

TCG1.dry = list() #TC greenness 
TCG2.dry = list() #TC greenness
TCG3.dry = list() #TC greenness

TCB1.dry = list() #TC brightness
TCB2.dry = list() #TC brightness
TCB3.dry = list() #TC brightness

TCA1.dry = list() #TC angle 
TCA2.dry = list() #TC angle 
TCA3.dry = list() #TC angle 

NDWI1.dry = list() #normalized difference moisture index
NDWI2.dry = list() #normalized difference moisture index
NDWI3.dry = list() #normalized difference moisture index

EVI1.dry = list() #normalized difference moisture index
EVI2.dry = list() #normalized difference moisture index
EVI3.dry = list() #normalized difference moisture index


for(i in 2001:2015){
  
  idx = which(EarlyDry.idx.M==i) #find the index for each year
  TCW1.dry1 <- subset(TCW1, idx, drop=FALSE)
  TCW1.dry[[i-2000]] <- calc(TCW1.dry1, median, na.rm = TRUE)
  TCB1.dry1 <- subset(TCB1, idx, drop=FALSE)
  TCB1.dry[[i-2000]] <- calc(TCB1.dry1, median, na.rm = TRUE)
  TCG1.dry1 <- subset(TCG1, idx, drop=FALSE)
  TCG1.dry[[i-2000]] <- calc(TCG1.dry1, median, na.rm = TRUE)
  TCA1.dry1 <- subset(TCA1, idx, drop=FALSE)
  TCA1.dry[[i-2000]] <- calc(TCA1.dry1, median, na.rm = TRUE)
  EVI1.dry1 <- subset(EVI1, idx, drop=FALSE)
  EVI1.dry[[i-2000]] <- calc(EVI1.dry1, median, na.rm = TRUE)
  NDWI1.dry1 <- subset(NDWI1, idx, drop=FALSE)
  NDWI1.dry[[i-2000]] <- calc(NDWI1.dry1, median, na.rm = TRUE)
  
  TCW2.dry1 <- subset(TCW2, idx, drop=FALSE)
  TCW2.dry[[i-2000]] <- calc(TCW2.dry1, median, na.rm = TRUE)
  TCB2.dry1 <- subset(TCB2, idx, drop=FALSE)
  TCB2.dry[[i-2000]] <- calc(TCB2.dry1, median, na.rm = TRUE)
  TCG2.dry1 <- subset(TCG2, idx, drop=FALSE)
  TCG2.dry[[i-2000]] <- calc(TCG2.dry1, median, na.rm = TRUE)
  TCA2.dry1 <- subset(TCA2, idx, drop=FALSE)
  TCA2.dry[[i-2000]] <- calc(TCA2.dry1, median, na.rm = TRUE)
  EVI2.dry1 <- subset(EVI2, idx, drop=FALSE)
  EVI2.dry[[i-2000]] <- calc(EVI2.dry1, median, na.rm = TRUE)
  NDWI2.dry1 <- subset(NDWI2, idx, drop=FALSE)
  NDWI2.dry[[i-2000]] <- calc(NDWI2.dry1, median, na.rm = TRUE)
  
  TCW3.dry1 <- subset(TCW3, idx, drop=FALSE)
  TCW3.dry[[i-2000]] <- calc(TCW3.dry1, median, na.rm = TRUE)
  TCB3.dry1 <- subset(TCB3, idx, drop=FALSE)
  TCB3.dry[[i-2000]] <- calc(TCB3.dry1, median, na.rm = TRUE)
  TCG3.dry1 <- subset(TCG3, idx, drop=FALSE)
  TCG3.dry[[i-2000]] <- calc(TCG3.dry1, median, na.rm = TRUE)
  TCA3.dry1 <- subset(TCA3, idx, drop=FALSE)
  TCA3.dry[[i-2000]] <- calc(TCA3.dry1, median, na.rm = TRUE)
  EVI3.dry1 <- subset(EVI3, idx, drop=FALSE)
  EVI3.dry[[i-2000]] <- calc(EVI3.dry1, median, na.rm = TRUE)
  NDWI3.dry1 <- subset(NDWI3, idx, drop=FALSE)
  NDWI3.dry[[i-2000]] <- calc(NDWI3.dry1, median, na.rm = TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

TCW1.dry = stack(TCW1.dry)
TCW2.dry = stack(TCW2.dry)
TCW3.dry = stack(TCW3.dry)

TCB1.dry = stack(TCB1.dry)
TCB2.dry = stack(TCB2.dry)
TCB3.dry = stack(TCB3.dry)

TCG1.dry = stack(TCG1.dry)
TCG2.dry = stack(TCG2.dry)
TCG3.dry = stack(TCG3.dry)

TCA1.dry = stack(TCA1.dry)
TCA2.dry = stack(TCA2.dry)
TCA3.dry = stack(TCA3.dry)

EVI1.dry = stack(EVI1.dry)
EVI2.dry = stack(EVI2.dry)
EVI3.dry = stack(EVI3.dry)

NDWI1.dry = stack(NDWI1.dry)
NDWI2.dry = stack(NDWI2.dry)
NDWI3.dry = stack(NDWI3.dry)


#write.raster
writeRaster(TCW1.dry,".\\NBAR_results3\\TCW1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW2.dry,".\\NBAR_results3\\TCW2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW3.dry,".\\NBAR_results3\\TCW3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCB1.dry,".\\NBAR_results3\\TCB1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCB2.dry,".\\NBAR_results3\\TCB2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCB3.dry,".\\NBAR_results3\\TCB3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCG1.dry,".\\NBAR_results3\\TCG1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG2.dry,".\\NBAR_results3\\TCG2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG3.dry,".\\NBAR_results3\\TCG3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCA1.dry,".\\NBAR_results3\\TCA1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA2.dry,".\\NBAR_results3\\TCA2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA3.dry,".\\NBAR_results3\\TCA3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(EVI1.dry,".\\NBAR_results3\\EVI1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI2.dry,".\\NBAR_results3\\EVI2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI3.dry,".\\NBAR_results3\\EVI3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(NDWI1.dry,".\\NBAR_results3\\NDWI1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI2.dry,".\\NBAR_results3\\NDWI2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI3.dry,".\\NBAR_results3\\NDWI3.dry.wa.grd", overwrite=TRUE) #write a raster stack files


#calculate trends
# test.pt = extract(TCW1.dry, 1) #test function
TCW1.dry.trd <- calc(TCW1.dry, mk1)
TCW2.dry.trd <- calc(TCW2.dry, mk1)
TCW3.dry.trd <- calc(TCW3.dry, mk1)

TCB1.dry.trd <- calc(TCB1.dry, mk1)
TCB2.dry.trd <- calc(TCB2.dry, mk1)
TCB3.dry.trd <- calc(TCB3.dry, mk1)

TCG1.dry.trd <- calc(TCG1.dry, mk1)
TCG2.dry.trd <- calc(TCG2.dry, mk1)
TCG3.dry.trd <- calc(TCG3.dry, mk1)

TCA1.dry.trd <- calc(TCA1.dry, mk1)
TCA2.dry.trd <- calc(TCA2.dry, mk1)
TCA3.dry.trd <- calc(TCA3.dry, mk1)

EVI1.dry.trd <- calc(EVI1.dry, mk1)
EVI2.dry.trd <- calc(EVI2.dry, mk1)
EVI3.dry.trd <- calc(EVI3.dry, mk1)

NDWI1.dry.trd <- calc(NDWI1.dry, mk1)
NDWI2.dry.trd <- calc(NDWI2.dry, mk1)
NDWI3.dry.trd <- calc(NDWI3.dry, mk1)

#save the results
writeRaster(TCW1.dry.trd,".\\NBAR_results3\\TCW1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW2.dry.trd,".\\NBAR_results3\\TCW2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW3.dry.trd,".\\NBAR_results3\\TCW3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCB1.dry.trd,".\\NBAR_results3\\TCB1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCB2.dry.trd,".\\NBAR_results3\\TCB2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCB3.dry.trd,".\\NBAR_results3\\TCB3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCG1.dry.trd,".\\NBAR_results3\\TCG1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG2.dry.trd,".\\NBAR_results3\\TCG2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG3.dry.trd,".\\NBAR_results3\\TCG3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(TCA1.dry.trd,".\\NBAR_results3\\TCA1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA2.dry.trd,".\\NBAR_results3\\TCA2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA3.dry.trd,".\\NBAR_results3\\TCA3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(EVI1.dry.trd,".\\NBAR_results3\\EVI1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI2.dry.trd,".\\NBAR_results3\\EVI2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI3.dry.trd,".\\NBAR_results3\\EVI3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files

writeRaster(NDWI1.dry.trd,".\\NBAR_results3\\NDWI1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI2.dry.trd,".\\NBAR_results3\\NDWI2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI3.dry.trd,".\\NBAR_results3\\NDWI3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files


#classified into different 
#TCW
names(TCW1.dry.trd) <- c("tau", "p", "obs")
names(TCW2.dry.trd) <- c("tau", "p", "obs")
names(TCW3.dry.trd) <- c("tau", "p", "obs")

TCW.trd1.grd = TCW1.dry.trd$tau
TCW.trd1.grd[TCW1.dry.trd$tau < 0 & TCW1.dry.trd$p <= 0.1] = 1
TCW.trd1.grd[TCW1.dry.trd$tau >= 0 & TCW1.dry.trd$p <= 0.1] = 2
TCW.trd1.grd[TCW1.dry.trd$p > 0.1] = 3
TCW.trd1.grd[is.na(TCW1.dry.trd$p)] = 4

TCW.trd2.grd = TCW1.dry.trd$tau
TCW.trd2.grd[TCW2.dry.trd$tau < 0 & TCW2.dry.trd$p <= 0.1] = 1
TCW.trd2.grd[TCW2.dry.trd$tau >= 0 & TCW2.dry.trd$p <= 0.1] = 2
TCW.trd2.grd[TCW2.dry.trd$p > 0.1] = 3
TCW.trd2.grd[is.na(TCW2.dry.trd$p)] = 4

TCW.trd3.grd = TCW3.dry.trd$tau
TCW.trd3.grd[TCW3.dry.trd$tau < 0 & TCW3.dry.trd$p <= 0.1] = 1
TCW.trd3.grd[TCW3.dry.trd$tau >= 0 & TCW3.dry.trd$p <= 0.1] = 2
TCW.trd3.grd[TCW3.dry.trd$p > 0.1] = 3
TCW.trd3.grd[is.na(TCW3.dry.trd$p)] = 4

writeRaster(TCW.trd1.grd,".\\NBAR_results3\\TCW1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCW.trd2.grd,".\\NBAR_results3\\TCW2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCW.trd3.grd,".\\NBAR_results3\\TCW3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

#TCG
names(TCG1.dry.trd) <- c("tau", "p", "obs")
names(TCG2.dry.trd) <- c("tau", "p", "obs")
names(TCG3.dry.trd) <- c("tau", "p", "obs")

TCG.trd1.grd = TCG1.dry.trd$tau
TCG.trd1.grd[TCG1.dry.trd$tau < 0 & TCG1.dry.trd$p <= 0.1] = 1
TCG.trd1.grd[TCG1.dry.trd$tau >= 0 & TCG1.dry.trd$p <= 0.1] = 2
TCG.trd1.grd[TCG1.dry.trd$p > 0.1] = 3
TCG.trd1.grd[is.na(TCG1.dry.trd$p)] = 4

TCG.trd2.grd = TCG1.dry.trd$tau
TCG.trd2.grd[TCG2.dry.trd$tau < 0 & TCG2.dry.trd$p <= 0.1] = 1
TCG.trd2.grd[TCG2.dry.trd$tau >= 0 & TCG2.dry.trd$p <= 0.1] = 2
TCG.trd2.grd[TCG2.dry.trd$p > 0.1] = 3
TCG.trd2.grd[is.na(TCG2.dry.trd$p)] = 4

TCG.trd3.grd = TCG3.dry.trd$tau
TCG.trd3.grd[TCG3.dry.trd$tau < 0 & TCG3.dry.trd$p <= 0.1] = 1
TCG.trd3.grd[TCG3.dry.trd$tau >= 0 & TCG3.dry.trd$p <= 0.1] = 2
TCG.trd3.grd[TCG3.dry.trd$p > 0.1] = 3
TCG.trd3.grd[is.na(TCG3.dry.trd$p)] = 4

writeRaster(TCG.trd1.grd,".\\NBAR_results3\\TCG1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCG.trd2.grd,".\\NBAR_results3\\TCG2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCG.trd3.grd,".\\NBAR_results3\\TCG3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

#TCB
names(TCB1.dry.trd) <- c("tau", "p", "obs")
names(TCB2.dry.trd) <- c("tau", "p", "obs")
names(TCB3.dry.trd) <- c("tau", "p", "obs")

TCB.trd1.grd = TCB1.dry.trd$tau
TCB.trd1.grd[TCB1.dry.trd$tau < 0 & TCB1.dry.trd$p <= 0.1] = 1
TCB.trd1.grd[TCB1.dry.trd$tau >= 0 & TCB1.dry.trd$p <= 0.1] = 2
TCB.trd1.grd[TCB1.dry.trd$p > 0.1] = 3
TCB.trd1.grd[is.na(TCB1.dry.trd$p)] = 4

TCB.trd2.grd = TCB1.dry.trd$tau
TCB.trd2.grd[TCB2.dry.trd$tau < 0 & TCB2.dry.trd$p <= 0.1] = 1
TCB.trd2.grd[TCB2.dry.trd$tau >= 0 & TCB2.dry.trd$p <= 0.1] = 2
TCB.trd2.grd[TCB2.dry.trd$p > 0.1] = 3
TCB.trd2.grd[is.na(TCB2.dry.trd$p)] = 4

TCB.trd3.grd = TCB3.dry.trd$tau
TCB.trd3.grd[TCB3.dry.trd$tau < 0 & TCB3.dry.trd$p <= 0.1] = 1
TCB.trd3.grd[TCB3.dry.trd$tau >= 0 & TCB3.dry.trd$p <= 0.1] = 2
TCB.trd3.grd[TCB3.dry.trd$p > 0.1] = 3
TCB.trd3.grd[is.na(TCB3.dry.trd$p)] = 4

writeRaster(TCB.trd1.grd,".\\NBAR_results3\\TCB1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCB.trd2.grd,".\\NBAR_results3\\TCB2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCB.trd3.grd,".\\NBAR_results3\\TCB3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

#TCA
names(TCA1.dry.trd) <- c("tau", "p", "obs")
names(TCA2.dry.trd) <- c("tau", "p", "obs")
names(TCA3.dry.trd) <- c("tau", "p", "obs")

TCA.trd1.grd = TCA1.dry.trd$tau
TCA.trd1.grd[TCA1.dry.trd$tau < 0 & TCA1.dry.trd$p <= 0.1] = 1
TCA.trd1.grd[TCA1.dry.trd$tau >= 0 & TCA1.dry.trd$p <= 0.1] = 2
TCA.trd1.grd[TCA1.dry.trd$p > 0.1] = 3
TCA.trd1.grd[is.na(TCA1.dry.trd$p)] = 4

TCA.trd2.grd = TCA1.dry.trd$tau
TCA.trd2.grd[TCA2.dry.trd$tau < 0 & TCA2.dry.trd$p <= 0.1] = 1
TCA.trd2.grd[TCA2.dry.trd$tau >= 0 & TCA2.dry.trd$p <= 0.1] = 2
TCA.trd2.grd[TCA2.dry.trd$p > 0.1] = 3
TCA.trd2.grd[is.na(TCA2.dry.trd$p)] = 4

TCA.trd3.grd = TCA3.dry.trd$tau
TCA.trd3.grd[TCA3.dry.trd$tau < 0 & TCA3.dry.trd$p <= 0.1] = 1
TCA.trd3.grd[TCA3.dry.trd$tau >= 0 & TCA3.dry.trd$p <= 0.1] = 2
TCA.trd3.grd[TCA3.dry.trd$p > 0.1] = 3
TCA.trd3.grd[is.na(TCA3.dry.trd$p)] = 4

writeRaster(TCA.trd1.grd,".\\NBAR_results3\\TCA1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCA.trd2.grd,".\\NBAR_results3\\TCA2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCA.trd3.grd,".\\NBAR_results3\\TCA3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

#EVI
names(EVI1.dry.trd) <- c("tau", "p", "obs")
names(EVI2.dry.trd) <- c("tau", "p", "obs")
names(EVI3.dry.trd) <- c("tau", "p", "obs")

EVI.trd1.grd = EVI1.dry.trd$tau
EVI.trd1.grd[EVI1.dry.trd$tau < 0 & EVI1.dry.trd$p <= 0.1] = 1
EVI.trd1.grd[EVI1.dry.trd$tau >= 0 & EVI1.dry.trd$p <= 0.1] = 2
EVI.trd1.grd[EVI1.dry.trd$p > 0.1] = 3
EVI.trd1.grd[is.na(EVI1.dry.trd$p)] = 4

EVI.trd2.grd = EVI1.dry.trd$tau
EVI.trd2.grd[EVI2.dry.trd$tau < 0 & EVI2.dry.trd$p <= 0.1] = 1
EVI.trd2.grd[EVI2.dry.trd$tau >= 0 & EVI2.dry.trd$p <= 0.1] = 2
EVI.trd2.grd[EVI2.dry.trd$p > 0.1] = 3
EVI.trd2.grd[is.na(EVI2.dry.trd$p)] = 4

EVI.trd3.grd = EVI3.dry.trd$tau
EVI.trd3.grd[EVI3.dry.trd$tau < 0 & EVI3.dry.trd$p <= 0.1] = 1
EVI.trd3.grd[EVI3.dry.trd$tau >= 0 & EVI3.dry.trd$p <= 0.1] = 2
EVI.trd3.grd[EVI3.dry.trd$p > 0.1] = 3
EVI.trd3.grd[is.na(EVI3.dry.trd$p)] = 4

writeRaster(EVI.trd1.grd,".\\NBAR_results3\\EVI1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(EVI.trd2.grd,".\\NBAR_results3\\EVI2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(EVI.trd3.grd,".\\NBAR_results3\\EVI3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

#NDWI
names(NDWI1.dry.trd) <- c("tau", "p", "obs")
names(NDWI2.dry.trd) <- c("tau", "p", "obs")
names(NDWI3.dry.trd) <- c("tau", "p", "obs")

NDWI.trd1.grd = NDWI1.dry.trd$tau
NDWI.trd1.grd[NDWI1.dry.trd$tau < 0 & NDWI1.dry.trd$p <= 0.1] = 1
NDWI.trd1.grd[NDWI1.dry.trd$tau >= 0 & NDWI1.dry.trd$p <= 0.1] = 2
NDWI.trd1.grd[NDWI1.dry.trd$p > 0.1] = 3
NDWI.trd1.grd[is.na(NDWI1.dry.trd$p)] = 4

NDWI.trd2.grd = NDWI1.dry.trd$tau
NDWI.trd2.grd[NDWI2.dry.trd$tau < 0 & NDWI2.dry.trd$p <= 0.1] = 1
NDWI.trd2.grd[NDWI2.dry.trd$tau >= 0 & NDWI2.dry.trd$p <= 0.1] = 2
NDWI.trd2.grd[NDWI2.dry.trd$p > 0.1] = 3
NDWI.trd2.grd[is.na(NDWI2.dry.trd$p)] = 4

NDWI.trd3.grd = NDWI3.dry.trd$tau
NDWI.trd3.grd[NDWI3.dry.trd$tau < 0 & NDWI3.dry.trd$p <= 0.1] = 1
NDWI.trd3.grd[NDWI3.dry.trd$tau >= 0 & NDWI3.dry.trd$p <= 0.1] = 2
NDWI.trd3.grd[NDWI3.dry.trd$p > 0.1] = 3
NDWI.trd3.grd[is.na(NDWI3.dry.trd$p)] = 4

writeRaster(NDWI.trd1.grd,".\\NBAR_results3\\NDWI1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(NDWI.trd2.grd,".\\NBAR_results3\\NDWI2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(NDWI.trd3.grd,".\\NBAR_results3\\NDWI3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
