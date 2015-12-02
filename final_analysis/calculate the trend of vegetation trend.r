# 11/26/2015
# reanalysis based on the WWF regions
# only for the level 2

setwd("D:\\users\\Zhihua\\MODIS")
library("MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)
library(Kendall)

#3. read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
proj.sin = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#highlight some small regions
regions = readOGR(dsn="D:/users/Zhihua/MODIS/NBAR_results3", layer="Val_regions")
projection(regions) <- proj.geo 


#to speed up the processing, devided to study area into four subregions
ext = extent(county_b_ghana)
ext.subregion = list(c(ext@xmin, (ext@xmin+ext@xmax)/2, ext@ymin, (ext@ymin+ext@ymax)/2),
                     c((ext@xmin+ext@xmax)/2, ext@xmax, ext@ymin, (ext@ymin+ext@ymax)/2),
                     c(ext@xmin, (ext@xmin+ext@xmax)/2, (ext@ymin+ext@ymax)/2, ext@ymax),
                     c((ext@xmin+ext@xmax)/2, ext@xmax, (ext@ymin+ext@ymax)/2, ext@ymax))

###################################################################################
## Section 1: calculate the trend of vegetation trend
## And also the number of clear obs
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

#for (subregion in 1:length(ext.subregion)){  #start of region loop
  subregion = 1
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
  
  band1.1 = crop(band1[[i]], ext.subregion[[subregion]])
  band2.1 = crop(band2[[i]], ext.subregion[[subregion]])
  band3.1 = crop(band3[[i]], ext.subregion[[subregion]])
  band4.1 = crop(band4[[i]], ext.subregion[[subregion]])
  band5.1 = crop(band5[[i]], ext.subregion[[subregion]])
  band6.1 = crop(band6[[i]], ext.subregion[[subregion]])
  band7.1 = crop(band7[[i]], ext.subregion[[subregion]])
  band.qc.1 = crop(band.qc[[i]], ext.subregion[[subregion]])
  
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

TCW2 = list() # select 0,1,2
TCG2 = list() #TC greenness
TCB2 = list() #TC brightness
TCA2 = list() #TC angle 
NDWI2 = list() #normalized difference moisture index
EVI2 = list() #normalized difference moisture index

list2 = list() # for calculating the number of clear obs

for (i in 1:nlayers(band11) ){
  
  #calculate indices
  brightness = 0.4395*band11[[i]]+0.5945*band21[[i]]+0.246*band31[[i]]+0.3918*band41[[i]]+0.3506*band51[[i]]+0.2136*band61[[i]]+0.2678*band71[[i]]
  greenness = -0.4064*band11[[i]]+0.5129*band21[[i]]-0.2744*band31[[i]]-0.2893*band41[[i]]+0.4882*band51[[i]]-0.0036*band61[[i]]-0.4169*band71[[i]]
  wetness = 0.1147*band11[[i]]+0.2489*band21[[i]]+0.2408*band31[[i]]+0.3132*band41[[i]]-0.3122*band51[[i]]-0.6416*band61[[i]]-0.5087*band71[[i]]
  EVI.1 = 2.5*(band21[[i]]*0.0001-band11[[i]]*0.0001)/(1+band21[[i]]*0.0001-band11[[i]]*0.0001)
  TCA.1 = atan(greenness/brightness)
  NDWI.1 = (band21[[i]]-band61[[i]])/(band21[[i]]+band61[[i]])
  
  qc.r = QAfind.mt2(band.qc1[[i]])
  
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
  list2[[i]] = qc.r1
  
  
  print(paste("Finish calculating TC for ", i, "of", nlayers(band11), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

rm(list = c("brightness","greenness","wetness","EVI.1","TCA.1","NDWI.1","qc.r","qc.r1"))
rm(list = c("brightness1","greenness1","wetness1","EVI.11","TCA11","NDWI11"))

TCW2 = stack(TCW2)
TCG2 = stack(TCG2)
TCB2 = stack(TCB2)
TCA2 = stack(TCA2)
EVI2 = stack(EVI2)
NDWI2 = stack(NDWI2)
list2 = stack(list2)

# calculate annual median TCW
TCW2.dry = list()
TCG2.dry = list() #TC greenness
TCB2.dry = list() #TC brightness
TCA2.dry = list() #TC angle 
NDWI2.dry = list() #normalized difference moisture index
EVI2.dry = list() #normalized difference moisture index


for(i in 2001:2015){
  
  idx = which(EarlyDry.idx.M==i) #find the index for each year
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
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

TCW2.dry = stack(TCW2.dry)
TCB2.dry = stack(TCB2.dry)
TCG2.dry = stack(TCG2.dry)
TCA2.dry = stack(TCA2.dry)
EVI2.dry = stack(EVI2.dry)
NDWI2.dry = stack(NDWI2.dry)


#write.raster
writeRaster(TCW2.dry,paste(".\\NBAR_results4\\TCW2.dry.wa.subregion.", subregion, ".grd", sep = ""),overwrite=TRUE) #write a raster stack files
writeRaster(TCB2.dry,paste(".\\NBAR_results4\\TCB2.dry.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(TCG2.dry,paste(".\\NBAR_results4\\TCG2.dry.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(TCA2.dry,paste(".\\NBAR_results4\\TCA2.dry.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(EVI2.dry,paste(".\\NBAR_results4\\EVI2.dry.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(NDWI2.dry,paste(".\\NBAR_results4\\NDWI2.dry.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files

#### calculate number of clear obs for each year;  ############
# calculate number of non-NA 

length.nonna = function(x){length(which(!is.na(x)))}

list2.dry = list()
for(i in 2001:2015){
  
  idx = which(EarlyDry.idx.M==i) #find the index for each year

  list2.dry1 <- subset(list2, idx, drop=FALSE)
  list2.dry[[i-2000]] <- calc(list2.dry1,length.nonna)
  

  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

list2.dry = stack(list2.dry)

writeRaster(list2.dry,paste(".\\NBAR_results4\\Obs.num.level2.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files


####calculate number of clear obs for each day of year  ############
obs.date2 = list()
for(i in 1:18){ #a total of 18 date, 15 years
 
  idx = (0:14)*18+i
  
  list2.dry1 <- subset(list2, idx, drop=FALSE)
  obs.date2[[i]] <- calc(list2.dry1,length.nonna)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

obs.date2 = stack(obs.date2)

writeRaster(obs.date2, paste(".\\NBAR_results4\\Obs.num.date.level2.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files

#calculate trends
# test.pt = extract(TCW1.dry, 1) #test function
TCW2.dry.trd <- calc(TCW2.dry, mk1)
TCB2.dry.trd <- calc(TCB2.dry, mk1)
TCG2.dry.trd <- calc(TCG2.dry, mk1)
TCA2.dry.trd <- calc(TCA2.dry, mk1)
EVI2.dry.trd <- calc(EVI2.dry, mk1)
NDWI2.dry.trd <- calc(NDWI2.dry, mk1)

#save the results
writeRaster(TCW2.dry.trd,paste(".\\NBAR_results4\\TCW2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(TCB2.dry.trd,paste(".\\NBAR_results4\\TCB2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(TCG2.dry.trd,paste(".\\NBAR_results4\\TCG2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(TCA2.dry.trd,paste(".\\NBAR_results4\\TCA2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(EVI2.dry.trd,paste(".\\NBAR_results4\\EVI2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
writeRaster(NDWI2.dry.trd,paste(".\\NBAR_results4\\NDWI2.dry.trend.wa.subregion.", subregion, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files


#classified into different 
#TCW
names(TCW2.dry.trd) <- c("tau", "p", "obs")

TCW.trd2.grd = TCW2.dry.trd$tau
TCW.trd2.grd[TCW2.dry.trd$tau < 0 & TCW2.dry.trd$p <= 0.1] = 1
TCW.trd2.grd[TCW2.dry.trd$tau >= 0 & TCW2.dry.trd$p <= 0.1] = 2
TCW.trd2.grd[TCW2.dry.trd$p > 0.1] = 3
TCW.trd2.grd[is.na(TCW2.dry.trd$p)] = 4

writeRaster(TCW.trd2.grd,paste(".\\NBAR_results4\\TCW2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#TCG
names(TCG2.dry.trd) <- c("tau", "p", "obs")

TCG.trd2.grd = TCG2.dry.trd$tau
TCG.trd2.grd[TCG2.dry.trd$tau < 0 & TCG2.dry.trd$p <= 0.1] = 1
TCG.trd2.grd[TCG2.dry.trd$tau >= 0 & TCG2.dry.trd$p <= 0.1] = 2
TCG.trd2.grd[TCG2.dry.trd$p > 0.1] = 3
TCG.trd2.grd[is.na(TCG2.dry.trd$p)] = 4

writeRaster(TCG.trd2.grd,paste(".\\NBAR_results4\\TCG2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#TCB
names(TCB2.dry.trd) <- c("tau", "p", "obs")

TCB.trd2.grd = TCB2.dry.trd$tau
TCB.trd2.grd[TCB2.dry.trd$tau < 0 & TCB2.dry.trd$p <= 0.1] = 1
TCB.trd2.grd[TCB2.dry.trd$tau >= 0 & TCB2.dry.trd$p <= 0.1] = 2
TCB.trd2.grd[TCB2.dry.trd$p > 0.1] = 3
TCB.trd2.grd[is.na(TCB2.dry.trd$p)] = 4

writeRaster(TCB.trd2.grd,paste(".\\NBAR_results4\\TCB2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#TCA
names(TCA2.dry.trd) <- c("tau", "p", "obs")

TCA.trd2.grd = TCA2.dry.trd$tau
TCA.trd2.grd[TCA2.dry.trd$tau < 0 & TCA2.dry.trd$p <= 0.1] = 1
TCA.trd2.grd[TCA2.dry.trd$tau >= 0 & TCA2.dry.trd$p <= 0.1] = 2
TCA.trd2.grd[TCA2.dry.trd$p > 0.1] = 3
TCA.trd2.grd[is.na(TCA2.dry.trd$p)] = 4

writeRaster(TCA.trd2.grd,paste(".\\NBAR_results4\\TCA2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#EVI
names(EVI2.dry.trd) <- c("tau", "p", "obs")

EVI.trd2.grd = EVI2.dry.trd$tau
EVI.trd2.grd[EVI2.dry.trd$tau < 0 & EVI2.dry.trd$p <= 0.1] = 1
EVI.trd2.grd[EVI2.dry.trd$tau >= 0 & EVI2.dry.trd$p <= 0.1] = 2
EVI.trd2.grd[EVI2.dry.trd$p > 0.1] = 3
EVI.trd2.grd[is.na(EVI2.dry.trd$p)] = 4

writeRaster(EVI.trd2.grd,paste(".\\NBAR_results4\\EVI2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#NDWI
names(NDWI2.dry.trd) <- c("tau", "p", "obs")

NDWI.trd2.grd = NDWI2.dry.trd$tau
NDWI.trd2.grd[NDWI2.dry.trd$tau < 0 & NDWI2.dry.trd$p <= 0.1] = 1
NDWI.trd2.grd[NDWI2.dry.trd$tau >= 0 & NDWI2.dry.trd$p <= 0.1] = 2
NDWI.trd2.grd[NDWI2.dry.trd$p > 0.1] = 3
NDWI.trd2.grd[is.na(NDWI2.dry.trd$p)] = 4

writeRaster(NDWI.trd2.grd,paste(".\\NBAR_results4\\NDWI2.dry.trend.wa.subregion.", subregion, ".tif", sep = ""), format="GTiff", overwrite=TRUE)

#sort( sapply(ls(),function(x){object.size(get(x))}))

#free up some memory
rm(list = c("TCW2","TCG2","TCB2","TCA2","EVI2","NDWI2","list2"))
rm(list = c("band11","band21","band31","band41","band51","band61","band71","band.qc1"))
rm(list = c("TCA2.dry1","TCW2.dry1","TCG2.dry1","TCB2.dry1","EVI2.dry1","NDWI2.dry1")) 

#moasic the results::::trend
list.indices = c("TCW", "EVI", "TCA","TCB","TCG","NDWI")
for (i in 1:length(list.indices)){
  tmp1 = raster(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.trend.wa.subregion.1.tif", sep = ""))
  tmp2 = raster(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.trend.wa.subregion.2.tif", sep = ""))
  tmp3 = raster(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.trend.wa.subregion.3.tif", sep = ""))
  tmp4 = raster(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.trend.wa.subregion.4.tif", sep = ""))
  
  tmp = mosaic(tmp1, tmp2, tmp3, tmp4,fun=mean)
  
  writeRaster(tmp,paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.trend.wa.tif", sep = ""), format="GTiff", overwrite=TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#moasic the results::::annual indices

for (i in 1:length(list.indices)){
  tmp1 = stack(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.wa.subregion.1.grd", sep = ""))
  tmp2 = stack(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.wa.subregion.2.grd", sep = ""))
  tmp3 = stack(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.wa.subregion.3.grd", sep = ""))
  tmp4 = stack(paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.wa.subregion.4.grd", sep = ""))
  
  tmp = mosaic(tmp1, tmp2, tmp3, tmp4,fun=mean)
  
  writeRaster(tmp,paste(".\\NBAR_results4\\", list.indices[[i]],"2.dry.wa.grd", sep = ""), overwrite=TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}



