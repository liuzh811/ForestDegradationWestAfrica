
setwd("D:\\users\\Zhihua\\MODIS")
library("MODIS")
library(rgdal)
library(raster)
library(rasterVis)

#3. read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

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

# band albedo quality
band.qc = preStack(path=".\\NBAR_QC", pattern="*Albedo_Band_Quality.tif$")
YearDOY.band.qc = substr(band.qc , 20, 26)
band.qc = stack(band.qc[EarlyDry.idx])

#crop the dataset to study area
band.qc1 = list()

for(i in 1:nlayers(band.qc)){
  
  band.qc.1 = crop(band.qc[[i]], county_b_ghana)
  band.qc1[[i]] = band.qc.1
  
  print(paste("Finish Cropping data for ", i, "of", nlayers(band.qc), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

band.qc1 = stack(band.qc1)

# calculate wetness indices, using different level of QC values

list1 = list() # select 0,1
list2 = list() # select 0,1,2
list3 = list() # select 0,1,2,3

for (i in 1:nlayers(band.qc1) ){
  
  qc.r = QAfind.mt2(band.qc1[[i]])
  
  #select 0,1
  qc.r1 = qc.r[[2]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  list1[[i]] = qc.r1
  
  #select 0,1,2
  qc.r1 = qc.r[[3]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  list2[[i]] = qc.r1
  #select 0,1,2,3
  qc.r1 = qc.r[[4]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  list3[[i]] = qc.r1
  
  
  print(paste("Finish calculating TC for ", i, "of", nlayers(band.qc1), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

list1 = stack(list1)
list2 = stack(list2)
list3 = stack(list3)

# calculate number of non-NA 

length.nonna = function(x){length(which(!is.na(x)))}

list1.dry = list() 
list2.dry = list()
list3.dry = list() 
for(i in 2001:2015){
  
  idx = which(EarlyDry.idx.M==i) #find the index for each year
  list1.dry1 <- subset(list1, idx, drop=FALSE)
  list1.dry[[i-2000]] <- calc(list1.dry1,length.nonna)
  
  list2.dry1 <- subset(list2, idx, drop=FALSE)
  list2.dry[[i-2000]] <- calc(list2.dry1,length.nonna)
  
  list3.dry1 <- subset(list3, idx, drop=FALSE)
  list3.dry[[i-2000]] <- calc(list3.dry1,length.nonna)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

list1.dry = stack(list1.dry)
list2.dry = stack(list2.dry)
list3.dry = stack(list3.dry)


writeRaster(list1.dry,".\\NBAR_results2\\Obs.num.level1.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(list2.dry,".\\NBAR_results2\\Obs.num.level2.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(list3.dry,".\\NBAR_results2\\Obs.num.level3.wa.grd", overwrite=TRUE) #write a raster stack files

lcc.grd = raster(".\\NBAR_results2\\lc.forest.tif")
list1.dry[lcc.grd != 1] = NA
list2.dry[lcc.grd != 1] = NA
list3.dry[lcc.grd != 1] = NA

names(list1.dry) <- paste("Y", 2001:2015, sep = "")
names(list2.dry) <- paste("Y", 2001:2015, sep = "")
names(list3.dry) <- paste("Y", 2001:2015, sep = "")

list1.dry.mn = calc(list1.dry, mean, na.rm = TRUE)
list2.dry.mn = calc(list2.dry, mean, na.rm = TRUE)
list3.dry.mn = calc(list3.dry, mean, na.rm = TRUE)

list.dry.mn = stack(list1.dry.mn, list2.dry.mn, list3.dry.mn)
names(list.dry.mn) <- c("QC1","QC2","QC3")

png(file = ".\\NBAR_results2\\Obs.num.mean.png", width = 3000, height = 5000, units = "px", res = 300)

recls2plot(list.dry.mn, threshold = c(2, 4, 6, 8, 10, 12), color = "Reds")

dev.off()

png(file = ".\\NBAR_results2\\Obs.num.qc.level1.png", width = 5000, height = 5000, units = "px", res = 300)

recls2plot(list1.dry, threshold = c(2, 4, 6, 8, 10, 12), color = "Reds")

dev.off()

png(file = ".\\NBAR_results2\\Obs.num.qc.level2.png", width = 5000, height = 5000, units = "px", res = 300)

recls2plot(list2.dry, threshold = c(2, 4, 6, 8, 10, 12), color = "Reds")

dev.off()

png(file = ".\\NBAR_results2\\Obs.num.qc.level3.png", width = 5000, height = 5000, units = "px", res = 300)

recls2plot(list3.dry, threshold = c(2, 4, 6, 8, 10, 12), color = "Reds")

dev.off()



