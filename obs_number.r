
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


##calculate the frequency of clear obs for each date
obs.date1 = list()
obs.date2 = list()
obs.date3 = list()

for(i in 1:18){ #a total of 18 date, 15 years
  
  idx = (0:14)*18+i
  
  list1.dry1 <- subset(list1, idx, drop=FALSE)
  obs.date1[[i]] <- calc(list1.dry1,length.nonna)
  
  list2.dry1 <- subset(list2, idx, drop=FALSE)
  obs.date2[[i]] <- calc(list2.dry1,length.nonna)
  
  list3.dry1 <- subset(list3, idx, drop=FALSE)
  obs.date3[[i]] <- calc(list3.dry1,length.nonna)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

obs.date1 = stack(obs.date1)
obs.date2 = stack(obs.date2)
obs.date3 = stack(obs.date3)

writeRaster(obs.date1,".\\NBAR_results2\\Obs.num.date.level1.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(obs.date2,".\\NBAR_results2\\Obs.num.date.level2.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(obs.date3,".\\NBAR_results2\\Obs.num.date.level3.wa.grd", overwrite=TRUE) #write a raster stack files

obs.date1[lcc.grd != 1] = NA
obs.date2[lcc.grd != 1] = NA
obs.date3[lcc.grd != 1] = NA

names(obs.date1) <- paste("DOY", substr(names(band.qc1)[c(1:18)], 14, 16), sep = "")
names(obs.date2) <- paste("DOY", substr(names(band.qc1)[c(1:18)], 14, 16), sep = "")
names(obs.date3) <- paste("DOY", substr(names(band.qc1)[c(1:18)], 14, 16), sep = "")


p.strip <- list(cex=1.1, lines=2)

png(file = ".\\NBAR_results2\\Obs.num.DOY.qc.level1.png", width = 5000, height = 5000, units = "px", res = 300)

levelplot(obs.date1, par.settings=BuRdTheme(),
          #maxpixels = nrow(obs.date1)*ncol(obs.date1),
          colorkey= list(labels= list(cex = 1.5)),
          scales=list(x=list(cex=1.1),y=list(cex=1.1)), 
          par.strip.text=p.strip) +
  layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()

png(file = ".\\NBAR_results2\\Obs.num.DOY.qc.level2.png", width = 5000, height = 5000, units = "px", res = 300)

levelplot(obs.date2, par.settings=BuRdTheme(),
          #maxpixels = nrow(obs.date1)*ncol(obs.date1),
          colorkey= list(labels= list(cex = 1.5)),
          scales=list(x=list(cex=1.1),y=list(cex=1.1)), 
          par.strip.text=p.strip) +
  layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()

png(file = ".\\NBAR_results2\\Obs.num.DOY.qc.level3.png", width = 5000, height = 5000, units = "px", res = 300)

levelplot(obs.date3, par.settings=BuRdTheme(),
          #maxpixels = nrow(obs.date1)*ncol(obs.date1),
          colorkey= list(labels= list(cex = 1.5)),
          scales=list(x=list(cex=1.1),y=list(cex=1.1)), 
          par.strip.text=p.strip) +
  layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()

#############plot number of year to calculate the trend
TCW1.dry.trd = stack(".\\NBAR_results2\\TCW1.dry.trend.wa.grd")
TCW2.dry.trd = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.grd")
TCW3.dry.trd = stack(".\\NBAR_results2\\TCW3.dry.trend.wa.grd")
lcc.grd = raster(".\\NBAR_results2\\lc.forest.tif")

Numb = stack(TCW1.dry.trd[[3]],TCW2.dry.trd[[3]],TCW3.dry.trd[[3]])
Numb.recls = recls2(Numb,threshold = c(4,6,8,10,12))
Numb.recls[lcc.grd != 1] = NA

#set color scheme
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

names(Numb.recls) <-c("Better.than.1","Better.than.2","Better.than.3")

png(file = ".\\NBAR_results2\\trend.numb.wa.png", width = 2000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.3, lines=2)
levelplot(Numb.recls,
          maxpixels = nrow(Numb.recls)*ncol(Numb.recls),
          at= breaks2_change, margin=FALSE,
          col.regions= clscolor_change,
          colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          par.strip.text=p.strip) +
  layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()


