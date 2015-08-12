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

forest_reserve= readOGR("D:\\users\\Zhihua\\MODIS\\Reserves_corrected_Ghana", layer = "Forest_Reserves")
forest_reserve.geo= spTransform(forest_reserve, CRS(proj.geo))

protect = readOGR(dsn="D:\\users\\Zhihua\\GeoData_West Africa\\protectarea",layer="protected_areas_west_africa")

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

for (i in 1:nlayers(band11) ){
  
 wetness = 0.1147*band11[[i]]+0.2489*band21[[i]]+0.2408*band31[[i]]+0.3132*band41[[i]]-0.3122*band51[[i]]-0.6416*band61[[i]]-0.5087*band71[[i]]
  qc.r = QAfind.mt2(band.qc1[[i]])
 
 #select 0,1
  qc.r1 = qc.r[[2]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  wetness1 = wetness*qc.r1
  TCW1[[i]] = wetness1

 #select 0,1,2
  qc.r1 = qc.r[[3]] 
  qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
  qc.r1[qc.r1<=2] = 1
  wetness1 = wetness*qc.r1
  TCW2[[i]] = wetness1
 #select 0,1,2,3
 qc.r1 = qc.r[[4]] 
 qc.r1[qc.r1>2] = NA #if 5 bands have clear data, then, use it 
 qc.r1[qc.r1<=2] = 1
 wetness1 = wetness*qc.r1
 TCW3[[i]] = wetness1
 
  print(paste("Finish calculating TC for ", i, "of", nlayers(band11), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

TCW1 = stack(TCW1)
TCW2 = stack(TCW2)
TCW3 = stack(TCW3)

# calculate annual median TCW
TCW1.dry = list() 
TCW2.dry = list()
TCW3.dry = list() 
for(i in 2001:2015){
  
  idx = which(EarlyDry.idx.M==i) #find the index for each year
  TCW1.dry1 <- subset(TCW1, idx, drop=FALSE)
  TCW1.dry[[i-2000]] <- calc(TCW1.dry1, median, na.rm = TRUE)

  TCW2.dry1 <- subset(TCW2, idx, drop=FALSE)
  TCW2.dry[[i-2000]] <- calc(TCW2.dry1, median, na.rm = TRUE)
  
  TCW3.dry1 <- subset(TCW3, idx, drop=FALSE)
  TCW3.dry[[i-2000]] <- calc(TCW3.dry1, median, na.rm = TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

TCW1.dry = stack(TCW1.dry)
TCW2.dry = stack(TCW2.dry)
TCW3.dry = stack(TCW3.dry)

#write.raster
writeRaster(TCW1.dry,".\\NBAR_results2\\TCW1.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW2.dry,".\\NBAR_results2\\TCW2.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW3.dry,".\\NBAR_results2\\TCW3.dry.wa.grd", overwrite=TRUE) #write a raster stack files

#calculate trends
# test.pt = extract(TCW1.dry, 1) #test function
TCW1.dry.trd <- calc(TCW1.dry, mk1)
TCW2.dry.trd <- calc(TCW2.dry, mk1)
TCW3.dry.trd <- calc(TCW3.dry, mk1)

#save the results
writeRaster(TCW1.dry.trd,".\\NBAR_results2\\TCW1.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW2.dry.trd,".\\NBAR_results2\\TCW2.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW3.dry.trd,".\\NBAR_results2\\TCW3.dry.trend.wa.grd", overwrite=TRUE) #write a raster stack files


#only retain significant trend area
TCW1.dry.trd[[1]][TCW1.dry.trd[[2]] > 0.1] = NA
TCW2.dry.trd[[1]][TCW2.dry.trd[[2]] > 0.1] = NA
TCW3.dry.trd[[1]][TCW3.dry.trd[[2]] > 0.1] = NA

writeRaster(TCW1.dry.trd[[1]],".\\NBAR_results2\\TCW1.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCW2.dry.trd[[1]],".\\NBAR_results2\\TCW2.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)
writeRaster(TCW3.dry.trd[[1]],".\\NBAR_results2\\TCW3.dry.trend.wa.tif", format="GTiff", overwrite=TRUE)

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

#clasnames_change <- c("< 0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6",">0.6")
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
