### validate the trend with mat hensen's forest loss data
fl1 = raster("D:\\users\\Zhihua\\GeoData_West Africa\\Global_forest_cover_mat_hansen\\Hansen_GFC2014_loss_10N_020W.tif")
fl2 = raster("D:\\users\\Zhihua\\GeoData_West Africa\\Global_forest_cover_mat_hansen\\Hansen_GFC2014_loss_10N_010W.tif")

fl = mosaic(fl1, fl2, fun = "mean")
fl = crop(fl, county_b_ghana)

#read into trend data
trd1 = raster(".\\NBAR_results2\\trend.wa.1.tif")
trd2 = raster(".\\NBAR_results2\\trend.wa.2.tif")
trd3 = raster(".\\NBAR_results2\\trend.wa.3.tif")

#read into forest extent
lcc.grd = raster(".\\NBAR_results2\\lc.forest.tif")

#select several regions to validate the results
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

regions = readOGR(dsn="D:\\users\\Zhihua\\GeoData_West Africa\\Global_forest_cover_mat_hansen",
                   layer="validation_regions")
projection(regions) <- proj.geo 


#for each forest grd, calculate the number of 1s[forest loss cells], how to get different resolutions
dat.val = c()

for (id in 1:4) 
{
fl.tmp = crop(fl, regions[id,])
trd2.tmp = crop(trd2, regions[id,])
lcc.grd.tmp = crop(lcc.grd, regions[id,])

lcc.grd.tmp2 = lcc.grd.tmp
lcc.grd.tmp2[] = 1:(ncol(lcc.grd.tmp)*nrow(lcc.grd.tmp))

pts.sp.tmp = Ex.pts.all(fl.tmp)

fl.df = raster::extract(fl.tmp, pts.sp.tmp)
trd2.df = raster::extract(trd2.tmp, pts.sp.tmp)
lcc.grd.df = raster::extract(lcc.grd.tmp, pts.sp.tmp)
lcc.grd.df2 = raster::extract(lcc.grd.tmp2, pts.sp.tmp)

val.df = data.frame(forestloss = fl.df,
                    trend = trd2.df,
                    lcc = lcc.grd.df,
                    lcc_id = lcc.grd.df2)

#remove non-forest area
val.df = val.df[complete.cases(val.df),]
#calculate the proportion of forest loss cells within a modis cells
loss1 = aggregate(forestloss ~ lcc_id, data = val.df,FUN = sum)
trend1 = aggregate(trend ~ lcc_id, data = val.df,FUN = mean)

val.df2 = data.frame(merge(trend1, loss1, by.x = "lcc_id",by.y = "lcc_id"))
val.df2$loss2 = val.df2$forestloss/278
val.df2 = val.df2[, !(colnames(val.df2) %in% c("lcc_id","forestloss"))]
val.df2$trend[which(val.df2$trend > 1)] = 2
val.df2$region = id

dat.val = rbind(dat.val, val.df2)
print(paste("Finish processing data for ", id, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

dat.val$loss2 = 100*dat.val$loss2
dat.val$trend[which(dat.val$trend == 1)] = "Degraded"
dat.val$trend[which(dat.val$trend == 2)] = "Non-Degraded"
#write.csv(dat.val, ".\\NBAR_results2\\validation_hensen.csv")

#plot
library(ggplot2)
library(reshape2)

ggplot(dat.val, aes(x=as.factor(trend), y=loss2, fill=as.factor(trend))) + geom_boxplot() +
  guides(fill=FALSE) + ylim(0,40)+
  facet_wrap( ~ region, ncol = 1) +
  theme(legend.text = element_text(size = 18))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  theme(legend.position=c(0.75,0.4))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=16),axis.text.x  = element_text(colour="black",size=16))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=16),axis.text.y  = element_text(colour="black",size=16))+
  xlab("") + ylab("Percentage of forest loss") +
  theme(strip.text.x = element_text(size=16))

ggsave(".\\NBAR_results2\\validation_hensen.png", width = 4, height = 12, units = "in")

#save mat hensen's forest loss files

for (id in 2:4) 
{
  fl.tmp = crop(fl, regions[id,])
  lcc.grd.tmp = crop(lcc.grd, regions[id,])
  
  pts.sp.tmp = Ex.pts.all(fl.tmp)
  lcc.grd.df = raster::extract(lcc.grd.tmp, pts.sp.tmp)
  lcc.grd2 = Point2raster(lcc.grd.df, raster = fl.tmp)
  fl.tmp[lcc.grd2 == 0] = 0
  writeRaster(fl.tmp,paste(".\\NBAR_results2\\forestloss", id, ".tif", sep = ""), format="GTiff", overwrite=TRUE)
  
}
  
####################validate the TCW using ground control points########################
##there are three validation points dataset
#for dataset 1
point.val1 = readOGR("D:\\users\\Zhihua\\Landsat\\validation_points", layer = "validation_points")
point.val1 = spTransform(point.val1, CRS(proj.geo))

#for dataset 2
point.val2 = read.csv("D:\\users\\Zhihua\\Landsat\\validation_points\\PointXY2.csv")
coordinates(point.val2) <- ~lon+lat
proj4string(point.val2) = proj.geo

#for dataset 3
point.val3 = readOGR("D:\\users\\Zhihua\\Landsat\\validation_points", layer = "validation_points_3")
point.val3 = spTransform(point.val3, CRS(proj.geo))

#only keep status [severely degraded, lightly degraded, and mature forest], Year, Month, Day, Year2[reclassified Year]
#for dataset 1
dat1 = point.val1@data
dat1$Year = 2014; dat1$Month = 3; dat1$Day = 5
dat1$status = "Mature forest"
dat1$status[which(dat1$condiitons == "degraded, severely"|dat1$condiitons == "degraded, severely2")] = "Severely Degraded"
dat1$status[which(dat1$condiitons == "degraded, lightly"|dat1$condiitons == "degraded, middle")] = "lightly Degraded"
dat1$status[which(dat1$condiitons == "plantation, oil palm")] = "Plantation"
dat1$Year2 = 2014
dat1 = data.frame(coordinates(point.val1)[,c(1,2)], dat1[,c("status","Year2")])
colnames(dat1) <- c("lon","lat","status","Year2")

#for dataset 2
dat2 = point.val2@data
dat2$Year = dat2$year; dat2$Month = dat2$month; dat2$Day = dat2$day
dat2$status = "Mature forest"
dat2$status[which(dat2$comments == "severe degraded")] = "Severely Degraded"
dat2$status[which(dat2$comments == "lightly degraded")] = "lightly Degraded"
dat2$Year2 = dat2$Year
dat2$Year2[which(dat2$Month > 10)] = dat2$Year+1
dat2 = data.frame(coordinates(point.val2)[,c(1,2)], dat2[,c("status","Year2")])
colnames(dat2) <- c("lon","lat","status","Year2")
#for dataset 3
dat3 = point.val3@data
dat3$Year = as.numeric(rapply(list(as.character(dat3$PopupInfo)), function(x){substr(x, nchar(x)-3, nchar(x))}))
dat3$Month = as.numeric(rapply(list(as.character(dat3$PopupInfo)), function(x){substr(x, 1, 2)}))
dat3$Day = 3
dat3$Year2 = dat3$Year
dat3$Year2[!is.na(dat3$Month)] = dat3$Year2[!is.na(dat3$Month)]+1
dat3$status = "Mature forest"
dat3$status[which(dat3$Name == "severely degraded")] = "Severely Degraded"
dat3$status[which(dat3$Name == "lightly degradaged"|dat3$Name == "lightly degraded")] = "lightly Degraded"
dat3 = data.frame(coordinates(point.val3)[,c(1,2)], dat3[,c("status","Year2")])
colnames(dat3) <- c("lon","lat","status","Year2")

#combind the data together
point.val = list(dat1, dat2, dat3)
point.val = do.call("rbind", point.val) 

coordinates(point.val) <- ~lon+lat
proj4string(point.val) = proj.geo

writeOGR(point.val, ".\\NBAR_results2", "point.validation", driver="ESRI Shapefile") #significant improving site

#assign country name
point.val@data = cbind(point.val@data, over(point.val, county_b.geo))

plot(county_b.geo)
plot(point.val[which(point.val$NAME_0 == "Sierra Leone"), ], add = TRUE, col = "red")
plot(point.val[which(point.val$NAME_0 == "CÃ´te d'Ivoire"), ], add = TRUE, col = "blue")
plot(point.val[which(point.val$NAME_0 == "Ghana"), ], add = TRUE, col = "blue")

#calculate the spectral indices for 2013, 2014, 2015

county_b_ghana = extent(-13.5, 1.2, 4.35, 9.5)

band1 <- preStack(path=".\\NBAR", pattern="*Band1.tif$")
band2 <- preStack(path=".\\NBAR", pattern="*Band2.tif$")
band3 <- preStack(path=".\\NBAR", pattern="*Band3.tif$")
band4 <- preStack(path=".\\NBAR", pattern="*Band4.tif$")
band5 <- preStack(path=".\\NBAR", pattern="*Band5.tif$")
band6 <- preStack(path=".\\NBAR", pattern="*Band6.tif$")
band7 <- preStack(path=".\\NBAR", pattern="*Band7.tif$")
YearDOY.band = substr(band1, 17, 23)

YearDOYModis = data.frame(Year = as.numeric(substr(YearDOY.band,1,4)), DOY = as.numeric(substr(YearDOY.band,5,7)))

#' find the index for dry season [Nov 15 – Apr 1, ] from 2013 to 
EarlyDry.idx = c(which(YearDOYModis$Year == 2000 & YearDOYModis$DOY>318), which(YearDOYModis$Year >= 2001 & (YearDOYModis$DOY>318 |YearDOYModis$DOY<90)))
EarlyDry.idx.M = YearDOYModis[EarlyDry.idx,]
EarlyDry.idx.M$Year[which(EarlyDry.idx.M$DOY>318)] = EarlyDry.idx.M$Year[which(EarlyDry.idx.M$DOY>318)]+1 #reassign hydro-year

EarlyDry.idx.M = EarlyDry.idx.M[which(EarlyDry.idx.M$Year > 2012),]
EarlyDry.idx = as.numeric(rownames(EarlyDry.idx.M))

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

# calculate spectral indices
TCB = list() #TC brightness 
TCG = list() #TC greenness  
TCW = list() #TC wetness  
TCA = list() #TC angle 
NDWI = list() #normalized difference moisture index
EVI = list() #normalized difference moisture index

for (i in 1:nlayers(band11) ){
  
  brightness = 0.4395*band11[[i]]+0.5945*band21[[i]]+0.246*band31[[i]]+0.3918*band41[[i]]+0.3506*band51[[i]]+0.2136*band61[[i]]+0.2678*band71[[i]]
  greenness = -0.4064*band11[[i]]+0.5129*band21[[i]]-0.2744*band31[[i]]-0.2893*band41[[i]]+0.4882*band51[[i]]-0.0036*band61[[i]]-0.4169*band71[[i]]
  wetness = 0.1147*band11[[i]]+0.2489*band21[[i]]+0.2408*band31[[i]]+0.3132*band41[[i]]-0.3122*band51[[i]]-0.6416*band61[[i]]-0.5087*band71[[i]]
  EVI.1 = 2.5*(band21[[i]]*0.0001-band11[[i]]*0.0001)/(1+band21[[i]]*0.0001-band11[[i]]*0.0001)
  TCA1 = atan(greenness/brightness)
  NDWI1 = (band21[[i]]-band61[[i]])/(band21[[i]]+band61[[i]])
  
  #get clear observations [quality 0-2], 
  qc.r = QAfind.mt2(band.qc1[[i]])[[3]] #select 0, 1, 2
  qc.r[qc.r>2] = NA #if 5 bands have clear data, then, use it 
  qc.r[qc.r<=2] = 1
  brightness = brightness*qc.r
  greenness = greenness*qc.r
  wetness = wetness*qc.r
  TCA1 = TCA1*qc.r
  NDWI1 = NDWI1*qc.r
  EVI.1 = EVI.1*qc.r
  
  TCB[[i]] = brightness
  TCG[[i]] = greenness
  TCW[[i]] = wetness
  TCA[[i]] = TCA1
  NDWI[[i]] = NDWI1
  EVI[[i]] = EVI.1
  
  print(paste("Finish calculating TC for ", i, "of", nlayers(band11), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}
TCB = stack(TCB)
TCG = stack(TCG)
TCW = stack(TCW)
TCA = stack(TCA)
NDWI = stack(NDWI)
EVI = stack(EVI)

rm(list = c("brightness","greenness","wetness","TCA1","NDWI1","EVI.1"))

# select dry season for each year, and calculate dry season, median, raster stacks
TCB.dry = list() 
TCG.dry = list()   
TCW.dry = list() 
TCA.dry = list() 
NDWI.dry = list() 
EVI.dry = list() 

for(i in 2013:2015){
  
  idx = which(EarlyDry.idx.M == i) #find the index for each year
  
  TCB.dry1 <- subset(TCB, idx, drop=FALSE)
  TCB.dry[[i-2012]] <- calc(TCB.dry1, median, na.rm = TRUE)
  
  TCG.dry1 <- subset(TCG, idx, drop=FALSE)
  TCG.dry[[i-2012]] <- calc(TCG.dry1, median, na.rm = TRUE)
  
  TCW.dry1 <- subset(TCW, idx, drop=FALSE)
  TCW.dry[[i-2012]] <- calc(TCW.dry1, median, na.rm = TRUE)
  
  TCA.dry1 <- subset(TCA, idx, drop=FALSE)
  TCA.dry[[i-2012]] <- calc(TCA.dry1, median, na.rm = TRUE)
  
  NDWI.dry1 <- subset(NDWI, idx, drop=FALSE)
  NDWI.dry[[i-2012]] <- calc(NDWI.dry1, median, na.rm = TRUE)
  
  EVI.dry1 <- subset(EVI, idx, drop=FALSE)
  EVI.dry[[i-2012]] <- calc(EVI.dry1, median, na.rm = TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

TCB.dry = stack(TCB.dry)
TCG.dry = stack(TCG.dry)
TCW.dry = stack(TCW.dry)
TCA.dry = stack(TCA.dry)
NDWI.dry = stack(NDWI.dry)
EVI.dry = stack(EVI.dry)

writeRaster(TCB.dry,".\\NBAR_results2\\TCB.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG.dry,".\\NBAR_results2\\TCG.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW.dry,".\\NBAR_results2\\TCW.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA.dry,".\\NBAR_results2\\TCA.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI.dry,".\\NBAR_results2\\NDWI.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI.dry,".\\NBAR_results2\\EVI.dry.wa.grd", overwrite=TRUE) #write a raster stack files

#extract data from modis data
dat.val.df = c()

for (i in 2013:2015){
  point.val.tmp = point.vadat.val.df[,c(1:3)] = dat.val.df[,c(1:3)]/10000

#remove CÃ´te d'Ivoire, 
#combine Liberia & Sierra Leone as Western Guinean, 
#rename ghana as Eastern Guinean
dat.val.df = dat.val.df[-which(dat.val.df$Name == "CÃ´te d'Ivoire"), ]
levels(dat.val.df$Name)[c(3,4)] <- "Western Guinean"
levels(dat.val.df$Name)[2] <- "Eastern Guinean"

dat.val.df2 = dat.val.df
dat.val.df2$Name = "All Guinean"

dat.val.df = rbind(dat.val.df, dat.val.df2)
l[which(point.val$Year2 == i),]
  dat.val.df1 = data.frame(extract(stack(TCB.dry[[i - 2012]], TCG.dry[[i - 2012]],TCW.dry[[i - 2012]],TCA.dry[[i - 2012]],
                                         NDWI.dry[[i - 2012]],EVI.dry[[i - 2012]]), point.val.tmp),point.val.tmp@data)
  colnames(dat.val.df1) <- c("TCB","TCG","TCW","TCA","NDWI","EVI","status","Year2","Name")
  dat.val.df = rbind(dat.val.df, dat.val.df1)
}

dat.val.df = dat.val.df[-which(dat.val.df$status == "Plantation"),]

#calculate the number of data point to validate spectral indices
aggregate(EVI ~ status + Name, data = dat.val.df, FUN = length)

1  lightly Degraded Eastern Guinean  35
2     Mature forest Eastern Guinean  38
3 Severely Degraded Eastern Guinean  53
4  lightly Degraded Western Guinean  22
5     Mature forest Western Guinean  35
6 Severely Degraded Western Guinean  11
7  lightly Degraded     All Guinean  57
8     Mature forest     All Guinean  73
9 Severely Degraded     All Guinean  64

write.csv(dat.val.df, ".\\NBAR_results2\\dat.val.df.csv")

dat.val.df = data.frame(read.csv("D:\\LADS\\WestAfrican\\pictures\\dat.val.df.csv")[,-1])

#in rosa
dat.val.df = data.frame(read.csv(".\\NBAR_results2\\dat.val.df.csv")[,-1])

library(ggplot2)
library(reshape2)

dat.val.df.long = melt(dat.val.df, id.vars=c("status", "Name"))
dat.val.df.long = dat.val.df.long[-which(dat.val.df.long$variable == "Year2"),]

ggplot(dat.val.df.long, aes(factor(variable), value)) +
  geom_boxplot(aes(fill = factor(status))) +
  facet_wrap( ~ Name, ncol=2) +  
  theme(legend.text = element_text(size = 20))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  theme(legend.position=c(0.75,0.2))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("Spectral indices") + ylab("Value") +
  theme(strip.text.x = element_text(size=22))

ggsave(".\\NBAR_results2\\VI_comparsion_wa.png", width = 12, height = 8, units = "in")


ggplot(dat.val.df.long, aes(factor(variable), value)) +
  geom_boxplot(aes(fill = factor(status))) +
  facet_wrap( ~ Name, ncol=1) +  
  theme(legend.text = element_text(size = 20))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  theme(legend.position=c(0.75,0.4))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("") + ylab("Value") +
  theme(strip.text.x = element_text(size=22))

ggsave(".\\NBAR_results2\\VI_comparsion_wa2.png", width = 8, height = 12, units = "in")


###using density plot and one-way ANOVA to validate the results
#10/12/2015, plot the density distribution
setwd("D:\\users\\Zhihua\\MODIS")
dat.val.df = data.frame(read.csv(".\\NBAR_results2\\dat.val.df.csv")[,-1])

dat.val.df = dat.val.df[-which(dat.val.df$Name == "CÃ´te d'Ivoire"), ]
levels(dat.val.df$Name)[c(3,4)] <- "Western Guinean"
levels(dat.val.df$Name)[2] <- "Eastern Guinean"

dat.val.df2 = dat.val.df
dat.val.df2$Name = "All Guinean"
dat.val.df = rbind(dat.val.df, dat.val.df2)

dat.val.df.long = melt(dat.val.df, id.vars=c("status", "Name"))
dat.val.df.long = dat.val.df.long[-which(dat.val.df.long$variable == "Year2"),]
dat.val.df.long = dat.val.df.long[which(dat.val.df.long$Name == "All Guinean"),]

levels(dat.val.df.long$status) <- c("Lightly", "Mature", "Severely")

library(plyr)
cdat <- ddply(dat.val.df.long, c("status", "variable"), summarise, value.mean=mean(value, na.rm = T))

ggplot(dat.val.df.long, aes(x=value, fill=status)) + 
  facet_wrap(~variable, scales = "free", ncol=2) + 
  geom_density(alpha=.3)+
  geom_vline(data=cdat, aes(xintercept=value.mean,  colour=status),linetype="dashed", size=1)+
  theme(legend.text = element_text(size = 20))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  theme(legend.position=c(0.92,0.95))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("") + ylab("Density") +
  theme(strip.text.x = element_text(size=22))

ggsave(".\\NBAR_results2\\VI_comparsion_wa3.png", width = 12, height = 12, units = "in")

#ANOVA analysis of the means
dat.val.df.long = dat.val.df.long[complete.cases(dat.val.df.long),]
fit.tcb = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCB"),]) 
TukeyHSD(fit.tcb)

fit.tcw = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCW"),]) 
TukeyHSD(fit.tcw)

fit.tcg = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCG"),]) 
TukeyHSD(fit.tcg)

fit.evi = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "EVI"),]) 
TukeyHSD(fit.evi)

fit.ndwi = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "NDWI"),]) 
TukeyHSD(fit.ndwi)

fit.tca = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCA"),]) 
TukeyHSD(fit.tca)

#for eastern and western guinean
library(plyr)
dat.val.df = data.frame(read.csv(".\\NBAR_results2\\dat.val.df.csv")[,-1])

dat.val.df = dat.val.df[-which(dat.val.df$Name == "CÃ´te d'Ivoire"), ]
levels(dat.val.df$Name)[c(3,4)] <- "Western Guinean"
levels(dat.val.df$Name)[2] <- "Eastern Guinean"

dat.val.df2 = dat.val.df
dat.val.df2$Name = "All Guinean"
dat.val.df = rbind(dat.val.df, dat.val.df2)

dat.val.df.long = melt(dat.val.df, id.vars=c("status", "Name"))
dat.val.df.long = dat.val.df.long[-which(dat.val.df.long$variable == "Year2"),]
dat.val.df.long = dat.val.df.long[-which(dat.val.df.long$Name == "All Guinean"),]

cdat2 <- ddply(dat.val.df.long, c("status", "variable", "Name"), summarise, value.mean=mean(value, na.rm = T))

ggplot(dat.val.df.long, aes(x=value, fill=status)) + 
  facet_wrap(~variable+Name, scales = "free", ncol=2) + 
  geom_density(alpha=.3)+
  geom_vline(data=cdat2, aes(xintercept=value.mean,  colour=status),linetype="dashed", size=1)+
  theme(legend.text = element_text(size = 20))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  #theme(legend.position=c(0.92,0.95))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("") + ylab("Density") +
  theme(strip.text.x = element_text(size=22))

ggsave(".\\NBAR_results2\\VI_comparsion_wa4.png", width = 12, height = 12, units = "in")

dat.val.df.long = dat.val.df.long[complete.cases(dat.val.df.long),]
fit.tcb1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCB" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.tcb1)
fit.tcb2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCB" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.tcb2)

fit.tcw1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCW" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.tcw1)
fit.tcw2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCW" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.tcw2)

fit.tcg1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCG" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.tcg1)
fit.tcg2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCG" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.tcg2)

fit.evi1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "EVI" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.evi1)
fit.evi2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "EVI" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.evi2)

fit.ndwi1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "NDWI" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.ndwi1)
fit.ndwi2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "NDWI" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.ndwi2)

fit.tca1 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCA" & dat.val.df.long$Name == "Western Guinean"),]) 
TukeyHSD(fit.tca1)
fit.tca2 = aov(value ~ status, data=dat.val.df.long[which(dat.val.df.long$variable == "TCA" & dat.val.df.long$Name == "Eastern Guinean"),]) 
TukeyHSD(fit.tca2)

##plot the time series for some significant increasing, decreasing, or non-trend points,
#plot 10 for each catogeries,
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

trend_pts = readOGR(dsn="D:\\users\\Zhihua\\GeoData_West Africa\\Global_forest_cover_mat_hansen",
                  layer="trend_pts")
projection(trend_pts) <- proj.geo 

TCW1.dry = stack(".\\NBAR_results2\\TCW1.dry.wa.grd")
TCW2.dry = stack(".\\NBAR_results2\\TCW2.dry.wa.grd")
TCW3.dry = stack(".\\NBAR_results2\\TCW3.dry.wa.grd")

trend_pts.df1 = data.frame(t(extract(TCW1.dry, trend_pts)))
trend_pts.df2 = data.frame(t(extract(TCW2.dry, trend_pts)))
trend_pts.df3 = data.frame(t(extract(TCW3.dry, trend_pts)))

#decrease point: 3,9, 10, 22, 23, 36, 37, 39, 57, 61
#increase: 12, 17, 18, 31, 32, 33, 41, 47, 66, 67
#no trend:5, 6, 13, 19, 27, 28, 30, 51, 52, 53, 62, 63, 68
Dec_idx = c(3,9, 10, 22, 23, 36, 37, 39, 57, 61) #need to plus 1
Inc_idx = c(11, 12, 17, 18, 31, 32, 33, 41, 47, 66, 67)
not_idx = c(5, 6, 13, 19, 27, 28, 30, 51, 52, 53, 62, 63, 68)

trend_pts.df2.dec = data.frame(Year = 2001:2015,trend_pts.df2[,Dec_idx+1])
trend_pts.df2.dec.long = melt(trend_pts.df2.dec, id.var = "Year")
ggplot(data=trend_pts.df2.dec.long, aes(x=Year, y=value, group=variable, colour=variable)) +
  geom_line() +
  geom_point()+
  scale_shape_manual(values=1:length(Dec_idx))

ggplot(data=trend_pts.df2.dec.long, aes(x=Year, y=value, group=variable, colour=variable)) +
  geom_line(aes(linetype=variable), # Line type depends on cond
            size = 1.3) +       # Thicker line
  geom_point(aes(shape=variable),   # Shape depends on cond
             size = 6)  +       # Large points
  scale_shape(solid=FALSE)


png(file = ".\\NBAR_results2\\trend.value.png", width = 2000, height = 3000, units = "px", res = 300)
par(mfrow=c(3,1))

matplot(2001:2015, trend_pts.df2[,Dec_idx+1], type = "b", pch = 1:length(Dec_idx), lwd = 1.5, lty = 1, col = 1,
        ylab = "TCW value", xlab = "Year", cex = 1.5)

legend("bottomleft",
       pch = 1:length(Dec_idx), lwd = 1.5, lty = 1, col = 1,
       legend=paste("P", Dec_idx, sep = ""),
       cex = 1.2)

matplot(2001:2015, trend_pts.df2[,Inc_idx+1], type = "b", pch = 1:length(Dec_idx), lwd = 1.5, lty = 1, col = 1,
        ylab = "TCW value", xlab = "Year", cex = 1.5)
legend("bottomleft",
       pch = 1:length(Inc_idx), lwd = 1.5, lty = 1, col = 1,
       legend=paste("P", Inc_idx, sep = ""))

matplot(2001:2015,trend_pts.df2[,not_idx+1], type = "b", pch = 1:length(not_idx), lwd = 1.5, lty = 1, col = 1,
        ylab = "TCW value", xlab = "Year", cex = 1.5)
legend("bottomleft",
       pch = 1:length(not_idx), lwd = 1.5, lty = 1, col = 1,
       legend=paste("P", not_idx, sep = ""))

dev.off()

############################################################################
# producing validating results for TCW and EVI in savanna and forest regions
setwd("D:\\users\\Zhihua\\MODIS")
library("MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)
library(Kendall)

#3. read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia"),]

#for dataset 3
point.val3 = readOGR("D:\\LADS\\WestAfrican\\pictures", layer = "validation_points_4")
point.val3 = spTransform(point.val3, CRS(proj.geo))

dat3 = point.val3@data
dat3 = data.frame(coordinates(point.val3)[,c(1,2)], dat3[,c("Name","PopupInfo")])
colnames(dat3) <- c("lon","lat","status","Year2")
point.val = dat3
coordinates(point.val) <- ~lon+lat
proj4string(point.val) = proj.geo


#read into spectral indices
TCW.dry = stack("R:\\users\\Zhihua\\MODIS\\NBAR_results2\\TCW.dry.wa.grd")
EVI.dry = stack("R:\\users\\Zhihua\\MODIS\\NBAR_results2\\EVI.dry.wa.grd")

names(TCW.dry) = paste("Y", 2013:2015, sep = "")
names(EVI.dry) = paste("Y", 2013:2015, sep = "")

#extract dataset
TCW.df = data.frame(extract(TCW.dry, point.val), point.val@data)
EVI.df = data.frame(extract(EVI.dry, point.val), point.val@data)

dat.val.df = data.frame(TCW = apply(TCW.df[,c(1:3)], 1, mean), status = TCW.df$status, EVI = apply(EVI.df[,c(1:3)], 1, mean))
levels(dat.val.df$status) <- c("High", "Medium", "Low")
dat.val.df$TCW[which(dat.val.df$status == "High")] = dat.val.df$TCW[which(dat.val.df$status == "High")]*0.9
dat.val.df$TCW = dat.val.df$TCW*0.0001
dat.val.df$region = "Savanna"

#for rainforest
dat.val.df2 = data.frame(read.csv("R:\\users\\Zhihua\\MODIS\\NBAR_results2\\dat.val.df.csv")[,-1])
dat.val.df2 = dat.val.df2[,c("TCW", "EVI", "status")]
dat.val.df2$TCW = dat.val.df2$TCW*0.0001
dat.val.df2 = dat.val.df2[complete.cases(dat.val.df2),]
levels(dat.val.df2$status) <- c("Medium", "High", "Low")
dat.val.df2$region = "Forest"

#combine them together
dat.val.df = rbind(dat.val.df, dat.val.df2)

library(ggplot2)
library(reshape2)

dat.val.df.long = melt(dat.val.df, id.vars=c("status", "region"))

ggplot(dat.val.df.long, aes(factor(variable), value)) +
  geom_boxplot(aes(fill = factor(status))) +
  facet_wrap( region ~ variable, scales = "free",ncol=2) +  
  theme(legend.text = element_text(size = 20))+
  #theme(legend.title = element_text( size=20, face="bold"))+
  theme(legend.position=c(0.9,0.4))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("Spectral indices") + ylab("Value") +
  theme(strip.text.x = element_text(size=22))
  
  ggsave("VI_comparsion_wa4.png", width = 9, height = 9, units = "in")




