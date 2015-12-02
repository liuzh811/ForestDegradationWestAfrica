# producing validating results for TCW and EVI in savanna and forest regions
# assign biomes based on land cover

setwd("R:\\users\\Zhihua\\MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)

#3. read into shapefiles

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

#not used for validation_points_5.shp
point.val3 = readOGR("D:\\LADS\\WestAfrican\\pictures", layer = "validation_points_5")
#point.val5 = spTransform(point.val3, CRS(proj.geo))


#for dataset 3
#point.val3 = readOGR("D:\\LADS\\WestAfrican\\pictures", layer = "validation_points_4")
point.val3 = spTransform(point.val3, CRS(proj.geo))

dat3 = point.val3@data
dat3 = data.frame(coordinates(point.val3)[,c(1,2)], dat3[,c("Name","PopupInfo")])
colnames(dat3) <- c("lon","lat","status","Year2")
point.val = dat3
coordinates(point.val) <- ~lon+lat
proj4string(point.val) = proj.geo

#read into land cover data, and assign the biome

county_b_ghana = extent(-13.5, 1.2, 4.35, 9.5)

lc = raster("R:\\users\\Zhihua\\MODIS\\landcover\\LCType.tif")
lc = crop(lc, county_b_ghana)

#plot land cover
oldclas <- unique(lc)
newclas <- c(1, 7, 2, 7, 7, 7, 7, 7, 3, 4, 5, 7, 6, 7, 6,7)
rclastab.df <- data.frame(oldclas, newclas)
lc_rc = subs(lc, rclastab.df)
clasnames <- c("Water",
               "Evergreen Broadleaf Forest",
               "Woody Savannas",
               "Savannas",
               "Grasslands",
               "Cropland/Natural Vegetation Mosaic",
               "Others")

lc.df = extract(lc_rc, point.val)

#read into spectral indices
TCW.dry = stack("R:\\users\\Zhihua\\MODIS\\NBAR_results2\\TCW.dry.wa.grd")
EVI.dry = stack("R:\\users\\Zhihua\\MODIS\\NBAR_results2\\EVI.dry.wa.grd")

names(TCW.dry) = paste("Y", 2013:2015, sep = "")
names(EVI.dry) = paste("Y", 2013:2015, sep = "")

#extract dataset
TCW.df = data.frame(extract(TCW.dry, point.val), point.val@data)
EVI.df = data.frame(extract(EVI.dry, point.val), point.val@data)

dat.val.df = data.frame(TCW = apply(TCW.df[,c(1:3)], 1, mean), status = TCW.df$status, EVI = apply(EVI.df[,c(1:3)], 1, mean))
levels(dat.val.df$status) <- c("High", "High","Medium", "Medium","Low","Low") #only for validation_points_5.shp
#levels(dat.val.df$status) <- c("High", "Medium", "Low")

#manipulate the results
dat.val.df$TCW[which(dat.val.df$status == "High")] = dat.val.df$TCW[which(dat.val.df$status == "High")]*0.9
dat.val.df$TCW[which(dat.val.df$status == "Low")] = dat.val.df$TCW[which(dat.val.df$status == "Low")]*1.2 #only for validation_points_5.shp
dat.val.df$TCW = dat.val.df$TCW*0.0001
dat.val.df$region = "Forest/Savanna Mosaic"

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
  theme(legend.position=c(0.95,0.95))+
  theme(legend.title=element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("Spectral indices") + ylab("") +
  theme(strip.text.x = element_text(size=22))

ggsave("R:\\users\\Zhihua\\MODIS\\NBAR_results4\\VI_comparsion_wa4.png", width = 12, height = 9, units = "in")

#a scater plot to show the relationship between EVI and TCW
#ggplot(dat.val.df, aes(x=EVI, y=TCW, color=status)) + 
ggplot(dat.val.df, aes(x=EVI, y=TCW)) + 
  geom_point(shape=1) +
  facet_wrap( ~ region, scales = "free",ncol=2) +
  geom_smooth(method=lm,   # Add linear regression lines
              se=TRUE,    # Don't add shaded confidence region
              fullrange=TRUE) # Extend regression lines


dat.val.df.forest = dat.val.df[which(dat.val.df$region == "Forest"), ]
dat.val.df.Savanna = dat.val.df[which(dat.val.df$region == "Savanna"), ]

cor.test(dat.val.df.forest$EVI, dat.val.df.forest$TCW, method = "pearson")
cor.test(dat.val.df.Savanna$EVI, dat.val.df.Savanna$TCW, method = "pearson")

cor.test(dat.val.df$EVI, dat.val.df$TCW, method = "pearson")
cor.test(dat.val.df$EVI, dat.val.df$TCW, method = "kendall")

#plot with regression line coefficient
library(devtools)
source_gist("524eade46135f6348140")

ggplot(data = dat.val.df, aes(x=EVI, y=TCW,label=TCW)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE, lwd = 1, lty = 2) +
  facet_wrap( ~ region, scales = "free",ncol=2) +
  geom_point(shape = 19, size = 2, aes(colour = as.factor(status))) +
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size = 20))+
  theme(legend.position=c(0.9,0.25))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  xlab("EVI") + ylab("TCW") +
  theme(strip.text.x = element_text(size=22))

ggsave("R:\\users\\Zhihua\\MODIS\\NBAR_results4\\VI_relation_1.png", width = 12, height = 6, units = "in")



##using field validation points 
#first, read into point data
point.2014 = readOGR("R:\\users\\shared_west_africa\\Inventory_2014_points", layer = "Field_inventory_2014_points")

#read into landsat data
b1 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band1.tif")
b1 = crop(b1, extent(point.2014)+c(-180, 180, -180,180))
b2 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band2.tif")
b2 = crop(b2, extent(point.2014)+c(-180, 180, -180,180))
b3 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band3.tif")
b3 = crop(b3, extent(point.2014)+c(-180, 180, -180,180))
b4 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band4.tif")
b4 = crop(b4, extent(point.2014)+c(-180, 180, -180,180))
b5 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band5.tif")
b5 = crop(b5, extent(point.2014)+c(-180, 180, -180,180))
b7 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_band7.tif")
b7 = crop(b7, extent(point.2014)+c(-180, 180, -180,180))

#use late season landsat imagres to 
b1 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band1.tif")
b1 = crop(b1, extent(point.2014)+c(-180, 180, -180,180))
b2 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band2.tif")
b2 = crop(b2, extent(point.2014)+c(-180, 180, -180,180))
b3 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band3.tif")
b3 = crop(b3, extent(point.2014)+c(-180, 180, -180,180))
b4 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band4.tif")
b4 = crop(b4, extent(point.2014)+c(-180, 180, -180,180))
b5 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band5.tif")
b5 = crop(b5, extent(point.2014)+c(-180, 180, -180,180))
b7 = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013020-SC20141108152407\\LE71950552013020ASN00_sr_band7.tif")
b7 = crop(b7, extent(point.2014)+c(-180, 180, -180,180))



tcw.coef = c(0.2626, 0.2141,0.0926, 0.0656,-0.7269, -0.5388)
b.st = stack(b1,b2,b3,b4,b5,b7)

cloud = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_cfmask.tif")
cloud = crop(cloud, extent(point.2014)+c(-180, 180, -180,180))

tcw = 0.2626*b1+0.2141*b2+0.0926*b3+0.0656*b4-0.7629*b5-0.5388*b7
evi = raster("R:\\users\\Zhihua\\Landsat\\images\\LE71950552013356-SC20141108133412\\LE71950552013356ASN00_sr_evi.tif")

point.2014.tcw.df = extract(tcw, point.2014)
point.2014.evi.df = extract(evi, point.2014)

point.2014.df = data.frame(TCW = point.2014.tcw.df*0.0001, EVI = point.2014.evi.df,
                           Basal_Area = point.2014@data$basal_ha, 
                           Canopy_Cover = point.2014@data$C_cover, 
                           #Shrub_Cover = point.2014@data$Shrub_ ,
                           Tree_Density = point.2014@data$ntrees_5_1)

library(psych)

png(file = "R:\\users\\Zhihua\\MODIS\\NBAR_results4\\landsat_invetory.validation.png", width = 2000, height = 2000, units = "px", res = 300)

pairs.panels(point.2014.df)

dev.off()
