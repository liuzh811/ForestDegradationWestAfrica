# spatial variability of trend
# analyze trend on differnt countries, wwf ecoregions, land cover, and protected status

#load libraries
library("raster")
library("dismo")
library("sp")
library("rgdal")
library("reshape2")
library("ggplot2")

setwd("R:/users/Zhihua/MODIS")

TCW.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/TCW2.dry.trend.wa2.tif") #trend after moving window analysis
EVI.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/EVI2.dry.trend.wa2.tif")
lc_rc = raster(".\\NBAR_results4\\lc_rc.wa.tif")
county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")

#read into shapefiles, countries
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                    county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

county_b_ghana.grd = rasterize(county_b_ghana, county_b_ghana.r)
county_b_ghana.grd = county_b_ghana.grd*county_b_ghana.r #1:"CÃ´te d'Ivoire", 2:"Ghana", 3:"Guinea", 4:Liberia, 5:Sierra Leone

#protected area #1: protected-Reserve, 2: protected-Eco-Reserve 3:non-protected
protected = readOGR(dsn="R:\\users\\Zhihua\\GeoData_West Africa\\protectarea",layer="protected_areas_west_africa2")

# protected@data$status = 1 #
# protected@data$status[which(protected@data$IUCN_CAT == "Ia" | 
#                               protected@data$IUCN_CAT == "II" | protected@data$IUCN_CAT == "IV")] = 2

# protected.grd = rasterize(protected, county_b_ghana.r, field = "status")
protected.grd = rasterize(protected, county_b_ghana.r, field = "status2")
protected.grd = protected.grd*county_b_ghana.r #protected
protected.grd[is.na(protected.grd)] = 3  #3 non-protected
protected.grd = protected.grd*county_b_ghana.r 

#read into ecoregion; ECO_NAME2: 1: west guinea forest; 2: east guinea foreste; 3: mosaic, 4: savvannas,5:other
eco.sp = readOGR(dsn="R:\\users\\Zhihua\\WWF_ecoregion\\official_teow\\official",layer="wwf_terr_ecos1")
eco.sp.grd = rasterize(eco.sp, county_b_ghana.r, field = eco.sp@data$ECO_NAME2)
eco.sp.grd = eco.sp.grd*county_b_ghana.r

#plot spatial variability of trend
# by country
lc_rc1 = county_b_ghana.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-16,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = county_b_ghana.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-16,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                 data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))))

lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("CÃ´te d'Ivoire", "Ghana", "Guinea", "Liberia", "Sierra Leone")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover") + 
  #theme(legend.position="none")+
  theme(legend.position=c(0.075,0.4))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.text.y = element_text(size=18)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend))

ggsave(".\\NBAR_results4\\trend_EVI&TCW_country.png", width = 10, height = 7.5, units = "in")

# by ecoregion
lc_rc1 = eco.sp.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-16,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = eco.sp.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-16,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                 data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 5), ] #remove class 5

lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Western Guinean\n Lowland Forests", "Eastern Guinean\n Forests", 
                               "Guinean Forest\n-Savanna Mosaic", "West Sudanian\n Savanna")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover") + 
  #theme(legend.position="none")+
  theme(legend.position=c(0.075,0.4))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.text.y = element_text(size=18)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend))

ggsave(".\\NBAR_results4\\trend_EVI&TCW_ecoregion.png", width = 10, height = 7.5, units = "in")

# by protected status
lc_rc1 = protected.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-10,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = protected.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-10,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                 data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))))

lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Reserve", "Eco-Reserve","Non-Protected")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover") + 
  #theme(legend.position="none")+
  theme(legend.position=c(0.075,0.4))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.text.y = element_text(size=18)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend))

ggsave(".\\NBAR_results4\\trend_EVI&TCW_protect.png", width = 10, height = 7.5, units = "in")

# by protected status * country
lc_rc1 = county_b_ghana.grd*10000 + protected.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-43,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/10000)
lc_rc1.df$protect = floor((lc_rc1.df$value - lc_rc1.df$lc*10000)/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*10000-lc_rc1.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = county_b_ghana.grd*10000 + protected.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-43,]
lc_rc2.df$lc = floor(lc_rc2.df$value/10000)
lc_rc2.df$protect = floor((lc_rc2.df$value - lc_rc2.df$lc*10000)/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*10000-lc_rc2.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = "TCW"), 
                 data.frame(lc_rc2.df, VI = "EVI"))

lc_rc1.df.long = melt(lc_rc.df[,c("lc","protect","trend","prop","VI")], id.vars=c("lc", "protect","trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("CÃ´te d'Ivoire", "Ghana", "Guinea", "Liberia", "Sierra Leone")
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  theme(legend.position=c(0.75,0.1))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=22))+ 
  theme(strip.text.y = element_text(size=22))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=1))

ggsave(".\\NBAR_results4\\trend_EVI&TCW_Protect&country.png", width = 7.5, height = 9, units = "in")

# by protected status * ecoregions
lc_rc1 = eco.sp.grd*10000 + protected.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-46,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/10000)
lc_rc1.df$protect = floor((lc_rc1.df$value - lc_rc1.df$lc*10000)/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*10000-lc_rc1.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = eco.sp.grd*10000 + protected.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-46,]
lc_rc2.df$lc = floor(lc_rc2.df$value/10000)
lc_rc2.df$protect = floor((lc_rc2.df$value - lc_rc2.df$lc*10000)/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*10000-lc_rc2.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = "TCW"), 
                 data.frame(lc_rc2.df, VI = "EVI"))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 5), ] #remove class 5

lc_rc1.df.long = melt(lc_rc.df[,c("lc","protect","trend","prop","VI")], id.vars=c("lc", "protect","trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Western Guinean\n Lowland Forests", "Eastern Guinean\n Forests", 
                               "Guinean Forest\n-Savanna Mosaic", "West Sudanian\n Savanna")
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  theme(legend.position=c(0.75,0.1))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=22))+ 
  theme(strip.text.y = element_text(size=22))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=1))

ggsave(".\\NBAR_results4\\trend_EVI&TCW_Protect&Ecoregion.png", width = 9, height = 9, units = "in")

























# 3.1 read into land cover and classify
setwd("R:\\users\\Zhihua\\MODIS")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

lc = raster("D:\\users\\Zhihua\\MODIS\\landcover\\LCType.tif")
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
clscolor = c("#4575b4", "#005a32","#01665e","#5ab4ac","#7fbf7b","#a6761d","#fc8d59")    		
breaks2 <- 0:7				
legendbrks2 <- 1:7 - 0.5


png(".\\NBAR_results4\\LC_type.png",height = 3000, width = 5000, res = 300, units = "px")

levelplot(lc_rc, main="Land Cover Type",
          maxpixels = nrow(lc_rc)*ncol(lc_rc),
          at= breaks2, margin=FALSE,
          col.regions= clscolor,
          colorkey= list(labels= list(labels= clasnames,at= legendbrks2, cex = 1.5)),space = "bottom") +
  layer(sp.polygons(county_b.geo, col = "black", lwd = 2)) +
  layer(sp.polygons(regions, col = "red", lwd = 2))

dev.off()

writeRaster(lc_rc,".\\NBAR_results4\\lc_rc.wa.tif", format="GTiff", overwrite=TRUE)

#############  plot trend for each biome ##############################

TCW.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/TCW2.dry.trend.wa.tif")
EVI.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/EVI2.dry.trend.wa.tif")
TCG.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/TCG2.dry.trend.wa.tif")
TCA.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/TCA2.dry.trend.wa.tif")
TCB.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/TCB2.dry.trend.wa.tif")
NDWI.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/NDWI2.dry.trend.wa.tif")

#county_b_ghana.r = rasterize(county_b_ghana, TCW.trd2.grd)
#county_b_ghana.r = county_b_ghana.r >= 1
#county_b_ghana.r[county_b_ghana.r != 1] = NA
#county_b_ghana.r[lc_rc == 1] = NA

#writeRaster(county_b_ghana.r,".\\NBAR_results4\\studyarea.mask.tif", format="GTiff", overwrite=TRUE)
county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")
lc_rc = raster(".\\NBAR_results4\\lc_rc.wa.tif")

TCW.trd2.grd = TCW.trd2.grd*county_b_ghana.r
EVI.trd2.grd = EVI.trd2.grd*county_b_ghana.r
TCG.trd2.grd = TCG.trd2.grd*county_b_ghana.r
TCA.trd2.grd = TCA.trd2.grd*county_b_ghana.r
TCB.trd2.grd = TCB.trd2.grd*county_b_ghana.r
NDWI.trd2.grd = NDWI.trd2.grd*county_b_ghana.r

lc_rc1 = lc_rc*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-25,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 4)

lc_rc2 = lc_rc*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-25,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 4)

lc_rc3 = lc_rc*100 + NDWI.trd2.grd 
lc_rc3.df = data.frame(freq(lc_rc3))[-25,]
lc_rc3.df$lc = floor(lc_rc3.df$value/100)
lc_rc3.df$trend = lc_rc3.df$value - lc_rc3.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc3.df, FUN = "sum")
lc_rc3.df$prop = 100*lc_rc3.df$count/rep(lc_sum$count, each = 4)

lc_rc4 = lc_rc*100 + TCG.trd2.grd 
lc_rc4.df = data.frame(freq(lc_rc4))[-25,]
lc_rc4.df$lc = floor(lc_rc4.df$value/100)
lc_rc4.df$trend = lc_rc4.df$value - lc_rc4.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc4.df, FUN = "sum")
lc_rc4.df$prop = 100*lc_rc4.df$count/rep(lc_sum$count, each = 4)

lc_rc5 = lc_rc*100 + TCA.trd2.grd 
lc_rc5.df = data.frame(freq(lc_rc5))[-25,]
lc_rc5.df$lc = floor(lc_rc5.df$value/100)
lc_rc5.df$trend = lc_rc5.df$value - lc_rc5.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc5.df, FUN = "sum")
lc_rc5.df$prop = 100*lc_rc5.df$count/rep(lc_sum$count, each = 4)


lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                 data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))),
                 data.frame(lc_rc3.df, VI = rep("NDWI", nrow(lc_rc2.df))),
                 data.frame(lc_rc4.df, VI = rep("TCG", nrow(lc_rc2.df))),
                 data.frame(lc_rc5.df, VI = rep("TCA", nrow(lc_rc2.df))))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 1 | lc_rc.df$lc == 7| lc_rc.df$lc == 5), ] #remove class 1: water; and 7: other class
lc_rc.df = lc_rc.df[-which(lc_rc.df$trend == 4), ] #not calculated classes

library(reshape2)
library(ggplot2)

lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- clasnames[c(2,3,4,6)]


ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover") + 
  theme(legend.position="none")+
theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.text.y = element_text(size=18))


ggsave(".\\NBAR_results4\\trend_biomes.png", width = 10, height = 7.5, units = "in")

#only plot TCW and EVI
ggplot(data=lc_rc1.df.long[which(lc_rc1.df.long$VI == "TCW"|lc_rc1.df.long$VI == "EVI"),], 
       aes(x=trend, y=value, fill=trend)) +
       facet_grid(VI ~ lc) +
       geom_bar(stat="identity", position=position_dodge(), colour="black") +
       xlab("") + ylab("Percentage of Land Cover") + 
  theme(legend.position="none")+
  theme(axis.title.x = element_text(colour="black", size=14),axis.text.x  = element_text(colour="black",size=12))+
  theme(axis.title.y = element_text(colour="black", size=14),axis.text.y  = element_text(colour="black",size=12))+
  theme(strip.text.x = element_text(size=14))+
  theme(strip.text.y = element_text(size=14))


ggsave(".\\NBAR_results4\\trend_biomes_EVI&TCW.png", width = 10, height = 7.5, units = "in")


#PLOT combination of EVI and TCW, TCW2: 1:negative, 2:positive, 3:no trend  		   
trd = TCW.trd2.grd
trd[] = NA
trd[TCW.trd2.grd == 2 & EVI.trd2.grd ==1] = 1
trd[TCW.trd2.grd == 2 & EVI.trd2.grd ==3] = 2
trd[TCW.trd2.grd == 2 & EVI.trd2.grd ==2] = 3
trd[TCW.trd2.grd == 3 & EVI.trd2.grd ==2] = 4
trd[TCW.trd2.grd == 1 & EVI.trd2.grd ==2] = 5
trd[TCW.trd2.grd == 1 & EVI.trd2.grd ==3] = 6
trd[TCW.trd2.grd == 1 & EVI.trd2.grd ==1] = 7
trd[TCW.trd2.grd == 3 & EVI.trd2.grd ==1] = 8
trd[is.na(trd)&TCW.trd2.grd < 4] = 9

#writeRaster(trd,".\\NBAR_results4\\Trend_com_data.tif", format="GTiff", overwrite=TRUE)

trd = raster(".\\NBAR_results4\\Trend_com_data.tif")
trd.arg = c("HLD", "WDI","WHLG","HDI","HLG", "WDD","WHLD","HDD","NoTrend")

lc_rc1 = lc_rc*100 + trd
lc_rc1.df = data.frame(freq(lc_rc1))[-55,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 9)

lc_rc.df = lc_rc1.df
lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 1 | lc_rc.df$lc == 7| lc_rc.df$lc == 5), ] #remove class 1: water; and 7: other class


library(reshape2)
library(ggplot2)
require("RColorBrewer")

color1 <- c(brewer.pal(9,"Reds")[c(3,6,9)], brewer.pal(9,"Greens")[c(3,6,9)],brewer.pal(9,"Blues")[c(3,6)],"#cccccc")

lc_rc.df.long = melt(lc_rc.df[,c("lc","trend","prop")], id.vars=c("lc", "trend"))

lc_rc.df.long$trend = factor(lc_rc.df.long$trend)
lc_rc.df.long$lc = factor(lc_rc.df.long$lc)
#levels(lc_rc.df.long$lc) <- clasnames[c(2,3,4,6)]
levels(lc_rc.df.long$lc) <- c("Forest", "Woody Savannas", "Savannas", "Crop/Natural Mosaic")
levels(lc_rc.df.long$trend) <- trd.arg

ggplot(data=lc_rc.df.long, aes(x=trend, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  facet_wrap( ~ lc,ncol=2) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.75,0.3))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=22))+ 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc.df.long$trend),
                    labels=levels(lc_rc.df.long$trend))

ggsave(".\\NBAR_results4\\Trend_com_stat.png", width = 9, height = 6, units = "in")

# plot
require("RColorBrewer")

breaks2_change <- 0:9        
legendbrks2_change <- 1:9 - 0.5

arg1 <- list(at=seq(1,9,1), labels=trd.arg) #these are the class names
labels = trd.arg

trd2 = trd
trd2[lc_rc == 1|lc_rc == 5|lc_rc == 7] = NA

png(".\\NBAR_results4\\Trend_com.png",height = 3000, width = 5000, res = 300, units = "px")

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(trd2,
          maxpixels = nrow(trd)*ncol(trd),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) + 
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5)) 

dev.off()
