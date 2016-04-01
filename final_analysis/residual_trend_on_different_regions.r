########## section 1: calculate 3 by 3 moving windos trends ################
# define a function first
Ratio.MW = function(inraster, size = 3){
  
  inraster.neg = inraster == 1
  inraster.pos = inraster == 2
  inraster.neg.3 <- focal(inraster.neg, w=matrix(1, nc=size, nr=size), fun=sum)
  inraster.pos.3 <- focal(inraster.pos, w=matrix(1, nc=size, nr=size), fun=sum)
  
  ratio = inraster
  ratio[] = 0
  ratio[inraster.neg.3 > inraster.pos.3 ] = 1 #Negative Dominant
  ratio[inraster.pos.3 > inraster.neg.3 ] = 2 #Positive Dominant
  ratio[inraster.neg.3 == inraster.pos.3] = 3   #no/equal change
  ratio[ratio == 0] = NA
  
  return(ratio)
  
}

setwd("R:\\users\\Zhihua\\MODIS")
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)

# read into shapefiles
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))
county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")

# read into residual tif
TCW.trd2.grd.res = raster(".\\NBAR_results4\\tcw.df.residual.trd.wa2.tif")
EVI.trd2.grd.res = raster(".\\NBAR_results4\\evi.df.residual.trd.wa2.tif")

#calculate ratio
TCW.trd2.grd.res.ratio = Ratio.MW(TCW.trd2.grd.res, size = 3)
EVI.trd2.grd.res.ratio = Ratio.MW(EVI.trd2.grd.res, size = 3)

writeRaster(TCW.trd2.grd.res.ratio,".\\NBAR_results4\\Residual.TCW2.dry.trend.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(EVI.trd2.grd.res.ratio,".\\NBAR_results4\\Residual.EVI2.dry.trend.wa2.tif", format="GTiff", overwrite=TRUE)

########## section 2: plot trend ################
ratio = stack(TCW.trd2.grd.res.ratio, EVI.trd2.grd.res.ratio)
ratio = ratio*county_b_ghana.r
names(ratio) <- c("TCW","EVI")

breaks2_change <- 0:3        
legendbrks2_change <- 1:3 - 0.5

labels = c("Negative Dominant","Positive Dominant","No/Equal Change")
labels = c("-","+","NT")

arg1 <- list(at=seq(1,3,1), labels=labels) #these are the class names
color1=c("#cccccc", "#6baed6","#fb6a4a","#74c476")
color1=c("#fb6a4a","#74c476","#cccccc")


png(file = ".\\NBAR_results4\\residual.trend.wa.TCW&EVI_cliped-window3.png", width = 6000, height = 5000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(ratio,
          maxpixels = nrow(ratio)*ncol(ratio),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) 
  # +  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5)) 

dev.off()


########## section 3: calculate trend on different regions ################
library("raster")
library("dismo")
library("sp")
library("rgdal")
library("reshape2")
library("ggplot2")


setwd("R:/users/Zhihua/MODIS")

TCW.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/Residual.TCW2.dry.trend.wa2.tif")
EVI.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/Residual.EVI2.dry.trend.wa2.tif")
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

#protected area 
#1: protected, 2: strictly protect (Eco-protected) 3:non-protected
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

#statistics
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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_country.png", width = 10, height = 7.5, units = "in")

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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion.png", width = 10, height = 7.5, units = "in")

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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_protect.png", width = 10, height = 7.5, units = "in")

#by land cover
lc_rc1 = lc_rc*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-22,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = lc_rc*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-22,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                 data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 1 | lc_rc.df$lc == 7| lc_rc.df$lc == 5), ] #remove class 1: water; and 7: other class
#lc_rc.df = lc_rc.df[-which(lc_rc.df$trend == 4), ] #not calculated classes
lc_rc1.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- clasnames[c(2,3,4,6)]
levels(lc_rc1.df.long$lc) <- c("Tropical Forest", "Woody Savannas","Savannas","Crop/Natural\n Mosaic")

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

ggsave(".\\NBAR_results4\\residual.trend_biomes_EVI&TCW-2.png", width = 10, height = 7.5, units = "in")


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
levels(lc_rc1.df.long$lc) <- c("CIV", "GHA", "GIN", "LBR", "SLE")
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.6,0.45))+
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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_Protect&country.png", width = 12, height = 9, units = "in")

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

levels(lc_rc1.df.long$lc) <- c("WGLF", "EGF", "GFSM", "WSS")                               
                               
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=1))

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_Protect&Ecoregion.png", width = 12, height = 9, units = "in")


# by protected status * land cover
lc_rc1 = lc_rc*10000 + protected.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-55,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/10000)
lc_rc1.df$protect = floor((lc_rc1.df$value - lc_rc1.df$lc*10000)/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*10000-lc_rc1.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

lc_rc2 = lc_rc*10000 + protected.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-55,]
lc_rc2.df$lc = floor(lc_rc2.df$value/10000)
lc_rc2.df$protect = floor((lc_rc2.df$value - lc_rc2.df$lc*10000)/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*10000-lc_rc2.df$protect*100
lc_sum = aggregate(count~protect+lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = "TCW"), 
                 data.frame(lc_rc2.df, VI = "EVI"))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 1 | lc_rc.df$lc == 7| lc_rc.df$lc == 5), ] #remove class 1: water; and 7: other class
lc_rc1.df.long = melt(lc_rc.df[,c("lc","protect","trend","prop","VI")], id.vars=c("lc", "protect","trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Tropical Forest", "Woody Savannas","Savannas","Crop/Natural\n Mosaic")
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP")

color1 = c("#fb6a4a", "#67a9cf", "#cccccc")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(VI ~ lc) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  theme(legend.position=c(0.75,0.45))+
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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_Protect&landcover.png", width = 12, height = 9, units = "in")



# by protected status * ecoregions, and adding mean for all the region in 1/21/2016
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

#add change in all the regions
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

lc_rc2.df = lc_rc.df[-which(lc_rc.df$lc == 5), ] #remove class 5
lc_rc2.df.long = melt(lc_rc2.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc2.df.long = data.frame(cbind(lc = lc_rc2.df.long[,1], protect = 4, lc_rc2.df.long[,c(2:5)]))

#combine them together
lc_rc1.df.long = rbind(lc_rc1.df.long, lc_rc2.df.long)

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
#change the order of the trend
lc_rc1.df.long$trend <- ordered(lc_rc1.df.long$trend, levels = c("Negative", "No Trend", "Positive"))

lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Western Guinean\n Lowland Forests", "Eastern Guinean\n Forests", 
                               "Guinean Forest\n-Savanna Mosaic", "West Sudanian\n Savanna")
levels(lc_rc1.df.long$lc) <- c("WGLF", "EGF", "GFSM", "WSS")                               
#change the order of the ecoregion
lc_rc1.df.long$lc <- ordered(lc_rc1.df.long$lc, levels = c("WSS", "GFSM", "EGF","WGLF"))
                           
lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected", "All")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP", "ALL")

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.position="top")+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=3))

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion2.png", width = 10, height = 7.5, units = "in")

# by protected status * countries, and adding mean for all the countries in 1/21/2016
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

#add change in all the regions, without consideration of projected status
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

lc_rc2.df.long = melt(lc_rc.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))
lc_rc2.df.long = data.frame(cbind(lc = lc_rc2.df.long[,1], protect = 4, lc_rc2.df.long[,c(2:5)]))

#combine them together
lc_rc1.df.long = rbind(lc_rc1.df.long, lc_rc2.df.long)

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
#change the order of the trend
lc_rc1.df.long$trend <- ordered(lc_rc1.df.long$trend, levels = c("Negative", "No Trend", "Positive"))

lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("CÃ´te d'Ivoire", "Ghana", "Guinea", "Liberia", "Sierra Leone")
levels(lc_rc1.df.long$lc) <- c("CIV", "GHA", "GIN", "LBR", "SLE")

lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected", "All")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP", "ALL")

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")

ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.position="top")+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=3))

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_country2.png", width = 10, height = 7.5, units = "in")

#plot the comparsion between Ghana and Ivory coast
# for TCW
lc_rc1 = county_b_ghana.grd*1000 + protected.grd*100+eco.sp.grd*10 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-138,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/1000)
lc_rc1.df$protect = floor((lc_rc1.df$value - lc_rc1.df$lc*1000)/100)
lc_rc1.df$ecoregion = floor((lc_rc1.df$value - lc_rc1.df$lc*1000-lc_rc1.df$protect*100)/10)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*1000-lc_rc1.df$protect*100-lc_rc1.df$ecoregion*10
#only keep cote devior [1] and ghana [2], and {2-4}
lc_rc1.df = lc_rc1.df[which(lc_rc1.df$lc < 3 & lc_rc1.df$ecoregion > 1 & lc_rc1.df$ecoregion < 5), ]
lc_sum = aggregate(count~ecoregion+lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

#add change in all the regions, remove the protected status

lc_rc2 = county_b_ghana.grd*1000 + eco.sp.grd*10 + TCW.trd2.grd
lc_rc2.df = data.frame(freq(lc_rc2))[-57,]
lc_rc2.df$lc = floor(lc_rc2.df$value/1000)
lc_rc2.df$ecoregion = floor((lc_rc2.df$value - lc_rc2.df$lc*1000)/10)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*1000 - lc_rc2.df$ecoregion*10
#only keep cote devior [1] and ghana [2], and {2-4}
lc_rc2.df = lc_rc2.df[which(lc_rc2.df$lc < 3 & lc_rc2.df$ecoregion > 1 & lc_rc2.df$ecoregion < 5), ]
lc_sum = aggregate(count~ecoregion+lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc2.df = data.frame(lc_rc2.df[,c("value", "count","lc")], protect = 4, lc_rc2.df[,c("ecoregion", "trend","prop")])

lc_rc.df = rbind(lc_rc1.df, lc_rc2.df)

lc_rc1.df.long = melt(lc_rc.df[,c("lc","protect","ecoregion","trend","prop")], id.vars=c("lc", "protect","ecoregion","trend"))


lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")

#change the order of the trend
lc_rc1.df.long$trend <- ordered(lc_rc1.df.long$trend, levels = c("Negative", "No Trend", "Positive"))

lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("CIV", "GHA")                               

lc_rc1.df.long$ecoregion = factor(lc_rc1.df.long$ecoregion)
levels(lc_rc1.df.long$ecoregion) <- c("EGF", "GFSM", "WSS")                            
                             
#change the order of the ecoregion
lc_rc1.df.long$ecoregion <- ordered(lc_rc1.df.long$ecoregion, levels = c("WSS", "GFSM", "EGF"))

lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected", "All")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP", "ALL")

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")


ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(ecoregion~lc) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.position="top")+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=3))

ggsave(".\\NBAR_results4\\residual.trend_ghana_cold_tcw.png", width = 10, height = 7.5, units = "in")

# for EVI
lc_rc1 = county_b_ghana.grd*1000 + protected.grd*100+eco.sp.grd*10 + EVI.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-138,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/1000)
lc_rc1.df$protect = floor((lc_rc1.df$value - lc_rc1.df$lc*1000)/100)
lc_rc1.df$ecoregion = floor((lc_rc1.df$value - lc_rc1.df$lc*1000-lc_rc1.df$protect*100)/10)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*1000-lc_rc1.df$protect*100-lc_rc1.df$ecoregion*10
#only keep cote devior [1] and ghana [2], and {2-4}
lc_rc1.df = lc_rc1.df[which(lc_rc1.df$lc < 3 & lc_rc1.df$ecoregion > 1 & lc_rc1.df$ecoregion < 5), ]
lc_sum = aggregate(count~ecoregion+lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 3)

#add change in all the regions, remove the protected status

lc_rc2 = county_b_ghana.grd*1000 + eco.sp.grd*10 + EVI.trd2.grd
lc_rc2.df = data.frame(freq(lc_rc2))[-57,]
lc_rc2.df$lc = floor(lc_rc2.df$value/1000)
lc_rc2.df$ecoregion = floor((lc_rc2.df$value - lc_rc2.df$lc*1000)/10)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*1000 - lc_rc2.df$ecoregion*10
#only keep cote devior [1] and ghana [2], and {2-4}
lc_rc2.df = lc_rc2.df[which(lc_rc2.df$lc < 3 & lc_rc2.df$ecoregion > 1 & lc_rc2.df$ecoregion < 5), ]
lc_sum = aggregate(count~ecoregion+lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 3)

lc_rc2.df = data.frame(lc_rc2.df[,c("value", "count","lc")], protect = 4, lc_rc2.df[,c("ecoregion", "trend","prop")])

lc_rc.df = rbind(lc_rc1.df, lc_rc2.df)

lc_rc1.df.long = melt(lc_rc.df[,c("lc","protect","ecoregion","trend","prop")], id.vars=c("lc", "protect","ecoregion","trend"))


lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")

#change the order of the trend
lc_rc1.df.long$trend <- ordered(lc_rc1.df.long$trend, levels = c("Negative", "No Trend", "Positive"))

lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("CIV", "GHA")                               

lc_rc1.df.long$ecoregion = factor(lc_rc1.df.long$ecoregion)
levels(lc_rc1.df.long$ecoregion) <- c("EGF", "GFSM", "WSS")                            
                             
#change the order of the ecoregion
lc_rc1.df.long$ecoregion <- ordered(lc_rc1.df.long$ecoregion, levels = c("WSS", "GFSM", "EGF"))

lc_rc1.df.long$protect = factor(lc_rc1.df.long$protect)
levels(lc_rc1.df.long$protect) <- c("Reserve", "Eco-Reserve","Non-Protected", "All")
levels(lc_rc1.df.long$protect) <- c("R", "ER","NP", "ALL")

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")


ggplot(data=lc_rc1.df.long, aes(x=protect, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(ecoregion~lc) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.position="top")+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=3))

ggsave(".\\NBAR_results4\\residual.trend_ghana_cold_evi.png", width = 10, height = 7.5, units = "in")


# by country * ecoregions, and adding mean for all the region in 3/31/2016
lc_rc1 = eco.sp.grd*10000 + county_b_ghana.grd*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-57,] #
lc_rc1.df$lc = floor(lc_rc1.df$value/10000)
lc_rc1.df$country = floor((lc_rc1.df$value - lc_rc1.df$lc*10000)/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*10000-lc_rc1.df$country*100
lc_sum = data.frame(aggregate(count~country+lc, data = lc_rc1.df, FUN = "sum"))
test = merge.data.frame(lc_rc1.df, lc_sum, by = c("lc", "country"))
test$prop = 100*test$count.x/test$count.y
lc_rc1.df = test

lc_rc2 = eco.sp.grd*10000 + county_b_ghana.grd*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-57,]
lc_rc2.df$lc = floor(lc_rc2.df$value/10000)
lc_rc2.df$country = floor((lc_rc2.df$value - lc_rc2.df$lc*10000)/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*10000-lc_rc2.df$country*100
lc_sum = data.frame(aggregate(count~country+lc, data = lc_rc2.df, FUN = "sum"))
test = merge.data.frame(lc_rc2.df, lc_sum, by = c("lc", "country"))
test$prop = 100*test$count.x/test$count.y
lc_rc2.df = test

lc_rc.df = rbind(data.frame(lc_rc1.df, VI = "TCW"), 
                 data.frame(lc_rc2.df, VI = "EVI"))

lc_rc.df = lc_rc.df[-which(lc_rc.df$lc == 5), ] #remove class 5

lc_rc1.df.long = melt(lc_rc.df[,c("lc","country","trend","prop","VI")], id.vars=c("lc", "country","trend","VI"))

#
lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Negative", "Positive", "No Trend")
#change the order of the trend
lc_rc1.df.long$trend <- ordered(lc_rc1.df.long$trend, levels = c("Negative", "No Trend", "Positive"))

lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- c("Western Guinean\n Lowland Forests", "Eastern Guinean\n Forests", 
                               "Guinean Forest\n-Savanna Mosaic", "West Sudanian\n Savanna")
levels(lc_rc1.df.long$lc) <- c("WGLF", "EGF", "GFSM", "WSS")                               
#change the order of the ecoregion
lc_rc1.df.long$lc <- ordered(lc_rc1.df.long$lc, levels = c("WSS", "GFSM", "EGF","WGLF"))

lc_rc1.df.long$country = factor(lc_rc1.df.long$country)
levels(lc_rc1.df.long$country) <- c("CÃ´te d'Ivoire", "Ghana", "Guinea", "Liberia", "Sierra Leone")
levels(lc_rc1.df.long$country) <- c("CIV", "GHA", "GIN", "LBR", "SLE")

#re-order country based annual mean rainfall
#lc_rc1.df.long$country <- ordered(lc_rc1.df.long$country, levels = c("GIN", "SLE", "LBR","CIV","GHA"))
lc_rc1.df.long$country <- ordered(lc_rc1.df.long$country, levels = c("SLE", "LBR","GIN", "CIV","GHA"))

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")

#remove conutry == LBR & lc == GFSM
lc_rc1.df.long = lc_rc1.df.long[-which(lc_rc1.df.long$country == "LBR" & lc_rc1.df.long$lc == "GFSM"), ]


ggplot(data=lc_rc1.df.long, aes(x=country, y=value, fill=trend)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  facet_grid(lc ~ VI) +
  xlab("") + ylab("Percentage of Land Cover") +
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_text(face="bold", colour="black", size=22),axis.text.x  = element_text(colour="black",size=20))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=22),axis.text.y  = element_text(colour="black",size=20))+
  #theme(legend.position=c(0.7,0.45))+
  theme(legend.position="top")+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(strip.text.x = element_text(size=20))+ 
  theme(strip.text.y = element_text(size=20))+ 
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(lc_rc1.df.long$trend),
                    labels=levels(lc_rc1.df.long$trend)) +
  guides(fill=guide_legend(ncol=3))

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion-country.png", width = 10, height = 7.5, units = "in")

#plot some highlighted area
regions = readOGR(dsn="R:/users/Zhihua/MODIS/NBAR_results3", layer="Val_regions")
projection(regions) <- proj.geo 
regions.sin <- spTransform(regions, CRS(proj.sin)) 

TCW.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/Residual.TCW2.dry.trend.wa2.tif")
EVI.trd2.grd = raster("R:/users/Zhihua/MODIS/NBAR_results4/Residual.EVI2.dry.trend.wa2.tif")
#read into mat hansen's landsat forest loss and reclassify
landsat_netloss = raster("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid.tif")

landsat_netloss_rcls = recls2(landsat_netloss,threshold = c(0, 0.05, 0.1, 0.2, 0.5))
breaks2_change_hansen <- 0:6        
legendbrks2_change_hansen <- 1:6 - 0.5
color1_hansen=c('#66c2a4','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')
labels_hansen = c("No Loss","0-0.05","0.05-0.1","0.1-0.2","0.2-0.5",">0.5")

#trend color scheme
breaks2_change_trd <- 0:4        
legendbrks2_change_trd <- 1:4 - 0.5
arg1_trd <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels_trd = c("Negative","Positive","No Trend","Not Calculated")
color1_trd = c("#e66101", "#1a9641","#ffffbf","#2b83ba")
color1_trd=c("#fb6a4a","#67a9cf","#cccccc","#74c476")

regionid = c(2,4,7)
# plot in one figure


png(file = paste("NBAR_results4/residual.Highlighted_region", "combined.png", sep = ""), width = 3000, height = 3000, units = "px", res = 300)

par(mfrow=c(3,3),mar=c(0.02, 0.02, 0.02, 0.02))

for (i in regionid){
  TCW.trd2.grd1 = crop(TCW.trd2.grd, regions[i,])
  EVI.trd2.grd1 = crop(EVI.trd2.grd, regions[i,])
  landsat_netloss_rcls1 = crop(landsat_netloss_rcls, regions[i,])

  #tcw
  plot(TCW.trd2.grd1,col = color1_trd[c(1:3)], 
       legend=FALSE,
       axes=FALSE,
       box=FALSE,
  )
  plot(county_b_ghana, add = TRUE)
  #evi
  plot(EVI.trd2.grd1,col = color1_trd[c(1:3)], 
       legend=FALSE,
       axes=FALSE,
       box=FALSE,
  )
  plot(county_b_ghana, add = TRUE)
 
#landsat
  
  plot(landsat_netloss_rcls1,col = color1_hansen[sort(unique(landsat_netloss_rcls1))], 
       legend=FALSE,
       axes=FALSE,
       box=FALSE,
  )
  plot(county_b_ghana, add = TRUE)
}

dev.off()

tcw.trd.df = data.frame(freq(TCW.trd2.grd))[-4,]; tcw.trd.df$prop = tcw.trd.df$count/sum(tcw.trd.df$count)
evi.trd.df = data.frame(freq(EVI.trd2.grd))[-4,]; evi.trd.df$prop = evi.trd.df$count/sum(evi.trd.df$count)
