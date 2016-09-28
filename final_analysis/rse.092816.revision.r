
	 
# revise on 9/27/2016
# make figures and calculate propotations based on residual trend, before the 3 by 3 moving box
Figure 4 (trend map) is in C:\zhihua\dataset\rs_revision

########## section 3: calculate trend on different regions ################
library("raster")
library("dismo")
library("sp")
library("rgdal")
library("reshape2")
library("ggplot2")


setwd("F:/Rosa/MODIS")

TCW.trd2.grd = raster(".\\NBAR_results4\\tcw.df.residual.trd.wa2.tif")
EVI.trd2.grd = raster(".\\NBAR_results4\\evi.df.residual.trd.wa2.tif")

TCW.trd2.grd[TCW.trd2.grd == 4] = 3
EVI.trd2.grd[EVI.trd2.grd == 4] = 3
#read into land cover, and remove changed area
lc = raster("C:/zhihua/dataset/rs_revision/data/LCTypeWAF2.tif")
TCW.trd2.grd[lc == 0| lc == 11 | lc ==13 | lc == 16] = NA
EVI.trd2.grd[lc == 0| lc == 11 | lc ==13 | lc == 16] = NA

lc_rc = raster(".\\NBAR_results4\\lc_rc.wa.tif")
county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")

#read into shapefiles, countries
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
county_b = readOGR(dsn="F:/Rosa\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

county_b_ghana.grd = rasterize(county_b_ghana, county_b_ghana.r)
county_b_ghana.grd = county_b_ghana.grd*county_b_ghana.r #1:"CÃ´te d'Ivoire", 2:"Ghana", 3:"Guinea", 4:Liberia, 5:Sierra Leone

#protected area 
#1: protected, 2: strictly protect (Eco-protected) 3:non-protected
protected = readOGR(dsn="F:/Rosa\\GeoData_West Africa\\protectarea",layer="protected_areas_west_africa2")

# protected@data$status = 1 #
# protected@data$status[which(protected@data$IUCN_CAT == "Ia" | 
#                               protected@data$IUCN_CAT == "II" | protected@data$IUCN_CAT == "IV")] = 2

# protected.grd = rasterize(protected, county_b_ghana.r, field = "status")
protected.grd = rasterize(protected, county_b_ghana.r, field = "status2")
protected.grd = protected.grd*county_b_ghana.r #protected
protected.grd[is.na(protected.grd)] = 3  #3 non-protected
protected.grd = protected.grd*county_b_ghana.r 

#read into ecoregion; ECO_NAME2: 1: west guinea forest; 2: east guinea foreste; 3: mosaic, 4: savvannas,5:other
eco.sp = readOGR(dsn="F:/Rosa\\WWF_ecoregion\\official_teow\\official",layer="wwf_terr_ecos1")
eco.sp.grd = rasterize(eco.sp, county_b_ghana.r, field = eco.sp@data$ECO_NAME2)
eco.sp.grd = eco.sp.grd*county_b_ghana.r

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

write.csv(lc_rc1.df.long, ".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion2-092816.csv")

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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion2-092816.png", width = 10, height = 7.5, units = "in")

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
lc_rc1.df.long$country <- ordered(lc_rc1.df.long$country, levels = c("GHA", "CIV","GIN", "LBR","SLE"))

#color1 = c("#fb6a4a", "#67a9cf", "#cccccc"): color order for c("Negative", "Positive", "No Trend")
color1 = c("#fb6a4a", "#cccccc", "#67a9cf") # color order for c("Negative", "No Trend","Positive")

#remove conutry == LBR & lc == GFSM
lc_rc1.df.long = lc_rc1.df.long[-which(lc_rc1.df.long$country == "LBR" & lc_rc1.df.long$lc == "GFSM"), ]

write.csv(lc_rc1.df.long, ".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion-country-092816.csv")

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

ggsave(".\\NBAR_results4\\residual.trend_EVI&TCW_ecoregion-country-092816.png", width = 10, height = 7.5, units = "in")



