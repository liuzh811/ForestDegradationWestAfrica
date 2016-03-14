#3/13/2016
#plot annual indices for EVI and TCW only

setwd("D:\\users\\Zhihua\\MODIS")
library(rgdal)
library(raster)
library(rasterVis)

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
proj.sin = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))
county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#read into mask files
county_b_ghana.r = raster(".\\NBAR_results4\\studyarea.mask.tif")

# read into annual indices
TCW.trd2.grd = stack(".\\NBAR_results4\\tcw2.dry.wa.grd")
EVI.trd2.grd = stack(".\\NBAR_results4\\evi2.dry.wa.grd")

#clip to shape
TCW.trd2.grd = TCW.trd2.grd*county_b_ghana.r
EVI.trd2.grd = EVI.trd2.grd*county_b_ghana.r

names(TCW.trd2.grd) <- paste("Y", 2001:2015, sep = "")
names(EVI.trd2.grd) <- paste("Y", 2001:2015, sep = "")

#plot TCW
library(colorRamps)

color_tc3 = rev(rainbow(99, start=0,end=1))
color_tc32 = rev(blue2green2red(99))

qu.val = quantile(as.vector(as.matrix(TCW.trd2.grd)), prob = c(0.01, 0.99), na.rm = T)

breaks_tc3 <- round(seq(min(minValue(TCW.trd2.grd)),max(maxValue(TCW.trd2.grd)),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(minValue(TCW.trd2.grd)),max(maxValue(TCW.trd2.grd)),length.out=10),0)

breaks_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=10),0)

png(".\\NBAR_results4\\annual_tcw.png", height = 5000, width = 5000, res = 300, units = "px")

levelplot(TCW.trd2.grd, main="Annual TC Wetness Map",
          at= breaks_tc3, margin=FALSE,
          maxpixels = nrow(TCW.trd2.grd)*ncol(TCW.trd2.grd),
          col.regions=color_tc32,
          colorkey= list(labels= list(labels= legendbrks2_tc3,at= legendbrks2_tc3, cex = 1.5), space = "bottom"),
          layout=c(3, 5))+
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 2))

dev.off()


#plot EVI
library(colorRamps)

EVI.trd2.grd = EVI.trd2.grd*1000

color_tc3 = rev(rainbow(99, start=0,end=1))
color_tc32 = rev(blue2green2red(99))

qu.val = quantile(as.vector(as.matrix(EVI.trd2.grd)), prob = c(0.01, 0.99), na.rm = T)

breaks_tc3 <- round(seq(min(minValue(EVI.trd2.grd)),max(maxValue(EVI.trd2.grd)),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(minValue(EVI.trd2.grd)),max(maxValue(EVI.trd2.grd)),length.out=10),0)

breaks_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=10),0)

png(".\\NBAR_results4\\annual_evi.png", height = 5000, width = 5000, res = 300, units = "px")

levelplot(EVI.trd2.grd, main="Annual EVI Map",
          at= breaks_tc3, margin=FALSE,
          maxpixels = nrow(EVI.trd2.grd)*ncol(EVI.trd2.grd),
          col.regions=color_tc32,
          colorkey= list(labels= list(labels= legendbrks2_tc3,at= legendbrks2_tc3, cex = 1.5), space = "bottom"),
          layout=c(3, 5))+
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 2))

dev.off()

#calculate the number of acculumative NAs at the begining and the end of the time series
# define a function first
ZeroAcc = function(x){ #x is a vector with 15 values
  
  if (( is.na(x[1]) & is.na(x[2]) & is.na(x[3]) & is.na(x[4]) & is.na(x[5]) ) |
        ( is.na(x[11]) & is.na(x[12]) & is.na(x[13]) & is.na(x[14]) & is.na(x[15]) )) 
    
  y = 1
  else (y = 0)
  
  return(y)
  
}

t1 <- calc(EVI.trd2.grd, ZeroAcc)

t1 = t1*county_b_ghana.r
t1 = t1+1

#plot
clscolor_change = c("#4d9221","#c51b7d")

breaks2_change <- 0:2        
legendbrks2_change <- 1:2 - 0.5
clasnames_change = c("No", "Yes")

png(".\\NBAR_results4\\Consecutive5nas.png", height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1.3, lines=2)
levelplot(t1,
          maxpixels = nrow(t1)*ncol(t1),
          at= breaks2_change, margin=FALSE,
          col.regions= clscolor_change,
          colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5)) 

dev.off()

#precent of yes: 1.4%, most of these are in the coastal area, the trend is not calculated at all
freq(t1)
