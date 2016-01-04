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


