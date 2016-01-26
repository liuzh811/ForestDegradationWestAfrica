setwd("D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2")

library(raster)
library(rgdal)

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#read into data
fn = list.files(path = ".", pattern = "*.nc$")
r = stack(fn)
r = crop(r, county_b_ghana) #crop the dataset to western Africa boundary
names(r) <- fn

R_year = calc(r, mean)  

plot(R_year*365)
plot(county_b_ghana, add = TRUE)

#plot
recls2 = function(x, threshold = threshold){
  #reclassify
  x.list = list()
  for(layer in 1:nlayers(x))
  {
    x2 = x[[layer]]
    x2[!is.na(x[[layer]])] = 1
    for(i in 1:length(threshold)){
      x2[x[[layer]]>threshold[i]] = i+1 
    }
    x2[is.na(x[[layer]])] = NA
    x.list[[layer]] <- x2
  }
  x.list = stack(x.list)
  return(x.list)
}

R_year = R_year*365
threshold = c(1000,1250,1500,1750,2000)
R_year2 = recls2(R_year, threshold=threshold)

county_b_ghana.grd = rasterize(county_b_ghana, R_year,field = 1)
R_year2 = R_year2*county_b_ghana.grd


#plot

library("rasterVis")

breaks2_change <- 0:6        
legendbrks2_change <- 1:6 - 0.5

arg1 <- list(at=seq(1,6,1), labels=c("<1000","1000-1250","1250-1500","1500-1750","1750-2000",">2000")) #these are the class names
labels = c("<1000","1000-1250","1250-1500","1500-1750","1750-2000",">2000")
color1=c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac')

png(file = "D:\\users\\Zhihua\\MODIS\\NBAR_results4\\TRMM.map.png", width = 3000, height = 2000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(R_year2,
          maxpixels = nrow(R_year2)*ncol(R_year2),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) 

dev.off()
