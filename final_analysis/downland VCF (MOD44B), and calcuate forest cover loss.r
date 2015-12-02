# section 5: downland VCF (MOD44B), and calcuate forest cover loss
setwd("D:\\users\\Zhihua\\MODIS")
library(XML)
library(RCurl)

#construct date first
mod.date = c("2000.03.05", "2001.03.06","2002.03.06","2003.03.06","2004.03.05","2005.03.06","2006.03.06","2007.03.06",
             "2008.03.05","2009.03.06","2010.03.06","2011.03.06","2012.03.05","2013.03.06","2014.03.06")
tileh = c("16","17","18")
tilev = c("07","08")

#get url
url <- "http://e4ftl01.cr.usgs.gov/MOLT/"

for (i in 5:length(mod.date)){
  url1 <- paste(url, "MOD44B.051", "/",mod.date[i],"/", sep = "")
  doc <- htmlParse(url1)
  links <- xpathSApply(doc, "//a/@href")
  free(doc)
  for(j in 1:length(tileh)){      
    for(k in 1:length(tilev)){
      #only select h17v07/h17v08
      fn = links[which(substr(links, 17, 22) == paste("h",tileh[j],"v",tilev[k], sep = ""))]
      download.file(paste(url1, fn[1], sep = ""), 
                    destfile = paste(getwd(),"/","MOD44B","/", fn[1], sep = ""),
                    mode = "wb")
      
    } # end of k
  } #end if j
} #end of i

#change the hdf into raster files using IDL

#calculate trend and plot
# in local machine
setwd("R:\\users\\Zhihua\\MODIS\\MOD44B")
require("raster")
library("rgdal")

modis.grid.sp = readOGR(dsn="R:\\users\\Zhihua\\MODIS\\MODIS_tile_shp\\modis_sinusoidal",layer="modis_sinusoidal_grid_world")
modis.grid.sp = modis.grid.sp[which(modis.grid.sp@data$h >= 16 & modis.grid.sp@data$h <= 18 & modis.grid.sp@data$v <= 8 & modis.grid.sp@data$v >= 7),]

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b_ghana.sinu = spTransform(county_b_ghana, CRS(proj.sin))
regions.sin <- spTransform(regions, CRS(proj.sin)) 

tileh.info = paste("h", tileh, sep = "") 
tilev.info = paste("v", tilev, sep = "")

tile.info = c(paste(tileh.info, tilev.info[1], sep = ""), paste(tileh.info, tilev.info[2], sep = ""))

tif.fn = list.files(path = "R:\\users\\Zhihua\\MODIS\\MOD44B\\hdf2tif", pattern = "*.Percent_tree_cover.tif$")
tif.fn = tif.fn[which(nchar(tif.fn) < 40)]
tif.stack = stack(paste("R:\\users\\Zhihua\\MODIS\\MOD44B\\hdf2tif\\", tif.fn, sep = ""))

projection(tif.stack) <- projection(modis.grid.sp)
extent(tif.stack) <- extent(modis.grid.sp)

tif.stack = crop(tif.stack , county_b_ghana.sinu)
tif.stack[tif.stack >= 100] = NA

#calculate trend
tif.trd = calc(tif.stack, mk1)

names(tif.trd) <- c("tau", "p", "obs")

tif.trd1 = tif.trd$tau
tif.trd1[tif.trd$tau < 0 & tif.trd$p <= 0.1] = 1
tif.trd1[tif.trd$tau >= 0 & tif.trd$p <= 0.1] = 2
tif.trd1[tif.trd$p > 0.1] = 3
tif.trd1[is.na(tif.trd$p)] = 4

#calculate tree cover difference

tif.trd2 = calc(tif.stack[[c(11:15)]], mean, na.rm = TRUE) - calc(tif.stack[[c(1:5)]], mean, na.rm = TRUE)

writeRaster(tif.trd1,paste("R:/users/Zhihua/MODIS/MOD44B/results/", "MOD44B.VCF.trend.tif", sep = ""), format="GTiff", overwrite=TRUE)
writeRaster(tif.trd2,paste("R:/users/Zhihua/MODIS/MOD44B/results/", "MOD44B.VCF.dif.tif", sep = ""), format="GTiff", overwrite=TRUE)

#plot trend
tif.trd1.geo = raster("./MOD44B/results/MOD44B.VCF.trend.geo.tif")
tif.trd1.geo = tif.trd1.geo*county_b_ghana.r
#tif.trd1.geo[lc_rc == 1|lc_rc == 5|lc_rc == 7] = NA
  
breaks2_change <- 0:4        
legendbrks2_change <- 1:4 - 0.5

arg1 <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels = c("Negative","Positive","No Trend","Not Calculated")
color1=c("#e66101", "#1a9641","#ffffbf","#2b83ba")

#plot MOD44B

png(file = "D:\\users\\Zhihua\\MODIS\\NBAR_results4\\MOD44B_vcf_trend_geo_clip.png", width = 4000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(tif.trd1.geo,
          maxpixels = nrow(tif.trd1.geo)*ncol(tif.trd1.geo),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)),   
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5))

dev.off()

#plot change between last 5 years and first 5 year
tif.trd2.geo = raster("./MOD44B/results/MOD44B.VCF.dif.geo.tif")
tif.trd2_rcls = recls2(tif.trd2.geo,threshold = c(-50, -20, -10, -3, 3, 10, 20, 50))
tif.trd2_rcls = tif.trd2_rcls*county_b_ghana.r

breaks2_change <- 0:9        
legendbrks2_change <- 1:9 - 0.5

labels = c("< -50","- (50-20)","-(20-10)","-(10-5)","-3-3", "5-10","10-20","20-50","> 50")
color1=c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')

png(file = "D:\\users\\Zhihua\\MODIS\\NBAR_results4\\MOD44B_vcf_dif_geo_clip.png", width = 4000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(tif.trd2_rcls,
          maxpixels = nrow(tif.trd2_rcls)*ncol(tif.trd2_rcls),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5))

dev.off()

