# section iv: calculate the mat hensen's Landsat forest loss
library(raster)
setwd("R:\\users\\Zhihua\\MODIS")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
regions = readOGR(dsn="R:/users/Zhihua/MODIS/NBAR_results3", layer="Val_regions")
projection(regions) <- proj.geo 

#read into trend data
TCW.trd2.grd = raster("./NBAR_results4/TCW2.dry.trend.wa.tif")
modis.ext = extent(TCW.trd2.grd)

#read into hansen's data
lt_chips = c("10N_000E", "10N_010W", "10N_020W","20N_020W", "20N_010W","20N_000E")

for (i in 1:6){
  
  fl = raster(paste("./GlobalForestCover/Hansen_GFC2015_loss_", lt_chips[i], ".tif", sep = ""))
  fg = raster(paste("./GlobalForestCover/Hansen_GFC2015_gain_", lt_chips[i], ".tif", sep = ""))
  
  landsat.ext = extent(fl)
  intersect.ext = intersect(modis.ext, landsat.ext)
  fl1 = crop(fl, intersect.ext)
  fg1 = crop(fg, intersect.ext)
  
  TCW.trd2.grd1 = crop(TCW.trd2.grd, intersect.ext)
  
  nr = nrow(TCW.trd2.grd1)
  nc = ncol(TCW.trd2.grd1)
  
  TCW.trd2.grd2 = TCW.trd2.grd1
  TCW.trd2.grd2[] <- 1:(nr*nc)
  
  #change resolutions
  TCW.trd2.grd3 <- resample(TCW.trd2.grd2, fl1, method='ngb') 
  
  # calculate loss and gain for each zone
  fl1.1 = zonal(fl1, TCW.trd2.grd3, fun='sum')
  fg1.1 = zonal(fg1, TCW.trd2.grd3, fun='sum')
  
  fl1.grd <- TCW.trd2.grd2
  fl1.grd[] <- fl1.1[,2]
  
  fg1.grd <- TCW.trd2.grd2
  fg1.grd[] <- fg1.1[,2]
  
  net.grd = fl1.grd - fg1.grd
  
  writeRaster(fl1.grd,paste("./GlobalForestCover/Hansen_GFC2015_loss_modis_grid", lt_chips[i], ".tif", sep = ""), format="GTiff", overwrite=TRUE)
  writeRaster(fg1.grd,paste("./GlobalForestCover/Hansen_GFC2015_gain_modis_grid", lt_chips[i], ".tif", sep = ""), format="GTiff", overwrite=TRUE)
  writeRaster(net.grd,paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[i], ".tif", sep = ""), format="GTiff", overwrite=TRUE)
  
  print(paste("Finish calculating region ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#combine the results to together

net1 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[1], ".tif", sep = ""))
net2 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[2], ".tif", sep = ""))
net3 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[3], ".tif", sep = ""))
net4 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[4], ".tif", sep = ""))
net5 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[5], ".tif", sep = ""))
net6 = raster(paste("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid", lt_chips[6], ".tif", sep = ""))

net = mosaic(net1, net2, net3, net4,net5,net6,fun=mean)
net = net*county_b.geo.r
net = net/289
net[net < 0] = NA
writeRaster(net,"./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid.tif", format="GTiff", overwrite=TRUE)

#plot 
net = raster("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid.tif")
net_rcls = recls2(net,threshold = c(0, 0.05, 0.1, 0.2, 0.5))

net_rcls = net_rcls*county_b_ghana.r
net_rcls[lc_rc == 1|lc_rc == 5|lc_rc == 7] = NA

breaks2_change <- 0:6        
legendbrks2_change <- 1:6 - 0.5
color1=c('#66c2a4','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')
labels = c("No Loss","0-0.05","0.05-0.1","0.1-0.2","0.2-0.5",">0.5")


png(file = "D:\\users\\Zhihua\\MODIS\\NBAR_results4\\Landsat_forest_netloss_cliped.png", width = 4000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(net_rcls,
          maxpixels = nrow(net_rcls)*ncol(net_rcls),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5)) 

dev.off()
