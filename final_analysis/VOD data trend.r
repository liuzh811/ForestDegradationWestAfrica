##VOD data to do the analysis
#download VOD from: http://www.wenfo.org/wald/global-biomass/
#refs: Liu, YY, van Dijk, AI, de Jeu, RA, Canadell, JG, McCabe, MF, Evans, JP & Wang, G (2015) 
# Recent reversal in loss of global terrestrial biomass. Nature Climate Change,
# in local machine
library(rgdal)
library(raster)
library(dismo)
library(rasterVis)
library(Kendall)
library("ncdf") ##only for netcdf version 4 or eailer

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
county_b = readOGR(dsn="R:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))
county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#vod data
vod = brick("R:\\users\\Zhihua\\GeoData_West Africa\\VOD\\Global_annual_mean_VOD_1993_2012_20150331.nc")
#vod = vod[[c(9:20)]]

vod2 = list()
for(i in 1:nlayers(vod)){vod2[[i]] <- flip(raster(t(as.matrix(vod[[i]]))), direction='x')}

vod2 = stack(vod2)
extent(vod2) <- c(extent(vod)@ymin,extent(vod)@ymax, extent(vod)@xmin, extent(vod)@xmax)
projection(vod2) <- projection(vod)
#names(vod2) <- paste("Y", 2001:2012, sep="")
names(vod2) <- paste("Y", 1993:2012, sep="")

vod2 = crop(vod2, county_b_ghana)

vod_trd = calc(vod2, mk1)
names(vod_trd) <- c("tau", "p", "obs")

vod_trd.grd = vod_trd$tau
vod_trd.grd[vod_trd$tau < 0 & vod_trd$p <= 0.1] = 1
vod_trd.grd[vod_trd$tau >= 0 & vod_trd$p <= 0.1] = 2
vod_trd.grd[vod_trd$p > 0.1] = 3
vod_trd.grd[is.na(vod_trd$p)] = 4


#biomass data
bio = brick("R:\\users\\Zhihua\\GeoData_West Africa\\VOD\\Global_annual_mean_ABC_lc2001_1993_2012_20150331.nc")
bio = bio[[c(9:20)]]

bio2 = list()
for(i in 1:nlayers(bio)){bio2[[i]] <- flip(raster(t(as.matrix(bio[[i]]))), direction='x')}

bio2 = stack(bio2)
extent(bio2) <- c(extent(bio)@ymin,extent(bio)@ymax, extent(bio)@xmin, extent(bio)@xmax)
projection(bio2) <- projection(bio)
names(bio2) <- paste("Y", 2001:2012, sep="")


#county_b_ghana = extent(-13.5, 1.2, 4.35, 9.5)
bio2 = crop(bio2, county_b_ghana)

bio_trd = calc(bio2, mk1)
names(bio_trd) <- c("tau", "p", "obs")

bio_trd.grd = bio_trd$tau
bio_trd.grd[bio_trd$tau < 0 & bio_trd$p <= 0.1] = 1
bio_trd.grd[bio_trd$tau >= 0 & bio_trd$p <= 0.1] = 2
bio_trd.grd[bio_trd$p > 0.1] = 3
bio_trd.grd[is.na(bio_trd$p)] = 4

#plot
county_b_ghana.r2 = rasterize(county_b_ghana, vod_trd.grd)
county_b_ghana.r2 = county_b_ghana.r2 >= 1
county_b_ghana.r2[county_b_ghana.r2 != 1] = NA
bio_trd.grd = bio_trd.grd*county_b_ghana.r2
vod_trd.grd = vod_trd.grd*county_b_ghana.r2

breaks2_change <- 0:4        
legendbrks2_change <- 1:4 - 0.5

arg1 <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels = c("Negative","Positive","No Trend","Not Calculated")
color1=c("#e66101", "#1a9641","#ffffbf","#2b83ba")

trd = stack(vod_trd.grd, bio_trd.grd)
names(trd) <- c("VOD", "Biomass")
#using levelplot methods
png(file = "R:\\users\\Zhihua\\MODIS\\NBAR_results4\\VOD.trend_1993-2012.png", width = 3000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(vod_trd.grd,
          maxpixels = nrow(vod_trd.grd)*ncol(vod_trd.grd),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) 

dev.off()
