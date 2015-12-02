## section 6, plot start from here
## 6.1 plot vegetation trend
TCW.trd2.grd = raster(".\\NBAR_results4\\TCW2.dry.trend.wa.tif")
EVI.trd2.grd = raster(".\\NBAR_results4\\EVI2.dry.trend.wa.tif")

#get boundry
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

protected = readOGR(dsn="D:\\users\\Zhihua\\GeoData_West Africa\\protectarea",layer="protected_areas_west_africa")

county_b.geo.r = rasterize(county_b.geo, TCW.trd2.grd)
county_b.geo.r = county_b.geo.r >= 1

#highlight some regions
#regions = readOGR(dsn=".\\NBAR_results3", layer="Val_regions")
#projection(regions) <- proj.geo 

#clip to shape
TCW.trd2.grd = TCW.trd2.grd*county_b_ghana.r
TCW.trd2.grd[lc_rc == 1|lc_rc == 5|lc_rc == 7] = NA

EVI.trd2.grd = EVI.trd2.grd*county_b_ghana.r
EVI.trd2.grd[lc_rc == 1|lc_rc == 5|lc_rc == 7] = NA

trd = stack(TCW.trd2.grd, EVI.trd2.grd)
names(trd) <- c("TCW","EVI")

breaks2_change <- 0:4        
legendbrks2_change <- 1:4 - 0.5

arg1 <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels = c("Negative","Positive","No Trend","Not Calculated")
color1=c("#e66101", "#1a9641","#ffffbf","#2b83ba")

#using levelplot methods
png(file = ".\\NBAR_results4\\trend.wa.TCW&EVI_cliped.png", width = 3000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.5, lines=2, fontface='bold')
levelplot(trd,
          maxpixels = nrow(trd)*ncol(trd),
          at= breaks2_change, margin=FALSE,
          col.regions= color1,
          colorkey= list(labels= list(labels= labels,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5)) +
  latticeExtra::layer(sp.polygons(regions[c(2,4,7), ], col = "red", lwd = 2.5)) 
#latticeExtra::layer(sp.polygons(protected, col = "blue", lwd = 1))

dev.off()

#6.2 plot highlighted regions

regions = readOGR(dsn="D:/users/Zhihua/MODIS/NBAR_results3", layer="Val_regions")
projection(regions) <- proj.geo 
regions.sin <- spTransform(regions, CRS(proj.sin)) 


TCW.trd2.grd = raster(".\\NBAR_results4\\TCW2.dry.trend.wa.tif")
EVI.trd2.grd = raster(".\\NBAR_results4\\EVI2.dry.trend.wa.tif")

#read into mat hansen's landsat forest loss and reclassify
landsat_netloss = raster("./GlobalForestCover/Hansen_GFC2015_net_loss_modis_grid.tif")

landsat_netloss_rcls = recls2(landsat_netloss,threshold = c(0, 0.05, 0.1, 0.2, 0.5))
breaks2_change_hansen <- 0:6        
legendbrks2_change_hansen <- 1:6 - 0.5
color1_hansen=c('#66c2a4','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')
labels_hansen = c("No Loss","0-0.05","0.05-0.1","0.1-0.2","0.2-0.5",">0.5")

#read into MOD44B VCF and reclassify
vcf_trend = raster("./MOD44B/results/MOD44B.VCF.trend.geo.tif")

vcf_dif = raster("./MOD44B/results/MOD44B.VCF.dif.geo.tif")
vcf_dif_rcls = recls2(vcf_dif,threshold = c(-50, -20, -10, -3, 3, 10, 20, 50))

breaks2_change_vcf <- 0:9        
legendbrks2_change_vcf <- 1:9 - 0.5
labels_vcf = c("< -50","- (50-20)","-(20-10)","-(10-5)","-3-3", "5-10","10-20","20-50","> 50")
color1_vcf=c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')

#trend color scheme
breaks2_change_trd <- 0:4        
legendbrks2_change_trd <- 1:4 - 0.5
arg1_trd <- list(at=seq(1,4,1), labels=c("Negative","Positive","No Trend","Not Calculated")) #these are the class names
labels_trd = c("Negative","Positive","No Trend","Not Calculated")
color1_trd = c("#e66101", "#1a9641","#ffffbf","#2b83ba")

#trend combination
trd_com = raster(".\\NBAR_results4\\Trend_com_data.tif")
breaks2_change_trd.com <- 0:9        
legendbrks2_change_trd.com <- 1:9 - 0.5
trd.arg_trd.com = c("HLD", "WDI","WHLG","HDI","HLG", "WDD","WHLD","HDD","NoTrend")
arg1_trd.com <- list(at=seq(1,9,1), labels=trd.arg_trd.com) #these are the class names
labels_trd.com = trd.arg_trd.com
color1_trd.com <- c(brewer.pal(9,"Reds")[c(3,6,9)], brewer.pal(9,"Greens")[c(3,6,9)],brewer.pal(9,"Blues")[c(3,6)],"#cccccc")

#plot start from here
regionid = c(2,4,7)

for (i in regionid){

  TCW.trd2.grd1 = crop(TCW.trd2.grd, regions[i,])
  EVI.trd2.grd1 = crop(EVI.trd2.grd, regions[i,])
  landsat_netloss_rcls1 = crop(landsat_netloss_rcls, regions[i,])
  
  vcf_trend1 = crop(vcf_trend, regions[i,])
  vcf_dif_rcls1 = crop(vcf_dif_rcls, regions[i,])
  trd_com1 = crop(trd_com, regions[i,])
  
  png(file = paste("NBAR_results4/Highlighted_region", i, ".png", sep = ""), width = 4500, height = 3000, units = "px", res = 300)
  
par(mfrow=c(2,3),mar=c(0.2, 0.2, 0.2, 0.2))
#tcw
plot(TCW.trd2.grd1,col = color1_trd[c(1:3)], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
ext = extent(TCW.trd2.grd1)
legend("topright",
  #x = ext@xmax - 0.16, y = ext@ymin+0.1, 
  legend = labels_trd[c(1:3)], 
  fill = color1_trd[c(1:3)],  cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white",title="a:) TCW trend")

#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "a:) TCW trend", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")


#evi
plot(EVI.trd2.grd1,col = color1_trd[c(1:3)], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
legend(#x = ext@xmax - 0.16, y = ext@ymin+0.1, 
        "topright",
       legend = labels_trd[c(1:3)], 
       fill = color1_trd[c(1:3)],  cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white",title="b:) EVI trend")
#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "b:) EVI trend", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")

#vcf_trend
plot(vcf_trend1,col = color1_trd[c(1:3)], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
ext = extent(vcf_trend1)
legend("topright",
  #x = ext@xmax - 0.16, y = ext@ymin+0.1, 
  legend = labels_trd[c(1:3)], 
       fill = color1_trd[c(1:3)],  cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white",title="c:) VCF trend")
#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "c:) VCF trend", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")

#vcf_dif_rcls1
plot(vcf_dif_rcls1,col = color1_vcf[sort(unique(vcf_dif_rcls1))], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
legend("topright",
  #x = ext@xmax - 0.16, y = ext@ymax-0.1, 
  legend = labels_vcf[sort(unique(vcf_dif_rcls1))], 
       fill = color1_vcf[sort(unique(vcf_dif_rcls1))],  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "white",title="c:) VCF Change")
#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "d:) VCF Change", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")

#plot.new()
#legend(x = 0, y = 0.9, legend = "d): VCF change",  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "transparent")
#legend(x = 0.05, y = 0.83, legend = labels_vcf, fill = color1_vcf,  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "transparent")

#legend(x = 0.5, y = 0.9, legend = "e): Percent\n of foret loss",  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "transparent")
#legend(x = 0.55, y = 0.73, legend = labels_hansen, fill = color1_hansen,  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "transparent")

#landsat

plot(landsat_netloss_rcls1,col = color1_hansen[sort(unique(landsat_netloss_rcls1))], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
legend("topright",
  #x = ext@xmax - 0.16, y = ext@ymax-0.1, 
  legend = labels_hansen[sort(unique(landsat_netloss_rcls1))], 
       fill = color1_hansen[sort(unique(landsat_netloss_rcls1))],  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "white",title="e:) Forest Loss")
#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "e:) Forest Loss", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")

#trend combination
plot(trd_com1,col = color1_trd.com[sort(unique(trd_com1))], 
     legend=FALSE,
     axes=FALSE,
     box=FALSE,
)
plot(county_b_ghana, add = TRUE)
legend("topright",
  #x = ext@xmax - 0.16, y = ext@ymax-0.1, 
  legend = labels_trd.com[sort(unique(trd_com1))], 
       fill = color1_trd.com[sort(unique(trd_com1))],  cex = 1.75,  box.lwd = 0,box.col = "transparent",bg = "white",title="f:) Change Attribution")
#legend(x = (ext@xmin + ext@xmax)/2, y = ext@ymax-0.01, "f:) Change Attribution", cex = 1.75,  box.lwd = 0,box.col = "white",bg = "white")


dev.off()

}

#6.3 plot number of clear obs
#for number of clear obs to calculate indice for each year

  tmp1 = stack(paste(".\\NBAR_results4\\Obs.num.level2.wa.subregion.1.grd", sep = ""))
  tmp2 = stack(paste(".\\NBAR_results4\\Obs.num.level2.wa.subregion.2.grd", sep = ""))
  tmp3 = stack(paste(".\\NBAR_results4\\Obs.num.level2.wa.subregion.3.grd", sep = ""))
  tmp4 = stack(paste(".\\NBAR_results4\\Obs.num.level2.wa.subregion.4.grd", sep = ""))
  
  list2.dry = mosaic(tmp1, tmp2, tmp3, tmp4,fun=mean)
  list2.dry = list2.dry*county_b_ghana.r
  
  names(list2.dry) <- paste("Y", 2001:2015, sep = "")

png(file = ".\\NBAR_results4\\Obs.num.qc.level2.png", width = 5000, height = 5000, units = "px", res = 300)
  
recls2plot(list2.dry, threshold = c(2, 4, 6, 8, 10, 12), color = "Spectral")

dev.off()

#for number of clear obs for each day of year cross 18 years
tmp1 = stack(paste(".\\NBAR_results4\\Obs.num.date.level2.wa.subregion.1.grd", sep = ""))
tmp2 = stack(paste(".\\NBAR_results4\\Obs.num.date.level2.wa.subregion.2.grd", sep = ""))
tmp3 = stack(paste(".\\NBAR_results4\\Obs.num.date.level2.wa.subregion.3.grd", sep = ""))
tmp4 = stack(paste(".\\NBAR_results4\\Obs.num.date.level2.wa.subregion.4.grd", sep = ""))

obs.date2 = mosaic(tmp1, tmp2, tmp3, tmp4,fun=mean)
obs.date2 = obs.date2*county_b_ghana.r

names(obs.date2) <- paste("DOY", c(seq.int(321, 361, by = 8), "001", "009", paste("0", seq.int(17, 89, by = 8), sep = "")), sep = "")

png(file = ".\\NBAR_results4\\Obs.num.DOY.qc.level2.png", width = 5000, height = 5000, units = "px", res = 300)

recls2plot(obs.date2, threshold = c(2, 4, 6, 8, 10, 12), color = "Spectral")

dev.off()

# for number of clear obs to calculate thre trend
tmp1 = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.subregion.1.grd")
tmp2 = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.subregion.2.grd")
tmp3 = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.subregion.3.grd")
tmp4 = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.subregion.4.grd")

Numb = mosaic(tmp1[[3]], tmp2[[3]], tmp3[[3]], tmp4[[3]],fun=mean)
Numb.recls = recls2(Numb,threshold = c(4,6,8,10,12))
# Numb.recls[lcc.grd != 1] = NA
Numb.recls = Numb.recls*county_b_ghana.r

#set color scheme
threshold = c(4,6,8,10,12)
color = "Spectral"

clscolor_change = brewer.pal(length(threshold)+2,color)[c(2:(length(threshold)+2))]
clscolor_change = c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837")

breaks2_change <- 0:(length(threshold)+1)        
legendbrks2_change <- 1:(length(threshold)+1) - 0.5

clasnames_change = c(paste("<=", threshold[1], sep = " "))
for(j in 1:(length(threshold)-1)){
  
  clasnames_change = c(clasnames_change, paste(threshold[j], "-", threshold[j+1], sep = " "))
  
}
clasnames_change[length(threshold)+1] <- c(paste(">", threshold[length(threshold)], sep = " "))

png(file = ".\\NBAR_results4\\trend.numb.wa.png", width = 3000, height = 3000, units = "px", res = 300)

p.strip <- list(cex=1.3, lines=2)
levelplot(Numb.recls,
          maxpixels = nrow(Numb.recls)*ncol(Numb.recls),
          at= breaks2_change, margin=FALSE,
          col.regions= clscolor_change,
          colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5)),
          scales=list(x=list(cex=1.3),y=list(cex=1.3)), 
          xlab=list(label = "Longtitude", cex=1.3),ylab=list(label = "Latitude", cex=1.3),
          par.strip.text=p.strip) +
  latticeExtra::layer(sp.polygons(county_b_ghana, col = "black", lwd = 1.5))
  
dev.off()
