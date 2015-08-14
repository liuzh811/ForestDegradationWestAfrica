######################################################################################
# find forest cells
#1 read data 
# MODIS land cover climatotology
lc.modis2 = raster("D:\\users\\Zhihua\\MODIS\\landcover\\LCType.tif")
lc.modis2 = crop(lc.modis2,county_b_ghana)

# ESA GlobCover 2009 data
lc.esa = raster("D:\\users\\Zhihua\\GeoData_West Africa\\GlobCover2009\\Globcover2009_V2.3_Global_2009\\GLOBCOVER_L4_200901_200912_V2.3.tif")
lc.esa = crop(lc.esa, county_b_ghana)

# ESA GlobCover 2000 data
lc.esa.2000 = readGDAL("D:\\users\\Zhihua\\GeoData_West Africa\\GLC2000_Africa_Grid\\africa_v5")
lc.esa.2000 = raster(lc.esa.2000)
lc.esa.2000 = crop(lc.esa.2000, county_b_ghana)

#2 extract data
lc.modis2.df = raster::extract(lc.modis2, pts.sp) #based on global land cover climatology
lc.esa.df = raster::extract(lc.esa, pts.sp)
lc.esa.2000.df = raster::extract(lc.esa.2000, pts.sp)

#3 find forest data
lcc = data.frame(lc.modis = lc.modis2.df, lc.esa =  lc.esa.df,lc.esa.2000 = lc.esa.2000.df)
lcc$lc = 0
lcc$lc[which((lcc$lc.modis == 2 & lcc$lc.esa == 40) | (lcc$lc.esa.2000 == 1 & lcc$lc.modis == 2))] = 1
lcc.grd = Point2raster(lcc$lc, raster = TCW1.dry[[1]]) #use this one as final tropical forest extent
writeRaster(lcc.grd,".\\NBAR_results2\\lc.forest.tif", format="GTiff", overwrite=TRUE)

#4 plot extend of change
#get degradation (=1), increasing(=2), No trend(=3), NOT calculated because of not enough data (=4), NON-forest (=NA)

TCW1.dry.trd1 = stack(".\\NBAR_results2\\TCW1.dry.trend.wa.grd") 
TCW2.dry.trd1 = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.grd") 
TCW3.dry.trd1 = stack(".\\NBAR_results2\\TCW3.dry.trend.wa.grd") 

trend.vi11.df = raster::extract(TCW1.dry.trd1, pts.sp) 
trend.vi21.df = raster::extract(TCW2.dry.trd1, pts.sp) 
trend.vi31.df = raster::extract(TCW3.dry.trd1, pts.sp) 

#for quality 1
lcc.df = data.frame(lcc, trend.vi11.df)
colnames(lcc.df)[c(5:7)] <- c("tau","p","length")
lcc.df$trend = NA
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau < 0 & lcc.df$p <= 0.1)] = 1
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau >= 0 & lcc.df$p <= 0.1)] = 2
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$p > 0.1)] = 3
lcc.df$trend[which(lcc.df$lc == 1 & is.na(lcc.df$p))] = 4

trd.grd = Point2raster(lcc.df$trend, raster = TCW1.dry[[1]]) #use this one as final tropical forest extent
writeRaster(trd.grd,".\\NBAR_results2\\trend.wa.1.tif", format="GTiff", overwrite=TRUE)

table(lcc.df$trend)
> table(lcc.df$trend) #0.06068643

1      2      3      4 
50523  17078 357754  27481
#degradatoion rate: 50523/(452836-27481) = 0.1187784
#increasing rate: 17078/(452836-27481) = 0.04014999
#cloud #0.06068643 27481/452836

#for quality 2
lcc.df = data.frame(lcc, trend.vi21.df)
colnames(lcc.df)[c(5:7)] <- c("tau","p","length")
lcc.df$trend = NA
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau < 0 & lcc.df$p <= 0.1)] = 1
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau >= 0 & lcc.df$p <= 0.1)] = 2
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$p > 0.1)] = 3
lcc.df$trend[which(lcc.df$lc == 1 & is.na(lcc.df$p))] = 4

trd.grd = Point2raster(lcc.df$trend, raster = TCW1.dry[[1]]) #use this one as final tropical forest extent
writeRaster(trd.grd,".\\NBAR_results2\\trend.wa.2.tif", format="GTiff", overwrite=TRUE)

table(lcc.df$trend)
> table(lcc.df$trend)

1      2      3      4 
54508  21284 369296   7748

54508/(452836-7748)
21284/(452836-7748)
7748/452836
#for quality 3
lcc.df = data.frame(lcc, trend.vi31.df)
colnames(lcc.df)[c(5:7)] <- c("tau","p","length")
lcc.df$trend = NA
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau < 0 & lcc.df$p <= 0.1)] = 1
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$tau >= 0 & lcc.df$p <= 0.1)] = 2
lcc.df$trend[which(lcc.df$lc == 1 & lcc.df$p > 0.1)] = 3
lcc.df$trend[which(lcc.df$lc == 1 & is.na(lcc.df$p))] = 4

trd.grd = Point2raster(lcc.df$trend, raster = TCW1.dry[[1]]) #use this one as final tropical forest extent
writeRaster(trd.grd,".\\NBAR_results2\\trend.wa.3.tif", format="GTiff", overwrite=TRUE)

table(lcc.df$trend)
> table(lcc.df$trend)

1      2      3 
72862  26839 353135
72862/452836
26839/452836

###########
