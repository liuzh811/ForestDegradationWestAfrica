#1. read into TRMM data
fn = list.files(path = "D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2", pattern = "*.nc$")
r = stack(paste("D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2\\", fn, sep = ""))
#r = crop(r, county_b.geo) #crop the dataset to western Africa boundary
r = crop(r, county_b_ghana) #crop the dataset to western Africa boundary
names(r) <- fn

#2. calulate rainfall 1, 2,3,4,5,6 month before dry season, and during the dry season
rf.1mon = list() #1 month before
rf.2mon = list() #2 month before
rf.3mon = list() #3 month before
rf.8mon = list() #rainy season rainfall

rf.4mon = list()
rf.5mon = list()
rf.6mon = list()
rf.7mon = list()
rf.9mon = list()
rf.10mon = list()
rf.dry = list()


for(i in 2000:2014){
  idx.1 = which(as.numeric(substr(names(r), 13,16))==i & as.numeric(substr(names(r), 18,19))==11 & as.numeric(substr(names(r), 21,22))==15)
  idx.1mon = (idx.1-30):(idx.1+60)
  idx.2mon = (idx.1-60):(idx.1+30)
  idx.3mon = (idx.1-90):idx.1
  idx.8mon = (idx.1-240):idx.1
  
  rf.1mon[[i-1999]] <- calc(r[[idx.1mon]], sum) 
  rf.2mon[[i-1999]] <- calc(r[[idx.2mon]], sum) 
  rf.3mon[[i-1999]] <- calc(r[[idx.3mon]], sum) 
  rf.8mon[[i-1999]] <- calc(r[[idx.8mon]], sum) 
  
  
  idx.4mon = (idx.1-120):idx.1
  idx.5mon = (idx.1-150):idx.1
  idx.6mon = (idx.1-180):idx.1
  idx.7mon = (idx.1-210):idx.1
  idx.9mon = (idx.1-270):idx.1
  idx.10mon = (idx.1-300):idx.1
  
  if(i < 2014){idx.dry = idx.1:(idx.1+135)} else (idx.dry = idx.1:(idx.1+70))
  rf.4mon[[i-1999]] <- calc(r[[idx.4mon]], sum)
  rf.5mon[[i-1999]] <- calc(r[[idx.5mon]], sum)
  rf.6mon[[i-1999]] <- calc(r[[idx.6mon]], sum)
  rf.7mon[[i-1999]] <- calc(r[[idx.7mon]], sum)
  rf.9mon[[i-1999]] <- calc(r[[idx.9mon]], sum)
  rf.10mon[[i-1999]] <- calc(r[[idx.10mon]], sum)
  rf.dry[[i-1999]] <- calc(r[[idx.dry]], sum)
  
}

rf.1mon = stack(rf.1mon)
rf.2mon = stack(rf.2mon)
rf.3mon = stack(rf.3mon)
rf.8mon = stack(rf.8mon)

rf.4mon = stack(rf.4mon)
rf.5mon = stack(rf.5mon)
rf.6mon = stack(rf.6mon)
rf.7mon = stack(rf.7mon)
rf.9mon = stack(rf.9mon)
rf.10mon = stack(rf.10mon)
rf.dry = stack(rf.dry)


#calculate trends
rf.1mon.trd <- calc(rf.1mon, mk1)
rf.2mon.trd <- calc(rf.2mon, mk1)
rf.3mon.trd <- calc(rf.3mon, mk1)
rf.8mon.trd <- calc(rf.8mon, mk1)

#only retain significant trend area
rf.1mon.trd[[1]][rf.1mon.trd[[2]] > 0.1] = NA
rf.2mon.trd[[1]][rf.2mon.trd[[2]] > 0.1] = NA
rf.3mon.trd[[1]][rf.3mon.trd[[2]] > 0.1] = NA
rf.8mon.trd[[1]][rf.8mon.trd[[2]] > 0.1] = NA

trend.rain = stack(rf.1mon.trd[[1]],rf.2mon.trd[[1]],rf.3mon.trd[[1]],rf.8mon.trd[[1]])

names(trend.rain) <- c(paste(1:3, ".month.before.dry.season", sep = ""),"Rainy.Season")

png(file = ".\\NBAR_results2\\rainfall.trend.png", width = 4000, height = 3000, units = "px", res = 300)

levelplot(trend.rain)+
  layer(sp.polygons(county_b.geo, col = "black", lwd = 1)) 

dev.off()

# 3. residual analysis
# 3. find the strongest relationship between rainfall and VI
# 3.1. get points
TCW2.dry = stack(".\\NBAR_results2\\TCW2.dry.wa.grd")
pts.sp = Ex.pts.all(TCW2.dry[[1]])

#extract rainfall
rf.1mon.df = raster::extract(rf.1mon, pts.sp)
rf.2mon.df = raster::extract(rf.2mon, pts.sp)
rf.3mon.df = raster::extract(rf.3mon, pts.sp)
rf.8mon.df = raster::extract(rf.8mon, pts.sp)

rf.4mon.df = raster::extract(rf.4mon, pts.sp)
rf.5mon.df = raster::extract(rf.5mon, pts.sp)
rf.6mon.df = raster::extract(rf.6mon, pts.sp)
rf.7mon.df = raster::extract(rf.7mon, pts.sp)
rf.9mon.df = raster::extract(rf.9mon, pts.sp)
rf.10mon.df = raster::extract(rf.10mon, pts.sp)
rf.dry.df = raster::extract(rf.dry, pts.sp)

#extract vi
tcw.df = raster::extract(TCW2.dry, pts.sp)

#extract forest points
lcc.grd = raster(".\\NBAR_results2\\lc.forest.tif")
lcc.df = raster::extract(lcc.grd, pts.sp) 

#compute the trend for forest pixels
lcc.forest.idx = which(lcc.df == 1)
lcc.forest.idx1 = which(lcc.df == 0) #non foret

#replace non-foret cells as NA
rf.1mon.df[lcc.forest.idx1,] = NA
rf.2mon.df[lcc.forest.idx1,] = NA
rf.3mon.df[lcc.forest.idx1,] = NA
rf.8mon.df[lcc.forest.idx1,] = NA
tcw.df[lcc.forest.idx1,] = NA

rf.4mon.df[lcc.forest.idx1,] = NA
rf.5mon.df[lcc.forest.idx1,] = NA
rf.6mon.df[lcc.forest.idx1,] = NA
rf.7mon.df[lcc.forest.idx1,] = NA
rf.9mon.df[lcc.forest.idx1,] = NA
rf.10mon.df[lcc.forest.idx1,] = NA
rf.dry.df[lcc.forest.idx1,] = NA


rf.vi = data.frame(rf.dry = apply(rf.dry.df, 2,median, na.rm = T),
                   rf.1mon = apply(rf.1mon.df, 2,median, na.rm = T),
                   rf.2mon = apply(rf.2mon.df, 2,median, na.rm = T),
                   rf.3mon = apply(rf.3mon.df, 2,median, na.rm = T),
                   rf.4mon = apply(rf.4mon.df, 2,median, na.rm = T),
                   rf.5mon = apply(rf.5mon.df, 2,median, na.rm = T),
                   rf.6mon = apply(rf.6mon.df, 2,median, na.rm = T),
                   rf.7mon = apply(rf.7mon.df, 2,median, na.rm = T),
                   rf.8mon = apply(rf.8mon.df, 2,median, na.rm = T),
                   rf.9mon = apply(rf.9mon.df, 2,median, na.rm = T),
                   rf.10mon = apply(rf.10mon.df, 2,median, na.rm = T),
                   tcw = apply(tcw.df, 2,median, na.rm = T))   
colnames(rf.vi) <- c("No.lag","1.Month.lag","2.Month.lag","3.Month.lag",
                              "4.Month.lag","5.Month.lag","6.Month.lag",
                              "7.Month.lag","8.Month.lag","9.Month.lag","10.Month.lag","TCW")

write.csv(rf.vi,"rainfall.TCW.relation.csv")

library(psych)
png(file = ".\\NBAR_results2\\rainfall.TCW.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.vi)

dev.off()

png(file = ".\\NBAR_results2\\rainfall.TCW.relation2.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(rf.vi[,c(1,2,3,4,12)])

dev.off()

