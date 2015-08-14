#####################################################################################################
#####   based on one month lag: from de Wasseige et al (2003) and Zhou et al (2014)    ##############
#####################################################################################################

## Section 2: residual trend analysis
## 1. select the suitable temporal scale for rainfall summarization
# 1.1 de Wasseige et al (2003 Agricultural and Forest Meteorology) described that LAI lag one month after dry season starts
#     in moist deciduous Congo–Guinean forest. [cited by Zhou et al 2014 nature congo forest]
# 1.2 Zhou et al. (2014) find April-May-June EVI correlated strongly with April-May-June and March-Apr-May rainfall.
#     [no lag & one month before]
#1.3 Enquist, B. J. & Enquist, C. A. F. Long-term change within a neotropical forest: assessing differential 
#    functional and floristic responses to disturbance and drought. Glob. Change Biol. 17, 1408–1424 (2011) 
#   Enquist and Enquist note that the tropical forest had different sensistivity to rainfall change, therefore we 
# examined up to 3 month before the dry season

# 3.1. get points
TCW2.dry = stack(".\\NBAR_results2\\TCW2.dry.wa.grd")
pts.sp = Ex.pts.all(TCW2.dry[[1]])

# 3.2 extract data
#extract rainfall
rf.1mon.df = raster::extract(rf.1mon, pts.sp)

# extract TCW
tcw.df = raster::extract(TCW2.dry, pts.sp)

# extract forest points
lcc.grd = raster(".\\NBAR_results2\\lc.forest.tif")
lcc.df = raster::extract(lcc.grd, pts.sp) 
lcc.forest.idx = which(lcc.df == 1)
lcc.forest.idx1 = which(lcc.df == 0) #non foret

#replace non-foret cells as NA to speedup calculation
rf.1mon.df[lcc.forest.idx1,] = NA
tcw.df[lcc.forest.idx1,] = NA

# 3.3 residual analysis, using 1 month before season rainfall
tcw.1.df = data.frame(tcw.df, rf.1mon.df)

tcw.rf.cor = apply(tcw.1.df, 1, fun.cor2)                                #based on strongest relationship
tcw.df.pred = sweep(rf.1mon.df, MARGIN = 1, STATS = tcw.rf.cor[1,], FUN="*") #predict from rainfall data
tcw.df.pred = sweep(tcw.df.pred, MARGIN = 1, STATS = tcw.rf.cor[2,], FUN="+") #add the intercept
tcw.df.residual = tcw.df - tcw.df.pred                                        #calculate residual
tcw.df.residual.trd <- apply(tcw.df.residual, 1, mk1)                             #calculate residual trend

# 3.4 convert to raster, and classify
tcw.df.residual.trd.value = as.numeric(tcw.df.residual.trd[1,])                                    
tcw.df.residual.trd.p = as.numeric(tcw.df.residual.trd[2,])                                    
# change to raster
corr = Point2raster(tcw.rf.cor[1,], raster = TCW2.dry[[1]])
corr[lcc.grd != 1] = NA

tcw.df.residual.trd.value = Point2raster(tcw.df.residual.trd.value, raster = TCW2.dry[[1]])
tcw.df.residual.trd.p = Point2raster(tcw.df.residual.trd.p, raster = TCW2.dry[[1]])

tcw.df.residual.trd.value2 = tcw.df.residual.trd.value
tcw.df.residual.trd.value2[tcw.df.residual.trd.value < 0 & tcw.df.residual.trd.p <= 0.1] = 1
tcw.df.residual.trd.value2[tcw.df.residual.trd.value >= 0 & tcw.df.residual.trd.p <= 0.1] = 2
tcw.df.residual.trd.value2[tcw.df.residual.trd.p > 0.1] = 3
tcw.df.residual.trd.value2[lcc.grd == 1&is.na(tcw.df.residual.trd.p)] = 4

#3.5 save results
writeRaster(tcw.df.residual.trd.value2,".\\NBAR_results2\\tcw.df.residual.trd.wa2.tif", format="GTiff", overwrite=TRUE)
writeRaster(corr,".\\NBAR_results2\\tcw.rainfall.cor.wa3.tif", format="GTiff", overwrite=TRUE)

##compare residual trend and raw trend
trd.grd2 = raster(".\\NBAR_results2\\trend.wa.2.tif")
trd.comparsion = trd.grd2 == tcw.df.residual.trd.value2
freq(trd.comparsion)
