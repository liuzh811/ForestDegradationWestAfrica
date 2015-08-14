# calculate spectral indices
TCB = list() #TC brightness 
TCG = list() #TC greenness  
TCW = list() #TC wetness  
TCA = list() #TC angle 
NDWI = list() #normalized difference moisture index
EVI = list() #normalized difference moisture index

for (i in 1:nlayers(band11) ){
  
  brightness = 0.4395*band11[[i]]+0.5945*band21[[i]]+0.246*band31[[i]]+0.3918*band41[[i]]+0.3506*band51[[i]]+0.2136*band61[[i]]+0.2678*band71[[i]]
  greenness = -0.4064*band11[[i]]+0.5129*band21[[i]]-0.2744*band31[[i]]-0.2893*band41[[i]]+0.4882*band51[[i]]-0.0036*band61[[i]]-0.4169*band71[[i]]
  wetness = 0.1147*band11[[i]]+0.2489*band21[[i]]+0.2408*band31[[i]]+0.3132*band41[[i]]-0.3122*band51[[i]]-0.6416*band61[[i]]-0.5087*band71[[i]]
  EVI.1 = 2.5*(band21[[i]]*0.0001-band11[[i]]*0.0001)/(1+band21[[i]]*0.0001-band11[[i]]*0.0001)
  TCA1 = atan(greenness/brightness)
  NDWI1 = (band21[[i]]-band61[[i]])/(band21[[i]]+band61[[i]])
  
  #get clear observations [quality 0-2], 
  qc.r = QAfind.mt2(band.qc1[[i]])[[3]] #select 0, 1, 2
  qc.r[qc.r>2] = NA #if 5 bands have clear data, then, use it 
  qc.r[qc.r<=2] = 1
  brightness = brightness*qc.r
  greenness = greenness*qc.r
  wetness = wetness*qc.r
  TCA1 = TCA1*qc.r
  NDWI1 = NDWI1*qc.r
  EVI.1 = EVI.1*qc.r
  
  TCB[[i]] = brightness
  TCG[[i]] = greenness
  TCW[[i]] = wetness
  TCA[[i]] = TCA1
  NDWI[[i]] = NDWI1
  EVI[[i]] = EVI.1
  
  print(paste("Finish calculating TC for ", i, "of", nlayers(band11), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}
TCB = stack(TCB)
TCG = stack(TCG)
TCW = stack(TCW)
TCA = stack(TCA)
NDWI = stack(NDWI)
EVI = stack(EVI)

rm(list = c("brightness","greenness","wetness","TCA1","NDWI1","EVI.1"))

# select dry season for each year, and calculate dry season, median, raster stacks
TCB.dry = list() 
TCG.dry = list()   
TCW.dry = list() 
TCA.dry = list() 
NDWI.dry = list() 
EVI.dry = list() 

for(i in 2013:2015){
  
  idx = which(EarlyDry.idx.M == i) #find the index for each year
  
  TCB.dry1 <- subset(TCB, idx, drop=FALSE)
  TCB.dry[[i-2012]] <- calc(TCB.dry1, median, na.rm = TRUE)
  
  TCG.dry1 <- subset(TCG, idx, drop=FALSE)
  TCG.dry[[i-2012]] <- calc(TCG.dry1, median, na.rm = TRUE)
  
  TCW.dry1 <- subset(TCW, idx, drop=FALSE)
  TCW.dry[[i-2012]] <- calc(TCW.dry1, median, na.rm = TRUE)
  
  TCA.dry1 <- subset(TCA, idx, drop=FALSE)
  TCA.dry[[i-2012]] <- calc(TCA.dry1, median, na.rm = TRUE)
  
  NDWI.dry1 <- subset(NDWI, idx, drop=FALSE)
  NDWI.dry[[i-2012]] <- calc(NDWI.dry1, median, na.rm = TRUE)
  
  EVI.dry1 <- subset(EVI, idx, drop=FALSE)
  EVI.dry[[i-2012]] <- calc(EVI.dry1, median, na.rm = TRUE)
  
  print(paste("Finish Calculating Year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

TCB.dry = stack(TCB.dry)
TCG.dry = stack(TCG.dry)
TCW.dry = stack(TCW.dry)
TCA.dry = stack(TCA.dry)
NDWI.dry = stack(NDWI.dry)
EVI.dry = stack(EVI.dry)

writeRaster(TCB.dry,".\\NBAR_results2\\TCB.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCG.dry,".\\NBAR_results2\\TCG.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCW.dry,".\\NBAR_results2\\TCW.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(TCA.dry,".\\NBAR_results2\\TCA.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(NDWI.dry,".\\NBAR_results2\\NDWI.dry.wa.grd", overwrite=TRUE) #write a raster stack files
writeRaster(EVI.dry,".\\NBAR_results2\\EVI.dry.wa.grd", overwrite=TRUE) #write a raster stack files
