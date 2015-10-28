#define a function to get significant change points
Ex.pts = function(x, sig.level = 0.1){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  dry.evi.sig2 = x <= sig.level
  dry.evi.sig2[dry.evi.sig2 ==0] = NA
  dry.evi.sig2.xy = rasterToPoints(dry.evi.sig2) #get raster coordinate, from left to right, top to bottom
  dry.evi.sig2.xy = data.frame(dry.evi.sig2.xy)
  dry.evi.sig2.sp <- SpatialPoints(coords = cbind(dry.evi.sig2.xy$x,dry.evi.sig2.xy$y),proj4string = CRS(proj.geo))
  return(dry.evi.sig2.sp)
}

#define a function to get all points for a raster
Ex.pts.all = function(x){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

#define a function to get all points for a raster
Ex.pts.all.nonNA = function(x){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  #pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}
#define a function to filter the data
#m is the data quality raster
#return a raster stack, layer 1-4 represent different quality levels
QAfind.mt2 = function(m){ 
  ext = extent(m)
  proj.geo = projection(m)
  
  m = as.matrix(m)
  nc = ncol(m);nr = nrow(m)
  
  bit.start <- c(1,2,seq(5, 32, 4))
  bit.end <- c(1,4,seq(8, 32, 4))
  
  #change value into bit
  qc.bit = apply(m,1,function(x){paste(rev(as.integer(intToBits(rev(x)))), collapse="")})
  bit.start.band1 = bit.start[9]+32*(0:(nc-1))
  bit.end.band1 = bit.end[9]+32*(0:(nc-1))
  
  bit.start.band2 = bit.start[8]+32*(0:(nc-1))
  bit.end.band2 = bit.end[8]+32*(0:(nc-1))
  #band 3
  bit.start.band3 = bit.start[7]+32*(0:(nc-1))
  bit.end.band3 = bit.end[7]+32*(0:(nc-1))
  
  #band 4
  bit.start.band4 = bit.start[6]+32*(0:(nc-1))
  bit.end.band4 = bit.end[6]+32*(0:(nc-1))
  
  #band 5
  bit.start.band5 = bit.start[5]+32*(0:(nc-1))
  bit.end.band5 = bit.end[5]+32*(0:(nc-1))
  
  #band 6
  bit.start.band6 = bit.start[4]+32*(0:(nc-1))
  bit.end.band6 = bit.end[4]+32*(0:(nc-1))
  
  #band 7
  bit.start.band7 = bit.start[3]+32*(0:(nc-1))
  bit.end.band7 = bit.end[3]+32*(0:(nc-1))
  
  #return a list, each element is a points
  qc.bit.band1 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band1, bit.end.band1),2)}) 
  qc.bit.band2 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band2, bit.end.band2),2)})
  qc.bit.band3 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band3, bit.end.band3),2)})
  qc.bit.band4 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band4, bit.end.band4),2)})
  qc.bit.band5 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band5, bit.end.band5),2)})
  qc.bit.band6 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band6, bit.end.band6),2)})
  qc.bit.band7 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start.band7, bit.end.band7),2)})
  
  qc = data.frame(band1 = unlist(qc.bit.band1), 
                  band2 = unlist(qc.bit.band2),
                  band3 = unlist(qc.bit.band3),
                  band4 = unlist(qc.bit.band4),
                  band5 = unlist(qc.bit.band5),
                  band6 = unlist(qc.bit.band6),
                  band7= unlist(qc.bit.band7))
  
  qc0 = apply(qc, 1, function(x){length(which(x > 0))})  #only select quality as 0
  qc1 = apply(qc, 1, function(x){length(which(x > 1))})  #only select quality as 0, 1
  qc2 = apply(qc, 1, function(x){length(which(x > 2))})  #only select quality as 0, 1,2
  qc3 = apply(qc, 1, function(x){length(which(x > 3))})  #only select quality as 0, 1,2
    
  #change back to raster
  q0 = raster(matrix(qc0, nrow = nr, ncol = nc, byrow = TRUE),    
                xmn=ext@xmin, xmx=ext@xmax,
                ymn=ext@ymin, ymx=ext@ymax, 
                crs=CRS(proj.geo))
  
  q1 = raster(matrix(qc1, nrow = nr, ncol = nc, byrow = TRUE),    
              xmn=ext@xmin, xmx=ext@xmax,
              ymn=ext@ymin, ymx=ext@ymax, 
              crs=CRS(proj.geo))
  
  q2 = raster(matrix(qc2, nrow = nr, ncol = nc, byrow = TRUE),    
              xmn=ext@xmin, xmx=ext@xmax,
              ymn=ext@ymin, ymx=ext@ymax, 
              crs=CRS(proj.geo))
  
  q3 = raster(matrix(qc3, nrow = nr, ncol = nc, byrow = TRUE),    
              xmn=ext@xmin, xmx=ext@xmax,
              ymn=ext@ymin, ymx=ext@ymax, 
              crs=CRS(proj.geo))
  
  return(stack(q0,q1,q2,q3))
}

m = raster(matrix(1:8, nrow = 2, byrow = TRUE))
plot(QAfind.mt2(m))


##calculate trends, reture reulst: first = tau, second = sign, third = number of obs
#m is a time seies of VI for one point, 
#if the data point >6, and consective NA values < 8, then calculate trend 
cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the 
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}

mk1 = function(m,Year=c(2001:2015), data.length = 6, na.length = 8) { 
  require(Kendall)
  m1 = m; m1[is.na(m)] = 0
  m.length = length(which(!is.na(m)))
  na.length1 = max(cumul_zeros(m1))
  if (m.length>=data.length & na.length1 < na.length){
    mk = MannKendall(m)
    return(c(mk[[1]][1],mk[[2]][1],m.length))
  } else {return(c(NA,NA,m.length))}
}

# define a function to find the correlationship between vi and rainfall
# for residual analysis, 
# first value is coefficient, 2nd is intercept, 3rd is p value
fun.cor <- function(x) {
  if (length(which(is.na(x[1:15]))) > 8 | length(which(is.na(x[16:30])))>8) 
  {return(c(NA,NA,NA))} 
  else 
  {
    lm1 = lm(as.numeric(x[1:15]) ~ as.numeric(x[16:30]))
    return(c(lm1$coefficients[2],lm1$coefficients[1], anova(lm1)$'Pr(>F)'[1]))
  }
}

#find correlations
COR.test <- function(x) {
  x = as.numeric(x)
  if (length(which(!is.na(x[1:15]))) < 4 | length(which(!is.na(x[16:30]))) < 4) 
  {
    c(NA,NA)
  } else 
  {
    test = cor.test(as.numeric(x[1:15]), as.numeric(x[16:30]),na.action = na.omit)
    c(test$estimate,test$p.value)
  }
}

x = rbind(1:30,1:30,31:60)
apply(x,1,COR.test)

plot(1:15, 16:30)
plot(16:30, 1:15)

#calculate the slope and significant level
fun.slp <- function(x) {
  x1 = 1:length(x)
  lm = lm(x ~ x1)
  #return(c(lm$coefficients[2], summary(lm)$coefficients[,4][2])) #return coef and sig
  return(lm$coefficients[2])
}

#reclassification
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

#reclassification and plot a raster or raster stack
recls2plot = function(x, 
                      threshold = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                      color = "Greens")
{
  require(rasterVis)
  require(raster)
  require(RColorBrewer)
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
  clscolor_change = brewer.pal(length(threshold)+2,color)[c(2:(length(threshold)+2))]
  breaks2_change <- 0:(length(threshold)+1)        
  legendbrks2_change <- 1:(length(threshold)+1) - 0.5
  
  #clasnames_change <- c("< 0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6",">0.6")
  clasnames_change = c(paste("<=", threshold[1], sep = " "))
  for(j in 1:(length(threshold)-1)){
    
    clasnames_change = c(clasnames_change, paste(threshold[j], "-", threshold[j+1], sep = " "))
    
  }
  clasnames_change[length(threshold)+1] <- c(paste(">", threshold[length(threshold)], sep = " "))
  
  levelplot(x.list,
            maxpixels = nrow(x.list)*ncol(x.list),
            at= breaks2_change, margin=FALSE,
            col.regions= clscolor_change,
            colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5))) +  
    latticeExtra::layer(sp.polygons(county_b.geo, col = "black", lwd = 1.5))
  
}

# define a function to change points to raster
Point2raster = function(points, raster){ #points is vector, raster the template
  nr = nrow(raster)
  nc = ncol(raster)
  ext = extent(raster)
  proj.geo = projection(raster)
  r = raster(matrix(points, nrow = nr, ncol = nc, byrow = TRUE),
             xmn=ext@xmin, xmx=ext@xmax,
             ymn=ext@ymin, ymx=ext@ymax, 
             crs=CRS(proj.geo))
  return(r)
}


#plot correlation between rainfall and VI
plot.cor = function(cor,                 #correlation table, first row is correlationship, second row is significant level
                    inraster,           #raster template    
                    threshold = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), #change to categroes    
                    title.txt)               #figure title
  
{
  
  require(rasterVis)
  require(raster)
  require(RColorBrewer)
  
  clscolor_change = brewer.pal(9,"Greens")[c(2:9)]
  breaks2_change <- 0:8        
  legendbrks2_change <- 1:8 - 0.5
  
  clasnames_change <- c("< 0","0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6",">0.6")
  
  tcw.1.df.r = Point2raster(cor[1,], inraster)
  tcw.1.df.sig = Point2raster(cor[2,], inraster)
  
  threshold =  c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  tcw.1.df.r2 = recls2(tcw.1.df.r, threshold = threshold)
  tcw.1.df.sig.sp = Ex.pts(tcw.1.df.sig)
  
  levelplot(tcw.1.df.r2,
            maxpixels = nrow(tcw.1.df.r2)*ncol(tcw.1.df.r2),
            main=title.txt,
            at= breaks2_change, margin=FALSE,
            col.regions= clscolor_change,
            colorkey= list(labels= list(labels= clasnames_change,at= legendbrks2_change, cex = 1.5))) +
    layer(sp.polygons(county_b.geo, col = "black", lwd = 4)) 
  
}

