
TCW1.dry.trd = stack(".\\NBAR_results2\\TCW1.dry.trend.wa.grd")
TCW2.dry.trd = stack(".\\NBAR_results2\\TCW2.dry.trend.wa.grd")
TCW3.dry.trd = stack(".\\NBAR_results2\\TCW3.dry.trend.wa.grd")

TCW1.dry.trd[[1]][TCW1.dry.trd[[2]] > 0.1] = NA
TCW2.dry.trd[[1]][TCW2.dry.trd[[2]] > 0.1] = NA
TCW3.dry.trd[[1]][TCW3.dry.trd[[2]] > 0.1] = NA


trend.vi1.df = raster::extract(TCW1.dry.trd[[1]], pts.sp) #TCW.dry.trd[[1]] contains only the significant change
trend.vi2.df = raster::extract(TCW2.dry.trd[[1]], pts.sp) #TCW.dry.trd[[1]] contains only the significant change
trend.vi3.df = raster::extract(TCW3.dry.trd[[1]], pts.sp) #TCW.dry.trd[[1]] contains only the significant change

lc.modis.df = raster::extract(lc.modis.forest, pts.sp) #based on yearly MODIS land cover data
lc.modis2.df = raster::extract(lc.modis2, pts.sp) #based on global land cover climatology
lc.esa.df = raster::extract(lc.esa, pts.sp)
lc.esa.2000.df = raster::extract(lc.esa.2000, pts.sp)

rain_an.df = raster::extract(f.mean, pts.sp)
rain_slp.df = raster::extract(rf.1mon.slp, pts.sp)

pop_den_change.df = raster::extract(pop_den_dif.GPW4, pts.sp)
pop_den.df = raster::extract(pop_den.afripop, pts.sp)

d2crop.modis.df = raster::extract(lc.modis.crop.dis, pts.sp)
d2crop.modis2.df = raster::extract(lc.modis2.crop.dis, pts.sp)
d2crop.global.ag.df = raster::extract(lc.global.ag.dis, pts.sp)

protect.df = raster::extract(protect.grd, pts.sp)
protect.df[which(protect.df>0)] = 1 #protected
protect.df[is.na(protect.df>0)] = 0 #non-protected

topo.df = raster::extract(stack(dem, slp, tpi5), pts.sp)
colnames(topo.df) <- c("dem","slope","tpi")

d2rd.df = raster::extract(d2rd, pts.sp.alb)
d2set.df = raster::extract(d2set, pts.sp.alb)

data = data.frame(vi1 = trend.vi1.df, 
                  vi2 = trend.vi1.df,
                  vi3 = trend.vi1.df,
                  lc_modis = lc.modis.df, # MODIS yearly land cover
                  lc_modis2 = lc.modis2.df,  # MODIS land cover climatology
                  lc_esa = lc.esa.df,   #ESA land cover 2009
                  lc.esa.2000 = lc.esa.2000.df,  #ESA land cover 2000
                  rain_an = rain_an.df,
                  rain_slp = rain_slp.df,
                  topo.df,
                  pop_den_change = pop_den_change.df,
                  pop_den_2000 = pop_den.df,                  
                  d2rd = d2rd.df,
                  d2set = d2set.df,
                  d2crop_modis = d2crop.modis.df,
                  d2crop_modis2 = d2crop.modis2.df,
                  d2crop_global.ag = d2crop.global.ag.df,
                  protect = protect.df)

data = data[which((data$lc_esa == 40 & data$lc_modis2 == 2) | (data$lc.esa.2000 == 1 & data$lc_modis2 == 2)),] #only select forest pixels

write.csv(data, ".\\NBAR_results2\\rf.data.080715.csv")


################## boosted regression tree analysis ############################
#prepare data, 
# select forest pixels first based on MODIS and ESA GlobalCover 2009
# if both products shows tropical forest, then, it is tropical forest
data2 = data[,c("vi2","lc_modis","d2crop_modis2","rain_an","rain_slp",
                "pop_den_change",
                "dem","slope","tpi",
                      "pop_den_2000","d2rd","d2set","protect")]
colnames(data2)[c(1,2,3)] = c("vi", "lc","d2crop")
data4 = data2[which(data2$lc == 1),] #only select forest pixels

#select ESA data
var.names = colnames(data)
ref.inf = list()
#for (i in 2:3){
  i = 1 #only used 1 
  data2 = data[,c(var.names[i],"lc_modis2","lc_esa","lc.esa.2000","d2crop_modis2","rain_an","rain_slp",
                  #"pop_den_change",
                  #"dem",
                  "slope",
                  "tpi",
                  "pop_den_2000",
                  "d2rd","d2set","protect")]
  
  colnames(data2)[c(1,2,3,4,5)] = c("vi", "lc_modis","lc_esa","lc.esa.2000","d2crop")
  data4 = data2[which((data2$lc_esa == 40 & data2$lc_modis == 2) | (data2$lc.esa.2000 == 1 & data2$lc_modis == 2)),] #only select forest pixels
  
  #remove lc column
  data4 = data4[, !(colnames(data4) %in% c("lc_esa","lc_modis","lc.esa.2000"))]
  
  data4.1 = data4[which(data4$vi < 0),] #select decreasing wetness, degradation
  data4.1$vi = 1
  
  data4.2 = data4[c(which(data4$vi >= 0), is.na(data4$vi)),] #select non degradation cells
  data4.2$vi = 0
  #data4.2 = data4.2[sample(1:nrow(data4.2), nrow(data4.1)), ]
  
  data4 = rbind(data4.1, data4.2)
  table(data4$vi)
  
  data4 = data4[complete.cases(data4),] #remove all na rows
  
  #data4 = data4[-which(data4$vi == "3"),] #remove no trend points
  data4$protect = factor(data4$protect)
  
  data4$d2crop = log10(data4$d2crop+1)
  data4$d2set = log10(data4$d2set+1)
  data4$d2rd = log10(data4$d2rd+1)
  
  #data4 = data4[which(data4$pop_den_change >= 0),]
  data4 = data4[which(data4$pop_den_2000 >= 0),]
  
  #data4$pop_den_change = log10(data4$pop_den_change+1)
  data4$pop_den_2000 = log10(data4$pop_den_2000+1)
  
  #select pixels only > 1200 annual rainfall
  data4 = data4[which(data4$rain_an > 1200), ] 
  #cor(data3$dem, data3$rain_an)
    
  data4.1 = data4[which(data4$vi ==1 ),] #degradation cells
  data4.2 = data4[which(data4$vi ==0 ),] #non-degradation cells

#repeat 10 times, and get the average results

  rel.imp = data.frame(var = colnames(data4)[c(2:ncol(data4))])
  auc = c()
  gbm.list = list()
  gbm.int.list = list()

  N = 0
  
repeat {
  N = N+1
  data3 = rbind(data4.1[sample(1:nrow(data4.1), 2500), ],
                data4.2[sample(1:nrow(data4.2), 2500), ])
  #pairs.panels(data3)
  
  #fit a BRT model
  gbm1=gbm.step(data = data3,                                           # your dataset
                gbm.x = 2:(ncol(data3)),       # predictors index
                gbm.y = 1,             # response variable index
                family = "bernoulli",                             # fit a binominal distribution
                tree.complexity = 3, 
                var.monotone = c(-1,0,1,0,0,1,-1,-1,0), 
                n.trees = 40, 
                learning.rate = 0.01, 
                bag.fraction = 0.5,
                silent = FALSE)
  
  gbm.list[[i]] <- gbm1
  
  find.int <- gbm.interactions(gbm1)
  gbm.int.list[[i]] <- find.int
  #calculate AUC
  preds = predict(gbm1, data3, n.trees = gbm1$gbm.call$best.trees, type="response")
  dev = calc.deviance(obs=data3$vi, pred=preds, calc.mean=TRUE)
  
  d <- cbind(data3$vi, preds)
  pres <- d[d[,1]==1, 2]
  abse <- d[d[,1]==0, 2]
  e <- dismo::evaluate(p=pres, a=abse)
  plot(e,"ROC")
  auc = c(auc, e@auc)
  
  #plot the relative importance of each variables
  rel.imp1 = summary(gbm1,n.trees=gbm1$gbm.call$best.trees,cBars=28,las=1,cex.axis=1.5,cex.lab=1.6,mar=c(5,8,4,2)+0.1,
                     axes=TRUE,axisnames=TRUE, plotit=T)
  rel.imp1 = data.frame(rel.imp1)
  
  rel.imp = merge(rel.imp, rel.imp1, by.x = "var",by.y = "var")
  print(paste("final fitting ", N, "th BRT modelling at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
  if(N>10){break}
}

  rel.imp.mn = data.frame(var = rel.imp[,1], rel.mn = apply(rel.imp[,c(2:11)], 1, mean))
  rel.imp.std = data.frame(var = rel.imp[,1], rel.std = apply(rel.imp[,c(2:11)], 1, sd))
  rel.imp1 = merge(rel.imp.mn, rel.imp.std, by.x = "var",by.y = "var")
  a = rel.imp1[with(rel.imp1, order(-rel.mn)), ]
# write.csv(a, ".\\NBAR_results2\\brt.relative.inf.csv")

# using the following gmb1 to plot partial plot 
# save(gbm1, file = ".\\NBAR_results2\\gbm1.partial.plot.Rdata")
ColumnNames <- colnames(data3)[vars]
Labels1 <- c("D2Crop","MeanAnnulRainfall","Rate_Rainfall_change", 
            #"Pop Den Change",
            #"Elevation",
            "Slope", "TPI",
            "PopDen",
            "D2Rd", "D2Set","Prot")

labels1 <- data.frame(ColumnNames,Labels1)
labels1 <- labels1[match(colnames(data3)[vars], labels1$ColumnNames),]
a = merge(a, labels1, by.x = "var", by.y = "ColumnNames")
a = a[with(a, order(-rel.mn)), ]

  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }

png(file = ".\\NBAR_results2\\relative.inf.png", width = 2000, height = 2000, units = "px", res = 300)
mar=c(3,0,0,0)

barx <- barplot(height=a$rel.mn,horiz=FALSE,xlab="",
                ylim = c(-4.9, 60),
                cex.axis=1.5,
                cex.lab=1.5)
error.bar(barx,a$rel.mn, a$rel.std,length=0.03)
for(i in 1:10){text(0.7+1.2*(i-1),-2.5, labels=a$Labels1[i],cex=1, srt=15)}
mtext("Relative importance (%)", side = 2, cex = 1.5,outer=F, padj = -3,adj = 0.5) #http://stat.ethz.ch/R-manual/R-devel/library/graphics/html/mtext.html

dev.off()

  # partial plot for each variable
  source("D:\\users\\Zhihua\\MODIS\\NBAR_results2\\gbm.plot2.R")
  vars <- 2:(ncol(data3))
  
  ColumnNames <- colnames(data3)[vars]
  Labels <- c("Distance to Cropland","Annual Rainfall","Rate of Rainfall Change", 
              #"Pop Den Change",
              #"Elevation",
              "Slope", "TPI",
              "Population Density",
              "Distance to Road", "Distance to Settlement","Protection Status")
  
  labels <- data.frame(ColumnNames,Labels)
  labels <- labels[match(colnames(data3)[vars], labels$ColumnNames),]
  
  smry1<- summary(gbm1)
  
  labels <- labels[order(match(colnames(data3)[vars],smry1$var)),2]
  labels <- as.vector(labels)

#plot population density on degradations
png(file = ".\\NBAR_results2\\partial.plot.pop_den.png", width = 1500, height = 1500, units = "px", res = 300)
  
mar = c(0.1,0.1,0.1,0.1)
 gbm.plot(gbm1,
         rug=TRUE,
         smooth=TRUE,
         #n.plots=1,
         variable.no = 6,
         x.label = labels[5],
         y.label = "Marginal Effects on Degradation",
         common.scale=TRUE,
         write.title = FALSE,
         show.contrib=FALSE,
         plot.layout=c(1,1),
         cex.axis = 1.3,
         cex.lab=1.5)

dev.off()

#plot distance to cropland on degradations
png(file = ".\\NBAR_results2\\partial.plot.d2crop.png", width = 1500, height = 1500, units = "px", res = 300)

mar = c(0.1,0.1,0.1,0.1)
gbm.plot(gbm1,
         rug=TRUE,
         smooth=TRUE,
         #n.plots=1,
         variable.no = 1,
         x.label = labels[1],
         y.label = "Marginal Effects on Degradation",
         common.scale=TRUE,
         write.title = FALSE,
         show.contrib=FALSE,
         plot.layout=c(1,1),
         cex.axis = 1.3,
         cex.lab=1.5)

dev.off()


  out.fig.names = paste(".\\NBAR_results2\\partial.plot.png", sep = "")
  png(file = out.fig.names, width = 4000, height = 3000, units = "px", res = 300)

  gbm.plot2(gbm1,
            rug=TRUE,
            smooth=TRUE,
            n.plots=9,
            x.label = labels[5],
            y.label = "Marginal Effects on Degradation",
            common.scale=TRUE,
            write.title = FALSE,
            show.contrib=TRUE,
            plot.layout=c(3,3),
            cex.lab=1.5)
  
  dev.off()
  

# using this one for variable.interaction plot
  save(gbm1, find.int, file = ".\\NBAR_results2\\gbm1.variable.interaction.Rdata")
#load(".\\NBAR_results2\\gbm1.variable.interaction.Rdata")
find.int <- gbm.interactions(gbm1)
find.int$interactions
find.int$rank.list

  out.fig.names = paste(".\\NBAR_results2\\variable.interaction.png", sep = "")
  png(file = out.fig.names, width = 2500, height = 2500, units = "px", res = 300)
  
  par(mfrow=c(2,2), mar = c(0.5,0.5,0.1,0.1))
  gbm.perspec(gbm1, find.int$rank.list[1,1], find.int$rank.list[1,3],
              x.label = Labels[find.int$rank.list[1,1]],y.label = Labels[find.int$rank.list[1,3]],z.label = "Degradation Probality",
              #y.range=c(15,20),
              z.range=c(0.4,1), cex.axis = 1.3, cex.lab = 1.3,smooth = "average")
  
gbm.perspec(gbm1, find.int$rank.list[2,1], find.int$rank.list[2,3],
            x.label = Labels[find.int$rank.list[2,1]],y.label = Labels[find.int$rank.list[2,3]],z.label = "Degradation Probality",
            #y.range=c(15,20),
            z.range=c(0.4,1), cex.axis = 1.3, cex.lab = 1.3,smooth = "average")

gbm.perspec(gbm1, find.int$rank.list[3,1], find.int$rank.list[3,3],
            x.label = Labels[find.int$rank.list[3,1]],y.label = Labels[find.int$rank.list[3,3]],z.label = "Degradation Probality",
            #y.range=c(15,20),
            z.range=c(0.4,1), cex.axis = 1.3, cex.lab = 1.3,smooth = "average")

gbm.perspec(gbm1, find.int$rank.list[4,1], find.int$rank.list[4,3],
            x.label = Labels[find.int$rank.list[4,1]],y.label = Labels[find.int$rank.list[4,3]],z.label = "Degradation Probality",
            #y.range=c(15,20),
            z.range=c(0.4,1), cex.axis = 1.3, cex.lab = 1.3,smooth = "average")


  dev.off() 
  
#}
