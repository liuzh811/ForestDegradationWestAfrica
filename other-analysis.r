random forest
write.csv(data3, ".\\NBAR_results2\\random_forest_data.csv")

library(randomForest)
rf <- randomForest(vi~., data=data3, na.action=na.omit, ntree = 200,mtry = ncol(data3)-1)
rf
plot(rf) 
hist(treesize(rf))
summary(rf)
names(rf)
rf$err.rate
rf$importance
rf$ntree

library(plotmo)
plotmo(rf, type="prob")

png(file = ".\\NBAR_results2\\rf.png", width = 3000, height = 3000, units = "px", res = 300)

par(mfrow=c(3,2))
varImpPlot(rf)

partialPlot(rf, data3,  x.var = "pop_den_2000",xlab = "Population density in 2000")
partialPlot(rf, data3,  x.var = "d2crop",xlab = "Distance to Cropland")
partialPlot(rf, data3,  x.var = "d2set",xlab = "Distance to Settlement")
partialPlot(rf, data3,  x.var = "rain_an",xlim = c(1300, 3000),xlab = "Annual Rainfall (mm)")
partialPlot(rf, data3,  x.var = "d2rd",xlab = "Distance to Road")

dev.off()

#################### rpart #################################################
library(rpart)
fit = rpart(vi~., 
            #method="class", 
            data=data3)

printcp(fit) # display the results
plotcp(fit) # visualize cross-validation results
summary(fit) # detailed summary of splits

# plot tree
plot(fit, uniform=TRUE, main="Classification Tree for Trend")
text(fit, use.n=TRUE, all=TRUE, cex=0.75)

############### PARTY package ###################################################
##Conditional random forest is better than random forest if you catogeries have multiple levels 
#from http://www.listendata.com/2014/11/random-forest-with-r.html
#https://dinsdalelab.sdsu.edu/metag.stats/code/randomforest.html
library(party)

ct = ctree(vi~., data = data3, 
           controls = ctree_control(mtry = 0,maxdepth = 3)) #no variable random,
plot(ct, main="Conditional Inference Tree")

# random forest using Conditional Inference Tree DO NOT USE THIS, LONG TIME TO WAIT
cf = cforest(vi~., data = data3,
             control = cforest_unbiased(mtry = 0, ntree = 200,maxdepth = 3))
varimp = varimp(cf)

barplot(varimp)
slotNames(cf)
slotNames(cf@ensemble)

#plot a first tree http://stackoverflow.com/questions/19924402/cforest-prints-empty-tree
pt <- party:::prettytree(cf@ensemble[[1]], names(cf@data@get("input")))
nt <- new("BinaryTree")
nt@tree <- pt
nt@data <- cf@data
nt@responses <- cf@responses
nt
plot(nt)

#plot a tree
tr <- party:::prettytree(cf@ensemble[[1]], names(cf@data@get("input")))
plot(new("BinaryTree", tree=tr, data=cf@data, responses=cf@responses))



