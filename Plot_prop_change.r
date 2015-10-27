
#after run calculating all indices and plot trend

lc_rc1 = lc_rc*100 + TCW.trd2.grd 
lc_rc1.df = data.frame(freq(lc_rc1))[-29,]
lc_rc1.df$lc = floor(lc_rc1.df$value/100)
lc_rc1.df$trend = lc_rc1.df$value - lc_rc1.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc1.df, FUN = "sum")
lc_rc1.df$prop = 100*lc_rc1.df$count/rep(lc_sum$count, each = 4)
lc_rc1.df = lc_rc1.df[-which(lc_rc1.df$lc == 1 | lc_rc1.df$lc == 7), ] #remove class 1: water; and 7: other class
lc_rc1.df = lc_rc1.df[-which(lc_rc1.df$trend == 4), ] #not calculated classes

lc_rc2 = lc_rc*100 + EVI.trd2.grd 
lc_rc2.df = data.frame(freq(lc_rc2))[-29,]
lc_rc2.df$lc = floor(lc_rc2.df$value/100)
lc_rc2.df$trend = lc_rc2.df$value - lc_rc2.df$lc*100
lc_sum = aggregate(count~lc, data = lc_rc2.df, FUN = "sum")
lc_rc2.df$prop = 100*lc_rc2.df$count/rep(lc_sum$count, each = 4)
lc_rc2.df = lc_rc2.df[-which(lc_rc2.df$lc == 1 | lc_rc2.df$lc == 7), ] #remove class 1: water; and 7: other class
lc_rc2.df = lc_rc2.df[-which(lc_rc2.df$trend == 4), ] #not calculated classes

lc_rc1.df = rbind(data.frame(lc_rc1.df, VI = rep("TCW", nrow(lc_rc1.df))), 
                  data.frame(lc_rc2.df, VI = rep("EVI", nrow(lc_rc2.df))))

library(reshape2)
library(ggplot2)

lc_rc1.df.long = melt(lc_rc1.df[,c("lc","trend","prop","VI")], id.vars=c("lc", "trend","VI"))

lc_rc1.df.long$trend = factor(lc_rc1.df.long$trend)
levels(lc_rc1.df.long$trend) <- c("Positive", "Negative", "No Trend")
lc_rc1.df.long$lc = factor(lc_rc1.df.long$lc)
levels(lc_rc1.df.long$lc) <- clasnames[2:6]


ggplot(data=lc_rc1.df.long, aes(x=trend, y=value, fill=trend)) +
  facet_grid(VI ~ lc) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("") + ylab("Percentage of Land Cover")
  
ggsave(".\\NBAR_results3\\trend_biomes.png", width = 12, height = 8, units = "in")
