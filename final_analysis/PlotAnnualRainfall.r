setwd("D:\\users\\Zhihua\\TRMM\\TRMM_3B42_daily_v7_2")

library(raster)
library(rgdal)

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

county_b = readOGR(dsn="D:\\users\\Zhihua\\TRMM\\gadm_v2_shp",layer="gadm2_westernafrican_dis")
county_b.geo = spTransform(county_b, CRS(proj.geo))

county_b_ghana = county_b.geo[which(county_b.geo$NAME_0=="Ghana"|county_b.geo$NAME_0=="CÃ´te d'Ivoire"|
                                      county_b.geo$NAME_0=="Sierra Leone"|county_b.geo$NAME_0=="Liberia" 
                                    |county_b.geo$NAME_0=="Guinea"),]

#read into data
fn = list.files(path = ".", pattern = "*.nc$")
r = stack(fn)
r = crop(r, county_b_ghana) #crop the dataset to western Africa boundary
names(r) <- fn

R_year = calc(r, mean)  

plot(R_year*365)
plot(county_b_ghana, add = TRUE)

#plot
