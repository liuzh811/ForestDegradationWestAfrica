##############################################################################
# downlaod MODIS MCD43A4 and MCD43A2 data
#load library
library("MODIS")
library(rgdal)
library(raster)

MODISoptions(localArcPath="D:\\users\\Zhihua\\MODIS",
             outDirPath="D:\\users\\Zhihua\\MODIS",
             gdalPath='c:/OSGeo4W64/bin')

getProduct() # list available products

dates <- as.POSIXct( as.Date(c("1/1/2000","30/7/2015"),format = "%d/%m/%Y") )
dates2 <- transDate(dates[1],dates[2]) # Transform input dates from before
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

#download reflectance
runGdal(product="MCD43A4",  #Nadir BRDF-Adjusted Reflectance, 1000 m reso
        begin=dates2$beginDOY,
        end = dates2$endDOY,
        tileH = 16:18,tileV = 7:8,
        SDSstring = "1", #only extract the first layers
        outProj=proj.geo,
        job = "NBAR")

#download QC
runGdal(product="MCD43A4",  #Nadir BRDF-Adjusted Reflectance, 1000 m reso
        begin=dates2$beginDOY,
        end = dates2$endDOY,
        tileH = 16:18,tileV = 7:8,
        SDSstring = "1", #only extract the first layers
        outProj=proj.geo,
        job = "NBAR_QC")

#################################################################################
# download TRMM 3B42 V7 daily rainfall data
#download website: http://pmm.nasa.gov/data-access/downloads/trmm
#go to: 3B42 Research Derived Daily Product
#click: NetCD/Simple Subset Wizard (GES DISC), go to: http://disc.sci.gsfc.nasa.gov/SSW/#keywords=TRMM_3B42_daily%207
#select geographic region and time period
#need to downland wget: http://www.gnu.org/software/wget/

#download
wget --content-disposition -i 3b42v72015.txt

###################################################################################
## downloading version 6
#write a function to downloa MCD43 v6 data
dlmcd43 <- function(product, #e.g., MCD43A2, MCD43A4
                    version, #e.g., 6
                    start_date, #yyyymmdd, e.g., "20051101"
                    end_date, #yyyymmdd, e.g., "20060401"
                    tileh, #e.g., c("17","18")
                    tilev, #e.g., c("07","08")
                    output_loc) #e.g., "MCD43A2V006"
  {
  
  require(XML)
  library(RCurl)
  
  #construct date first
  d31 = c(paste("0", 1:9, sep = ""), as.character(10:31))
  mon.leap = c(rep("01", 31),rep("02",29), rep("03",31),rep("04",30),rep("05",31),rep("06",30),
               rep("07",31),rep("08",31),rep("09",30),rep("10",31),rep("11",30),rep("12",31))
  day.leap = c(d31, d31[-c(30,31)], d31, d31[-30],d31,d31[-31],d31,d31,d31[-31],d31,d31[-31],d31)
  mod.date = data.frame(year = as.character(c(rep(2000, 366), rep(2001, 365),rep(2002, 365),rep(2003, 365),rep(2004, 366),
                                              rep(2005, 365),rep(2006, 365),rep(2007, 365),rep(2008, 366),
                                              rep(2009, 365),rep(2010, 365),rep(2011, 365),rep(2012, 366),
                                              rep(2013, 365),rep(2014, 365),rep(2015, 365),rep(2016, 366))),
                        month = c(mon.leap, mon.leap[-60],mon.leap[-60],mon.leap[-60],mon.leap,
                                  mon.leap[-60],mon.leap[-60],mon.leap[-60],mon.leap,
                                  mon.leap[-60],mon.leap[-60],mon.leap[-60],mon.leap,
                                  mon.leap[-60],mon.leap[-60],mon.leap[-60],mon.leap),
                        day = c(day.leap, day.leap[-60],day.leap[-60],day.leap[-60],day.leap,
                                day.leap[-60],day.leap[-60],day.leap[-60],day.leap,
                                day.leap[-60],day.leap[-60],day.leap[-60],day.leap,
                                day.leap[-60],day.leap[-60],day.leap[-60],day.leap))
  #get url
  url <- "http://e4ftl01.cr.usgs.gov/MOTA/"
  start_idx = which(mod.date$year == as.character(substr(start_date, 1, 4)) &
                    mod.date$month == as.character(substr(start_date, 5, 6)) &
                    mod.date$day == as.character(substr(start_date, 7, 8)))
 
  end_idx = which(mod.date$year == as.character(substr(end_date, 1, 4)) &
                      mod.date$month == as.character(substr(end_date, 5, 6)) &
                      mod.date$day == as.character(substr(end_date, 7, 8)))

  #create a folder to store hdf data
  dir.create(paste(getwd(),"/",output_loc, sep = ""))
  for (i in start_idx:end_idx){
    url1 <- paste(url, product, ".00",version,"/",
                  paste(mod.date$year[i],mod.date$month[i],mod.date$day[i], sep = "."),
                  "/", sep = "")
    doc <- htmlParse(url1)
    links <- xpathSApply(doc, "//a/@href")
    free(doc)
    for(j in 1:length(tileh)){      
      for(k in 1:length(tilev)){
        #only select h17v07/h17v08
        fn = links[which(substr(links, 18, 23) == paste("h",tileh[j],"v",tilev[k], sep = ""))]
        download.file(paste(url1, fn[1], sep = ""), 
                      destfile = paste(getwd(),"/",output_loc,"/", fn[1], sep = ""),
                      mode = "wb")
        
      } # end of k
    } #end if j
    print(paste("Finish downloading ", i - start_idx, "of", end_idx-start_idx, "of", product, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
    } #end of i
}

setwd("D:/users/Zhihua/MODIS/NBARV006")
dlmcd43(product = "MCD43A4", 
        version = 6,
        start_date="20051101",
        end_date="20060401",
        tileh = "17",
        tilev = c("07","08"),
        output_loc = "MCD43A4")

#"MCD43A2"; data quality layers
dlmcd43(product = "MCD43A2", 
        version = 6,
        start_date="20051101",
        end_date="20060401",
        tileh = "17",
        tilev = c("07","08"),
        output_loc = "MCD43A2")

##############################################################
#change hdf to tif files
#create a MCD43A2tif folder first
require(rgdal)
require(raster)
library(gdalUtils)

gdal_setInstallation(search_path = "C:\\OSGeo4W64\\bin", rescan = TRUE,
                     ignore.full_scan = FALSE, verbose = FALSE)

dir = "D:/users/Zhihua/MODIS/NBARV006"
setwd(dir)

hdf.fn = list.files(path = "./MCD43A2", pattern = "*.hdf$")

#gdalinfo(paste(dir, "/","MCD43A2/", hdf.fn[1], sep = ""))
#gdalsrsinfo(paste(dir, "/","MCD43A2/", hdf.fn[1], sep = ""))
#data layers 12:18 are data quality for band1 - band7

for(i in 1:length(hdf.fn)){
  for (j in 12:18){
gdal_translate(paste("./MCD43A2/", hdf.fn[i],sep = ""),
               paste("./MCD43A2tif/", substr(hdf.fn[i], 1, 24), "quality_band", j-11, ".tif",sep = ""),
               of="GTiff",
               output_Raster=TRUE,verbose=TRUE,sd_index=j)
      }
      print(paste(" HDF to TIF converting for ", i, "of", length(hdf.fn), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
      }

#################################################################
#mosaic files
#create a MCD43A2tif_mosaic folder first

dir = "D:/users/Zhihua/MODIS/NBARV006"
setwd(dir)

band = paste("band",1:7, sep = "")

qa.list1 = list.files(path = "./MCD43A2tif", pattern = "*.tif$")

for (i in 1:length(band)){
qa.list = qa.list1[which(substr(qa.list1, 33, 37) == band[i])]
  
h17v07.b1 = qa.list[which(substr(qa.list, 18,23) == "h17v07")]
h17v08.b1 = qa.list[which(substr(qa.list, 18,23) == "h17v08")]

for (j in 1:length(h17v07.b1)){
tmp1 = paste("./MCD43A2tif/", h17v07.b1[j], sep = "")
tmp2 = paste("./MCD43A2tif/", h17v08.b1[j], sep = "")

mosaic_rasters(gdalfile=c(tmp1,tmp2),
               dst_dataset=paste("./MCD43A2tif_mosaic/", "test_mosaic.tif",sep = ""),
               separate=TRUE,
               of="GTiff",
               verbose=TRUE)

tmp = stack("./MCD43A2tif_mosaic/test_mosaic.tif")
tmp = calc(tmp, sum, na.rm = TRUE)

writeRaster(tmp,paste("./MCD43A2tif_mosaic/", substr(h17v07.b1[j], 1,17), band[i],".tif", sep = ""), format="GTiff", overwrite=TRUE)
file.remove("./MCD43A2tif_mosaic/test_mosaic.tif")

print(paste(" Finish mosaicing for", band[i], "of", length(h17v07.b1), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}
}

############################################################################
#caluclate percent of non-NA for level 2, [0,1,2 were used] 
b1 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band1.tif$"), sep = ""))
b2 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band2.tif$"), sep = ""))
b3 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band3.tif$"), sep = ""))
b4 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band4.tif$"), sep = ""))
b5 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band5.tif$"), sep = ""))
b6 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band6.tif$"), sep = ""))
b7 = stack(paste("./MCD43A2tif_mosaic/", 
                 list.files(path = "./MCD43A2tif_mosaic", pattern = "*band7.tif$"), sep = ""))


#replace 
b1[b1 > 2] = NA
b1[!is.na(b1)] = 1

b2[b2 > 2] = NA
b2[!is.na(b2)] = 1

b3[b3 > 2] = NA
b3[!is.na(b3)] = 1

b4[b4 > 2] = NA
b4[!is.na(b4)] = 1

b5[b5 > 2] = NA
b5[!is.na(b5)] = 1

b6[b6 > 2] = NA
b6[!is.na(b6)] = 1

b7[b7 > 2] = NA
b7[!is.na(b7)] = 1

#calculate sum
b.list = list()
for (i in 1:length(b1)){
  tmp = calc(stack(b1[i], b2[i],b3[i],b4[i],b5[i],b6[i],b7[i]), sum, na.rm = TRUE)
  tmp = tmp >= 5
  b.list[[i]] <- tmp
}

b.list = stack(b.list)
b.list = calc(b.list, sum, na.rm = TRUE)
prop = b.list/length(b1)

