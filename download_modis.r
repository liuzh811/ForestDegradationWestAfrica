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

