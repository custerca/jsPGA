
#### Import the hydrography file
# Goedatabase retrieved from: https://gisdata.mn.gov/dataset/water-dnr-hydrography
# Retrieved on 28 June 2023

library(sf)
library(sp)
library(data.table)
library(rgeos)
library(rgdal)

file.choose()

# edit for proper file path
mn_lakeinformation <- st_read(dsn = "filepath\\fgdb_water_dnr_hydrography\\water_dnr_hydrography_uncompressed.gdb", 
                              layer = "dnr_hydro_features_all")


setDT(mn_lakeinformation)
mn_lakeinformation[ , .N , is.na(fw_id)] # lakes assigned fisheries and wildife division ID numbers:

#keep only a subset of the data
mn_lakeinformation[ , .N  , wb_class ]
mn_lakeinformation[ wb_class %in% c('Island or Land', 'Lake or Pond', 'Riverine polygon', 'Reservoir'), .N , ] #subset matching the DNR's openwater defined dataset

mn_lakeinformation <- mn_lakeinformation[ wb_class %in% c('Lake or Pond', 'Reservoir'), , ]

#chris wants lat/long
#make table:
# ref:https://sites.google.com/a/lakeheadu.ca/yong-luo/blog/convert-utm-to-longlat-in-r
# prepare UTM coordinates matrix
utmcoor<-SpatialPoints(cbind(mn_lakeinformation$center_utm_x,mn_lakeinformation$center_utm_y), proj4string=CRS("+proj=utm +zone=15")) # UTM zone MN: 15
#utmdata$X and utmdata$Y are corresponding to UTM Easting and Northing, respectively.
#zone= UTM zone
# converting
longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))


unlist(longlatcoor)$coords.x1

mn_lakeinformation[ , `:=` (Lon = unlist(longlatcoor)$coords.x1, Lat = unlist(longlatcoor)$coords.x2 )]


#add NHD ids:
#https://www.sciencebase.gov/catalog/item/6206d3c2d34ec05caca53071
nhds <- fread("filepath\\lake_id_crosswalk.csv") 


#I want to cast this thing wider so that each MN DOW only has one row
molten <- melt(nhds[!is.na(MNDOW_ID), ], id.vars = "MNDOW_ID", na.rm = TRUE)

any(duplicated(molten))

molten <- molten[!duplicated(molten)]

molten[ , case := seq_len(.N) , .(MNDOW_ID,variable) ]

molten[ , variable.v := paste(variable,case, sep = ".")]

nhds.mnkey <- dcast(molten, MNDOW_ID ~ variable.v, value.var = "value")

names(nhds.mnkey)[names(nhds.mnkey)== "site_id.1"] <- "NHD_ID"


#add a col in the NHDs that strips out the colnames
nhds.mnkey[ , MNDOW_ID_c := gsub("MNDOW_", "", MNDOW_ID)  ]

mn_lakeinformation[nhds.mnkey, on = .(dowlknum = MNDOW_ID_c) , NHD_ID := NHD_ID  ]

#check coverage
mn_lakeinformation[ , .N , is.na(NHD_ID) ]

#things look okay?
mn_lakeinformation[ !is.na(NHD_ID) , hist(log(acres)) , ] #dist of sizes

mn_lakeinformation[ , .N , pw_parent_name][order(-N)] #dist of names
mn_lakeinformation[ , .N , pw_parent_name][ , .(nameocc  =  .N) , N ][ order(-N)]

# drop primary geoms
mn_lakeinformation[ ,shape := NULL]

fwrite(mn_lakeinformation[acres >= 1], file = "filepath\\mn_lakeinformation_1acplus.csv")

fwrite(mn_lakeinformation, file = "filepath\\mn_lakeinformation_allopenwat.csv")


############################# Metadata for outputs ########################

# This script generates two files: 
#1) mn_lakeinformation_allopenwat.csv
#2) n_lakeinformation_1acplus.csv

# these files grab the attributes for lake in MN from https://gisdata.mn.gov/dataset/water-dnr-hydrography
# the hydrography dataset is then thinned to only wb_class of 'Lake or Pond' or 'Reservoir'
# next I add a centroid Lat & Long column by converting the existing UTM centroids
# finally I match the dow numbers to a USGS temp data crosswalk from https://www.sciencebase.gov/catalog/item/6206d3c2d34ec05caca53071
# the data are exported in totality (1)
# or limited only to waterbodies 1ac opr greater in size 

################### crosswalk problems ##################


# late in the project we realized that the USGS dat release contained a slightly reduced version of the crosswalk. We think that it's because all of the following lakes are Canadian border waters and maybe the USGS team was hesitant to touch canadian waters (possible that they feared the Canadian Mounted Police)
# the following lakes were Identified as present in a different croswalk that we'd used, so we add them manually here to the publicly available one.


# MNDOW_ID                                      site_id
# <char>                                       <char>
# mndow_16003400  nhdhr_{275873AC-3A1F-4B64-A4F4-C8F9447596D9}
# mndow_16003600  nhdhr_{28921AC1-0FB6-4B3A-90B3-E3B04B90EC4E}
# mndow_16023000                              nhdhr_114539715
# mndow_16032900   nhdhr_6430882C-FC9A-424D-8CFF-3F07909A6C31
# mndow_16033100   nhdhr_DDAEF384-10CE-4109-8346-D0B93ACC11A4
# mndow_16046300   nhdhr_ACFB1F46-E1C4-4257-A1BC-6D0FCE084964
# mndow_16058000   nhdhr_84B4F117-CE2C-420F-AAB0-DC6B72EDE7E5
# mndow_16058100                               nhdhr_80994029
# mndow_16061000   nhdhr_C8B5BE7D-E148-4118-B1B7-D6B0B869021D
# mndow_16061700   nhdhr_8DAC2F2F-D634-4BA6-A20E-49C587FA92A5
# mndow_16063300   nhdhr_C93489F3-69A4-4272-BE3A-E4FD1600AF54
# mndow_38001200   nhdhr_78D51277-53AC-459C-B067-285DEE38CD7A
# mndow_38021100   nhdhr_23508A03-E0C7-4C49-82A9-C5334F287754
# mndow_38052100   nhdhr_346025D4-3DEF-4F72-A1EF-B65F7095ACF1
# mndow_38052300   nhdhr_0B4C4C75-92D7-4173-A828-20954437BBD8
# mndow_38053200   nhdhr_B62A3EBC-3812-446B-8FFC-1BEB9A006D12
# mndow_69060800   nhdhr_26F12E6A-07D1-4DC6-8BEA-980AD41B687A
# mndow_69069300   nhdhr_0E445706-1F3E-47FA-8A83-D38E2C953307
# mndow_69069400   nhdhr_EBD8B5AF-7EBC-43A3-82BA-AF41830C144C
# mndow_69106400   nhdhr_891627FD-D1DF-4F1F-A7F9-907DEDE8239B
# MNDOW_ID                                      site_id
#render these as a table and append to nhds.mnkey

# Creating the data
antiset <- data.table(
  MNDOW_ID = c("mndow_16003400", "mndow_16003600", "mndow_16023000", "mndow_16032900", "mndow_16033100", "mndow_16046300", "mndow_16058000", "mndow_16058100", "mndow_16061000", "mndow_16061700", "mndow_16063300", "mndow_38001200", "mndow_38021100", "mndow_38052100", "mndow_38052300", "mndow_38053200", "mndow_69060800", "mndow_69069300", "mndow_69069400", "mndow_69106400"),
  NHD_ID = c("nhdhr_{275873AC-3A1F-4B64-A4F4-C8F9447596D9}", "nhdhr_{28921AC1-0FB6-4B3A-90B3-E3B04B90EC4E}", "nhdhr_114539715", "nhdhr_6430882C-FC9A-424D-8CFF-3F07909A6C31", "nhdhr_DDAEF384-10CE-4109-8346-D0B93ACC11A4", "nhdhr_ACFB1F46-E1C4-4257-A1BC-6D0FCE084964", "nhdhr_84B4F117-CE2C-420F-AAB0-DC6B72EDE7E5", "nhdhr_80994029", "nhdhr_C8B5BE7D-E148-4118-B1B7-D6B0B869021D", "nhdhr_8DAC2F2F-D634-4BA6-A20E-49C587FA92A5", "nhdhr_C93489F3-69A4-4272-BE3A-E4FD1600AF54", "nhdhr_78D51277-53AC-459C-B067-285DEE38CD7A", "nhdhr_23508A03-E0C7-4C49-82A9-C5334F287754", "nhdhr_346025D4-3DEF-4F72-A1EF-B65F7095ACF1", "nhdhr_0B4C4C75-92D7-4173-A828-20954437BBD8", "nhdhr_B62A3EBC-3812-446B-8FFC-1BEB9A006D12", "nhdhr_26F12E6A-07D1-4DC6-8BEA-980AD41B687A", "nhdhr_0E445706-1F3E-47FA-8A83-D38E2C953307", "nhdhr_EBD8B5AF-7EBC-43A3-82BA-AF41830C144C", "nhdhr_891627FD-D1DF-4F1F-A7F9-907DEDE8239B")
)

mndow_nhdhr_xwalk <- rbind(nhds.mnkey[ , .(MNDOW_ID, NHD_ID) ,],antiset)

write_rds(mndow_nhdhr_xwalk, file = "filepath\\mndow_nhdhr_xwalk.rds")


















