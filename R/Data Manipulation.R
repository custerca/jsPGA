library(tidyverse)
library(sf)
library(ncdf4)
library(data.table)

# Setting working directory
# setwd()

# MNDOW (state ID) to NHD ID crosswalk
# Created in the MN_lakes_info.R script
mnxwalk <- readRDS("data/mndow_nhdhr_xwalk.rds")
setDT(mnxwalk)
setnames(mnxwalk,"MNDOW_ID","DOW")
setkey(mnxwalk,DOW)

# Raw MN fish sample data
# Available on ScienceBase DOI:
mnfish <- unique(read.csv("data/MN_Fish_Data_21JUN2023.csv"))
setDT(mnfish)
mnfish$SURVEYDATE <- as.Date(mnfish$SURVEYDATE,format="%m/%d/%Y")
mnfish$SURVEY_ID <- as.factor(mnfish$SURVEY_ID)
mnfish$DOW <- paste0("mndow_",mnfish$DOW)

# Merging MN fish data with NHD ID crosswalk
# summarizing into one sample per species - gear - lake - sample date
MNfish <- data.table::merge.data.table(
                mnfish, 
                mnxwalk, 
                by="DOW")[
                  # filtering out shortjaw cisco, selecting columns
                  COMMON_NAME != "shortjaw cisco",.(DOW,SURVEYDATE,YEAR,EFFORT,COMMON_NAME,TOTAL_CATCH,GEAR_CATEGORY,site_id)
                  # combining different namings of cisco
                ][,COMMON_NAME:=ifelse(COMMON_NAME %in% c("tullibee (cisco)","cisco species"),"cisco",COMMON_NAME)
                  # filtering years
                ][(YEAR <= 2022) & (YEAR >=1998)
                  # only samples between June and September
                ][(month(SURVEYDATE)>=6) & (month(SURVEYDATE)<=9)
                  # summing total catch for each lake - date - gear - effort - species
                ][,.(TOTAL_CATCH=sum(TOTAL_CATCH)),by=.(DOW,site_id,SURVEYDATE,YEAR,EFFORT,COMMON_NAME,GEAR_CATEGORY)
                  # summing catch and effort for each lake - date - gear - species
                ][,.(TOTAL_CATCH=sum(TOTAL_CATCH),EFFORT=sum(EFFORT)),by=.(DOW,site_id,SURVEYDATE,YEAR,COMMON_NAME,GEAR_CATEGORY)]

setnames(MNfish,"GEAR_CATEGORY","GEAR")

# Creating data table to include column for sample year and previous 4 years to create 5 year rolling means for lake temperature and secchi (water clarity)
MNyears <- unique(MNfish[,.(DOW,site_id,SURVEYDATE,YEAR)])[,.(YEAR5=seq(YEAR-4,YEAR,by=1)),by=.(DOW,site_id,SURVEYDATE)]

# reading in NLDAS lake temperature data, only retaining columns of interest
# Corson-Dosch, H.R., Mcaliley, W.A., Platt, L.R.C., Padilla, J.A., and Read, J.S., 2023, 
# Daily water column temperature predictions for thousands of Midwest U.S. lakes between 1979-2022 
# and under future climate scenarios: U.S. Geological Survey data release, 
# https://doi.org/10.5066/P9EQQER7.

laketemp <- arrow::read_feather("Data/lake_temperature_metrics_GLM_NLDAS.feather",
                                col_select = c("site_id","year","mean_surf_jul"))
setDT(laketemp)
# renaming column for merging
setnames(laketemp,"year","YEAR5")

# merging 5-year mean DT with temperature data
MNyear_temp <- data.table::merge.data.table(MNyears,
                                            laketemp,
                                            all.x=TRUE,
                                            by=c("site_id","YEAR5"))

# checking to make sure every lake has all 5 years of data
mnchkr1 <- MNyear_temp[,.(chkr=sum(!is.na(mean_surf_jul))),by=c("site_id","DOW","SURVEYDATE")]
# count(mnchkr1,chkr)
# # some only have 0 or 4

# Cleanest data.  All lakes/surveys with missing data along the way (i.e., no site ID or surface temp) were filtered out
MNfish_temp <- data.table::merge.data.table(MNfish,
                                            # merging with a filtered DT to only lakes with all 5 years of lake temperatures
                                            MNyear_temp[site_id %in% mnchkr1[chkr==5,site_id] &
                                                          DOW %in% mnchkr1[chkr==5,DOW],
                                                        # Mean temperature across all 5 years
                                                        .(july5yr=mean(mean_surf_jul)),
                                                        by=.(site_id,DOW,SURVEYDATE)],
                                            by=c("DOW","site_id","SURVEYDATE"))[
                                              # Removing NAs
                                              !is.na(july5yr)]

# Water clarity data
# Max Gilnes, Rensselaer Polytechnic Institute,Troy, NY, United States, 08/2023, written communication)

secchi <- read.csv("data/annual_median_remote_sensing_secchi.csv")
setDT(secchi)
setnames(secchi,"year","YEAR5")
setnames(secchi,"annual.median.rs","secchi")

# Merging secchi with 5 year DT
MNyear_secchi <- data.table::merge.data.table(MNyears,
                                              unique(secchi[,.(site_id,YEAR5,secchi)]),
                                              by=c("site_id","YEAR5"))

# checking to make sure every lake has all 5 years of data
mnchkr2 <- MNyear_secchi[,.(chkr=sum(!is.na(secchi))),by=c("site_id","DOW","SURVEYDATE")]
# count(mnchkr2,chkr)
# # some only have 0 or 3 or 4

MNfish_secchi <- data.table::merge.data.table(MNfish_temp,
                                              # merging with a filtered DT to only lakes with all 5 years of water clarity
                                              MNyear_secchi[site_id %in% mnchkr2[chkr==5,site_id] &
                                                              DOW %in% mnchkr2[chkr==5,DOW],
                                                            # Mean secchi across all 5 years
                                                            .(secchi5=mean(secchi)),
                                                            by=.(site_id,DOW,SURVEYDATE)],
                                              by=c("DOW","site_id","SURVEYDATE"))[
                                                # Removing NAs
                                                !is.na(secchi5)]

# All LAGOS data obtained from 
# LAGOS-US Geo https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1136.3
# LAGOS-US Locus https://portal.edirepository.org/nis/mapbrowse?packageid=edi.854.1
# LAGOS land use data
lagos_lc <- read.csv("Data/lagos geo/zone_landuse.csv")
setDT(lagos_lc)
# watershed spatial division
lagos_ws <- lagos_lc[spatial_division=="ws"]
rm(lagos_lc)

lagos_info <- read.csv("Data/lagos locus/lake_information.csv")
setDT(lagos_info)

lagos_nhd <- lagos_info[,site_id:=paste0("nhdhr_",lake_nhdid)
                        # filtering to sites within MN fish data
                        ][site_id %in% MNfish_temp$site_id
                          # selecting columns and mutating data types
                          ][,.(lagoslakeid=as.character(lagoslakeid),
                               zoneid=as.character(ws_zoneid),
                               lake_elevation_m,
                               lon=lake_lon_decdeg,
                               lat=lake_lat_decdeg),
                            by=site_id]

lagos_char <- read.csv("Data/lagos locus/lake_characteristics.csv")
setDT(lagos_char)

# selecting columns and mutating data types
lagos_char_2 <- lagos_char[,.(lagoslakeid=as.character(lagoslakeid),lake_waterarea_ha)]

lagos_nhd_2 <- data.table::merge.data.table(lagos_nhd,
                                            lagos_char_2,
                                            by="lagoslakeid")


# merging all LAGOS data
lagos_landuse <- merge(lagos_nhd_2,
                       lagos_ws,
                       by="zoneid")[,
                                    .(lon,lat,
                                      elevation=lake_elevation_m,
                                      lakearea=lake_waterarea_ha,
                                      total_dev = sum(c(nlcd_devopen21_pct, nlcd_devlow22_pct,
                                                        nlcd_devmed23_pct, nlcd_devhi24_pct), na.rm=T),
                                      total_ag = sum(c(nlcd_past81_pct, nlcd_cultcrop82_pct), na.rm=T),
                                      total_for = sum(nlcd_fordec41_pct,nlcd_forcon42_pct,nlcd_formix43_pct,na.rm=T)),
                                    by=c("site_id","year")]

# function to find closest matching NLCD year
nlcdyr <- function(x){
  lcyear <- unique(lagos_landuse$year)
  y <- lcyear[min(which(abs(lcyear - x) == min(abs(lcyear - x))))]
  
  return(y)
}

# adding a (temporarily poorly) named column for closest NLCD year
MNfish_secchi_nlcd <- MNfish_secchi[,year := nlcdyr(YEAR),by=seq_len(nrow(MNfish_secchi))]

MNfish_full <- data.table::merge.data.table(MNfish_secchi_nlcd,
                                            lagos_landuse,
                                            by=c("site_id","year"))

setnames(MNfish_full,"year","NLCDyear")
saveRDS(MNfish_full,file="data/MNfishfull.rds")

# selecting environmental covariates then transforming and standardizing them
zmat <- unique(MNfish_full[,.(DOW,site_id,SURVEYDATE,
                              secchi5,elevation,lakearea,
                              total_dev,total_ag,total_for)])[,.(DOW,site_id,SURVEYDATE,
                                                                 secchi.z=as.numeric(scale(secchi5)),
                                                                 elevation.z=as.numeric(scale(elevation)),
                                                                 lakearea.z=as.numeric(scale(log(lakearea))),
                                                                 total.dev.z=as.numeric(scale(car::logit(total_dev/100))),
                                                                 total.ag.z=as.numeric(scale(car::logit(total_ag/100))),
                                                                 total.for.z=as.numeric(scale(car::logit(total_for/100))))]

MNfish_z <- data.table::merge.data.table(MNfish_full[,.(DOW,site_id,SURVEYDATE,YEAR,
                                                        COMMON_NAME,GEAR,TOTAL_CATCH,EFFORT,
                                                        july5yr,lon,lat)],
                                         zmat,
                                         by=c("DOW","site_id","SURVEYDATE"))
  
saveRDS(MNfish_z,file="data/MNfishz.rds")

# this data set is eventually used for future predictions where we treat land cover as most recently available
# merging all LAGOS data and selecting most recent year (2016) of land use data
lagos_landuse_2016 <- data.table::merge.data.table(lagos_nhd_2,
                                                   lagos_ws,
                                                   by="zoneid")[year==2016,
                                                                .(lon,lat,
                                                                  elevation=lake_elevation_m,
                                                                  lakearea=lake_waterarea_ha,
                                                                  total_dev = sum(c(nlcd_devopen21_pct, nlcd_devlow22_pct,
                                                                                    nlcd_devmed23_pct, nlcd_devhi24_pct), na.rm=T),
                                                                  total_ag = sum(c(nlcd_past81_pct, nlcd_cultcrop82_pct), na.rm=T),
                                                                  total_for = sum(nlcd_fordec41_pct,nlcd_forcon42_pct,nlcd_formix43_pct,na.rm=T)),
                                                                by=site_id]

MNfish_full_2016 <- data.table::merge.data.table(MNfish_secchi,
                                                 lagos_landuse_2016,
                                                 by="site_id")

saveRDS(MNfish_full_2016,file="data/MNfishfull_2016.rds")

# selecting environmental covariates then transforming and standardizing them
zmat_2016 <- unique(MNfish_full_2016[,.(DOW,site_id,SURVEYDATE,
                                      secchi5,elevation,lakearea,
                                      total_dev,total_ag,total_for)])[,.(DOW,site_id,SURVEYDATE,
                                                                         secchi.z=as.numeric(scale(secchi5)),
                                                                         elevation.z=as.numeric(scale(elevation)),
                                                                         lakearea.z=as.numeric(scale(log(lakearea))),
                                                                         total.dev.z=as.numeric(scale(car::logit(total_dev/100))),
                                                                         total.ag.z=as.numeric(scale(car::logit(total_ag/100))),
                                                                         total.for.z=as.numeric(scale(car::logit(total_for/100))))]


MNfish_z_2016 <- data.table::merge.data.table(MNfish_full_2016[,.(DOW,site_id,SURVEYDATE,YEAR,
                                                             COMMON_NAME,GEAR,TOTAL_CATCH,EFFORT,
                                                             july5yr,lon,lat)],
                                         zmat_2016,
                                         by=c("DOW","site_id","SURVEYDATE"))

saveRDS(MNfish_z_2016,file="data/MNfishz_2016.rds")

############################################  Data manipulation for predictions  #############################################

library(rstan)
library(tidyverse)
library(data.table)

#setwd()



TRC = function(temp, CTmax, Topt, sigma){
  trc = case_when(
    temp <= Topt ~ exp(-((temp - Topt)/(2*sigma))^2), 
    CTmax >= temp & temp > Topt ~ 1 - ((temp - Topt)/(Topt - CTmax))^2,
    temp > CTmax ~ 0)
  return(trc)
}

# Full MN lake data
# Created in the MN_lakes_info.R script
MNlakeraw <- read.csv("data/mn_lakeinformation_allopenwat/mn_lakeinformation_1acplus.csv")
MNlakes <- MNlakeraw %>%
  filter(!is.na(dowlknum)) %>%
  mutate(DOW=paste0("mndow_",dowlknum)) %>%
  select(DOW,site_id=NHD_ID,lon=Lon,lat=Lat) %>%
  distinct()
rm(MNlakeraw)
MNlakes$site_id[MNlakes$site_id == ""] <- NA

# Correcting an incorrect coordinate
# RAW - lat: 34.55589, lon: -43.99598
# Replacing with coordinates from fishfull (from LAGOS)
# unique(fishfull[fishfull$DOW=="mndow_11030500"]$lon)
# unique(fishfull[fishfull$DOW=="mndow_11030500"]$lat)
MNlakes[MNlakes$DOW=="mndow_11030500",]$lon <- -94.35107
MNlakes[MNlakes$DOW=="mndow_11030500",]$lat <- 46.44632

setDT(MNlakes)

dow_order <- MNlakes %>%
  arrange(DOW) %>%
  distinct(DOW, .keep_all = T) %>%
  select(DOW,site_id)


# secchi (just using the most recent year for future predictions)
# Max Gilnes, Rensselaer Polytechnic Institute,Troy, NY, United States, 08/2023, written communication)
secchi <- read.csv("data/annual_median_remote_sensing_secchi.csv")
setDT(secchi)
setnames(secchi,"annual.median.rs","secchi")

secchi_rec <- secchi[,.(my=max(year),secchi,year),by=site_id][year==my]
secchi_rec$my <- NULL
secchi_rec$year <- NULL

MN_work <- data.table::merge.data.table(dow_order,secchi_rec,by="site_id") %>%
  distinct()

# LAGOS land use data
# All LAGOS data obtained from 
# LAGOS-US Geo https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1136.3
# LAGOS-US Locus https://portal.edirepository.org/nis/mapbrowse?packageid=edi.854.1
lagos_lc <- read.csv("Data/lagos geo/zone_landuse.csv")
setDT(lagos_lc)
lagos_ws <- lagos_lc[spatial_division=="ws"]
rm(lagos_lc)

lagos_info <- read.csv("Data/lagos locus/lake_information.csv")
setDT(lagos_info)

lagos_nhd <- lagos_info[,site_id:=paste0("nhdhr_",lake_nhdid)
][site_id %in% MN_work$site_id
][,.(lagoslakeid=as.character(lagoslakeid),
     zoneid=as.character(ws_zoneid),
     lake_elevation_m),
  by=site_id]

lagos_char <- read.csv("Data/lagos locus/lake_characteristics.csv")
setDT(lagos_char)

lagos_char_2 <- lagos_char[,.(lagoslakeid=as.character(lagoslakeid),
                              lake_waterarea_ha)]

lagos_nhd_2 <- data.table::merge.data.table(lagos_nhd,
                                            lagos_char_2,
                                            by="lagoslakeid")

lagos_landuse <- data.table::merge.data.table(lagos_nhd_2,
                                              lagos_ws,
                                              by="zoneid")[year==2016,
                                                           .(elevation=lake_elevation_m,
                                                             lakearea=lake_waterarea_ha,
                                                             total_dev = sum(c(nlcd_devopen21_pct, nlcd_devlow22_pct,
                                                                               nlcd_devmed23_pct, nlcd_devhi24_pct), na.rm=T),
                                                             total_ag = sum(c(nlcd_past81_pct, nlcd_cultcrop82_pct), na.rm=T),
                                                             total_for = sum(nlcd_fordec41_pct,nlcd_forcon42_pct,nlcd_formix43_pct,na.rm=T)),
                                                           by=site_id]

MN_work <- data.table::merge.data.table(MN_work,lagos_landuse,by="site_id")

# NLDAS Lake temperature
# Corson-Dosch, H.R., Mcaliley, W.A., Platt, L.R.C., Padilla, J.A., and Read, J.S., 2023, 
# Daily water column temperature predictions for thousands of Midwest U.S. lakes between 1979-2022 
# and under future climate scenarios: U.S. Geological Survey data release, 
# https://doi.org/10.5066/P9EQQER7.
laketemp <- arrow::read_feather("Data/lake_temperature_metrics_GLM_NLDAS.feather",
                                col_select = c("site_id","year","mean_surf_jul")) %>%
  # Grabbing year 2000 to match "historical" GCM and year 2021 as "current" temperature data
  filter(year %in% c(2000,2021)) %>%
  pivot_wider(names_from = "year",values_from = "mean_surf_jul",names_prefix = "NLDAS_")

MN_work <- data.table::merge.data.table(MN_work,laketemp,by="site_id")

# GCM temperature data
# Corson-Dosch, H.R., Mcaliley, W.A., Platt, L.R.C., Padilla, J.A., and Read, J.S., 2023, 
# Daily water column temperature predictions for thousands of Midwest U.S. lakes between 1979-2022 
# and under future climate scenarios: U.S. Geological Survey data release, 
# https://doi.org/10.5066/P9EQQER7.
temp_gcm <- read.csv("data/gcm_july_temps_by_lake_time_period.csv")

# Difference between historical temps and future temps
gcm_diff <- temp_gcm %>%
  select(site_id,GCM,time.period,mean.july.surf.temp) %>%
  pivot_wider(names_from = "time.period",values_from = "mean.july.surf.temp") %>%
  mutate(LC=`late-century`-`historic`,
         MC=`mid-century`-`historic`) %>%
  select(site_id,GCM,LC,MC)

# MN_work <- data.table::merge.data.table(MN_work,gcm_mean,by="site_id")
MN_work <- data.table::merge.data.table(MN_work,gcm_diff,by="site_id")


# Full dataset for predictions
MN_full <- left_join(dow_order,MN_work,by=c("DOW","site_id")) %>%
  select("DOW","site_id","secchi","elevation","lakearea","total_dev","total_ag","total_for",
         "NLDAS_2000","NLDAS_2021","GCM","LC","MC",starts_with("comp"))

MN_full <- MN_full[complete.cases(MN_full),]

##### create spatial basis functions using svd

locs <- MNlakes %>%
  filter(DOW %in% MN_full$DOW) %>%
  arrange(DOW) %>%
  distinct(DOW, .keep_all = T) %>%
  select(lon, lat) %>%
  as.matrix()

dmat <- fields::rdist.earth(locs, miles = F)
diag(dmat) <- 0

C = fields::Matern(dmat,range = 100,smoothness = 2.5)

dcomp <- svd(C)
saveRDS(dcomp,"data/dcomp_fullMN.rds")

dow_order <- MNlakes %>%
  filter(DOW %in% MN_full$DOW) %>%
  arrange(DOW) %>%
  distinct(DOW, .keep_all = T) %>%
  select(DOW,site_id)

decomp_df <- as_tibble(dcomp$u[,1:16] %*% diag(sqrt(dcomp$d)[1:16]), .name_repair = ~paste0("comp", 1:16)) %>%
  mutate(lon = locs[,1],
         lat = locs[,2],
         DOW = dow_order$DOW)

MN_full <- left_join(MN_full,decomp_df,by="DOW")

MN_sbv <- select(MN_full,DOW,starts_with("comp"))
saveRDS(MN_sbv,"data/MNsbv_M16.rds")

# Standardizing with original mean and scale used for model fit
og_fishdat <- readRDS("data/MNfishfull.rds")
Zsecchi=scale(og_fishdat$secchi5)
Zelevation=scale(og_fishdat$elevation)
Zlakearea=scale(log(og_fishdat$lakearea))
Ztotaldev=scale(car::logit(og_fishdat$total_dev/100))
Ztotalag=scale(car::logit(og_fishdat$total_ag/100))

setDT(MN_full)
MN_fullz <- MN_full[,`:=`(secchi.z=as.numeric(scale(secchi,attr(Zsecchi, "scaled:center"), attr(Zsecchi, "scaled:scale"))),
                          elevation.z=as.numeric(scale(elevation,attr(Zelevation, "scaled:center"), attr(Zelevation, "scaled:scale"))),
                          lakearea.z=as.numeric(scale(log(lakearea),attr(Zlakearea, "scaled:center"), attr(Zlakearea, "scaled:scale"))),
                          total.dev.z=as.numeric(scale(car::logit(total_dev/100), attr(Ztotaldev, "scaled:center"), attr(Ztotaldev, "scaled:scale"))),
                          total.ag.z=as.numeric(scale(car::logit(total_ag/100),attr(Ztotalag, "scaled:center"), attr(Ztotalag, "scaled:scale")))),
                    by=c("site_id","DOW")]

MN_fullcc <- MN_fullz[complete.cases(MN_fullz),]
saveRDS(MN_fullcc,file="data/MNfullcc.rds")




