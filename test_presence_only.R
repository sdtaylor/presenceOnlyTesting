library(dplyr)
library(tidyr)
library(doParallel)
library(magrittr)
library(stringr)
library(geosphere) #for regularCoordinates & areaPolygon
library(raster)
library(gbm)
library(Metrics) #for auc

dataFolder='~/data/bbs/'

studyYears=2000:2014

counts=read.csv(paste(dataFolder, 'BBS_counts.csv', sep=''))
routes=read.csv(paste(dataFolder, 'BBS_routes.csv', sep='')) %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID, lat=lati, lon=loni) %>%
  dplyr::arrange(siteID, lon, lat)
species=read.csv(paste(dataFolder, 'BBS_species.csv', sep=''))
weather=read.csv(paste(dataFolder, 'BBS_weather.csv', sep='')) %>%
  mutate(siteID=paste(countrynum, statenum, route,sep='-')) %>%
  dplyr::select(siteID,Year=year, RPID=rpid,runtype)

#Some records are of genus only and "unidentified". Get rid of those.
#Or maybe find someway to incorperate them, because now some sites of false-negs. 
unidSpp=species %>%
  dplyr::select(Aou=AOU, name=english_common_name) %>%
  mutate(unID=ifelse(str_sub(name, 1,5)=='unid.', 1,0)  ) %>%
  filter(unID==1) %>%
  extract2('Aou')

#Filter weather to years of study so it can be used to calculate occData next.
weather=weather %>%
  filter(Year %in% studyYears)

#data frame for data just in the study years and without unidentified spp. 
#and only from runs where the BBS standards =1
occData=counts %>%
  filter(Year %in% studyYears) %>%
  mutate(siteID=paste(countrynum, statenum, Route,sep='-')) %>%
  filter(!Aou %in% unidSpp) %>%
  dplyr::select(Aou, siteID,Year, RPID) %>%
  left_join(weather, by=c('siteID','Year','RPID')) %>%
  filter(runtype==1)  %>%
  dplyr::select(-runtype, -RPID) %>%
  distinct()

#Only keep sites that have been sampled 10 or more years in the past 14.
sitesToKeep=occData %>%
  dplyr::select(siteID, Year) %>%
  distinct() %>%
  group_by(siteID) %>%
  summarize(yearsSampled=n()) %>%
  ungroup() %>%
  filter(yearsSampled >=10) %>%
  extract2('siteID')

occData = occData %>%
  dplyr::select(-Year) %>%
  distinct() %>%
  filter(siteID %in% sitesToKeep)
#  mutate(presence=1) %>%
#  complete(Aou, siteID, fill=list(presence=0))

routes = routes %>%
  filter(siteID %in% sitesToKeep)

################################################################################
#Subset sites into training and testing ah la Harris 2015. 
# All the points within inner.radius of a center point will be in the test set.
# Everything more than outer.radius away will be in the training set.
# radii are in meters. 
north_america_extent=extent(-130, -60, 25, 50)

centers = as.data.frame(regularCoordinates(12)) %>%
  filter(lon >= -130, lon<=-60, lat>=25, lat<=50)

inner.radius = 1.5E5
outer.radius = 3E5

dists=pointDistance(as.matrix(dplyr::select(routes, lon, lat)),as.matrix(centers), lonlat=TRUE, allpairs = TRUE)
dists=as.data.frame(dists)
dists$siteID=routes$siteID
dists=dists %>%
  gather(centerPoint, dist, -siteID)

in_train=dists %>%
  group_by(siteID) %>%
  top_n(1, -dist) %>%
  ungroup() %>%
  filter(dist>=outer.radius) %>%
  extract2('siteID')
in_test=dists %>%
  group_by(siteID) %>%
  top_n(1, -dist) %>%
  filter(dist<=inner.radius) %>%
  extract2('siteID') 

routes=routes %>%
  mutate(type=ifelse(siteID %in% in_train, 'train', ifelse(siteID %in% in_test, 'test','neither'))) %>%
  filter(type != 'neither')

#########################################################################
#Background points for presence only modeling
north_america_extent=extent(-130, -60, 25, 50)
background_points=spsample(as(north_america_extent, 'SpatialPolygons') , 5000, 'random')

#########################################################################
#Get bioclim data for all points
bioclim_files=list.files('./bioclim_cropped/', pattern='*bil', full.names = TRUE)
bioclim_files=bioclim_files[!grepl('xml', bioclim_files)]
bioclim_stack=stack(bioclim_files)

route_locations=routes[,c('lon','lat')]
coordinates(route_locations) = c('lon','lat')
route_data = as.data.frame(extract(bioclim_stack, route_locations)) 
route_data$siteID=routes$siteID
route_data$type=routes$type

route_data = route_data %>%
  filter(!is.na(bio1))

#Flags for background data will always be 
background_data= as.data.frame(extract(bioclim_stack, background_points)) %>%
  filter(!is.na(bio1)) %>%
  mutate(presence=0, type='train')

########################################################################
#Get the area of a species distrubtion by drawing a convex hull over presence points
#75% sure this is calculating as expected. 
get_area=function(p){
  hull=grDevices::chull(p$lon, p$lat)
  hull=c(hull, hull[1])
  hull_coord=p[hull,]
  
  hull_polygon=SpatialPolygons(list(Polygons(list(Polygon(hull_coord)), ID=1)))
  
  #convert area to square km
  geosphere::areaPolygon(hull_polygon)/(1000^2)
}

########################################################################
#Wrapper for model
modelFormula=as.formula('presence ~ bio1+bio2+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio16+bio17+bio18+bio19')

sdm_model=function(trainData, testData){
  model=gbm(modelFormula, n.trees=5000, distribution = 'bernoulli', interaction.depth = 4, shrinkage=0.001, 
            data= trainData)
  perf=gbm.perf(model, plot.it=FALSE)
  x=predict(model, n.trees=perf, newdata=testData, type='response')
  return(x)
}

########################################################################
#Species to use in analysis. hand picked by looking at presence/absence data over north america
spp_list=c(60,1840,3370,3290,3090,2971,2970,2890,2882,4430,4340,4300,3620,7610,7260,7070,6883,6882,6710,7510)

results=data.frame()

for(this_sp in spp_list){
  
  print(paste('Species:', this_sp))
  #Setup datasets for this species
  this_sp_data=occData %>%
    filter(Aou==this_sp) %>%
    mutate(presence=1) %>%
    right_join(route_data, 'siteID') %>%
    mutate(presence=ifelse(is.na(presence), 0, 1))
  
  #Presence/absence training data.
  pa_train=this_sp_data %>%
    filter(type=='train')
  
  #Presence only data with background points used for abscences
  po_train=this_sp_data %>%
    filter(type=='train', presence==1) %>%
    bind_rows(background_data)
  
  #evaluation used for both is true presence/absence from test set holdout
  evaluation=this_sp_data %>%
    filter(type=='test')
  
  #Train p/a model on route only training subset
  pa_predictions=sdm_model(trainData=pa_train, testData=evaluation)
  #evaulate p/a model on route testing subset (presence and absence)
  pa_auc=auc(evaluation$presence, pa_predictions)

  #train p/o model on route only traiing subset (without absences) and all background data
  po_predictions=sdm_model(trainData=po_train, testData=evaluation)
  #evaluate p/o model on route testing subset (presence and absence) 
  po_auc=auc(evaluation$presence, po_predictions)
  
  #Get the area of the species distribution.
  all_presence_sites=this_sp_data %>%
    filter(presence==1) %>%
    extract2('siteID')
  
  #Get lat long of presence only sites to calculate area
  present_sites=this_sp_data %>%
    filter(presence==1) %>%
    extract2('siteID')
  
  present_sites=routes %>%
    filter(siteID %in% present_sites) %>%
    dplyr::select(lon,lat)
  
  sp_area=get_area(present_sites)
  
  results_this_sp=data.frame(Aou=this_sp, po_auc=po_auc, pa_auc=pa_auc, sp_area=sp_area)
  
  results=rbind(results, results_this_sp)
}





