#library(devtools)
#devtools::install_github("sbashevkin/spacetools")

library(tidyverse)
library(deltareportr)
library(lubridate)
library(hms)
library(sf)
library(stars)
library(AICcmodavg)
library(scales)
library(cvTools)
require(patchwork)
require(geofacet)
require(gamm4)
require(dtplyr)
library(deltareportr)

library(rgdal)
library(rmapshaper)
library(broom)
library(rgeos)
library(mgcv)

source("soap_checker/soap_check.R")

# No need to specify absolute file paths in an Rstudio project (if you open it as a project).
data_root<-file.path("data-raw")

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline_UTM10", "Bay_Delta_Poly_Outline_UTM10.shp"))

###################################################
################# Setting Boundary ############
###################################################

Delta.xy.aut <- tidy(Delta.aut)
head(Delta.xy.aut)

deltacoords <- Delta.xy.aut %>% dplyr::select(long,lat,piece)
names(deltacoords) <- c("x", "y", "piece")
borderlist <- split(deltacoords, deltacoords$piece)
names(borderlist)

Delta.xy.aut <- Delta.xy.aut %>% rename(x = long, y = lat)
#Delta.xy.aut$piece

border.aut <- lapply(borderlist, "[", c(1,2))
nr <- seq(1,length(borderlist))

border.aut <- lapply(nr, function(n) as.list.data.frame(border.aut[[n]]))

###################################################
################# Setting knots using 10x10 points ############
###################################################

N <- 20
gx <- seq(min(Delta.xy.aut[,1]), max(Delta.xy.aut[,1]), len = N)
gy <- seq(min(Delta.xy.aut[,2]), max(Delta.xy.aut[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knots <- gp[with(gp, inSide(border.aut, x, y)), ]
row.names(knots)<-1:nrow(knots)


plot(Delta.aut, col="grey")
points(knots, pch=21, bg="orange")
text(knots, labels=rownames(knots))

#Or alternatively you can use the custom knots that Mike created and has been tested to work
#knots <- read.csv(file.path(data_root,"Test_Knots.csv"))
#knots <- knots[,-1]
#plot(Delta.xy.aut, col="grey")
#points(knots, pch=21, bg="orange")
#text(knots, labels=rownames(knots))

####################################################################################################################################
#Load integrated temperature dataset and region shape files-------------------------------
#Same steps from Sam's Temperature QAQC R script

# Load Delta Shapefile from Brian
Delta<-st_read(file.path(data_root,"Delta Subregions"))%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

# Load data
Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USBR", "USGS"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  mutate(Temperature_bottom=if_else(Temperature_bottom>30, NA_real_, Temperature_bottom))%>% #Remove bad bottom temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data betwen 5AM and 8PM
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% #BM: Changed from Sam's original code to make it non-ordered
  mutate(Date_num = as.numeric(Date), # Create numeric version of date for models
         Time = as_hms(Datetime))%>% # Create variable for time-of-day, not date. 
  mutate(Time_num=as.numeric(Time))%>% # Create numeric version of time for models (=seconds since midnight)
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source)%>%
  summarise(N=n())%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_join(Delta) # Add subregions

# Remove any subregions that do not contain at least one of these >50 samples stations from the major monitoring programs
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion) | SubRegion=="Georgiana Slough") # Retain Georgiana Slough because it's surrounded by well-sampled regions
# Visualize sampling regions of major surveys

# Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >50 samples stations from the major monitoring programs
Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)



#Create new data frame as to not affect Sam's original dataset
#Filter just those that have bottom temperature measurements

#Also subset data to 2011 and on-------------------------------------------------------
#As they have the best spatiotemporal coverage

Data_subset<- Data %>% filter(!is.na(Temperature_bottom))
summary(Data_subset$Temperature_bottom)
Data_subset <- Data_subset %>% filter(Year>2010)
str(Data_subset)

#Remove temperatures that seem unreasonable for bottom and surface (those below 5 C and above 30 C)
Data_subset<- Data_subset %>% filter(Temperature_bottom>=5&Temperature_bottom<=30)
Data_subset<- Data_subset %>% filter(Temperature>=5&Temperature<=30)
hist(Data_subset$Temperature_bottom)

#Convert data to UTM since lat and long don't seem to work well with soap-film
cord.dec <- SpatialPoints(cbind(Data_subset$Longitude, -Data_subset$Latitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=10 +datum=WGS84"))
cord.UTM.data.frame <-as.data.frame(cord.UTM)
cord.UTM.data.frame$coords.x2<-abs(cord.UTM.data.frame$coords.x2)

Data_subset$x<-cord.UTM.data.frame$coords.x1
Data_subset$y<-cord.UTM.data.frame$coords.x2

#Keep only data points that are inside the boundaries
Data_subset_inside <- Data_subset[with(Data_subset, inSide(bnd = border.aut, x, y)), ]
#Keep data points outside the boundaries just to see
Data_subset_outside <- Data_subset[!with(Data_subset, inSide(bnd = border.aut, x, y)), ]

#Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_inside$x,Data_subset_inside$y, pch=21, bg="purple")
points(Data_subset_outside$x,Data_subset_outside$y, pch=21, bg="red")


################ Soap Checker ################

soap_check(bnd = border.aut, knots = knots)



# Editing out knots too close to the border with sf ----------------------------------------------------------
Delta.aut_sf<-st_as_sf(Delta.aut)
knots_sf<-st_as_sf(knots, coords=c("x","y"), crs=st_crs(Delta.aut_sf), remove=F)

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_sf, color="red")+
  theme_bw()

distances<-Delta.aut_sf %>%
  st_cast(to = 'LINESTRING') %>% #TUrn polygon into linestring
  st_distance(y = knots_sf) # Get distance of each point from that perimeter linestring

#remove knots within 400 m of boundary
knots_sf_edited<-knots_sf[-which(distances<units::set_units(400, "m")),]

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_sf, color="red")+
  geom_sf(data=knots_sf_edited, color="blue")+
  theme_bw()

knots_edited<-knots_sf_edited%>%
  st_drop_geometry()

############### Fit test model  ######################
m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut)),
              data = Data_subset_inside, family = tw(), method = "REML", knots = knots_edited)

########################################################################################################################################################################
########Test temperature stratification model

# Calculate difference of bottom temperature from surface temperature as response variable
Data_subset_inside$Temperature_difference <- Data_subset_inside$Temperature_bottom-Data_subset_inside$Temperature
hist(Data_subset_inside$Temperature_difference)

Data_subset_inside<- Data_subset_inside %>% filter(Temperature_difference>-5)
hist(Data_subset_inside$Temperature_difference)

#Standardize temperature covariate
Data_subset_inside$Temperature_s<-(Data_subset_inside$Temperature-mean(Data_subset_inside$Temperature))/sd(Data_subset_inside$Temperature)

#Temperature Anomaly Model---------------------------------------------
##Julian day and temperature model (not adjusted for space for now)
temperature_anomaly_GAM<- gam(Temperature ~ s(Julian_day_s,bs="cc",k=5),data=Data_subset_inside)
summary(temperature_anomaly_GAM)
gam.check(temperature_anomaly_GAM)
plot(temperature_anomaly_GAM)

#Add term to the dataset
Data_subset_inside$Temperature_prediction <-predict(temperature_anomaly_GAM,Data_subset_inside)
Data_subset_inside$Temperature_anomaly <-Data_subset_inside$Temperature - Data_subset_inside$Temperature_prediction 

#Run test model with tensor product smooth per previous attempt without soap-film
model_bottom_03_r <- bam(Temperature_difference ~  te(x, y, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("sf", "tp","cc"), k=c(20,7,7),xt = list(list(bnd = border.aut),NULL,NULL))+
                           te(x, y, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("sw", "tp","cc"), k=c(20,7,7),xt = list(list(bnd = border.aut),NULL,NULL)),
                         data = Data_subset_inside, method="fREML", discrete=T, nthreads=4, knots =knots_edited)

gam.check(model_bottom_03_r)
vis.gam(model_bottom_03_r,view=c("x","y"),cond=list(Julian_day_s=0, Temperature_anomaly=-1),plot.type="contour")

# Save model to make predict work
saveRDS(model_bottom_03_r, "model_bottom_03_r.Rds")

############################################################################################################################################
####################ANYTHING BELOW IS EXPLORATORY CODE 

### TO MAKE PREDICT WORK
# 1) Clear environment
# 2) restart R
# 3) Load mgcv
library(mgcv)
# 4) load model
model_bottom_03_r<-readRDS("model_bottom_03_r.Rds")

# 5) Predict
predict(model_bottom_03_r,newdata=Data_subset_inside)

# 6) Load any other packages you need

#Error in if (!setting_geom) { : missing value where TRUE/FALSE needed
str(Data_subset_inside)

#Cross validation test ---------------------------------------------------------------
#with K=5
k <- 5 #the number of folds

set.seed(2020)
folds <- cvFolds(NROW(Data_subset_inside), K=k)

Data_subset_inside$holdoutpred <- rep(0,nrow(Data_subset_inside))

for(i in 1:k){
  train <- Data_subset_inside[folds$subsets[folds$which != i], ] #Set the training set
  validation <- Data_subset_inside[folds$subsets[folds$which == i], ] #Set the validation set
  
  newlm <- bam(Temperature_difference ~  te(x, y, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("sf", "tp","cc"), k=c(20,7,7),xt = list(list(bnd = border.aut,nmax=nmax),NULL,NULL))+
                 te(x, y, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("sw", "tp","cc"), k=c(20,7,7),xt = list(list(bnd = border.aut,nmax=nmax),NULL,NULL)),
               data = train, method="fREML", discrete=T, nthreads=4, knots =knots_edit_v2) #Get your new linear model (just fit on the train data)
  newpred <- predict(newlm,newdata=validation) #Get the predicitons for the validation set (from the model just fit on the train data)
  
  Data_subset_inside[folds$subsets[folds$which == i], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
}

Data_subset_inside$holdoutpred_CorrectSign<-ifelse(Data_subset_inside$Temperature_difference==0,NA,(ifelse(Data_subset_inside$Temperature_difference>0&Data_subset_inside$holdoutpred>0|Data_subset_inside$Temperature_difference<0&Data_subset_inside$holdoutpred<0,1,0)))

summary(Data_subset_inside$holdoutpred_CorrectSign)
