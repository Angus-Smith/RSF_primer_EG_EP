## Angus Smith MSc thesis

#####################################################################
#### [RSF_primer_EG_EP] 02 RSF SSF lesson
#####################################################################

###########
## TO DOs
# - 



##########################
#### Overview
##########################
## Load point data and landscape data. Prep for use in the RSF cours:
## https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html
## 


##########################
#### Settings
##########################

rm(list = ls())



##########################
#### Setup
##########################

setwd("~/Google_Drive/projects/MSc/RSF_primer")

source("RSF_primer_EG_EP/[RSF_primer_EG_EP] functions.R")

pcks <- list("sp", "sf", "raster", "rgdal", "rgeos", "ggmap", "adehabitatHR", "sjPlot")

sapply(pcks, require, char = TRUE)


##########################
#### Load data
##########################

elevation <- raster("data/RSF_layers/elev.tif")
dist.roads <- raster("data/RSF_layers/distroads.tif")
time.since.fire <- raster("data/RSF_layers/tsfire.tif")

time.since.fire[is.na(time.since.fire)] <- 100

dist.roads <- crop(dist.roads, elevation)
time.since.fire <- crop(time.since.fire, elevation)

dist.roads <- resample(dist.roads, elevation, "bilinear")
time.since.fire <- resample(time.since.fire, elevation, "bilinear")

fmch.habitat <- raster::stack(elevation, time.since.fire, dist.roads)

load("data/RSF_layers/elk.rda")
caribou.sp <- elk
caribou.df <- as.data.frame(caribou.sp)
caribou.df$x <- caribou.df$coords.x1
caribou.df$y <- caribou.df$coords.x2
caribou.df$time <- caribou.df$dttm

##########################
#### Lesson
##########################


# Loading habitat RasterStack and plotting it, just to remind us of the beauty
plot(fmch.habitat)

# Plotting the caribou data over elevation
plot(fmch.habitat[[1]])
points(caribou.sp, pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.2))


# RSF
# Compute MCP and plot it
caribou.mcp <- mcp(caribou.sp, 100)
lines(caribou.mcp, col = "red", lwd = 2)

# Generating available (random) points within MCP and plot them
xy.obs <- caribou.df[sample(1:nrow(caribou.df), 200), c("x", "y")]
xy.random <- spsample(caribou.mcp, 1000, "random")@coords
plot(xy.random, asp = 1, col = "darkblue", pch = 19, cex = 0.5)
points(xy.obs, pch = 19, col = "orange", cex = 0.5)

# Creating a data frame of used and available locations and label them as 'TRUE' (used) or 'FALSE' (available)
Data.rsf <- data.frame(X = c(xy.obs[, 1], xy.random[, 1]), 
                       Y = c(xy.obs[, 2], xy.random[, 2]), 
                       Used = c(rep(TRUE, nrow(xy.obs)), rep(FALSE, nrow(xy.random))))

# Extract the habitat variables from the raster stack and add them to this data frame
Data.rsf$elevation <- raster::extract(fmch.habitat[[1]], Data.rsf[, 1:2])
Data.rsf$time.since.fire <- raster::extract(fmch.habitat[[2]], Data.rsf[, 1:2])
Data.rsf$dist.roads <- raster::extract(fmch.habitat[[3]], Data.rsf[, 1:2])

# Boxplots comparing used and available locations with respect to the habitat variables
par(mfrow=c(1,3), bty = "l", tck = 0.01, mgp = c(1.5,.25,0))
boxplot(elevation ~ Used, data = Data.rsf, main = "elevation")
boxplot(time.since.fire ~ Used, data = Data.rsf, main = "time since fire")
boxplot(dist.roads ~ Used, data = Data.rsf, main = "distance to roads")

# Fitting an RSF!!  And looking at coefficients
RSF.fit <- glm(Used ~ scale(elevation) + scale(time.since.fire) + scale(dist.roads), data = Data.rsf, family = "binomial")
summary(RSF.fit)$coef


# Plotting coefficients using 'sjPlot' package
plot_model(RSF.fit)

# Looking at correlation between elevation and distance to roads
with(Data.rsf, plot(elevation, dist.roads, col = factor(Used)))

# A model with interaction between elevation and distance to roads
RSF.fit2 <- glm(Used ~ scale(elevation) * scale(dist.roads) + scale(time.since.fire), data = Data.rsf, family = "binomial")
plot_model(RSF.fit2)

# Comparing AIC values between two models
AIC(RSF.fit, RSF.fit2)

# Adding a squared term to elevation for elevation, fitting the new model
# Then plotting the coefficients and comparing AICs for all three models

Data.rsf$e2 <- Data.rsf$elevation^2
RSF.fit3 <-  glm(Used ~ (scale(elevation) + scale(e2)) * scale(dist.roads) + scale(time.since.fire), data = Data.rsf, family = "binomial")
plot_model(RSF.fit3)
AIC(RSF.fit, RSF.fit2, RSF.fit3)

# Removing the non significant interaction between elevation2 and distance to roads
RSF.fit4 <- update(RSF.fit3, formula = ~ . - scale(e2):scale(dist.roads))
AIC(RSF.fit, RSF.fit2, RSF.fit3, RSF.fit4)


# RSF prediction map

# Make sure we have prediction scopes that match, and creating a new raster called 'fmch.rsf'
# with the same extent and resolution; we will fill this raster with predicted values

fmch.habitat.cropped <- crop(fmch.habitat, extent(caribou.sp))
fmch.rsf <- fmch.habitat.cropped[[1]]

# Extract the approporiate vectors for all habitat layers:
fmch.elevation.vector <- fmch.habitat.cropped[[1]]@data@values
fmch.fire.vector <- fmch.habitat.cropped[[2]]@data@values
fmch.dist.roads.vector <- fmch.habitat.cropped[[3]]@data@values
fmch.elevation2.vector <- (fmch.habitat.cropped[[1]]@data@values)^2
Predict.df <- data.frame(elevation = fmch.elevation.vector, 
                         time.since.fire=fmch.fire.vector, 
                         dist.roads=fmch.dist.roads.vector, 
                         e2=fmch.elevation2.vector)

# use the predict function to fill in the newly created habitat raster with predicted values
fmch.rsf@data@values <- predict(RSF.fit3, newdata = Predict.df)

# Plot the RSF with caribou points on top!
dev.off()
plot(fmch.rsf, main = "RSF", col = grey.colors(100, start = 0, end = 1))
points(caribou.sp, col = alpha("darkblue", 0.3), cex = 0.5)



# SSFs
# Loading more packages
require(amt)
require(sp)

# Picking one animal to run a quick SSF
fa1403<- caribou.df %>%
  filter(ANIMAL_ID == "1010")

# making a 'track_xyt' object required by 'amt' and looking at the sampling rate
track1 <- mk_track(fa1403, x, y, time, crs=CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
summarize_sampling_rate(track1)

# calculating step lengths and plotting a histogram 
steps.1 <- step_lengths(track1)
hist(steps.1, breaks = 30, main="Histogram of step lengths", xlab="Step lengths")

# creating steps from the used (real) points
true.steps1 <- steps(track1)
true.steps1

# creating 15 available (random) steps per used step
true.random.steps1 <- random_steps(true.steps1, n = 15)

# extracting covariates at the endpoints of each step (both used and available)
ssf1 <- true.random.steps1 %>% extract_covariates(fmch.habitat) 







##########################
#### CODE BREAKS HERE
##########################









# running a quick model with each of the covariates, including the step length
ssf.1 <- clogit(case_ ~ strata(step_id_) + scale(elevation) + scale(time.since.fire) + scale(dist.roads) +
                  scale(sl_), method = 'approximate', data = ssf1)

plot_model(ssf.1, title="SSF Coefficients")


# Now doing the same process but with multiple animals rather than just one
nested <- 
  caribou.df %>% 
  nest(-id) %>% 
  mutate(track = map(data, ~ mk_track(., x, y, time, crs= CRS("+proj=utm +zone=6 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))))

names(nested)


# calculating step lengths,  plotting a histogram and looking at sampling rates
step.lengths <- unlist(map(nested$track, step_lengths))
hist(step.lengths, breaks = 20, main="Histogram of step lengths", xlab="Step lengths")
map(nested$track, summarize_sampling_rate)

# creating steps from the used points
true.steps <- map(nested$track, steps)
head(true.steps[[1]])

# creating 10 random (available) steps per used step
true.random.steps <- map(true.steps, random_steps, n=10)

# combining all lists together to get normal data frame with all individuals 
all.steps <- bind_rows(true.random.steps, .id="id")
all.steps

# plotting a histogram of turning angles
with(all.steps, hist(ta_[case_==T], breaks=20, main = "Histogram of turning angles", xlab="turning angle"))

# converting steps to SpatialPointsDataFrame and projecting it
steps.sp <- all.steps
coordinates(steps.sp) <- ~x2_+y2_
proj4string(steps.sp) <- CRS("+proj=utm +zone=6 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Extracting habitat covariates and annotating GPS locations with habitat data
covariates <- raster::extract(fmch.habitat, steps.sp)
steps.sp$elevation <- covariates[,1]
steps.sp$tsf.years <- covariates[,2]
steps.sp$roads <- covariates[,3]
head(steps.sp)

# Looking at used versus available points with respect to elevation
boxplot(elevation ~ case_, data = steps.sp, ylab="elevation")

# Fitting a conditional logistic regression model
clogit(case_ ~ strata(step_id_) + scale(elevation) + scale(tsf.years) +
         scale(roads) + cluster(id) + scale(sl_), method = 'approximate', data = steps.sp)


