# 1. Setting up the R environment -----------------------------------------

# install and load gen3sis
require(gen3sis)
# get package version
print(paste("gen3sis version:", packageVersion("gen3sis")))

# package here robust simple path management
require(here)
# we will also use terra for landscape visualization
require(terra)
# package ape used for plotting phylogenies
require(ape)

# set path to directory containing the landscape
landscape_dir <- here("space")
# look at folder structure
list.files(landscape_dir)

# set path to config_file
config_file <- here("config/config_simple.R")


# 2. Understanding the basics ---------------------------------------------


# load landscapes
lc <- readRDS(file.path(landscape_dir,"landscapes.rds"))
# class(lc) "list"
# get names of landscape variables
names(lc)

# get first time step
first_step_pos <- ncol(lc$mean_temp)
first_step <- colnames(lc$mean_temp)[first_step_pos]
# get first 10 sites of mean temperature for the 2 last time steps and the first (oldest) time step
df_mean_temp <- lc$mean_temp
lc$mean_temp[100:110, c(1:4, first_step_pos)]

# plot mean_temp for first and last time step
par(mfrow = c(1,2))
plot(r_first <- rast(lc$mean_temp[ ,c("x", "y", first_step)], type = "xyz"), main = "first")
plot(r_last <- rast(lc$mean_temp[ ,c("x", "y", "0")], type = "xyz"), main = "last")
par(mfrow = c(1,1)) # set it back

## juliana plot ##
df_mean_temp_selected <- lc$mean_temp[, c(1:3, first_step_pos)]
# or
df_mean_temp_selected <- lc$mean_temp[, c("x", "y", "0", first_step)]

raster_obj <- rast(df_mean_temp_selected)

# plot
plot(raster_obj)

## end ##

# Create a way to visualize the last 100 time steps of this gen3sis input. If possible think of abstracting for any x,y,z1,z2,z3… temporal data-frame. If possible, try to make a function.

# solution 1
for (ti in 100:0) {
  ri <- rast(lc$elevation[,c("x", "y", as.character(ti))], type="xyz")
  plot(ri, main=paste(ti/100, "Ma")) # divide by 100 to get Ma since 1 time-step =10 kyr
  Sys.sleep(0.1)
}

# solution 2
# define animation function
plot_landenv <-  function(df, times=100:0, reverse=F, speed=0.1){
  # df is a data frame with x, y coordinates and time steps as columns with environmental variables
  # reverse is a boolean. Reverse the order of time steps?
  # times is either a vector of time steps or a character "all"
  # speed in seconds
  #df <- lc$temp
  if (times[1] == "all"){
    times <- names(df)[!names(df)%in%c("x", "y")]
    # or
    # times <- setdiff(names(df), c("x", "y"))

  }
  if (reverse){
    times <- rev(times)
  }
  for (ti in times){
    ri <- terra::rast(df[,c("x", "y", ti)], type="xyz")
    plot(ri, main = paste(as.numeric(ti)/100, "Ma"))
    Sys.sleep(speed)
  }
}
# call animation function
plot_landenv(lc$elevation, times=c(100:0), speed=0.1)
# example of plotting entire time series
# plot_landenv(lc$min_temp, times="all", reverse=T, speed=0.03)

# load config
cf <- create_input_config(config_file = here("config/config_simple.R"))
# list all main elements of the config file
names(cf$gen3sis)

# list all elements of the general section, i.e. the main settings and not so much on the eco-evolutionary processes
names(cf$gen3sis$general)

# 3. Run a basic simulation -----------------------------------------------

# run simulation
sim <- run_simulation(config = here("config/config_simple.R"), 
                      landscape = here("space"), 
                      output_directory = here("output"))

sim <- readRDS(here("output/config_simple/sgen3sis.rds"))

#check elements inside the sim object
names(sim)

# visualize the outputs
plot_summary(sim)

# plot richness from summary in custom fashion
na_mask <- is.na(lc$elevation[,"0"])
rich <- sim$summary$`richness-final`
rich[na_mask,3] <- NA
plot(rast(rich, type="xyz"), col = c("grey", gen3sis::color_richness_non_CVDCBP(max(rich, na.rm=T))), main = "Richness")

# plot richness at time step 32 using saved data
sps32 <- readRDS(here("output/config_simple/species/species_t_32.rds"))
lc32 <- readRDS(here("output/config_simple/landscapes/landscape_t_32.rds"))
plot_richness(sps32, lc32)

phy <- read.nexus(file.path(here("output/config_simple/phy.nex")))
plot(phy)

# 4. Customize simulations ------------------------------------------------

conf <- create_input_config(here("config/config_M2_TH.R"))
conf$gen3sis$speciation$divergence_threshold = 30

conf$gen3sis$general$end_of_timestep_observer = function(data, vars, config){
  plot_richness(data$all_species, data$landscape)
  # make p/a matrices
  out_dir <- config$directories$output
  if(!file.exists(file.path(out_dir,"occs"))){
    dir.create(file.path(out_dir, "occs"))
  }
  # cell names
  all_cells <- rownames(data$landscape$coordinates)
  # get 0 for absence and 1 for presence in each grid cell
  asp <- do.call(cbind,
                 lapply(data$all_species, FUN = function(x) {
                   ifelse(all_cells %in% names(x$abundance), 1, 0)
                 }))
  # colnames are species names
  colnames(asp ) <- unlist(lapply(data$all_species, function(x){x$id}))
  # column bind with x/y coordinates
  presence_absence_matrix <- cbind(data$landscape$coordinates, asp)
  saveRDS(presence_absence_matrix, 
          file = file.path(out_dir,"occs", paste0("pa_t_",vars$ti, ".rds")))
}

simod <- run_simulation(config = conf, 
                        landscape = here("space"), 
                        output_directory=tempdir())
plot_summary(simod)

# 5. Troubleshoot ---------------------------------------------------------

# Use 'stop_time' to halt execution at a specific timestep in the landscape object:
stop_time <- "62"

get_dispersal_values <- function(n, species, landscape, config) {
  if (landscape$timestep == stop_time) {
    browser()
  }
  
  # You can also check the 'vars' object for the current timestep:
  vars <- dynGet("vars", inherits = TRUE)
  if (vars$ti == stop_time) {
    browser()
  }
}
