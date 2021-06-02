if (!require("pacman")) install.packages("pacman")

pacman::p_load(cowplot,
               data.table,
               dynatopmodel,
               foreach,
               ggplot2,
               grid,
               gridExtra,
               lightgbm,
               lubridate,
               MASS,
               pbapply,
               plyr,
               purrr,
               raster,
               RColorBrewer,
               readxl,
               psych,
               reshape2,
               rgdal,
               scales,
               sf,
               sp,
               stringr,
               SpaDES,
               stringi,
               tidyr,
               tidyverse,
               ggrepel,
               vegan,
               R.devices)

p_load_gh("hunzikp/velox")

setwd("/mnt/shared/analysis/Andrews_simpleGBM/") # set to project directory
years = "All" # can change to year(s) -- for example, "2013" or "c(2013, 2017)" or "2013:2015"

dir.create("data_processed")
dir.create("output")
dir.create("temp")

start = Sys.time()
source("scripts/utility_functions.R")
source("scripts/001 - raster setup.R")
source("scripts/002 - set up temperature.R")
source("scripts/003 - spatial prediction.R")
Sys.time() - start

