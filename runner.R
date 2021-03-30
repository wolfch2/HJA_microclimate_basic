require(cowplot)
require(data.table)
require(doMC)
require(doSNOW)
require(dynatopmodel)
require(foreach)
require(ggplot2)
require(grid)
require(gridExtra)
require(lightgbm)
require(lubridate)
require(MASS)
require(pbapply)
require(plyr)
require(png)
require(purrr)
require(raster)
require(RColorBrewer)
require(readxl)
require(reshape2)
require(rgdal)
require(scales)
require(sf)
require(sp)
require(stringr)
require(stringi)
require(tidyr)
require(tidyverse)
require(vegan)
require(velox) # devtools::install_github("hunzikp/velox")
require(R.devices) # https://www.jottr.org/2018/07/21/suppressgraphics/

N_CORE_LARGE = 24 # max. number of cores to use
N_CORE_SMALL = 6 # number of cores for memmory intensive tasks

setwd("/home/chrisgraywolf/shared/analysis/Andrews_simpleGBM/") # set to project directory

set.tempdir("temp")
setPaths(cachePath="temp", inputPath="temp", modulePath="temp", outputPath="temp", silent = FALSE)

dir.create("data_processed")
dir.create("output")

registerDoMC(N_CORE_LARGE)

start = Sys.time()
source("scripts/utility_functions.R")
source("scripts/001 - raster setup.R")
source("scripts/002 - set up temperature.R")
source("scripts - Andrews2/003 - spatial prediction.R")
Sys.time() - start

