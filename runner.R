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
source("scripts/001 - raster setup.R") # 20 min
source("scripts/002 - set up temperature.R") # 27 sec

#source("scripts - Andrews2/000 - conceptual figure.R")

#source("scripts - Andrews2/002 - gridMET comparison.R", encoding = "Latin1") # 7 sec
#source("scripts - Andrews2/003 - set up temperature.R") # 27 sec
source("scripts - Andrews2/004 - main models.R") # 45 sec
#source("scripts - Andrews2/005 - sd models.R") # 1 sec
#source("scripts - Andrews2/006 - spatial prediction.R") # 2 min
#source("scripts - Andrews2/007 - ALE plot.R") # 4 sec
Sys.time() - start

