# Simple HJA microclimate models

Code for microclimate temperature modeling in the HJ Andrews forest. Before running the code, you should set up your directory structure as follows:

```
your project directory
├── data
│    └── data_input
└── scripts
```

Where the R scripts in this repository should be placed in the "scripts" folder and the following datasets should be placed in the "data\_input" folder:
1. HJA variables\_final.xlsx [available in this repository]
2. site\_locations.xlsx [available in this repository]
3. temperature\_cleaned [folder with temperature data csv's]
4. rasters listed in "HJA variables\_final.xlsx" [see HJA data webpage]
5. harvest, hjaveg8, shapefiles [for PC1/PC2 only; see HJA data webpage, I think]

Before running the code, you will need to have the R packages listed in "runner.R" installed and to set the project directory in "runner.R" (note "velox" is no longer on CRAN). You might also want to lower N\_CORE\_SMALL and N\_CORE\_LARGE in this script if you have limited cpu/memory. The doMC package is used for parallelizing the code, which works in Linux (and possibly MacOS), but may not work in Windows. If you are using Windows, you can skip this package and parallelism will be  disabled (with some warnings).

The output of the code, which can be obtained by executing "runner.R", is a set of predicted temperature rasters (with a corresponding paneled image) corresponding to the temperature metrics in Wolf et al. (2021) along with mean min/mean/max daily temperature for each month. All relevant methods (e.g., raster processing) are described in Wolf et al. (2021) with two differences:

1. monthly temperature metrics are included here
2. all models are for the mean response only (no quantile regression is used)

Reference:

Christopher Wolf, David M. Bell, Hankyu Kim, Michael Paul Nelson, Mark Schulze, and Matthew G. Betts. Temporal consistency of undercanopy thermal refugia in old-growth forest. *Agricultural and Forest Meteorology*, In Review.
