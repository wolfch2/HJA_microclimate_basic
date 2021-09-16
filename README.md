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

Before running the code, you will need to set the project directory in "runner.R". To obtain predictions for specific year(s), you can change the "years" variable in this script. Note that memory use peaks at ~5.5 GB and runtime is approximately 3.5 minutes (on my desktop computer). Optionally, you can change the "framework" variable from "lightgbm" to "xgboost" to use the xgboost gradient boosting framework instead of lightgbm.

The output of the code, which can be obtained by executing "runner.R", is a set of predicted temperature rasters (with a paneled figure) corresponding to the temperature metrics in Wolf et al. (2021), Frey et al. (2016), and mean min/mean/max daily temperature for each month. All relevant methods (e.g., raster processing) are described in Wolf et al. (2021) with three differences:

1. additional temperature metrics are included here
2. all models are for the mean response only (no quantile regression is used)
3. as in Frey et al. (2016), only 25 and 250 meter spatial scales are considered (predictions are made at 25 m resolution)

If you want to model other temperature metrics, you can add them to the "T\_custom\_metrics" function in the "002 - set up temperature.R" script. Optionally, formatted names can be added to the "format\_names" function in "utility\_functions.R".

References:

Christopher Wolf, David M. Bell, Hankyu Kim, Michael Paul Nelson, Mark Schulze, and Matthew G. Betts. Temporal consistency of undercanopy thermal refugia in old-growth forest. *Agricultural and Forest Meteorology* 307 (2021): 108520.

Frey, Sarah JK, et al. "Spatial models reveal the microclimatic buffering capacity of old-growth forests." *Science advances* 2.4 (2016): e1501392.
