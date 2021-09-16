############################## set up data

rast_reduce = readRDS("data_processed/rast_reduce.RDS")
rast_reduce = rast_reduce[! rast_reduce$variable %in% c("PC1","PC2"),] # not using PC vars for main models
rast_reduce$var_scale = paste0(rast_reduce$variable, "XX", rast_reduce$scale)
rast_spread = tidyr::spread(rast_reduce[,c("site","var_scale","value"),],key="var_scale",value="value")

temp_merged = readRDS("data_processed/temperature_metrics.RDS")
temp_merged$site_year = paste(temp_merged$LOCATION_CODE, temp_merged$Year)
if(years != "All"){
        temp_merged = temp_merged[temp_merged$Year %in% years,]
}

data = merge(temp_merged, rast_spread, by="site")

############################## build predictor matrix for newdata

predictors = readRDS("data_processed/predictors.RDS")
predictor_mat = as.matrix(predictors)
pred_mask = readRDS("data_processed/pred_mask.RDS")

############################## fit models!

pred_rast_list = foreach(var=sort(unique(data$variable))) %do% {  # startup is slow since predictor_mat big
        print(var)
        gc() # seems to help w/ memory allocation errors
        covar_data = data[data$variable==var,]

        if(framework == "lightgbm"){
                # https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html
                # can likely predict faster w/ treelite https://github.com/Microsoft/LightGBM/issues/2094
                mod = lightgbm(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                        label = covar_data$value,
                        num_leaves = round(0.75*2^5),
                        learning_rate = 1,
                        objective = "regression",
                        nrounds = 25,
                        nthread = 1) # nthread > 1 also enables parallel prediction!
                pred = predict(mod, data = predictor_mat[,unique(rast_reduce$var_scale)])
        } else{ # framework == "xgboost"
                mod = xgboost(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                        label = covar_data$value,
                        booster = "gbtree",
                        objective = "reg:linear",
                        nrounds=25,
                        max_depth=5,
                        nthread=1)
                pred = predict(mod, newdata = predictor_mat[,unique(rast_reduce$var_scale)])
        }

        out = pred_mask
        out[] = pred
        out[pred_mask[] == 0] = NA
        return(out)
}

pred_rast_stack = stack(pred_rast_list)
names(pred_rast_stack) = sort(unique(data$variable))
writeRaster(pred_rast_stack, filename=paste0("output/",names(pred_rast_stack)), bylayer=TRUE,format="GTiff", overwrite=TRUE)

############################## plot prediction rasters

temps = readRDS("data_processed/temperature_metrics.RDS")
sites = readRDS("data_processed/sites.RDS")
sites= merge(temps, sites, by.x="SITECODE", by.y="LOCATION_CODE")

map_theme = theme(axis.ticks=element_blank(),
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_line(colour="transparent"),
                      panel.grid.minor = element_line(colour="transparent"),
                      plot.margin=margin(0,1,0,0),
                      legend.position=c(0.02,0.98),
                      legend.justification=0:1,
                      legend.title=element_text(size=6.5),
                      legend.text=element_text(size=6.5),
                      legend.key.height=unit(0.4,"lines"),
                      legend.key.width=unit(0.6,"lines"),
                      legend.background=element_rect(fill=NA)) 

plot_list = lapply(names(pred_rast_stack), function(var){
        print(var)
        rast = pred_rast_stack[[var]]
        rast = aggregate(rast, 4) # optional - aggregate to reduce file size

        rast_pts <- data.frame(rasterToPoints(rast))
        colnames(rast_pts) <- c('x','y','value')

        sites_var = sites[sites$variable == var,]
        sites_var = sites_var[! duplicated(sites_var$LOCATION_CODE),]
        print(nrow(sites_var))

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           map_theme

        return(p)
})

pred_plot = plot_grid(plotlist=plot_list, nrow=7)

png("output/pred_map_all.png", width=13, height=11, units="in", res=200)
plot(pred_plot)
dev.off()
