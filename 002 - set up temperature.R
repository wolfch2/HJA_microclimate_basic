# https://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r
monthdays = Vectorize(function(m, y){
  y = as.numeric(y); m=as.numeric(m);
  # Quick check for leap year
  leap = 0
  if (y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0))
    leap = 1

  # Return the number of days in the month
  return(switch(m,
                31,
                28 + leap,  # adds 1 if leap year
                31,
                30,
                31,
                30,
                31,
                31,
                30,
                31,
                30,
                31))
})

# compute mean/max/min provided at least prop_obs proportion of obs. are present
missing_mean = function(x, n_days, prop_obs=4/5){
        out = mean(x)
        if(length(x) < prop_obs * n_days) out = NA
        return(out)
}
missing_max = function(x, n_days, prop_obs=4/5){
        out = max(x)
        if(length(x) < prop_obs * n_days) out = NA
        return(out)
}
missing_min = function(x, n_days, prop_obs=4/5){
        out = min(x)
        if(length(x) < prop_obs * n_days) out = NA
        return(out)
}

# sitedata contains daily temperature summaries for a given site and year
# (SiteInfo, Date, MaxT, MinT, MeanT)
T_custom_metrics = function(sitedata, prop_obs = 0.8){ # see StreamThermal package for input format
        sitedata = separate(sitedata, "Date", c("Year", "Month", "Day"), sep = "-")
        sitedata$Month = as.numeric(sitedata$Month)
        y = sitedata$Year[1]
        out = data.frame(SiteInfo=sitedata$SiteInfo[1], # max daily mean looks better than mean daily max (maybe mean daily max is getting thrown off more by errors?)
                         Jul_Sep_mean_max=missing_mean(sitedata$MaxT[sitedata$Month %in% 7:9],sum(monthdays(7:9,y)),prop_obs),
                         Jul_Sep_mean_mean=missing_mean(sitedata$MeanT[sitedata$Month %in% 7:9],sum(monthdays(7:9,y)),prop_obs),
                         GDD_winter_5=ifelse(sum(sitedata$Month <= 3) == sum(monthdays(1:3,y)), sum(pmax((sitedata$MaxT + sitedata$MinT)/2 - 5, 0)[sitedata$Month <= 3]), NA),
                         Apr_Jun_mean_max=missing_mean(sitedata$MaxT[sitedata$Month %in% 4:6],sum(monthdays(4:6,y)),prop_obs),
                         Apr_Jun_mean_mean=missing_mean(sitedata$MeanT[sitedata$Month %in% 4:6],sum(monthdays(4:6,y)),prop_obs),                         
                         Apr_Jun_mean_min=missing_mean(sitedata$MinT[sitedata$Month %in% 4:6],sum(monthdays(4:6,y)),prop_obs),
                         stringsAsFactors=FALSE)
        for(m in 1:12){
                out[,paste0("m_",m,"_min")] = missing_mean(sitedata$MinT[sitedata$Month %in% m],sum(monthdays(m,y)),prop_obs)
                out[,paste0("m_",m,"_mean")] = missing_mean(sitedata$MeanT[sitedata$Month %in% m],sum(monthdays(m,y)),prop_obs)
                out[,paste0("m_",m,"_max")] = missing_mean(sitedata$MaxT[sitedata$Month %in% m],sum(monthdays(m,y)),prop_obs)
        }
        return(out)
}

############################## set up temperature data

temp_files = list.files("data_input/temperature_cleaned",full.names=TRUE)
temperature = rbindlist(pblapply(temp_files, function(temp_file){
        name = as.numeric(gsub(".*BIRD_|015_cleaned.*","",temp_file))
        temp = fread(temp_file)        
        temp[, "LOCATION_CODE" := name] # https://stackoverflow.com/questions/19072053/adding-columns-to-a-data-table
        return(temp)
}))
date_time_sep = stri_split_fixed(temperature$DATE_TIME, " ", 2, simplify=TRUE)
temperature[, "Date" := date_time_sep[,1]]
temperature[, "Time" := date_time_sep[,2]]

# https://stackoverflow.com/questions/12064202/apply-several-summary-functions-on-several-variables-by-group-in-one-call
temperature = temperature[ , .(AIRTEMP_MIN_DAY = min(TEMP_C),
                                AIRTEMP_MEAN_DAY = mean(TEMP_C),
                                AIRTEMP_MAX_DAY = max(TEMP_C)), by = .(LOCATION_CODE, Date)] # aggregate to day level

temperature = separate(temperature, "Date", c("Year","Month","Day"), "-", convert=TRUE, remove=FALSE)
temperature$Date = ymd(temperature$Date)

############################## compute annual metrics

temp_split = split(temperature,paste(temperature$Year,temperature$LOCATION_CODE))
temperature_metrics = foreach(i=1:length(temp_split)) %dopar% {
        print(i)
        df = temp_split[[i]]
        sitedata_micro = data.frame(SiteInfo=df$LOCATION_CODE,
                      Date=df$Date,
                      MaxT=df$AIRTEMP_MAX_DAY,
                      MinT=df$AIRTEMP_MIN_DAY,
                      MeanT=df$AIRTEMP_MEAN_DAY,
                      stringsAsFactors=FALSE)
        temps = t(T_custom_metrics(sitedata_micro)[,-1])
        out = data.frame(site=df$LOCATION_CODE[1],
                   Year=df$Year[1],
                   variable=rownames(temps),
                   value=temps,
                   stringsAsFactors=FALSE)
        return(out)
}
temperature_metrics = data.frame(rbindlist(temperature_metrics))
temperature_metrics = temperature_metrics[! is.na(temperature_metrics$value),]
table(temperature_metrics$variable, temperature_metrics$Year) # looks fine..

temperature_metrics$POINT = temperature_metrics$LOCATION_CODE = temperature_metrics$SITECODE = temperature_metrics$site

saveRDS(temperature_metrics, "data_processed/temperature_metrics.RDS")
