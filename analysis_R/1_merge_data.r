rm(list = ls())

library(tidyverse)
library(stringr)
library(lubridate)

#####################
## 1. OpenET
#####################
OpenET <- read.csv("./OpenET_PhaseII_model_ET_dataset/daily_data.csv")

# site name
colnames(OpenET)[1] <- "site"
site <- unique(OpenET$site)

# Date as r date
OpenET$DATE <- as.Date(OpenET$DATE,format = "%m/%d/%Y")

#####################
## 2. Flux data
#####################
getwd()
setwd("./flux_ET_dataset/daily_data_files")
getwd()

for(i in c(1:length(site))){
	file_name <- paste0(site[i],"_daily_data.csv")
	df <- read.csv(file_name)
	
	# Date as r Date
	colnames(df)[1] <- "DATE"
	df$DATE <- as.Date(df$DATE)
	
	# select required variables
	
	if(sum(colnames(df) == 'Rn') == 0){df$Rn <- NA}
	if(sum(colnames(df) == 'G') == 0){df$G <- NA}
	if(sum(colnames(df) == 'ebr') == 0){df$ebr <- NA}
	
	if(sum(colnames(df) == 'LE') == 0){df$LE <- NA}
	if(sum(colnames(df) == 'H') == 0){df$H <- NA}
	if(sum(colnames(df) == 'LE_corr') == 0){df$LE_corr <- NA}
	if(sum(colnames(df) == 'H_corr') == 0){df$H_corr <- NA}
	
	if(sum(colnames(df) == 'LE_subday_gaps') == 0){df$LE_subday_gaps <- NA}
	if(sum(colnames(df) == 'ET_gap') == 0){df$ET_gap <- NA}
	
	if(sum(colnames(df) == 'ET') == 0){df$ET <- NA}
	if(sum(colnames(df) == 'ET_corr') == 0){df$ET_corr <- NA}

	if(sum(colnames(df) == 'ws') == 0){df$ws <- NA}
	if(sum(colnames(df) == 't_avg') == 0){df$t_avg <- NA}
	if(sum(colnames(df) == 't_min') == 0){df$t_min <- NA}
	if(sum(colnames(df) == 't_max') == 0){df$t_max <- NA}
	if(sum(colnames(df) == 'rh') == 0){df$rh <- NA}
	if(sum(colnames(df) == 't_dew') == 0){df$t_dew <- NA}
	if(sum(colnames(df) == 'vpd') == 0){df$vpd <- NA}
	if(sum(colnames(df) == 'vp') == 0){df$vp <- NA}
	if(sum(colnames(df) == 'rso') == 0){df$rso <- NA}
	if(sum(colnames(df) == 'sw_in') == 0){df$sw_in <- NA}
	if(sum(colnames(df) == 'sw_out') == 0){df$sw_out <- NA}
	if(sum(colnames(df) == 'lw_in') == 0){df$lw_in <- NA}
	if(sum(colnames(df) == 'lw_out') == 0){df$lw_out <- NA}
	if(sum(colnames(df) == 'ppt') == 0){df$ppt <- NA}
	
	if(sum(colnames(df) == 'ASCE_ETo') == 0){df$ASCE_ETo <- NA}
	
	df <- df %>% select(DATE, Rn, G, H, LE, ebr, H_corr, LE_corr, LE_subday_gaps,
											ET_gap, ET, ET_corr, ws, t_avg, t_min, t_max, rh,
											t_dew, vpd, vp, rso, sw_in, sw_out, lw_in, lw_out,ppt, ASCE_ETo,
											gridMET_ETr, gridMET_ETo, gridMET_prcp)
	
	df <- df %>% filter(!is.na(ET_corr))
	
	# add site column
	df$site <- site[i]
	
	#merge OpenET data
	merged <- left_join(df,OpenET,by = c("site","DATE"))

	#merge all sites
	if(i==1){
		Flux_OpenET <- merged
	} else {
		Flux_OpenET <- full_join(Flux_OpenET,merged)
	}
}

setwd("..")
setwd("..")

#####################
## 3. RADET data
#####################
RADET <- read.csv("./RADET_dataset/daily_RADET.csv")
head(RADET)
RADET$DATE <- as.Date(RADET$DATE)

combined <- left_join(Flux_OpenET,RADET)

#####################
## 4. Save
#####################
write.csv(combined,"combined_data.csv",row.names = F)

