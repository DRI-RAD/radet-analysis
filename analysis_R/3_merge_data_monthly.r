rm(list = ls())

library(tidyverse)
library(stringr)
library(lubridate)

#####################
## 1. OpenET
#####################
OpenET <- read.csv("./OpenET_PhaseII_model_ET_dataset/monthly_data.csv")
head(OpenET)

# site name
colnames(OpenET)[1] <- "site"
site <- unique(OpenET$site)

# Date as r date
OpenET$DATE <- as.Date(OpenET$DATE,format = "%m/%d/%Y")

OpenET <- OpenET %>%
	mutate(YEAR = year(DATE),
				 MON = month(DATE)) %>%
	select(!DATE)

#####################
## 2. Flux data
#####################

getwd()
setwd("./flux_ET_dataset/monthly_data_files")
getwd()

ls <- list.files()

for(i in c(1:length(ls))){
	file_name <- paste0(ls[i])
	df <- read.csv(file_name)
	
	# Date as r Date
	colnames(df)[1] <- "DATE"
	df$DATE <- as.Date(df$DATE)
	
	# select required variables
	
	if(sum(colnames(df) == 'ET') == 0){df$ET <- NA}
	if(sum(colnames(df) == 'ET_corr') == 0){df$ET_corr <- NA}
	
	
	df <- df %>% select(DATE, ET, ET_corr)
	df <- df %>% filter(!is.na(ET_corr))
	
	# add site column
	df$site <- str_remove(ls[i], "_monthly.*$")
	
	#merge all sites
	if(i==1){
		Flux_mon <- df
	} else {
		Flux_mon <- full_join(Flux_mon,df)
	}
}

setwd("..")
setwd("..")

Flux_mon <- Flux_mon %>%
	mutate(YEAR = year(DATE),
				 MON = month(DATE)) %>%
	select(!DATE)

######################################
## 2. Flux data (daily to integrate)
######################################
daily <- read.csv("combined_data.csv")
daily$DATE <- as.Date(daily$DATE)
head(daily)

daily <- daily %>%
	mutate(YEAR = year(DATE),
				 MON = month(DATE))

FLUX_mon_from_daily <-
	daily %>% 
	group_by(site,YEAR, MON) %>%
	summarise(ET_corr_relaxed = sum(ET_corr),
						ET_no_gap_count = sum(!is.na(ET)),
						)
head(FLUX_mon_from_daily)

Flux_OpenET <- full_join(Flux_mon,FLUX_mon_from_daily)
Flux_OpenET <- left_join(Flux_OpenET,OpenET)
head(Flux_OpenET)

#####################
## 3. RADET data
#####################
RADET <- read.csv("./RADET_dataset/monthly_RADET.csv")
head(RADET)

combined <- left_join(Flux_OpenET,RADET)
head(combined)

#####################
## 4. Save
#####################
write.csv(combined,"combined_data_monthly.csv",row.names = F)
