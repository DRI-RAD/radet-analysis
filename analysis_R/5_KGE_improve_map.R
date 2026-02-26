rm(list = ls())
library(tidyverse)
library(readr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)  
library(tigris)

## data
df <- read.csv("combined_data.csv")
pts <- df %>% group_by(site) %>% tally()

###
meta <- read.csv("./flux_ET_dataset/station_metadata.csv",skip = 1)
colnames(meta)[1] <- "site"

meta_use <- meta %>% 
	select(site,General.classification,Latitude,Longitude) %>%
	rename(site=site, 
				 class0 = General.classification,
				 latitude = Latitude,
				 longitude = Longitude
	) %>%
	mutate(class1 = ifelse(class0 %in% c("Croplands", "Wetland/Riparian"),"Croplands/Wetland",class0))

pts <- inner_join(pts,meta_use)

########### Figure 2
# to sf (EPSG:4326)
pts_sf <- st_as_sf(pts,
									 coords = c("longitude", "latitude"),
									 crs    = 4326,
									 remove = FALSE)

# North America map
na_map <- ne_countries(continent = "North America",
											 scale = "medium",
											 returnclass = "sf")

# State boundaries
options(tigris_use_cache = TRUE)
states <- states(cb = TRUE, resolution = "20m") |> st_as_sf()

# Shape and fill values
shape_vals <- c(
	Croplands = 22,  
	Shrublands = 23,  
	Grasslands = 21,  
	`Evergreen Forests` = 24,  
	`Mixed Forests` = 25,  
	`Wetland/Riparian` = 21
)

fill_vals <- c(
	`Evergreen Forests` = "#1b9e77",
	`Mixed Forests` = "#d95f02",
	Shrublands = "#e7298a",
	Grasslands = "#66a61e",
	`Wetland/Riparian` = "#1f78b4",
	Croplands = "#e6ab02"
)

size_vals <- c(Croplands = 3,  
							 Shrublands = 3,  
							 Grasslands = 3,  
							 `Evergreen Forests` = 3,  
							 `Mixed Forests` = 3,  
							 `Wetland/Riparian` = 2
) 

# Slight jitter for point separation
pts_sf_jit <- st_jitter(pts_sf, factor = 0.01)

# visualization (Figure 2)
ggplot() +
	# geom_sf(data = na_map,   fill = "floralwhite", color = "gray60", size = 0.3) +
	geom_sf(data = states,   fill = "floralwhite",  color = "gray40", size = 0.3) +
	geom_sf(data = pts_sf_jit,
					aes(shape = class0, fill = class0, size = class0),
					color  = "black",stroke = 0.3, alpha = 0.8) +
	scale_shape_manual(values = shape_vals, name = "Land cover") +
	scale_size_manual(values = size_vals, name = "Land cover") +
	scale_fill_manual(values = fill_vals,  name = "Land cover") +
	# ✅ Curved projection, but limits stay in degrees (EPSG:4326)
	coord_sf(
		crs = "+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-95",   # or "+proj=ortho +lat_0=40 +lon_0=-95"
		default_crs = st_crs(4326),                         # data & limits in lon/lat
		xlim = c(-125, -69),
		ylim = c(25, 52),
		expand = FALSE
	) +
	
	# labs(x = "Longitude", y = "Latitude") +
	theme_minimal() +
	# annotation_scale(location = "bl",
	# 								 width_hint = 0.2,
	# 								 pad_x = unit(0.25, "cm"),
	# 								 pad_y = unit(0.25, "cm")) +
	theme(legend.position = "right") 


########### Figure 11
### KGE improvement
KGE_improvement <- read.csv("KEG_improvement.csv")
pts2 <- inner_join(KGE_improvement,pts)
pts2 <- pts2 %>% filter(!is.na(dif))

# to sf (EPSG:4326)
pts_sf <- st_as_sf(pts2,
									 coords = c("longitude", "latitude"),
									 crs    = 4326,
									 remove = FALSE)

shape_vals <- c(Croplands = 22,  
								Shrublands = 23,  
								Grasslands = 21,  
								`Evergreen Forests` = 24,  
								`Mixed Forests` = 25,
								`Wetland/Riparian` = 21
)  

size_vals <- c(Croplands = 3,  
							 Shrublands = 3,  
							 Grasslands = 3,  
							 `Evergreen Forests` = 3,  
							 `Mixed Forests` = 3,
							 `Wetland/Riparian` = 2
)  

pts_sf_jit <- st_jitter(pts_sf, factor = 0.01) 


# Figure 11
ggplot() +
	# geom_sf(data = na_map,   fill = "floralwhite", color = "gray60", size = 0.3) +
	geom_sf(data = states,   fill = "floralwhite",            color = "gray40", size = 0.3) +
	geom_sf(data = pts_sf_jit %>% arrange(abs(dif)),
					aes(shape = class0, fill = dif, size = class0),
					color  = "black", stroke = 0.3,alpha = 0.8) +
	scale_shape_manual(values = shape_vals, name = "Land cover") +
	scale_size_manual(values = size_vals, name = "Land cover") +
	scale_fill_gradient2(
		limits = c(-0.401, 0.401),           # set scale range
		oob = scales::squish,            # force out-of-range values to boundary color
		low = "red", mid = "floralwhite", high = "blue", midpoint = 0,
		name = "Δ KGE"                     # optional label
	) +
	# ✅ Curved projection, but limits stay in degrees (EPSG:4326)
	coord_sf(
		crs = "+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-95",   # or "+proj=ortho +lat_0=40 +lon_0=-95"
		default_crs = st_crs(4326),                         # data & limits in lon/lat
		xlim = c(-125, -69),
		ylim = c(25, 52),
		expand = FALSE
	) +
	
	# labs(x = "Longitude", y = "Latitude") +
	# annotation_scale(location = "bl",
	# 								 width_hint = 0.2,
	# 								 pad_x = unit(0.25, "cm"),
	# 								 pad_y = unit(0.25, "cm")) +
	theme_minimal() +
	theme(legend.position = "right") 



############################### daily (Figure 6)
KGE_improvement <- read.csv("KEG_improvement_daily.csv")
pts2 <- inner_join(KGE_improvement,pts)
pts2 <- pts2 %>% filter(!is.na(dif))

# to sf (EPSG:4326)
pts_sf <- st_as_sf(pts2,
									 coords = c("longitude", "latitude"),
									 crs    = 4326,
									 remove = FALSE)

shape_vals <- c(Croplands = 22,  
								Shrublands = 23,  
								Grasslands = 21,  
								`Evergreen Forests` = 24,  
								`Mixed Forests` = 25,
								`Wetland/Riparian` = 21
)  

size_vals <- c(Croplands = 3,  
							 Shrublands = 3,  
							 Grasslands = 3,  
							 `Evergreen Forests` = 3,  
							 `Mixed Forests` = 3,
							 `Wetland/Riparian` = 2
)  

pts_sf_jit <- st_jitter(pts_sf, factor = 0.01)   


# Figure 6
ggplot() +
	# geom_sf(data = na_map,   fill = "floralwhite", color = "gray60", size = 0.3) +
	geom_sf(data = states,   fill = "floralwhite",            color = "gray40", size = 0.3) +
	geom_sf(data = pts_sf_jit %>% arrange(abs(dif)),
					aes(shape = class0, fill = dif, size = class0),
					color  = "black", stroke = 0.3,alpha = 0.8) +
	scale_shape_manual(values = shape_vals, name = "Land cover") +
	scale_size_manual(values = size_vals, name = "Land cover") +
	scale_fill_gradient2(
		limits = c(-0.401, 0.401),           # set scale range
		oob = scales::squish,            # force out-of-range values to boundary color
		low = "red", mid = "floralwhite", high = "blue", midpoint = 0,
		name = "Δ KGE"                     # optional label
	) +
	# ✅ Curved projection, but limits stay in degrees (EPSG:4326)
	coord_sf(
		crs = "+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-95",   # or "+proj=ortho +lat_0=40 +lon_0=-95"
		default_crs = st_crs(4326),                         # data & limits in lon/lat
		xlim = c(-125, -69),
		ylim = c(25, 52),
		expand = FALSE
	) +
	
	# labs(x = "Longitude", y = "Latitude") +
	# annotation_scale(location = "bl",
	# 								 width_hint = 0.2,
	# 								 pad_x = unit(0.25, "cm"),
	# 								 pad_y = unit(0.25, "cm")) +
	theme_minimal() +
	theme(legend.position = "right") 
