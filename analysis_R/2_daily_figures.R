rm(list = ls())

library(tidyverse)
library(Metrics)
library(hydroGOF)
library(ggpubr)
library(stringr)
library(cowplot)
library(ggpmisc)

# loading data
df <- read.csv("combined_data.csv")
df$DATE <- as.Date(df$DATE)
head(df)

### meta data
meta <- read.csv("./flux_ET_dataset/station_metadata.csv",skip = 1)
colnames(meta)[1] <- "site"
head(meta)

meta_use <- meta %>% 
	select(site,General.classification) %>%
	rename(site=site, 
				 class0 = General.classification)

df0 <- left_join(df,meta_use)

#################################
######### Figures
#################################

# Note: Figures 3 and S3 are not included here because the SFE calculation requires 
# satellite-derived available energy and gridMET-driven meteorological variables.

### Figure 4 and S4
df0 %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	ggplot(aes(ET_corr,RADET)) +
	geom_point(alpha = 0.2) +
	geom_abline(lty=2) +
	facet_wrap(~class0, nrow = 2) +
	theme_bw() +
	geom_smooth(method = "lm",formula = 'y~x+0') +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label))),
		parse = TRUE,
		coef.digits = 2
	) +
	stat_regline_equation(aes(label =  paste(..rr.label..)),
												label.y = 10,label.x = 0) +
	labs(x = expression("EBR-corrected in situ ET ("*mm~d^-1*")"),
			 y = expression("RADET ("*mm~d^-1*")"),
	) +
	coord_cartesian(xlim = c(0,11.8),ylim = c(0,11.8))

df0 %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET > 0) %>%
	ggplot(aes(ET,RADET)) +
	geom_point(alpha = 0.2) +
	geom_abline(lty=2) +
	facet_wrap(~class0) +
	theme_bw() +
	geom_smooth(method = "lm",formula = 'y~x+0') +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label))),
		parse = TRUE,
		coef.digits = 2
	) +
	stat_regline_equation(aes(label =  paste(..rr.label..)),
												label.y = 10,label.x = 0) +
	labs(x = expression("EBR-uncorrected in situ ET ("*mm~d^-1*")"),
			 y = expression("RADET ("*mm~d^-1*")"),
	) +
	coord_cartesian(xlim = c(0,11.8),ylim = c(0,11.8))

### Figure S7
# US-Me sites
df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	filter(site %in% unique(site[grepl("US-Me",site)]))%>%
	group_by(class0) %>%
	summarise(Ensemble = KGE(Ensemble,ET_corr),
						RADET = KGE(RADET,ET_corr)) 
df0 %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	filter(site %in% unique(site[grepl("US-Me",site)]))%>%
	ggplot(aes(ET_corr,Ensemble)) +
	geom_point(alpha = 0.8,aes(shape = site,color ="Ensemble")) +
	geom_point(alpha = 0.8,aes(shape = site,color ="RADET", y = RADET)) +
	geom_abline(lty=2) +
	theme_bw() +
	geom_smooth(method = "lm",formula = 'y~x+0',aes(color ="Ensemble")) +
	geom_smooth(method = "lm",formula = 'y~x+0',aes(color ="RADET", y = RADET)) +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label)),color ="Ensemble"),
		parse = TRUE,
		coef.digits = 2.5, label.y = 0.95
	) +
	geom_label(aes(x=3,y=7.7, label = "KGE = 0.32"), label.size = 0, color = "#D55E00") +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label)),color ="RADET", y = RADET),
		parse = TRUE,
		coef.digits = 2.5, label.y = 0.87
	) +
	geom_label(aes(x=3,y=7, label = "KGE = 0.64"), label.size = 0, color = "#0072B2") +
	labs(x = expression("EBR-corrected in situ ET ("*mm~d^-1*")"),
			 y = expression("Model ("*mm~d^-1*")"),
			 title = "a. Evergreen forests in Oregon"
	) + 
	scale_color_manual(
		name = "model",
		values = c("Ensemble" = "#D55E00", "RADET" = "#0072B2")) 


### Nevada shrubland sites (spring valley)
df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	filter(site %in% c("SPV_3","SV_5","SV_6"))%>%
	summarise(Ensemble = KGE(Ensemble,ET_corr),
						RADET = KGE(RADET,ET_corr)) 

df0 %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	filter(site %in% c("SPV_3","SV_5","SV_6"))%>%
	ggplot(aes(ET_corr,Ensemble)) +
	geom_point(alpha = 0.8,aes(shape = site,color ="Ensemble")) +
	geom_point(alpha = 0.8,aes(shape = site,color ="RADET", y = RADET)) +
	geom_abline(lty=2) +
	theme_bw() +
	geom_smooth(method = "lm",formula = 'y~x+0',aes(color ="Ensemble")) +
	geom_smooth(method = "lm",formula = 'y~x+0',aes(color ="RADET", y = RADET)) +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label)),color ="Ensemble"),
		parse = TRUE,
		coef.digits = 2, label.y = 0.95
	) +
	geom_label(aes(x=3.5,y=5.7, label = "KGE = 0.43"), label.size = 0, color = "#D55E00") +
	stat_poly_eq(
		formula = y ~ x + 0,
		aes(label = paste(
			after_stat(eq.label)),color ="RADET", y = RADET),
		parse = TRUE,
		coef.digits = 2, label.y = 0.87
	) +
	geom_label(aes(x=3.5,y=5.19, label = "KGE = 0.71"), label.size = 0, color = "#0072B2") +
	labs(x = expression("EBR-corrected in situ ET ("*mm~d^-1*")"),
			 y = expression("Model ("*mm~d^-1*")"),
			 title = "b. Shrublands & grasslands in Nevada"
	)  + 
	scale_color_manual(
		name = "model",
		values = c("Ensemble" = "#D55E00", "RADET" = "#0072B2")
	)


### Figure 5 (replace ET_corr to ET for Figure S5)
stat_table <- df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	group_by(class0, site) %>%
	summarise(Ensemble = max(-1,KGE(Ensemble,ET_corr)),
						RADET = max(-1,KGE(RADET,ET_corr)),
						PT.JPL = max(-1,KGE(PT.JPL,ET_corr)),
						geeSEBAL = max(-1,KGE(geeSEBAL,ET_corr)),
						eeMETRIC = max(-1,KGE(eeMETRIC,ET_corr)),
						SSEBop = max(-1,KGE(SSEBop,ET_corr)),
						DisALEXI = max(-1,KGE(DisALEXI,ET_corr)),
						length = n(),
						weight = sqrt(length),
						
	) %>%
	filter(length > 5) %>%
	group_by(class0) %>%
	summarise(Ensemble = weighted.mean(Ensemble,w = weight),
						RADET = weighted.mean(RADET,w = weight),
						PT.JPL = weighted.mean(PT.JPL,w = weight),
						geeSEBAL = weighted.mean(geeSEBAL,w = weight),
						eeMETRIC = weighted.mean(eeMETRIC,w = weight),
						SSEBop = weighted.mean(SSEBop,w = weight),
						DisALEXI = weighted.mean(DisALEXI,w = weight),
	)

stat_table_long <- stat_table %>% gather(key = "model",value = "KGE", -class0)
stat_table_long$model <- factor(stat_table_long$model,
																levels = c("geeSEBAL","PT.JPL","SSEBop","eeMETRIC",
																					 "DisALEXI","Ensemble","RADET"))
a <- stat_table_long %>%
	ggplot(aes(class0,KGE,fill = model)) +
	geom_bar(stat = "identity",
					 position = position_dodge(width = 0.8),
					 width = 0.8,
					 color = "black") +
	scale_fill_manual(values = c("gray100", "gray80","gray60","gray40","gray20",
															 "wheat","red")) +
	theme_bw() +
	labs(y=expression("KGE"),
			 x = "")  + theme(axis.title.x = element_blank())


stat_table <- df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	group_by(class0, site) %>%
	summarise(Ensemble = NSE(Ensemble,ET_corr),
						RADET = NSE(RADET,ET_corr),
						PT.JPL = NSE(PT.JPL,ET_corr),
						geeSEBAL = NSE(geeSEBAL,ET_corr),
						eeMETRIC = NSE(eeMETRIC,ET_corr),
						SSEBop = NSE(SSEBop,ET_corr),
						DisALEXI = NSE(DisALEXI,ET_corr),
						length = n(),
						weight = sqrt(length),
						
	) %>%
	mutate(Ensemble = ifelse(Ensemble < -1, -1, Ensemble),
				 RADET = ifelse(RADET < -1, -1, RADET),
				 PT.JPL = ifelse(PT.JPL < -1, -1, PT.JPL),
				 geeSEBAL = ifelse(geeSEBAL < -1, -1, geeSEBAL),
				 eeMETRIC = ifelse(eeMETRIC < -1, -1, eeMETRIC),
				 SSEBop = ifelse(SSEBop < -1, -1, SSEBop),
				 DisALEXI = ifelse(DisALEXI < -1, -1, DisALEXI)
	) %>%
	filter(length > 5) %>%
	group_by(class0) %>%
	summarise(Ensemble = weighted.mean(Ensemble,w = weight),
						RADET = weighted.mean(RADET,w = weight),
						PT.JPL = weighted.mean(PT.JPL,w = weight),
						geeSEBAL = weighted.mean(geeSEBAL,w = weight),
						eeMETRIC = weighted.mean(eeMETRIC,w = weight),
						SSEBop = weighted.mean(SSEBop,w = weight),
						DisALEXI = weighted.mean(DisALEXI,w = weight),
	)

stat_table_long <- stat_table %>% gather(key = "model",value = "NSE", -class0)

stat_table_long$model <- factor(stat_table_long$model,
																levels = c("geeSEBAL","PT.JPL","SSEBop","eeMETRIC",
																					 "DisALEXI","Ensemble","RADET"))
b <- stat_table_long %>%
	ggplot(aes(class0,NSE,fill = model)) +
	geom_bar(stat = "identity",
					 position = position_dodge(width = 0.8),
					 width = 0.8,
					 color = "black") +
	scale_fill_manual(values = c("gray100", "gray80","gray60","gray40","gray20",
															 "wheat","red")) +
	theme_bw() +
	labs(y=expression("NSE"),
			 x = "")  +
	coord_cartesian(ylim = c(-0.75,0.75)) + theme(axis.title.x = element_blank())

stat_table <- df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	group_by(class0, site) %>%
	summarise(Ensemble = rmse(Ensemble,ET_corr),
						RADET = rmse(RADET,ET_corr),
						PT.JPL = rmse(PT.JPL,ET_corr),
						geeSEBAL = rmse(geeSEBAL,ET_corr),
						eeMETRIC = rmse(eeMETRIC,ET_corr),
						SSEBop = rmse(SSEBop,ET_corr),
						DisALEXI = rmse(DisALEXI,ET_corr),
						length = n(),
						weight = sqrt(length),
						
	) %>%
	filter(length > 3) %>%
	group_by(class0) %>%
	summarise(Ensemble = weighted.mean(Ensemble,w = weight),
						RADET = weighted.mean(RADET,w = weight),
						PT.JPL = weighted.mean(PT.JPL,w = weight),
						geeSEBAL = weighted.mean(geeSEBAL,w = weight),
						eeMETRIC = weighted.mean(eeMETRIC,w = weight),
						SSEBop = weighted.mean(SSEBop,w = weight),
						DisALEXI = weighted.mean(DisALEXI,w = weight),
	)

stat_table_long <- stat_table %>% gather(key = "model",value = "NSE", -class0)

stat_table_long$model <- factor(stat_table_long$model,
																levels = c("geeSEBAL","PT.JPL","SSEBop","eeMETRIC",
																					 "DisALEXI","Ensemble","RADET"))
c <- stat_table_long %>%
	ggplot(aes(class0,NSE,fill = model)) +
	geom_bar(stat = "identity",
					 position = position_dodge(width = 0.8),
					 width = 0.8,
					 color = "black") +
	scale_fill_manual(values = c("gray100", "gray80","gray60","gray40","gray20",
															 "wheat","red")) +
	theme_bw() +
	labs(y=expression("RMSE (mm "*d^-1*")"),
			 x = "")   + theme(axis.title.x = element_blank())

stat_table <- df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	group_by(class0, site) %>%
	summarise(Ensemble = mae(Ensemble,ET_corr),
						RADET = mae(RADET,ET_corr),
						PT.JPL = mae(PT.JPL,ET_corr),
						geeSEBAL = mae(geeSEBAL,ET_corr),
						eeMETRIC = mae(eeMETRIC,ET_corr),
						SSEBop = mae(SSEBop,ET_corr),
						DisALEXI = mae(DisALEXI,ET_corr),
						length = n(),
						weight = sqrt(length),
						
	) %>%
	filter(length > 3) %>%
	group_by(class0) %>%
	summarise(Ensemble = weighted.mean(Ensemble,w = weight),
						RADET = weighted.mean(RADET,w = weight),
						PT.JPL = weighted.mean(PT.JPL,w = weight),
						geeSEBAL = weighted.mean(geeSEBAL,w = weight),
						eeMETRIC = weighted.mean(eeMETRIC,w = weight),
						SSEBop = weighted.mean(SSEBop,w = weight),
						DisALEXI = weighted.mean(DisALEXI,w = weight),
	)

stat_table_long <- stat_table %>% gather(key = "model",value = "NSE", -class0)

stat_table_long$model <- factor(stat_table_long$model,
																levels = c("geeSEBAL","PT.JPL","SSEBop","eeMETRIC",
																					 "DisALEXI","Ensemble","RADET"))
d <- stat_table_long %>%
	ggplot(aes(class0,NSE,fill = model)) +
	geom_bar(stat = "identity",
					 position = position_dodge(width = 0.8),
					 width = 0.8,
					 color = "black") +
	scale_fill_manual(values = c("gray100", "gray80","gray60","gray40","gray20",
															 "wheat","red")) +
	theme_bw() +
	labs(y=expression("MAE (mm "*d^-1*")"),
			 x = "")   + theme(axis.title.x = element_blank())


mbe <- function (predicted,actual) {return(mean(predicted - actual))}
stat_table <- df0 %>%
	ungroup() %>%
	filter(!is.na(Ensemble) & !is.na(RADET) & ET_corr > 0) %>%
	group_by(class0, site) %>%
	summarise(Ensemble = mbe(Ensemble,ET_corr),
						RADET = mbe(RADET,ET_corr),
						PT.JPL = mbe(PT.JPL,ET_corr),
						geeSEBAL = mbe(geeSEBAL,ET_corr),
						eeMETRIC = mbe(eeMETRIC,ET_corr),
						SSEBop = mbe(SSEBop,ET_corr),
						DisALEXI = mbe(DisALEXI,ET_corr),
						length = n(),
						weight = sqrt(length),
						
	) %>%
	filter(length > 3) %>%
	group_by(class0) %>%
	summarise(Ensemble = weighted.mean(Ensemble,w = weight),
						RADET = weighted.mean(RADET,w = weight),
						PT.JPL = weighted.mean(PT.JPL,w = weight),
						geeSEBAL = weighted.mean(geeSEBAL,w = weight),
						eeMETRIC = weighted.mean(eeMETRIC,w = weight),
						SSEBop = weighted.mean(SSEBop,w = weight),
						DisALEXI = weighted.mean(DisALEXI,w = weight),
	)

stat_table_long <- stat_table %>% gather(key = "model",value = "NSE", -class0)

stat_table_long$model <- factor(stat_table_long$model,
																levels = c("geeSEBAL","PT.JPL","SSEBop","eeMETRIC",
																					 "DisALEXI","Ensemble","RADET"))
e <- stat_table_long %>%
	ggplot(aes(class0,NSE,fill = model)) +
	geom_bar(stat = "identity",
					 position = position_dodge(width = 0.8),
					 width = 0.8,
					 color = "black") +
	scale_fill_manual(values = c("gray100", "gray80","gray60","gray40","gray20",
															 "wheat","red")) +
	theme_bw() +
	labs(y=expression("MBE (mm "*d^-1*")"),
			 x = "")   + theme(axis.title.x = element_blank())

plot_grid(plot_grid(a+theme(legend.position = "none"),
										b+theme(legend.position = "none"),
										c+theme(legend.position = "none"),
										d+theme(legend.position = "none"),
										e+theme(legend.position = "none"),
										ncol=1,
										labels = "auto"),
					get_legend(a + labs(fill = "ET model")),
					ncol = 2,rel_widths = c(1,0.2))


#### KGE improvement
KEG_improvement <- df0 %>%
	group_by(site) %>%
	summarise(Ensemble = KGE(Ensemble,ET_corr),
						RADET = KGE(RADET,ET_corr),
						dif = RADET-Ensemble
	)

write.csv(KEG_improvement,"KEG_improvement_daily.csv",row.names = F)
