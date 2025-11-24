library(tidyverse)
library(mgcv)
library(gratia)
library(zoo)
library(cowplot)
library(ggpubr)

load(file="./input/Zerofill_Sets_20251107.RData") # zero filled sets for different precision thresholds (no threshold, 0.90, 0.95, and 0.99)

# GOAL: plot vocal activity & occupied sites at each different threshold; then run each through the gam to see whether phenometrics are impacted


# FOCAL SPECIES -----------------------------------------------------------

# list of species
focal_list <- c("pawr", "towa")


focal_spp_none <- trials.0fill.0 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, captured_local_date, common_name, aou4) %>%
  summarise(nVoc=sum(nhits))
focal_spp_90 <- trials.0fill.09 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, captured_local_date, common_name, aou4) %>%
  summarise(nVoc=sum(nhits))
focal_spp_95 <- trials.0fill.095 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, captured_local_date, common_name, aou4) %>%
  summarise(nVoc=sum(nhits))
focal_spp_99 <- trials.0fill.099 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, captured_local_date, common_name, aou4) %>%
  summarise(nVoc=sum(nhits))

focal_spp_none$common_name <- factor(focal_spp_none$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))
focal_spp_90$common_name <- factor(focal_spp_90$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))
focal_spp_95$common_name <- factor(focal_spp_95$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))
focal_spp_99$common_name <- factor(focal_spp_99$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))

# color values for discrete viridis
v_cols_4 <- c("#fde725", "#21918c", "#3b528b", "#440154")

thresholds_plot <- ggplot() +
  geom_line(data=focal_spp_none, aes(x=as.Date(captured_local_date), y=nVoc), color=v_cols_4[1], alpha=1, linewidth=0.7) +
  geom_line(data=focal_spp_90, aes(x=as.Date(captured_local_date), y=nVoc), color=v_cols_4[2], alpha=0.8, linewidth=0.7) +
  geom_line(data=focal_spp_95, aes(x=as.Date(captured_local_date), y=nVoc), color=v_cols_4[3],  alpha=0.8, linewidth=0.7) +
  geom_line(data=focal_spp_99, aes(x=as.Date(captured_local_date), y=nVoc), color=v_cols_4[4], alpha=0.8, linewidth=0.7) +
  theme_bw() +
  labs(x="Date", y="Daily counts of BirdNET detections", color="Precision threshold") +
  theme(legend.position = "bottom", text = element_text(size=12),
        panel.grid = element_blank(), strip.background = element_rect(fill="white") ,strip.text = element_text(face = "bold")) +
  facet_wrap(~common_name)

towa_inset <- ggplot(focal_spp_99[focal_spp_99$aou4=="towa",]) +
  geom_line(aes(x=as.Date(captured_local_date), y=nVoc), color=v_cols_4[4], alpha=0.8, linewidth=0.7) +
  scale_x_date(limits = c(as.Date("2021-04-01"), as.Date("2021-09-08"))) +
  theme_bw() +
  labs(x="", y="", color="", title="pr(TP)=0.99") +
  theme(plot.title = element_text(size=10), 
        legend.position = "none", text = element_text(size=10),
        panel.grid = element_blank(),
        panel.background = element_blank())

combined_plot <- ggdraw() + draw_plot(thresholds_plot) + draw_plot(towa_inset, x = 0.56, y = 0.55, width = 0.38, height = 0.35)

# Effect of thresholding on naive occupancy

naiveOcc_none <- trials.0fill.0 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, location, common_name, aou4) %>%
  summarise(nVoc=sum(nhits)) %>%
  group_by(threshold, common_name, aou4) %>%
  summarise(nstat = n_distinct(location[nVoc>0]))
naiveOcc_90 <- trials.0fill.09 %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, location, common_name, aou4) %>%
  summarise(nVoc=sum(nhits)) %>%
  group_by(threshold, common_name, aou4) %>%
  summarise(nstat = n_distinct(location[nVoc>0]))
naiveOcc_95 <- trials.0fill.095  %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, location, common_name, aou4) %>%
  summarise(nVoc=sum(nhits)) %>%
  group_by(threshold, common_name, aou4) %>%
  summarise(nstat = n_distinct(location[nVoc>0]))
naiveOcc_99 <- trials.0fill.099  %>% filter(aou4 %in% focal_list) %>% 
  group_by(threshold, location, common_name, aou4) %>%
  summarise(nVoc=sum(nhits)) %>%
  group_by(threshold, common_name, aou4) %>%
  summarise(nstat = n_distinct(location[nVoc>0]))

naiveOcc_none$common_name <- factor(naiveOcc_none$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))
naiveOcc_90$common_name <- factor(naiveOcc_90$common_name, levels =c("Pacific Wren", "Townsend's Warbler"))
naiveOcc_95$common_name <- factor(naiveOcc_95$common_name, levels =c("Pacific Wren", "Townsend's Warbler"))
naiveOcc_99$common_name <- factor(naiveOcc_99$common_name, levels = c("Pacific Wren", "Townsend's Warbler"))

nstat_all <- rbind(naiveOcc_none, naiveOcc_90, naiveOcc_95, naiveOcc_99)
head(nstat_all)
#nstat_all$threshold <- as.character(nstat_all$threshold)

naiveOcc_plot <- ggplot(nstat_all) +
  geom_line(aes(x=threshold, y=nstat, group=aou4)) +
  geom_point(aes(x=threshold, y=nstat), size=3) +
  geom_point(aes(x=threshold, y=nstat, color=factor(threshold)), size=1) +
  theme_bw() +
  labs(x = "Precision threshold", y="# ARUs with detections \n above threshold", color="Precision threshold") +
  facet_wrap(~common_name) +
  scale_color_manual(values=v_cols_4) +
#  scale_color_viridis_c(direction = -1, begin = 0.1, end = 0.8) +
  theme(legend.position = "bottom", text = element_text(size=12),
        panel.grid = element_blank(), strip.background = element_rect(fill="white"), strip.text = element_text(face="bold"))
  
thresh_exp_plot <- ggarrange(combined_plot, naiveOcc_plot,
                             ncol=1, heights=c(3,2), common.legend = T)

#ggsave(thresh_exp_plot, filename="output/figures/threshold_expt_plot_pawr_towa.tif",
 #      width = 6.5, height = 7, dpi=600)

# Run GAMs ----------------------------------------------------------------

# Run each species for each of 4 pr(TP) thresholds (no threshold, 0.90, 0.95, 0.99)


# all in a loop -----------------------------------------------------------

threshold_expt_dat <- rbind(trials.0fill.0, trials.0fill.09, trials.0fill.095, trials.0fill.099) %>% filter(aou4 %in% focal_list, ElevBin=="Mid")

threshold_expt_dat$ElevBin <- "Mid"

unique(threshold_expt_dat$threshold)
threshold_expt_dat$threshlab <- ifelse(threshold_expt_dat$threshold==0, "no", 
                                       ifelse(threshold_expt_dat$threshold==0.90, "pr(TP)>0.9",
                                              ifelse(threshold_expt_dat$threshold==0.95, "pr(TP)>0.95", "pr(TP)>0.99")))

# mgcv requires characters to be converted to factors
threshold_expt_dat$location <- as.factor(threshold_expt_dat$location)
threshold_expt_dat$aou4 <- as.factor(threshold_expt_dat$aou4)
threshold_expt_dat$ElevBin <- as.factor(threshold_expt_dat$ElevBin)

species.list <- focal_list
elev.bin <- unique(threshold_expt_dat$ElevBin)
spline.type <- "tp"
thresh.list <- unique(threshold_expt_dat$threshlab)
region<-"OLYM"

# Storage object
export.stats <- data.frame(
  species = NULL,
  elev.bin = NULL,
  threshold = NULL,
  metric = NULL,
  jday = NULL
)


# Empty list to store a list of models
mod.list <- list()

# Start plotting PDF
#pdf(file = paste0("./output/model_output/Threshold_Expt_SpeciesPlotsElev_",region,"_", spline.type, "_FocalSpp_=7_maxperc=0.3_pred_by_spp_", Sys.Date(),".pdf"), width = 8.5, height = 9.5)
tiff(file = paste0("./output/model_output/Threshold_Expt_SpeciesPlotsElev_",region,"_", spline.type, "_FocalSpp_=7_maxperc=0.3_pred_by_spp_", Sys.Date(),".tif"), width = 8.5, height = 9.5, units = "in", res = 300)
layout(matrix(1:8, nrow = 4, ncol = 2, byrow = FALSE))


for(i in 1:length(species.list)) {
  for(j in 1:length(thresh.list)) {
    # Model for species i, thershold j, Mid elevations
    data.ijk <- threshold_expt_dat[threshold_expt_dat$aou4 == tolower(species.list[i]) & threshold_expt_dat$ElevBin == elev.bin & threshold_expt_dat$threshlab==thresh.list[j], ]
    # Make sure we have enough data to continue
    if(dim(data.ijk)[1] > 0 & length(unique(data.ijk$location)) > 4) {
      # Fit GAM model
      model.ijk <- gam(cbind(nhits, ntrials) ~ s(JDay, bs = spline.type, k = 7) + s(location, bs = "re"),
        data = data.ijk,
        method = "REML",
        family = binomial,
        knots = list(JDay = c(50,300)) # make starting and ending knots way before and after breeding szn (50=20 Feb; 300=28 Oct)
      )
      # store model in mod.list
      mod.list[[paste(toupper(data.ijk$aou4[1]),data.ijk$ElevBin[1],sep="_")]] <- model.ijk
      
      # Extract phenometrics
      new_data <- data.frame(JDay = seq(summary(threshold_expt_dat$JDay[threshold_expt_dat$ElevBin==elev.bin & threshold_expt_dat$threshlab == thresh.list[j] & threshold_expt_dat$aou4==species.list[i]])[[1]], summary(threshold_expt_dat$JDay[threshold_expt_dat$ElevBin==elev.bin & threshold_expt_dat$threshlab == thresh.list[j]& threshold_expt_dat$aou4==species.list[i]])[[6]]), 
                             location = data.ijk$location[1]) # dummy location to satisfy predict()
      predictions <- predict(model.ijk, newdata = new_data,
                             type="link",
                             se.fit = T,
                             exclude = 's(location)')
      # standard errors
      critval <- 1.96 ## approx 95% CI
      
      predictions$fit.trans <- model.ijk$family$linkinv(predictions$fit)
      predictions$upr.trans <- model.ijk$family$linkinv(predictions$fit + (critval * predictions$se.fit))
      predictions$lwr.trans <- model.ijk$family$linkinv(predictions$fit - (critval * predictions$se.fit))
      
      # Finding mins / maxes
      pred.zoo <- as.zoo(predictions$fit.trans)
      min.days <- new_data$JDay[which(rollapply(pred.zoo, width = 9, FUN = function(x) which.min(x)==5)) + 4]
      max.days <- new_data$JDay[which(rollapply(pred.zoo, width = 9, FUN = function(x) which.max(x)==5)) + 4]
      
      # Test for a max at the start and/or end if starting/ending is > 30% of max probability elsewhere
      start.max <- as.vector(predictions$fit.trans[1] > 0.3*max(predictions$fit.trans))
      end.max <- as.vector(predictions$fit.trans[length(predictions$fit.trans)] > 0.3*max(predictions$fit.trans))
      # Eliminate max values that are not > 30% of the global max
      max.days <- max.days[predictions$fit.trans[new_data$JDay %in% max.days] > 0.3*max(predictions$fit.trans)]
      # vector including boundaries if applicable
      all.max <- c(if(start.max == T) {min(new_data$JDay)} else {NULL},
                   max.days,
                   if(end.max == T) {max(new_data$JDay)} else {NULL})
      # vector that only includes min if min is between maxes
      all.min <- min.days[min.days < max(all.max) & min.days > min(all.max)]
      # Only calculate half-max if there are non-boundary maxes  
      half.rise <- NULL
      half.fall <- NULL
      if(length(max.days) > 0) {
        for(k in 1:length(max.days)) {
          # calculate half-rise for every non-boundary max
          # calculated as half-y between max and preceding min
          max.k <- max.days[k]
          if(min(min.days) < min(max.k)) {
            min.k <- min.days[max(which(min.days < max.k))]
          } else {
            min.k <- min(new_data$JDay) # if min at previous edge, select edge
          }
          
          half.y <- mean(c(predictions$fit.trans[new_data$JDay == max.k],
                           predictions$fit.trans[new_data$JDay == min.k]))
          half.rise.k <- new_data$JDay[as.numeric(attr(which.min(abs((predictions$fit.trans - half.y)[which(new_data$JDay == min.k):which(new_data$JDay == max.k)])), "names"))]
          
          # calculate half-fall for every non-boundary max
          # calculated as half-y between max and following min
          if(max(min.days) > max.k) {
            min.k <- min.days[min(which(min.days > max.k))]
          } else {
            min.k <- max(new_data$JDay) # if min at next edge, select edge
          }
          half.y <- mean(c(predictions$fit.trans[new_data$JDay == max.k],
                           predictions$fit.trans[new_data$JDay == min.k]))
          half.fall.k <- new_data$JDay[as.numeric(attr(which.min(abs((predictions$fit.trans - half.y)[which(new_data$JDay == max.k):which(new_data$JDay == min.k)])), "names"))]
          
          # save in object
          half.rise <- c(half.rise, half.rise.k)
          half.fall <- c(half.fall, half.fall.k)
        }
      }
      
      # Export data.frame
      df.ijk <- rbind(if(length(min.days) > 0) {data.frame(species = toupper(species.list[i]),
                                                           threshold = thresh.list[j],
                                                           elev.bin = elev.bin,
                                                           metric = "Minimum", 
                                                           jday = min.days)},
                      if(start.max == TRUE) {data.frame(species = toupper(species.list[i]),
                                                        threshold = thresh.list[j],
                                                        elev.bin = elev.bin,
                                                        metric = "Start Boundary", 
                                                        jday = min(new_data$JDay))},
                      if(end.max == TRUE) {data.frame(species = toupper(species.list[i]),
                                                      threshold = thresh.list[j],
                                                      elev.bin = elev.bin,
                                                      metric = "End Boundary", 
                                                      jday = max(new_data$JDay))},
                      if(length(max.days) > 0) {rbind(
                        data.frame(species = toupper(species.list[i]),
                                   threshold = thresh.list[j],
                                   elev.bin = elev.bin,
                                   metric = "Maximum", 
                                   jday = max.days),
                        data.frame(species = toupper(species.list[i]),
                                   threshold = thresh.list[j],
                                   elev.bin = elev.bin,
                                   metric = "Half Rise", 
                                   jday = half.rise),
                        data.frame(species = toupper(species.list[i]),
                                   threshold = thresh.list[j],
                                   elev.bin = elev.bin,
                                   metric = "Half Fall", 
                                   jday = half.fall))
                      }
      )
      # Store this species/bin into export data.frame
      export.stats <- rbind(export.stats, df.ijk)
      
      data.summary <- data.ijk %>% group_by(JDay) %>% summarise(nhits =sum(nhits), ntrials= sum(ntrials), vocalrate = nhits/ntrials)
      
      # Create plot
      plot(new_data$JDay, predictions$fit.trans, type = "l", xlab = "Day of Year", ylab = "Predicted Vocal Rate",
           ylim = c(0, max(predictions$upr.trans)*1.1),
           xlim = range(threshold_expt_dat$JDay),
           main = paste(toupper(species.list[i]), "in Mid -", thresh.list[j], "threshold"), lwd = 2)
      lines(new_data$JDay, predictions$upr.trans)
      lines(new_data$JDay, predictions$lwr.trans)
      points(data.ijk$JDay, ((data.ijk$nhits / data.ijk$ntrials)/ max((data.ijk$nhits / data.ijk$ntrials)))*max(predictions$fit.trans), pch = 16, col = "#00000020")
      abline(v = all.max, col = "royalblue")
      abline(v = half.rise, col = "#edc948")
      abline(v = half.fall, col = "#e15759")
    }
  }
}

dev.off()

## Export CSV of metrics
#write.csv(x = export.stats, file = paste0("./output/model_output/Threshold_Expt_Phenometrics_",region,"_",spline.type, "_FocalSpp_k=7_maxperc=0.3_pred_by_spp_thresh=",thresh, "_", Sys.Date(), ".csv"))

#saveRDS(mod.list, file = paste0("./output/model_output/Threshold_Expt_mod.list", region,spline.type, "_FocalSpp_k=7_maxperc=0.3_pred_by_spp_", Sys.Date(), ".rds"))




ggplot(export.stats) +
  geom_point(aes(x=jday, y=threshold, color=species)) + 
  facet_grid(metric~species, scales = "free")

