# This code creates the GAM figures in Appendix for "Passive acoustic data as phenological distributions..." manuscript


library(tidyverse)
library(mgcv)
library(gratia)
library(zoo)


# LOOP FOR FINAL PLOTS ----------------------------------------------------
load("input/Zerofill_Sets_20251107.RData") # loads in zero-filled thresholded datasets

bin0 <- trials.0fill.095[!is.na(trials.0fill.095$threshold),]
thresh <- "0.95"
#bin0$ElevBin[bin0$ElevBin=="Med"] <- "Mid"
# mcgv doesn't like characters, convert to factors
bin0$location <- as.factor(bin0$location)
bin0$aou4 <- as.factor(bin0$aou4)
bin0$ElevBin <- factor(bin0$ElevBin, levels=c("Low","Mid"))

# truncate Mid data to day 120
bin0.occ <- bin0[!(bin0$ElevBin=="Mid" & bin0$JDay<120),] # ".occ" is now a meaningless name, but keeping out of laziness

# read in results from model run (script 4)
# model list
mod.list <- readRDS("output/model_output/mod.listOLYMtp_ALLDATA_k=7_maxperc=0.3_pred_by_spp_thresh=0.95_2025-11-04.rds") #or whatever your file is named

# calculated phenometrics
phenometrics <- read_csv("output/model_output/Phenometrics_OLYM_tp_ALLDATA_k=7_maxperc=0.3_pred_by_spp_thresh=0.95_2025-11-04.csv") %>% select(-1) %>% mutate(ElevBin = elev.bin)

phenometrics$key <- paste(phenometrics$species, phenometrics$elev.bin, sep="_")
phenometrics$metric <- factor(phenometrics$metric, levels = c("Start Boundary", "Half Rise", "Maximum", "Half Fall", "End Boundary", "Minimum"))

bin0.occ$aou4 <- toupper(bin0.occ$aou4)
species.list <- unique(bin0.occ$aou4)
elev.bins <- unique(bin0.occ$ElevBin)

# load in list of species x elev combos to keep
nhits.results <- trials.0fill.095 %>% left_join(loc_data) %>% filter(ElevBin !="High") %>% mutate(species=toupper(aou4), elev.bin=ElevBin, keep="unused") %>% group_by(species, elev.bin) %>% summarise(nhits = sum(nhits))



keeplist$key <- paste(keeplist$species, keeplist$elev.bin, sep="_")

color.key <- data.frame(metric = c("Start Boundary", "Half Rise", "Maximum", "Half Fall", "End Boundary"), color = c("#edc948", "#edc948", "royalblue", "#e15759", "#e15759"))
color.vec <- setNames(color.key$color, color.key$metric)



for (species in species.list) {
  plot.list <- list()
  metric.list <- list()
  dat.list <- list()
  
  for (elev.bin in elev.bins) {
    model.name <- paste(species, elev.bin, sep = "_")
    model.ijk <- mod.list[[model.name]]
    dat <- bin0.occ %>% filter(aou4 == species, ElevBin == elev.bin)
    phenoms <- phenometrics %>% filter(key == model.name) %>% filter(metric!="Minimum")
    
    # Generate prediction data
    new_data <- data.frame(
      JDay = seq(min(dat$JDay), max(dat$JDay)),
      location = dat$location[1]
    )
    
    predictions <- predict(model.ijk, newdata = new_data,
                           type = "link", se.fit = TRUE,
                           exclude = "s(location)")
    
    critval <- 1.96
    new_data <- new_data %>%
      mutate(
        fit = model.ijk$family$linkinv(predictions$fit),
        lwr = model.ijk$family$linkinv(predictions$fit - critval * predictions$se.fit),
        upr = model.ijk$family$linkinv(predictions$fit + critval * predictions$se.fit),
        ElevBin = elev.bin
      )
    
    # Raw data summary
    dat_summary <- dat %>%
      group_by(JDay) %>%
      summarise(rate = sum(nhits) / sum(ntrials), ntrials = sum(ntrials), .groups = "drop") %>%
      mutate(ElevBin = elev.bin)
    
    dat$rate.s <- ((dat$nhits/dat$ntrials) / max(dat$nhits/dat$ntrials)) * max(new_data$fit)
    
    # Combine for plotting
    plot.df <- left_join(new_data, dat_summary, by = c("JDay", "ElevBin")) %>%
      mutate(rate.s = (rate/max(rate)) * max(fit))
    
    plot.list[[elev.bin]] <- plot.df
    metric.list[[elev.bin]] <- phenoms
    dat.list[[elev.bin]] <- dat
  }
  
  plot.df.all <- bind_rows(plot.list)
  plot.metric.all <- bind_rows(metric.list)
  raw.all <- bind_rows(dat.list)
  
  # Check include flag
  keeplist.annot <- keeplist %>%
    filter(species == !!species) %>%
    mutate(
      JDay = -Inf,
      y = -Inf,
      label = ifelse(Include==0,"X",""),
      ElevBin = elev.bin
    )
  
  # Build plot
  p <- ggplot(plot.df.all, aes(x = JDay)) +
  #  geom_point(aes(y = rate), alpha = 0.3) +
    geom_point(data = raw.all, aes(y=rate.s, alpha=0.01), color="gray80") +
    geom_line(aes(y = fit), size=1) +
    geom_line(aes(y = lwr)) +
    geom_line(aes(y = upr)) +
 #   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_vline(data = plot.metric.all, aes(xintercept = jday, color = metric)) +
    geom_text(data = keeplist.annot, aes(x = JDay, y = y, label = label),
              inherit.aes = FALSE,
              color = "red", size = 20, alpha = 0.5,
              hjust = -0.5, vjust = -0.5)+
    facet_wrap(~ElevBin, nrow = 1, scales = "free_y") +
    scale_color_manual(values=color.vec) +
    labs(title = paste("Species:", species), y = "Predicted Daily Vocal Rate", x = "Day of Year") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"), panel.grid = element_blank())
  
 
  
  # add heatmap 
  
  q <- ggplot(plot.df.all, aes(x = JDay, y = 1, fill = ntrials)) +
    geom_tile() +
    facet_wrap(~ElevBin) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_void() +
    theme(strip.text = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "bottom", legend.key.width = unit(1.5, "cm"), legend.key.height =unit(0.2, "cm"))  +
    labs(fill = "Trials per Day")
  
  both <- ggarrange(p, q, ncol=1, heights = c(1,0.2), align="v")
  
  # Save plot
  ggsave(filename = paste0("./output/figures/revision_figs/gam_plots/", species, "_", Sys.Date(), ".tif"),
         plot = both, width = 8.5, height = 3.25, units = "in", dpi = 300)
}

