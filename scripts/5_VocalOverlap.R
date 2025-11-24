# This code generates visuals for singing period overlap, as well as summarizes vocal phenoperiods by elevation and migratory strategy.

# Run after script 4 so all dependencies are already loaded


# Load libraries
library(tidyverse)
library(mgcv)
library(gratia)
library(zoo) 
library(ggpubr)


cbFriendly <- c("#edc948","#e15759","#b07aa1","#004949", "#ff6db6", "#f28e2b", "#4e79a7")

# Run this file directly after script 4, or load in results from GAM run:
load("./output/model_output/multispecies_GAM_results.RData") # loads in loc_data, nhits.results, phenometrics ("export.stats"), and mod.results (table of GAM fit results)

pheno <- export.stats

# Read in bird table of species info
birdtable <- read.csv("./input/Siegel_et_al_2012_Elevation_Ranges_of_Birds_NCCN_non0s.csv")
head(birdtable)

bt2 <- birdtable %>% 
  separate(Species, into = c("Common Name", "Scientific Name"), sep="[()]") %>%
  mutate(Low.Fit = ifelse(Band_Code %in% pheno$species & Low.Bin==1, 1, 0), 
         Mid.Fit = ifelse(Band_Code %in% pheno$species & Mid.Bin==1, 1, 0))  %>%
  filter(Band_Code %in% toupper(species.list)) %>% 
  left_join(nhits.results, by=c("Band_Code" = "species"))

bt2$Mid.Bin[bt2$Band_Code=="WETA"] <- 1 # unexpected # of hits for this elevation; we will review the model
setdiff(toupper(bt2$Band_Code), export.stats$species) # make sure all species are in both tables (should read 0)
setdiff(export.stats$species,toupper(bt2$Band_Code))


### CREATE KEEP LIST - BIRDS WHOSE PHENOMETRICS SHOULD BE ANALYZED ###

# first, filter results to birds whose published ranges overlap the elev stratum definitions (mean +/- SD).
keeplist <- bt2 %>% select(Band_Code, Low.Bin, Mid.Bin) %>% mutate(species = Band_Code, Low=Low.Bin, Mid=Mid.Bin, .keep="unused") %>% pivot_longer(!species, names_to = "elev.bin", values_to = "Include")

pheno.join <- left_join(pheno, keeplist,by=c("species", "elev.bin"), relationship = "many-to-many") %>% left_join(nhits.results, by=c("species", "elev.bin"))

#n_distinct(verif_bind$common_name)
n_distinct(bin0.occ$aou4)
n_distinct(pheno.join$species) # does not include BTYW

head(pheno.join)

# then exclude data-deficient species: fewer than 1000 detections overall, resulted in no model fit (ruhu, swth, wewp)

nhits.results %>% filter(nhits < 1000)

pheno.join$Include[pheno.join$species=="RUHU" | pheno.join$species=="SWTH" | pheno.join$species=="WEWP"] <- 0
pheno.join$Include[pheno.join$species=="YRWA" & pheno.join$elev.bin=="Low"] <- 0
pheno.join$Include[pheno.join$species=="WAVI" & pheno.join$elev.bin=="Mid"] <- 0
pheno.join$Include[pheno.join$species=="WETA" & pheno.join$elev.bin=="Mid"] <- 1


# Who is vocalizing?------------------------------------------------------------------------------------------------------------------------------------

pheno_match <- pheno.join[pheno.join$Include==1,]

n_distinct(pheno_match$species)


# For each elevation bin, for each species, it is 'singing' between the first half-rise (or left boundary) and last half-fall (or right boundary)
# below code will generate warnings for spp x elev combos that cannot be estimated; that is OK
pheno.period <- data.frame(species = NULL, elev.bin = NULL, start.date = NULL, end.date = NULL)
for(i in 1:length(species.list)) {
  frame.i <- data.frame(species = toupper(species.list[i]), elev.bin = c("Low", "Mid"), start.date = NA, end.date = NA)
  pheno.i <- pheno_match[pheno_match$species == toupper(species.list[i]), ]
  frame.i[1, "start.date"] <- min(pheno.i[pheno.i$metric %in% c("Half Rise", "Start Boundary") & pheno.i$elev.bin == "Low", "jday"])
  frame.i[1, "end.date"] <- max(pheno.i[pheno.i$metric %in% c("Half Fall", "End Boundary") & pheno.i$elev.bin == "Low", "jday"])
  frame.i[2, "start.date"] <- min(pheno.i[pheno.i$metric %in% c("Half Rise", "Start Boundary") & pheno.i$elev.bin == "Mid", "jday"])
  frame.i[2, "end.date"] <- max(pheno.i[pheno.i$metric %in% c("Half Fall", "End Boundary") & pheno.i$elev.bin == "Mid", "jday"])
  pheno.period <- rbind(pheno.period, frame.i)
}

# If no pheno-phases or unable to assign period, eliminate
pheno.period <- pheno.period %>% filter(!grepl("Inf", start.date) & !grepl("Inf", end.date))

pheno.period

pheno.period %>% group_by(elev.bin) %>% summarise(nspecies = n_distinct(species))


# compare phenoperiods by spp, elev bins, migratory strategy... -----------

head(pheno.period)
unique(pheno.period$species)
pheno.period$duration <- pheno.period$end.date - pheno.period$start.date

pheno_table <- left_join(pheno.period, bt2, by=c("species"="Band_Code", "elev.bin"="elev.bin"))
#write.csv(pheno_table, "./output/pheno_table_25spp.csv") # generates Table S1

pheno_table$MigStrat <- factor(pheno_table$X, levels = c("R", "SDM", "LDM", "IRR"))
pheno_table$start.date.ymd <- as.Date(pheno_table$start.date, origin=as.Date("2021-01-01"))
pheno_table$end.date.ymd <- as.Date(pheno_table$end.date, origin=as.Date("2021-01-01"))

#pheno_table %>% group_by(MigStrat, elev.bin) %>% summarise(mean_dur = mean(duration), sd_dur = sd(duration)) %>% write_csv("./output/phenopd_by_elev_migrate_25spp.csv")

# peaks
# retrieve maxima from all species combos
# calculate difference between maxima for spp with >1 maximum

maxlist <- pheno_match %>% filter(metric=="Maximum" | metric=="Start Boundary") %>%
  group_by(species, elev.bin)

pheno_match %>% filter(metric=="Maximum"|metric=="Start Boundary") %>%
  group_by(species, elev.bin) %>% 
  summarise(diffs = max(jday)-min(jday)) %>%
  filter(diffs > 0) 

maxlist$Date <- as.Date(maxlist$jday, origin="2021-01-01")
maxlist <- left_join(maxlist, bt2, by=c("species"="Band_Code", "elev.bin"="elev.bin", "nhits"="nhits"))
maxlist$MigStrat <-ifelse(maxlist$X=="R", "Resident", 
                          ifelse(maxlist$X=="SDM", "Short-\ndistance", 
                                 ifelse(maxlist$X=="LDM", "Long-\ndistance", "Irruptive")))
maxlist$MigStrat <- factor(maxlist$MigStrat, levels = c("Resident", "Short-\ndistance", "Long-\ndistance", "Irruptive"))

pheno_table$MigStrat <- ifelse(pheno_table$X=="R", "Resident", 
                         ifelse(pheno_table$X=="SDM", "Short-\ndistance", 
                                ifelse(pheno_table$X=="LDM", "Long-\ndistance", "Irruptive")))
pheno_table$MigStrat <- factor(pheno_table$MigStrat, levels = c("Resident", "Short-\ndistance", "Long-\ndistance", "Irruptive"))

pheno_table$ElevLab <- ifelse(pheno_table$elev.bin=="Low", "LOW-ELEVATION", "MID-ELEVATION")
maxlist$ElevLab <- ifelse(maxlist$elev.bin=="Low", "LOW-ELEVATION", "MID-ELEVATION")

phenopd_spp_plot <- pheno_table %>% ggplot() +
  geom_linerange(aes(x=reorder(species, -start.date), ymin=start.date.ymd, ymax=end.date.ymd, color=MigStrat), alpha=0.7, linewidth=3) +
  geom_point(data=maxlist, aes(x=species, y=Date, group=MigStrat, shape=metric)) +
  coord_flip() +
  facet_grid(MigStrat~ElevLab, scales="free_y", space = "free_y") +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_rect(fill = "white"), panel.grid.major.x = element_line(linetype = 3), panel.grid.minor.x=element_blank(), panel.grid.major.y = element_line(linetype = 3)) +
  scale_color_manual(values = cbFriendly) +
  labs(color="Migratory Strategy", x="Species") +
  scale_y_date(date_labels = "%b %d", breaks = as.Date(c("2021-04-01","2021-05-01", "2021-06-01", "2021-07-01", "2021-08-01"))) 


ggsave(file=paste0("output/figures/25spp_phenopd_grouped_", Sys.Date(), ".tif"), plot = phenopd_spp_plot, width=6.5, height = 7, units = "in", dpi=600)

by_elevbin_migstrat <- pheno_table %>% group_by(elev.bin, ElevLab, MigStrat) %>% 
  summarise(n.species = n_distinct(species), mean.dur = mean(duration), sd.dur = sd(duration),
            mean.start = mean(start.date), sd.start = sd(start.date), 
            med.start = median(start.date), med.dur = median(duration))



# peaks

head(maxlist)

first_peaks <- maxlist %>% group_by(species, elev.bin, ElevLab) %>% summarise(first_peak = min(jday), first_peak_Date = as.Date(first_peak, origin="2021-01-01"))
pheno_table <- left_join(pheno_table, first_peaks)

pheno_table %>% group_by(X) %>% summarise(mean_peak1 = mean(first_peak))


# Fig 4 combined boxplot/dotplot ------------------------------------------

startplot.c <- ggplot(pheno_table) +
  geom_violin(outliers=FALSE,aes(x=MigStrat, y=start.date.ymd, fill=MigStrat), alpha=0.6) +
  geom_point(position=position_jitter(width = 0.2, height=0), aes(x=MigStrat, y=start.date.ymd, shape=MigStrat, fill=MigStrat), alpha=0.8, size=2, stroke=0.8) +
  #   geom_text(data=by_elevbin_migstrat, aes(label=paste("n =",n.species), x=MigStrat, y = 70, group=MigStrat)) +
  facet_wrap(~ElevLab, scales = "free_x") +
  theme_bw() +
  labs(x="", y="Start date of vocal phenoperiod", color="Migratory Strategy") +
  theme(legend.position = "none",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.major.y = element_line(linetype = 3)) +
  scale_color_manual(values = cbFriendly) +
  scale_fill_manual(values =  cbFriendly) +
  scale_shape_manual(values = c(22,23,24,25)) +
  scale_y_date(date_labels = "%b %d", breaks = as.Date(c("2021-04-01","2021-05-01", "2021-06-01", "2021-07-01", "2021-08-01"))) 


peak_plot.c <- ggplot(pheno_table) +
  geom_violin(outliers=FALSE,aes(x=MigStrat, y=first_peak_Date, fill=MigStrat), alpha=0.6) +
  geom_point(position=position_jitter(width = 0.2, height=0), aes(x=MigStrat, y=first_peak_Date, shape=MigStrat, fill=MigStrat), alpha=0.8, size=2, stroke=0.8) +
  # geom_text(data=by_elevbin_migstrat, aes(label=paste("n =",n.species), x=MigStrat, y = 70, group=MigStrat)) +
  facet_wrap(~ElevLab, ncol=2) +
  theme_bw() + 
  labs(y="Date of first peak", x="") +
  theme(legend.position = "none",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.major.y = element_line(linetype = 3)) +
  scale_color_manual(values = cbFriendly) +
  scale_fill_manual(values =  cbFriendly) +
  scale_shape_manual(values = c(22,23,24,25)) +
  scale_y_date(date_labels = "%b %d", breaks = as.Date(c("2021-04-01","2021-05-01", "2021-06-01", "2021-07-01", "2021-08-01"))) 

durplot.c <- ggplot(pheno_table) +
  geom_violin(aes(x=MigStrat, y=duration, fill=MigStrat), alpha=0.6) +
  geom_point(position=position_jitter(width = 0.2, height=0.1), aes(x=MigStrat, y=duration, shape=MigStrat, fill=MigStrat),  alpha=0.8, size=2, stroke=0.8) +
#  geom_text(data=by_elevbin_migstrat, aes(label=paste("n =",n.species), x=MigStrat, y = 10, group=MigStrat)) +
  facet_wrap(~ElevLab, ncol=2) +
  theme_bw() +
  labs(x="", y="Duration (days)", fill="Migratory Strategy") +
  theme(legend.position = "none",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x=element_blank(), panel.grid.major.y = element_line(linetype = 3)) +
  #  scale_color_manual(values = wes_palette("BottleRocket2", n=4)) +
  scale_fill_manual(values =  cbFriendly) +
  scale_shape_manual(values = c(22,23,24,25)) 

startpeakdur.hybrid <- ggarrange(startplot.c, peak_plot.c, durplot.c,
                            ncol=1, common.legend = T, legend = "none", align = "v") 
#  labels = c("a) Duration", "b) Start Day"))

pheno_table %>% group_by(X, MigStrat) %>% summarise(meandur = mean(duration), sddur = sd(duration))

ggsave(filename = paste0("./output/figures/start_peak_dur_migstrat_25spp_", Sys.Date(), ".tif"), plot = startpeakdur.hybrid,
       dpi=600, width = 6.5, height = 8.5, units = "in")  

# Double Peaks

dbl_peak <- pheno_match %>% filter(metric=="Maximum") %>%
  left_join(bt2,  by=c("species"="Band_Code", "elev.bin"="elev.bin")) %>%
  group_by(species, elev.bin, X) %>% 
  summarise(firstpeak = min(jday), secondpeak=max(jday), diffs = max(jday)-min(jday)) %>%
  filter(diffs > 0)

dbl_peak %>% group_by(elev.bin) %>% 
  summarise(nspec = n_distinct(species),
            meanfirst = mean(firstpeak),
            sdfirst = sd(firstpeak),
            meansecond = mean(secondpeak),
            sdsecond = sd(secondpeak),
            avgdiff = mean(diffs))

dbl_peak %>% group_by(X) %>% 
  summarise(nspec = n_distinct(species),
            meanfirst = mean(firstpeak),
            sdfirst = sd(firstpeak),
            meansecond = mean(secondpeak),
            sdsecond = sd(secondpeak),
            avgdiff = mean(diffs))

mean(dbl_peak$firstpeak)
mean(dbl_peak$secondpeak)

# half-decline

head(pheno)

declines <- pheno_match %>% filter(metric=="Half Fall")
declines %>% group_by(elev.bin, species) %>% arrange(desc(jday)) %>% slice(1:1) %>% group_by(elev.bin) %>% summarise(mean_fall = mean(jday))

