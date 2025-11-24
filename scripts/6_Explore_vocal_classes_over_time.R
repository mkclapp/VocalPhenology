# to be run after script 5
load("~/RProjects/Pheno_MS/input/OLYM_verifications_models_thresholds_20250204.RData")
head(all_verifs)
head(loc_data)


cbFriendly <- c("#edc948","#e15759","#b07aa1","#004949", "#ff6db6", "#f28e2b", "#4e79a7")

library(tidyverse)
# create minimal location data frame for merging ElevBins

locdat_min <- loc_data %>% select(location, ElevM, ElevBin)
locdat_min$ElevBin <- as.character(locdat_min$ElevBin)
locdat_min$ElevBin[locdat_min$ElevBin=="Med"] <- "Mid"

# pull in migratory strategy 
head(bt2)

# keeplist <- bt2 %>% select(`Common Name`, X, Band_Code, Low.Bin, Mid.Bin) %>% mutate(common_name = `Common Name`, species = Band_Code, Low=Low.Bin, Mid=Mid.Bin, .keep="unused") %>% pivot_longer(cols = c("Low", "Mid"), names_to = "ElevBin", values_to = "Include")
# #keeplist$common_name <- substr(keeplist$common_name,1,nchar(keeplist$common_name)-1)

activebirds <- bt2[bt2$Band_Code %in% pheno_match$species,]
activebirds$aou4 <- tolower(activebirds$Band_Code)

tot_vocs <- all_verifs %>% filter(is_species_present==1 & needs_review!=1) %>% 
  mutate(Song = ifelse(grepl("Song", focal_sound_class), "Song", "Call")) %>%
  filter(!is.na(Song)) %>%
  left_join(locdat_min) %>% 
  group_by(aou4, common_name, ElevBin, Week=week(as.Date(captured_local_date)), Week_date = as.Date(paste(2021, Week, "1", sep = "-"), format = "%Y-%W-%u")) %>%
  summarise(tot_vocs = n()) %>%
  left_join(activebirds, by=c("aou4"="aou4", "ElevBin"="elev.bin"))

songs_calls <- all_verifs %>% filter(is_species_present==1) %>% 
  mutate(Song = ifelse(grepl("Drum", focal_sound_class), "Song", ifelse(grepl("Song", focal_sound_class), "Song", "Call"))) %>%
  filter(!is.na(Song)) %>%
  left_join(locdat_min) %>% 
  group_by(common_name, ElevBin, Song, Week=week(as.Date(captured_local_date)),Week_date = as.Date(paste(2021, Week, "1", sep = "-"), format = "%Y-%W-%u")) %>% 
  summarise(n_vocs = n()) %>%
  left_join(tot_vocs) %>%
  mutate(prop_vocs = n_vocs/tot_vocs)


# migratory strategies individually
# Residents
songs_calls %>% filter(ElevBin !="High", X=="R") %>%
  ggplot() +
  geom_col(aes(x=Week_date, y=n_vocs, fill=Song), alpha=0.5) +
#  geom_line(aes(x=Week_date, y=n_vocs, color=Song)) +
  geom_vline(aes(xintercept=as.Date("2021-05-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-06-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-07-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-08-01"))) +
  facet_wrap(common_name~ElevBin, scales="free_y", ncol=4) +
  theme_pubclean() +
  scale_fill_manual(values=cbFriendly)  +
  labs(title = "Resident", x="Date", y="verified BirdNET samples", fill="Vocal Class")+
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), filename = "./output/figures/SuppFig_SoundClasses_OLYM_Resident.png",
       height = 8, width=6.5, units = "in", dpi=600)


# Migratory
songs_calls %>% filter(ElevBin !="High", X=="SDM" | X=="LDM") %>%
  ggplot() +
  geom_col(aes(x=Week_date, y=n_vocs, fill=Song), alpha=0.5) +
#  geom_line(aes(x=Week_date, y=n_vocs, color=Song)) +
  geom_vline(aes(xintercept=as.Date("2021-05-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-06-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-07-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-08-01"))) +
  facet_wrap(common_name~ElevBin, scales="free_y", ncol=4) +
  theme_pubclean() +
  scale_fill_manual(values=cbFriendly)  +
  labs(title = "Migratory", x="Date",y="verified BirdNET samples", fill="Vocal Class")+
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), filename = "./output/figures/SuppFig_SoundClasses_OLYM_Migratory.png",
       height = 8, width=6.5, units = "in", dpi=600)

# IRR
songs_calls %>% filter(ElevBin !="High", X=="IRR") %>%
  ggplot() +
  geom_point(aes(x=Week_date, y=n_vocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week_date, y=n_vocs, color=Song)) +
  geom_vline(aes(xintercept=as.Date("2021-05-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-06-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-07-01"))) +
  geom_vline(aes(xintercept=as.Date("2021-08-01"))) +
  facet_wrap(common_name~ElevBin, scales="free_y", ncol=2) +
  theme_bw() +
  scale_color_manual(values=cbFriendly)  +
  labs(title = "Irruptive", x="Date")
  
songs_calls %>% filter(ElevBin == "Low") %>%
ggplot() +
  geom_point(aes(x=Week_date, y=prop_vocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week_date, y=prop_vocs, color=Song)) +
  facet_wrap(~common_name, scales="free") +
  theme_bw() +
  scale_color_manual(values=cbFriendly) +
  labs(title = "Low Stratum")

songs_calls %>% filter(ElevBin == "Mid") %>%
  ggplot() +
  geom_point(aes(x=Week_date, y=prop_vocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week_date, y=prop_vocs, color=Song)) +
  facet_wrap(~common_name, scales="free") +
  theme_bw() +
  scale_color_manual(values=cbFriendly) +
  labs(title="Mid Stratum")

songs_calls %>% filter(ElevBin == "Low") %>%
  ggplot() +
  geom_point(aes(x=Week_date, y=n_vocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week_date, y=n_vocs, color=Song)) +
  facet_wrap(~common_name, scales="free") +
  theme_bw() +
  scale_color_manual(values=cbFriendly) +
  labs(title = "Low Stratum")



ggplot(songs_calls) +
  geom_point(aes(x=Week_date, y=n_vocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week_date, y=n_vocs, color=Song)) +
  facet_wrap(~common_name, scales="free") +
  theme_bw() +
  scale_color_manual(values=cbFriendly)


ggplot(vocal_class) +
  geom_point(aes(x=Week, y=nvocs, color=Song), alpha=0.5) +
  geom_line(aes(x=Week, y=nvocs, color=Song)) +
  facet_wrap(~common_name, scales="free") +
  theme_bw() +
  scale_color_manual(values=cbFriendly)



# Vocal Classes and CS ----------------------------------------------------

verifs4cs <- all_verifs %>% filter(is_species_present==1 & needs_review!=1) %>% 
  mutate(Song = ifelse(grepl("Drum", focal_sound_class), "Song", ifelse(grepl("Song", focal_sound_class), "Song", "Call"))) %>%
  filter(!is.na(Song)) %>%
  left_join(locdat_min) #%>%
 # left_join(keeplist)

ggplot(verifs4cs) +
  geom_histogram(aes(x=confidence, fill=Song), bins = 50, alpha=0.5) + 
  facet_wrap(~common_name, scales = "free_y", ncol=5) +
  scale_fill_manual(values=cbFriendly)  +
  theme_pubclean() +
  labs(x="Confidence Score",y="verified BirdNET samples", fill="Vocal Class")+
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), filename = "./output/figures/SuppFig_SongsCalls_CS.png",
       width = 11, height =8, units = "in", dpi=600)




songs_calls_all <- all_verifs %>%
  mutate(Song = ifelse(grepl("Drum", focal_sound_class), "Song", ifelse(grepl("Song", focal_sound_class), "Song", "Call"))) %>%
  filter(!is.na(Song)) %>%
  left_join(locdat_min) %>% 
  filter(ElevBin != "High") %>%
  mutate(Week=week(as.Date(captured_local_date)),Week_date = as.Date(paste(2021, Week, "1", sep = "-"), format = "%Y-%W-%u")) 

ggplot() +
  geom_histogram(data=songs_calls_all[songs_calls_all$is_species_present==1,], aes(x=Week_date), color="black", fill="forestgreen",bins = 50, alpha=0.4) +
  geom_histogram(data=songs_calls_all[songs_calls_all$is_species_present==0,], aes(x=Week_date), color="black", fill="red", bins = 50, alpha=0.4) +
  facet_grid(aou4~Song, scales = "free_y") +
  scale_fill_manual(values=cbFriendly)  +
  theme_pubclean() +
  labs(x="Confidence Score",y="verified BirdNET samples", fill="TP/FP", title="Migratory")+
  theme(legend.position = "bottom")

ggplot() +
  geom_histogram(data=songs_calls_all[songs_calls_all$is_species_present==1,], aes(x=Week_date), color="black", fill="#4e79a7",bins = 50, alpha=0.7) +
  geom_histogram(data=songs_calls_all[songs_calls_all$is_species_present==0,], aes(x=Week_date), color="black", fill="#f28e2b", bins = 50, alpha=0.7) +
  facet_wrap(~aou4, scales = "free_y", ncol=4) +
  scale_fill_manual(values=cbFriendly)  +
  theme_pubclean() +
  labs(x="Time of Season",y="verified BirdNET samples", fill="TP/FP", title="True and False Positive Rates by Season")+
  theme(legend.position = "bottom")

ggsave("./output/figures/revision_figs/tpfp_rates_season_portrait.png", width = 8.5, height = 11, units = "in", dpi=300)
# TP/FP by season, THRESHOLDED data

df_long
df_95 <- df_long %>% filter(threshold==0.95)

filtered_verifs <- songs_calls_all %>% left_join(df_95, by=c("aou4"="species")) %>%
  filter(confidence > cs_at_threshold)

ggplot() +
  geom_histogram(data=filtered_verifs[filtered_verifs$is_species_present==1,], aes(x=Week_date), color="black", fill="#4e79a7",bins = 50, alpha=0.7) +
  geom_histogram(data=filtered_verifs[filtered_verifs$is_species_present==0,], aes(x=Week_date), color="black", fill="#f28e2b", bins = 50, alpha=0.7) +
  facet_wrap(~aou4, scales = "free_y", ncol=4) +
#  scale_fill_manual(values=cbFriendly)  +
  theme_pubclean() +
  labs(x="Time of Season",y="verified BirdNET samples", fill="TP/FP", title="True and False Positive Rates by Season, 95% pr(TP)")+
  theme(legend.position = "bottom")

ggsave("./output/figures/revision_figs/tpfp_rates_season_prTP95_portrait.png", width = 8.5, height = 11, units = "in", dpi=300)

# we find that thresholding flattens the few false positive "blips" that happen in the unthresholded data. so thresholding is more likely to cause Type II error (failing to detect a pattern when it exists) rather than Type I error. This, combined with the fact that 
