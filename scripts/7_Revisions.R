# revisions
library(tidyverse)
library(lme4)

# Does precision change seasonally?
load("input/OLYM_verifications_models_thresholds2025-11-07.RData")
hist(yday(all_verifs$captured_local_date))

loc_data$ElevBin[loc_data$ElevBin=="Med"] <- "Mid"

all_verifs

breakdown <- all_verifs %>% group_by(common_name, alternative_sound_class, focal_sound_class) %>% summarise(n())

raw_prec <- all_verifs %>% group_by(common_name) %>% summarise(tot = n(), tp = sum(is_species_present), fp = tot-tp, raw_precision = round(sum(is_species_present)/tot,3))
breakdown <- left_join(breakdown, raw_prec)

all_verifs$yday <- yday(all_verifs$captured_local_date)
all_verifs$species <- factor(all_verifs$common_name)


# does jday (1) explain variation in the data and (2) result in a better model fit?

all_verifs[is.na(all_verifs$aou4),] # check that they all got codes
all_verifs <- all_verifs %>% left_join(loc_data) %>% filter(ElevBin != "High") # add relevant env. vars

all_verifs %>% group_by(common_name, ElevBin) %>% summarise(n()) %>% View() # check to see how uneven


ggplot(all_verifs) +
  geom_histogram(aes(ElevM, fill=ElevBin))

# create a list of each species' data
list_df <- split(all_verifs, all_verifs$aou4)

# Initialize a list to store models
modlist_logit <- list()

# Iterate over each unique value of 'aou4'
for (aou4_val in names(list_df)) {
  
  subset_df <- list_df[[aou4_val]]
  
  model1 <- glm(is_species_present ~ logit, data = subset_df, family = binomial())
  
  modlist_logit[[aou4_val]] <- model1
}


modlist_jday <- list()

# Iterate over each unique value of 'aou4'
for (aou4_val in names(list_df)) {
  
  subset_df <- list_df[[aou4_val]]
  
  model1 <- glm(is_species_present ~ logit + yday, data = subset_df, family = binomial())
  
  modlist_jday[[aou4_val]] <- model1
}

# finally, elevation


modlist_elevbin <- list()

# Iterate over each unique value of 'aou4'
for (aou4_val in names(list_df)) {
  
  subset_df <- list_df[[aou4_val]]
  
  model1 <- glm(is_species_present ~ logit + as.factor(ElevBin), data = subset_df, family = binomial())
  
  modlist_elevbin[[aou4_val]] <- model1
}

 
# solvefor <- function(p, mod) {
#   beta0 = coef(mod)[1]
#   beta1 = coef(mod)[2]
#   logodds = (log(p/(1-p)) - beta0)/beta1
#   Cscore = plogis(logodds)
#   return(Cscore)
# }

# so we have 3 lists of models.
# Assuming modlist1 and modlist2 are your two model lists
sig_results <- data.frame(dataset = character(), p_logit = numeric(), p_jday = numeric(), p_elev = numeric(), n_low = numeric(), n_mid = numeric(), AIC_mod1 = numeric(),AIC_mod2 = numeric(),AIC_mod3 = numeric(),better_model = character(),stringsAsFactors = FALSE)

#jday first
for (name in names(modlist_logit)) {
  m1 <- modlist_logit[[name]]
  m2 <- modlist_jday[[name]]
  m3 <- modlist_elevbin[[name]]
  
  # Likelihood ratio test
  pval_logit <- summary(m1)$coefficients["logit", "Pr(>|z|)"] # p-val for logit (does logit explain variance?)
  test_jday <- anova(m1, m2, test = "Chisq")
  pval_jday <- test_jday$`Pr(>Chi)`[2]  # p-value for the added predictor
  test_elev <- anova(m1, m3, test = "Chisq")
  pval_elev <- test_elev$`Pr(>Chi)`[2]  # p-value for the added predictor
  n_low <- nrow(all_verifs[all_verifs$aou4==name & all_verifs$ElevBin=="Low",])
  n_mid <- nrow(all_verifs[all_verifs$aou4==name & all_verifs$ElevBin=="Mid",])
  # AICs
  aic1 <- AIC(m1)
  aic2 <- AIC(m2)
  aic3 <- AIC(m3)
  aics <- c(mod1 = aic1, mod2 = aic2, mod3 = aic3)
  sorted <- sort(aics)
  delta <- sorted[2] - sorted[1]
  
  best_model <- if (delta > 3) names(sorted)[1] else "none"

  
  sig_results <- rbind(sig_results, data.frame(dataset = toupper(name), p_logit = round(pval_logit,3), p_jday = round(pval_jday,3), p_elev = round(pval_elev,3), n_low = n_low, n_mid = n_mid, AIC_mod1 = aic1, AIC_mod2 = aic2, AIC_mod3 = aic3, better_model = best_model
  ))
}
sig_results

write_csv(sig_results, "./output/figures/revision_figs/covars_on_precision_compare.csv")

### PLOT PREDICTIONS ###

# create list of dataframes for each spp

data_list <- all_verifs %>%
  group_by(aou4) %>%
  group_split() %>%
  setNames(unique(all_verifs$aou4))


generate_predictions <- function(mod, df, model_type, dataset_name, factor_var) {
  # Create grid using raw jday
  score_seq <- seq(min(df$logit, na.rm = TRUE), max(df$logit, na.rm = TRUE), length.out = 100)
  jday_seq  <- seq(min(df$yday, na.rm = TRUE), max(df$yday, na.rm = TRUE), length.out = max(df$yday)-min(df$yday))
  # If factor_var is provided, get its levels
  if (!is.null(factor_var)) {
    factor_levels <- unique(df[[factor_var]])
    grid <- expand.grid(logit = score_seq, yday = jday_seq, factor_level = factor_levels)
    names(grid)[names(grid) == "factor_level"] <- factor_var  # rename to match model
  } else {
    grid <- expand.grid(logit = score_seq, yday = jday_seq)
  }
  
  
  # Predict directly using raw jday
  grid$predicted_prob <- predict(mod, newdata = grid, type = "response")
  
  # Annotate
  grid$model_type <- model_type
  grid$dataset <- dataset_name
  
  return(grid)
}


results_logit <- imap_dfr(modlist_logit, ~ generate_predictions(.x, data_list[[.y]], "Logit only", .y, factor_var = NULL))
results_logit_jday <- imap_dfr(modlist_jday, ~ generate_predictions(.x, data_list[[.y]], "Logit + JDay", .y, factor_var = NULL))
results_logit_elevbin <- imap_dfr(modlist_elevbin, ~ generate_predictions(.x, data_list[[.y]], "Logit + ElevBin", .y, factor_var = "ElevBin"))

# Combine
all_results <- bind_rows(results_logit, results_logit_jday, results_logit_elevbin)
head(df_long)
df_long$dataset <- df_long$species
df95 <- df_long[df_long$threshold==0.95,]

results_logit_jday <- left_join(results_logit_jday, df95, by="dataset")

ggplot(results_logit) +
  geom_raster(aes(x = logit, y = yday, fill = predicted_prob), alpha=0.8) +
  scale_fill_viridis_c(name = "Predicted\npr(TP)", option="D", direction=-1) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8)) +
  facet_wrap(~toupper(dataset)) +
  geom_vline(data=df95, aes(xintercept = logit_threshold)) +
  labs(title = "Predicted Probabilities Across Species, logit-only model", x = "Score", y = "Julian Day") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"))
#ggsave(filename = "output/figures/revision_figs/predicted_prtp_logit_only_allspp.png", width = 8.5, height=8, units="in", dpi=300)

ggplot(data = results_logit_jday) +
  geom_raster(aes(x = logit, y = yday, fill = predicted_prob), alpha=0.8) +
  facet_wrap(~ toupper(dataset) )+
  geom_vline(data=df95, aes(xintercept = logit_threshold)) +
  scale_fill_viridis_c(name = "Predicted\npr(TP)", direction = -1) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8)) +
  labs(x = "Score", y = "Day of Year") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"))

#ggsave(filename = "output/figures/revision_figs/predicted_prtp_logit_jday_allspp.tif", width = 8.5, height=8, units="in", dpi=300)

ggplot(data = results_logit_elevbin) +
  geom_raster(aes(x = logit, y = ElevBin, fill = predicted_prob), alpha=0.8) +
  facet_wrap(~ toupper(dataset) )+
  geom_vline(data=df95, aes(xintercept = logit_threshold)) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8)) +
  scale_fill_viridis_c(name = "Predicted\npr(TP)", direction = -1) +
  labs(x = "Score", y = "Day of Year") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"))
#ggsave(filename = "output/figures/revision_figs/predicted_prtp_logit_elev_allspp.tif", width = 8.5, height=8, units="in", dpi=300)

