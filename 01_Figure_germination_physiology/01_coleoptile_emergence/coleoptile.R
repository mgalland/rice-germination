# Coleoptile emergence
library("tidyverse")

#############
# Import data
#############
df = read.csv("Figure1_germination_physiology/01_coleoptile_emergence/values_coleoptile_emergence.csv", 
              stringsAsFactors = F, 
              na.strings = "NA") %>% 
  as_tibble()


#################################################
# Calculate mean and SD of germination percentage
#################################################
n_germ <- 
  df %>%
  filter(phenotype == "germinated") %>% 
  pivot_longer(- c("petri_id", "phenotype"), names_to = "time", values_to = "value") %>% 
  group_by(time) %>% 
  summarise(mean_germ = mean(value, na.rm = TRUE), sd_germ = sd(value, na.rm = TRUE))

total <- 
  df %>%
  filter(phenotype == "total_seeds") %>% 
  pivot_longer(- c("petri_id", "phenotype"), names_to = "time", values_to = "value") %>% 
  group_by(time) %>% 
  summarise(mean_total = mean(value, na.rm = TRUE))

seed_perc_germ <- inner_join(n_germ, total) %>% 
  mutate(mean_percentage = mean_germ / mean_total * 100) %>% 
  mutate(sd_percentage = sd_germ / mean_total * 100)
seed_perc_germ

write.csv(x = seed_perc_germ, file = "Figure1_germination_physiology/01_coleoptile_emergence/mean_sd_values_coleoptile_emergence.csv", 
          row.names = F, 
          quote = F)

######
# Plot
######

ggplot(seed_perc_germ, aes(x= time, y = mean_percentage)) +
  geom_point() +
  geom_line(group = 1, size = 1, colour = "darkblue") +
  geom_errorbar(data = seed_perc_germ, 
                mapping = aes(ymin = mean_percentage - sd_percentage, 
                              ymax = mean_percentage + sd_percentage), 
                width = 0.5) +
  labs(x = "Time (h)", y = "Mean germination percentage (n = 7)")
ggsave(filename = "Figure1_germination_physiology/01_coleoptile_emergence/Fig1A_coleoptile_emergence.pdf")
ggsave(filename = "Figure1_germination_physiology/01_coleoptile_emergence/Fig1A_coleoptile_emergence.png")
ggsave(filename = "Figure1_germination_physiology/01_coleoptile_emergence/Fig1A_coleoptile_emergence.svg")


