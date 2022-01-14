library("tidyverse")
library("patchwork")

############
# Embryo AA
############

emb <- read.csv("00_data/embryo_metabolome.csv", stringsAsFactors = F, check.names = F) %>% 
  filter(class == "Amino Acid") %>% 
  select(- class) %>% 
  pivot_longer(- compound, names_to = "sample", values_to = "abundance") %>% 
  separate(col = sample, into = c("time","rep")) %>% 
  mutate(time = gsub(pattern = "E",replacement = "", x = time)) %>% 
  mutate(time = as.numeric(time)) %>% 
  group_by(compound, time) %>% 
  summarise(median_abundance = median(abundance)) %>% 
  ungroup() 
head(emb)  

p_emb <- ggplot(emb, aes(x = time, y = median_abundance, color = compound)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ compound, scales = "free_y") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,4,8,12,16,24)) +
  ggtitle("Amino acid profiles in the embryo") +
  labs(x = "Time after imbibition", y = "Normalised abundance in the embryo (AU)")
p_emb

############
# Endosperm AA
############
endo <- read.csv("00_data/endosperm_metabolome_germination.csv", stringsAsFactors = F, check.names = F) %>% 
  filter(class == "Amino Acid") %>% 
  select(- class) %>% 
  pivot_longer(- compound, names_to = "sample", values_to = "abundance") %>% 
  separate(col = sample, into = c("time","rep")) %>% 
  mutate(time = gsub(pattern = "A",replacement = "", x = time)) %>% 
  mutate(time = as.numeric(time)) %>% 
  group_by(compound, time) %>% 
  summarise(median_abundance = median(abundance)) %>% 
  ungroup() 
head(endo)  

p_endo <- 
  ggplot(endo, aes(x = time, y = median_abundance, color = compound)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ compound, scales = "free_y") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,4,8,12,16,24)) + 
  ggtitle("Amino acid profiles in the endosperm") +
  labs(x = "Time after imbibition", y = "Normalised abundance in the endosperm (AU)")
p_endo

############
# Save plots
############
p_emb + p_endo
ggsave(filename = "05_amino_acids_profiles/aa_profiles.pdf", width = 20, height = 7)
ggsave(filename = "05_amino_acids_profiles/aa_profiles.png", width = 20, height = 7)
