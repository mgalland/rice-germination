library("tidyverse")
library("patchwork")
source("02_Figure_PCA/mypca.R")

############
# Embryo PCA
############

emb <- read.csv("00_data/embryo_metabolome.csv", stringsAsFactors = F, check.names = F) %>% 
  select(- class) %>% 
  column_to_rownames("compound") %>% 
  t(.)

emb_pca <- mypca(x = emb, center = T, scale = T)


### score plot from PC1 and PC2
scores <- emb_pca$scores[,1:2] %>% 
  rownames_to_column("sample") %>% 
  separate(sample, into = c("time", "rep"), sep = "-")

scores$time <- factor(scores$time, levels = c("E0", "E4", "E8", "E12", "E16","E24"))

# explained variance
explained_var = emb_pca$explained_var$exp_var[1:2]

p_emb <- ggplot(scores) + 
  geom_point(aes(x = PC1, y = PC2, col = time), size = 4) + 
  xlab(paste0('PC1(',explained_var[1],'%)')) + 
  ylab(paste0('PC2(',explained_var[2],'%)')) + 
  ggtitle('PCA score plot from the embryo metabolome') +
  scale_color_brewer(type = "qual", palette = 3)
p_emb
ggsave(filename = "02_Figure_PCA/embryo_metabolome.pdf", )


##################
# 02 Endosperm PCA
############W######

endo <- read.csv("00_data/endosperm_metabolome_germination.csv", stringsAsFactors = F, check.names = F) %>% 
  select(- class) %>% 
  column_to_rownames("compound") %>% 
  t(.)

endo_pca <- mypca(x = endo, center = T, scale = T)


### score plot from PC1 and PC2
scores <- endo_pca$scores[,1:2] %>% 
  rownames_to_column("sample") %>% 
  separate(sample, into = c("time", "rep"), sep = "-")

scores$time <- factor(scores$time, levels = c("A0", "A4", "A8", "A12", "A16","A24"))

# explained variance
explained_var = endo_pca$explained_var$exp_var[1:2]

p_endo <- ggplot(scores) + 
  geom_point(aes(x = PC1, y = PC2, col = time), size = 4) + 
  xlab(paste0('PC1(',explained_var[1],'%)')) + 
  ylab(paste0('PC2(',explained_var[2],'%)')) + 
  ggtitle('PCA score plot from the endosperm metabolome') +
  scale_color_brewer(type = "qual", palette = 2)
p_endo
ggsave(filename = "02_Figure_PCA/endosperm_metabolome.pdf", )

p_emb + p_endo
ggsave(filename = "02_Figure_PCA/Figure_02_pca_metabolites.pdf", width = 16, height = 10)
ggsave(filename = "02_Figure_PCA/Figure_02_pca_metabolites.svg", width = 16, height = 10)