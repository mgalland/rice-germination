library("tidyverse")
library("patchwork")
source("02_Figure_PCA/mypca.R")

############
# Embryo PCA
############

emb <- read.csv("00_data/embryo_transcriptome_germination.csv", stringsAsFactors = F) %>% 
  select(- locus, - gene_model, - annotation) %>% 
  column_to_rownames("probe") %>% 
  t(.)

emb_pca <- mypca(x = emb, center = T, scale = T)


### score plot from PC1 and PC2
scores <- emb_pca$scores[,1:2] %>% 
  rownames_to_column("sample") %>% 
  separate(sample, into = c("time", "rep"), sep = "_")

scores$time <- factor(scores$time, levels = c("E0", "E4", "E8", "E12", "E16","E24"))

# explained variance
explained_var = emb_pca$explained_var$exp_var[1:2]

p_emb <- ggplot(scores) + 
  geom_point(aes(x = PC1, y = PC2, col = time), size = 4) + 
  xlab(paste0('PC1(',explained_var[1],'%)')) + 
  ylab(paste0('PC2(',explained_var[2],'%)')) + 
  ggtitle('PCA score plot from the embryo transcriptome') +
  scale_color_brewer(type = "qual", palette = 2)
p_emb
ggsave(filename = "02_Figure_PCA/embryo_transcriptome.pdf", )


### Loadings
emb_loadings <- emb_pca$loadings %>% select(PC1, PC2) %>% rownames_to_column("probe")
emb_loadings_top10_PC1 <- 
  emb_loadings %>% 
  select(probe, PC1) %>% 
  mutate(absPC1 = abs(PC1)) %>% 
  top_n(n = 10, wt = absPC1)
emb_loadings_top10_PC2 <- 
  emb_loadings %>% 
  select(probe, PC2) %>% 
  mutate(absPC2 = abs(PC2)) %>% 
  top_n(n = 10, wt = absPC2)

probes2locus <- read.csv("00_data/embryo_transcriptome_germination.csv", stringsAsFactors = F) %>% 
  select(probe, locus, annotation)

emb_loadings_top10_PC1 %>% inner_join(., probes2locus, by = "probe") %>% write.csv(file = "02_Figure_PCA/emb_top10_PC1.csv", row.names = F)
emb_loadings_top10_PC2 %>% inner_join(., probes2locus, by = "probe") %>% write.csv(file = "02_Figure_PCA/emb_top10_PC2.csv", row.names = F)

##################
# 02 Endosperm PCA
############W######

endo <- read.csv("00_data/endosperm_transcriptome_germination.csv", stringsAsFactors = F) %>% 
  select(- locus, - gene_model, - annotation) %>% 
  column_to_rownames("probe") %>% 
  t(.)

endo_pca <- mypca(x = endo, center = T, scale = T)


### score plot from PC1 and PC2
scores <- endo_pca$scores[,1:2] %>% 
  rownames_to_column("sample") %>% 
  separate(sample, into = c("time", "rep"), sep = "_")

scores$time <- factor(scores$time, levels = c("A0", "A4", "A8", "A12", "A16","A24"))

# explained variance
explained_var = endo_pca$explained_var$exp_var[1:2]

p_endo <- ggplot(scores) + 
  geom_point(aes(x = PC1, y = PC2, col = time), size = 4) + 
  xlab(paste0('PC1(',explained_var[1],'%)')) + 
  ylab(paste0('PC2(',explained_var[2],'%)')) + 
  ggtitle('PCA score plot from the endosperm transcriptome') +
  scale_color_brewer(type = "qual", palette = 2)
p_endo
ggsave(filename = "02_Figure_PCA/endosperm_transcriptome.pdf", )


p_emb + p_endo
ggsave(filename = "02_Figure_PCA/Figure_02_pca_RNA.pdf", width = 16, height = 10)
ggsave(filename = "02_Figure_PCA/Figure_02_pca_RNA.svg", width = 16, height = 10)

### Endosperm loadings
endo_loadings <- endo_pca$loadings %>% select(PC1, PC2) %>% rownames_to_column("probe")
endo_loadings_top10_PC1 <- 
  endo_loadings %>% 
  select(probe, PC1) %>% 
  mutate(absPC1 = abs(PC1)) %>% 
  top_n(n = 10, wt = absPC1)
endo_loadings_top10_PC2 <- 
  endo_loadings %>% 
  select(probe, PC2) %>% 
  mutate(absPC2 = abs(PC2)) %>% 
  top_n(n = 10, wt = absPC2)

probes2locus <- read.csv("00_data/endosperm_transcriptome_germination.csv", stringsAsFactors = F) %>% 
  select(probe, locus, annotation)

endo_loadings_top10_PC1 %>% inner_join(., probes2locus, by = "probe") %>% write.csv(file = "02_Figure_PCA/endo_top10_PC1.csv", row.names = F)
endo_loadings_top10_PC2 %>% inner_join(., probes2locus, by = "probe") %>% write.csv(file = "02_Figure_PCA/endo_top10_PC2.csv", row.names = F)
