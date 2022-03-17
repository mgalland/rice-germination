library("tidyverse")
library("RColorBrewer")
library("rlang")
library("matrixStats")
library("GGally")
library("patchwork")

###########################################
# Import cluster information for the embryo
##########################################

clusters <- read.csv("03_Figure_rna_prot_regulation/clusters.csv",stringsAsFactors = FALSE) %>% 
  rename("gene_model" = "protein") %>% 
  filter(tissue == "endosperm")

##################################
# Import embryo transcriptome data
##################################

rna <- read.csv("00_data/endosperm_transcriptome_germination.csv",
                check.names = F, 
                stringsAsFactors = F) %>% 
  select(- probe, -locus, - annotation) 

### Average value per gene (several probes per gene sometimes)
rna_averaged <- 
  rna %>% 
  pivot_longer(- gene_model, names_to = "sample", values_to = "expr") %>% 
  group_by(gene_model, sample) %>% 
  summarise(avg_expr = mean(expr)) %>% 
  ungroup() %>% 
  separate(sample, into = c("time","rep"), sep = "_") %>% 
  group_by(gene_model, time) %>% 
  summarise(avg_expr = mean(avg_expr))

rna_wide <- pivot_wider(rna_averaged, 
                        id_cols = "gene_model", 
                        names_from = "time", 
                        values_from = "avg_expr") %>% 
  column_to_rownames("gene_model") %>% 
  relocate(A4, .after = A0) %>% 
  relocate(A8, .after = A4) %>% 
  as.matrix()

## Scale matrix by maximum
gene_maximums <- rowMaxs(rna_wide)
rna_scaled <- rna_wide / gene_maximums
rna_scaled <- as.data.frame(rna_scaled) %>% rownames_to_column("gene_model")

### Clean up
rm(rna)
rm(rna_averaged)
rm(rna_wide)

####################
# Add cluster number
####################

df4plot <- inner_join(rna_scaled, clusters, by = "gene_model") 

#########################
# Parallel coordinate plot
#########################
my_theme <- function(base_size = 12){
  theme_bw(base_size = base_size)
}

p96 <- df4plot %>% 
  filter(table == "S2H_96mRNAs") %>% 
  ggparcoord(., columns = 2:7, groupColumn = 8, scale = "globalminmax") + 
  facet_wrap(~cluster) + 
  my_theme() +
  labs(x = "Time after imbibition", y = "Normalised transcript abundance (AU)") +
  ggtitle("96 transcripts corresponding to 96 proteins down-accumulated in the endosperm") 
p96

p11 <- df4plot %>% 
  filter(table == "S2F_12mRNAs") %>% 
  ggparcoord(., columns = 2:7, groupColumn = 8, scale = "globalminmax") + 
  facet_wrap(~cluster) +
  labs(x = "Time after imbibition", y = "Normalised transcript abundance (AU)") +
  ggtitle("10 transcripts corresponding to 10 proteins up-accumulated in the endosperm") + 
  my_theme()
p11

p96 + p11 
ggsave("04_RNA_clusters/endosperm_clusters.png", width = 20, height = 8)
ggsave("04_RNA_clusters/endosperm_clusters.pdf", width = 20, height = 8)
ggsave("04_RNA_clusters/endosperm_clusters.svg", width = 20, height = 8)





