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
  filter(tissue == "embryo")

##################################
# Import embryo transcriptome data
##################################

rna <- read.csv("00_data/embryo_transcriptome_germination.csv",
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
  relocate(E4, .after = E0) %>% 
  relocate(E8, .after = E4) %>% 
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

p229 <- df4plot %>% 
  filter(table == "S2B_229mRNAs") %>% 
  ggparcoord(., columns = 2:7, groupColumn = 8, scale = "globalminmax") + 
  facet_wrap(~cluster) + 
  my_theme() +
  labs(x = "Time after imbibition", y = "Normalised transcript abundance (AU)") +
  ggtitle("229 transcripts corresponding to 229 proteins up-accumulated in the embryo") 
p229

p133 <- df4plot %>% 
  filter(table == "S2D_133mRNAs") %>% 
  ggparcoord(., columns = 2:7, groupColumn = 8, scale = "globalminmax") + 
  facet_wrap(~cluster) +
  labs(x = "Time after imbibition", y = "Normalised transcript abundance (AU)") +
  ggtitle("133 transcripts corresponding to 133 proteins down-accumulated in the embryo") + 
  my_theme()
p133

p229 + p133 
ggsave("04_RNA_clusters/embryo_clusters.png", width = 20, height = 8)
ggsave("04_RNA_clusters/embryo_clusters.pdf", width = 20, height = 8)
ggsave("04_RNA_clusters/embryo_clusters.svg", width = 20, height = 8)





