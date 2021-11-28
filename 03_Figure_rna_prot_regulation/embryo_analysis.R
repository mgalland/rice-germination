library("tidyverse")
library("pheatmap")
library("RColorBrewer")


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

rna_wide <- pivot_wider(rna_averaged, id_cols = "gene_model", names_from = "time", values_from = "avg_expr")


#######################################
# Import embryo prot data (up proteins)
#######################################

prot_up <- read.csv("03_Figure_rna_prot_regulation/265_emb_prot_up_0_vs_24HAI.csv", 
                    check.names = F, 
                    stringsAsFactors = F) %>% 
  rename("gene_model" = "protein") %>% 
  mutate(profile = "up") %>% 
  mutate(t_test = as.numeric(t_test))

prot_down <- read.csv("03_Figure_rna_prot_regulation/152_emb_prot_down_0_vs_24HAI.csv", 
                    check.names = F, 
                    stringsAsFactors = F) %>% 
  rename("gene_model" = "protein") %>% 
  mutate(profile = "down")

prot <- bind_rows(prot_down, prot_up) %>% 
  select(gene_model, profile, log2ratio_E24vsE0)

############################
# Merge RNA and protein data
############################

rna_filtered <- inner_join(rna_wide, prot) 

rm(rna)
rm(rna_averaged)
rm(rna_wide)

#######################################
# Heatmap number 01 = embryo protein up
#######################################

## Create matrix for up-regulated embryo proteins
gene_mat <- 
  rna_filtered %>% 
  filter(profile == "up") %>% 
  column_to_rownames("gene_model") %>% 
  select(- profile, - log2ratio_E24vsE0) %>% 
  relocate(E4, .after = E0) %>% 
  relocate(E8, .after = E4)


# sample annotation
my_sample_col <- data.frame(sample = rep(c("phase_I", "phase_II"), c(3,3)))
row.names(my_sample_col) <- colnames(gene_mat)
my_colours <- list(sample = c("phase_I" = "#7fc97f", "phase_II" = "#beaed4"))

# gene annotation
my_gene_col <- 
  rna_filtered %>% 
  filter(profile == "up") %>% 
  column_to_rownames("gene_model") %>% 
  select(log2ratio_E24vsE0)

pheatmap(gene_mat, 
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         cluster_cols = FALSE, 
         fontsize_row = 5, 
         angle_col = 0, 
         scale = "none",
         annotation_col = my_sample_col,
         annotation_colors = my_colours,
         annotation_row = my_gene_col)

#########################################
# Heatmap number 02 = embryo protein down
#########################################

## Create matrix for down-regulated embryo proteins
gene_mat_down_prot <- 
  rna_filtered %>% 
  filter(profile == "down") %>% 
  column_to_rownames("gene_model") %>% 
  select(- profile, - log2ratio_E24vsE0) %>% 
  relocate(E4, .after = E0) %>% 
  relocate(E8, .after = E4)

# gene annotation
my_gene_col <- 
  rna_filtered %>% 
  filter(profile == "down") %>% 
  column_to_rownames("gene_model") %>% 
  select(log2ratio_E24vsE0) %>% 
  mutate(log2ratio_E24vsE0 = - log2ratio_E24vsE0) ## Takes the invert to display

my_colours <- list(sample = c("phase_I" = "#7fc97f", "phase_II" = "#beaed4"), 
                   log2ratio_E24vsE0 = brewer.pal(name = "Blues"))

pheatmap(gene_mat_down_prot, 
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         cluster_cols = FALSE, 
         fontsize_row = 5, 
         angle_col = 0, 
         scale = "none",
         annotation_col = my_sample_col,
         annotation_colors = my_colours,
         annotation_row = my_gene_col)
