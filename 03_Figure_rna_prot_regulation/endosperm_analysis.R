library("tidyverse")
library("pheatmap")
library("RColorBrewer")


#####################################
# Import endosperm transcriptome data
#####################################

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

rna_wide <- pivot_wider(rna_averaged, id_cols = "gene_model", names_from = "time", values_from = "avg_expr")


############################
# Import endosperm prot data 
############################

prot <- read.csv("03_Figure_rna_prot_regulation/109_endosperm_prot_down_0_vs_24HAI.csv", 
                    check.names = F, 
                    stringsAsFactors = F) %>% 
  rename("gene_model" = "protein") %>% 
  select(gene_model, log2ratio_A24vsA0)

############################
# Merge RNA and protein data
############################

rna_filtered <- inner_join(rna_wide, prot) 

rm(rna)
rm(rna_averaged)
rm(rna_wide)


#########################################
# Heatmap = endosperm protein down
#########################################

## Create matrix for down-regulated endosperm proteins
gene_mat_down_prot <- 
  rna_filtered %>% 
  column_to_rownames("gene_model") %>% 
  select(- log2ratio_A24vsA0) %>% 
  relocate(A4, .after = A0) %>% 
  relocate(A8, .after = A4)

my_sample_col <- data.frame(sample = rep(c("phase_I", "phase_II"), c(3,3)))
row.names(my_sample_col) <- colnames(gene_mat_down_prot)
my_colours <- list(sample = c("phase_I" = "#7fc97f", "phase_II" = "#beaed4"))

# gene annotation
my_gene_col <- 
  rna_filtered %>% 
  column_to_rownames("gene_model") %>% 
  select(log2ratio_A24vsA0) %>% 
  mutate(log2ratio_A24vsA0 = - log2ratio_A24vsA0) ## Takes the invert to display

my_colours <- list(sample = c("phase_I" = "#7fc97f", "phase_II" = "#beaed4"))

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
