library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("rlang")
library("matrixStats")


#######################
# heatmap color palette
#######################
my_palette <- colorRampPalette(colors = c("#fff5eb","#d94801"))

###########################
# Function to make heatmaps
###########################

create_heatmap <- function(selected_table = "S2F_12mRNAs") {
  # get genes of interest corresponding to clusters in suppl. table
  selected_table <- enquo(selected_table)
  genes_from_clusters <- filter(clusters, table == !!selected_table)
  
  # Retrieve scaled rna abundance for genes from clusters
  scaled_rna <- inner_join(x = rna_scaled, 
                           y = genes_from_clusters, 
                           by = "gene_model") %>% 
    arrange(cluster) 
  gene_mat <- scaled_rna %>% 
    column_to_rownames("gene_model") %>% 
    select(- table, - cluster, - tissue)

  my_gene_col <- genes_from_clusters %>% 
    select(- table, - tissue) %>% 
    column_to_rownames("gene_model")

  
  pheatmap(gene_mat, 
           color = my_palette(10),
           cluster_rows = F, 
           cluster_cols = FALSE, 
           fontsize_row = 5, 
           angle_col = 0,
           annotation_row = my_gene_col,  
           scale = "none")
}

############################
# Import cluster information
############################

clusters <- read.csv("03_Figure_rna_prot_regulation/clusters.csv",stringsAsFactors = FALSE) %>% 
  rename("gene_model" = "protein")

######################################
# Import endosperm transcriptome data
######################################

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

#####################
# Endosperm heatmaps
####################

create_heatmap(selected_table = "S2H_96mRNAs") # Figure 9
create_heatmap(selected_table = "S2F_12mRNAs") # Figure 10

