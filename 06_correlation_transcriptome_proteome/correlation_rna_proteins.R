suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("patchwork"))

#####################
# Embryo RNA log 2 FC
#####################

emb_rna <- read.csv("00_data/embryo_transcriptome_germination.csv", 
                check.names = F, 
                stringsAsFactors = F) %>% 
  select(- probe, -locus, - annotation) %>% 
  ### Average value per gene (several probes per gene sometimes)
  pivot_longer(- gene_model, names_to = "sample", values_to = "expr") %>% 
  group_by(gene_model, sample) %>% 
  summarise(avg_expr = mean(expr)) %>% 
  ungroup() %>% 
  separate(sample, into = c("time","rep"), sep = "_") %>% 
  group_by(gene_model, time) %>% 
  summarise(avg_expr = mean(avg_expr)) %>% 
  mutate(time = factor(x = time, levels = c("E0","E4","E8","E12","E16","E24"))) %>% 
  ## Keep only 0 and 24 HAI and calculate log2 fold change
  filter(time == "E0" | time == "E24") %>% 
  pivot_wider(id_cols = gene_model, names_from = "time", values_from = "avg_expr") %>% 
  mutate(emb_rna_log2fc = log2(E24) - log2(E0)) %>% 
  rename("gene" = "gene_model") %>% 
  select(- E0, -E24)

##########################
# Embryo protein log 2 FC
#########################

emb_protein <- read.csv("00_data/embryo_proteins_germination.csv", stringsAsFactors = F) %>% 
  rename("gene" = "protein") %>% 
  select(- locus, - t_test)


################################
# Embryo RNA Protein correlation
################################
emb_cor <- inner_join(emb_rna, emb_protein, by = "gene") %>% 
  with(., cor(emb_rna_log2fc, emb_protein_log2fc, method = "spearman"))

N_genes_embryo <- inner_join(emb_rna, emb_protein, by = "gene") %>% nrow()

p_emb <- inner_join(emb_rna, emb_protein, by = "gene") %>% 
  ggplot(., aes(x = emb_rna_log2fc, y = emb_protein_log2fc)) +
  geom_point() + 
  xlim(-6, +6) +
  ylim(-6,+6) +
  theme_classic() + 
  labs(x = "Embryo mRNA log2 fold change (24 vs 0 HAI)", 
       y = "Embryo protein log2 fold change (24 vs 0 HAI)") +
  annotate(geom = "text", x = -3, y = 5, 
    label = paste("Spearman correlation coefficient:", as.character(round(emb_cor, 2)))) +
  annotate(geom = "text", x = -4, y = 4, 
           label = paste("Number of genes: ", as.character(N_genes_embryo))) +
  geom_abline(intercept = 0, slope = 1, colour = "blue") +
  ggtitle("Embryo")
p_emb

#####################
# Endosperm RNA log 2 FC
#####################

endo_rna <- read.csv("00_data/endosperm_transcriptome_germination.csv", 
                    check.names = F, 
                    stringsAsFactors = F) %>% 
  select(- probe, -locus, - annotation) %>% 
  ### Average value per gene (several probes per gene sometimes)
  pivot_longer(- gene_model, names_to = "sample", values_to = "expr") %>% 
  group_by(gene_model, sample) %>% 
  summarise(avg_expr = mean(expr)) %>% 
  ungroup() %>% 
  separate(sample, into = c("time","rep"), sep = "_") %>% 
  group_by(gene_model, time) %>% 
  summarise(avg_expr = mean(avg_expr)) %>% 
  ## Keep only 0 and 24 HAI and calculate log2 fold change
  filter(time == "A0" | time == "A24") %>% 
  pivot_wider(id_cols = gene_model, names_from = "time", values_from = "avg_expr") %>% 
  mutate(endosperm_rna_log2fc = log2(A24) - log2(A0)) %>% 
  rename("gene" = "gene_model") %>% 
  select(- A0, -A24)

##########################
# Endosperm protein log 2 FC
#########################

endosperm_protein <- read.csv("00_data/endosperm_proteins_germination.csv", 
                              stringsAsFactors = F) %>% 
  rename("gene" = "protein") %>% 
  select(- locus)

##################
# Endosperm correlation
###################

endosperm_cor <- inner_join(endo_rna, endosperm_protein, by = "gene") %>% 
  with(., cor(endosperm_rna_log2fc, endosperm_protein_log2fc, method = "spearman"))

N_genes_endosperm <- inner_join(endo_rna, endosperm_protein, by = "gene") %>% nrow()

p_endosperm <- inner_join(endo_rna, endosperm_protein, by = "gene") %>% 
  ggplot(., aes(x = endosperm_rna_log2fc, y = endosperm_protein_log2fc)) +
  geom_point() + 
  xlim(-2, +2) +
  ylim(-2,+2) +
  theme_classic() + 
  labs(x = "Endosperm mRNA log2 fold change (24 vs 0 HAI)", 
       y = "Endosperm protein log2 fold change (24 vs 0 HAI)") +
  annotate(geom = "text", x = -1.05, y = 2, 
           label = paste("Spearman correlation coefficient:", as.character(round(endosperm_cor, 2)))) +
  annotate(geom = "text", x = -1.5, y = 1.5, 
           label = paste("Number of genes: ", as.character(N_genes_endosperm))) +
  geom_abline(intercept = 0, slope = 1, colour = "orange") +
  ggtitle("Endosperm")
p_endosperm

##########################
# Final correlation figure
##########################

p_combined <- p_emb + p_endosperm
p_combined
ggsave(p_combined, 
       filename = "06_correlation_transcriptome_proteome/FigureSX_correlation_RNA_protein.pdf",
       width = 12, 
       height = 8)
