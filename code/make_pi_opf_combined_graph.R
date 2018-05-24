### Combined graph of select SCFA analysis
### picrust and metagenomes
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

# Store the genes associated with each pathway
butyrate_genes <- c("K00929", "K01034", "K01035")
propionate_genes <- c("K01908", "K01895", "K00932", "K19697")
acetate_genes <- c("K18372", "K00467", "K00156", "K00925", "K01512", "K01067", "K18118", "K01026", 
                   "K01905", "K22224", "K01895", "K00128", "K14085", "K00149", "K00138")


# Notes on what the butyrate kegg orthologs are:
  # K00929 - butyrate kinase
  # K01034 - acetate CoA (alpha)
  # K01035 - acetate CoA (beta)

# Load in needed metadata and rename sample column
picrust_data <- read_csv("data/process/tables/selected_scfa_gene_data.csv") %>% 
  filter(kegg_ortholog %in% butyrate_genes)

opf_data <- read_csv("data/process/tables/select_scfa_opf_data.csv") %>% 
  filter(kegg_id %in% butyrate_genes)

imputed_fig <- picrust_data %>% 
  mutate(dx = factor(dx, 
                levels = c("normal", "adenoma", "cancer"), 
                labels = c("Control", "Adenoma", "Carcinoma")), 
         kegg_ortholog = factor(kegg_ortholog, 
                                levels = c("K00929", "K01034", "K01035"), 
                                labels = c("Butyrate Kinase", "Acetate CoA\n(alpha subunit)", 
                                           "Acetate CoA\n(beta subunit)"))) %>% 
  ggplot(aes(kegg_ortholog, log10(rel_abund), fill = dx)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_vline(xintercept=seq(1.5, length(unique(picrust_data$kegg_ortholog))-0.5, 1), 
             lwd=1, colour="gray") + 
  theme_bw() + ggtitle("A") + 
  scale_fill_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  labs(x = "", y = expression(Log["10"]~Imputed~Relative~Abundance)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


metagenome_fig <- opf_data %>% 
  mutate(disease = factor(disease, 
                levels = c("Healthy", "Adenoma", "Cancer"), 
                labels = c("Control", "Adenoma", "Carcinoma")), 
         kegg_id = factor(kegg_id, 
                          levels = c("K00929"), labels = c("Butyrate Kinase"))) %>% 
  ggplot(aes(kegg_id, combined_count, fill = disease)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  theme_bw() + ggtitle("B") + 
  scale_fill_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  labs(x = "", y = "Counts per Million") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Combine the graphs together
gene_plot <- grid.arrange(imputed_fig, metagenome_fig, 
                          layout_matrix = rbind(c(1, 1, 2)))
# Write out to specific directory
ggsave("results/figures/butyrate_specific_genes_graph.pdf", gene_plot, width = 9, height = 5)

