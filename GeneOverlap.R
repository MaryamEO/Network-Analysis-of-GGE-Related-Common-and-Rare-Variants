
# loading the libraries 
if (!require ('pacman')){
  install.packages('pacman', dependencies = TRUE)
  library(pacman)
}


p_load(readr,
       dplyr,
       ggplot2,
       data.table,
       tidyverse,
       lmtest,
       ggtext,
       patchwork,
       ggpubr,
       janitor,
       ggdist,
       gghalves,
       colorspace,
       scales,
       knitr,
       VarfromPDB,
       openxlsx,
       biomaRt,
       GeneOverlap,
       limma,
       clusterProfiler,
       org.Hs.eg.db,
       HGNChelper)


# loading the lasso genes and performing 
### https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.25/
### https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
### https://www.biostars.org/p/168203/
### retrieving the protein coding genes from Genome assembly GRCh37.p13


#subnetworks="/mnt/data/epi25_rare/data/epipgx_cogie_CENet/Gene_Set/C226_DMV.tsv"

subnetworks <- as.character(commandArgs(trailingOnly = TRUE)[1])

phenotype <-  as.character(commandArgs(trailingOnly = TRUE)[2])

# the gene IDs as the universe database 

genes <- fread ("/mnt/isilon/projects/isbsequencing/epi_rare/net_06032025/subnetworks/hgnc_complete_set.txt")

# Loading the Known Epilepsy genes

Epilepsy <- fread ("/mnt/data/epi25_rare/Network/Analysis/Epilepsy_Gene/genes4epilepsy/EpilepsyGenes_v2024-09.tsv")

if (phenotype == "Epilepsy") {
  
  print("Epilepsy")
  
  gene_list_input <- as.vector(Epilepsy$Gene)
  

  } else {
  
    print("GGE")
  
    GGE <- Epilepsy %>%
      dplyr::filter(str_detect(`Phenotype(s)`, "GGE"))
    
    gene_list_input <- as.vector(GGE$Gene)
  }


Gene_Sets <- fread (base::paste(subnetworks,
                                sep = ""),
                    header = TRUE)

currentHumanMap  <- getCurrentHumanMap()

subnetwork_list_input <- HGNChelper::checkGeneSymbols(Gene_Sets$gene, 
                                                      map = currentHumanMap )

subnetwork_list_input <- subnetwork_list_input %>%
  dplyr::mutate(Suggested.Symbol = ifelse (is.na(Suggested.Symbol), x, Suggested.Symbol))


# convergence across the databases 

Epilepsy_genes <- all ( gene_list_input %in% genes$symbol)

gene_list_input <- gene_list_input[gene_list_input %in% genes$symbol]

subnetwork_genes <- all (subnetwork_list_input$Suggested.Symbol %in% genes$symbol)

subnetwork_list_input <- subnetwork_list_input[subnetwork_list_input$Suggested.Symbol %in% genes$symbol, ]

# The hypergeometric distribution to measure the statistical significance 

go.obj <- newGeneOverlap(subnetwork_list_input$Suggested.Symbol,
                          gene_list_input,
                          genome.size=base::length(c(genes$symbol))
                         )
                         
go.obj <- testGeneOverlap(go.obj)
        
OR <- getOddsRatio(go.obj)

Overlapping_p_value <- go.obj@pval

Jaccard_Index <- go.obj@Jaccard

intersection_size <- base::length(go.obj@intersection)

phenotype_class <- phenotype

category_class <- unique(Gene_Sets$variant)

subnetwork_code <- unique(Gene_Sets$subnetwork)
 
Number_subnetwork_Gene <- base::length(go.obj@listA)
        
Number_Epilepsy_Gene <- base::length(go.obj@listB)
        

table <- data.frame(Epilepsy_genes,
                    subnetwork_genes,
                    OR,
                    Overlapping_p_value,
                    Jaccard_Index,
                    intersection_size,
                    phenotype_class,
                    category_class,
                    subnetwork_code,
                    Number_subnetwork_Gene,
                    Number_Epilepsy_Gene)
      
write_tsv(table, 
          base::paste("/mnt/isilon/projects/isbsequencing/epi_rare/net_06032025/subnetworks/GeneOverlap/Gene_overlap_",
                      phenotype_class,
                      "_",
                      category_class,
                      "_",
                      subnetwork_code,
                      ".tsv", 
                      sep = ""),
          col_names = TRUE)
        
  
    


# Heat map
# 
# gene_overlap <- fread ("/mnt/isilon/projects/isbsequencing/epi_rare/IRS/Table/Gene_overlap.tsv",
#                        header = TRUE)
# 
# 
# 
# 
# heatmap_plot <- ggplot(gene_overlap, aes(x = AF, y = gene_list_class, fill = -log10(Overlapping_p_value))) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   geom_text(aes(label = round(OR, 3)), 
#             color = "black", 
#             vjust = 0.5, 
#             hjust = 0.5,
#             size = 2,
#             angle = 90) +
#   labs(title = "Lasso Regression Retrieved Gene Set Enrichment Anlysis", x = "AF", y = "Gene Set") +
#   theme_void() +
#   theme(legend.position = "bottom",
#         #title = element_text(size = 25, colour = "black"),
#         #legend.text = element_text(size=25),
#         #legend.title = element_text(size=25),
#         axis.line = element_line(colour = "black"),
#         axis.text.y = element_text(size = 15, colour = "black"),
#         axis.text.x = element_text(size = 7, colour = "black", angle = 90),
#         axis.text.x.bottom = element_text(size = 7, colour = "black"),
#         axis.title.x = element_text(size = 15, colour = "black"),
#         axis.title.y = element_text(size = 15, colour = "black"),
#         strip.text.y = element_text(size = 15,
#                                     face = "bold"),
#         strip.text.x = element_text(size = 15,
#                                     face = "bold")) +
#   facet_grid(vars(phenotype_class), vars(category_class))
# 
# print(heatmap_plot)










  







        