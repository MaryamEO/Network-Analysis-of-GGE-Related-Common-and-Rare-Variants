
# loading the libraries 
if (!require ('pacman')){
  install.packages('pacman', dependencies = TRUE)
  library(pacman)
}


p_load(readr,
       igraph,
       dplyr,
       ggplot2,
       data.table,
       tidyverse,
       scales,
       knitr,
       orthogene,
       limma,
       enrichR,
       EnrichIntersect,
       ggraph,
       tidygraph
       ,signatureSearch
       ,networkD3,
       gplots,
       simplifyEnrichment
       )


set.seed(123)

# Mining the Drugable Genes for Personalized Medicine

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") 
}

if (websiteLive) {
  dbs <- listEnrichrDbs()
}

dbs <- c("DGIdb_Drug_Targets_2024") #  the dataset 

# Reading the file with the subnetworks information 

subNetwork <- fread ("hierarchy_filter_GGE_ptv.tsv")

# Making a list of genes within each subnetwork 

list_subNetwork <- lapply(split(subNetwork, by = "represents"), function(x) unlist(strsplit(x$CD_MemberList, " ")))

clusters <- names(list_subNetwork)

# Retrieving the lists of drugs and their corresponding go terms associated with the genes within each cluster 

data <- tibble::tibble()

for(cluster in clusters) {
  
  print (cluster)
  
  gene_list_input <- as.data.frame(list_subNetwork[[cluster]]) 
  
  # Making sure the name of genes are based on the latest version 
  
  gene_list_input <- gene_list_input %>%
    dplyr::mutate(Symbols = alias2SymbolTable(gene_list_input$`list_subNetwork[[cluster]]`)) %>%
    dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), `list_subNetwork[[cluster]]`, Symbols))
  
  # Mining the druggable genes 
  
  drugs <- tibble::tibble()
  
  if (websiteLive) {
    enriched <- enrichr(gene_list_input$Symbols
                        , dbs)
  }
  
  DGIdb <- if (websiteLive) enriched[["DGIdb_Drug_Targets_2024"]]
  
  DGIdb <- DGIdb[order(DGIdb$Adjusted.P.value),]
  
  DGIdb <- DGIdb %>%
    dplyr::mutate(dbs = "DGIdb_Drug_Targets_2024") %>%
    dplyr::mutate(subnetwork = cluster)
  
  drugs <- DGIdb
  
  # Listing the genes in the parent network

  if (cluster == "C201") {
    
    write_tsv(drugs,
              "Drug_network_PTV_C201.tsv",
              col_names = TRUE)
  }
  
  # Drugs 

  if (nrow(drugs) > 0) {
    
    # The Drug Set Enrichment Analysis (DSEA) with GSEA algorithm (dsea_GSEA function)
    
    drugList <- c(unique(drugs$Term))
    
    tryCatch( 
      {
        dsea_hyperG(
          drugList,
          type = "GO",
          ont = "BP",
          pvalueCutoff = 0.05,
          pAdjustMethod = "bonferroni",
          qvalueCutoff = 0.2,
          minGSSize = 10,
          maxGSSize = 500
          ) -> drug_GO 
        
        drug_go_result <- drug_GO@result
        
        drug_go_result <- drug_go_result %>%
          dplyr::mutate(subnetwork = cluster)
        data <- data %>%
          dplyr::bind_rows(drug_go_result)
        }, 
      error = function(e){
        print(e)
        }
      )
    }
  }
  
write_tsv(data,
          "Drug_GO_Cluster_PTV.tsv",
          col_names = TRUE)


  
# Clustering the semantic similarity matrix of the identified GO terms related to the set of drugs 

networkData <-  fread ("Drug_GO_Cluster_PTV.tsv")

parent_subnetwork <- networkData %>%
  dplyr::filter(subnetwork == "C201")

matrix = GO_similarity(unique(parent_subnetwork$ID))

cl = binary_cut(matrix)

export_to_shiny_app(matrix, cl)

# Network visulization of the drug-gene interactions

drugs <- fread ("Drug_network_PTV_C201.tsv")

drugs <- drugs %>% 
  tidyr::separate_rows(Genes, sep = ";")

gene <- length(unique(drugs$Genes))

count <- drugs%>%
  dplyr::group_by(Term) %>%
  dplyr::summarise(n = n())

drugs <- drugs %>%
  dplyr::mutate(Symbols = alias2SymbolTable(drugs$Genes)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), Genes, Symbols))

networkData <- drugs %>%
  dplyr::select(Symbols,
                Term)

drug.net <-  as_tbl_graph(networkData)

drug.net %>%
  activate(nodes) %>%
  mutate(centrality = centrality_alpha()) %>%
  mutate (degree = centrality_degree())-> drug_net

drug_net <- as_tibble(drug_net)


node_gene <- c(networkData$Symbols)
node_drug <- c(networkData$Term)

new_min <- 0.01
new_max <- 3

# common variant genes

common <- fread("sig_zscore_Common_2023_GGE.tsv")

common <- common %>%
  dplyr::mutate(Symbols = alias2SymbolTable(common$gene)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), gene, Symbols))

# rare variant genes

rare <- fread("sig_zscore_GGE_PTV.tsv")

rare <- rare %>%
  dplyr::mutate(Symbols = alias2SymbolTable(rare$gene)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), gene, Symbols))

nodeFactors <- factor(sort(unique(c(node_gene, node_drug))))

nodes <- data.frame(name = nodeFactors)

nodes <- nodes%>%
  dplyr::mutate(group = ifelse(nodes$name %in% node_gene, "Network Proximal Genes", "Drugs")) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% common$Symbols & !(nodes$name %in% rare$Symbols), "Common Variant Genes", group)) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% rare$Symbols & !(nodes$name %in% common$Symbols), "Damaging Missense URVs Genes", group)) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% rare$Symbols & nodes$name %in% common$Symbols, "Common and Damaging Missense URVs Genes", group)) %>%
  dplyr::mutate(group = as.factor(group)) %>%
  dplyr::inner_join(drug_net) %>%
  dplyr::mutate(scaled_sizes = ((degree * (new_max - new_min)) + new_min))

write_tsv(nodes,
          "nodes_PTV.tsv",
          col_names = TRUE)

  
  # Scaling centrality values to the desired size range

node_gene <- match(node_gene, levels(nodeFactors)) - 1

node_drug <- match(node_drug, levels(nodeFactors)) - 1

links <- data.frame(node_gene, node_drug)

links <- links%>%
  dplyr::mutate(n = 1) %>%
  dplyr::mutate(as.factor(n))

library(htmltools)
library(magrittr)
library(htmlwidgets)

forceNetwork(Links = links,
             Nodes = nodes,
             Source = 'node_gene',
             Target = 'node_drug',
             NodeID = 'name',
             Nodesize = "scaled_sizes",
             colourScale = JS("d3.scaleOrdinal().range(['#88db8b','#d6cec3', 'blue', 'red', '#bf00ff']);"),
             Value = "n",
             Group = 'group',
             zoom = TRUE,
             opacity = 1,
             charge = -2,
             fontSize = 20,
             legend = TRUE
             ) %>%
  onRender(  # Modify node styles to set text color
           "function(el, x) {
      d3.selectAll('.node text').style('fill', 'black');
    }"
  ) %>%
  htmlwidgets::prependContent(htmltools::tags$h1("Drug–Gene Interactions Due to Contribution of GGE-related Common Variants and Protein Truncating URVs")) %>%
  htmlwidgets::prependContent(htmltools::tags$footer("The size of the nodes corresponds to their degree")) %>%
  saveNetwork(file = 'GGE_Common_PTV.html')





