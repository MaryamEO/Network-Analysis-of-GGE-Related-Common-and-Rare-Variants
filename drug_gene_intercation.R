
# loading the libraries 
if (!require ('pacman')){
  install.packages('pacman', dependencies = TRUE)
  library(pacman)
}


p_load(data.table,
       tidyverse,
       openxlsx,
       limma,
       rDGIdb,
       tidygraph,
       networkD3,
       simplifyEnrichment,
       signatureSearch
       )


set.seed(123)

#  Mining the Drug able Genome for Personalized Medicine per gene in each subcluster ---- 

databases <- rDGIdb::sourceDatabases()


#  reading the file with the subnetworks information 
subNetwork <- fread ("hierarchy_filter_GGE_ptv.tsv")

# make a list of genes within each subnetwork 
list_subNetwork <- lapply(split(subNetwork, by = "represents"), function(x) unlist(strsplit(x$CD_MemberList, " ")))

clusters <- names(list_subNetwork)

data <- tibble::tibble()

for(cluster in clusters) {
  
  print (cluster)
  gene_list_input <- as.data.frame(list_subNetwork[[cluster]]) 
  
  # make sure the name of genes are based on the latest version 
  gene_list_input <- gene_list_input %>%
    dplyr::mutate(Symbols = alias2SymbolTable(gene_list_input$`list_subNetwork[[cluster]]`)) %>%
    dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), `list_subNetwork[[cluster]]`, Symbols))
  
  # Mining the Druggable Genome for Personalized Medicine per gene in each subcluster 
  drugs <- tibble::tibble()
  for (database in databases) {
    res <- rDGIdb::queryDGIdb(genes = gene_list_input$Symbols,
                              interactionTypes = NULL,
                              geneCategories = NULL,
                              sourceDatabases = database)
    
    drug_result <- detailedResults(res)
    drugs <- drugs %>%
      dplyr::bind_rows(drug_result)
    
  }
  
  if (cluster == "C314") {
    
    write_tsv(drugs,
              "Drug_network_PTV_C314.tsv",
              col_names = TRUE)
  }
  
# retrieving the lists of drugs and their corresponding go terms associated with the genes within each cluster 
# The Drug Set Enrichment Analysis (DSEA) with GSEA algorithm (dsea_GSEA function)

  if (nrow(drugs) > 0) {
  
    drugList <- c(unique(drugs$Drug))
    
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


  

# GO similarity matrix # ----

networkData <-  fread ("Drug_GO_Cluster_PTV.tsv")

test <- networkData %>%
  dplyr::filter(subnetwork == "C314")

# https://jokergoo.github.io/simplifyEnrichment/index.html

mat = GO_similarity(unique(test$ID))
cl = binary_cut(mat)
export_to_shiny_app(mat, cl)

# network visualization # ----

# https://stackoverflow.com/questions/42324064/directed-graph-using-networkd3-and-a-data-frame
# http://www.sthda.com/english/articles/33-social-network-analysis/136-network-analysis-and-manipulation-using-r/ 
# https://stackoverflow.com/questions/46899144/r-add-title-to-networkd3-plot-and-save
# https://tidygraph.data-imaginist.com/reference/centrality.html

drugs <- fread ("Drug_network_PTV_C314.tsv")

drugs <- drugs %>%
  dplyr::mutate(Symbols = alias2SymbolTable(drugs$Gene)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), Gene, Symbols))

networkData <- drugs %>%
  dplyr::select(Gene,
                Drug)

drug.net <-  as_tbl_graph(networkData)

drug.net %>%
  activate(nodes) %>%
  mutate(centrality = centrality_alpha()) %>%
  mutate (degree = centrality_degree())-> drug_net

drug_net <- as_tibble(drug_net)

node_gene <- c(networkData$Gene)
node_drug <- c(networkData$Drug)

# Define the desired new range
new_min <- 0.01
new_max <- 3

# common variant genes
common <- fread("GGE_common.tsv")

common <- common %>%
  dplyr::mutate(Symbols = alias2SymbolTable(common$gene)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), gene, Symbols))
                          
# rare variant genes
rare <- fread("GGE_PTV.tsv")
rare <- rare %>%
  dplyr::mutate(Symbols = alias2SymbolTable(rare$gene)) %>%
  dplyr::mutate (Symbols = dplyr::if_else(is.na(Symbols), gene, Symbols))

nodeFactors <- factor(sort(unique(c(node_gene, node_drug))))
                          
nodes <- data.frame(name = nodeFactors)
                          
nodes <- nodes%>%
  dplyr::mutate(group = ifelse(nodes$name %in% node_gene, "Network Proximal Genes", "Drugs")) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% common$Symbols & !(nodes$name %in% rare$Symbols), "Common Variant Genes", group)) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% rare$Symbols & !(nodes$name %in% common$Symbols), "Protein-Truncating URVs Genes", group)) %>%
  dplyr::mutate(group = ifelse(nodes$name %in% rare$Symbols & nodes$name %in% common$Symbols, "Common and Protein-Truncating URVs Genes", group)) %>%
  dplyr::mutate(group = as.factor(group)) %>%
  dplyr::inner_join(drug_net) %>%
  dplyr::mutate(scaled_sizes = ((degree * (new_max - new_min)) + new_min))

write_tsv(nodes,
          "nodes_PTV.tsv",
          col_names = TRUE)

# Scale centrality values to the desired size range

node_gene <- match(node_gene, levels(nodeFactors)) - 1
node_drug <- match(node_drug, levels(nodeFactors)) - 1

links <- data.frame(node_gene, node_drug)
                          
links <- links%>%
  dplyr::mutate(n = 1) %>%
  dplyr::mutate(as.factor(n))

# https://www.rdocumentation.org/packages/networkD3/versions/0.4/topics/forceNetwork
# https://d3js.org/d3-scale/ordinal#scaleOrdinal
# HTML generation
                          
library(htmltools)
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
  htmlwidgets::prependContent(htmltools::tags$h1("Drug–Gene Interactions Due to Contribution of GGE-related Common Variants and Protein-Truncating URVs ")) %>%
  htmlwidgets::prependContent(htmltools::tags$footer("The size of the nodes corresponds to their degree")) %>%
  saveNetwork(file = 'GGE_Common_PTV.html')
