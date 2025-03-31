
# loading the libraries 
if (!require ('pacman')){
  install.packages('pacman', dependencies = TRUE)
  library(pacman)
}


p_load(data.table,
       tidyverse,
       enrichR)

# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}
if (websiteLive) dbs <- listEnrichrDbs()

dbs <- c("Allen_Brain_Atlas_10x_scRNA_2021") #  the dataset 

# loading the network 
# This is the datatable mentioned in the supplementary as "Subnetworks Hierarchy"

subNetwork <- fread ("hierarchy_filter_GGE_ptv.tsv")

# splitting the network to the identified subnetworks 

list_subNetwork <- lapply(split(subNetwork, by = "represents"), function(x) unlist(strsplit(x$CD_MemberList, " ")))

clusters <- names(list_subNetwork)

Table <- tibble::tibble()

empty <- c()

for(cluster in clusters) {
  
  print (cluster)
  
  gene_list_input <- as.vector(list_subNetwork[[cluster]])
  
  if (websiteLive) {
    enriched <- enrichr(gene_list_input
                        , dbs)
  }
  
  Allen_Brain_Atlas <- if (websiteLive) enriched[["Allen_Brain_Atlas_10x_scRNA_2021"]]
  
  Allen_Brain_Atlas <- Allen_Brain_Atlas[order(Allen_Brain_Atlas$Adjusted.P.value),]
  
  Allen_Brain_Atlas <- Allen_Brain_Atlas %>%
    dplyr::mutate(dbs = "Allen_Brain_Atlas_10x_scRNA_2021") %>%
    dplyr::mutate(subnetwork = cluster)
  
  if (nrow(Allen_Brain_Atlas)> 0) {
    
    print("save")
    
    Table <- Table %>%
      dplyr::bind_rows(Allen_Brain_Atlas)
    
  } else {
    
    print("no enrichment")
    
    empty <- append(empty, cluster)
    
  }
}


write_tsv(Table, "Allen_Brain_Atlas_PTV.tsv",
          col_names = TRUE)


# Map the data to this taxonomy that has been provided in this paper "https://www.nature.com/articles/s41586-021-03465-8"

Data <- Table 

ID <- readxl::read_excel (base::paste("Supplementary Table 1.xlsx",
                                      sep = ""))
                          
ID <- ID %>%
  dplyr::select(`pCL_name (or CL_name)`,
                `pCL_id (or CL_id)`,
                neuron_type)

Data_cluster <- Data %>%
  dplyr::mutate(Overlap = sprintf('"%s"', Data$Overlap)) %>%
  filter(grepl('Human', Term)) %>%
  mutate(Term = gsub('Human', '', Term),  # removeing 'Human' from the terms
         matches = str_match(Term, "^\\s*(.*?)\\s+(\\S+)$"),  # extracting matches
         Cell = matches[, 2],  # Extract first part (Cell)
         Expression = matches[, 3]) %>%
  dplyr::inner_join(ID,
                    by = c("Cell" = "pCL_name (or CL_name)")) %>%
  dplyr::select(- matches,
                - Cell,
                - Expression)

write_tsv(Data_cluster, 
          "Human_Allen_Brain_Atlas_PTV.tsv",
          col_names = TRUE)
  

Data_cluster_sub <- Data_cluster %>%
  dplyr::group_by(subnetwork) %>%
  dplyr::filter (Adjusted.P.value < 0.05) %>%
  dplyr::filter (Odds.Ratio == max(Odds.Ratio)) %>%
  mutate(min_Adjusted_P_value = Term) %>%
  ungroup()

data <- as.data.frame(table(Data_cluster_sub$neuron_type 
                            ,useNA = "ifany"
                            ))

# percentages

data$fraction <- data$Freq / sum(data$Freq)

# the cumulative percentages (top of each rectangle)

data$ymax <- cumsum(data$fraction)

# the bottom of each rectangle

data$ymin <- c(0, head(data$ymax, n=-1))

# label position

data$labelPosition <- (data$ymax + data$ymin) / 2

# label

data$label <- paste0(data$Var1, "\n value: ", round(data$fraction * 100), "%")

# A donut chart

ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_hue() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  labs(title = base::paste("Neuronal Subtype-Specific Enrichment Analysis Using GGE Subnetworks:", 
                           "\n",
                           "Common and Protein-Truncating URVs"))+
  theme_void() +
  theme(legend.position = "none",
        title =element_text(size=18, face='bold')) -> fig


ggsave(paste("Allen_Brain_PTV",
             ".png",
             sep = ""),
       fig,
       device = "png",
       limitsize = FALSE,
       dpi = 300,
       width = 40,
       height = 40,
       units = "cm",
       bg = "white")
