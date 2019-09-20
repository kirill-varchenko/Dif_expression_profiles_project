library(KEGGgraph)
library(KEGG.db)
library(KEGGREST)
library(gage)
library(readxl)
# library(dplyr)
library(ggplot2)

### GLOBAL ###

use_local_kegg <- TRUE
local_kegg_dir <- 'kegg_data'

min_genes <- 5 # min number of genes in the pathway to cutoff

### FUNCTIONS ###

# parsing KEGG graphs
# using local files from local_kegg_dir when use_local_kegg is true
# tmpfiles otherwise
load_hsa <- function(hsa_name) {
  filehsa <- ifelse(use_local_kegg, paste0(local_kegg_dir, '/hsa', hsa_name, '.xml'), tempfile())
  if ( !(use_local_kegg & file.exists(filehsa)) ) {
    downloadable <- tryCatch(
      expr = {
        retrieveKGML(hsa_name, organism='hsa', destfile = filehsa, method='internal', quiet=TRUE)
        print(paste(hsa_name, 'downloaded to:', filehsa))
        TRUE
      },
      error = function(e) {
        print(paste0('Unable to download: hsa', hsa_name))
        FALSE
      },
      warning = function(e) {
        print(paste0('Unable to download: hsa', hsa_name))
        FALSE
      }
    )
    if (!downloadable) {
      if (file.exists(filehsa)) {
        file.remove(filehsa)
      }
      return(NA)
    }
  } else {
    print(filehsa)
  }
  mapkG <- parseKGML2Graph(filehsa, expandGenes=T, genesOnly=T)
  return(mapkG)
}

# main function that finds all descendants for each gene 
gene_iterator_forward <- function(goal_gene, 
                                  unique_symbol_names,
                                  forward_edges_list,
                                  diff_expressed_gene_symb){
  
  flag_visited <- rep(0,length(unique_symbol_names))
  names(flag_visited) <- unique_symbol_names
  
  all_children <- forward_edges_list[[goal_gene]]
  
  flag_visited[all_children] <- 1
  counter <- 1
  
  while (counter <= length(all_children)) {
    current_child <- all_children[counter]
    children_of_this_child <<- forward_edges_list[[current_child]]
    
    for_children_to_add <<- names(which(flag_visited[children_of_this_child] == 0))
    
    all_children <- c(all_children, for_children_to_add)
    
    flag_visited[for_children_to_add] <- 1
    counter <- counter + 1
    
  }
  #intersection between all gene descendants and differentially expressed ones
  all_children <- intersect(all_children, diff_expressed_gene_symb) 
  all_children <- union(all_children, goal_gene)
  return (all_children) #returns only differentially expressed descendants for each gene
}

### START HERE ###

data(kegg.gs)
data(kegg.gs.dise)

# the ability to choose between two datasets, depending on p_value
# you can download tast dataframe and read it from your computer (with your own path)
diff_expressed <- read_excel('gse_16357_001.xlsx') 
diff_expressed <- diff_expressed$gse16357_entrez
diff_expressed <- diff_expressed[-grep('/+', diff_expressed)]

all_KEGG_pathways <- append(kegg.gs, kegg.gs.dise) #browse pathways from KEGG database
list_intersected_genes <- lapply(all_KEGG_pathways, function(x) intersect(diff_expressed, x)) #intersection between genes in datasets and pathway genes from prev step
short_list <- list_intersected_genes[lengths(list_intersected_genes) >= min_genes] #cutoff by min 5 genes in the pathway
pathids <- substr(names(short_list), start=4, stop=8) #extracting IDs

#downloading pathways that contains genes from dataset
# load('graph_list.RData')
graph_list <- lapply(pathids, load_hsa)
graph_list <- graph_list[!is.na(graph_list)]
#save(graph_list, file='graph_list.RData')
#ok_graphs <- sapply(graph_list, function(x) typeof(x) == 'S4')
#graph_list <- graph_list[ok_graphs]
#mega_graph <- mergeKEGGgraphs(graph_list) #merge all pathways in one big graph


#for each graph-pathway in megagraph: change hsa IDs to UniProt IDs (aka 'real names'), 
#find all descendants for each gene in the given pathway, get the list of genes 
#with the largest number of differentially expressed descendants             
all_mega_ancestorz <- c()
for (next_graph in graph_list){
  print(next_graph)
  mega_nodes <- nodes(next_graph)
  mega_edges <- edges(next_graph)
  if( length(next_graph@edgeData@data) < length(mega_nodes)) {
    next
  }
  
  symbol_names <- sapply(mega_nodes, function(x) getKEGGnodeData(next_graph,x)@graphics@name)
  hsa_names <- sapply(mega_nodes, function(x) getKEGGnodeData(next_graph,x)@name)
  
  unique_symbol_names <- unique(unname(symbol_names))

  dict_hsa_to_symb = c()
  
  for(k in 1:length(symbol_names)){
    current_hsa <- sapply(hsa_names[[k]], function(x) substr(x,start = 5, stop = nchar(x)))
    current_symbol <- rep(symbol_names[k],length(current_hsa))
    names(current_symbol) <- current_hsa
    dict_hsa_to_symb <- c(dict_hsa_to_symb, current_symbol)
  }
  
  diff_expressed_gene_symb <- dict_hsa_to_symb[diff_expressed]
  diff_expressed_gene_symb<-unique(unname(diff_expressed_gene_symb[!is.na(diff_expressed_gene_symb)]))

  forward_edges_list <- vector('list', length(unique_symbol_names))
  names(forward_edges_list) <- unique_symbol_names
  
  counter <- 0
  for (node_name in attributes(mega_edges)$`names`){
    counter <- counter + 1
  
    disp_node_name <- getKEGGnodeData(next_graph,node_name)@graphics@name

    current_isoform_children <- unname(sapply(mega_edges[[node_name]], 
                                              function(x) getKEGGnodeData(next_graph,x)@graphics@name))
    
    if(length(current_isoform_children)!=0){
      forward_edges_list[[disp_node_name]] <- union(forward_edges_list[[disp_node_name]],
                                                    current_isoform_children)  
    }
  }
  forward_edges_list <- forward_edges_list[lapply(forward_edges_list,length)>0]
  
  all_children <- lapply(diff_expressed_gene_symb, function(goal_gene){
                                      gene_iterator_forward(goal_gene, 
                                                            unique_symbol_names,
                                                            forward_edges_list,
                                                            diff_expressed_gene_symb)})  
  
  names(all_children) <- diff_expressed_gene_symb
  
  pure_all_children <- all_children
  
  whomax <- names( which.max( sapply(all_children, length) ) )
  main_ancestorz <- whomax
  whomax_descendants <- length(all_children[[whomax]]) 
  while( whomax_descendants > 0 ){
    all_children <- sapply(all_children, function(x) setdiff(x,all_children[[whomax]]))
    whomax <- names( which.max( sapply(all_children, length) ) )
    main_ancestorz <- c(main_ancestorz, whomax)
    whomax_descendants <- length(all_children[[whomax]]) 
  }
  
  print(main_ancestorz)
  all_mega_ancestorz <- c(all_mega_ancestorz,main_ancestorz)
}
#all_mega_ancestorz <- unique(all_mega_ancestorz)

sorted_table <- sort(table(all_mega_ancestorz), decreasing = T) #get the most abundant genes (within all pathways)

#choose the top 5 most abundant genes from previous resuit (all_mega_ancestorz)
freq_anz <- as.data.frame(sorted_table[1:5])
colnames(freq_anz) <- c('Gene_name', 'n_path')

ggplot(freq_anz, aes(reorder(Gene_name, -n_path), n_path, fill = Gene_name)) + 
  geom_bar(stat = 'identity') + ggtitle('Most abundant genes') 


