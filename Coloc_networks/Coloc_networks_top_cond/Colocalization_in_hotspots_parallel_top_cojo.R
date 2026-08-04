rm(list=ls(all=TRUE))
time.start <- Sys.time()

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

library(readr)
library(dplyr)
library(igraph)
library(ggplot2)
library(scales)
library(ggrepel)
library(Matrix)
library(ggraph)
library(tidygraph)
library(data.table)


#################################################
## 0. USER DEFINED PARAMETERS (MODIFY THIS PART)
################################################
colocalization_results_filtered_path ="Colocalization_results_filtered.csv"
hotspots_df_path = "Hotspots_df.csv"
LB_results_path = "LB_in_hotsposts.csv"
type_of_clustering = "connected" # either "greedy", "connected"
##############################################

#################################################
## Include parameter info in the logs
cat("Coloc results: ", colocalization_results_filtered_path,"\n")
cat("Hotspost df: ", hotspots_df_path,"\n")
cat("LB results in hotspot: ", LB_results_path,"\n")
#################################################


#############################################
###### 1. Upload all the data needed
#############################################
args <- commandArgs(trailingOnly = TRUE)
k = as.numeric(args[[1]])
colocalization_results_filtered <- fread(colocalization_results_filtered_path)
hotspot_df <- read_csv(hotspots_df_path)[,-1]
LB_in_hotsposts <- read_csv(LB_results_path)[,-1]

##############################################
##### 2.Computation
##########################################
cat("Processing row", k, "\n")
# Create a graph that has as node the seqid linked to a specific hotspot
sel_hotspot <- as.character(hotspot_df[k, 1])
chr_hotspot <- as.character(hotspot_df[k, 2])
LB_in_hotspost <- LB_in_hotsposts %>%
  filter(full_hotspot_gene_window == sel_hotspot & chr == chr_hotspot)
colocalization_results_filtered_in_hotspot <- colocalization_results_filtered %>%
  filter(full_hotspot_gene_window_locus_a == sel_hotspot & chr_locus_a == chr_hotspot &
           full_hotspot_gene_window_locus_b == sel_hotspot & chr_locus_b == chr_hotspot)

cat("Starting computation of the adjacency matrix \n")

colocalization_results_filtered_in_hotspot <- colocalization_results_filtered_in_hotspot %>%
  mutate(
    node_name_a = paste(trait_a, top_cond_a, sep=":"),
    node_name_b = paste(trait_b, top_cond_b, sep=":")
  )
nodes_name <- unique(c(colocalization_results_filtered_in_hotspot$node_name_a, colocalization_results_filtered_in_hotspot$node_name_b))
# Keep only the nodes that at least colocalize with another signal

node_index <- setNames(seq_along(nodes_name), nodes_name)
row_idx <- node_index[colocalization_results_filtered_in_hotspot$node_name_a]
col_idx <- node_index[colocalization_results_filtered_in_hotspot$node_name_b]

adj_matrix <- sparseMatrix(
  i = c(row_idx, col_idx),
  j = c(col_idx, row_idx),
  x = 1,
  dims = c(length(nodes_name), length(nodes_name)),
  dimnames = list(nodes_name, nodes_name)
)
adj_matrix@x[] <- 1

cat("Computation of the adjacency matrix finished \n")

g_est <- graph_from_adjacency_matrix(adj_matrix , mode = "undirected", diag = F,add.colnames = NA, add.rownames = NULL )
foldname <- paste(chr_hotspot,"_", sel_hotspot, sep ="")
dir.create(foldname)
filename <- paste(foldname, "/Estimated_graph.txt", sep="")
summary_adj_matrix <- summary(adj_matrix)
summary_adj_matrix$i_Name <- rownames(adj_matrix)[summary_adj_matrix$i]
summary_adj_matrix$j_Name <- colnames(adj_matrix)[summary_adj_matrix$j]
write.table(summary_adj_matrix, file=filename)
nodes_plot_label <- match(V(g_est)$name, nodes_name)
V(g_est)$plot_label <- nodes_plot_label
nod_legend <- data.frame(nodes_names = V(g_est)$name, label= match(V(g_est)$name, nodes_name))
filename <- paste(foldname, "/Node_legend.txt", sep="")
write.table(nod_legend, file=filename)

cat("Starting community computation \n")
cfg <- cluster_fast_greedy(g_est)
comp <- components(g_est)
cat("Community computation done \n")

#Annotate results
communities <-data.frame(nodes = V(g_est)$name,
                         group_conncetd_component = comp$membership,
                         group_greedy_clustering = cfg$membership )

communities$phenotype_id <- rep(NA, nrow(communities))
communities$chr <- rep(NA, nrow(communities))
communities$start <- rep(NA, nrow(communities))
communities$end <- rep(NA, nrow(communities))

for(i in 1: nrow(communities)){
  node_name <- communities$nodes[i]
  temp_a <- colocalization_results_filtered_in_hotspot %>% 
    filter(node_name_a == node_name) %>%
    select(trait_a, chr_locus_a, start_locus_a, end_locus_a ) %>% 
    rename(phenotype_id = trait_a, chr=chr_locus_a,
           start= start_locus_a, end = end_locus_a)
  temp_b <- colocalization_results_filtered_in_hotspot %>% 
    filter(node_name_b == node_name) %>%
    select(trait_b, chr_locus_b, start_locus_b, end_locus_b ) %>% 
    rename(phenotype_id = trait_b, chr=chr_locus_b,
           start= start_locus_b, end = end_locus_b)
  temp <- rbind(temp_a,temp_b)
  if(nrow(temp)>0){
    communities$phenotype_id[i] <- temp$phenotype_id[1]
    communities$chr[i] <- temp$chr[1]
    communities$start[i] <- temp$start[1]
    communities$end[i] <- temp$end[1]
  }
}

if(type_of_clustering == "connected"){
  communities$group = communities$group_conncetd_component
  V(g_est)$group <- communities$group_conncetd_component
} else {
  communities$group = communities$group_greedy_clustering
  V(g_est)$group <- communities$group_greedy_clustering
}


communities <-  full_join(LB_in_hotspost, communities, by= c("chr", "start","end", "phenotype_id"))


communities$community_type <-"all_trans"
for (group_id in unique(communities$group)) {
  temp <- communities %>% filter(group == group_id)
  
  if (nrow(temp) == 1) {
    communities$community_type[communities$group == group_id] <- temp$cis_or_trans
  } else {
    cis_nodes <- sum(temp$cis_or_trans == "cis")
    trans_nodes <- sum(temp$cis_or_trans == "trans")
    
    # Assign community types based on the conditions
    communities$community_type[communities$group == group_id] <- dplyr::case_when(
      cis_nodes == 0 ~ "all_trans",
      cis_nodes == 1 & trans_nodes >= 1 ~ "one_cis_multi_trans",
      cis_nodes > 1 & trans_nodes > 1 ~ "multi_cis_multi_trans",
      cis_nodes > 1 & trans_nodes == 0 ~ "multi_cis",
      TRUE ~ "unknown" 
    )
  }
}

filename <- paste(foldname, "/Hotspot_communities_",type_of_clustering,".csv", sep="")
write.csv(communities, file = filename)


cat("Saving the plots \n")
num_colors <- length((unique(V(g_est)$group)))
colrs <- scales::hue_pal()(num_colors)

weight.community=function(row,weigth.within,weight.between){
  if(as.numeric(V(g_est)$group[V(g_est)$name == row[1]])==as.numeric(V(g_est)$group[V(g_est)$name == row[2]])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}

E(g_est)$weight=apply(get.edgelist(g_est),1,weight.community,2,1)
V(g_est)$cis_trans <- NA
for(i in  1: length(V(g_est)$name)){
  temp_name <- V(g_est)$name[i]
  temp <- communities %>% filter(nodes == temp_name)
  if(nrow(temp)>0){
  V(g_est)$cis_trans[i] <- temp$cis_or_trans[1]}
}

layout=layout.fruchterman.reingold(g_est,weights=E(g_est)$weight)
vertex_shapes <- ifelse(V(g_est)$cis_trans == "cis", "circle", "csquare") 
plot_name <- paste(foldname, "/Plots_graph_",type_of_clustering,".jpg", sep="")

jpeg(plot_name, width = 485, height = 485)
plot(g_est,
     layout = layout, 
     edge.arrow.size = 0.1, 
     vertex.color = colrs[V(g_est)$group],  
     vertex.size = 5,
     vertex.shape = vertex_shapes,
     vertex.frame.color = colrs[V(g_est)$group],
     vertex.label=NA,
     vertex.label.color = "black",
     vertex.label.cex = 0.3,
     vertex.label.dist = 0,
     edge.curved = 0.2,
     edge.width = 0.2,
     margin = c(0, 0, 0, 0))

legend("bottomright",
       legend = c("Cis", "Trans"),
       pch = c(21, 22),
       pt.bg = c("black", "white"),
       pt.cex = 1.5,
       col = "black",
       bty = "n",
       title = "Node Shapes")
dev.off()



degrees <- degree(g_est)
degree_df <- data.frame(Degree = degrees)

ggplot(degree_df, aes(x = Degree)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  ggtitle("Node Degree Distribution") +
  xlab("Degree") +
  ylab("Frequency") +
  #theme_minimal()+
  theme_bw()
plot_name <- paste(foldname, "/Degree_plot_",type_of_clustering, ".png", sep="")
ggsave(plot_name, width = 25, height = 20, units = "cm")


mapping <- read.delim("somascan_tss_ncbi_grch37_version_20241212.txt", sep=",")
mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")

mapping <- mapping[, c("target","TSS", "UniProt_ID") ]
colnames(mapping) <- c("phenotype_id","phenotype_pos","UniProt_ID")

communities <- left_join(communities, mapping, by=c("phenotype_id", "UniProt_ID"), relationship = "many-to-many")

communities$phenotype_cis_end<-(communities$phenotype_pos+500000)
communities$phenotype_cis_start<-(communities$phenotype_pos-500000)

communities$phenotype_chr<-ifelse(communities$phenotype_chr=="X","23",communities$phenotype_chr)
communities$phenotype_chr<-ifelse(communities$phenotype_chr=="Y","23",communities$phenotype_chr)
table(communities$phenotype_chr)
communities$phenotype_chr<-as.numeric(communities$phenotype_chr)

data_map_cum <- communities %>%
  group_by(phenotype_chr) %>%
  dplyr::summarise(max_bp_map = max(phenotype_cis_start)) %>%
  mutate(bp_add_map = lag(cumsum(as.numeric(max_bp_map)), default = 0)) %>%
  select(phenotype_chr, bp_add_map)

communities <- communities %>%
  inner_join(data_map_cum, by = "phenotype_chr") %>%
  mutate(bp_cum_map = phenotype_cis_start+ bp_add_map)

save(communities, file=paste(foldname, "/communities_adam_plot_",type_of_clustering, ".RData", sep="") )

communities_adam_plot <- communities %>%
  select(POS, bp_cum_map, group, cis_or_trans, UniProt_ID)

ggplot(communities_adam_plot, aes(x = POS, y = bp_cum_map, color = as.factor(group), shape = cis_or_trans)) +
  geom_point(size = 3) + 
  geom_text_repel(
    data = subset(communities, cis_or_trans == "cis"),  # Only label "cis" points
    aes(label = UniProt_ID),                             
    size = 3,                                            
    max.overlaps = 80                                   
  ) +
  scale_shape_manual(values = c("cis" = 16, "trans" = 22)) +  # 16 = filled circle, 22 =  filled triangle
  labs(
    title = "Coding gene by pQTL with cis/trans",
    x = "pQTL positions",
    y = "Coding gene position",
    color = "Group",
    shape = "Cis or Trans"
  ) +
  theme_bw()

plot_name <- paste(foldname, "/Regional_plot_",type_of_clustering,".png", sep="")
ggsave(plot_name, width = 25, height = 20, units = "cm")

ccat('Compute genomics distance \n')
com <- as.numeric(names(table(communities$group)))
dist_matrix <- matrix(0, length(com) , length(com))
colnames(dist_matrix) <- com
rownames(dist_matrix) <- com
for(i in 1:(length(com)-1)){
  temp_t1 <- communities %>%
    filter(group == as.numeric(rownames(dist_matrix)[i])) %>% 
    mutate(SNP_pos =  as.numeric(sub("^\\d+:(\\d+):.*$", "\\1", SNPID)))
  median_pos_t1 <- median(temp_t1$SNP_pos,na.rm = TRUE)
  for (j in (i+1):length(com)){
    temp_t2 <- communities %>%
      filter(group == as.numeric(colnames(dist_matrix)[j])) %>% 
      mutate(SNP_pos = as.numeric(sub("^\\d+:(\\d+):.*$", "\\1", SNPID)))
    median_pos_t2 <- median(temp_t2$SNP_pos, na.rm = TRUE)
    dist_matrix[i,j]= abs(median_pos_t2-median_pos_t1)
    dist_matrix[j,i]=abs(median_pos_t2-median_pos_t1)
  }
}
filename <- paste(foldname, "/Gen_distance_communities_", type_of_clustering,".csv", sep="")
write.csv(dist_matrix, file = filename)

communities_to_retain <- communities_adam_plot %>%
  group_by(group) %>%
  summarise(n = n()) %>%
  filter(n > 4)

communities_to_retain <- na.omit(communities_to_retain)

communities_res <- communities_adam_plot %>%
  filter(group %in% communities_to_retain$group)

ggplot(communities_res, aes(x = POS, y = bp_cum_map, color = as.factor(group), shape = cis_or_trans)) +
  geom_point(size = 3) +  
  geom_text_repel(
    data = subset(communities_res, cis_or_trans == "cis"), 
    aes(label = UniProt_ID),                            
    size = 3,                                            
    max.overlaps = 80                                     
  ) +
  scale_shape_manual(values = c("cis" = 16, "trans" = 22)) + 
  labs(
    title = "Coding gene by pQTL with cis/trans",
    x = "pQTL positions",
    y = "Coding gene position",
    color = "Group",
    shape = "Cis or Trans"
  ) +
  #theme_minimal() + 
  theme_bw()
plot_name <- paste(foldname, "/Regional_plot_",type_of_clustering,"_restrict_to_large_communities.png", sep="")
ggsave(plot_name, width = 25, height = 20, units = "cm")
cat('All plots are saved')

cat('Computing the correlation between the communities in INTERVAL \n')
prot_res_df <- read_delim("INTERVAL_NonImp_residuals_final.txt", delim="\t")

communities_df <- communities %>% arrange(group)
seq_in_hot <- unique(communities_df$phenotype_id)
prot_res_df_in_hot <- prot_res_df[seq_in_hot]
corr_matrix <- cor(prot_res_df_in_hot)
com <- as.numeric(names(table(communities$group)))
median_corr_matrix <- matrix(0, length(com) , length(com))
colnames(median_corr_matrix) <- com
rownames(median_corr_matrix) <- com
for(i in 1:length(com)){
  temp_t1 <- communities_df %>%
    filter(group == rownames(median_corr_matrix)[i])
  for (j in i:length(com)){
    temp_t2 <- communities_df %>%
      filter(group == colnames(median_corr_matrix)[j])
    temp_corr_matrix <- corr_matrix[temp_t1$phenotype_id, temp_t2$phenotype_id]
    if(nrow(temp_t1)>1 & nrow(temp_t1)>1){
      median_corr_matrix[i,j] <- median(temp_corr_matrix[upper.tri(temp_corr_matrix)])
      median_corr_matrix[j,i] <- median(temp_corr_matrix[upper.tri(temp_corr_matrix)])}
    else{
      median_corr_matrix[i,j] <- median(temp_corr_matrix)
      median_corr_matrix[j,i] <- median(temp_corr_matrix)
    }
  }
}
filename <- paste(foldname, "/Corr_communities_",type_of_clustering,".csv", sep="")
write_csv(as.data.frame(median_corr_matrix), file = filename)

num_colors <- length(unique(V(g_est)$group))
palette <- scales::hue_pal()(num_colors)

tg <- as_tbl_graph(g_est)
layout <- create_layout(tg, layout = 'graphopt') 

p <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(aes(
  color = as.factor(V(g_est)$group),
  shape = as.factor(V(g_est)$cis_trans),
  fill = ifelse(V(g_est)$cis_trans == "cis", as.factor(V(g_est)$group), "white")
  ), size = 4) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = c(palette, white = "white")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) +
  guides(
    color = guide_legend(title = "Community Group"),
    shape = guide_legend(title = "Cis or Trans"),
    fill = "none"  
    )
ggsave(p, filename=paste(foldname, "/HighQuality_Network_",type_of_clustering, ".png", sep=""), width=8, height=8, dpi=300)

layout <- create_layout(tg, layout = 'stress') 
p <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(aes(
    color = as.factor(V(g_est)$group),
    shape = as.factor(V(g_est)$cis_trans),
    fill = ifelse(V(g_est)$cis_trans == "cis", as.factor(V(g_est)$group), "white")
  ), size = 4) +
  geom_node_text(aes(label = V(g_est)$plot_label), size = 3, repel = TRUE) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = c(palette, white = "white")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) +
  guides(
    color = guide_legend(title = "Community Group"),
    shape = guide_legend(title = "Cis or Trans"),
    fill = "none"  # Optional: hide fill legend if not needed 
  )
ggsave(p, filename=paste0(foldname, "/HighQuality_Network_investigation_",type_of_clustering,".png", sep=""), width=8, height=8, dpi=300)


cat("\n Computational time of: ")
cat( Sys.time() - time.start )
