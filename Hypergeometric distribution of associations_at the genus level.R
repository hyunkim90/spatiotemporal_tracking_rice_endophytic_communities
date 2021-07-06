#### Hypergeometric distribution of associations

otu.list.5 <- read.csv('otu_table_genus_merged.csv')
all.associate.net <- read.csv('Associations_of_all_networks_final_2.csv')
node.all.net <- read.csv('All_network_node_2.csv')

###Optimized
hypergeometric.analysis<-function(net.node, net.asso,keyword){

vertex.genus <- otu.list.5$Genus[which(otu.list.5$OTU_id %in% unique(net.node$id))]

G0_cor_df_padj <- net.asso

for (i in as.character(G0_cor_df_padj$Source)){
  G0_cor_df_padj$Source_Genus[which(G0_cor_df_padj$Source == i)] <- as.character(otu.list.5$Genus[which(otu.list.5$OTU_id ==i)])

}

G0_cor_df_padj$Target <- as.character(G0_cor_df_padj$Target)
G0_cor_df_padj$Target_Genus <- 0
for (i in as.character(G0_cor_df_padj$Target)){
  G0_cor_df_padj$Target_Genus[which(G0_cor_df_padj$Target == i)] <- as.character(otu.list.5$Genus[which(otu.list.5$OTU_id ==i)])

}

G0_cor_df_padj$Genus_link <- paste(G0_cor_df_padj$Source_Genus,"--",G0_cor_df_padj$Target_Genus)


overrep <- function(p1 = phyl1,
                    p2 = phyl2){
  if (p1 != p2){
    v1 = table(vertex.genus)[p1]
    v2 = table(vertex.genus)[p2]
    m <- v1*v2*(length(unique(net.node$id))-1)/length(unique(net.node$id))
    n <- length(unique(net.node$id))*(length(unique(net.node$id))-1) - m
    k <- length(G0_cor_df_padj$Source)
    x <- length(G0_cor_df_padj$Genus_link[which(G0_cor_df_padj$Genus_link == paste(p1,"--",p2))])
    p = phyper(x, m, n, k)}
  else {
    v1 = table(vertex.genus)[p1]
    m <- v1*c(v1-1)*(length(unique(net.node$id))-1)/length(unique(net.node$id))
    n <- length(unique(net.node$id))*(length(unique(net.node$id))-1) - m
    k <- length(G0_cor_df_padj$Source)
    x <- length(G0_cor_df_padj$Genus_link[which(G0_cor_df_padj$Genus_link == paste(p1,"--",p2))])
    p = phyper(x, m, n, k)
  }
}

p.over <- matrix(nrow = length(unique(vertex.genus)),ncol =  length(unique(vertex.genus)))
for(i in 1: length(unique(vertex.genus))){
  for(j in 1: length(unique(vertex.genus))){
    po <- overrep(p1 = names(sort(table(vertex.genus),decreasing = TRUE))[i],
                  p2 = names(sort(table(vertex.genus),decreasing = TRUE))[j])
    p.over[i,j] <- po
  }
  print(i)
}

p.over.adj <- matrix(p.adjust(p.over,method = "fdr"),nrow =  length(unique(vertex.genus)))
require(corrplot)

row.names(p.over.adj) <- names(sort(table(vertex.genus)[table(vertex.genus)>0],decreasing = TRUE))
colnames(p.over.adj) <- names(sort(table(vertex.genus)[table(vertex.genus)>0],decreasing = TRUE))

write.csv(p.over.adj,paste0("P value_Overrepresentation_meta_network_",keyword,".csv"))

return(p.over.adj)}


hypergeometric.analysis(net.node = node.all.net, net.asso = all.associate.net, "All_net")
