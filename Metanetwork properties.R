##### Properties of metanetwork
otu.list <- read.csv("otu_list_of_all_communities.csv")
### unidentified 

otu.list.2<- tidyr::unite(otu.list, "Class2", Phylum, Class, sep = "_")
otu.list.3<- tidyr::unite(otu.list.2, "Order", Class2, Order, sep = "_")
otu.list.4<- tidyr::unite(otu.list.3, "Family", Order, Family, sep = "_")
otu.list.5<- tidyr::unite(otu.list.4, "Genus", Family, Genus, sep = "_")

all.asso <- read.csv("Associations_of_all_networks_final_2.csv")

node.all.net <- read.csv('All_network_node_2.csv')


node.all.net <- read.csv('Metanetwork_node_V2.csv')
head(node.all.net)
as.character(all.asso$Source)[grepl("^F",as.character(all.asso$Source))]
class(all.asso$Source)
###igraph object
G0_cor_df_padj <- all.asso

nodeattrib_G0_combine <- data.frame(node=union(G0_cor_df_padj$Source,G0_cor_df_padj$Target))
nodeattrib_G0_combine$kingdom <- 0

for (i in as.character(nodeattrib_G0_combine$node))
{
  if (i %in% as.character(nodeattrib_G0_combine$node)[grepl("^B",as.character(nodeattrib_G0_combine$node))] == TRUE)
  {nodeattrib_G0_combine[nodeattrib_G0_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G0_combine[nodeattrib_G0_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G0_combine) <- as.character(nodeattrib_G0_combine$node)
nodeattrib_G0_combine$kingdom

head(nodeattrib_G0_combine)


all_G0_net <- graph_from_data_frame(G0_cor_df_padj,direct=F, vertices=nodeattrib_G0_combine)

## Number of nodes
length(V(all_G0_net)) # 80

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G0_net)))) #4423
length(grep("^F",names(V(all_G0_net)))) #2092


## Connections 
bb_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G0) #75129

ff_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G0) #44625

fb_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G0) #78994



## Network properties
meta_degree <- sort(igraph::degree(all_G0_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "2275"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "61.012"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G0_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  2.86590479879909"
print(paste("mean clustering coefficient = ", transitivity(all_G0_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.237793793145947"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G0_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality =  0.000285786149406026"
print(paste("mean closeness centrality = ", mean(closeness(all_G0_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0920390704424148"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G0_net, V(all_G0_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors =  61.0124328472755"
##

net <- all_G0_net
G0_all_deg <- igraph::degree(net,mode="all")
G0_all_betweenness <- betweenness(net, normalized = TRUE)
G0_all_closeness <- closeness(net, normalized = TRUE)
G0_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G0_all_transitivity)<- V(net)$name
G0_all_transitivity[is.na(G0_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
G0_all_deg.1percent <- G0_all_deg[G0_all_deg >= quantile(G0_all_deg,prob=1-n/100)]
length(G0_all_deg.1percent) #8

G0_all_betweenness.1percent <- G0_all_betweenness[G0_all_betweenness >= quantile(G0_all_betweenness,prob=1-n/100)]
length(G0_all_betweenness.1percent) #7

G0_all_closeness.1percent <- G0_all_closeness[G0_all_closeness >= quantile(G0_all_closeness,prob=1-n/100)]
length(G0_all_closeness.1percent) #7

intersect(names(G0_all_deg.1percent), names(G0_all_betweenness.1percent))
intersect(names(G0_all_deg.1percent), names(G0_all_closeness.1percent))
intersect(names(G0_all_betweenness.1percent), names(G0_all_closeness.1percent))

### network properties of bacteria and fungi
df.G0.degree<-data.frame(G0_all_deg)
head(df.G0.degree)

names(df.G0.degree)[1] <- c("Degree")

df.G0.closeness<-data.frame(G0_all_closeness)
head(df.G0.closeness)

names(df.G0.closeness)[1] <- c("Closeness")

df.G0.betweenness<-data.frame(G0_all_betweenness)
head(df.G0.betweenness)

names(df.G0.betweenness)[1] <- c("Betweenness")

df.G0.degree$kingdom <- ifelse(grepl("^B",rownames(df.G0.degree)),'Bacteria', 'Fungi')
df.G0.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G0.closeness)),'Bacteria', 'Fungi')
df.G0.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G0.betweenness)),'Bacteria', 'Fungi')
head(df.G0.betweenness)



#### Network hub
#### Cluster using Gephi
library(brainGraph)
library(igraph)
####Assigning community membership (module)
  ## Calling cluster defined in Gephi
  cluster_gephi <- node.all.net
  cluster_gephi <- cluster_gephi[c("Id","modularity_class")]
  rownames(cluster_gephi)<-cluster_gephi$Id
  cluster_gephi<-cluster_gephi[-c(1)]
  cluster_gephi$modularity_class <- cluster_gephi$modularity_class+1 ## Because of 0
  vec.cluster<-unlist(cluster_gephi, recursive = F, use.names = T)
  names(vec.cluster) <- rownames(cluster_gephi)
  
  ### Among-module connectivity
  Pi<-part_coeff(all_G0_net, vec.cluster)
  z_p_table<-data.frame(Pi)
  head(data.frame(Pi))
  ### within-module connectivity
  Zi<-within_module_deg_z_score(all_G0_net, vec.cluster)
  z_p_table$zi <- Zi
  
  
  ###Assigning Kingdom
  z_p_table$kingdom <- 0
  z_p_table$kingdom[grep("^B",rownames(z_p_table))] <- "Bacteria"
  z_p_table$kingdom[grep("^F",rownames(z_p_table))] <- "Fungi"
  
  z_p_table$Role <- 0
  z_p_table$Role[z_p_table$zi >= 2.5 & z_p_table$Pi >= 0.62] <- "Network hub"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi < 0.62] <- "Peripheral"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi > 0.62] <- "Connector"
  z_p_table$Role[z_p_table$zi > 2.5 & z_p_table$Pi < 0.62] <- "Module hub"
  
  write.csv(z_p_table, "z_p_table_gephi_metanetwork.csv")
  
  ##plotting
  p<-ggplot(z_p_table, aes(x=Pi, y=zi, color=kingdom)) +
    xlab('\n Among-module connectivity (Pi)')+
    ylab("Within-module connectivity (Zi) \n") +
    geom_point(size=3, alpha=0.7) +
    theme(aspect.ratio = 1)+
    # ggtitle("Volcano Plot \n") +
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size=FALSE) +
    #scale_x_continuous(breaks=seq(0,1,0.2))+
    #scale_y_continuous(breaks=seq(-20,0,-5))+
    scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
    geom_hline(yintercept=2.5, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', size = 0.75)
p

length(rownames(z_p_table)[which(z_p_table$Role == "Network hub")]) #135/6515
length(rownames(z_p_table)[which(z_p_table$Role == "Module hub")]) #39/6515
length(rownames(z_p_table)[which(z_p_table$Role == "Connector")]) #2766/6515
length(rownames(z_p_table)[which(z_p_table$Role == "Peripheral")]) #
