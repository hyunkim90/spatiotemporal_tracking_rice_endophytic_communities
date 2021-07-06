#### Core jaccard distance
core.list.bac<-list(core.bac.75.BS, core.bac.75.RS,core.bac.75.R,
                   core.bac.75.S1, core.bac.75.S2,core.bac.75.S3,
                   core.bac.75.S4, core.bac.75.S5,core.bac.75.S6,
                   core.bac.75.S7, core.bac.75.S8,core.bac.75.S9,
                   core.bac.75.L1,core.bac.75.L2,core.bac.75.L3,
                   core.bac.75.FL,core.bac.75.G)

edgelist.trim.p <- core.list.bac

jac.mat <- NULL
for(i in 1:17){
  jac.i <- NULL
   for(j in 1:17){
    a <- table(edgelist.trim.p[[i]] %in% edgelist.trim.p[[j]])[2]
    b <- length(unique(c(edgelist.trim.p[[i]], edgelist.trim.p[[j]])))
    jac <- a/b
    jac.i <- c(jac.i, jac)
  }
  jac.mat <- rbind(jac.mat,jac.i) 
}
diag(jac.mat) <- 0

jac.mat

env.names_2 <- c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL","G")

row.names(jac.mat) <- env.names_2
colnames(jac.mat) <- env.names_2
jac.mat[is.na(jac.mat)] <- 0

jac.mat.pos<-melt(jac.mat)
write.csv(jac.mat.pos,"jac.mat.core.csv")





#### Hub


top_10_degree <- function(keyword){
  
  edge.file <- read.table(paste0('C:/Users/Hyun Kim/Desktop/SNU/Endophytes/Analysis/Analysis with 2018 samples/Network for SparCC/Stem network/FastSpar output/Input for gephi/0.3/1_0.3_edge_',keyword,'.tsv'), sep = '\t', header = T)
  edge.file$Source<- gsub(" ", ".", edge.file$Source)
  edge.file$Source <- gsub("-", ".", edge.file$Source)
  edge.file$Source <- gsub("\\(", ".", edge.file$Source)
  edge.file$Source <- gsub("\\)", ".", edge.file$Source)
  
  edge.file$Target<- gsub(" ", ".", edge.file$Target)
  edge.file$Target <- gsub("-", ".", edge.file$Target)
  edge.file$Target <- gsub("\\(", ".", edge.file$Target)
  edge.file$Target <- gsub("\\)", ".", edge.file$Target)
  
  
  
  G0_cor_df_padj <- edge.file
  
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
  
  net <- all_G0_net
  G0_all_deg <- igraph::degree(net,mode="all")
  
  df.G106_all_deg <- data.frame(G0_all_deg)
  df.G106_all_deg$Node <- rownames(df.G106_all_deg)
  df.G106_all_deg.top10 <- df.G106_all_deg %>% arrange(desc(G0_all_deg)) %>% head(n=10)
  df.G106_all_deg.top10$Network <- keyword
  names(df.G106_all_deg.top10)[1] <- "Degree"
  return(df.G106_all_deg.top10)
  
}

top.10.deg.BS.48<-top_10_degree("BS_48")
top.10.deg.BS.62<-top_10_degree("BS_62")
top.10.deg.BS.76<-top_10_degree("BS_76")
top.10.deg.BS.90<-top_10_degree("BS_90")
top.10.deg.BS.106<-top_10_degree("BS_106")
top.10.deg.BS.120<-top_10_degree("BS_120")
top.10.deg.BS.141<-top_10_degree("BS_141")

top.10.deg.RS.48<-top_10_degree("RS_48")
top.10.deg.RS.62<-top_10_degree("RS_62")
top.10.deg.RS.76<-top_10_degree("RS_76")
top.10.deg.RS.90<-top_10_degree("RS_90")
top.10.deg.RS.106<-top_10_degree("RS_106")
top.10.deg.RS.120<-top_10_degree("RS_120")
top.10.deg.RS.141<-top_10_degree("RS_141")

top.10.deg.R.48<-top_10_degree("R_48")
top.10.deg.R.62<-top_10_degree("R_62")
top.10.deg.R.76<-top_10_degree("R_76")
top.10.deg.R.90<-top_10_degree("R_90")
top.10.deg.R.106<-top_10_degree("R_106")
top.10.deg.R.120<-top_10_degree("R_120")
top.10.deg.R.141<-top_10_degree("R_141")

top.10.deg.S1.48<-top_10_degree("S1_48")
top.10.deg.S1.62<-top_10_degree("S1_62")
top.10.deg.S1.76<-top_10_degree("S1_76")
top.10.deg.S1.90<-top_10_degree("S1_90")
top.10.deg.S1.106<-top_10_degree("S1_106")
top.10.deg.S1.120<-top_10_degree("S1_120")
top.10.deg.S1.141<-top_10_degree("S1_141")

top.10.deg.S2.48<-top_10_degree("S2_48")
top.10.deg.S2.62<-top_10_degree("S2_62")
top.10.deg.S2.76<-top_10_degree("S2_76")
top.10.deg.S2.90<-top_10_degree("S2_90")
top.10.deg.S2.106<-top_10_degree("S2_106")
top.10.deg.S2.120<-top_10_degree("S2_120")
top.10.deg.S2.141<-top_10_degree("S2_141")

top.10.deg.S3.48<-top_10_degree("S3_48")
top.10.deg.S3.62<-top_10_degree("S3_62")
top.10.deg.S3.76<-top_10_degree("S3_76")
top.10.deg.S3.90<-top_10_degree("S3_90")
top.10.deg.S3.106<-top_10_degree("S3_106")
top.10.deg.S3.120<-top_10_degree("S3_120")
top.10.deg.S3.141<-top_10_degree("S3_141")

top.10.deg.S4.62<-top_10_degree("S4_62")
top.10.deg.S4.76<-top_10_degree("S4_76")
top.10.deg.S4.90<-top_10_degree("S4_90")
top.10.deg.S4.106<-top_10_degree("S4_106")
top.10.deg.S4.120<-top_10_degree("S4_120")
top.10.deg.S4.141<-top_10_degree("S4_141")

top.10.deg.S5.76<-top_10_degree("S5_76")
top.10.deg.S5.90<-top_10_degree("S5_90")
top.10.deg.S5.106<-top_10_degree("S5_106")
top.10.deg.S5.120<-top_10_degree("S5_120")
top.10.deg.S5.141<-top_10_degree("S5_141")

top.10.deg.S6.76<-top_10_degree("S6_76")
top.10.deg.S6.90<-top_10_degree("S6_90")
top.10.deg.S6.106<-top_10_degree("S6_106")
top.10.deg.S6.120<-top_10_degree("S6_120")
top.10.deg.S6.141<-top_10_degree("S6_141")

top.10.deg.S7.90<-top_10_degree("S7_90")
top.10.deg.S7.106<-top_10_degree("S7_106")
top.10.deg.S7.120<-top_10_degree("S7_120")
top.10.deg.S7.141<-top_10_degree("S7_141")

top.10.deg.S8.90<-top_10_degree("S8_90")
top.10.deg.S8.106<-top_10_degree("S8_106")
top.10.deg.S8.120<-top_10_degree("S8_120")
top.10.deg.S8.141<-top_10_degree("S8_141")

top.10.deg.S9.90<-top_10_degree("S9_90")
top.10.deg.S9.106<-top_10_degree("S9_106")
top.10.deg.S9.120<-top_10_degree("S9_120")
top.10.deg.S9.141<-top_10_degree("S9_141")

top.10.deg.L1.48<-top_10_degree("L1_48")
top.10.deg.L1.62<-top_10_degree("L1_62")
top.10.deg.L1.76<-top_10_degree("L1_76")
top.10.deg.L1.90<-top_10_degree("L1_90")
top.10.deg.L1.106<-top_10_degree("L1_106")
top.10.deg.L1.120<-top_10_degree("L1_120")
top.10.deg.L1.141<-top_10_degree("L1_141")

top.10.deg.L2.48<-top_10_degree("L2_48")
top.10.deg.L2.62<-top_10_degree("L2_62")
top.10.deg.L2.76<-top_10_degree("L2_76")
top.10.deg.L2.90<-top_10_degree("L2_90")
top.10.deg.L2.106<-top_10_degree("L2_106")
top.10.deg.L2.120<-top_10_degree("L2_120")
top.10.deg.L2.141<-top_10_degree("L2_141")

top.10.deg.L3.62<-top_10_degree("L3_62")
top.10.deg.L3.76<-top_10_degree("L3_76")
top.10.deg.L3.90<-top_10_degree("L3_90")
top.10.deg.L3.106<-top_10_degree("L3_106")
top.10.deg.L3.120<-top_10_degree("L3_120")
top.10.deg.L3.141<-top_10_degree("L3_141")

top.10.deg.FL.76<-top_10_degree("FL_76")
top.10.deg.FL.90<-top_10_degree("FL_90")
top.10.deg.FL.106<-top_10_degree("FL_106")
top.10.deg.FL.120<-top_10_degree("FL_120")
top.10.deg.FL.141<-top_10_degree("FL_141")

top.10.deg.G.0<-top_10_degree("G0")
top.10.deg.G.76<-top_10_degree("G76")
top.10.deg.G.90<-top_10_degree("G90")
top.10.deg.G.106<-top_10_degree("G106")
top.10.deg.G.120<-top_10_degree("G120")
top.10.deg.G.141<-top_10_degree("G141")


top10.deg.comp.age <- rbind(top.10.deg.BS.48, top.10.deg.BS.62, top.10.deg.BS.76, top.10.deg.BS.90,top.10.deg.BS.106,top.10.deg.BS.120,top.10.deg.BS.141,
                            top.10.deg.RS.48, top.10.deg.RS.62, top.10.deg.RS.76, top.10.deg.RS.90,top.10.deg.RS.106,top.10.deg.RS.120,top.10.deg.RS.141,
                            top.10.deg.R.48, top.10.deg.R.62, top.10.deg.R.76, top.10.deg.R.90,top.10.deg.R.106,top.10.deg.R.120,top.10.deg.R.141,
                            top.10.deg.S1.48, top.10.deg.S1.62, top.10.deg.S1.76, top.10.deg.S1.90,top.10.deg.S1.106,top.10.deg.S1.120,top.10.deg.S1.141,
                            top.10.deg.S2.48, top.10.deg.S2.62, top.10.deg.S2.76, top.10.deg.S2.90,top.10.deg.S2.106,top.10.deg.S2.120,top.10.deg.S2.141,
                            top.10.deg.S3.48, top.10.deg.S3.62, top.10.deg.S3.76, top.10.deg.S3.90,top.10.deg.S3.106,top.10.deg.S3.120,top.10.deg.S3.141,
                            top.10.deg.S4.62, top.10.deg.S4.76, top.10.deg.S4.90,top.10.deg.S4.106,top.10.deg.S4.120,top.10.deg.S4.141,
                            top.10.deg.S5.76, top.10.deg.S5.90,top.10.deg.S5.106,top.10.deg.S5.120,top.10.deg.S5.141,
                            top.10.deg.S6.76, top.10.deg.S6.90,top.10.deg.S6.106,top.10.deg.S6.120,top.10.deg.S6.141,
                            top.10.deg.S7.90,top.10.deg.S7.106,top.10.deg.S7.120,top.10.deg.S7.141,
                            top.10.deg.S8.90,top.10.deg.S8.106,top.10.deg.S8.120,top.10.deg.S8.141,
                            top.10.deg.S9.90,top.10.deg.S9.106,top.10.deg.S9.120,top.10.deg.S9.141,
                            top.10.deg.L1.48, top.10.deg.L1.62, top.10.deg.L1.76, top.10.deg.L1.90,top.10.deg.L1.106,top.10.deg.L1.120,top.10.deg.L1.141,
                            top.10.deg.L2.48, top.10.deg.L2.62, top.10.deg.L2.76, top.10.deg.L2.90,top.10.deg.L2.106,top.10.deg.L2.120,top.10.deg.L2.141,
                            top.10.deg.L3.62, top.10.deg.L3.76, top.10.deg.L3.90,top.10.deg.L3.106,top.10.deg.L3.120,top.10.deg.L3.141,
                            top.10.deg.FL.76, top.10.deg.FL.90,top.10.deg.FL.106,top.10.deg.FL.120,top.10.deg.FL.141,
                            top.10.deg.G.0, top.10.deg.G.76, top.10.deg.G.90,top.10.deg.G.106,top.10.deg.G.120,top.10.deg.G.141)

unique(top10.deg.comp.age$Node)

### Order by kingdom and phylum, class
otu.list <- read.csv("otu_list_of_all_communities.csv")
otu.list.hub <- subset(otu.list, OTU_id %in% unique(top10.deg.comp.age$Node))
length(otu.list.hub$OTU_id) #521

otu.list.hub$kingdom <- ifelse(grepl("^B",otu.list.hub$OTU_id) == T, "Bacteria","Fungi")
write.csv(otu.list.hub, "taxonomy of hub.csv")

otu.list.hub<-read.csv("taxonomy of hub.csv")

otu.list.hub$Phylum2 <- otu.list.hub$Phylum
otu.list.hub$Phylum2 <- as.character(otu.list.hub$Phylum2)
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
otu.list.hub$Phylum2[which(otu.list.hub$Kingdom == "Bacteria" &otu.list.hub$Phylum == "unidentified" & otu.list.hub$Class == "unidentified")] <- "Bacteria_unidientified"

otu.list.hub$Phylum2[which(otu.list.hub$Class == "Dothideomycetes")] <- "Dothideomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Eurotiomycetes")] <- "Eurotiomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Leotiomycetes")] <- "Leotiomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Pezizomycetes")] <- "Pezizomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Saccharomycetes")] <- "Saccharomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Sordariomycetes")] <- "Sordariomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Saccharomycetes")] <- "Saccharomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Saccharomycetes")] <- "Saccharomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Phylum == "Ascomycota" & otu.list.hub$Class == "unidentified")] <- "Asco_unidientified"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Agaricomycetes")] <- "Agaricomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Exobasidiomycetes")] <- "Exobasidiomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Malasseziomycetes")] <- "Malasseziomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Microbotryomycetes")] <- "Microbotryomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Tremellomycetes")] <- "Tremellomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Phylum == "Basidiomycota" & otu.list.hub$Class == "unidentified")] <- "Basidio_unidientified"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Ustilaginomycetes")] <- "Ustilaginomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Class == "Wallemiomycetes")] <- "Wallemiomycetes"
otu.list.hub$Phylum2[which(otu.list.hub$Kingdom == "Fungi" &otu.list.hub$Phylum == "unidentified" & otu.list.hub$Class == "unidentified")] <- "Fungi_unidientified"
unique(otu.list.hub$Phylum2)

otu.list.hub$Phylum2

write.csv(otu.list.gen.hub,"hub table gen.csv")

top10.deg.comp.age$Phylum2 <- 0
for (i in as.character(unique(top10.deg.comp.age$Node))){
  top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Node == i)] <-  unique(otu.list.hub$Phylum2[which(otu.list.hub$OTU_id == i)])
  }

### OTU level
hub.count<-read.csv("Hub count.csv")


hub.count$BS_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_48")] == T, 1, 0)
hub.count$BS_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_62")] == T, 1, 0)
hub.count$BS_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_76")] == T, 1, 0)
hub.count$BS_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_90")] == T, 1, 0)
hub.count$BS_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_106")] == T, 1, 0)
hub.count$BS_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_120")] == T, 1, 0)
hub.count$BS_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "BS_141")] == T, 1, 0)

hub.count$RS_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_48")] == T, 1, 0)
hub.count$RS_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_62")] == T, 1, 0)
hub.count$RS_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_76")] == T, 1, 0)
hub.count$RS_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_90")] == T, 1, 0)
hub.count$RS_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_106")] == T, 1, 0)
hub.count$RS_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_120")] == T, 1, 0)
hub.count$RS_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "RS_141")] == T, 1, 0)

hub.count$R_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_48")] == T, 1, 0)
hub.count$R_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_62")] == T, 1, 0)
hub.count$R_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_76")] == T, 1, 0)
hub.count$R_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_90")] == T, 1, 0)
hub.count$R_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_106")] == T, 1, 0)
hub.count$R_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_120")] == T, 1, 0)
hub.count$R_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "R_141")] == T, 1, 0)

hub.count$S1_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_48")] == T, 1, 0)
hub.count$S1_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_62")] == T, 1, 0)
hub.count$S1_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_76")] == T, 1, 0)
hub.count$S1_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_90")] == T, 1, 0)
hub.count$S1_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_106")] == T, 1, 0)
hub.count$S1_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_120")] == T, 1, 0)
hub.count$S1_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S1_141")] == T, 1, 0)

hub.count$S2_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_48")] == T, 1, 0)
hub.count$S2_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_62")] == T, 1, 0)
hub.count$S2_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_76")] == T, 1, 0)
hub.count$S2_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_90")] == T, 1, 0)
hub.count$S2_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_106")] == T, 1, 0)
hub.count$S2_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_120")] == T, 1, 0)
hub.count$S2_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S2_141")] == T, 1, 0)

hub.count$S3_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_48")] == T, 1, 0)
hub.count$S3_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_62")] == T, 1, 0)
hub.count$S3_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_76")] == T, 1, 0)
hub.count$S3_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_90")] == T, 1, 0)
hub.count$S3_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_106")] == T, 1, 0)
hub.count$S3_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_120")] == T, 1, 0)
hub.count$S3_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S3_141")] == T, 1, 0)

hub.count$S4_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_62")] == T, 1, 0)
hub.count$S4_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_76")] == T, 1, 0)
hub.count$S4_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_90")] == T, 1, 0)
hub.count$S4_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_106")] == T, 1, 0)
hub.count$S4_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_120")] == T, 1, 0)
hub.count$S4_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S4_141")] == T, 1, 0)

hub.count$S5_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S5_76")] == T, 1, 0)
hub.count$S5_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S5_90")] == T, 1, 0)
hub.count$S5_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S5_106")] == T, 1, 0)
hub.count$S5_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S5_120")] == T, 1, 0)
hub.count$S5_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S5_141")] == T, 1, 0)

hub.count$S6_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S6_76")] == T, 1, 0)
hub.count$S6_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S6_90")] == T, 1, 0)
hub.count$S6_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S6_106")] == T, 1, 0)
hub.count$S6_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S6_120")] == T, 1, 0)
hub.count$S6_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S6_141")] == T, 1, 0)

hub.count$S7_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S7_90")] == T, 1, 0)
hub.count$S7_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S7_106")] == T, 1, 0)
hub.count$S7_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S7_120")] == T, 1, 0)
hub.count$S7_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S7_141")] == T, 1, 0)

hub.count$S8_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S8_90")] == T, 1, 0)
hub.count$S8_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S8_106")] == T, 1, 0)
hub.count$S8_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S8_120")] == T, 1, 0)
hub.count$S8_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S8_141")] == T, 1, 0)

hub.count$S9_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S9_90")] == T, 1, 0)
hub.count$S9_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S9_106")] == T, 1, 0)
hub.count$S9_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S9_120")] == T, 1, 0)
hub.count$S9_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "S9_141")] == T, 1, 0)

hub.count$L1_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_48")] == T, 1, 0)
hub.count$L1_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_62")] == T, 1, 0)
hub.count$L1_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_76")] == T, 1, 0)
hub.count$L1_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_90")] == T, 1, 0)
hub.count$L1_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_106")] == T, 1, 0)
hub.count$L1_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_120")] == T, 1, 0)
hub.count$L1_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L1_141")] == T, 1, 0)

hub.count$L2_48 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_48")] == T, 1, 0)
hub.count$L2_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_62")] == T, 1, 0)
hub.count$L2_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_76")] == T, 1, 0)
hub.count$L2_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_90")] == T, 1, 0)
hub.count$L2_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_106")] == T, 1, 0)
hub.count$L2_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_120")] == T, 1, 0)
hub.count$L2_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L2_141")] == T, 1, 0)

hub.count$L3_62 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_62")] == T, 1, 0)
hub.count$L3_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_76")] == T, 1, 0)
hub.count$L3_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_90")] == T, 1, 0)
hub.count$L3_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_106")] == T, 1, 0)
hub.count$L3_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_120")] == T, 1, 0)
hub.count$L3_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "L3_141")] == T, 1, 0)

hub.count$FL_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "FL_76")] == T, 1, 0)
hub.count$FL_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "FL_90")] == T, 1, 0)
hub.count$FL_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "FL_106")] == T, 1, 0)
hub.count$FL_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "FL_120")] == T, 1, 0)
hub.count$FL_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "FL_141")] == T, 1, 0)

hub.count$G_0 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G0")] == T, 1, 0)
hub.count$G_76 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G76")] == T, 1, 0)
hub.count$G_90 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G90")] == T, 1, 0)
hub.count$G_106 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G106")] == T, 1, 0)
hub.count$G_120 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G120")] == T, 1, 0)
hub.count$G_141 <- ifelse(hub.count$Hub %in% top10.deg.comp.age$Node[which(top10.deg.comp.age$Network == "G141")] == T, 1, 0)

head(hub.count)

write.csv(hub.count, "hub count final.csv")



### Phyla and class level
hub.count.phyla <- read.csv('Hub count phyla.csv')
hub.count.phyla$BS_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_48")] == T, 1, 0)
hub.count.phyla$BS_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_62")] == T, 1, 0)
hub.count.phyla$BS_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_76")] == T, 1, 0)
hub.count.phyla$BS_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_90")] == T, 1, 0)
hub.count.phyla$BS_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_106")] == T, 1, 0)
hub.count.phyla$BS_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_120")] == T, 1, 0)
hub.count.phyla$BS_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "BS_141")] == T, 1, 0)

hub.count.phyla$RS_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_48")] == T, 1, 0)
hub.count.phyla$RS_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_62")] == T, 1, 0)
hub.count.phyla$RS_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_76")] == T, 1, 0)
hub.count.phyla$RS_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_90")] == T, 1, 0)
hub.count.phyla$RS_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_106")] == T, 1, 0)
hub.count.phyla$RS_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_120")] == T, 1, 0)
hub.count.phyla$RS_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "RS_141")] == T, 1, 0)

hub.count.phyla$R_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_48")] == T, 1, 0)
hub.count.phyla$R_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_62")] == T, 1, 0)
hub.count.phyla$R_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_76")] == T, 1, 0)
hub.count.phyla$R_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_90")] == T, 1, 0)
hub.count.phyla$R_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_106")] == T, 1, 0)
hub.count.phyla$R_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_120")] == T, 1, 0)
hub.count.phyla$R_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "R_141")] == T, 1, 0)

hub.count.phyla$S1_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_48")] == T, 1, 0)
hub.count.phyla$S1_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_62")] == T, 1, 0)
hub.count.phyla$S1_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_76")] == T, 1, 0)
hub.count.phyla$S1_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_90")] == T, 1, 0)
hub.count.phyla$S1_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_106")] == T, 1, 0)
hub.count.phyla$S1_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_120")] == T, 1, 0)
hub.count.phyla$S1_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S1_141")] == T, 1, 0)

hub.count.phyla$S2_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_48")] == T, 1, 0)
hub.count.phyla$S2_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_62")] == T, 1, 0)
hub.count.phyla$S2_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_76")] == T, 1, 0)
hub.count.phyla$S2_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_90")] == T, 1, 0)
hub.count.phyla$S2_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_106")] == T, 1, 0)
hub.count.phyla$S2_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_120")] == T, 1, 0)
hub.count.phyla$S2_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S2_141")] == T, 1, 0)

hub.count.phyla$S3_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_48")] == T, 1, 0)
hub.count.phyla$S3_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_62")] == T, 1, 0)
hub.count.phyla$S3_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_76")] == T, 1, 0)
hub.count.phyla$S3_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_90")] == T, 1, 0)
hub.count.phyla$S3_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_106")] == T, 1, 0)
hub.count.phyla$S3_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_120")] == T, 1, 0)
hub.count.phyla$S3_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S3_141")] == T, 1, 0)

hub.count.phyla$S4_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_62")] == T, 1, 0)
hub.count.phyla$S4_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_76")] == T, 1, 0)
hub.count.phyla$S4_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_90")] == T, 1, 0)
hub.count.phyla$S4_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_106")] == T, 1, 0)
hub.count.phyla$S4_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_120")] == T, 1, 0)
hub.count.phyla$S4_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S4_141")] == T, 1, 0)

hub.count.phyla$S5_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S5_76")] == T, 1, 0)
hub.count.phyla$S5_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S5_90")] == T, 1, 0)
hub.count.phyla$S5_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S5_106")] == T, 1, 0)
hub.count.phyla$S5_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S5_120")] == T, 1, 0)
hub.count.phyla$S5_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S5_141")] == T, 1, 0)

hub.count.phyla$S6_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S6_76")] == T, 1, 0)
hub.count.phyla$S6_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S6_90")] == T, 1, 0)
hub.count.phyla$S6_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S6_106")] == T, 1, 0)
hub.count.phyla$S6_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S6_120")] == T, 1, 0)
hub.count.phyla$S6_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S6_141")] == T, 1, 0)

hub.count.phyla$S7_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S7_90")] == T, 1, 0)
hub.count.phyla$S7_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S7_106")] == T, 1, 0)
hub.count.phyla$S7_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S7_120")] == T, 1, 0)
hub.count.phyla$S7_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S7_141")] == T, 1, 0)

hub.count.phyla$S8_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S8_90")] == T, 1, 0)
hub.count.phyla$S8_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S8_106")] == T, 1, 0)
hub.count.phyla$S8_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S8_120")] == T, 1, 0)
hub.count.phyla$S8_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S8_141")] == T, 1, 0)

hub.count.phyla$S9_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S9_90")] == T, 1, 0)
hub.count.phyla$S9_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S9_106")] == T, 1, 0)
hub.count.phyla$S9_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S9_120")] == T, 1, 0)
hub.count.phyla$S9_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "S9_141")] == T, 1, 0)

hub.count.phyla$L1_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_48")] == T, 1, 0)
hub.count.phyla$L1_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_62")] == T, 1, 0)
hub.count.phyla$L1_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_76")] == T, 1, 0)
hub.count.phyla$L1_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_90")] == T, 1, 0)
hub.count.phyla$L1_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_106")] == T, 1, 0)
hub.count.phyla$L1_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_120")] == T, 1, 0)
hub.count.phyla$L1_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L1_141")] == T, 1, 0)

hub.count.phyla$L2_48 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_48")] == T, 1, 0)
hub.count.phyla$L2_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_62")] == T, 1, 0)
hub.count.phyla$L2_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_76")] == T, 1, 0)
hub.count.phyla$L2_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_90")] == T, 1, 0)
hub.count.phyla$L2_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_106")] == T, 1, 0)
hub.count.phyla$L2_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_120")] == T, 1, 0)
hub.count.phyla$L2_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L2_141")] == T, 1, 0)

hub.count.phyla$L3_62 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_62")] == T, 1, 0)
hub.count.phyla$L3_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_76")] == T, 1, 0)
hub.count.phyla$L3_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_90")] == T, 1, 0)
hub.count.phyla$L3_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_106")] == T, 1, 0)
hub.count.phyla$L3_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_120")] == T, 1, 0)
hub.count.phyla$L3_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "L3_141")] == T, 1, 0)

hub.count.phyla$FL_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "FL_76")] == T, 1, 0)
hub.count.phyla$FL_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "FL_90")] == T, 1, 0)
hub.count.phyla$FL_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "FL_106")] == T, 1, 0)
hub.count.phyla$FL_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "FL_120")] == T, 1, 0)
hub.count.phyla$FL_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "FL_141")] == T, 1, 0)

hub.count.phyla$G_0 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G0")] == T, 1, 0)
hub.count.phyla$G_76 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G76")] == T, 1, 0)
hub.count.phyla$G_90 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G90")] == T, 1, 0)
hub.count.phyla$G_106 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G106")] == T, 1, 0)
hub.count.phyla$G_120 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G120")] == T, 1, 0)
hub.count.phyla$G_141 <- ifelse(hub.count.phyla$Hub %in% top10.deg.comp.age$Phylum2[which(top10.deg.comp.age$Network == "G141")] == T, 1, 0)

write.csv(hub.count.phyla,"hub.count.phyla.csv")



### Genus level
otu.list.2<- tidyr::unite(otu.list, "Class2", Phylum, Class, sep = "_")
otu.list.3<- tidyr::unite(otu.list.2, "Order", Class2, Order, sep = "_")
otu.list.4<- tidyr::unite(otu.list.3, "Family", Order, Family, sep = "_")
otu.list.5<- tidyr::unite(otu.list.4, "Genus", Family, Genus, sep = "_")


otu.list.gen.hub <- subset(otu.list.5, OTU_id %in% unique(top10.deg.comp.age$Node))
otu.list.gen.hub <- read.csv('hub table gen.csv')
top10.deg.comp.age$Genus2 <- 0
for (i in as.character(unique(top10.deg.comp.age$Node))){
  top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Node == i)] <-  unique(as.character(otu.list.gen.hub$Genus[which(otu.list.gen.hub$OTU_id == i)]))
}

hub.count.genus <- read.csv('Hub count genus.csv')
hub.count.genus$BS_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_48")] == T, 1, 0)
hub.count.genus$BS_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_62")] == T, 1, 0)
hub.count.genus$BS_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_76")] == T, 1, 0)
hub.count.genus$BS_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_90")] == T, 1, 0)
hub.count.genus$BS_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_106")] == T, 1, 0)
hub.count.genus$BS_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_120")] == T, 1, 0)
hub.count.genus$BS_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "BS_141")] == T, 1, 0)

hub.count.genus$RS_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_48")] == T, 1, 0)
hub.count.genus$RS_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_62")] == T, 1, 0)
hub.count.genus$RS_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_76")] == T, 1, 0)
hub.count.genus$RS_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_90")] == T, 1, 0)
hub.count.genus$RS_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_106")] == T, 1, 0)
hub.count.genus$RS_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_120")] == T, 1, 0)
hub.count.genus$RS_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "RS_141")] == T, 1, 0)

hub.count.genus$R_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_48")] == T, 1, 0)
hub.count.genus$R_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_62")] == T, 1, 0)
hub.count.genus$R_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_76")] == T, 1, 0)
hub.count.genus$R_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_90")] == T, 1, 0)
hub.count.genus$R_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_106")] == T, 1, 0)
hub.count.genus$R_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_120")] == T, 1, 0)
hub.count.genus$R_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "R_141")] == T, 1, 0)

hub.count.genus$S1_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_48")] == T, 1, 0)
hub.count.genus$S1_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_62")] == T, 1, 0)
hub.count.genus$S1_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_76")] == T, 1, 0)
hub.count.genus$S1_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_90")] == T, 1, 0)
hub.count.genus$S1_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_106")] == T, 1, 0)
hub.count.genus$S1_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_120")] == T, 1, 0)
hub.count.genus$S1_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S1_141")] == T, 1, 0)

hub.count.genus$S2_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_48")] == T, 1, 0)
hub.count.genus$S2_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_62")] == T, 1, 0)
hub.count.genus$S2_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_76")] == T, 1, 0)
hub.count.genus$S2_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_90")] == T, 1, 0)
hub.count.genus$S2_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_106")] == T, 1, 0)
hub.count.genus$S2_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_120")] == T, 1, 0)
hub.count.genus$S2_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S2_141")] == T, 1, 0)

hub.count.genus$S3_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_48")] == T, 1, 0)
hub.count.genus$S3_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_62")] == T, 1, 0)
hub.count.genus$S3_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_76")] == T, 1, 0)
hub.count.genus$S3_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_90")] == T, 1, 0)
hub.count.genus$S3_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_106")] == T, 1, 0)
hub.count.genus$S3_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_120")] == T, 1, 0)
hub.count.genus$S3_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S3_141")] == T, 1, 0)

hub.count.genus$S4_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_62")] == T, 1, 0)
hub.count.genus$S4_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_76")] == T, 1, 0)
hub.count.genus$S4_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_90")] == T, 1, 0)
hub.count.genus$S4_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_106")] == T, 1, 0)
hub.count.genus$S4_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_120")] == T, 1, 0)
hub.count.genus$S4_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S4_141")] == T, 1, 0)

hub.count.genus$S5_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S5_76")] == T, 1, 0)
hub.count.genus$S5_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S5_90")] == T, 1, 0)
hub.count.genus$S5_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S5_106")] == T, 1, 0)
hub.count.genus$S5_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S5_120")] == T, 1, 0)
hub.count.genus$S5_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S5_141")] == T, 1, 0)

hub.count.genus$S6_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S6_76")] == T, 1, 0)
hub.count.genus$S6_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S6_90")] == T, 1, 0)
hub.count.genus$S6_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S6_106")] == T, 1, 0)
hub.count.genus$S6_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S6_120")] == T, 1, 0)
hub.count.genus$S6_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S6_141")] == T, 1, 0)

hub.count.genus$S7_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S7_90")] == T, 1, 0)
hub.count.genus$S7_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S7_106")] == T, 1, 0)
hub.count.genus$S7_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S7_120")] == T, 1, 0)
hub.count.genus$S7_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S7_141")] == T, 1, 0)

hub.count.genus$S8_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S8_90")] == T, 1, 0)
hub.count.genus$S8_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S8_106")] == T, 1, 0)
hub.count.genus$S8_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S8_120")] == T, 1, 0)
hub.count.genus$S8_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S8_141")] == T, 1, 0)

hub.count.genus$S9_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S9_90")] == T, 1, 0)
hub.count.genus$S9_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S9_106")] == T, 1, 0)
hub.count.genus$S9_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S9_120")] == T, 1, 0)
hub.count.genus$S9_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "S9_141")] == T, 1, 0)

hub.count.genus$L1_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_48")] == T, 1, 0)
hub.count.genus$L1_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_62")] == T, 1, 0)
hub.count.genus$L1_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_76")] == T, 1, 0)
hub.count.genus$L1_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_90")] == T, 1, 0)
hub.count.genus$L1_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_106")] == T, 1, 0)
hub.count.genus$L1_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_120")] == T, 1, 0)
hub.count.genus$L1_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L1_141")] == T, 1, 0)

hub.count.genus$L2_48 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_48")] == T, 1, 0)
hub.count.genus$L2_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_62")] == T, 1, 0)
hub.count.genus$L2_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_76")] == T, 1, 0)
hub.count.genus$L2_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_90")] == T, 1, 0)
hub.count.genus$L2_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_106")] == T, 1, 0)
hub.count.genus$L2_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_120")] == T, 1, 0)
hub.count.genus$L2_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L2_141")] == T, 1, 0)

hub.count.genus$L3_62 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_62")] == T, 1, 0)
hub.count.genus$L3_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_76")] == T, 1, 0)
hub.count.genus$L3_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_90")] == T, 1, 0)
hub.count.genus$L3_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_106")] == T, 1, 0)
hub.count.genus$L3_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_120")] == T, 1, 0)
hub.count.genus$L3_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "L3_141")] == T, 1, 0)

hub.count.genus$FL_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "FL_76")] == T, 1, 0)
hub.count.genus$FL_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "FL_90")] == T, 1, 0)
hub.count.genus$FL_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "FL_106")] == T, 1, 0)
hub.count.genus$FL_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "FL_120")] == T, 1, 0)
hub.count.genus$FL_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "FL_141")] == T, 1, 0)

hub.count.genus$G_0 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G0")] == T, 1, 0)
hub.count.genus$G_76 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G76")] == T, 1, 0)
hub.count.genus$G_90 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G90")] == T, 1, 0)
hub.count.genus$G_106 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G106")] == T, 1, 0)
hub.count.genus$G_120 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G120")] == T, 1, 0)
hub.count.genus$G_141 <- ifelse(hub.count.genus$Hub %in% top10.deg.comp.age$Genus2[which(top10.deg.comp.age$Network == "G141")] == T, 1, 0)

write.csv(hub.count.genus,"hub.count.genus.csv")
