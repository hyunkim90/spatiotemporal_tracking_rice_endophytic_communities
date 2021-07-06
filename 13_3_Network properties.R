### Network properties
## read cor and p

G0_cor_df_padj <- read.table("1_0.3_edge_G0.tsv", sep='\t', header =T)

nodeattrib_G0_combine <- data.frame(node=union(G0_cor_df_padj$Source,G0_cor_df_padj$Target))
nodeattrib_G0_combine$kingdom <- 0

for (i in as.character(nodeattrib_G0_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G0_combine[nodeattrib_G0_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G0_combine[nodeattrib_G0_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G0_combine) <- as.character(nodeattrib_G0_combine$node)
nodeattrib_G0_combine$kingdom



all_G0_net <- graph_from_data_frame(G0_cor_df_padj,direct=F, vertices=nodeattrib_G0_combine)

## Number of nodes
length(V(all_G0_net)) # 80

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G0_net)))) #54
length(grep("^F",names(V(all_G0_net)))) #26


## Connections 
bb_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G0) #100

ff_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G0) #10

fb_occur_G0 <- droplevels(G0_cor_df_padj[with(G0_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G0) #41



## Network properties
meta_degree <- sort(igraph::degree(all_G0_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G0_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G0_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G0_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G0_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G0_net, V(all_G0_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G0_net
G0_all_deg <- igraph::degree(net,mode="all")
G0_all_betweenness <- betweenness(net, normalized = TRUE)
G0_all_closeness <- closeness(net, normalized = TRUE)
G0_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G0_all_transitivity)<- V(net)$name
G0_all_transitivity[is.na(G0_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G0_all_deg.1percent <- G0_all_deg[G0_all_deg >= quantile(G0_all_deg,prob=1-n/100)]
length(G0_all_deg.1percent) #8

G0_all_betweenness.1percent <- G0_all_betweenness[G0_all_betweenness >= quantile(G0_all_betweenness,prob=1-n/100)]
length(G0_all_betweenness.1percent) #7

G0_all_closeness.1percent <- G0_all_closeness[G0_all_closeness >= quantile(G0_all_closeness,prob=1-n/100)]
length(G0_all_closeness.1percent) #7

intersect(names(G0_all_deg.1percent), names(G0_all_betweenness.1percent))
# [1] "B9_Methylobacterium"                                   
# [2] "B37_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
intersect(names(G0_all_deg.1percent), names(G0_all_closeness.1percent))
intersect(names(G0_all_betweenness.1percent), names(G0_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G0.degree<-data.frame(G0_all_deg)
head(df.G0.degree)
df.G0.degree$Group <- "G0"
names(df.G0.degree)[1] <- c("Degree")

df.G0.closeness<-data.frame(G0_all_closeness)
head(df.G0.closeness)
df.G0.closeness$Group <- "G0"
names(df.G0.closeness)[1] <- c("Closeness")

df.G0.betweenness<-data.frame(G0_all_betweenness)
head(df.G0.betweenness)
df.G0.betweenness$Group <- "G0"
names(df.G0.betweenness)[1] <- c("Betweenness")

df.G0.degree$kingdom <- ifelse(grepl("^B",rownames(df.G0.degree)),'Bacteria', 'Fungi')
df.G0.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G0.closeness)),'Bacteria', 'Fungi')
df.G0.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G0.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G0.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G0.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G0.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G0.degree)
head(df.G0.betweenness)
df.G0.degree$OTU_id <- rownames(df.G0.degree)
df.G0.betweenness$OTU_id <- rownames(df.G0.betweenness)
df.G0.closeness$OTU_id <- rownames(df.G0.closeness)
df.G0.net.properties <- merge(df.G0.degree, df.G0.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G0.net.properties <- merge(df.G0.net.properties, df.G0.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G0.net.properties$kingdom <- factor(df.G0.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-3
ggplot(df.G0.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G0.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G0.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G0.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G0.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G0.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

### G76
### Network properties
## read cor and p

G76_cor_df_padj <- read.table("1_0.3_edge_G76.tsv", sep='\t', header =T)

nodeattrib_G76_combine <- data.frame(node=union(G76_cor_df_padj$Source,G76_cor_df_padj$Target))
nodeattrib_G76_combine$kingdom <- 0

for (i in as.character(nodeattrib_G76_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G76_combine[nodeattrib_G76_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G76_combine[nodeattrib_G76_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G76_combine) <- as.character(nodeattrib_G76_combine$node)
nodeattrib_G76_combine$kingdom



all_G76_net <- graph_from_data_frame(G76_cor_df_padj,direct=F, vertices=nodeattrib_G76_combine)

## Number of nodes
length(V(all_G76_net)) # 363

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G76_net)))) #163
length(grep("^F",names(V(all_G76_net)))) #200


## Connections 
bb_occur_G76 <- droplevels(G76_cor_df_padj[with(G76_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G76) #484

ff_occur_G76 <- droplevels(G76_cor_df_padj[with(G76_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G76) #660

fb_occur_G76 <- droplevels(G76_cor_df_padj[with(G76_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G76) #889



## Network properties
meta_degree <- sort(igraph::degree(all_G76_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G76_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G76_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G76_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G76_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G76_net, V(all_G76_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G76_net
G76_all_deg <- igraph::degree(net,mode="all")
G76_all_betweenness <- betweenness(net, normalized = TRUE)
G76_all_closeness <- closeness(net, normalized = TRUE)
G76_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G76_all_transitivity)<- V(net)$name
G76_all_transitivity[is.na(G76_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G76_all_deg.1percent <- G76_all_deg[G76_all_deg >= quantile(G76_all_deg,prob=1-n/100)]
length(G76_all_deg.1percent) #8

G76_all_betweenness.1percent <- G76_all_betweenness[G76_all_betweenness >= quantile(G76_all_betweenness,prob=1-n/100)]
length(G76_all_betweenness.1percent) #7

G76_all_closeness.1percent <- G76_all_closeness[G76_all_closeness >= quantile(G76_all_closeness,prob=1-n/100)]
length(G76_all_closeness.1percent) #7

intersect(names(G76_all_deg.1percent), names(G76_all_betweenness.1percent))
#[1] "F131_unidentified"        "F590_o_Saccharomycetales"

intersect(names(G76_all_deg.1percent), names(G76_all_closeness.1percent))
intersect(names(G76_all_betweenness.1percent), names(G76_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G76.degree<-data.frame(G76_all_deg)
head(df.G76.degree)
df.G76.degree$Group <- "G76"
names(df.G76.degree)[1] <- c("Degree")

df.G76.closeness<-data.frame(G76_all_closeness)
head(df.G76.closeness)
df.G76.closeness$Group <- "G76"
names(df.G76.closeness)[1] <- c("Closeness")

df.G76.betweenness<-data.frame(G76_all_betweenness)
head(df.G76.betweenness)
df.G76.betweenness$Group <- "G76"
names(df.G76.betweenness)[1] <- c("Betweenness")

df.G76.degree$kingdom <- ifelse(grepl("^B",rownames(df.G76.degree)),'Bacteria', 'Fungi')
df.G76.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G76.closeness)),'Bacteria', 'Fungi')
df.G76.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G76.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G76.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G76.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G76.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G76.degree)
head(df.G76.betweenness)
df.G76.degree$OTU_id <- rownames(df.G76.degree)
df.G76.betweenness$OTU_id <- rownames(df.G76.betweenness)
df.G76.closeness$OTU_id <- rownames(df.G76.closeness)
df.G76.net.properties <- merge(df.G76.degree, df.G76.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G76.net.properties <- merge(df.G76.net.properties, df.G76.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G76.net.properties$kingdom <- factor(df.G76.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-1
ggplot(df.G76.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G76.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G76.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G76.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G76.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G76.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)



### G90
### Network properties
## read cor and p

G90_cor_df_padj <- read.table("1_0.3_edge_G90.tsv", sep='\t', header =T)

nodeattrib_G90_combine <- data.frame(node=union(G90_cor_df_padj$Source,G90_cor_df_padj$Target))
nodeattrib_G90_combine$kingdom <- 0

for (i in as.character(nodeattrib_G90_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G90_combine[nodeattrib_G90_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G90_combine[nodeattrib_G90_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G90_combine) <- as.character(nodeattrib_G90_combine$node)
nodeattrib_G90_combine$kingdom



all_G90_net <- graph_from_data_frame(G90_cor_df_padj,direct=F, vertices=nodeattrib_G90_combine)

## Number of nodes
length(V(all_G90_net)) # 192

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G90_net)))) #62
length(grep("^F",names(V(all_G90_net)))) #130


## Connections 
bb_occur_G90 <- droplevels(G90_cor_df_padj[with(G90_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G90) #73

ff_occur_G90 <- droplevels(G90_cor_df_padj[with(G90_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G90) #154

fb_occur_G90 <- droplevels(G90_cor_df_padj[with(G90_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G90) #115



## Network properties
meta_degree <- sort(igraph::degree(all_G90_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G90_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G90_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G90_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G90_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G90_net, V(all_G90_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G90_net
G90_all_deg <- igraph::degree(net,mode="all")
G90_all_betweenness <- betweenness(net, normalized = TRUE)
G90_all_closeness <- closeness(net, normalized = TRUE)
G90_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G90_all_transitivity)<- V(net)$name
G90_all_transitivity[is.na(G90_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G90_all_deg.1percent <- G90_all_deg[G90_all_deg >= quantile(G90_all_deg,prob=1-n/100)]
length(G90_all_deg.1percent) #8

G90_all_betweenness.1percent <- G90_all_betweenness[G90_all_betweenness >= quantile(G90_all_betweenness,prob=1-n/100)]
length(G90_all_betweenness.1percent) #7

G90_all_closeness.1percent <- G90_all_closeness[G90_all_closeness >= quantile(G90_all_closeness,prob=1-n/100)]
length(G90_all_closeness.1percent) #7

intersect(names(G90_all_deg.1percent), names(G90_all_betweenness.1percent))
intersect(names(G90_all_deg.1percent), names(G90_all_closeness.1percent))
intersect(names(G90_all_betweenness.1percent), names(G90_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G90.degree<-data.frame(G90_all_deg)
head(df.G90.degree)
df.G90.degree$Group <- "G90"
names(df.G90.degree)[1] <- c("Degree")

df.G90.closeness<-data.frame(G90_all_closeness)
head(df.G90.closeness)
df.G90.closeness$Group <- "G90"
names(df.G90.closeness)[1] <- c("Closeness")

df.G90.betweenness<-data.frame(G90_all_betweenness)
head(df.G90.betweenness)
df.G90.betweenness$Group <- "G90"
names(df.G90.betweenness)[1] <- c("Betweenness")

df.G90.degree$kingdom <- ifelse(grepl("^B",rownames(df.G90.degree)),'Bacteria', 'Fungi')
df.G90.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G90.closeness)),'Bacteria', 'Fungi')
df.G90.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G90.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G90.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G90.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G90.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G90.degree)
head(df.G90.betweenness)
df.G90.degree$OTU_id <- rownames(df.G90.degree)
df.G90.betweenness$OTU_id <- rownames(df.G90.betweenness)
df.G90.closeness$OTU_id <- rownames(df.G90.closeness)
df.G90.net.properties <- merge(df.G90.degree, df.G90.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G90.net.properties <- merge(df.G90.net.properties, df.G90.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G90.net.properties$kingdom <- factor(df.G90.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-3
ggplot(df.G90.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G90.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G90.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G90.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G90.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G90.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)


### G106
### Network properties
## read cor and p

G106_cor_df_padj <- read.table("1_0.3_edge_G106.tsv", sep='\t', header =T)

nodeattrib_G106_combine <- data.frame(node=union(G106_cor_df_padj$Source,G106_cor_df_padj$Target))
nodeattrib_G106_combine$kingdom <- 0

for (i in as.character(nodeattrib_G106_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G106_combine[nodeattrib_G106_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G106_combine[nodeattrib_G106_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G106_combine) <- as.character(nodeattrib_G106_combine$node)
nodeattrib_G106_combine$kingdom



all_G106_net <- graph_from_data_frame(G106_cor_df_padj,direct=F, vertices=nodeattrib_G106_combine)

## Number of nodes
length(V(all_G106_net)) # 192

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G106_net)))) #29
length(grep("^F",names(V(all_G106_net)))) #77


## Connections 
bb_occur_G106 <- droplevels(G106_cor_df_padj[with(G106_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G106) #29

ff_occur_G106 <- droplevels(G106_cor_df_padj[with(G106_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G106) #103

fb_occur_G106 <- droplevels(G106_cor_df_padj[with(G106_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G106) #51



## Network properties
meta_degree <- sort(igraph::degree(all_G106_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G106_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G106_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G106_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G106_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G106_net, V(all_G106_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G106_net
G106_all_deg <- igraph::degree(net,mode="all")
G106_all_betweenness <- betweenness(net, normalized = TRUE)
G106_all_closeness <- closeness(net, normalized = TRUE)
G106_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G106_all_transitivity)<- V(net)$name
G106_all_transitivity[is.na(G106_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G106_all_deg.1percent <- G106_all_deg[G106_all_deg >= quantile(G106_all_deg,prob=1-n/100)]
length(G106_all_deg.1percent) #8

G106_all_betweenness.1percent <- G106_all_betweenness[G106_all_betweenness >= quantile(G106_all_betweenness,prob=1-n/100)]
length(G106_all_betweenness.1percent) #7

G106_all_closeness.1percent <- G106_all_closeness[G106_all_closeness >= quantile(G106_all_closeness,prob=1-n/100)]
length(G106_all_closeness.1percent) #7

intersect(names(G106_all_deg.1percent), names(G106_all_betweenness.1percent))
intersect(names(G106_all_deg.1percent), names(G106_all_closeness.1percent))
intersect(names(G106_all_betweenness.1percent), names(G106_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G106.degree<-data.frame(G106_all_deg)
head(df.G106.degree)
df.G106.degree$Group <- "G106"
names(df.G106.degree)[1] <- c("Degree")

df.G106.closeness<-data.frame(G106_all_closeness)
head(df.G106.closeness)
df.G106.closeness$Group <- "G106"
names(df.G106.closeness)[1] <- c("Closeness")

df.G106.betweenness<-data.frame(G106_all_betweenness)
head(df.G106.betweenness)
df.G106.betweenness$Group <- "G106"
names(df.G106.betweenness)[1] <- c("Betweenness")

df.G106.degree$kingdom <- ifelse(grepl("^B",rownames(df.G106.degree)),'Bacteria', 'Fungi')
df.G106.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G106.closeness)),'Bacteria', 'Fungi')
df.G106.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G106.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G106.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G106.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G106.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G106.degree)
head(df.G106.betweenness)
df.G106.degree$OTU_id <- rownames(df.G106.degree)
df.G106.betweenness$OTU_id <- rownames(df.G106.betweenness)
df.G106.closeness$OTU_id <- rownames(df.G106.closeness)
df.G106.net.properties <- merge(df.G106.degree, df.G106.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G106.net.properties <- merge(df.G106.net.properties, df.G106.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G106.net.properties$kingdom <- factor(df.G106.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-1
ggplot(df.G106.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G106.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G106.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G106.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G106.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G106.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)



### G120
### Network properties
## read cor and p

G120_cor_df_padj <- read.table("1_0.3_edge_G120.tsv", sep='\t', header =T)

nodeattrib_G120_combine <- data.frame(node=union(G120_cor_df_padj$Source,G120_cor_df_padj$Target))
nodeattrib_G120_combine$kingdom <- 0

for (i in as.character(nodeattrib_G120_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G120_combine[nodeattrib_G120_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G120_combine[nodeattrib_G120_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G120_combine) <- as.character(nodeattrib_G120_combine$node)
nodeattrib_G120_combine$kingdom



all_G120_net <- graph_from_data_frame(G120_cor_df_padj,direct=F, vertices=nodeattrib_G120_combine)

## Number of nodes
length(V(all_G120_net)) # 72

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G120_net)))) #33
length(grep("^F",names(V(all_G120_net)))) #39


## Connections 
bb_occur_G120 <- droplevels(G120_cor_df_padj[with(G120_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G120) #56

ff_occur_G120 <- droplevels(G120_cor_df_padj[with(G120_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G120) #25

fb_occur_G120 <- droplevels(G120_cor_df_padj[with(G120_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G120) #47



## Network properties
meta_degree <- sort(igraph::degree(all_G120_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G120_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G120_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G120_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G120_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G120_net, V(all_G120_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G120_net
G120_all_deg <- igraph::degree(net,mode="all")
G120_all_betweenness <- betweenness(net, normalized = TRUE)
G120_all_closeness <- closeness(net, normalized = TRUE)
G120_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G120_all_transitivity)<- V(net)$name
G120_all_transitivity[is.na(G120_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G120_all_deg.1percent <- G120_all_deg[G120_all_deg >= quantile(G120_all_deg,prob=1-n/100)]
length(G120_all_deg.1percent) #8

G120_all_betweenness.1percent <- G120_all_betweenness[G120_all_betweenness >= quantile(G120_all_betweenness,prob=1-n/100)]
length(G120_all_betweenness.1percent) #7

G120_all_closeness.1percent <- G120_all_closeness[G120_all_closeness >= quantile(G120_all_closeness,prob=1-n/100)]
length(G120_all_closeness.1percent) #7

intersect(names(G120_all_deg.1percent), names(G120_all_betweenness.1percent))
intersect(names(G120_all_deg.1percent), names(G120_all_closeness.1percent))
intersect(names(G120_all_betweenness.1percent), names(G120_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G120.degree<-data.frame(G120_all_deg)
head(df.G120.degree)
df.G120.degree$Group <- "G120"
names(df.G120.degree)[1] <- c("Degree")

df.G120.closeness<-data.frame(G120_all_closeness)
head(df.G120.closeness)
df.G120.closeness$Group <- "G120"
names(df.G120.closeness)[1] <- c("Closeness")

df.G120.betweenness<-data.frame(G120_all_betweenness)
head(df.G120.betweenness)
df.G120.betweenness$Group <- "G120"
names(df.G120.betweenness)[1] <- c("Betweenness")

df.G120.degree$kingdom <- ifelse(grepl("^B",rownames(df.G120.degree)),'Bacteria', 'Fungi')
df.G120.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G120.closeness)),'Bacteria', 'Fungi')
df.G120.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G120.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G120.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G120.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G120.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G120.degree)
head(df.G120.betweenness)
df.G120.degree$OTU_id <- rownames(df.G120.degree)
df.G120.betweenness$OTU_id <- rownames(df.G120.betweenness)
df.G120.closeness$OTU_id <- rownames(df.G120.closeness)
df.G120.net.properties <- merge(df.G120.degree, df.G120.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G120.net.properties <- merge(df.G120.net.properties, df.G120.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G120.net.properties$kingdom <- factor(df.G120.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-1
ggplot(df.G120.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G120.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G120.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G120.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G120.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G120.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)


### G141
### Network properties
## read cor and p

G141_cor_df_padj <- read.table("1_0.3_edge_G141.tsv", sep='\t', header =T)

nodeattrib_G141_combine <- data.frame(node=union(G141_cor_df_padj$Source,G141_cor_df_padj$Target))
nodeattrib_G141_combine$kingdom <- 0

for (i in as.character(nodeattrib_G141_combine$node))
{
  if (i %in% bac.list$OTU_id == TRUE)
  {nodeattrib_G141_combine[nodeattrib_G141_combine$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattrib_G141_combine[nodeattrib_G141_combine$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_G141_combine) <- as.character(nodeattrib_G141_combine$node)
nodeattrib_G141_combine$kingdom



all_G141_net <- graph_from_data_frame(G141_cor_df_padj,direct=F, vertices=nodeattrib_G141_combine)

## Number of nodes
length(V(all_G141_net)) # 72

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_G141_net)))) #38
length(grep("^F",names(V(all_G141_net)))) #34


## Connections 
bb_occur_G141 <- droplevels(G141_cor_df_padj[with(G141_cor_df_padj, grepl("^B",Source) & grepl("^B",Target)),])
nrow(bb_occur_G141) #51

ff_occur_G141 <- droplevels(G141_cor_df_padj[with(G141_cor_df_padj, grepl("^F",Source) & grepl("^F",Target)),])
nrow(ff_occur_G141) #22

fb_occur_G141 <- droplevels(G141_cor_df_padj[with(G141_cor_df_padj, grepl("^F",Source) & grepl("^B",Target)),])
nrow(fb_occur_G141) #35



## Network properties
meta_degree <- sort(igraph::degree(all_G141_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "41"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "4.04850213980029"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_G141_net, directed = FALSE))) #"average shortest path length =  3.59523728700828"
#"average shortest path length =  5.78163545837797"
print(paste("mean clustering coefficient = ", transitivity(all_G141_net, "global"))) #"mean clustering coefficient =  0.949977913375762"
#"mean clustering coefficient =  0.270780363778604"
print(paste("mean betweenness centrality = ", mean(betweenness(all_G141_net, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
#"mean betweenness centrality = 0.00606465960717704"
print(paste("mean closeness centrality = ", mean(closeness(all_G141_net, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
#"mean closeness centrality =  0.0203369303504658"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_G141_net, V(all_G141_net)))))) #"mean number of neighbors =  223.182421227197"
#"mean number of neighbors = 4.04850213980029"
##

net <- all_G141_net
G141_all_deg <- igraph::degree(net,mode="all")
G141_all_betweenness <- betweenness(net, normalized = TRUE)
G141_all_closeness <- closeness(net, normalized = TRUE)
G141_all_transitivity <- transitivity(net, "local", vids = V(net))
names(G141_all_transitivity)<- V(net)$name
G141_all_transitivity[is.na(G141_all_transitivity)] <- 0


## Defining hub OTUs
n <- 3
G141_all_deg.1percent <- G141_all_deg[G141_all_deg >= quantile(G141_all_deg,prob=1-n/100)]
length(G141_all_deg.1percent) #8

G141_all_betweenness.1percent <- G141_all_betweenness[G141_all_betweenness >= quantile(G141_all_betweenness,prob=1-n/100)]
length(G141_all_betweenness.1percent) #7

G141_all_closeness.1percent <- G141_all_closeness[G141_all_closeness >= quantile(G141_all_closeness,prob=1-n/100)]
length(G141_all_closeness.1percent) #7

intersect(names(G141_all_deg.1percent), names(G141_all_betweenness.1percent))
intersect(names(G141_all_deg.1percent), names(G141_all_closeness.1percent))
intersect(names(G141_all_betweenness.1percent), names(G141_all_closeness.1percent))
#F87_Penicillium (ITS2 network)
#B17_f_Peptostreptococcaceae (ITS1 network)

### network properties of bacteria and fungi
df.G141.degree<-data.frame(G141_all_deg)
head(df.G141.degree)
df.G141.degree$Group <- "G141"
names(df.G141.degree)[1] <- c("Degree")

df.G141.closeness<-data.frame(G141_all_closeness)
head(df.G141.closeness)
df.G141.closeness$Group <- "G141"
names(df.G141.closeness)[1] <- c("Closeness")

df.G141.betweenness<-data.frame(G141_all_betweenness)
head(df.G141.betweenness)
df.G141.betweenness$Group <- "G141"
names(df.G141.betweenness)[1] <- c("Betweenness")

df.G141.degree$kingdom <- ifelse(grepl("^B",rownames(df.G141.degree)),'Bacteria', 'Fungi')
df.G141.closeness$kingdom <- ifelse(grepl("^B",rownames(df.G141.closeness)),'Bacteria', 'Fungi')
df.G141.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.G141.betweenness)),'Bacteria', 'Fungi')


print_degree <- function(deg){
  
  deg$kingdom <- factor(deg$kingdom, levels=c('Bacteria', 'Fungi'))  
  max.degree <- aggregate(deg$Degree, by = list(deg$kingdom), max)
  
  colnames(max.degree) <- c("kingdom", "maxdegree")
  
  # wilcoxon test
  x <- subset(deg, kingdom=='Bacteria')$Degree
  y <- subset(deg, kingdom=='Fungi')$Degree
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Degree \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}

p1 <- print_degree(df.G141.degree)
p1


## Betweenness centrality
print_betweenness <- function(betw){
  
  betw$kingdom <- factor(betw$kingdom, levels=c('Bacteria','Fungi'))  
  max.betweenness <- aggregate(betw$Betweenness, by = list(betw$kingdom), max)
  
  colnames(max.betweenness) <- c("kingdom", "maxbetweenness")
  
  
  # wilcoxon test
  x <- subset(betw, kingdom=='Bacteria')$Betweenness
  y <- subset(betw, kingdom=='Fungi')$Betweenness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Betweenness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p2 <- print_betweenness(df.G141.betweenness)
p2

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8, ncol=4)


## Closeness centrality
print_closeness <- function(closen){
  
  closen$kingdom <- factor(closen$kingdom, levels=c('Bacteria','Fungi'))  
  max.closeness <- aggregate(closen$Closeness, by = list(closen$kingdom), max)
  
  colnames(max.closeness) <- c("kingdom", "maxcloseness")
  
  
  x <- subset(closen, kingdom=='Bacteria')$Closeness
  y <- subset(closen, kingdom=='Fungi')$Closeness
  wilcox.test(x, y, conf.int = TRUE)
  
  p<-ggplot(data=closen, aes(x=kingdom, y=Closeness, fill=kingdom))+ geom_boxplot(width = 0.8) +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
    ylab("Closeness centrality \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
    scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(aspect.ratio = 2) + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  return(p)
}
p3 <- print_closeness(df.G141.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()



### Plotting hub
head(df.G141.degree)
head(df.G141.betweenness)
df.G141.degree$OTU_id <- rownames(df.G141.degree)
df.G141.betweenness$OTU_id <- rownames(df.G141.betweenness)
df.G141.closeness$OTU_id <- rownames(df.G141.closeness)
df.G141.net.properties <- merge(df.G141.degree, df.G141.betweenness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G141.net.properties <- merge(df.G141.net.properties, df.G141.closeness, by = c("OTU_id" = "OTU_id","Group" = "Group", "kingdom" = "kingdom"))
df.G141.net.properties$kingdom <- factor(df.G141.net.properties$kingdom, levels = c('Bacteria', 'Fungi'))
n<-3
ggplot(df.G141.net.properties, aes(x=Degree, y=Betweenness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Betweenness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G141.net.properties$Betweenness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G141.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)

ggplot(df.G141.net.properties, aes(x=Degree, y=Closeness, color=kingdom)) +
  xlab('\n Degree')+
  ylab("Closeness centrality \n") +
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
  scale_color_manual(values=c("Bacteria"="#c7bfab","Fungi"="#83a259")) + 
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
  geom_hline(yintercept=quantile(df.G141.net.properties$Closeness,prob=1-n/100), linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=quantile(df.G141.net.properties$Degree,prob=1-n/100), linetype='dashed', color='black', size = 0.75)



###
df.G0.net.properties
df.G76.net.properties
df.G90.net.properties
df.G106.net.properties
df.G120.net.properties
df.G141.net.properties


df.G.net.properties <- rbind(df.G0.net.properties,df.G76.net.properties, df.G90.net.properties,
                             df.G106.net.properties,df.G120.net.properties,df.G141.net.properties)


#### plotting
library(ggpubr)
p <- ggboxplot(data = df.G.net.properties, x="Group", y="Degree", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Degree\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = kingdom), label = "p.signif", method = "wilcox.test")
p

p <- ggboxplot(data = df.G.net.properties, x="Group", y="Betweenness", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Betweenness centrality\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = kingdom), label = "p.signif", method = "wilcox.test")
p

p <- ggboxplot(data = df.G.net.properties, x="Group", y="Betweenness", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Betweenness centrality\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = kingdom), label = "p.signif", method = "wilcox.test")
p


p <- ggboxplot(data = df.G.net.properties, x="Group", y="Closeness", fill = "kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Closeness centrality\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = kingdom), label = "p.signif", method = "wilcox.test")
p



### Network complexity
##### Measuring network complexity
install.packages("QuACN")
BiocManager::install("RBGL")
library(QuACN)

#### 
g.G0<-igraph.to.graphNEL(all_G0_net)
g.G76<-igraph.to.graphNEL(all_G76_net)
g.G90<-igraph.to.graphNEL(all_G90_net)
g.G106<-igraph.to.graphNEL(all_G106_net)
g.G120<-igraph.to.graphNEL(all_G120_net)
g.G141<-igraph.to.graphNEL(all_G141_net)

#### Bertz complexity index
bertz(g.G0) #976.7536
bertz(g.G76) #6151.023
bertz(g.G90) #2841.116
bertz(g.G106) #1370.054
bertz(g.G120) #855.7143
bertz(g.G141) #866.4692


#### Offdiagonal Complexity
offdiagonal(g.G0)#0.4584748
offdiagonal(g.G76) # 0.515314
offdiagonal(g.G90) #0.4091284
offdiagonal(g.G106) #0.3968748
offdiagonal(g.G120) #0.4528651
offdiagonal(g.G141) #0.3784092


##Atom-bond connectivity
atomBondConnectivity(g.G0) #88.1586
atomBondConnectivity(g.G76)#788.624
atomBondConnectivity(g.G90)#212.8144
atomBondConnectivity(g.G106)#112.6028
atomBondConnectivity(g.G120)#76.8729
atomBondConnectivity(g.G141)#69.10818


##Geometric-arithmetic Indices
geometricArithmetic1(g.G0)#141.7909
geometricArithmetic1(g.G76)#1918.323
geometricArithmetic1(g.G90)#314.289
geometricArithmetic1(g.G106)#173.886
geometricArithmetic1(g.G120)#120.803
geometricArithmetic1(g.G141)#102.7638



### Contribution of core OTUs in seed networks
bac.seed.id <- subset(bac.list, OTU %in% bac.seed)
fun.seed.id <- subset(fun.list, OTU %in% fun.seed)

bac.seed.id <- bac.seed.id$OTU_id
fun.seed.id <- fun.seed.id$OTU_id

core.seed.id <- c(bac.seed.id, fun.seed.id)

df.G.net.properties.seed.core <- subset(df.G.net.properties, OTU_id %in% bac.seed.id |OTU_id %in% fun.seed.id)


### Proportion of core OTUs in network nodes
length(intersect(core.seed.id,names(V(all_G0_net))))/length(V(all_G0_net)) #26.25%
length(intersect(core.seed.id,names(V(all_G76_net))))/length(V(all_G76_net))#7.16%
length(intersect(core.seed.id,names(V(all_G90_net))))/length(V(all_G90_net)) #13%
length(intersect(core.seed.id,names(V(all_G106_net))))/length(V(all_G106_net))#21.6%
length(intersect(core.seed.id,names(V(all_G120_net))))/length(V(all_G120_net))#30.5%
length(intersect(core.seed.id,names(V(all_G141_net))))/length(V(all_G141_net))#27.8%

### Composition of edges
draw_edge_proportion <- function(filename){
  edge.mer <- read.table(file = filename, sep = '\t', header = TRUE)
  head(edge.mer)
  
  edge.mer$BF <- ifelse(grepl("^B",edge.mer$Source) & grepl("^F",edge.mer$Target) | grepl("^F",edge.mer$Source) & grepl("^B",edge.mer$Target),'BF',ifelse(grepl("^B",edge.mer$Source) & grepl("^B",edge.mer$Target),'B','F'))
  edge.mer
  
  edge.mer$PN <- ifelse(edge.mer$positive,'positive','negative')
  edge.mer$count <- 1
  edge.mer %>% arrange(desc(Cor))
  edge.mer %>% arrange(desc(BF))
  
  df.edge.mer <- edge.mer %>% group_by(BF,PN) %>% summarise(count=sum(count))
  
  df.edge.mer$Proportion <- df.edge.mer$count / sum(df.edge.mer$count)
  df.edge.mer
  print(df.edge.mer %>% filter(BF == 'B') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(BF == 'BF') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(BF == 'F') %>% mutate(ratio = Proportion / sum(Proportion)))
  
  p.edge.mer <- ggplot(df.edge.mer, aes(x=BF, y = Proportion, fill = PN)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack', colour="black") +
    #scale_fill_discrete() +
    scale_fill_manual(values = c('red2','blue2')) +
    
    xlab('')+
    ylab("Proportion of edges \n") +
    #ggtitle("Phylum Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 2,reverse = T))+theme(aspect.ratio = 1.5)+
    theme(legend.position="bottom",legend.title=element_blank()) +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(),plot.background=element_blank())
  
  return(p.edge.mer)
}


draw_edge_proportion("1_0.3_edge_G0.tsv")
draw_edge_proportion("1_0.3_edge_G76.tsv")
draw_edge_proportion("1_0.3_edge_G90.tsv")
draw_edge_proportion("1_0.3_edge_G106.tsv")
draw_edge_proportion("1_0.3_edge_G120.tsv")
draw_edge_proportion("1_0.3_edge_G141.tsv")


### edge count of Core OTU
edge.G0 <- read.table(file = "1_0.3_edge_G0.tsv", sep = '\t', header = TRUE)
head(edge.G0)

edge.G0$core_association <- ifelse(edge.G0$Source %in% core.seed.id & edge.G0$Target %in% core.seed.id,'CC',ifelse(!(edge.G0$Target %in% core.seed.id) & !(edge.G0$Source %in% core.seed.id),'NN','CN'))
edge.G0$from_to <- paste(edge.G0$Source,"_to_",edge.G0$Target)
length(edge.G0$from_to[which(edge.G0$core_association=="CC")]) #15
length(edge.G0$from_to[which(edge.G0$core_association=="CN")]) #55
length(edge.G0$from_to[which(edge.G0$core_association=="NN")]) #81

edge.G76 <- read.table(file = "1_0.3_edge_G76.tsv", sep = '\t', header = TRUE)
head(edge.G76)

edge.G76$core_association <- ifelse(edge.G76$Source %in% core.seed.id & edge.G76$Target %in% core.seed.id,'CC',ifelse(!(edge.G76$Target %in% core.seed.id) & !(edge.G76$Source %in% core.seed.id),'NN','CN'))
edge.G76
edge.G76$from_to <- paste(edge.G76$Source,"_to_",edge.G76$Target)
length(edge.G76$from_to[which(edge.G76$core_association=="CC")]) #18
length(edge.G76$from_to[which(edge.G76$core_association=="CN")]) #259
length(edge.G76$from_to[which(edge.G76$core_association=="NN")]) #1756

edge.G90 <- read.table(file = "1_0.3_edge_G90.tsv", sep = '\t', header = TRUE)
head(edge.G90)

edge.G90$core_association <- ifelse(edge.G90$Source %in% core.seed.id & edge.G90$Target %in% core.seed.id,'CC',ifelse(!(edge.G90$Target %in% core.seed.id) & !(edge.G90$Source %in% core.seed.id),'NN','CN'))
edge.G90
edge.G90$from_to <- paste(edge.G90$Source,"_to_",edge.G90$Target)
length(edge.G90$from_to[which(edge.G90$core_association=="CC")]) #22
length(edge.G90$from_to[which(edge.G90$core_association=="CN")]) #62
length(edge.G90$from_to[which(edge.G90$core_association=="NN")]) #258

edge.G106 <- read.table(file = "1_0.3_edge_G106.tsv", sep = '\t', header = TRUE)
head(edge.G106)

edge.G106$core_association <- ifelse(edge.G106$Source %in% core.seed.id & edge.G106$Target %in% core.seed.id,'CC',ifelse(!(edge.G106$Target %in% core.seed.id) & !(edge.G106$Source %in% core.seed.id),'NN','CN'))
edge.G106
edge.G106$from_to <- paste(edge.G106$Source,"_to_",edge.G106$Target)
length(edge.G106$from_to[which(edge.G106$core_association=="CC")]) #15
length(edge.G106$from_to[which(edge.G106$core_association=="CN")]) #68
length(edge.G106$from_to[which(edge.G106$core_association=="NN")]) #100


edge.G120 <- read.table(file = "1_0.3_edge_G120.tsv", sep = '\t', header = TRUE)
head(edge.G120)

edge.G120$core_association <- ifelse(edge.G120$Source %in% core.seed.id & edge.G120$Target %in% core.seed.id,'CC',ifelse(!(edge.G120$Target %in% core.seed.id) & !(edge.G120$Source %in% core.seed.id),'NN','CN'))
edge.G120
edge.G120$from_to <- paste(edge.G120$Source,"_to_",edge.G120$Target)
length(edge.G120$from_to[which(edge.G120$core_association=="CC")]) #18
length(edge.G120$from_to[which(edge.G120$core_association=="CN")]) #52
length(edge.G120$from_to[which(edge.G120$core_association=="NN")]) #58

edge.G141 <- read.table(file = "1_0.3_edge_G141.tsv", sep = '\t', header = TRUE)
head(edge.G141)

edge.G141$core_association <- ifelse(edge.G141$Source %in% core.seed.id & edge.G141$Target %in% core.seed.id,'CC',ifelse(!(edge.G141$Target %in% core.seed.id) & !(edge.G141$Source %in% core.seed.id),'NN','CN'))
edge.G141
edge.G141$from_to <- paste(edge.G141$Source,"_to_",edge.G141$Target)
length(edge.G141$from_to[which(edge.G141$core_association=="CC")]) #13
length(edge.G141$from_to[which(edge.G141$core_association=="CN")]) #45
length(edge.G141$from_to[which(edge.G141$core_association=="NN")]) #50

cc.G141<-edge.G141$from_to[which(edge.G141$core_association=="CC")]
cc.G120<-edge.G120$from_to[which(edge.G120$core_association=="CC")]
intersect(cc.G141,cc.G120)


cc.G106<-edge.G106$from_to[which(edge.G106$core_association=="CC")]
cc.G90<-edge.G90$from_to[which(edge.G90$core_association=="CC")]
intersect(cc.G141,cc.G106)
intersect(cc.G141,cc.G90)

cc.G76<-edge.G76$from_to[which(edge.G76$core_association=="CC")]
cc.G0<-edge.G0$from_to[which(edge.G0$core_association=="CC")]
intersect(cc.G141,cc.G0)
intersect(cc.G141,cc.G76)
