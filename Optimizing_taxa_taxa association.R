### Taxonomic composition of associations
# Over-representation
require(ggplot2)
require(igraph)
require(corrplot)


# otu.list <- rbind(bac.list, fun.list)
# write.csv(otu.list, "otu list of all communities.csv")

otu.list <- read.csv("otu_list_of_all_communities.csv")
### unidentified 

otu.list.2<- tidyr::unite(otu.list, "Class2", Phylum, Class, sep = "_")
otu.list.3<- tidyr::unite(otu.list.2, "Order", Class2, Order, sep = "_")
otu.list.4<- tidyr::unite(otu.list.3, "Family", Order, Family, sep = "_")
otu.list.5<- tidyr::unite(otu.list.4, "Genus", Family, Genus, sep = "_")
write.csv(otu.list.5,"otu_table_genus_merged.csv")

all.asso <- read.csv("Associations_of_all_networks_final_2.csv")

node.all.net <- read.csv('All_network_node_2.csv')

vertex.phylum <- otu.list.2$Class2[which(otu.list.2$OTU_id %in% unique(node.all.net$id))]
length(unique(vertex.phylum))


##Defining dominant classes


###Bacteria
bac.clean.ss.f.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.f.18  <- phyloseq::filter_taxa(bac.clean.ss.f.18 , function(x) sum(x) != 0, TRUE)


df.phylum <- bac.clean.ss.f.18  %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.phylum$Age <- factor(df.phylum$Age, levels = order.age)



library(forcats) 
df.phylum %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
unique(df.phylum$Class)


df.phylum.2<- tidyr::unite(df.phylum, "Class2", Phylum, Class, sep = "_")


levels(df.phylum.2$Class2)
levels(df.phylum.2$Class2) = c(levels(df.phylum.2$Class2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 2,]$Class2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Class2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(desc(Abundance))
vec <- ord$Class2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]

dom.order.bac <- vec.order


###Fungi
fun.clean.ss.f.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.f.18  <- phyloseq::filter_taxa(fun.clean.ss.f.18 , function(x) sum(x) != 0, TRUE)


df.phylum <- fun.clean.ss.f.18  %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.phylum$Age <- factor(df.phylum$Age, levels = order.age)

df.phylum %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
df.phylum %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))

df.phylum.2 <-df.phylum
df.phylum.2$Class2<- paste(df.phylum.2$Phylum, df.phylum.2$Class, sep = "_")

levels(df.phylum.2$Class2)
levels(df.phylum.2$Class2) = c(levels(df.phylum.2$Class2), 'Low abundance')


# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 2,]$Class2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Class2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(desc(Abundance))
vec <- ord$Class2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]

dom.order.fun <- vec.order

dom.class <- c(dom.order.bac, dom.order.fun)

as.character(dom.class) -> dom.class
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


##all
corrplot(corr = 1-p.over.adj,
         #order = "hclust",
         #hclust.method = "single",
         col = getPalette(1),
         tl.col = "black",
         tl.cex=.6,
         #tl.pos = "ld",
         cl.pos = "n",
         method = "square",
         insig = "blank",
         sig.level = 0.05, 
         p.mat = p.over.adj,
         type = "lower"
)
dev.off()

#### Analysis in each compartment

all.asso <- read.csv("Associations_of_all_networks_final_2.csv")

### Association
BS.net.asso <- subset(all.asso, Compartment == "BS")
RS.net.asso <- subset(all.asso, Compartment == "RS")
R.net.asso <- subset(all.asso, Compartment == "R")
S1.net.asso <- subset(all.asso, Compartment == "S1")
S2.net.asso <- subset(all.asso, Compartment == "S2")
S3.net.asso <- subset(all.asso, Compartment == "S3")
S4.net.asso <- subset(all.asso, Compartment == "S4")
S5.net.asso <- subset(all.asso, Compartment == "S5")
S6.net.asso <- subset(all.asso, Compartment == "S6")
S7.net.asso <- subset(all.asso, Compartment == "S7")
S8.net.asso <- subset(all.asso, Compartment == "S8")
S9.net.asso <- subset(all.asso, Compartment == "S9")
L1.net.asso <- subset(all.asso, Compartment == "L1")
L2.net.asso <- subset(all.asso, Compartment == "L2")
L3.net.asso <- subset(all.asso, Compartment == "L3")
FL.net.asso <- subset(all.asso, Compartment == "FL")
G.net.asso <- subset(all.asso, Compartment == "G")


##node
node.all.net <- read.csv('All_network_node_2.csv')

BS.net.node <- subset(node.all.net, Compartment == "BS")
RS.net.node <- subset(node.all.net, Compartment == "RS")
R.net.node <- subset(node.all.net, Compartment == "R")
S1.net.node <- subset(node.all.net, Compartment == "S1")
S2.net.node <- subset(node.all.net, Compartment == "S2")
S3.net.node <- subset(node.all.net, Compartment == "S3")
S4.net.node <- subset(node.all.net, Compartment == "S4")
S5.net.node <- subset(node.all.net, Compartment == "S5")
S6.net.node <- subset(node.all.net, Compartment == "S6")
S7.net.node <- subset(node.all.net, Compartment == "S7")
S8.net.node <- subset(node.all.net, Compartment == "S8")
S9.net.node <- subset(node.all.net, Compartment == "S9")
L1.net.node <- subset(node.all.net, Compartment == "L1")
L2.net.node <- subset(node.all.net, Compartment == "L2")
L3.net.node <- subset(node.all.net, Compartment == "L3")
FL.net.node <- subset(node.all.net, Compartment == "FL")
G.net.node <- subset(node.all.net, Compartment == "G")



###Optimized
hypergeometric.analysis<-function(net.node, net.asso,keyword){
  vertex.class <- otu.list.2$Class2[which(as.character(otu.list.2$OTU_id) %in% as.character(unique(net.node$id)))]

G0_cor_df_padj <- net.asso

G0_cor_df_padj$Source_Class <- 0

for (i in as.character(G0_cor_df_padj$Source)){
  G0_cor_df_padj$Source_Class[which(G0_cor_df_padj$Source == i)] <- as.character(otu.list.2$Class2[which(otu.list.2$OTU_id ==i)])
  
}


G0_cor_df_padj$Target_Class <- 0
for (i in as.character(G0_cor_df_padj$Target)){
  G0_cor_df_padj$Target_Class[which(G0_cor_df_padj$Target == i)] <- as.character(otu.list.2$Class2[which(otu.list.2$OTU_id ==i)])
  
}

G0_cor_df_padj$Class_link <- paste(G0_cor_df_padj$Source_Class,"--",G0_cor_df_padj$Target_Class)

overrep <- function(p1 = phyl1,
                    p2 = phyl2){
  if (p1 != p2){
    v1 = table(vertex.class)[p1]
    v2 = table(vertex.class)[p2]
    m <- v1*v2*(length(unique(net.node$id))-1)/length(unique(net.node$id))
    n <- length(unique(net.node$id))*(length(unique(net.node$id))-1) - m
    k <- length(G0_cor_df_padj$Source)
    x <- length(G0_cor_df_padj$Class_link[which(G0_cor_df_padj$Class_link == paste(p1,"--",p2))])
    p = phyper(x, m, n, k)}
  else {
    v1 = table(vertex.class)[p1]
    m <- v1*c(v1-1)*(length(unique(net.node$id))-1)/length(unique(net.node$id))
    n <- length(unique(net.node$id))*(length(unique(net.node$id))-1) - m
    k <- length(G0_cor_df_padj$Source)
    x <- length(G0_cor_df_padj$Class_link[which(G0_cor_df_padj$Class_link == paste(p1,"--",p2))])
    p = phyper(x, m, n, k)
  }
}

p.over <- matrix(nrow = length(unique(vertex.class)),ncol =  length(unique(vertex.class)))
for(i in 1: length(unique(vertex.class))){
  for(j in 1: length(unique(vertex.class))){
    po <- overrep(p1 = names(sort(table(vertex.class),decreasing = TRUE))[i],
                  p2 = names(sort(table(vertex.class),decreasing = TRUE))[j])
    p.over[i,j] <- po
  }
  print(i)
}

p.over.adj <- matrix(p.adjust(p.over,method = "fdr"),nrow =  length(unique(vertex.class)))
require(corrplot)

row.names(p.over.adj) <- names(sort(table(vertex.class),decreasing = TRUE))
colnames(p.over.adj) <- names(sort(table(vertex.class),decreasing = TRUE))

write.csv(p.over.adj,paste0("P value_Overrepresentation_meta network_",keyword,".csv"))

return(p.over.adj)}


### All networks
result.all<-hypergeometric.analysis(net.node = node.all.net,net.asso =all.asso,keyword = "all")

### each compartment

result.BS <- hypergeometric.analysis(net.node = BS.net.node,net.asso =BS.net.asso,keyword = "BS")

result.RS <- hypergeometric.analysis(net.node = RS.net.node,net.asso =RS.net.asso,keyword = "RS")

result.R <- hypergeometric.analysis(net.node = R.net.node,net.asso =R.net.asso,keyword = "R")

result.S1 <- hypergeometric.analysis(net.node = S1.net.node,net.asso =S1.net.asso,keyword = "S1")

result.S2 <- hypergeometric.analysis(net.node = S2.net.node,net.asso =S2.net.asso,keyword = "S2")

result.S3 <- hypergeometric.analysis(net.node = S3.net.node,net.asso =S3.net.asso,keyword = "S3")

result.S4 <- hypergeometric.analysis(net.node = S4.net.node,net.asso =S4.net.asso,keyword = "S4")

result.S5 <- hypergeometric.analysis(net.node = S5.net.node,net.asso =S5.net.asso,keyword = "S5")

result.S6 <- hypergeometric.analysis(net.node = S6.net.node,net.asso =S6.net.asso,keyword = "S6")

result.S7 <- hypergeometric.analysis(net.node = S7.net.node,net.asso =S7.net.asso,keyword = "S7")

result.S8 <- hypergeometric.analysis(net.node = S8.net.node,net.asso =S8.net.asso,keyword = "S8")

result.S9 <- hypergeometric.analysis(net.node = S9.net.node,net.asso =S9.net.asso,keyword = "S9")

result.L1 <- hypergeometric.analysis(net.node = L1.net.node,net.asso =L1.net.asso,keyword = "L1")

result.L2 <- hypergeometric.analysis(net.node = L2.net.node,net.asso =L2.net.asso,keyword = "L2")

result.L3 <- hypergeometric.analysis(net.node = L3.net.node,net.asso =L3.net.asso,keyword = "L3")

result.FL <- hypergeometric.analysis(net.node = FL.net.node,net.asso =FL.net.asso,keyword = "FL")

result.G <- hypergeometric.analysis(net.node = G.net.node,net.asso =G.net.asso,keyword = "G")

##dominant


dominant.class.plot <- function(resultfile){ 
  
p.over.adj <- resultfile

common.class<-intersect(colnames(resultfile),dom.class)

corrplot(corr = 1-p.over.adj[colnames(p.over.adj) %in% common.class,],
         #order = "hclust",
         #hclust.method = "single",
         col = getPalette(1),
         tl.col = "black",
         tl.cex=.6,
         #tl.pos = "ld",
         cl.pos = "n",
         method = "square",
         insig = "blank",
         sig.level = 0.05, 
         p.mat = p.over.adj[colnames(p.over.adj) %in% common.class,],
         type = "lower"
)

}
dominant.class.plot(result.BS)
dominant.class.plot(result.RS)
dominant.class.plot(result.R)
dominant.class.plot(result.S1)
dominant.class.plot(result.S2)
dominant.class.plot(result.S3)
dominant.class.plot(result.S4)
dominant.class.plot(result.S5)
dominant.class.plot(result.S6)
dominant.class.plot(result.S7)
dominant.class.plot(result.S8)
dominant.class.plot(result.S9)
dominant.class.plot(result.L1)
dominant.class.plot(result.L2)
dominant.class.plot(result.L3)
dominant.class.plot(result.FL)
dominant.class.plot(result.G)
dominant.class.plot(result.all)


##all
all.class.plot <- function(resultfile){ 
  
  p.over.adj <- resultfile
  corrplot(corr = 1-p.over.adj,
         #order = "hclust",
         #hclust.method = "single",
         col = getPalette(1),
         tl.col = "black",
         tl.cex=.6,
         #tl.pos = "ld",
         cl.pos = "n",
         method = "square",
         insig = "blank",
         sig.level = 0.05, 
         p.mat = p.over.adj,
         type = "lower"
  )}


all.class.plot(result.all)


dev.off()


hyper.genus <-read.csv("Overrepresentation_meta network_at_the_genus_level.csv")


###Dominant genus

###Bacteria
bac.clean.ss.f.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.f.18  <- phyloseq::filter_taxa(bac.clean.ss.f.18 , function(x) sum(x) != 0, TRUE)


df.phylum <- bac.clean.ss.f.18  %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.phylum$Age <- factor(df.phylum$Age, levels = order.age)



library(forcats) 

df.phylum.2<- tidyr::unite(df.phylum, "Class", Phylum, Class, sep = "_")
df.phylum.2<- tidyr::unite(df.phylum.2, "Order", Class, Order, sep = "_")
df.phylum.2<- tidyr::unite(df.phylum.2, "Family", Order, Family, sep = "_")
df.phylum.2<- tidyr::unite(df.phylum.2, "Genus", Family, Genus, sep = "_")

levels(df.phylum.2$Genus)
levels(df.phylum.2$Genus) = c(levels(df.phylum.2$Genus), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 2,]$Genus <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(desc(Abundance))
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]

dom.order.bac <- vec.order


###Fungi
fun.clean.ss.f.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.f.18  <- phyloseq::filter_taxa(fun.clean.ss.f.18 , function(x) sum(x) != 0, TRUE)


df.phylum <- fun.clean.ss.f.18  %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.phylum$Age <- factor(df.phylum$Age, levels = order.age)

df.phylum %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
df.phylum %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))

df.phylum.2 <-df.phylum
df.phylum.2$Genus<- paste(df.phylum.2$Phylum, df.phylum.2$Class, sep = "_")

levels(df.phylum.2$Genus)
levels(df.phylum.2$Genus) = c(levels(df.phylum.2$Genus), 'Low abundance')


# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 2,]$Genus <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(desc(Abundance))
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]

dom.order.fun <- vec.order

dom.class <- c(dom.order.bac, dom.order.fun)

as.character(dom.class) -> dom.class
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


# Provence 

dom.class <- df1$region
df.tax.link <- function(g = g.p1){
  el <- get.edgelist(g)
  a <- V(g.link)$class
  names(a) <- V(g.link)$name
  a.rank <- a[el[,1]]
  b.rank <- a[el[,2]]
  ab.edge <- paste(a.rank,b.rank,sep=" ")
  ab.edge <- ab.edge[a.rank %in% dom.class & b.rank %in% dom.class]
  tax.df <- NULL
  for(i in 1:15){
    edge.spe <- ab.edge[grep(dom.class[i],ab.edge)]
    tax.link <- 1:15
    for(j in c(1:15)[-i]){
      tax.link[j] <-  length(grep(dom.class[j], edge.spe))
    }
    tax.link[i] <- length(grep(paste(dom.class[i],dom.class[i]),ab.edge))
    tax.df <- rbind(tax.df,tax.link)
  }
  return(tax.df)
} 

df.g.entire <- df.tax.link(g.link)
file.list <- list.files("graph_1/")

tiff(filename = "manuscript/overrep_subnet.tiff",
     type = "cairo",
     compression = "lzw",
     width = 4000,
     height = 2000,
     res = 500)
par(mfrow = c(2,7))
for(i in 1:14){
  g.p1 <- read.graph(paste("graph/",file.list[i],sep=""),"graphml")
  df.p1 <- df.tax.link(g.p1)
  overrep.p.mat <- matrix(nrow = 15,ncol = 15)
  for(j in 1:15){
    for(k in 1:15){
      p <- pbinom(q = df.p1[j,k], 
                  size = df.g.entire[i,j], 
                  prob = ecount(g.p1)/ecount(g.link),
                  lower.tail = FALSE)
      overrep.p.mat[j,k] <- p
    }}
  write.csv(x = overrep.p.mat,
            file = paste("overrepresentation/",
                         file.list[i],
                         sep=""))
  overrep.p.adj <- matrix(p.adjust(overrep.p.mat,method = "fdr"),nrow = 15)
  corrplot(corr = 1-overrep.p.adj,
           col = getPalette(1),
           tl.pos = "n",
           cl.pos = "n",
           #method = "square",
           insig = "blank",
           sig.level = 0.05,
           type = "lower",
           p.mat = overrep.p.adj
  )
  
  print(i)
}

dev.off()


vertex.sp

### Taxa-taxa combine
# Find genus association combination
BS.net.asso <- subset(all.asso, Compartment == "BS")
RS.net.asso <- subset(all.asso, Compartment == "RS")
R.net.asso <- subset(all.asso, Compartment == "R")
S1.net.asso <- subset(all.asso, Compartment == "S1")
S2.net.asso <- subset(all.asso, Compartment == "S2")
S3.net.asso <- subset(all.asso, Compartment == "S3")
S4.net.asso <- subset(all.asso, Compartment == "S4")
S5.net.asso <- subset(all.asso, Compartment == "S5")
S6.net.asso <- subset(all.asso, Compartment == "S6")
S7.net.asso <- subset(all.asso, Compartment == "S7")
S8.net.asso <- subset(all.asso, Compartment == "S8")
S9.net.asso <- subset(all.asso, Compartment == "S9")
L1.net.asso <- subset(all.asso, Compartment == "L1")
L2.net.asso <- subset(all.asso, Compartment == "L2")
L3.net.asso <- subset(all.asso, Compartment == "L3")
FL.net.asso <- subset(all.asso, Compartment == "FL")
G.net.asso <- subset(all.asso, Compartment == "G")



binomial.distribution.asso.compartment <- function(association.files, keyword, node.files){
  el.df <- subset(association.files, Compartment == keyword)
  head(el.df)
  
  G0_cor_df_padj <-el.df
  
  G0_cor_df_padj$Source_Genus <- 0
  
  for (i in as.character(G0_cor_df_padj$Source)){
    G0_cor_df_padj$Source_Genus[which(G0_cor_df_padj$Source == i)] <- as.character(otu.list.5$Genus[which(otu.list.5$OTU_id ==i)])
    
  }
  
  G0_cor_df_padj$Target <- as.character(G0_cor_df_padj$Target)
  G0_cor_df_padj$Target_Genus <- 0
  for (i in as.character(G0_cor_df_padj$Target)){
    G0_cor_df_padj$Target_Genus[which(G0_cor_df_padj$Target == i)] <- as.character(otu.list.5$Genus[which(otu.list.5$OTU_id ==i)])
    
  }
  
  G0_cor_df_padj$Genus_link <- paste(G0_cor_df_padj$Source_Genus,"--",G0_cor_df_padj$Target_Genus)
  
  
  net.genus<-subset(node.files, Compartment == keyword)
  
  
  vertex.sp <- otu.list.5$Genus[which(otu.list.5$OTU_id %in% unique(net.genus$id))]
  
  combine.1c <- NULL
  for (i in 1:length(G0_cor_df_padj$Genus_link)){
    a <- G0_cor_df_padj$Genus_link[i][which(G0_cor_df_padj$Source_Genus[i] != G0_cor_df_padj$Target_Genus[i])]
    combine.1c <- c(combine.1c,a)
  }
  
  combine <- sort(table(c(combine.1c)),decreasing = TRUE)
  combine.c <- combine[-grep("unidentified_unidentified_unidentified_unidentified_unidentified",names(combine))]
  
  # Overrepresentations
  overrep.sp <- function(p1 = sp1, 
                         p2 = sp2){
    vertex.sp <- otu.list.5$Genus[which(otu.list.5$OTU_id %in% unique(net.genus$id))]
    
    v1 = table(vertex.sp)[p1]
    v2 = table(vertex.sp)[p2]
    prob <- v1*v2/(length(unique(net.genus$id))*(length(unique(net.genus$id))-1))
    size <- length(G0_cor_df_padj$Source)
    q <- length(G0_cor_df_padj$Genus_link[which(G0_cor_df_padj$Genus_link == paste(p1,"--",p2))])
    p <- pbinom(q,size,prob,lower.tail = FALSE)
    return(p)
  }
  
  
  p.sp.combine <- NULL
  
  ###Without unidentified genus
  length(combine)
  length(combine.c)
  p1.sp <- unlist(base::strsplit(names(combine.c[1:length(combine)]),split = " -- "))[seq(1,(length(combine)*2),2)]
  p2.sp <- unlist(base::strsplit(names(combine.c[1:length(combine)]),split = " -- "))[seq(2,(length(combine)*2),2)]
  
  for(i in 1:length(combine)){
    p.sp <- overrep.sp(p1 = p1.sp[i],
                       p2 = p2.sp[i])
    p.sp.combine <- c(p.sp.combine,p.sp)
    print(i)
  }
  p.sp.adj <- p.adjust(p.sp.combine,method = "fdr")
  
  sp.combine <- names(combine[1:length(combine)])
  sp.combine.c <- sp.combine[p.sp.adj<0.05]
  
  sp.df <- data.frame(species1 = unlist(base::strsplit(sp.combine.c,split = " -- "))[seq(1,length(sp.combine.c)*2,2)],
                      species2 = unlist(base::strsplit(sp.combine.c,split = " -- "))[seq(2,length(sp.combine.c)*2,2)],
                      time = combine[1:length(combine)][p.sp.adj<0.05],
                      over.expression = p.sp.adj[p.sp.adj<0.05])
  head(sp.df )
  write.csv(sp.df, paste0("sp_combine_",keyword,".csv"))
  return(sp.df)
  
  
}



BS.combine<-binomial.distribution.asso.compartment(all.asso,"BS",node.all.net)
RS.combine<-binomial.distribution.asso.compartment(all.asso,"RS",node.all.net)
R.combine<-binomial.distribution.asso.compartment(all.asso,"R",node.all.net)

S1.combine<-binomial.distribution.asso.compartment(all.asso,"S1",node.all.net)
S2.combine<-binomial.distribution.asso.compartment(all.asso,"S2",node.all.net)
S3.combine<-binomial.distribution.asso.compartment(all.asso,"S3",node.all.net)
S4.combine<-binomial.distribution.asso.compartment(all.asso,"S4",node.all.net)
S5.combine<-binomial.distribution.asso.compartment(all.asso,"S5",node.all.net)
S6.combine<-binomial.distribution.asso.compartment(all.asso,"S6",node.all.net)
S7.combine<-binomial.distribution.asso.compartment(all.asso,"S7",node.all.net)
S8.combine<-binomial.distribution.asso.compartment(all.asso,"S8",node.all.net)
S9.combine<-binomial.distribution.asso.compartment(all.asso,"S9",node.all.net)

L1.combine<-binomial.distribution.asso.compartment(all.asso,"L1",node.all.net)
L2.combine<-binomial.distribution.asso.compartment(all.asso,"L2",node.all.net)
L3.combine<-binomial.distribution.asso.compartment(all.asso,"L3",node.all.net)
FL.combine<-binomial.distribution.asso.compartment(all.asso,"FL",node.all.net)

G.combine<-binomial.distribution.asso.compartment(all.asso,"G",node.all.net)


### Chord diagram
# Libraries
#install.packages(c('patchwork','hrbrthemes','circlize'))
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")

###Example
# # Load dataset from github
# data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)
# 
# # short names
# colnames(data) <- c("Africa", "East Asia", "Europe", "Latin Ame.",   "North Ame.",   "Oceania", "South Asia", "South East Asia", "Soviet Union", "West.Asia")
# rownames(data) <- colnames(data)
# 
# # I need a long format
# data_long <- data %>%
#   rownames_to_column %>%
#   gather(key = 'key', value = 'value', -rowname)
# 
# # parameters
# circos.clear()
# circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
# par(mar = rep(0, 4))
# 
# # color palette
# mycolor <- viridis(44, alpha = 1, begin = 0, end = 1, option = "D")
# mycolor <- mycolor[sample(1:44)]
# 
# # Base plot
# chordDiagram(
#   x = data_long,
#   grid.col = mycolor,
#   transparency = 0.25,
#   directional = 1,
#   direction.type = c("arrows", "diffHeight"),
#   diffHeight  = -0.04,
#   annotationTrack = "grid",
#   annotationTrackHeight = c(0.05, 0.1),
#   link.arr.type = "big.arrow",
#   link.sort = TRUE,
#   link.largest.ontop = TRUE)
# 
# # Add text and axis
# circos.trackPlotRegion(
#   track.index = 1,
#   bg.border = NA,
#   panel.fun = function(x, y) {
# 
#     xlim = get.cell.meta.data("xlim")
#     sector.index = get.cell.meta.data("sector.index")
# 
#     # Add names to the sector.
#     circos.text(
#       x = mean(xlim),
#       y = 3.2,
#       labels = sector.index,
#       facing = "bending",
#       cex = 0.8
#     )
# 
#     # Add graduation on axis
#     circos.axis(
#       h = "top",
#       major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)),
#       minor.ticks = 1,
#       major.tick.percentage = 0.5,
#       labels.niceFacing = FALSE)
#   }
# )


### Our data
dev.off()

microhabitat.list<- c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL","G")

for (i in as.character(microhabitat.list)){
  data_long <- read.csv(paste0("data_long_",i,".csv"), header=TRUE)
  
  fun.list$Phylum[is.na(fun.list$Phylum)] <- "unidentified"
  data_long$species1_kingdom <- ifelse(data_long$species1 %in% fun.list$Phylum, "Fungi","Bacteria")
  
  write.csv(data_long,paste0("data_long_",i,"_edit.csv"))
}

length(unique(c(data_long$species1,data_long$species2)))

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 3, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
# mycolor <- viridis(44, alpha = 1, begin = 0, end = 1, option = "D")
# mycolor <- mycolor[sample(1:44)]

mycolor <- c("Proteobacteria" = "#556B2F", "Ascomycota" = "#1E63AF", "Basidiomycota" = "#BE4146", "Chytridiomycota" = "#CCCC99",
             "Glomeromycota" = "#34608d","Monoblepharomycota" = "#CC0000","Mortierellomycota" = "#4E734E","Mucoromycota" = "#E4AF2C",
             "Rozellomycota" = "#1f9e89", "unidentified" = "#D3D3D3", "WPS-2" = "#2fb47c", "Tenericutes" = "#91d742",
             "WS2" = "#423f85", "TA06" = "#5ac864", "WS4" = "#46317e", "BRC1" = "#4ec36b","Zixibacteria" = "#efe51c",
             "Verrucomicrobia" = "#CD69C9", "Acidobacteria" = "#8B5742","Actinobacteria" = "#EE6363","Armatimonadetes" = "#FF7F50",
             "Bacteroidetes" = "#63B8FF", "Chlamydiae" = "#218e8d", "Chloroflexi" = "#FFD700", "Cyanobacteria" = "#458B00",
             "Deinococcus-Thermus" = "#2c718e", "Elusimicrobia" = "#00C5CD", "Epsilonbacteraeota" = "#528B8B","FBP" = "#287d8e",
             "Fibrobacteres" = "#EE7942", "Firmicutes" = "#FFA54F","Fusobacteria" = "#443883", "Gemmatimonadetes" = "#CDAF95",
             "Latescibacteria" = "#404688","Nitrospinae" = "#20a386","Nitrospirae" = "#EE799F","Omnitrophicaeota" ="#fde725",
             "Patescibacteria" ="#BFEFFF","Planctomycetes" ="#43CD80","Rokubacteria" = "#dfe318", "Spirochaetes" = "#8B7D6B",
             "Dependentiae" = "#ff7f0a", "Kiritimatiellaeota" = "#c9000b", "Hydrogenedentes" = "#ff917f")
# Base plot
dev.off()

data_long <- read.csv("data_long_BS_edit.csv")

chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.4,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = F, 
  link.largest.ontop = F)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # # Add graduation on axis
    # circos.axis(
    #   h = "top", 
    #   major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
    #   minor.ticks = 1, 
    #   major.tick.percentage = 0.5,
    #   labels.niceFacing = FALSE)
  }
)
