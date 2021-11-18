### Temporal changes in bray curtis distance of seeds with G0
bac.clean.log
fun.clean.log


### abundance table of shared seed communities
persistent.bac
persistent.fun

seed.sample.list <- map$SampleID[which(map$Microhabitat == "G")] 
seed.sample.list.2 <- f.meta$SampleID[which(f.meta$Microhabitat == "G")] 

#### Bray-Curtis dissimilarity - Bacteria
bac.clean.ss.rel <- transform(bac.clean.ss.f, transform = "compositional")

b.otu.lognorm <- otu_table(bac.clean.log) 

b.otu.lognorm <- otu_table(bac.clean.ss.rel) 

b.otu.lognorm.persist <- subset(b.otu.lognorm, rownames(b.otu.lognorm) %in% persistent.bac)
b.otu.lognorm.persist.t <- t(b.otu.lognorm.persist)
b.otu.lognorm.persist.t <- subset(b.otu.lognorm.persist.t,rownames(b.otu.lognorm.persist.t) %in% seed.sample.list)

bray.dist.bac<-vegdist(b.otu.lognorm.persist.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "Dissimilarity"
bray.dist.bac.melt <- subset(bray.dist.bac.melt, Dissimilarity != 0)

map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)


b.meta.seed <- subset(map, Compartment == "Seed" & Year == "year_2018")


f.meta<-sample_data(fun.clean.ss.f)
f.meta <- data.frame(f.meta)
f.meta.seed <- subset(f.meta, Compartment == "Grain" & Year == "year_2018")


G.0 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_0"))]
G.76 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_76"))]
G.90 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_90"))]
G.106 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_106"))]
G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_120"))]
G.141 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance between seed 76 and leaves
bray.dist.bac.76.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.76 & Var2 %in% G.0)
bray.dist.bac.76.G.2 <- subset(bray.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.76)
bray.dist.bac.76.G <- rbind(bray.dist.bac.76.G.1, bray.dist.bac.76.G.2)
length(bray.dist.bac.76.G$Dissimilarity)
bray.dist.bac.76.G$Comparison <- "G_76"
mean(bray.dist.bac.76.G$Dissimilarity) #0.8236615 #0.8418308 (RA)

bray.dist.bac.90.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.90 & Var2 %in% G.0)
bray.dist.bac.90.G.2 <- subset(bray.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.90)
bray.dist.bac.90.G <- rbind(bray.dist.bac.90.G.1, bray.dist.bac.90.G.2)
length(bray.dist.bac.90.G$Dissimilarity)
bray.dist.bac.90.G$Comparison <- "G_90"
mean(bray.dist.bac.90.G$Dissimilarity) #0.5700637 #0.5858067 (RA)

bray.dist.bac.106.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.106 & Var2 %in% G.0)
bray.dist.bac.106.G.2 <- subset(bray.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.106)
bray.dist.bac.106.G <- rbind(bray.dist.bac.106.G.1, bray.dist.bac.106.G.2)
length(bray.dist.bac.106.G$Dissimilarity)
bray.dist.bac.106.G$Comparison <- "G_106"
mean(bray.dist.bac.106.G$Dissimilarity) # 0.4586626 #0.4844428 (RA)

bray.dist.bac.120.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.120 & Var2 %in% G.0)
bray.dist.bac.120.G.2 <- subset(bray.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.120)
bray.dist.bac.120.G <- rbind(bray.dist.bac.120.G.1, bray.dist.bac.120.G.2)
length(bray.dist.bac.120.G$Dissimilarity)
bray.dist.bac.120.G$Comparison <- "G_120"
mean(bray.dist.bac.120.G$Dissimilarity) #0.4168994 #0.4338359 (RA)

bray.dist.bac.141.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.141 & Var2 %in% G.0)
bray.dist.bac.141.G.2 <- subset(bray.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.141)
bray.dist.bac.141.G <- rbind(bray.dist.bac.141.G.1, bray.dist.bac.141.G.2)
length(bray.dist.bac.141.G$Dissimilarity)
bray.dist.bac.141.G$Comparison <- "G_141"
mean(bray.dist.bac.141.G$Dissimilarity) #0.3919073 # 0.4121608 (RA)
 

bray.dist.bac.G <- rbind(bray.dist.bac.76.G,bray.dist.bac.90.G,bray.dist.bac.106.G,
                            bray.dist.bac.120.G,bray.dist.bac.141.G)



###Fungi
fun.clean.ss.rel <- transform(fun.clean.ss.f, transform = "compositional")

f.otu.lognorm <- otu_table(fun.clean.log) 

f.otu.lognorm <- otu_table(fun.clean.ss.rel) 

f.otu.lognorm.persist <- subset(f.otu.lognorm, rownames(f.otu.lognorm) %in% persistent.fun)
f.otu.lognorm.persist.t <- t(f.otu.lognorm.persist)
f.otu.lognorm.persist.t <- subset(f.otu.lognorm.persist.t,rownames(f.otu.lognorm.persist.t) %in% seed.sample.list.2)

bray.dist.fun<-vegdist(f.otu.lognorm.persist.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "Dissimilarity"
bray.dist.fun.melt <- subset(bray.dist.fun.melt, Dissimilarity != 0)




G.0 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_0"))]
G.76 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_76"))]
G.90 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_90"))]
G.106 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_106"))]
G.120 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_120"))]
G.141 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance between seed 76 and leaves
bray.dist.fun.76.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.76 & Var2 %in% G.0)
bray.dist.fun.76.G.2 <- subset(bray.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.76)
bray.dist.fun.76.G <- rbind(bray.dist.fun.76.G.1, bray.dist.fun.76.G.2)
length(bray.dist.fun.76.G$Dissimilarity)
bray.dist.fun.76.G$Comparison <- "G_76"
mean(bray.dist.fun.76.G$Dissimilarity) #0.7052215 #0.8744395 (RA)

bray.dist.fun.90.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.90 & Var2 %in% G.0)
bray.dist.fun.90.G.2 <- subset(bray.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.90)
bray.dist.fun.90.G <- rbind(bray.dist.fun.90.G.1, bray.dist.fun.90.G.2)
length(bray.dist.fun.90.G$Dissimilarity)
bray.dist.fun.90.G$Comparison <- "G_90"
mean(bray.dist.fun.90.G$Dissimilarity) #0.6340816 # 0.855687 (RA)

bray.dist.fun.106.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.106 & Var2 %in% G.0)
bray.dist.fun.106.G.2 <- subset(bray.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.106)
bray.dist.fun.106.G <- rbind(bray.dist.fun.106.G.1, bray.dist.fun.106.G.2)
length(bray.dist.fun.106.G$Dissimilarity)
bray.dist.fun.106.G$Comparison <- "G_106"
mean(bray.dist.fun.106.G$Dissimilarity) # 0.6211113 #0.8156643 (RA)

bray.dist.fun.120.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.120 & Var2 %in% G.0)
bray.dist.fun.120.G.2 <- subset(bray.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.120)
bray.dist.fun.120.G <- rbind(bray.dist.fun.120.G.1, bray.dist.fun.120.G.2)
length(bray.dist.fun.120.G$Dissimilarity)
bray.dist.fun.120.G$Comparison <- "G_120" 
mean(bray.dist.fun.120.G$Dissimilarity) #0.3933539 #0.6279877 (RA)

bray.dist.fun.141.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.141 & Var2 %in% G.0)
bray.dist.fun.141.G.2 <- subset(bray.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.141)
bray.dist.fun.141.G <- rbind(bray.dist.fun.141.G.1, bray.dist.fun.141.G.2)
length(bray.dist.fun.141.G$Dissimilarity)
bray.dist.fun.141.G$Comparison <- "G_141"
mean(bray.dist.fun.141.G$Dissimilarity) #0.4005389 #0.6654811 (RA)


bray.dist.fun.G <- rbind(bray.dist.fun.76.G,bray.dist.fun.90.G,bray.dist.fun.106.G,
                         bray.dist.fun.120.G,bray.dist.fun.141.G)


#### Kruskal-Wallis test
max.bray <- aggregate(bray.dist.bac.G$Dissimilarity, by = list(bray.dist.bac.G$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.bac.G$Comparison <- as.factor(bray.dist.bac.G$Comparison)
bray.dist.bac.G$Comparison <- factor(bray.dist.bac.G$Comparison, levels = c("G_76","G_90","G_106","G_120","G_141"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.G)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.G,
              method="bh")
PT = DT$res

#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- c("G_106", "G_120","G_141","G_76","G_90")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.bac.G, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p


##Fungi
#### Kruskal-Wallis test
max.bray <- aggregate(bray.dist.fun.G$Dissimilarity, by = list(bray.dist.fun.G$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.fun.G$Comparison <- as.factor(bray.dist.fun.G$Comparison)
bray.dist.fun.G$Comparison <- factor(bray.dist.fun.G$Comparison, levels = c("G_76","G_90","G_106","G_120","G_141"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.G)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.G,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn$Group <- c("G_106", "G_120","G_141","G_76","G_90")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.fun.G, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p


write.csv(bray.dist.bac.G,"Bray-curtis dissimilarity_bacterial_persistent_OTUs.csv")
write.csv(bray.dist.fun.G,"Bray-curtis dissimilarity_fungal_persistent_OTUs.csv")



#### Heritability of persistent OTUs
h2.bac.specificity.2 <-read.csv("Bac_h2_taxa_reassigned.csv")
h2.fun.specificity.2 <-read.csv("Fun_h2_taxa_reassigned.csv")

head(h2.bac.specificity.2)
h2.bac.G.persist <- subset(h2.bac.specificity.2, OTU %in% persistent.bac & Section == "G")
nrow(h2.bac.G.persist)

h2.fun.G.persist <- subset(h2.fun.specificity.2, OTU %in% persistent.fun & Section == "G")
nrow(h2.fun.G.persist)



#### Environmental effects or host physiology?
### persistent OTUs in seeds grown in Chuncheon (different site and weather condition)
seed.sample.list.3 <- map$SampleID[which(map$Compartment == "Seed")] 
seed.sample.list.4 <- f.meta$SampleID[which(f.meta$Compartment %in% c("Seed","Grain"))] 

#### Bray-Curtis dissimilarity - Bacteria
b.otu.lognorm <- otu_table(bac.clean.log)
b.otu.lognorm.persist <- subset(b.otu.lognorm, rownames(b.otu.lognorm) %in% persistent.bac)
b.otu.lognorm.persist.t <- t(b.otu.lognorm.persist)
b.otu.lognorm.persist.t <- subset(b.otu.lognorm.persist.t,rownames(b.otu.lognorm.persist.t) %in% seed.sample.list.3)

bray.dist.bac<-vegdist(b.otu.lognorm.persist.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac <-as.matrix(bray.dist.bac)

bray.dist.bac.lower<-get_lower_tri(bray.dist.bac)
bray.dist.bac.melt <- melt(as.matrix(bray.dist.bac.lower), na.rm = T)
head(bray.dist.bac.melt)
names(bray.dist.bac.melt)[3] <- "Dissimilarity"
bray.dist.bac.melt <- subset(bray.dist.bac.melt, Dissimilarity != 0)

map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)


b.meta.seed <- subset(map, Compartment == "Seed")

##2017
CC1.G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("CC1G120"))]
CC2.G.140 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("CC2G140"))]
UF1.G.140 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("UF1G140"))]
UF1.G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("UF1G120"))]
G.76 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_76"))]
G.90 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_90"))]
G.106 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_106"))]
G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_120"))]
G.141 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance
## Comparison with G.141
bray.dist.bac.CC1.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.141 & Var2 %in% CC1.G.120)
bray.dist.bac.CC1.G.2 <- subset(bray.dist.bac.melt, Var1 %in% CC1.G.120 & Var2 %in% G.141)
bray.dist.bac.CC1.G <- rbind(bray.dist.bac.CC1.G.1, bray.dist.bac.CC1.G.2)
length(bray.dist.bac.CC1.G$Dissimilarity)
bray.dist.bac.CC1.G$Comparison <- "CC1.120"
mean(bray.dist.bac.CC1.G$Dissimilarity) #0.480772

bray.dist.bac.CC2.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.141 & Var2 %in% CC2.G.140)
bray.dist.bac.CC2.G.2 <- subset(bray.dist.bac.melt, Var1 %in% CC2.G.140 & Var2 %in% G.141)
bray.dist.bac.CC2.G <- rbind(bray.dist.bac.CC2.G.1, bray.dist.bac.CC2.G.2)
length(bray.dist.bac.CC2.G$Dissimilarity)
bray.dist.bac.CC2.G$Comparison <- "CC2.140"
mean(bray.dist.bac.CC2.G$Dissimilarity) #0.5138156

bray.dist.bac.UF1.G.1 <- subset(bray.dist.bac.melt, Var1 %in% G.141 & Var2 %in% UF1.G.140 )
bray.dist.bac.UF1.G.2 <- subset(bray.dist.bac.melt, Var1 %in% UF1.G.140  & Var2 %in% G.141)
bray.dist.bac.UF1.G <- rbind(bray.dist.bac.UF1.G.1, bray.dist.bac.UF1.G.2)
length(bray.dist.bac.UF1.G$Dissimilarity)
bray.dist.bac.UF1.G$Comparison <- "UF1.140"
mean(bray.dist.bac.UF1.G$Dissimilarity) # 0.4677836


## Comparison within 2017
bray.dist.bac.CC12.G.1 <- subset(bray.dist.bac.melt, Var1 %in% CC2.G.140 & Var2 %in% CC1.G.120)
bray.dist.bac.CC12.G.2 <- subset(bray.dist.bac.melt, Var1 %in% CC1.G.120 & Var2 %in% CC2.G.140)
bray.dist.bac.CC12.G <- rbind(bray.dist.bac.CC12.G.1, bray.dist.bac.CC12.G.2)
length(bray.dist.bac.CC12.G$Dissimilarity)
bray.dist.bac.CC12.G$Comparison <- "CC1_2"
mean(bray.dist.bac.CC12.G$Dissimilarity) #0.3235775

bray.dist.bac.CC2UF1.G.1 <- subset(bray.dist.bac.melt, Var1 %in% UF1.G.140 & Var2 %in% CC2.G.140)
bray.dist.bac.CC2UF1.G.2 <- subset(bray.dist.bac.melt, Var1 %in% CC2.G.140 & Var2 %in% UF1.G.140)
bray.dist.bac.CC2UF1.G <- rbind(bray.dist.bac.CC2UF1.G.1, bray.dist.bac.CC2UF1.G.2)
length(bray.dist.bac.CC2UF1.G$Dissimilarity)
bray.dist.bac.CC2UF1.G$Comparison <- "CC2_UF1"
mean(bray.dist.bac.CC2UF1.G$Dissimilarity) #0.299885

bray.dist.bac.CC1UF1.G.1 <- subset(bray.dist.bac.melt, Var1 %in% UF1.G.140 & Var2 %in% CC1.G.120)
bray.dist.bac.CC1UF1.G.2 <- subset(bray.dist.bac.melt, Var1 %in% CC1.G.120 & Var2 %in% UF1.G.140)
bray.dist.bac.CC1UF1.G <- rbind(bray.dist.bac.CC1UF1.G.1, bray.dist.bac.CC1UF1.G.2)
length(bray.dist.bac.CC1UF1.G$Dissimilarity)
bray.dist.bac.CC1UF1.G$Comparison <- "CC1_UF1"
mean(bray.dist.bac.CC1UF1.G$Dissimilarity) #0.3468656

bray.dist.bac.1718 <- rbind(bray.dist.bac.CC1.G,bray.dist.bac.CC2.G,bray.dist.bac.UF1.G,
                            bray.dist.bac.CC12.G,bray.dist.bac.CC2UF1.G,bray.dist.bac.CC1UF1.G)


##Fungi
f.otu.lognorm <- otu_table(fun.clean.log)
f.otu.lognorm.persist <- subset(f.otu.lognorm, rownames(f.otu.lognorm) %in% persistent.fun)
f.otu.lognorm.persist.t <- t(f.otu.lognorm.persist)
f.otu.lognorm.persist.t <- subset(f.otu.lognorm.persist.t,rownames(f.otu.lognorm.persist.t) %in% seed.sample.list.4)

bray.dist.fun<-vegdist(f.otu.lognorm.persist.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun <-as.matrix(bray.dist.fun)

bray.dist.fun.lower<-get_lower_tri(bray.dist.fun)
bray.dist.fun.melt <- melt(as.matrix(bray.dist.fun.lower), na.rm = T)
head(bray.dist.fun.melt)
names(bray.dist.fun.melt)[3] <- "Dissimilarity"
bray.dist.fun.melt <- subset(bray.dist.fun.melt, Dissimilarity != 0)

f.meta<-sample_data(fun.clean.ss.f)
f.meta <- data.frame(f.meta)
f.meta.seed <- subset(f.meta, Compartment %in% c("Seed","Grain"))


##2017

CC1.G.120 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("CC1G120"))]
CC2.G.140 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("CC2G140"))]
UF1.G.140 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("UF1G140"))]
G.141 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance
bray.dist.fun.CC1.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.141 & Var2 %in% CC1.G.120)
bray.dist.fun.CC1.G.2 <- subset(bray.dist.fun.melt, Var1 %in% CC1.G.120 & Var2 %in% G.141)
bray.dist.fun.CC1.G <- rbind(bray.dist.fun.CC1.G.1, bray.dist.fun.CC1.G.2)
length(bray.dist.fun.CC1.G$Dissimilarity)
bray.dist.fun.CC1.G$Comparison <- "CC1.120"
mean(bray.dist.fun.CC1.G$Dissimilarity) #0.4324067

bray.dist.fun.CC2.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.141 & Var2 %in% CC2.G.140)
bray.dist.fun.CC2.G.2 <- subset(bray.dist.fun.melt, Var1 %in% CC2.G.140 & Var2 %in% G.141)
bray.dist.fun.CC2.G <- rbind(bray.dist.fun.CC2.G.1, bray.dist.fun.CC2.G.2)
length(bray.dist.fun.CC2.G$Dissimilarity)
bray.dist.fun.CC2.G$Comparison <- "CC2.140"
mean(bray.dist.fun.CC2.G$Dissimilarity) #0.3833042

bray.dist.fun.UF1.G.1 <- subset(bray.dist.fun.melt, Var1 %in% G.141 & Var2 %in% UF1.G.140 )
bray.dist.fun.UF1.G.2 <- subset(bray.dist.fun.melt, Var1 %in% UF1.G.140  & Var2 %in% G.141)
bray.dist.fun.UF1.G <- rbind(bray.dist.fun.UF1.G.1, bray.dist.fun.UF1.G.2)
length(bray.dist.fun.UF1.G$Dissimilarity)
bray.dist.fun.UF1.G$Comparison <- "UF1.140"
mean(bray.dist.fun.UF1.G$Dissimilarity) # 0.4211388


## Comparison within 2017
bray.dist.fun.CC12.G.1 <- subset(bray.dist.fun.melt, Var1 %in% CC2.G.140 & Var2 %in% CC1.G.120)
bray.dist.fun.CC12.G.2 <- subset(bray.dist.fun.melt, Var1 %in% CC1.G.120 & Var2 %in% CC2.G.140)
bray.dist.fun.CC12.G <- rbind(bray.dist.fun.CC12.G.1, bray.dist.fun.CC12.G.2)
length(bray.dist.fun.CC12.G$Dissimilarity)
bray.dist.fun.CC12.G$Comparison <- "CC1_2"
mean(bray.dist.fun.CC12.G$Dissimilarity) #0.3634917

bray.dist.fun.CC2UF1.G.1 <- subset(bray.dist.fun.melt, Var1 %in% UF1.G.140 & Var2 %in% CC2.G.140)
bray.dist.fun.CC2UF1.G.2 <- subset(bray.dist.fun.melt, Var1 %in% CC2.G.140 & Var2 %in% UF1.G.140)
bray.dist.fun.CC2UF1.G <- rbind(bray.dist.fun.CC2UF1.G.1, bray.dist.fun.CC2UF1.G.2)
length(bray.dist.fun.CC2UF1.G$Dissimilarity)
bray.dist.fun.CC2UF1.G$Comparison <- "CC2_UF1"
mean(bray.dist.fun.CC2UF1.G$Dissimilarity) #0.3727367

bray.dist.fun.CC1UF1.G.1 <- subset(bray.dist.fun.melt, Var1 %in% UF1.G.140 & Var2 %in% CC1.G.120)
bray.dist.fun.CC1UF1.G.2 <- subset(bray.dist.fun.melt, Var1 %in% CC1.G.120 & Var2 %in% UF1.G.140)
bray.dist.fun.CC1UF1.G <- rbind(bray.dist.fun.CC1UF1.G.1, bray.dist.fun.CC1UF1.G.2)
length(bray.dist.fun.CC1UF1.G$Dissimilarity)
bray.dist.fun.CC1UF1.G$Comparison <- "CC1_UF1"
mean(bray.dist.fun.CC1UF1.G$Dissimilarity) #0.4356085


bray.dist.fun.1718 <- rbind(bray.dist.fun.CC1.G,bray.dist.fun.CC2.G,bray.dist.fun.UF1.G,
                            bray.dist.fun.CC12.G,bray.dist.fun.CC2UF1.G,bray.dist.fun.CC1UF1.G)


write.csv(bray.dist.bac.1718,"Bray-curtis dissimilarity_bacterial_persistent_OTUs_20172018.csv")
write.csv(bray.dist.fun.1718,"Bray-curtis dissimilarity_fungal_persistent_OTUs_20172018.csv")


#### Kruskal-Wallis test
max.bray <- aggregate(bray.dist.fun.1718$Dissimilarity, by = list(bray.dist.fun.1718$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.fun.1718$Comparison <- as.factor(bray.dist.fun.1718$Comparison)
bray.dist.fun.1718$Comparison <- factor(bray.dist.fun.1718$Comparison, levels = c("CC1_2","CC1_UF1","CC2_UF1","CC1.120","CC2.140","UF1.140"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.1718)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.1718,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn$Group <- c("CC1.120", "CC1_2","CC1_UF1","CC2.140","CC2_UF1","UF1.140")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.fun.1718, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p



max.bray <- aggregate(bray.dist.bac.1718$Dissimilarity, by = list(bray.dist.bac.1718$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.bac.1718$Comparison <- as.factor(bray.dist.bac.1718$Comparison)
bray.dist.bac.1718$Comparison <- factor(bray.dist.bac.1718$Comparison, levels = c("CC1_2","CC1_UF1","CC2_UF1","CC1.120","CC2.140","UF1.140"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.1718)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.1718,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn$Group <- c("CC1.120", "CC1_2","CC1_UF1","CC2.140","CC2_UF1","UF1.140")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.bac.1718, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p


#### Core
bac.core<-read.xlsx("Bacterial core profile.xlsx",1)
fun.core<-read.xlsx("Fungal core profile.xlsx",1)

#### Bray-Curtis dissimilarity - Bacteria
b.otu.lognorm <- otu_table(bac.clean.log)
b.otu.lognorm.core <- subset(b.otu.lognorm, rownames(b.otu.lognorm) %in% bac.core$OTU)
b.otu.lognorm.core.t <- t(b.otu.lognorm.core)
b.otu.lognorm.core.t <- subset(b.otu.lognorm.core.t,rownames(b.otu.lognorm.core.t) %in% seed.sample.list.3)

bray.dist.bac.core<-vegdist(b.otu.lognorm.core.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.bac.core <-as.matrix(bray.dist.bac.core)

bray.dist.bac.core.lower<-get_lower_tri(bray.dist.bac.core)
bray.dist.bac.core.melt <- melt(as.matrix(bray.dist.bac.core.lower), na.rm = T)
head(bray.dist.bac.core.melt)
names(bray.dist.bac.core.melt)[3] <- "Dissimilarity"
bray.dist.bac.core.melt <- subset(bray.dist.bac.core.melt, Dissimilarity != 0)

map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)


b.meta.seed <- subset(map, Compartment == "Seed")

##2017
CC1.G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("CC1G120"))]
CC2.G.140 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("CC2G140"))]
UF1.G.140 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("UF1G140"))]
UF1.G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("UF1G120"))]
G.76 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_76"))]
G.90 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_90"))]
G.106 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_106"))]
G.120 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_120"))]
G.141 <- b.meta.seed$SampleID[which(b.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance
## Comparison with G.141
bray.dist.bac.CC1.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% G.141 & Var2 %in% CC1.G.120)
bray.dist.bac.CC1.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% CC1.G.120 & Var2 %in% G.141)
bray.dist.bac.CC1.G <- rbind(bray.dist.bac.CC1.G.1, bray.dist.bac.CC1.G.2)
length(bray.dist.bac.CC1.G$Dissimilarity)
bray.dist.bac.CC1.G$Comparison <- "CC1.120"
mean(bray.dist.bac.CC1.G$Dissimilarity) #0.480772 #0.4784973(core)

bray.dist.bac.CC2.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% G.141 & Var2 %in% CC2.G.140)
bray.dist.bac.CC2.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% CC2.G.140 & Var2 %in% G.141)
bray.dist.bac.CC2.G <- rbind(bray.dist.bac.CC2.G.1, bray.dist.bac.CC2.G.2)
length(bray.dist.bac.CC2.G$Dissimilarity)
bray.dist.bac.CC2.G$Comparison <- "CC2.140"
mean(bray.dist.bac.CC2.G$Dissimilarity) #0.5138156 #0.5116184(core)

bray.dist.bac.UF1.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% G.141 & Var2 %in% UF1.G.140 )
bray.dist.bac.UF1.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% UF1.G.140  & Var2 %in% G.141)
bray.dist.bac.UF1.G <- rbind(bray.dist.bac.UF1.G.1, bray.dist.bac.UF1.G.2)
length(bray.dist.bac.UF1.G$Dissimilarity)
bray.dist.bac.UF1.G$Comparison <- "UF1.140"
mean(bray.dist.bac.UF1.G$Dissimilarity) # 0.4677836 #0.4624183(core)


## Comparison within 2017
bray.dist.bac.CC12.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% CC2.G.140 & Var2 %in% CC1.G.120)
bray.dist.bac.CC12.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% CC1.G.120 & Var2 %in% CC2.G.140)
bray.dist.bac.CC12.G <- rbind(bray.dist.bac.CC12.G.1, bray.dist.bac.CC12.G.2)
length(bray.dist.bac.CC12.G$Dissimilarity)
bray.dist.bac.CC12.G$Comparison <- "CC1_2"
mean(bray.dist.bac.CC12.G$Dissimilarity) #0.3235775 #0.3235775(core)

bray.dist.bac.CC2UF1.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% UF1.G.140 & Var2 %in% CC2.G.140)
bray.dist.bac.CC2UF1.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% CC2.G.140 & Var2 %in% UF1.G.140)
bray.dist.bac.CC2UF1.G <- rbind(bray.dist.bac.CC2UF1.G.1, bray.dist.bac.CC2UF1.G.2)
length(bray.dist.bac.CC2UF1.G$Dissimilarity)
bray.dist.bac.CC2UF1.G$Comparison <- "CC2_UF1"
mean(bray.dist.bac.CC2UF1.G$Dissimilarity) #0.299885 # 0.2954689(core)

bray.dist.bac.CC1UF1.G.1 <- subset(bray.dist.bac.core.melt, Var1 %in% UF1.G.140 & Var2 %in% CC1.G.120)
bray.dist.bac.CC1UF1.G.2 <- subset(bray.dist.bac.core.melt, Var1 %in% CC1.G.120 & Var2 %in% UF1.G.140)
bray.dist.bac.CC1UF1.G <- rbind(bray.dist.bac.CC1UF1.G.1, bray.dist.bac.CC1UF1.G.2)
length(bray.dist.bac.CC1UF1.G$Dissimilarity)
bray.dist.bac.CC1UF1.G$Comparison <- "CC1_UF1"
mean(bray.dist.bac.CC1UF1.G$Dissimilarity) #0.3468656 #0.3425778(core)

bray.dist.bac.core.1718 <- rbind(bray.dist.bac.CC1.G,bray.dist.bac.CC2.G,bray.dist.bac.UF1.G,
                            bray.dist.bac.CC12.G,bray.dist.bac.CC2UF1.G,bray.dist.bac.CC1UF1.G)


##Fungi
f.otu.lognorm <- otu_table(fun.clean.log)
f.otu.lognorm.core <- subset(f.otu.lognorm, rownames(f.otu.lognorm) %in% fun.core$OTU)
f.otu.lognorm.core.t <- t(f.otu.lognorm.core)
f.otu.lognorm.core.t <- subset(f.otu.lognorm.core.t,rownames(f.otu.lognorm.core.t) %in% seed.sample.list.4)

bray.dist.fun.core<-vegdist(f.otu.lognorm.core.t, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

bray.dist.fun.core <-as.matrix(bray.dist.fun.core)

bray.dist.fun.core.lower<-get_lower_tri(bray.dist.fun.core)
bray.dist.fun.core.melt <- melt(as.matrix(bray.dist.fun.core.lower), na.rm = T)
head(bray.dist.fun.core.melt)
names(bray.dist.fun.core.melt)[3] <- "Dissimilarity"
bray.dist.fun.core.melt <- subset(bray.dist.fun.core.melt, Dissimilarity != 0)

f.meta<-sample_data(fun.clean.ss.f)
f.meta <- data.frame(f.meta)
f.meta.seed <- subset(f.meta, Compartment %in% c("Seed","Grain"))


##2017

CC1.G.120 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("CC1G120"))]
CC2.G.140 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("CC2G140"))]
UF1.G.140 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("UF1G140"))]
G.141 <- f.meta.seed$SampleID[which(f.meta.seed$Replication %in% c("G_141"))]


### Bray-Curtis distance
bray.dist.fun.CC1.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% G.141 & Var2 %in% CC1.G.120)
bray.dist.fun.CC1.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% CC1.G.120 & Var2 %in% G.141)
bray.dist.fun.CC1.G <- rbind(bray.dist.fun.CC1.G.1, bray.dist.fun.CC1.G.2)
length(bray.dist.fun.CC1.G$Dissimilarity)
bray.dist.fun.CC1.G$Comparison <- "CC1.120"
mean(bray.dist.fun.CC1.G$Dissimilarity) #0.4324067 #0.4439423

bray.dist.fun.CC2.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% G.141 & Var2 %in% CC2.G.140)
bray.dist.fun.CC2.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% CC2.G.140 & Var2 %in% G.141)
bray.dist.fun.CC2.G <- rbind(bray.dist.fun.CC2.G.1, bray.dist.fun.CC2.G.2)
length(bray.dist.fun.CC2.G$Dissimilarity)
bray.dist.fun.CC2.G$Comparison <- "CC2.140"
mean(bray.dist.fun.CC2.G$Dissimilarity) #0.3833042 #0.3874962

bray.dist.fun.UF1.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% G.141 & Var2 %in% UF1.G.140 )
bray.dist.fun.UF1.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% UF1.G.140  & Var2 %in% G.141)
bray.dist.fun.UF1.G <- rbind(bray.dist.fun.UF1.G.1, bray.dist.fun.UF1.G.2)
length(bray.dist.fun.UF1.G$Dissimilarity)
bray.dist.fun.UF1.G$Comparison <- "UF1.140"
mean(bray.dist.fun.UF1.G$Dissimilarity) # 0.4211388 #0.4062176


## Comparison within 2017
bray.dist.fun.CC12.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% CC2.G.140 & Var2 %in% CC1.G.120)
bray.dist.fun.CC12.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% CC1.G.120 & Var2 %in% CC2.G.140)
bray.dist.fun.CC12.G <- rbind(bray.dist.fun.CC12.G.1, bray.dist.fun.CC12.G.2)
length(bray.dist.fun.CC12.G$Dissimilarity)
bray.dist.fun.CC12.G$Comparison <- "CC1_2"
mean(bray.dist.fun.CC12.G$Dissimilarity) #0.3634917 #0.3868134

bray.dist.fun.CC2UF1.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% UF1.G.140 & Var2 %in% CC2.G.140)
bray.dist.fun.CC2UF1.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% CC2.G.140 & Var2 %in% UF1.G.140)
bray.dist.fun.CC2UF1.G <- rbind(bray.dist.fun.CC2UF1.G.1, bray.dist.fun.CC2UF1.G.2)
length(bray.dist.fun.CC2UF1.G$Dissimilarity)
bray.dist.fun.CC2UF1.G$Comparison <- "CC2_UF1"
mean(bray.dist.fun.CC2UF1.G$Dissimilarity) #0.3727367 #0.3634141

bray.dist.fun.CC1UF1.G.1 <- subset(bray.dist.fun.core.melt, Var1 %in% UF1.G.140 & Var2 %in% CC1.G.120)
bray.dist.fun.CC1UF1.G.2 <- subset(bray.dist.fun.core.melt, Var1 %in% CC1.G.120 & Var2 %in% UF1.G.140)
bray.dist.fun.CC1UF1.G <- rbind(bray.dist.fun.CC1UF1.G.1, bray.dist.fun.CC1UF1.G.2)
length(bray.dist.fun.CC1UF1.G$Dissimilarity)
bray.dist.fun.CC1UF1.G$Comparison <- "CC1_UF1"
mean(bray.dist.fun.CC1UF1.G$Dissimilarity) #0.4356085 #0.4435908


bray.dist.fun.core.1718 <- rbind(bray.dist.fun.CC1.G,bray.dist.fun.CC2.G,bray.dist.fun.UF1.G,
                            bray.dist.fun.CC12.G,bray.dist.fun.CC2UF1.G,bray.dist.fun.CC1UF1.G)


write.csv(bray.dist.bac.1718,"Bray-curtis dissimilarity_bacterial_persistent_OTUs_20172018.csv")
write.csv(bray.dist.fun.1718,"Bray-curtis dissimilarity_fungal_persistent_OTUs_20172018.csv")


#### Kruskal-Wallis test
max.bray <- aggregate(bray.dist.fun.core.1718$Dissimilarity, by = list(bray.dist.fun.1718$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.fun.core.1718$Comparison <- as.factor(bray.dist.fun.core.1718$Comparison)
bray.dist.fun.core.1718$Comparison <- factor(bray.dist.fun.core.1718$Comparison, levels = c("CC1_2","CC1_UF1","CC2_UF1","CC1.120","CC2.140","UF1.140"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.core.1718)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.core.1718,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn$Group <- c("CC1.120", "CC1_2","CC1_UF1","CC2.140","CC2_UF1","UF1.140")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.fun.core.1718, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p



max.bray <- aggregate(bray.dist.bac.core.1718$Dissimilarity, by = list(bray.dist.bac.1718$Comparison), max)
colnames(max.bray) <- c("Group", "Maxbray")

##Kruskal-Wallis test
bray.dist.bac.core.1718$Comparison <- as.factor(bray.dist.bac.core.1718$Comparison)
bray.dist.bac.core.1718$Comparison <- factor(bray.dist.bac.core.1718$Comparison, levels = c("CC1_2","CC1_UF1","CC2_UF1","CC1.120","CC2.140","UF1.140"))

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.core.1718)
kw$p.value
kw$p.value<- round(kw$p.value, 10)

#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.core.1718,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn$Group <- c("CC1.120", "CC1_2","CC1_UF1","CC2.140","CC2_UF1","UF1.140")

hsd1 <- merge(dunn,max.bray, by = 'Group')
p <- ggplot(data = bray.dist.bac.core.1718, aes(x=Comparison, y=Dissimilarity)) + geom_boxplot(aes(fill = Comparison), width = 0.8) +
  theme_bw() + theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=Maxbray, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Bray-Curtis dissimilarity\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))

p


### Abundance distribution of persistent OTUs
### Relative abundance of persistent, lost, and gained OTUs in developing seeds
##Bacteria
##2017 + 2018
bac.clean.ss.G <- subset_samples(bac.clean.ss.f, Compartment == "Seed")
bac.clean.ss.G <- phyloseq::filter_taxa(bac.clean.ss.G, function(x) sum(x) != 0, TRUE)


bac.G.comp.rel <- transform(bac.clean.ss.G, 'compositional')
melt.bac.G.comp.rel <-psmelt(bac.G.comp.rel)

melt.bac.G.comp.rel$Genus <- as.character(melt.bac.G.comp.rel$Genus)
melt.bac.G.comp.rel$Genus[is.na(melt.bac.G.comp.rel$Genus)] <- "Unidentified"
melt.otu.bac.G.comp.rel.persist <- subset(melt.bac.G.comp.rel, OTU %in% persistent.bac)

otu.bac.G.comp.rel.sum.tab <- melt.otu.bac.G.comp.rel.persist %>% group_by(Days, Year, OTU, Genus,Location) %>% summarise(Mean = mean(Abundance)) %>% group_by(Genus, Days, Year, Location) %>% summarise(TotalAbun = sum(Mean))


otu.bac.G.comp.rel.sum.tab$Days <- factor(otu.bac.G.comp.rel.sum.tab$Days, levels = c("0","76","80","90","106","120","140","141"))

##Ordering genus
ord <- melt.otu.bac.G.comp.rel.persist %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Unidentified"))]
vec.uniden.Low <- c("Unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)

otu.bac.G.comp.rel.sum.tab$Genus <- factor(otu.bac.G.comp.rel.sum.tab$Genus, levels =vec.reorder)

### plotting
p <- ggplot(otu.bac.G.comp.rel.sum.tab, aes(x=Days, y = TotalAbun, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Year+Location, scales = 'free', nrow =1)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 0.5)+
  xlab('')+
  ylab("Relative abundance\n") +coord_cartesian(ylim = c(0, 1))+
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p


## Fungi
##2017 + 2018
fun.clean.ss.G <- subset_samples(fun.clean.ss.f, Compartment %in% c("Seed","Grain"))
fun.clean.ss.G <- phyloseq::filter_taxa(fun.clean.ss.G, function(x) sum(x) != 0, TRUE)


fun.G.comp.rel <- transform(fun.clean.ss.G, 'compositional')

melt.fun.G.comp.rel <-psmelt(fun.G.comp.rel)

melt.fun.G.comp.rel$Genus <- as.character(melt.fun.G.comp.rel$Genus)
melt.fun.G.comp.rel$Genus[is.na(melt.fun.G.comp.rel$Genus)] <- "unidentified"
melt.otu.fun.G.comp.rel.persist <- subset(melt.fun.G.comp.rel, OTU %in% persistent.fun)

otu.fun.G.comp.rel.sum.tab <- melt.otu.fun.G.comp.rel.persist %>% group_by(Days, Year, OTU, Genus,Location) %>% summarise(Mean = mean(Abundance)) %>% group_by(Genus, Days, Year, Location) %>% summarise(TotalAbun = sum(Mean))


otu.fun.G.comp.rel.sum.tab$Days <- factor(otu.fun.G.comp.rel.sum.tab$Days, levels = c("0","76","80","90","106","120","140","141"))

##Ordering genus
ord <- melt.otu.fun.G.comp.rel.persist %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("unidentified"))]
vec.uniden.Low <- c("unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)

otu.fun.G.comp.rel.sum.tab$Genus <- factor(otu.fun.G.comp.rel.sum.tab$Genus, levels =vec.reorder)

### plotting
p <- ggplot(otu.fun.G.comp.rel.sum.tab, aes(x=Days, y = TotalAbun, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Year+Location, scales = 'free', nrow=1)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 0.5)+
  xlab('')+
  ylab("Relative abundance\n") +coord_cartesian(ylim = c(0, 1))+
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p


write.csv(otu.fun.G.comp.rel.sum.tab,"Raw data table for fungal genus_persistent OTU.csv")
write.csv(otu.bac.G.comp.rel.sum.tab,"Raw data table for bacterial genus_persistent OTU.csv")

write.csv(melt.otu.fun.G.comp.rel.persist,"Raw data table for RA of fungal persistent OTU.csv")
write.csv(melt.otu.bac.G.comp.rel.persist,"Raw data table for RA of bacterial persistent OTU.csv")
