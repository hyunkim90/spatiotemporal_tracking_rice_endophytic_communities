##### Jaccard distance with G_0
f.otu.lognorm <- otu_table(fun.clean.log.18)
jaccard.dist.fun<-vegdist(t(f.otu.lognorm), method="jaccard", binary=F, diag=FALSE, upper=FALSE, na.rm = T)
class(jaccard.dist.fun)

jaccard.dist.fun <-as.matrix(jaccard.dist.fun)

jaccard.dist.fun.lower<-get_lower_tri(jaccard.dist.fun)
jaccard.dist.fun.melt <- melt(as.matrix(jaccard.dist.fun.lower), na.rm = T)
head(jaccard.dist.fun.melt)
names(jaccard.dist.fun.melt)[3] <- "Dissimilarity"
jaccard.dist.fun.melt <- subset(jaccard.dist.fun.melt, Dissimilarity != 0)



#### Jaccard dissimiarity - Bacteria
b.otu.lognorm <- otu_table(bac.clean.log.18)
jaccard.dist.bac<-vegdist(t(b.otu.lognorm), method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)

jaccard.dist.bac <-as.matrix(jaccard.dist.bac)

jaccard.dist.bac.lower<-get_lower_tri(jaccard.dist.bac)
jaccard.dist.bac.melt <- melt(as.matrix(jaccard.dist.bac.lower), na.rm = T)
head(jaccard.dist.bac.melt)
names(jaccard.dist.bac.melt)[3] <- "Dissimilarity"
jaccard.dist.bac.melt <- subset(jaccard.dist.bac.melt, Dissimilarity != 0)

map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)


#### Distance with G_0
##Fungi

G.0<-f.meta.18$SampleID[which(f.meta.18$Days == 0 & f.meta.18$Microhabitat %in% c("G"))]

## 48
L1.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("L1"))]
L2.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("L2"))]

jaccard.dist.fun.48.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.48 & Var2 %in% G.0)
jaccard.dist.fun.48.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.48)
jaccard.dist.fun.48.L1 <- rbind(jaccard.dist.fun.48.L1.1, jaccard.dist.fun.48.L1.2)
length(jaccard.dist.fun.48.L1$Dissimilarity)
jaccard.dist.fun.48.L1$Comparison <- "L1"
mean(jaccard.dist.fun.48.L1$Dissimilarity)

jaccard.dist.fun.48.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.48 & Var2 %in% G.0)
jaccard.dist.fun.48.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.48)
jaccard.dist.fun.48.L2 <- rbind(jaccard.dist.fun.48.L2.1, jaccard.dist.fun.48.L2.2)
length(jaccard.dist.fun.48.L2$Dissimilarity)
jaccard.dist.fun.48.L2$Comparison <- "L2"
mean(jaccard.dist.fun.48.L2$Dissimilarity)


jaccard.dist.fun.48.L <- rbind(jaccard.dist.fun.48.L1,jaccard.dist.fun.48.L2)
jaccard.dist.fun.48.L$Days <- 48

##62
L1.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L1"))]
L2.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L2"))]
L3.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L3"))]

jaccard.dist.fun.62.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.62 & Var2 %in% G.0)
jaccard.dist.fun.62.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.62)
jaccard.dist.fun.62.L1 <- rbind(jaccard.dist.fun.62.L1.1, jaccard.dist.fun.62.L1.2)
length(jaccard.dist.fun.62.L1$Dissimilarity)
jaccard.dist.fun.62.L1$Comparison <- "L1"
mean(jaccard.dist.fun.62.L1$Dissimilarity)

jaccard.dist.fun.62.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.62 & Var2 %in% G.0)
jaccard.dist.fun.62.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.62)
jaccard.dist.fun.62.L2 <- rbind(jaccard.dist.fun.62.L2.1, jaccard.dist.fun.62.L2.2)
length(jaccard.dist.fun.62.L2$Dissimilarity)
jaccard.dist.fun.62.L2$Comparison <- "L2"
mean(jaccard.dist.fun.62.L2$Dissimilarity)

jaccard.dist.fun.62.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.62 & Var2 %in% G.0)
jaccard.dist.fun.62.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.62)
jaccard.dist.fun.62.L3 <- rbind(jaccard.dist.fun.62.L3.1, jaccard.dist.fun.62.L3.2)
length(jaccard.dist.fun.62.L3$Dissimilarity)
jaccard.dist.fun.62.L3$Comparison <- "L3"
mean(jaccard.dist.fun.62.L3$Dissimilarity)

jaccard.dist.fun.62.L <- rbind(jaccard.dist.fun.62.L1,jaccard.dist.fun.62.L2,jaccard.dist.fun.62.L3)
jaccard.dist.fun.62.L$Days <- 62


#76
L1.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L1"))]
L2.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L2"))]
L3.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L3"))]
FL.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.76.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.76 & Var2 %in% G.0)
jaccard.dist.fun.76.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.76)
jaccard.dist.fun.76.L1 <- rbind(jaccard.dist.fun.76.L1.1, jaccard.dist.fun.76.L1.2)
length(jaccard.dist.fun.76.L1$Dissimilarity)
jaccard.dist.fun.76.L1$Comparison <- "L1"
mean(jaccard.dist.fun.76.L1$Dissimilarity)

jaccard.dist.fun.76.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.76 & Var2 %in% G.0)
jaccard.dist.fun.76.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.76)
jaccard.dist.fun.76.L2 <- rbind(jaccard.dist.fun.76.L2.1, jaccard.dist.fun.76.L2.2)
length(jaccard.dist.fun.76.L2$Dissimilarity)
jaccard.dist.fun.76.L2$Comparison <- "L2"
mean(jaccard.dist.fun.76.L2$Dissimilarity)

jaccard.dist.fun.76.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.76 & Var2 %in% G.0)
jaccard.dist.fun.76.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.76)
jaccard.dist.fun.76.L3 <- rbind(jaccard.dist.fun.76.L3.1, jaccard.dist.fun.76.L3.2)
length(jaccard.dist.fun.76.L3$Dissimilarity)
jaccard.dist.fun.76.L3$Comparison <- "L3"
mean(jaccard.dist.fun.76.L3$Dissimilarity)

jaccard.dist.fun.76.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.76 & Var2 %in% G.0)
jaccard.dist.fun.76.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% FL.76)
jaccard.dist.fun.76.FL <- rbind(jaccard.dist.fun.76.FL.1, jaccard.dist.fun.76.FL.2)
length(jaccard.dist.fun.76.FL$Dissimilarity)
jaccard.dist.fun.76.FL$Comparison <- "FL"
mean(jaccard.dist.fun.76.FL$Dissimilarity)

jaccard.dist.fun.76.L <- rbind(jaccard.dist.fun.76.L1,jaccard.dist.fun.76.L2,jaccard.dist.fun.76.L3,
                               jaccard.dist.fun.76.FL)
jaccard.dist.fun.76.L$Days <- 76

###90 days
L1.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L1"))]
L2.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L2"))]
L3.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L3"))]
FL.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.90.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.90 & Var2 %in% G.0)
jaccard.dist.fun.90.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.90)
jaccard.dist.fun.90.L1 <- rbind(jaccard.dist.fun.90.L1.1, jaccard.dist.fun.90.L1.2)
length(jaccard.dist.fun.90.L1$Dissimilarity)
jaccard.dist.fun.90.L1$Comparison <- "L1"
mean(jaccard.dist.fun.90.L1$Dissimilarity)

jaccard.dist.fun.90.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.90 & Var2 %in% G.0)
jaccard.dist.fun.90.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.90)
jaccard.dist.fun.90.L2 <- rbind(jaccard.dist.fun.90.L2.1, jaccard.dist.fun.90.L2.2)
length(jaccard.dist.fun.90.L2$Dissimilarity)
jaccard.dist.fun.90.L2$Comparison <- "L2"
mean(jaccard.dist.fun.90.L2$Dissimilarity)

jaccard.dist.fun.90.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.90 & Var2 %in% G.0)
jaccard.dist.fun.90.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.90)
jaccard.dist.fun.90.L3 <- rbind(jaccard.dist.fun.90.L3.1, jaccard.dist.fun.90.L3.2)
length(jaccard.dist.fun.90.L3$Dissimilarity)
jaccard.dist.fun.90.L3$Comparison <- "L3"
mean(jaccard.dist.fun.90.L3$Dissimilarity)

jaccard.dist.fun.90.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.90 & Var2 %in% G.0)
jaccard.dist.fun.90.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% FL.90)
jaccard.dist.fun.90.FL <- rbind(jaccard.dist.fun.90.FL.1, jaccard.dist.fun.90.FL.2)
length(jaccard.dist.fun.90.FL$Dissimilarity)
jaccard.dist.fun.90.FL$Comparison <- "FL"
mean(jaccard.dist.fun.90.FL$Dissimilarity)


jaccard.dist.fun.90.L <- rbind(jaccard.dist.fun.90.L1,jaccard.dist.fun.90.L2,jaccard.dist.fun.90.L3,
                               jaccard.dist.fun.90.FL)

jaccard.dist.fun.90.L$Days <- 90

###106 days
L1.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L1"))]
L2.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L2"))]
L3.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L3"))]
FL.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("FL"))]


jaccard.dist.fun.106.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.106 & Var2 %in% G.0)
jaccard.dist.fun.106.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.106)
jaccard.dist.fun.106.L1 <- rbind(jaccard.dist.fun.106.L1.1, jaccard.dist.fun.106.L1.2)
length(jaccard.dist.fun.106.L1$Dissimilarity)
jaccard.dist.fun.106.L1$Comparison <- "L1"
mean(jaccard.dist.fun.106.L1$Dissimilarity)

jaccard.dist.fun.106.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.106 & Var2 %in% G.0)
jaccard.dist.fun.106.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.106)
jaccard.dist.fun.106.L2 <- rbind(jaccard.dist.fun.106.L2.1, jaccard.dist.fun.106.L2.2)
length(jaccard.dist.fun.106.L2$Dissimilarity)
jaccard.dist.fun.106.L2$Comparison <- "L2"
mean(jaccard.dist.fun.106.L2$Dissimilarity)

jaccard.dist.fun.106.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.106 & Var2 %in% G.0)
jaccard.dist.fun.106.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.106)
jaccard.dist.fun.106.L3 <- rbind(jaccard.dist.fun.106.L3.1, jaccard.dist.fun.106.L3.2)
length(jaccard.dist.fun.106.L3$Dissimilarity)
jaccard.dist.fun.106.L3$Comparison <- "L3"
mean(jaccard.dist.fun.106.L3$Dissimilarity)

jaccard.dist.fun.106.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.106 & Var2 %in% G.0)
jaccard.dist.fun.106.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% FL.106)
jaccard.dist.fun.106.FL <- rbind(jaccard.dist.fun.106.FL.1, jaccard.dist.fun.106.FL.2)
length(jaccard.dist.fun.106.FL$Dissimilarity)
jaccard.dist.fun.106.FL$Comparison <- "FL"
mean(jaccard.dist.fun.106.FL$Dissimilarity)

jaccard.dist.fun.106.L <- rbind(jaccard.dist.fun.106.L1,jaccard.dist.fun.106.L2,jaccard.dist.fun.106.L3,
                                jaccard.dist.fun.106.FL)
jaccard.dist.fun.106.L$Days <- 106
###120 days
L1.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L1"))]
L2.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L2"))]
L3.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L3"))]
FL.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.120.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.120 & Var2 %in% G.0)
jaccard.dist.fun.120.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.120)
jaccard.dist.fun.120.L1 <- rbind(jaccard.dist.fun.120.L1.1, jaccard.dist.fun.120.L1.2)
length(jaccard.dist.fun.120.L1$Dissimilarity)
jaccard.dist.fun.120.L1$Comparison <- "L1"
mean(jaccard.dist.fun.120.L1$Dissimilarity)

jaccard.dist.fun.120.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.120 & Var2 %in% G.0)
jaccard.dist.fun.120.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.120)
jaccard.dist.fun.120.L2 <- rbind(jaccard.dist.fun.120.L2.1, jaccard.dist.fun.120.L2.2)
length(jaccard.dist.fun.120.L2$Dissimilarity)
jaccard.dist.fun.120.L2$Comparison <- "L2"
mean(jaccard.dist.fun.120.L2$Dissimilarity)

jaccard.dist.fun.120.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.120 & Var2 %in% G.0)
jaccard.dist.fun.120.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.120)
jaccard.dist.fun.120.L3 <- rbind(jaccard.dist.fun.120.L3.1, jaccard.dist.fun.120.L3.2)
length(jaccard.dist.fun.120.L3$Dissimilarity)
jaccard.dist.fun.120.L3$Comparison <- "L3"
mean(jaccard.dist.fun.120.L3$Dissimilarity)

jaccard.dist.fun.120.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.120 & Var2 %in% G.0)
jaccard.dist.fun.120.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% FL.120)
jaccard.dist.fun.120.FL <- rbind(jaccard.dist.fun.120.FL.1, jaccard.dist.fun.120.FL.2)
length(jaccard.dist.fun.120.FL$Dissimilarity)
jaccard.dist.fun.120.FL$Comparison <- "FL"
mean(jaccard.dist.fun.120.FL$Dissimilarity)

jaccard.dist.fun.120.L <- rbind(jaccard.dist.fun.120.L1,jaccard.dist.fun.120.L2,jaccard.dist.fun.120.L3,
                                jaccard.dist.fun.120.FL)

jaccard.dist.fun.120.L$Days <- 120
###141 days
L1.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L1"))]
L2.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L2"))]
L3.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L3"))]
FL.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.141.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.141 & Var2 %in% G.0)
jaccard.dist.fun.141.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L1.141)
jaccard.dist.fun.141.L1 <- rbind(jaccard.dist.fun.141.L1.1, jaccard.dist.fun.141.L1.2)
length(jaccard.dist.fun.141.L1$Dissimilarity)
jaccard.dist.fun.141.L1$Comparison <- "L1"
mean(jaccard.dist.fun.141.L1$Dissimilarity)

jaccard.dist.fun.141.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.141 & Var2 %in% G.0)
jaccard.dist.fun.141.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L2.141)
jaccard.dist.fun.141.L2 <- rbind(jaccard.dist.fun.141.L2.1, jaccard.dist.fun.141.L2.2)
length(jaccard.dist.fun.141.L2$Dissimilarity)
jaccard.dist.fun.141.L2$Comparison <- "L2"
mean(jaccard.dist.fun.141.L2$Dissimilarity)

jaccard.dist.fun.141.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.141 & Var2 %in% G.0)
jaccard.dist.fun.141.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% L3.141)
jaccard.dist.fun.141.L3 <- rbind(jaccard.dist.fun.141.L3.1, jaccard.dist.fun.141.L3.2)
length(jaccard.dist.fun.141.L3$Dissimilarity)
jaccard.dist.fun.141.L3$Comparison <- "L3"
mean(jaccard.dist.fun.141.L3$Dissimilarity)

jaccard.dist.fun.141.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.141 & Var2 %in% G.0)
jaccard.dist.fun.141.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% FL.141)
jaccard.dist.fun.141.FL <- rbind(jaccard.dist.fun.141.FL.1, jaccard.dist.fun.141.FL.2)
length(jaccard.dist.fun.141.FL$Dissimilarity)
jaccard.dist.fun.141.FL$Comparison <- "FL"
mean(jaccard.dist.fun.141.FL$Dissimilarity)

jaccard.dist.fun.141.L <- rbind(jaccard.dist.fun.141.L1,jaccard.dist.fun.141.L2,jaccard.dist.fun.141.L3,
                                jaccard.dist.fun.141.FL)
jaccard.dist.fun.141.L$Days <- 141
##Stem

S1.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S1"))]
S2.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S2"))]
S3.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S3"))]

jaccard.dist.fun.48.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.48 & Var2 %in% G.0)
jaccard.dist.fun.48.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.48)
jaccard.dist.fun.48.S1 <- rbind(jaccard.dist.fun.48.S1.1, jaccard.dist.fun.48.S1.2)
length(jaccard.dist.fun.48.S1$Dissimilarity)
jaccard.dist.fun.48.S1$Comparison <- "S1"
mean(jaccard.dist.fun.48.S1$Dissimilarity)

jaccard.dist.fun.48.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.48 & Var2 %in% G.0)
jaccard.dist.fun.48.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.48)
jaccard.dist.fun.48.S2 <- rbind(jaccard.dist.fun.48.S2.1, jaccard.dist.fun.48.S2.2)
length(jaccard.dist.fun.48.S2$Dissimilarity)
jaccard.dist.fun.48.S2$Comparison <- "S2"
mean(jaccard.dist.fun.48.S2$Dissimilarity)

jaccard.dist.fun.48.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.48 & Var2 %in% G.0)
jaccard.dist.fun.48.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.48)
jaccard.dist.fun.48.S3 <- rbind(jaccard.dist.fun.48.S3.1, jaccard.dist.fun.48.S3.2)
length(jaccard.dist.fun.48.S3$Dissimilarity)
jaccard.dist.fun.48.S3$Comparison <- "S3"
mean(jaccard.dist.fun.48.S3$Dissimilarity)


jaccard.dist.fun.48.S <- rbind(jaccard.dist.fun.48.S1,jaccard.dist.fun.48.S2,jaccard.dist.fun.48.S3)
jaccard.dist.fun.48.S$Days <- 48


S1.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S1"))]
S2.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S2"))]
S3.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S3"))]
S4.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S4"))]

jaccard.dist.fun.62.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.62 & Var2 %in% G.0)
jaccard.dist.fun.62.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.62)
jaccard.dist.fun.62.S1 <- rbind(jaccard.dist.fun.62.S1.1, jaccard.dist.fun.62.S1.2)
length(jaccard.dist.fun.62.S1$Dissimilarity)
jaccard.dist.fun.62.S1$Comparison <- "S1"
mean(jaccard.dist.fun.62.S1$Dissimilarity)

jaccard.dist.fun.62.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.62 & Var2 %in% G.0)
jaccard.dist.fun.62.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.62)
jaccard.dist.fun.62.S2 <- rbind(jaccard.dist.fun.62.S2.1, jaccard.dist.fun.62.S2.2)
length(jaccard.dist.fun.62.S2$Dissimilarity)
jaccard.dist.fun.62.S2$Comparison <- "S2"
mean(jaccard.dist.fun.62.S2$Dissimilarity)

jaccard.dist.fun.62.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.62 & Var2 %in% G.0)
jaccard.dist.fun.62.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.62)
jaccard.dist.fun.62.S3 <- rbind(jaccard.dist.fun.62.S3.1, jaccard.dist.fun.62.S3.2)
length(jaccard.dist.fun.62.S3$Dissimilarity)
jaccard.dist.fun.62.S3$Comparison <- "S3"
mean(jaccard.dist.fun.62.S3$Dissimilarity)

jaccard.dist.fun.62.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.62 & Var2 %in% G.0)
jaccard.dist.fun.62.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.62)
jaccard.dist.fun.62.S4 <- rbind(jaccard.dist.fun.62.S4.1, jaccard.dist.fun.62.S4.2)
length(jaccard.dist.fun.62.S4$Dissimilarity)
jaccard.dist.fun.62.S4$Comparison <- "S4"
mean(jaccard.dist.fun.62.S4$Dissimilarity)

jaccard.dist.fun.62.S <- rbind(jaccard.dist.fun.62.S1,jaccard.dist.fun.62.S2,jaccard.dist.fun.62.S3,
                               jaccard.dist.fun.62.S4)
jaccard.dist.fun.62.S$Days <- 62



S1.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S1"))]
S2.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S2"))]
S3.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S3"))]
S4.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S4"))]
S5.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S5"))]
S6.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S6"))]

jaccard.dist.fun.76.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.76)
jaccard.dist.fun.76.S1 <- rbind(jaccard.dist.fun.76.S1.1, jaccard.dist.fun.76.S1.2)
length(jaccard.dist.fun.76.S1$Dissimilarity)
jaccard.dist.fun.76.S1$Comparison <- "S1"
mean(jaccard.dist.fun.76.S1$Dissimilarity)

jaccard.dist.fun.76.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.76)
jaccard.dist.fun.76.S2 <- rbind(jaccard.dist.fun.76.S2.1, jaccard.dist.fun.76.S2.2)
length(jaccard.dist.fun.76.S2$Dissimilarity)
jaccard.dist.fun.76.S2$Comparison <- "S2"
mean(jaccard.dist.fun.76.S2$Dissimilarity)

jaccard.dist.fun.76.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.76)
jaccard.dist.fun.76.S3 <- rbind(jaccard.dist.fun.76.S3.1, jaccard.dist.fun.76.S3.2)
length(jaccard.dist.fun.76.S3$Dissimilarity)
jaccard.dist.fun.76.S3$Comparison <- "S3"
mean(jaccard.dist.fun.76.S3$Dissimilarity)

jaccard.dist.fun.76.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.76)
jaccard.dist.fun.76.S4 <- rbind(jaccard.dist.fun.76.S4.1, jaccard.dist.fun.76.S4.2)
length(jaccard.dist.fun.76.S4$Dissimilarity)
jaccard.dist.fun.76.S4$Comparison <- "S4"
mean(jaccard.dist.fun.76.S4$Dissimilarity)

jaccard.dist.fun.76.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S5.76)
jaccard.dist.fun.76.S5 <- rbind(jaccard.dist.fun.76.S5.1, jaccard.dist.fun.76.S5.2)
length(jaccard.dist.fun.76.S5$Dissimilarity)
jaccard.dist.fun.76.S5$Comparison <- "S5"
mean(jaccard.dist.fun.76.S5$Dissimilarity)

jaccard.dist.fun.76.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.76 & Var2 %in% G.0)
jaccard.dist.fun.76.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S6.76)
jaccard.dist.fun.76.S6 <- rbind(jaccard.dist.fun.76.S6.1, jaccard.dist.fun.76.S6.2)
length(jaccard.dist.fun.76.S6$Dissimilarity)
jaccard.dist.fun.76.S6$Comparison <- "S6"
mean(jaccard.dist.fun.76.S6$Dissimilarity)

jaccard.dist.fun.76.S <- rbind(jaccard.dist.fun.76.S1,jaccard.dist.fun.76.S2,jaccard.dist.fun.76.S3,
                               jaccard.dist.fun.76.S4,jaccard.dist.fun.76.S5,jaccard.dist.fun.76.S6)
jaccard.dist.fun.76.S$Days <- 76
###90 days
S1.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S1"))]
S2.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S2"))]
S3.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S3"))]
S4.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S4"))]
S5.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S5"))]
S6.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S6"))]
S7.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S7"))]
S8.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S8"))]
S9.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.90.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.90)
jaccard.dist.fun.90.S1 <- rbind(jaccard.dist.fun.90.S1.1, jaccard.dist.fun.90.S1.2)
length(jaccard.dist.fun.90.S1$Dissimilarity)
jaccard.dist.fun.90.S1$Comparison <- "S1"
mean(jaccard.dist.fun.90.S1$Dissimilarity)

jaccard.dist.fun.90.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.90)
jaccard.dist.fun.90.S2 <- rbind(jaccard.dist.fun.90.S2.1, jaccard.dist.fun.90.S2.2)
length(jaccard.dist.fun.90.S2$Dissimilarity)
jaccard.dist.fun.90.S2$Comparison <- "S2"
mean(jaccard.dist.fun.90.S2$Dissimilarity)

jaccard.dist.fun.90.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.90)
jaccard.dist.fun.90.S3 <- rbind(jaccard.dist.fun.90.S3.1, jaccard.dist.fun.90.S3.2)
length(jaccard.dist.fun.90.S3$Dissimilarity)
jaccard.dist.fun.90.S3$Comparison <- "S3"
mean(jaccard.dist.fun.90.S3$Dissimilarity)

jaccard.dist.fun.90.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.90)
jaccard.dist.fun.90.S4 <- rbind(jaccard.dist.fun.90.S4.1, jaccard.dist.fun.90.S4.2)
length(jaccard.dist.fun.90.S4$Dissimilarity)
jaccard.dist.fun.90.S4$Comparison <- "S4"
mean(jaccard.dist.fun.90.S4$Dissimilarity)

jaccard.dist.fun.90.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S5.90)
jaccard.dist.fun.90.S5 <- rbind(jaccard.dist.fun.90.S5.1, jaccard.dist.fun.90.S5.2)
length(jaccard.dist.fun.90.S5$Dissimilarity)
jaccard.dist.fun.90.S5$Comparison <- "S5"
mean(jaccard.dist.fun.90.S5$Dissimilarity)

jaccard.dist.fun.90.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S6.90)
jaccard.dist.fun.90.S6 <- rbind(jaccard.dist.fun.90.S6.1, jaccard.dist.fun.90.S6.2)
length(jaccard.dist.fun.90.S6$Dissimilarity)
jaccard.dist.fun.90.S6$Comparison <- "S6"
mean(jaccard.dist.fun.90.S6$Dissimilarity)

jaccard.dist.fun.90.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S7.90)
jaccard.dist.fun.90.S7 <- rbind(jaccard.dist.fun.90.S7.1, jaccard.dist.fun.90.S7.2)
length(jaccard.dist.fun.90.S7$Dissimilarity)
jaccard.dist.fun.90.S7$Comparison <- "S7"
mean(jaccard.dist.fun.90.S7$Dissimilarity)

jaccard.dist.fun.90.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S8.90)
jaccard.dist.fun.90.S8 <- rbind(jaccard.dist.fun.90.S8.1, jaccard.dist.fun.90.S8.2)
length(jaccard.dist.fun.90.S8$Dissimilarity)
jaccard.dist.fun.90.S8$Comparison <- "S8"
mean(jaccard.dist.fun.90.S8$Dissimilarity)

jaccard.dist.fun.90.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.90 & Var2 %in% G.0)
jaccard.dist.fun.90.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S9.90)
jaccard.dist.fun.90.S9 <- rbind(jaccard.dist.fun.90.S9.1, jaccard.dist.fun.90.S9.2)
length(jaccard.dist.fun.90.S9$Dissimilarity)
jaccard.dist.fun.90.S9$Comparison <- "S9"
mean(jaccard.dist.fun.90.S9$Dissimilarity)


jaccard.dist.fun.90.S <- rbind(jaccard.dist.fun.90.S1,jaccard.dist.fun.90.S2,jaccard.dist.fun.90.S3,
                               jaccard.dist.fun.90.S4,jaccard.dist.fun.90.S5,jaccard.dist.fun.90.S6,
                               jaccard.dist.fun.90.S7,jaccard.dist.fun.90.S8,jaccard.dist.fun.90.S9)

jaccard.dist.fun.90.S$Days <- 90
###106 days
S1.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S1"))]
S2.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S2"))]
S3.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S3"))]
S4.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S4"))]
S5.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S5"))]
S6.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S6"))]
S7.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S7"))]
S8.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S8"))]
S9.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.106.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.106)
jaccard.dist.fun.106.S1 <- rbind(jaccard.dist.fun.106.S1.1, jaccard.dist.fun.106.S1.2)
length(jaccard.dist.fun.106.S1$Dissimilarity)
jaccard.dist.fun.106.S1$Comparison <- "S1"
mean(jaccard.dist.fun.106.S1$Dissimilarity)

jaccard.dist.fun.106.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.106)
jaccard.dist.fun.106.S2 <- rbind(jaccard.dist.fun.106.S2.1, jaccard.dist.fun.106.S2.2)
length(jaccard.dist.fun.106.S2$Dissimilarity)
jaccard.dist.fun.106.S2$Comparison <- "S2"
mean(jaccard.dist.fun.106.S2$Dissimilarity)

jaccard.dist.fun.106.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.106)
jaccard.dist.fun.106.S3 <- rbind(jaccard.dist.fun.106.S3.1, jaccard.dist.fun.106.S3.2)
length(jaccard.dist.fun.106.S3$Dissimilarity)
jaccard.dist.fun.106.S3$Comparison <- "S3"
mean(jaccard.dist.fun.106.S3$Dissimilarity)

jaccard.dist.fun.106.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.106)
jaccard.dist.fun.106.S4 <- rbind(jaccard.dist.fun.106.S4.1, jaccard.dist.fun.106.S4.2)
length(jaccard.dist.fun.106.S4$Dissimilarity)
jaccard.dist.fun.106.S4$Comparison <- "S4"
mean(jaccard.dist.fun.106.S4$Dissimilarity)

jaccard.dist.fun.106.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S5.106)
jaccard.dist.fun.106.S5 <- rbind(jaccard.dist.fun.106.S5.1, jaccard.dist.fun.106.S5.2)
length(jaccard.dist.fun.106.S5$Dissimilarity)
jaccard.dist.fun.106.S5$Comparison <- "S5"
mean(jaccard.dist.fun.106.S5$Dissimilarity)

jaccard.dist.fun.106.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S6.106)
jaccard.dist.fun.106.S6 <- rbind(jaccard.dist.fun.106.S6.1, jaccard.dist.fun.106.S6.2)
length(jaccard.dist.fun.106.S6$Dissimilarity)
jaccard.dist.fun.106.S6$Comparison <- "S6"
mean(jaccard.dist.fun.106.S6$Dissimilarity)

jaccard.dist.fun.106.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S7.106)
jaccard.dist.fun.106.S7 <- rbind(jaccard.dist.fun.106.S7.1, jaccard.dist.fun.106.S7.2)
length(jaccard.dist.fun.106.S7$Dissimilarity)
jaccard.dist.fun.106.S7$Comparison <- "S7"
mean(jaccard.dist.fun.106.S7$Dissimilarity)

jaccard.dist.fun.106.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S8.106)
jaccard.dist.fun.106.S8 <- rbind(jaccard.dist.fun.106.S8.1, jaccard.dist.fun.106.S8.2)
length(jaccard.dist.fun.106.S8$Dissimilarity)
jaccard.dist.fun.106.S8$Comparison <- "S8"
mean(jaccard.dist.fun.106.S8$Dissimilarity)

jaccard.dist.fun.106.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.106 & Var2 %in% G.0)
jaccard.dist.fun.106.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S9.106)
jaccard.dist.fun.106.S9 <- rbind(jaccard.dist.fun.106.S9.1, jaccard.dist.fun.106.S9.2)
length(jaccard.dist.fun.106.S9$Dissimilarity)
jaccard.dist.fun.106.S9$Comparison <- "S9"
mean(jaccard.dist.fun.106.S9$Dissimilarity)


jaccard.dist.fun.106.S <- rbind(jaccard.dist.fun.106.S1,jaccard.dist.fun.106.S2,jaccard.dist.fun.106.S3,
                                jaccard.dist.fun.106.S4,jaccard.dist.fun.106.S5,jaccard.dist.fun.106.S6,
                                jaccard.dist.fun.106.S7,jaccard.dist.fun.106.S8,jaccard.dist.fun.106.S9)
jaccard.dist.fun.106.S$Days <- 106
###120 days
S1.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S1"))]
S2.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S2"))]
S3.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S3"))]
S4.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S4"))]
S5.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S5"))]
S6.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S6"))]
S7.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S7"))]
S8.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S8"))]
S9.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.120.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.120)
jaccard.dist.fun.120.S1 <- rbind(jaccard.dist.fun.120.S1.1, jaccard.dist.fun.120.S1.2)
length(jaccard.dist.fun.120.S1$Dissimilarity)
jaccard.dist.fun.120.S1$Comparison <- "S1"
mean(jaccard.dist.fun.120.S1$Dissimilarity)

jaccard.dist.fun.120.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.120)
jaccard.dist.fun.120.S2 <- rbind(jaccard.dist.fun.120.S2.1, jaccard.dist.fun.120.S2.2)
length(jaccard.dist.fun.120.S2$Dissimilarity)
jaccard.dist.fun.120.S2$Comparison <- "S2"
mean(jaccard.dist.fun.120.S2$Dissimilarity)

jaccard.dist.fun.120.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.120)
jaccard.dist.fun.120.S3 <- rbind(jaccard.dist.fun.120.S3.1, jaccard.dist.fun.120.S3.2)
length(jaccard.dist.fun.120.S3$Dissimilarity)
jaccard.dist.fun.120.S3$Comparison <- "S3"
mean(jaccard.dist.fun.120.S3$Dissimilarity)

jaccard.dist.fun.120.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.120)
jaccard.dist.fun.120.S4 <- rbind(jaccard.dist.fun.120.S4.1, jaccard.dist.fun.120.S4.2)
length(jaccard.dist.fun.120.S4$Dissimilarity)
jaccard.dist.fun.120.S4$Comparison <- "S4"
mean(jaccard.dist.fun.120.S4$Dissimilarity)

jaccard.dist.fun.120.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S5.120)
jaccard.dist.fun.120.S5 <- rbind(jaccard.dist.fun.120.S5.1, jaccard.dist.fun.120.S5.2)
length(jaccard.dist.fun.120.S5$Dissimilarity)
jaccard.dist.fun.120.S5$Comparison <- "S5"
mean(jaccard.dist.fun.120.S5$Dissimilarity)

jaccard.dist.fun.120.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S6.120)
jaccard.dist.fun.120.S6 <- rbind(jaccard.dist.fun.120.S6.1, jaccard.dist.fun.120.S6.2)
length(jaccard.dist.fun.120.S6$Dissimilarity)
jaccard.dist.fun.120.S6$Comparison <- "S6"
mean(jaccard.dist.fun.120.S6$Dissimilarity)

jaccard.dist.fun.120.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S7.120)
jaccard.dist.fun.120.S7 <- rbind(jaccard.dist.fun.120.S7.1, jaccard.dist.fun.120.S7.2)
length(jaccard.dist.fun.120.S7$Dissimilarity)
jaccard.dist.fun.120.S7$Comparison <- "S7"
mean(jaccard.dist.fun.120.S7$Dissimilarity)

jaccard.dist.fun.120.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S8.120)
jaccard.dist.fun.120.S8 <- rbind(jaccard.dist.fun.120.S8.1, jaccard.dist.fun.120.S8.2)
length(jaccard.dist.fun.120.S8$Dissimilarity)
jaccard.dist.fun.120.S8$Comparison <- "S8"
mean(jaccard.dist.fun.120.S8$Dissimilarity)

jaccard.dist.fun.120.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.120 & Var2 %in% G.0)
jaccard.dist.fun.120.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S9.120)
jaccard.dist.fun.120.S9 <- rbind(jaccard.dist.fun.120.S9.1, jaccard.dist.fun.120.S9.2)
length(jaccard.dist.fun.120.S9$Dissimilarity)
jaccard.dist.fun.120.S9$Comparison <- "S9"
mean(jaccard.dist.fun.120.S9$Dissimilarity)


jaccard.dist.fun.120.S <- rbind(jaccard.dist.fun.120.S1,jaccard.dist.fun.120.S2,jaccard.dist.fun.120.S3,
                                jaccard.dist.fun.120.S4,jaccard.dist.fun.120.S5,jaccard.dist.fun.120.S6,
                                jaccard.dist.fun.120.S7,jaccard.dist.fun.120.S8,jaccard.dist.fun.120.S9)
jaccard.dist.fun.120.S$Days <- 120
###141 days
S1.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S1"))]
S2.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S2"))]
S3.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S3"))]
S4.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S4"))]
S5.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S5"))]
S6.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S6"))]
S7.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S7"))]
S8.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S8"))]
S9.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.141.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S1.141)
jaccard.dist.fun.141.S1 <- rbind(jaccard.dist.fun.141.S1.1, jaccard.dist.fun.141.S1.2)
length(jaccard.dist.fun.141.S1$Dissimilarity)
jaccard.dist.fun.141.S1$Comparison <- "S1"
mean(jaccard.dist.fun.141.S1$Dissimilarity)

jaccard.dist.fun.141.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S2.141)
jaccard.dist.fun.141.S2 <- rbind(jaccard.dist.fun.141.S2.1, jaccard.dist.fun.141.S2.2)
length(jaccard.dist.fun.141.S2$Dissimilarity)
jaccard.dist.fun.141.S2$Comparison <- "S2"
mean(jaccard.dist.fun.141.S2$Dissimilarity)

jaccard.dist.fun.141.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S3.141)
jaccard.dist.fun.141.S3 <- rbind(jaccard.dist.fun.141.S3.1, jaccard.dist.fun.141.S3.2)
length(jaccard.dist.fun.141.S3$Dissimilarity)
jaccard.dist.fun.141.S3$Comparison <- "S3"
mean(jaccard.dist.fun.141.S3$Dissimilarity)

jaccard.dist.fun.141.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S4.141)
jaccard.dist.fun.141.S4 <- rbind(jaccard.dist.fun.141.S4.1, jaccard.dist.fun.141.S4.2)
length(jaccard.dist.fun.141.S4$Dissimilarity)
jaccard.dist.fun.141.S4$Comparison <- "S4"
mean(jaccard.dist.fun.141.S4$Dissimilarity)

jaccard.dist.fun.141.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S5.141)
jaccard.dist.fun.141.S5 <- rbind(jaccard.dist.fun.141.S5.1, jaccard.dist.fun.141.S5.2)
length(jaccard.dist.fun.141.S5$Dissimilarity)
jaccard.dist.fun.141.S5$Comparison <- "S5"
mean(jaccard.dist.fun.141.S5$Dissimilarity)

jaccard.dist.fun.141.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S6.141)
jaccard.dist.fun.141.S6 <- rbind(jaccard.dist.fun.141.S6.1, jaccard.dist.fun.141.S6.2)
length(jaccard.dist.fun.141.S6$Dissimilarity)
jaccard.dist.fun.141.S6$Comparison <- "S6"
mean(jaccard.dist.fun.141.S6$Dissimilarity)

jaccard.dist.fun.141.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S7.141)
jaccard.dist.fun.141.S7 <- rbind(jaccard.dist.fun.141.S7.1, jaccard.dist.fun.141.S7.2)
length(jaccard.dist.fun.141.S7$Dissimilarity)
jaccard.dist.fun.141.S7$Comparison <- "S7"
mean(jaccard.dist.fun.141.S7$Dissimilarity)

jaccard.dist.fun.141.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S8.141)
jaccard.dist.fun.141.S8 <- rbind(jaccard.dist.fun.141.S8.1, jaccard.dist.fun.141.S8.2)
length(jaccard.dist.fun.141.S8$Dissimilarity)
jaccard.dist.fun.141.S8$Comparison <- "S8"
mean(jaccard.dist.fun.141.S8$Dissimilarity)

jaccard.dist.fun.141.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.141 & Var2 %in% G.0)
jaccard.dist.fun.141.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% S9.141)
jaccard.dist.fun.141.S9 <- rbind(jaccard.dist.fun.141.S9.1, jaccard.dist.fun.141.S9.2)
length(jaccard.dist.fun.141.S9$Dissimilarity)
jaccard.dist.fun.141.S9$Comparison <- "S9"
mean(jaccard.dist.fun.141.S9$Dissimilarity)


jaccard.dist.fun.141.S <- rbind(jaccard.dist.fun.141.S1,jaccard.dist.fun.141.S2,jaccard.dist.fun.141.S3,
                                jaccard.dist.fun.141.S4,jaccard.dist.fun.141.S5,jaccard.dist.fun.141.S6,
                                jaccard.dist.fun.141.S7,jaccard.dist.fun.141.S8,jaccard.dist.fun.141.S9)

jaccard.dist.fun.141.S$Days <- 141


##Root
R.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.48.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.48 & Var2 %in% G.0)
jaccard.dist.fun.48.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.48)
jaccard.dist.fun.48.R <- rbind(jaccard.dist.fun.48.R.1, jaccard.dist.fun.48.R.2)
length(jaccard.dist.fun.48.R$Dissimilarity)
jaccard.dist.fun.48.R$Comparison <- "R"
jaccard.dist.fun.48.R$Days <- 48
mean(jaccard.dist.fun.48.R$Dissimilarity)

R.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.62.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.62 & Var2 %in% G.0)
jaccard.dist.fun.62.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.62)
jaccard.dist.fun.62.R <- rbind(jaccard.dist.fun.62.R.1, jaccard.dist.fun.62.R.2)
length(jaccard.dist.fun.62.R$Dissimilarity)
jaccard.dist.fun.62.R$Comparison <- "R"
jaccard.dist.fun.62.R$Days <- 62
mean(jaccard.dist.fun.62.R$Dissimilarity)


R.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.76.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.76 & Var2 %in% G.0)
jaccard.dist.fun.76.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.76)
jaccard.dist.fun.76.R <- rbind(jaccard.dist.fun.76.R.1, jaccard.dist.fun.76.R.2)
length(jaccard.dist.fun.76.R$Dissimilarity)
jaccard.dist.fun.76.R$Comparison <- "R"
jaccard.dist.fun.76.R$Days <- 76
mean(jaccard.dist.fun.76.R$Dissimilarity)

R.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.90.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.90 & Var2 %in% G.0)
jaccard.dist.fun.90.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.90)
jaccard.dist.fun.90.R <- rbind(jaccard.dist.fun.90.R.1, jaccard.dist.fun.90.R.2)
length(jaccard.dist.fun.90.R$Dissimilarity)
jaccard.dist.fun.90.R$Comparison <- "R"
jaccard.dist.fun.90.R$Days <- 90
mean(jaccard.dist.fun.90.R$Dissimilarity)


R.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.106.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.106 & Var2 %in% G.0)
jaccard.dist.fun.106.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.106)
jaccard.dist.fun.106.R <- rbind(jaccard.dist.fun.106.R.1, jaccard.dist.fun.106.R.2)
length(jaccard.dist.fun.106.R$Dissimilarity)
jaccard.dist.fun.106.R$Comparison <- "R"
jaccard.dist.fun.106.R$Days <- 106
mean(jaccard.dist.fun.106.R$Dissimilarity)

R.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.120.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.120 & Var2 %in% G.0)
jaccard.dist.fun.120.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.120)
jaccard.dist.fun.120.R <- rbind(jaccard.dist.fun.120.R.1, jaccard.dist.fun.120.R.2)
length(jaccard.dist.fun.120.R$Dissimilarity)
jaccard.dist.fun.120.R$Comparison <- "R"
jaccard.dist.fun.120.R$Days <- 120
mean(jaccard.dist.fun.120.R$Dissimilarity)

R.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.141.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.141 & Var2 %in% G.0)
jaccard.dist.fun.141.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% R.141)
jaccard.dist.fun.141.R <- rbind(jaccard.dist.fun.141.R.1, jaccard.dist.fun.141.R.2)
length(jaccard.dist.fun.141.R$Dissimilarity)
jaccard.dist.fun.141.R$Comparison <- "R"
jaccard.dist.fun.141.R$Days <- 141
mean(jaccard.dist.fun.141.R$Dissimilarity)


##Rhizopshere
RS.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.48.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.48 & Var2 %in% G.0)
jaccard.dist.fun.48.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.48)
jaccard.dist.fun.48.RS <- rbind(jaccard.dist.fun.48.RS.1, jaccard.dist.fun.48.RS.2)
length(jaccard.dist.fun.48.RS$Dissimilarity)
jaccard.dist.fun.48.RS$Comparison <- "RS"
jaccard.dist.fun.48.RS$Days <- 48
mean(jaccard.dist.fun.48.RS$Dissimilarity)

RS.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.62.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.62 & Var2 %in% G.0)
jaccard.dist.fun.62.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.62)
jaccard.dist.fun.62.RS <- rbind(jaccard.dist.fun.62.RS.1, jaccard.dist.fun.62.RS.2)
length(jaccard.dist.fun.62.RS$Dissimilarity)
jaccard.dist.fun.62.RS$Comparison <- "RS"
jaccard.dist.fun.62.RS$Days <- 62
mean(jaccard.dist.fun.62.RS$Dissimilarity)


RS.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.76.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.76 & Var2 %in% G.0)
jaccard.dist.fun.76.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.76)
jaccard.dist.fun.76.RS <- rbind(jaccard.dist.fun.76.RS.1, jaccard.dist.fun.76.RS.2)
length(jaccard.dist.fun.76.RS$Dissimilarity)
jaccard.dist.fun.76.RS$Comparison <- "RS"
jaccard.dist.fun.76.RS$Days <- 76
mean(jaccard.dist.fun.76.RS$Dissimilarity)

RS.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.90.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.90 & Var2 %in% G.0)
jaccard.dist.fun.90.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.90)
jaccard.dist.fun.90.RS <- rbind(jaccard.dist.fun.90.RS.1, jaccard.dist.fun.90.RS.2)
length(jaccard.dist.fun.90.RS$Dissimilarity)
jaccard.dist.fun.90.RS$Comparison <- "RS"
jaccard.dist.fun.90.RS$Days <- 90
mean(jaccard.dist.fun.90.RS$Dissimilarity)


RS.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.106.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.106 & Var2 %in% G.0)
jaccard.dist.fun.106.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.106)
jaccard.dist.fun.106.RS <- rbind(jaccard.dist.fun.106.RS.1, jaccard.dist.fun.106.RS.2)
length(jaccard.dist.fun.106.RS$Dissimilarity)
jaccard.dist.fun.106.RS$Comparison <- "RS"
jaccard.dist.fun.106.RS$Days <- 106
mean(jaccard.dist.fun.106.RS$Dissimilarity)

RS.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.120.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.120 & Var2 %in% G.0)
jaccard.dist.fun.120.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.120)
jaccard.dist.fun.120.RS <- rbind(jaccard.dist.fun.120.RS.1, jaccard.dist.fun.120.RS.2)
length(jaccard.dist.fun.120.RS$Dissimilarity)
jaccard.dist.fun.120.RS$Comparison <- "RS"
jaccard.dist.fun.120.RS$Days <- 120
mean(jaccard.dist.fun.120.RS$Dissimilarity)

RS.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("RS"))]

jaccard.dist.fun.141.RS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% RS.141 & Var2 %in% G.0)
jaccard.dist.fun.141.RS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% RS.141)
jaccard.dist.fun.141.RS <- rbind(jaccard.dist.fun.141.RS.1, jaccard.dist.fun.141.RS.2)
length(jaccard.dist.fun.141.RS$Dissimilarity)
jaccard.dist.fun.141.RS$Comparison <- "RS"
jaccard.dist.fun.141.RS$Days <- 141
mean(jaccard.dist.fun.141.RS$Dissimilarity)

##Bulk soil
BS.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.48.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.48 & Var2 %in% G.0)
jaccard.dist.fun.48.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.48)
jaccard.dist.fun.48.BS <- rbind(jaccard.dist.fun.48.BS.1, jaccard.dist.fun.48.BS.2)
length(jaccard.dist.fun.48.BS$Dissimilarity)
jaccard.dist.fun.48.BS$Comparison <- "BS"
jaccard.dist.fun.48.BS$Days <- 48
mean(jaccard.dist.fun.48.BS$Dissimilarity)

BS.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.62.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.62 & Var2 %in% G.0)
jaccard.dist.fun.62.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.62)
jaccard.dist.fun.62.BS <- rbind(jaccard.dist.fun.62.BS.1, jaccard.dist.fun.62.BS.2)
length(jaccard.dist.fun.62.BS$Dissimilarity)
jaccard.dist.fun.62.BS$Comparison <- "BS"
jaccard.dist.fun.62.BS$Days <- 62
mean(jaccard.dist.fun.62.BS$Dissimilarity)


BS.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.76.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.76 & Var2 %in% G.0)
jaccard.dist.fun.76.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.76)
jaccard.dist.fun.76.BS <- rbind(jaccard.dist.fun.76.BS.1, jaccard.dist.fun.76.BS.2)
length(jaccard.dist.fun.76.BS$Dissimilarity)
jaccard.dist.fun.76.BS$Comparison <- "BS"
jaccard.dist.fun.76.BS$Days <- 76
mean(jaccard.dist.fun.76.BS$Dissimilarity)

BS.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.90.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.90 & Var2 %in% G.0)
jaccard.dist.fun.90.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.90)
jaccard.dist.fun.90.BS <- rbind(jaccard.dist.fun.90.BS.1, jaccard.dist.fun.90.BS.2)
length(jaccard.dist.fun.90.BS$Dissimilarity)
jaccard.dist.fun.90.BS$Comparison <- "BS"
jaccard.dist.fun.90.BS$Days <- 90
mean(jaccard.dist.fun.90.BS$Dissimilarity)


BS.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.106.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.106 & Var2 %in% G.0)
jaccard.dist.fun.106.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.106)
jaccard.dist.fun.106.BS <- rbind(jaccard.dist.fun.106.BS.1, jaccard.dist.fun.106.BS.2)
length(jaccard.dist.fun.106.BS$Dissimilarity)
jaccard.dist.fun.106.BS$Comparison <- "BS"
jaccard.dist.fun.106.BS$Days <- 106
mean(jaccard.dist.fun.106.BS$Dissimilarity)

BS.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.120.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.120 & Var2 %in% G.0)
jaccard.dist.fun.120.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.120)
jaccard.dist.fun.120.BS <- rbind(jaccard.dist.fun.120.BS.1, jaccard.dist.fun.120.BS.2)
length(jaccard.dist.fun.120.BS$Dissimilarity)
jaccard.dist.fun.120.BS$Comparison <- "BS"
jaccard.dist.fun.120.BS$Days <- 120
mean(jaccard.dist.fun.120.BS$Dissimilarity)

BS.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("BS"))]

jaccard.dist.fun.141.BS.1 <- subset(jaccard.dist.fun.melt, Var1 %in% BS.141 & Var2 %in% G.0)
jaccard.dist.fun.141.BS.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% BS.141)
jaccard.dist.fun.141.BS <- rbind(jaccard.dist.fun.141.BS.1, jaccard.dist.fun.141.BS.2)
length(jaccard.dist.fun.141.BS$Dissimilarity)
jaccard.dist.fun.141.BS$Comparison <- "BS"
jaccard.dist.fun.141.BS$Days <- 141
mean(jaccard.dist.fun.141.BS$Dissimilarity)

jaccard.dist.fun.w.G0 <- rbind(jaccard.dist.fun.48.L, jaccard.dist.fun.62.L,jaccard.dist.fun.76.L,
                               jaccard.dist.fun.90.L,jaccard.dist.fun.106.L,jaccard.dist.fun.120.L,
                               jaccard.dist.fun.141.L,jaccard.dist.fun.48.S,jaccard.dist.fun.62.S,
                               jaccard.dist.fun.76.S,jaccard.dist.fun.90.S,jaccard.dist.fun.106.S,
                               jaccard.dist.fun.120.S,jaccard.dist.fun.141.S,jaccard.dist.fun.48.R,
                               jaccard.dist.fun.62.R,jaccard.dist.fun.76.R,jaccard.dist.fun.90.R,
                               jaccard.dist.fun.106.R,jaccard.dist.fun.120.R,jaccard.dist.fun.141.R,
                               jaccard.dist.fun.48.RS,jaccard.dist.fun.62.RS,jaccard.dist.fun.76.RS,
                               jaccard.dist.fun.90.RS,jaccard.dist.fun.106.RS,jaccard.dist.fun.120.RS,
                               jaccard.dist.fun.141.RS,jaccard.dist.fun.48.BS,jaccard.dist.fun.62.BS,
                               jaccard.dist.fun.76.BS,jaccard.dist.fun.90.BS,jaccard.dist.fun.106.BS,
                               jaccard.dist.fun.120.BS,jaccard.dist.fun.141.BS)


###Bacteria

G.0<-map$SampleID[which(map$Days == 0 & map$Microhabitat %in% c("G"))]

## 48
L1.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("L1"))]
L2.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("L2"))]

jaccard.dist.bac.48.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.48 & Var2 %in% G.0)
jaccard.dist.bac.48.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.48)
jaccard.dist.bac.48.L1 <- rbind(jaccard.dist.bac.48.L1.1, jaccard.dist.bac.48.L1.2)
length(jaccard.dist.bac.48.L1$Dissimilarity)
jaccard.dist.bac.48.L1$Comparison <- "L1"
mean(jaccard.dist.bac.48.L1$Dissimilarity)

jaccard.dist.bac.48.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.48 & Var2 %in% G.0)
jaccard.dist.bac.48.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.48)
jaccard.dist.bac.48.L2 <- rbind(jaccard.dist.bac.48.L2.1, jaccard.dist.bac.48.L2.2)
length(jaccard.dist.bac.48.L2$Dissimilarity)
jaccard.dist.bac.48.L2$Comparison <- "L2"
mean(jaccard.dist.bac.48.L2$Dissimilarity)


jaccard.dist.bac.48.L <- rbind(jaccard.dist.bac.48.L1,jaccard.dist.bac.48.L2)
jaccard.dist.bac.48.L$Days <- 48

##62
L1.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L1"))]
L2.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L2"))]
L3.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L3"))]

jaccard.dist.bac.62.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.62 & Var2 %in% G.0)
jaccard.dist.bac.62.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.62)
jaccard.dist.bac.62.L1 <- rbind(jaccard.dist.bac.62.L1.1, jaccard.dist.bac.62.L1.2)
length(jaccard.dist.bac.62.L1$Dissimilarity)
jaccard.dist.bac.62.L1$Comparison <- "L1"
mean(jaccard.dist.bac.62.L1$Dissimilarity)

jaccard.dist.bac.62.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.62 & Var2 %in% G.0)
jaccard.dist.bac.62.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.62)
jaccard.dist.bac.62.L2 <- rbind(jaccard.dist.bac.62.L2.1, jaccard.dist.bac.62.L2.2)
length(jaccard.dist.bac.62.L2$Dissimilarity)
jaccard.dist.bac.62.L2$Comparison <- "L2"
mean(jaccard.dist.bac.62.L2$Dissimilarity)

jaccard.dist.bac.62.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.62 & Var2 %in% G.0)
jaccard.dist.bac.62.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.62)
jaccard.dist.bac.62.L3 <- rbind(jaccard.dist.bac.62.L3.1, jaccard.dist.bac.62.L3.2)
length(jaccard.dist.bac.62.L3$Dissimilarity)
jaccard.dist.bac.62.L3$Comparison <- "L3"
mean(jaccard.dist.bac.62.L3$Dissimilarity)

jaccard.dist.bac.62.L <- rbind(jaccard.dist.bac.62.L1,jaccard.dist.bac.62.L2,jaccard.dist.bac.62.L3)
jaccard.dist.bac.62.L$Days <- 62


#76
L1.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L1"))]
L2.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L2"))]
L3.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L3"))]
FL.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.76.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.76 & Var2 %in% G.0)
jaccard.dist.bac.76.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.76)
jaccard.dist.bac.76.L1 <- rbind(jaccard.dist.bac.76.L1.1, jaccard.dist.bac.76.L1.2)
length(jaccard.dist.bac.76.L1$Dissimilarity)
jaccard.dist.bac.76.L1$Comparison <- "L1"
mean(jaccard.dist.bac.76.L1$Dissimilarity)

jaccard.dist.bac.76.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.76 & Var2 %in% G.0)
jaccard.dist.bac.76.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.76)
jaccard.dist.bac.76.L2 <- rbind(jaccard.dist.bac.76.L2.1, jaccard.dist.bac.76.L2.2)
length(jaccard.dist.bac.76.L2$Dissimilarity)
jaccard.dist.bac.76.L2$Comparison <- "L2"
mean(jaccard.dist.bac.76.L2$Dissimilarity)

jaccard.dist.bac.76.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.76 & Var2 %in% G.0)
jaccard.dist.bac.76.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.76)
jaccard.dist.bac.76.L3 <- rbind(jaccard.dist.bac.76.L3.1, jaccard.dist.bac.76.L3.2)
length(jaccard.dist.bac.76.L3$Dissimilarity)
jaccard.dist.bac.76.L3$Comparison <- "L3"
mean(jaccard.dist.bac.76.L3$Dissimilarity)

jaccard.dist.bac.76.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.76 & Var2 %in% G.0)
jaccard.dist.bac.76.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% FL.76)
jaccard.dist.bac.76.FL <- rbind(jaccard.dist.bac.76.FL.1, jaccard.dist.bac.76.FL.2)
length(jaccard.dist.bac.76.FL$Dissimilarity)
jaccard.dist.bac.76.FL$Comparison <- "FL"
mean(jaccard.dist.bac.76.FL$Dissimilarity)

jaccard.dist.bac.76.L <- rbind(jaccard.dist.bac.76.L1,jaccard.dist.bac.76.L2,jaccard.dist.bac.76.L3,
                               jaccard.dist.bac.76.FL)
jaccard.dist.bac.76.L$Days <- 76

###90 days
L1.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L1"))]
L2.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L2"))]
L3.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L3"))]
FL.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.90.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.90 & Var2 %in% G.0)
jaccard.dist.bac.90.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.90)
jaccard.dist.bac.90.L1 <- rbind(jaccard.dist.bac.90.L1.1, jaccard.dist.bac.90.L1.2)
length(jaccard.dist.bac.90.L1$Dissimilarity)
jaccard.dist.bac.90.L1$Comparison <- "L1"
mean(jaccard.dist.bac.90.L1$Dissimilarity)

jaccard.dist.bac.90.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.90 & Var2 %in% G.0)
jaccard.dist.bac.90.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.90)
jaccard.dist.bac.90.L2 <- rbind(jaccard.dist.bac.90.L2.1, jaccard.dist.bac.90.L2.2)
length(jaccard.dist.bac.90.L2$Dissimilarity)
jaccard.dist.bac.90.L2$Comparison <- "L2"
mean(jaccard.dist.bac.90.L2$Dissimilarity)

jaccard.dist.bac.90.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.90 & Var2 %in% G.0)
jaccard.dist.bac.90.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.90)
jaccard.dist.bac.90.L3 <- rbind(jaccard.dist.bac.90.L3.1, jaccard.dist.bac.90.L3.2)
length(jaccard.dist.bac.90.L3$Dissimilarity)
jaccard.dist.bac.90.L3$Comparison <- "L3"
mean(jaccard.dist.bac.90.L3$Dissimilarity)

jaccard.dist.bac.90.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.90 & Var2 %in% G.0)
jaccard.dist.bac.90.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% FL.90)
jaccard.dist.bac.90.FL <- rbind(jaccard.dist.bac.90.FL.1, jaccard.dist.bac.90.FL.2)
length(jaccard.dist.bac.90.FL$Dissimilarity)
jaccard.dist.bac.90.FL$Comparison <- "FL"
mean(jaccard.dist.bac.90.FL$Dissimilarity)


jaccard.dist.bac.90.L <- rbind(jaccard.dist.bac.90.L1,jaccard.dist.bac.90.L2,jaccard.dist.bac.90.L3,
                               jaccard.dist.bac.90.FL)

jaccard.dist.bac.90.L$Days <- 90

###106 days
L1.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L1"))]
L2.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L2"))]
L3.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L3"))]
FL.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("FL"))]


jaccard.dist.bac.106.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.106 & Var2 %in% G.0)
jaccard.dist.bac.106.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.106)
jaccard.dist.bac.106.L1 <- rbind(jaccard.dist.bac.106.L1.1, jaccard.dist.bac.106.L1.2)
length(jaccard.dist.bac.106.L1$Dissimilarity)
jaccard.dist.bac.106.L1$Comparison <- "L1"
mean(jaccard.dist.bac.106.L1$Dissimilarity)

jaccard.dist.bac.106.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.106 & Var2 %in% G.0)
jaccard.dist.bac.106.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.106)
jaccard.dist.bac.106.L2 <- rbind(jaccard.dist.bac.106.L2.1, jaccard.dist.bac.106.L2.2)
length(jaccard.dist.bac.106.L2$Dissimilarity)
jaccard.dist.bac.106.L2$Comparison <- "L2"
mean(jaccard.dist.bac.106.L2$Dissimilarity)

jaccard.dist.bac.106.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.106 & Var2 %in% G.0)
jaccard.dist.bac.106.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.106)
jaccard.dist.bac.106.L3 <- rbind(jaccard.dist.bac.106.L3.1, jaccard.dist.bac.106.L3.2)
length(jaccard.dist.bac.106.L3$Dissimilarity)
jaccard.dist.bac.106.L3$Comparison <- "L3"
mean(jaccard.dist.bac.106.L3$Dissimilarity)

jaccard.dist.bac.106.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.106 & Var2 %in% G.0)
jaccard.dist.bac.106.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% FL.106)
jaccard.dist.bac.106.FL <- rbind(jaccard.dist.bac.106.FL.1, jaccard.dist.bac.106.FL.2)
length(jaccard.dist.bac.106.FL$Dissimilarity)
jaccard.dist.bac.106.FL$Comparison <- "FL"
mean(jaccard.dist.bac.106.FL$Dissimilarity)

jaccard.dist.bac.106.L <- rbind(jaccard.dist.bac.106.L1,jaccard.dist.bac.106.L2,jaccard.dist.bac.106.L3,
                                jaccard.dist.bac.106.FL)
jaccard.dist.bac.106.L$Days <- 106
###120 days
L1.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L1"))]
L2.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L2"))]
L3.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L3"))]
FL.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.120.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.120 & Var2 %in% G.0)
jaccard.dist.bac.120.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.120)
jaccard.dist.bac.120.L1 <- rbind(jaccard.dist.bac.120.L1.1, jaccard.dist.bac.120.L1.2)
length(jaccard.dist.bac.120.L1$Dissimilarity)
jaccard.dist.bac.120.L1$Comparison <- "L1"
mean(jaccard.dist.bac.120.L1$Dissimilarity)

jaccard.dist.bac.120.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.120 & Var2 %in% G.0)
jaccard.dist.bac.120.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.120)
jaccard.dist.bac.120.L2 <- rbind(jaccard.dist.bac.120.L2.1, jaccard.dist.bac.120.L2.2)
length(jaccard.dist.bac.120.L2$Dissimilarity)
jaccard.dist.bac.120.L2$Comparison <- "L2"
mean(jaccard.dist.bac.120.L2$Dissimilarity)

jaccard.dist.bac.120.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.120 & Var2 %in% G.0)
jaccard.dist.bac.120.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.120)
jaccard.dist.bac.120.L3 <- rbind(jaccard.dist.bac.120.L3.1, jaccard.dist.bac.120.L3.2)
length(jaccard.dist.bac.120.L3$Dissimilarity)
jaccard.dist.bac.120.L3$Comparison <- "L3"
mean(jaccard.dist.bac.120.L3$Dissimilarity)

jaccard.dist.bac.120.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.120 & Var2 %in% G.0)
jaccard.dist.bac.120.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% FL.120)
jaccard.dist.bac.120.FL <- rbind(jaccard.dist.bac.120.FL.1, jaccard.dist.bac.120.FL.2)
length(jaccard.dist.bac.120.FL$Dissimilarity)
jaccard.dist.bac.120.FL$Comparison <- "FL"
mean(jaccard.dist.bac.120.FL$Dissimilarity)

jaccard.dist.bac.120.L <- rbind(jaccard.dist.bac.120.L1,jaccard.dist.bac.120.L2,jaccard.dist.bac.120.L3,
                                jaccard.dist.bac.120.FL)

jaccard.dist.bac.120.L$Days <- 120
###141 days
L1.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L1"))]
L2.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L2"))]
L3.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L3"))]
FL.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.141.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.141 & Var2 %in% G.0)
jaccard.dist.bac.141.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L1.141)
jaccard.dist.bac.141.L1 <- rbind(jaccard.dist.bac.141.L1.1, jaccard.dist.bac.141.L1.2)
length(jaccard.dist.bac.141.L1$Dissimilarity)
jaccard.dist.bac.141.L1$Comparison <- "L1"
mean(jaccard.dist.bac.141.L1$Dissimilarity)

jaccard.dist.bac.141.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.141 & Var2 %in% G.0)
jaccard.dist.bac.141.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L2.141)
jaccard.dist.bac.141.L2 <- rbind(jaccard.dist.bac.141.L2.1, jaccard.dist.bac.141.L2.2)
length(jaccard.dist.bac.141.L2$Dissimilarity)
jaccard.dist.bac.141.L2$Comparison <- "L2"
mean(jaccard.dist.bac.141.L2$Dissimilarity)

jaccard.dist.bac.141.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.141 & Var2 %in% G.0)
jaccard.dist.bac.141.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% L3.141)
jaccard.dist.bac.141.L3 <- rbind(jaccard.dist.bac.141.L3.1, jaccard.dist.bac.141.L3.2)
length(jaccard.dist.bac.141.L3$Dissimilarity)
jaccard.dist.bac.141.L3$Comparison <- "L3"
mean(jaccard.dist.bac.141.L3$Dissimilarity)

jaccard.dist.bac.141.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.141 & Var2 %in% G.0)
jaccard.dist.bac.141.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% FL.141)
jaccard.dist.bac.141.FL <- rbind(jaccard.dist.bac.141.FL.1, jaccard.dist.bac.141.FL.2)
length(jaccard.dist.bac.141.FL$Dissimilarity)
jaccard.dist.bac.141.FL$Comparison <- "FL"
mean(jaccard.dist.bac.141.FL$Dissimilarity)

jaccard.dist.bac.141.L <- rbind(jaccard.dist.bac.141.L1,jaccard.dist.bac.141.L2,jaccard.dist.bac.141.L3,
                                jaccard.dist.bac.141.FL)
jaccard.dist.bac.141.L$Days <- 141
##Stem

S1.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S1"))]
S2.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S2"))]
S3.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S3"))]

jaccard.dist.bac.48.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.48 & Var2 %in% G.0)
jaccard.dist.bac.48.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.48)
jaccard.dist.bac.48.S1 <- rbind(jaccard.dist.bac.48.S1.1, jaccard.dist.bac.48.S1.2)
length(jaccard.dist.bac.48.S1$Dissimilarity)
jaccard.dist.bac.48.S1$Comparison <- "S1"
mean(jaccard.dist.bac.48.S1$Dissimilarity)

jaccard.dist.bac.48.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.48 & Var2 %in% G.0)
jaccard.dist.bac.48.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.48)
jaccard.dist.bac.48.S2 <- rbind(jaccard.dist.bac.48.S2.1, jaccard.dist.bac.48.S2.2)
length(jaccard.dist.bac.48.S2$Dissimilarity)
jaccard.dist.bac.48.S2$Comparison <- "S2"
mean(jaccard.dist.bac.48.S2$Dissimilarity)

jaccard.dist.bac.48.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.48 & Var2 %in% G.0)
jaccard.dist.bac.48.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.48)
jaccard.dist.bac.48.S3 <- rbind(jaccard.dist.bac.48.S3.1, jaccard.dist.bac.48.S3.2)
length(jaccard.dist.bac.48.S3$Dissimilarity)
jaccard.dist.bac.48.S3$Comparison <- "S3"
mean(jaccard.dist.bac.48.S3$Dissimilarity)


jaccard.dist.bac.48.S <- rbind(jaccard.dist.bac.48.S1,jaccard.dist.bac.48.S2,jaccard.dist.bac.48.S3)
jaccard.dist.bac.48.S$Days <- 48


S1.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S1"))]
S2.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S2"))]
S3.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S3"))]
S4.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S4"))]

jaccard.dist.bac.62.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.62 & Var2 %in% G.0)
jaccard.dist.bac.62.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.62)
jaccard.dist.bac.62.S1 <- rbind(jaccard.dist.bac.62.S1.1, jaccard.dist.bac.62.S1.2)
length(jaccard.dist.bac.62.S1$Dissimilarity)
jaccard.dist.bac.62.S1$Comparison <- "S1"
mean(jaccard.dist.bac.62.S1$Dissimilarity)

jaccard.dist.bac.62.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.62 & Var2 %in% G.0)
jaccard.dist.bac.62.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.62)
jaccard.dist.bac.62.S2 <- rbind(jaccard.dist.bac.62.S2.1, jaccard.dist.bac.62.S2.2)
length(jaccard.dist.bac.62.S2$Dissimilarity)
jaccard.dist.bac.62.S2$Comparison <- "S2"
mean(jaccard.dist.bac.62.S2$Dissimilarity)

jaccard.dist.bac.62.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.62 & Var2 %in% G.0)
jaccard.dist.bac.62.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.62)
jaccard.dist.bac.62.S3 <- rbind(jaccard.dist.bac.62.S3.1, jaccard.dist.bac.62.S3.2)
length(jaccard.dist.bac.62.S3$Dissimilarity)
jaccard.dist.bac.62.S3$Comparison <- "S3"
mean(jaccard.dist.bac.62.S3$Dissimilarity)

jaccard.dist.bac.62.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.62 & Var2 %in% G.0)
jaccard.dist.bac.62.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.62)
jaccard.dist.bac.62.S4 <- rbind(jaccard.dist.bac.62.S4.1, jaccard.dist.bac.62.S4.2)
length(jaccard.dist.bac.62.S4$Dissimilarity)
jaccard.dist.bac.62.S4$Comparison <- "S4"
mean(jaccard.dist.bac.62.S4$Dissimilarity)

jaccard.dist.bac.62.S <- rbind(jaccard.dist.bac.62.S1,jaccard.dist.bac.62.S2,jaccard.dist.bac.62.S3,
                               jaccard.dist.bac.62.S4)
jaccard.dist.bac.62.S$Days <- 62



S1.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S1"))]
S2.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S2"))]
S3.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S3"))]
S4.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S4"))]
S5.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S5"))]
S6.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S6"))]

jaccard.dist.bac.76.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.76)
jaccard.dist.bac.76.S1 <- rbind(jaccard.dist.bac.76.S1.1, jaccard.dist.bac.76.S1.2)
length(jaccard.dist.bac.76.S1$Dissimilarity)
jaccard.dist.bac.76.S1$Comparison <- "S1"
mean(jaccard.dist.bac.76.S1$Dissimilarity)

jaccard.dist.bac.76.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.76)
jaccard.dist.bac.76.S2 <- rbind(jaccard.dist.bac.76.S2.1, jaccard.dist.bac.76.S2.2)
length(jaccard.dist.bac.76.S2$Dissimilarity)
jaccard.dist.bac.76.S2$Comparison <- "S2"
mean(jaccard.dist.bac.76.S2$Dissimilarity)

jaccard.dist.bac.76.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.76)
jaccard.dist.bac.76.S3 <- rbind(jaccard.dist.bac.76.S3.1, jaccard.dist.bac.76.S3.2)
length(jaccard.dist.bac.76.S3$Dissimilarity)
jaccard.dist.bac.76.S3$Comparison <- "S3"
mean(jaccard.dist.bac.76.S3$Dissimilarity)

jaccard.dist.bac.76.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.76)
jaccard.dist.bac.76.S4 <- rbind(jaccard.dist.bac.76.S4.1, jaccard.dist.bac.76.S4.2)
length(jaccard.dist.bac.76.S4$Dissimilarity)
jaccard.dist.bac.76.S4$Comparison <- "S4"
mean(jaccard.dist.bac.76.S4$Dissimilarity)

jaccard.dist.bac.76.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S5.76)
jaccard.dist.bac.76.S5 <- rbind(jaccard.dist.bac.76.S5.1, jaccard.dist.bac.76.S5.2)
length(jaccard.dist.bac.76.S5$Dissimilarity)
jaccard.dist.bac.76.S5$Comparison <- "S5"
mean(jaccard.dist.bac.76.S5$Dissimilarity)

jaccard.dist.bac.76.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.76 & Var2 %in% G.0)
jaccard.dist.bac.76.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S6.76)
jaccard.dist.bac.76.S6 <- rbind(jaccard.dist.bac.76.S6.1, jaccard.dist.bac.76.S6.2)
length(jaccard.dist.bac.76.S6$Dissimilarity)
jaccard.dist.bac.76.S6$Comparison <- "S6"
mean(jaccard.dist.bac.76.S6$Dissimilarity)

jaccard.dist.bac.76.S <- rbind(jaccard.dist.bac.76.S1,jaccard.dist.bac.76.S2,jaccard.dist.bac.76.S3,
                               jaccard.dist.bac.76.S4,jaccard.dist.bac.76.S5,jaccard.dist.bac.76.S6)
jaccard.dist.bac.76.S$Days <- 76
###90 days
S1.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S1"))]
S2.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S2"))]
S3.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S3"))]
S4.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S4"))]
S5.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S5"))]
S6.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S6"))]
S7.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S7"))]
S8.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S8"))]
S9.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.90.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.90)
jaccard.dist.bac.90.S1 <- rbind(jaccard.dist.bac.90.S1.1, jaccard.dist.bac.90.S1.2)
length(jaccard.dist.bac.90.S1$Dissimilarity)
jaccard.dist.bac.90.S1$Comparison <- "S1"
mean(jaccard.dist.bac.90.S1$Dissimilarity)

jaccard.dist.bac.90.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.90)
jaccard.dist.bac.90.S2 <- rbind(jaccard.dist.bac.90.S2.1, jaccard.dist.bac.90.S2.2)
length(jaccard.dist.bac.90.S2$Dissimilarity)
jaccard.dist.bac.90.S2$Comparison <- "S2"
mean(jaccard.dist.bac.90.S2$Dissimilarity)

jaccard.dist.bac.90.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.90)
jaccard.dist.bac.90.S3 <- rbind(jaccard.dist.bac.90.S3.1, jaccard.dist.bac.90.S3.2)
length(jaccard.dist.bac.90.S3$Dissimilarity)
jaccard.dist.bac.90.S3$Comparison <- "S3"
mean(jaccard.dist.bac.90.S3$Dissimilarity)

jaccard.dist.bac.90.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.90)
jaccard.dist.bac.90.S4 <- rbind(jaccard.dist.bac.90.S4.1, jaccard.dist.bac.90.S4.2)
length(jaccard.dist.bac.90.S4$Dissimilarity)
jaccard.dist.bac.90.S4$Comparison <- "S4"
mean(jaccard.dist.bac.90.S4$Dissimilarity)

jaccard.dist.bac.90.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S5.90)
jaccard.dist.bac.90.S5 <- rbind(jaccard.dist.bac.90.S5.1, jaccard.dist.bac.90.S5.2)
length(jaccard.dist.bac.90.S5$Dissimilarity)
jaccard.dist.bac.90.S5$Comparison <- "S5"
mean(jaccard.dist.bac.90.S5$Dissimilarity)

jaccard.dist.bac.90.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S6.90)
jaccard.dist.bac.90.S6 <- rbind(jaccard.dist.bac.90.S6.1, jaccard.dist.bac.90.S6.2)
length(jaccard.dist.bac.90.S6$Dissimilarity)
jaccard.dist.bac.90.S6$Comparison <- "S6"
mean(jaccard.dist.bac.90.S6$Dissimilarity)

jaccard.dist.bac.90.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S7.90)
jaccard.dist.bac.90.S7 <- rbind(jaccard.dist.bac.90.S7.1, jaccard.dist.bac.90.S7.2)
length(jaccard.dist.bac.90.S7$Dissimilarity)
jaccard.dist.bac.90.S7$Comparison <- "S7"
mean(jaccard.dist.bac.90.S7$Dissimilarity)

jaccard.dist.bac.90.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S8.90)
jaccard.dist.bac.90.S8 <- rbind(jaccard.dist.bac.90.S8.1, jaccard.dist.bac.90.S8.2)
length(jaccard.dist.bac.90.S8$Dissimilarity)
jaccard.dist.bac.90.S8$Comparison <- "S8"
mean(jaccard.dist.bac.90.S8$Dissimilarity)

jaccard.dist.bac.90.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.90 & Var2 %in% G.0)
jaccard.dist.bac.90.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S9.90)
jaccard.dist.bac.90.S9 <- rbind(jaccard.dist.bac.90.S9.1, jaccard.dist.bac.90.S9.2)
length(jaccard.dist.bac.90.S9$Dissimilarity)
jaccard.dist.bac.90.S9$Comparison <- "S9"
mean(jaccard.dist.bac.90.S9$Dissimilarity)


jaccard.dist.bac.90.S <- rbind(jaccard.dist.bac.90.S1,jaccard.dist.bac.90.S2,jaccard.dist.bac.90.S3,
                               jaccard.dist.bac.90.S4,jaccard.dist.bac.90.S5,jaccard.dist.bac.90.S6,
                               jaccard.dist.bac.90.S7,jaccard.dist.bac.90.S8,jaccard.dist.bac.90.S9)

jaccard.dist.bac.90.S$Days <- 90
###106 days
S1.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S1"))]
S2.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S2"))]
S3.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S3"))]
S4.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S4"))]
S5.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S5"))]
S6.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S6"))]
S7.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S7"))]
S8.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S8"))]
S9.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.106.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.106)
jaccard.dist.bac.106.S1 <- rbind(jaccard.dist.bac.106.S1.1, jaccard.dist.bac.106.S1.2)
length(jaccard.dist.bac.106.S1$Dissimilarity)
jaccard.dist.bac.106.S1$Comparison <- "S1"
mean(jaccard.dist.bac.106.S1$Dissimilarity)

jaccard.dist.bac.106.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.106)
jaccard.dist.bac.106.S2 <- rbind(jaccard.dist.bac.106.S2.1, jaccard.dist.bac.106.S2.2)
length(jaccard.dist.bac.106.S2$Dissimilarity)
jaccard.dist.bac.106.S2$Comparison <- "S2"
mean(jaccard.dist.bac.106.S2$Dissimilarity)

jaccard.dist.bac.106.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.106)
jaccard.dist.bac.106.S3 <- rbind(jaccard.dist.bac.106.S3.1, jaccard.dist.bac.106.S3.2)
length(jaccard.dist.bac.106.S3$Dissimilarity)
jaccard.dist.bac.106.S3$Comparison <- "S3"
mean(jaccard.dist.bac.106.S3$Dissimilarity)

jaccard.dist.bac.106.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.106)
jaccard.dist.bac.106.S4 <- rbind(jaccard.dist.bac.106.S4.1, jaccard.dist.bac.106.S4.2)
length(jaccard.dist.bac.106.S4$Dissimilarity)
jaccard.dist.bac.106.S4$Comparison <- "S4"
mean(jaccard.dist.bac.106.S4$Dissimilarity)

jaccard.dist.bac.106.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S5.106)
jaccard.dist.bac.106.S5 <- rbind(jaccard.dist.bac.106.S5.1, jaccard.dist.bac.106.S5.2)
length(jaccard.dist.bac.106.S5$Dissimilarity)
jaccard.dist.bac.106.S5$Comparison <- "S5"
mean(jaccard.dist.bac.106.S5$Dissimilarity)

jaccard.dist.bac.106.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S6.106)
jaccard.dist.bac.106.S6 <- rbind(jaccard.dist.bac.106.S6.1, jaccard.dist.bac.106.S6.2)
length(jaccard.dist.bac.106.S6$Dissimilarity)
jaccard.dist.bac.106.S6$Comparison <- "S6"
mean(jaccard.dist.bac.106.S6$Dissimilarity)

jaccard.dist.bac.106.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S7.106)
jaccard.dist.bac.106.S7 <- rbind(jaccard.dist.bac.106.S7.1, jaccard.dist.bac.106.S7.2)
length(jaccard.dist.bac.106.S7$Dissimilarity)
jaccard.dist.bac.106.S7$Comparison <- "S7"
mean(jaccard.dist.bac.106.S7$Dissimilarity)

jaccard.dist.bac.106.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S8.106)
jaccard.dist.bac.106.S8 <- rbind(jaccard.dist.bac.106.S8.1, jaccard.dist.bac.106.S8.2)
length(jaccard.dist.bac.106.S8$Dissimilarity)
jaccard.dist.bac.106.S8$Comparison <- "S8"
mean(jaccard.dist.bac.106.S8$Dissimilarity)

jaccard.dist.bac.106.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.106 & Var2 %in% G.0)
jaccard.dist.bac.106.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S9.106)
jaccard.dist.bac.106.S9 <- rbind(jaccard.dist.bac.106.S9.1, jaccard.dist.bac.106.S9.2)
length(jaccard.dist.bac.106.S9$Dissimilarity)
jaccard.dist.bac.106.S9$Comparison <- "S9"
mean(jaccard.dist.bac.106.S9$Dissimilarity)


jaccard.dist.bac.106.S <- rbind(jaccard.dist.bac.106.S1,jaccard.dist.bac.106.S2,jaccard.dist.bac.106.S3,
                                jaccard.dist.bac.106.S4,jaccard.dist.bac.106.S5,jaccard.dist.bac.106.S6,
                                jaccard.dist.bac.106.S7,jaccard.dist.bac.106.S8,jaccard.dist.bac.106.S9)
jaccard.dist.bac.106.S$Days <- 106
###120 days
S1.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S1"))]
S2.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S2"))]
S3.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S3"))]
S4.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S4"))]
S5.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S5"))]
S6.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S6"))]
S7.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S7"))]
S8.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S8"))]
S9.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.120.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.120)
jaccard.dist.bac.120.S1 <- rbind(jaccard.dist.bac.120.S1.1, jaccard.dist.bac.120.S1.2)
length(jaccard.dist.bac.120.S1$Dissimilarity)
jaccard.dist.bac.120.S1$Comparison <- "S1"
mean(jaccard.dist.bac.120.S1$Dissimilarity)

jaccard.dist.bac.120.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.120)
jaccard.dist.bac.120.S2 <- rbind(jaccard.dist.bac.120.S2.1, jaccard.dist.bac.120.S2.2)
length(jaccard.dist.bac.120.S2$Dissimilarity)
jaccard.dist.bac.120.S2$Comparison <- "S2"
mean(jaccard.dist.bac.120.S2$Dissimilarity)

jaccard.dist.bac.120.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.120)
jaccard.dist.bac.120.S3 <- rbind(jaccard.dist.bac.120.S3.1, jaccard.dist.bac.120.S3.2)
length(jaccard.dist.bac.120.S3$Dissimilarity)
jaccard.dist.bac.120.S3$Comparison <- "S3"
mean(jaccard.dist.bac.120.S3$Dissimilarity)

jaccard.dist.bac.120.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.120)
jaccard.dist.bac.120.S4 <- rbind(jaccard.dist.bac.120.S4.1, jaccard.dist.bac.120.S4.2)
length(jaccard.dist.bac.120.S4$Dissimilarity)
jaccard.dist.bac.120.S4$Comparison <- "S4"
mean(jaccard.dist.bac.120.S4$Dissimilarity)

jaccard.dist.bac.120.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S5.120)
jaccard.dist.bac.120.S5 <- rbind(jaccard.dist.bac.120.S5.1, jaccard.dist.bac.120.S5.2)
length(jaccard.dist.bac.120.S5$Dissimilarity)
jaccard.dist.bac.120.S5$Comparison <- "S5"
mean(jaccard.dist.bac.120.S5$Dissimilarity)

jaccard.dist.bac.120.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S6.120)
jaccard.dist.bac.120.S6 <- rbind(jaccard.dist.bac.120.S6.1, jaccard.dist.bac.120.S6.2)
length(jaccard.dist.bac.120.S6$Dissimilarity)
jaccard.dist.bac.120.S6$Comparison <- "S6"
mean(jaccard.dist.bac.120.S6$Dissimilarity)

jaccard.dist.bac.120.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S7.120)
jaccard.dist.bac.120.S7 <- rbind(jaccard.dist.bac.120.S7.1, jaccard.dist.bac.120.S7.2)
length(jaccard.dist.bac.120.S7$Dissimilarity)
jaccard.dist.bac.120.S7$Comparison <- "S7"
mean(jaccard.dist.bac.120.S7$Dissimilarity)

jaccard.dist.bac.120.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S8.120)
jaccard.dist.bac.120.S8 <- rbind(jaccard.dist.bac.120.S8.1, jaccard.dist.bac.120.S8.2)
length(jaccard.dist.bac.120.S8$Dissimilarity)
jaccard.dist.bac.120.S8$Comparison <- "S8"
mean(jaccard.dist.bac.120.S8$Dissimilarity)

jaccard.dist.bac.120.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.120 & Var2 %in% G.0)
jaccard.dist.bac.120.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S9.120)
jaccard.dist.bac.120.S9 <- rbind(jaccard.dist.bac.120.S9.1, jaccard.dist.bac.120.S9.2)
length(jaccard.dist.bac.120.S9$Dissimilarity)
jaccard.dist.bac.120.S9$Comparison <- "S9"
mean(jaccard.dist.bac.120.S9$Dissimilarity)


jaccard.dist.bac.120.S <- rbind(jaccard.dist.bac.120.S1,jaccard.dist.bac.120.S2,jaccard.dist.bac.120.S3,
                                jaccard.dist.bac.120.S4,jaccard.dist.bac.120.S5,jaccard.dist.bac.120.S6,
                                jaccard.dist.bac.120.S7,jaccard.dist.bac.120.S8,jaccard.dist.bac.120.S9)
jaccard.dist.bac.120.S$Days <- 120
###141 days
S1.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S1"))]
S2.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S2"))]
S3.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S3"))]
S4.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S4"))]
S5.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S5"))]
S6.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S6"))]
S7.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S7"))]
S8.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S8"))]
S9.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.141.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S1.141)
jaccard.dist.bac.141.S1 <- rbind(jaccard.dist.bac.141.S1.1, jaccard.dist.bac.141.S1.2)
length(jaccard.dist.bac.141.S1$Dissimilarity)
jaccard.dist.bac.141.S1$Comparison <- "S1"
mean(jaccard.dist.bac.141.S1$Dissimilarity)

jaccard.dist.bac.141.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S2.141)
jaccard.dist.bac.141.S2 <- rbind(jaccard.dist.bac.141.S2.1, jaccard.dist.bac.141.S2.2)
length(jaccard.dist.bac.141.S2$Dissimilarity)
jaccard.dist.bac.141.S2$Comparison <- "S2"
mean(jaccard.dist.bac.141.S2$Dissimilarity)

jaccard.dist.bac.141.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S3.141)
jaccard.dist.bac.141.S3 <- rbind(jaccard.dist.bac.141.S3.1, jaccard.dist.bac.141.S3.2)
length(jaccard.dist.bac.141.S3$Dissimilarity)
jaccard.dist.bac.141.S3$Comparison <- "S3"
mean(jaccard.dist.bac.141.S3$Dissimilarity)

jaccard.dist.bac.141.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S4.141)
jaccard.dist.bac.141.S4 <- rbind(jaccard.dist.bac.141.S4.1, jaccard.dist.bac.141.S4.2)
length(jaccard.dist.bac.141.S4$Dissimilarity)
jaccard.dist.bac.141.S4$Comparison <- "S4"
mean(jaccard.dist.bac.141.S4$Dissimilarity)

jaccard.dist.bac.141.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S5.141)
jaccard.dist.bac.141.S5 <- rbind(jaccard.dist.bac.141.S5.1, jaccard.dist.bac.141.S5.2)
length(jaccard.dist.bac.141.S5$Dissimilarity)
jaccard.dist.bac.141.S5$Comparison <- "S5"
mean(jaccard.dist.bac.141.S5$Dissimilarity)

jaccard.dist.bac.141.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S6.141)
jaccard.dist.bac.141.S6 <- rbind(jaccard.dist.bac.141.S6.1, jaccard.dist.bac.141.S6.2)
length(jaccard.dist.bac.141.S6$Dissimilarity)
jaccard.dist.bac.141.S6$Comparison <- "S6"
mean(jaccard.dist.bac.141.S6$Dissimilarity)

jaccard.dist.bac.141.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S7.141)
jaccard.dist.bac.141.S7 <- rbind(jaccard.dist.bac.141.S7.1, jaccard.dist.bac.141.S7.2)
length(jaccard.dist.bac.141.S7$Dissimilarity)
jaccard.dist.bac.141.S7$Comparison <- "S7"
mean(jaccard.dist.bac.141.S7$Dissimilarity)

jaccard.dist.bac.141.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S8.141)
jaccard.dist.bac.141.S8 <- rbind(jaccard.dist.bac.141.S8.1, jaccard.dist.bac.141.S8.2)
length(jaccard.dist.bac.141.S8$Dissimilarity)
jaccard.dist.bac.141.S8$Comparison <- "S8"
mean(jaccard.dist.bac.141.S8$Dissimilarity)

jaccard.dist.bac.141.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.141 & Var2 %in% G.0)
jaccard.dist.bac.141.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% S9.141)
jaccard.dist.bac.141.S9 <- rbind(jaccard.dist.bac.141.S9.1, jaccard.dist.bac.141.S9.2)
length(jaccard.dist.bac.141.S9$Dissimilarity)
jaccard.dist.bac.141.S9$Comparison <- "S9"
mean(jaccard.dist.bac.141.S9$Dissimilarity)


jaccard.dist.bac.141.S <- rbind(jaccard.dist.bac.141.S1,jaccard.dist.bac.141.S2,jaccard.dist.bac.141.S3,
                                jaccard.dist.bac.141.S4,jaccard.dist.bac.141.S5,jaccard.dist.bac.141.S6,
                                jaccard.dist.bac.141.S7,jaccard.dist.bac.141.S8,jaccard.dist.bac.141.S9)

jaccard.dist.bac.141.S$Days <- 141


##Root
R.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.48.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.48 & Var2 %in% G.0)
jaccard.dist.bac.48.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.48)
jaccard.dist.bac.48.R <- rbind(jaccard.dist.bac.48.R.1, jaccard.dist.bac.48.R.2)
length(jaccard.dist.bac.48.R$Dissimilarity)
jaccard.dist.bac.48.R$Comparison <- "R"
jaccard.dist.bac.48.R$Days <- 48
mean(jaccard.dist.bac.48.R$Dissimilarity)

R.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.62.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.62 & Var2 %in% G.0)
jaccard.dist.bac.62.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.62)
jaccard.dist.bac.62.R <- rbind(jaccard.dist.bac.62.R.1, jaccard.dist.bac.62.R.2)
length(jaccard.dist.bac.62.R$Dissimilarity)
jaccard.dist.bac.62.R$Comparison <- "R"
jaccard.dist.bac.62.R$Days <- 62
mean(jaccard.dist.bac.62.R$Dissimilarity)


R.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.76.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.76 & Var2 %in% G.0)
jaccard.dist.bac.76.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.76)
jaccard.dist.bac.76.R <- rbind(jaccard.dist.bac.76.R.1, jaccard.dist.bac.76.R.2)
length(jaccard.dist.bac.76.R$Dissimilarity)
jaccard.dist.bac.76.R$Comparison <- "R"
jaccard.dist.bac.76.R$Days <- 76
mean(jaccard.dist.bac.76.R$Dissimilarity)

R.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.90.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.90 & Var2 %in% G.0)
jaccard.dist.bac.90.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.90)
jaccard.dist.bac.90.R <- rbind(jaccard.dist.bac.90.R.1, jaccard.dist.bac.90.R.2)
length(jaccard.dist.bac.90.R$Dissimilarity)
jaccard.dist.bac.90.R$Comparison <- "R"
jaccard.dist.bac.90.R$Days <- 90
mean(jaccard.dist.bac.90.R$Dissimilarity)


R.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.106.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.106 & Var2 %in% G.0)
jaccard.dist.bac.106.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.106)
jaccard.dist.bac.106.R <- rbind(jaccard.dist.bac.106.R.1, jaccard.dist.bac.106.R.2)
length(jaccard.dist.bac.106.R$Dissimilarity)
jaccard.dist.bac.106.R$Comparison <- "R"
jaccard.dist.bac.106.R$Days <- 106
mean(jaccard.dist.bac.106.R$Dissimilarity)

R.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.120.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.120 & Var2 %in% G.0)
jaccard.dist.bac.120.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.120)
jaccard.dist.bac.120.R <- rbind(jaccard.dist.bac.120.R.1, jaccard.dist.bac.120.R.2)
length(jaccard.dist.bac.120.R$Dissimilarity)
jaccard.dist.bac.120.R$Comparison <- "R"
jaccard.dist.bac.120.R$Days <- 120
mean(jaccard.dist.bac.120.R$Dissimilarity)

R.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.141.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.141 & Var2 %in% G.0)
jaccard.dist.bac.141.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% R.141)
jaccard.dist.bac.141.R <- rbind(jaccard.dist.bac.141.R.1, jaccard.dist.bac.141.R.2)
length(jaccard.dist.bac.141.R$Dissimilarity)
jaccard.dist.bac.141.R$Comparison <- "R"
jaccard.dist.bac.141.R$Days <- 141
mean(jaccard.dist.bac.141.R$Dissimilarity)


RS.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.48.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.48 & Var2 %in% G.0)
jaccard.dist.bac.48.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.48)
jaccard.dist.bac.48.RS <- rbind(jaccard.dist.bac.48.RS.1, jaccard.dist.bac.48.RS.2)
length(jaccard.dist.bac.48.RS$Dissimilarity)
jaccard.dist.bac.48.RS$Comparison <- "RS"
jaccard.dist.bac.48.RS$Days <- 48
mean(jaccard.dist.bac.48.RS$Dissimilarity)

RS.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.62.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.62 & Var2 %in% G.0)
jaccard.dist.bac.62.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.62)
jaccard.dist.bac.62.RS <- rbind(jaccard.dist.bac.62.RS.1, jaccard.dist.bac.62.RS.2)
length(jaccard.dist.bac.62.RS$Dissimilarity)
jaccard.dist.bac.62.RS$Comparison <- "RS"
jaccard.dist.bac.62.RS$Days <- 62
mean(jaccard.dist.bac.62.RS$Dissimilarity)


RS.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.76.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.76 & Var2 %in% G.0)
jaccard.dist.bac.76.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.76)
jaccard.dist.bac.76.RS <- rbind(jaccard.dist.bac.76.RS.1, jaccard.dist.bac.76.RS.2)
length(jaccard.dist.bac.76.RS$Dissimilarity)
jaccard.dist.bac.76.RS$Comparison <- "RS"
jaccard.dist.bac.76.RS$Days <- 76
mean(jaccard.dist.bac.76.RS$Dissimilarity)

RS.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.90.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.90 & Var2 %in% G.0)
jaccard.dist.bac.90.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.90)
jaccard.dist.bac.90.RS <- rbind(jaccard.dist.bac.90.RS.1, jaccard.dist.bac.90.RS.2)
length(jaccard.dist.bac.90.RS$Dissimilarity)
jaccard.dist.bac.90.RS$Comparison <- "RS"
jaccard.dist.bac.90.RS$Days <- 90
mean(jaccard.dist.bac.90.RS$Dissimilarity)


RS.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.106.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.106 & Var2 %in% G.0)
jaccard.dist.bac.106.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.106)
jaccard.dist.bac.106.RS <- rbind(jaccard.dist.bac.106.RS.1, jaccard.dist.bac.106.RS.2)
length(jaccard.dist.bac.106.RS$Dissimilarity)
jaccard.dist.bac.106.RS$Comparison <- "RS"
jaccard.dist.bac.106.RS$Days <- 106
mean(jaccard.dist.bac.106.RS$Dissimilarity)

RS.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.120.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.120 & Var2 %in% G.0)
jaccard.dist.bac.120.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.120)
jaccard.dist.bac.120.RS <- rbind(jaccard.dist.bac.120.RS.1, jaccard.dist.bac.120.RS.2)
length(jaccard.dist.bac.120.RS$Dissimilarity)
jaccard.dist.bac.120.RS$Comparison <- "RS"
jaccard.dist.bac.120.RS$Days <- 120
mean(jaccard.dist.bac.120.RS$Dissimilarity)

RS.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("RS"))]

jaccard.dist.bac.141.RS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% RS.141 & Var2 %in% G.0)
jaccard.dist.bac.141.RS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% RS.141)
jaccard.dist.bac.141.RS <- rbind(jaccard.dist.bac.141.RS.1, jaccard.dist.bac.141.RS.2)
length(jaccard.dist.bac.141.RS$Dissimilarity)
jaccard.dist.bac.141.RS$Comparison <- "RS"
jaccard.dist.bac.141.RS$Days <- 141
mean(jaccard.dist.bac.141.RS$Dissimilarity)

##Bulk soil
BS.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.48.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.48 & Var2 %in% G.0)
jaccard.dist.bac.48.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.48)
jaccard.dist.bac.48.BS <- rbind(jaccard.dist.bac.48.BS.1, jaccard.dist.bac.48.BS.2)
length(jaccard.dist.bac.48.BS$Dissimilarity)
jaccard.dist.bac.48.BS$Comparison <- "BS"
jaccard.dist.bac.48.BS$Days <- 48
mean(jaccard.dist.bac.48.BS$Dissimilarity)

BS.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.62.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.62 & Var2 %in% G.0)
jaccard.dist.bac.62.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.62)
jaccard.dist.bac.62.BS <- rbind(jaccard.dist.bac.62.BS.1, jaccard.dist.bac.62.BS.2)
length(jaccard.dist.bac.62.BS$Dissimilarity)
jaccard.dist.bac.62.BS$Comparison <- "BS"
jaccard.dist.bac.62.BS$Days <- 62
mean(jaccard.dist.bac.62.BS$Dissimilarity)


BS.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.76.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.76 & Var2 %in% G.0)
jaccard.dist.bac.76.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.76)
jaccard.dist.bac.76.BS <- rbind(jaccard.dist.bac.76.BS.1, jaccard.dist.bac.76.BS.2)
length(jaccard.dist.bac.76.BS$Dissimilarity)
jaccard.dist.bac.76.BS$Comparison <- "BS"
jaccard.dist.bac.76.BS$Days <- 76
mean(jaccard.dist.bac.76.BS$Dissimilarity)

BS.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.90.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.90 & Var2 %in% G.0)
jaccard.dist.bac.90.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.90)
jaccard.dist.bac.90.BS <- rbind(jaccard.dist.bac.90.BS.1, jaccard.dist.bac.90.BS.2)
length(jaccard.dist.bac.90.BS$Dissimilarity)
jaccard.dist.bac.90.BS$Comparison <- "BS"
jaccard.dist.bac.90.BS$Days <- 90
mean(jaccard.dist.bac.90.BS$Dissimilarity)


BS.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.106.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.106 & Var2 %in% G.0)
jaccard.dist.bac.106.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.106)
jaccard.dist.bac.106.BS <- rbind(jaccard.dist.bac.106.BS.1, jaccard.dist.bac.106.BS.2)
length(jaccard.dist.bac.106.BS$Dissimilarity)
jaccard.dist.bac.106.BS$Comparison <- "BS"
jaccard.dist.bac.106.BS$Days <- 106
mean(jaccard.dist.bac.106.BS$Dissimilarity)

BS.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.120.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.120 & Var2 %in% G.0)
jaccard.dist.bac.120.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.120)
jaccard.dist.bac.120.BS <- rbind(jaccard.dist.bac.120.BS.1, jaccard.dist.bac.120.BS.2)
length(jaccard.dist.bac.120.BS$Dissimilarity)
jaccard.dist.bac.120.BS$Comparison <- "BS"
jaccard.dist.bac.120.BS$Days <- 120
mean(jaccard.dist.bac.120.BS$Dissimilarity)

BS.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("BS"))]

jaccard.dist.bac.141.BS.1 <- subset(jaccard.dist.bac.melt, Var1 %in% BS.141 & Var2 %in% G.0)
jaccard.dist.bac.141.BS.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% BS.141)
jaccard.dist.bac.141.BS <- rbind(jaccard.dist.bac.141.BS.1, jaccard.dist.bac.141.BS.2)
length(jaccard.dist.bac.141.BS$Dissimilarity)
jaccard.dist.bac.141.BS$Comparison <- "BS"
jaccard.dist.bac.141.BS$Days <- 141
mean(jaccard.dist.bac.141.BS$Dissimilarity)



jaccard.dist.bac.w.G0 <- rbind(jaccard.dist.bac.48.L, jaccard.dist.bac.62.L,jaccard.dist.bac.76.L,
                               jaccard.dist.bac.90.L,jaccard.dist.bac.106.L,jaccard.dist.bac.120.L,
                               jaccard.dist.bac.141.L,jaccard.dist.bac.48.S,jaccard.dist.bac.62.S,
                               jaccard.dist.bac.76.S,jaccard.dist.bac.90.S,jaccard.dist.bac.106.S,
                               jaccard.dist.bac.120.S,jaccard.dist.bac.141.S,jaccard.dist.bac.48.R,
                               jaccard.dist.bac.62.R,jaccard.dist.bac.76.R,jaccard.dist.bac.90.R,
                               jaccard.dist.bac.106.R,jaccard.dist.bac.120.R,jaccard.dist.bac.141.R,
                               jaccard.dist.bac.48.RS,jaccard.dist.bac.62.RS,jaccard.dist.bac.76.RS,
                               jaccard.dist.bac.90.RS,jaccard.dist.bac.106.RS,jaccard.dist.bac.120.RS,
                               jaccard.dist.bac.141.RS,jaccard.dist.bac.48.BS,jaccard.dist.bac.62.BS,
                               jaccard.dist.bac.76.BS,jaccard.dist.bac.90.BS,jaccard.dist.bac.106.BS,
                               jaccard.dist.bac.120.BS,jaccard.dist.bac.141.BS)




write.csv(jaccard.dist.bac.w.G0,"jaccard.dist.bac.w.G0.csv")
write.csv(jaccard.dist.fun.w.G0,"jaccard.dist.fun.w.G0.csv")


#### Distance with G141
##Fungi

G.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("G"))]

## 48
L1.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("L1"))]
L2.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("L2"))]

jaccard.dist.fun.48.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.48 & Var2 %in% G.141)
jaccard.dist.fun.48.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.48)
jaccard.dist.fun.48.L1 <- rbind(jaccard.dist.fun.48.L1.1, jaccard.dist.fun.48.L1.2)
length(jaccard.dist.fun.48.L1$Dissimilarity)
jaccard.dist.fun.48.L1$Comparison <- "L1"
mean(jaccard.dist.fun.48.L1$Dissimilarity)

jaccard.dist.fun.48.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.48 & Var2 %in% G.141)
jaccard.dist.fun.48.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.48)
jaccard.dist.fun.48.L2 <- rbind(jaccard.dist.fun.48.L2.1, jaccard.dist.fun.48.L2.2)
length(jaccard.dist.fun.48.L2$Dissimilarity)
jaccard.dist.fun.48.L2$Comparison <- "L2"
mean(jaccard.dist.fun.48.L2$Dissimilarity)


jaccard.dist.fun.48.L <- rbind(jaccard.dist.fun.48.L1,jaccard.dist.fun.48.L2)
jaccard.dist.fun.48.L$Days <- 48

##62
L1.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L1"))]
L2.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L2"))]
L3.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("L3"))]

jaccard.dist.fun.62.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.62 & Var2 %in% G.141)
jaccard.dist.fun.62.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.62)
jaccard.dist.fun.62.L1 <- rbind(jaccard.dist.fun.62.L1.1, jaccard.dist.fun.62.L1.2)
length(jaccard.dist.fun.62.L1$Dissimilarity)
jaccard.dist.fun.62.L1$Comparison <- "L1"
mean(jaccard.dist.fun.62.L1$Dissimilarity)

jaccard.dist.fun.62.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.62 & Var2 %in% G.141)
jaccard.dist.fun.62.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.62)
jaccard.dist.fun.62.L2 <- rbind(jaccard.dist.fun.62.L2.1, jaccard.dist.fun.62.L2.2)
length(jaccard.dist.fun.62.L2$Dissimilarity)
jaccard.dist.fun.62.L2$Comparison <- "L2"
mean(jaccard.dist.fun.62.L2$Dissimilarity)

jaccard.dist.fun.62.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.62 & Var2 %in% G.141)
jaccard.dist.fun.62.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.62)
jaccard.dist.fun.62.L3 <- rbind(jaccard.dist.fun.62.L3.1, jaccard.dist.fun.62.L3.2)
length(jaccard.dist.fun.62.L3$Dissimilarity)
jaccard.dist.fun.62.L3$Comparison <- "L3"
mean(jaccard.dist.fun.62.L3$Dissimilarity)

jaccard.dist.fun.62.L <- rbind(jaccard.dist.fun.62.L1,jaccard.dist.fun.62.L2,jaccard.dist.fun.62.L3)
jaccard.dist.fun.62.L$Days <- 62


#76
L1.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L1"))]
L2.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L2"))]
L3.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("L3"))]
FL.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.76.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.76 & Var2 %in% G.141)
jaccard.dist.fun.76.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.76)
jaccard.dist.fun.76.L1 <- rbind(jaccard.dist.fun.76.L1.1, jaccard.dist.fun.76.L1.2)
length(jaccard.dist.fun.76.L1$Dissimilarity)
jaccard.dist.fun.76.L1$Comparison <- "L1"
mean(jaccard.dist.fun.76.L1$Dissimilarity)

jaccard.dist.fun.76.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.76 & Var2 %in% G.141)
jaccard.dist.fun.76.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.76)
jaccard.dist.fun.76.L2 <- rbind(jaccard.dist.fun.76.L2.1, jaccard.dist.fun.76.L2.2)
length(jaccard.dist.fun.76.L2$Dissimilarity)
jaccard.dist.fun.76.L2$Comparison <- "L2"
mean(jaccard.dist.fun.76.L2$Dissimilarity)

jaccard.dist.fun.76.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.76 & Var2 %in% G.141)
jaccard.dist.fun.76.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.76)
jaccard.dist.fun.76.L3 <- rbind(jaccard.dist.fun.76.L3.1, jaccard.dist.fun.76.L3.2)
length(jaccard.dist.fun.76.L3$Dissimilarity)
jaccard.dist.fun.76.L3$Comparison <- "L3"
mean(jaccard.dist.fun.76.L3$Dissimilarity)

jaccard.dist.fun.76.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.76 & Var2 %in% G.141)
jaccard.dist.fun.76.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% FL.76)
jaccard.dist.fun.76.FL <- rbind(jaccard.dist.fun.76.FL.1, jaccard.dist.fun.76.FL.2)
length(jaccard.dist.fun.76.FL$Dissimilarity)
jaccard.dist.fun.76.FL$Comparison <- "FL"
mean(jaccard.dist.fun.76.FL$Dissimilarity)

jaccard.dist.fun.76.L <- rbind(jaccard.dist.fun.76.L1,jaccard.dist.fun.76.L2,jaccard.dist.fun.76.L3,
                               jaccard.dist.fun.76.FL)
jaccard.dist.fun.76.L$Days <- 76

###90 days
L1.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L1"))]
L2.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L2"))]
L3.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("L3"))]
FL.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.90.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.90 & Var2 %in% G.141)
jaccard.dist.fun.90.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.90)
jaccard.dist.fun.90.L1 <- rbind(jaccard.dist.fun.90.L1.1, jaccard.dist.fun.90.L1.2)
length(jaccard.dist.fun.90.L1$Dissimilarity)
jaccard.dist.fun.90.L1$Comparison <- "L1"
mean(jaccard.dist.fun.90.L1$Dissimilarity)

jaccard.dist.fun.90.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.90 & Var2 %in% G.141)
jaccard.dist.fun.90.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.90)
jaccard.dist.fun.90.L2 <- rbind(jaccard.dist.fun.90.L2.1, jaccard.dist.fun.90.L2.2)
length(jaccard.dist.fun.90.L2$Dissimilarity)
jaccard.dist.fun.90.L2$Comparison <- "L2"
mean(jaccard.dist.fun.90.L2$Dissimilarity)

jaccard.dist.fun.90.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.90 & Var2 %in% G.141)
jaccard.dist.fun.90.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.90)
jaccard.dist.fun.90.L3 <- rbind(jaccard.dist.fun.90.L3.1, jaccard.dist.fun.90.L3.2)
length(jaccard.dist.fun.90.L3$Dissimilarity)
jaccard.dist.fun.90.L3$Comparison <- "L3"
mean(jaccard.dist.fun.90.L3$Dissimilarity)

jaccard.dist.fun.90.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.90 & Var2 %in% G.141)
jaccard.dist.fun.90.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% FL.90)
jaccard.dist.fun.90.FL <- rbind(jaccard.dist.fun.90.FL.1, jaccard.dist.fun.90.FL.2)
length(jaccard.dist.fun.90.FL$Dissimilarity)
jaccard.dist.fun.90.FL$Comparison <- "FL"
mean(jaccard.dist.fun.90.FL$Dissimilarity)


jaccard.dist.fun.90.L <- rbind(jaccard.dist.fun.90.L1,jaccard.dist.fun.90.L2,jaccard.dist.fun.90.L3,
                               jaccard.dist.fun.90.FL)

jaccard.dist.fun.90.L$Days <- 90

###106 days
L1.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L1"))]
L2.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L2"))]
L3.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("L3"))]
FL.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("FL"))]


jaccard.dist.fun.106.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.106 & Var2 %in% G.141)
jaccard.dist.fun.106.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.106)
jaccard.dist.fun.106.L1 <- rbind(jaccard.dist.fun.106.L1.1, jaccard.dist.fun.106.L1.2)
length(jaccard.dist.fun.106.L1$Dissimilarity)
jaccard.dist.fun.106.L1$Comparison <- "L1"
mean(jaccard.dist.fun.106.L1$Dissimilarity)

jaccard.dist.fun.106.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.106 & Var2 %in% G.141)
jaccard.dist.fun.106.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.106)
jaccard.dist.fun.106.L2 <- rbind(jaccard.dist.fun.106.L2.1, jaccard.dist.fun.106.L2.2)
length(jaccard.dist.fun.106.L2$Dissimilarity)
jaccard.dist.fun.106.L2$Comparison <- "L2"
mean(jaccard.dist.fun.106.L2$Dissimilarity)

jaccard.dist.fun.106.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.106 & Var2 %in% G.141)
jaccard.dist.fun.106.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.106)
jaccard.dist.fun.106.L3 <- rbind(jaccard.dist.fun.106.L3.1, jaccard.dist.fun.106.L3.2)
length(jaccard.dist.fun.106.L3$Dissimilarity)
jaccard.dist.fun.106.L3$Comparison <- "L3"
mean(jaccard.dist.fun.106.L3$Dissimilarity)

jaccard.dist.fun.106.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.106 & Var2 %in% G.141)
jaccard.dist.fun.106.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% FL.106)
jaccard.dist.fun.106.FL <- rbind(jaccard.dist.fun.106.FL.1, jaccard.dist.fun.106.FL.2)
length(jaccard.dist.fun.106.FL$Dissimilarity)
jaccard.dist.fun.106.FL$Comparison <- "FL"
mean(jaccard.dist.fun.106.FL$Dissimilarity)

jaccard.dist.fun.106.L <- rbind(jaccard.dist.fun.106.L1,jaccard.dist.fun.106.L2,jaccard.dist.fun.106.L3,
                                jaccard.dist.fun.106.FL)
jaccard.dist.fun.106.L$Days <- 106
###120 days
L1.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L1"))]
L2.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L2"))]
L3.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("L3"))]
FL.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.120.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.120 & Var2 %in% G.141)
jaccard.dist.fun.120.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.120)
jaccard.dist.fun.120.L1 <- rbind(jaccard.dist.fun.120.L1.1, jaccard.dist.fun.120.L1.2)
length(jaccard.dist.fun.120.L1$Dissimilarity)
jaccard.dist.fun.120.L1$Comparison <- "L1"
mean(jaccard.dist.fun.120.L1$Dissimilarity)

jaccard.dist.fun.120.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.120 & Var2 %in% G.141)
jaccard.dist.fun.120.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.120)
jaccard.dist.fun.120.L2 <- rbind(jaccard.dist.fun.120.L2.1, jaccard.dist.fun.120.L2.2)
length(jaccard.dist.fun.120.L2$Dissimilarity)
jaccard.dist.fun.120.L2$Comparison <- "L2"
mean(jaccard.dist.fun.120.L2$Dissimilarity)

jaccard.dist.fun.120.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.120 & Var2 %in% G.141)
jaccard.dist.fun.120.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.120)
jaccard.dist.fun.120.L3 <- rbind(jaccard.dist.fun.120.L3.1, jaccard.dist.fun.120.L3.2)
length(jaccard.dist.fun.120.L3$Dissimilarity)
jaccard.dist.fun.120.L3$Comparison <- "L3"
mean(jaccard.dist.fun.120.L3$Dissimilarity)

jaccard.dist.fun.120.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.120 & Var2 %in% G.141)
jaccard.dist.fun.120.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% FL.120)
jaccard.dist.fun.120.FL <- rbind(jaccard.dist.fun.120.FL.1, jaccard.dist.fun.120.FL.2)
length(jaccard.dist.fun.120.FL$Dissimilarity)
jaccard.dist.fun.120.FL$Comparison <- "FL"
mean(jaccard.dist.fun.120.FL$Dissimilarity)

jaccard.dist.fun.120.L <- rbind(jaccard.dist.fun.120.L1,jaccard.dist.fun.120.L2,jaccard.dist.fun.120.L3,
                                jaccard.dist.fun.120.FL)

jaccard.dist.fun.120.L$Days <- 120
###141 days
L1.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L1"))]
L2.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L2"))]
L3.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("L3"))]
FL.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("FL"))]

jaccard.dist.fun.141.L1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L1.141 & Var2 %in% G.141)
jaccard.dist.fun.141.L1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L1.141)
jaccard.dist.fun.141.L1 <- rbind(jaccard.dist.fun.141.L1.1, jaccard.dist.fun.141.L1.2)
length(jaccard.dist.fun.141.L1$Dissimilarity)
jaccard.dist.fun.141.L1$Comparison <- "L1"
mean(jaccard.dist.fun.141.L1$Dissimilarity)

jaccard.dist.fun.141.L2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L2.141 & Var2 %in% G.141)
jaccard.dist.fun.141.L2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L2.141)
jaccard.dist.fun.141.L2 <- rbind(jaccard.dist.fun.141.L2.1, jaccard.dist.fun.141.L2.2)
length(jaccard.dist.fun.141.L2$Dissimilarity)
jaccard.dist.fun.141.L2$Comparison <- "L2"
mean(jaccard.dist.fun.141.L2$Dissimilarity)

jaccard.dist.fun.141.L3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% L3.141 & Var2 %in% G.141)
jaccard.dist.fun.141.L3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% L3.141)
jaccard.dist.fun.141.L3 <- rbind(jaccard.dist.fun.141.L3.1, jaccard.dist.fun.141.L3.2)
length(jaccard.dist.fun.141.L3$Dissimilarity)
jaccard.dist.fun.141.L3$Comparison <- "L3"
mean(jaccard.dist.fun.141.L3$Dissimilarity)

jaccard.dist.fun.141.FL.1 <- subset(jaccard.dist.fun.melt, Var1 %in% FL.141 & Var2 %in% G.141)
jaccard.dist.fun.141.FL.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% FL.141)
jaccard.dist.fun.141.FL <- rbind(jaccard.dist.fun.141.FL.1, jaccard.dist.fun.141.FL.2)
length(jaccard.dist.fun.141.FL$Dissimilarity)
jaccard.dist.fun.141.FL$Comparison <- "FL"
mean(jaccard.dist.fun.141.FL$Dissimilarity)

jaccard.dist.fun.141.L <- rbind(jaccard.dist.fun.141.L1,jaccard.dist.fun.141.L2,jaccard.dist.fun.141.L3,
                                jaccard.dist.fun.141.FL)
jaccard.dist.fun.141.L$Days <- 141
##Stem

S1.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S1"))]
S2.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S2"))]
S3.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("S3"))]

jaccard.dist.fun.48.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.48 & Var2 %in% G.141)
jaccard.dist.fun.48.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.48)
jaccard.dist.fun.48.S1 <- rbind(jaccard.dist.fun.48.S1.1, jaccard.dist.fun.48.S1.2)
length(jaccard.dist.fun.48.S1$Dissimilarity)
jaccard.dist.fun.48.S1$Comparison <- "S1"
mean(jaccard.dist.fun.48.S1$Dissimilarity)

jaccard.dist.fun.48.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.48 & Var2 %in% G.141)
jaccard.dist.fun.48.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.48)
jaccard.dist.fun.48.S2 <- rbind(jaccard.dist.fun.48.S2.1, jaccard.dist.fun.48.S2.2)
length(jaccard.dist.fun.48.S2$Dissimilarity)
jaccard.dist.fun.48.S2$Comparison <- "S2"
mean(jaccard.dist.fun.48.S2$Dissimilarity)

jaccard.dist.fun.48.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.48 & Var2 %in% G.141)
jaccard.dist.fun.48.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.48)
jaccard.dist.fun.48.S3 <- rbind(jaccard.dist.fun.48.S3.1, jaccard.dist.fun.48.S3.2)
length(jaccard.dist.fun.48.S3$Dissimilarity)
jaccard.dist.fun.48.S3$Comparison <- "S3"
mean(jaccard.dist.fun.48.S3$Dissimilarity)


jaccard.dist.fun.48.S <- rbind(jaccard.dist.fun.48.S1,jaccard.dist.fun.48.S2,jaccard.dist.fun.48.S3)
jaccard.dist.fun.48.S$Days <- 48


S1.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S1"))]
S2.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S2"))]
S3.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S3"))]
S4.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("S4"))]

jaccard.dist.fun.62.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.62 & Var2 %in% G.141)
jaccard.dist.fun.62.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.62)
jaccard.dist.fun.62.S1 <- rbind(jaccard.dist.fun.62.S1.1, jaccard.dist.fun.62.S1.2)
length(jaccard.dist.fun.62.S1$Dissimilarity)
jaccard.dist.fun.62.S1$Comparison <- "S1"
mean(jaccard.dist.fun.62.S1$Dissimilarity)

jaccard.dist.fun.62.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.62 & Var2 %in% G.141)
jaccard.dist.fun.62.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.62)
jaccard.dist.fun.62.S2 <- rbind(jaccard.dist.fun.62.S2.1, jaccard.dist.fun.62.S2.2)
length(jaccard.dist.fun.62.S2$Dissimilarity)
jaccard.dist.fun.62.S2$Comparison <- "S2"
mean(jaccard.dist.fun.62.S2$Dissimilarity)

jaccard.dist.fun.62.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.62 & Var2 %in% G.141)
jaccard.dist.fun.62.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.62)
jaccard.dist.fun.62.S3 <- rbind(jaccard.dist.fun.62.S3.1, jaccard.dist.fun.62.S3.2)
length(jaccard.dist.fun.62.S3$Dissimilarity)
jaccard.dist.fun.62.S3$Comparison <- "S3"
mean(jaccard.dist.fun.62.S3$Dissimilarity)

jaccard.dist.fun.62.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.62 & Var2 %in% G.141)
jaccard.dist.fun.62.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.62)
jaccard.dist.fun.62.S4 <- rbind(jaccard.dist.fun.62.S4.1, jaccard.dist.fun.62.S4.2)
length(jaccard.dist.fun.62.S4$Dissimilarity)
jaccard.dist.fun.62.S4$Comparison <- "S4"
mean(jaccard.dist.fun.62.S4$Dissimilarity)

jaccard.dist.fun.62.S <- rbind(jaccard.dist.fun.62.S1,jaccard.dist.fun.62.S2,jaccard.dist.fun.62.S3,
                               jaccard.dist.fun.62.S4)
jaccard.dist.fun.62.S$Days <- 62



S1.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S1"))]
S2.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S2"))]
S3.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S3"))]
S4.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S4"))]
S5.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S5"))]
S6.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("S6"))]

jaccard.dist.fun.76.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.76)
jaccard.dist.fun.76.S1 <- rbind(jaccard.dist.fun.76.S1.1, jaccard.dist.fun.76.S1.2)
length(jaccard.dist.fun.76.S1$Dissimilarity)
jaccard.dist.fun.76.S1$Comparison <- "S1"
mean(jaccard.dist.fun.76.S1$Dissimilarity)

jaccard.dist.fun.76.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.76)
jaccard.dist.fun.76.S2 <- rbind(jaccard.dist.fun.76.S2.1, jaccard.dist.fun.76.S2.2)
length(jaccard.dist.fun.76.S2$Dissimilarity)
jaccard.dist.fun.76.S2$Comparison <- "S2"
mean(jaccard.dist.fun.76.S2$Dissimilarity)

jaccard.dist.fun.76.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.76)
jaccard.dist.fun.76.S3 <- rbind(jaccard.dist.fun.76.S3.1, jaccard.dist.fun.76.S3.2)
length(jaccard.dist.fun.76.S3$Dissimilarity)
jaccard.dist.fun.76.S3$Comparison <- "S3"
mean(jaccard.dist.fun.76.S3$Dissimilarity)

jaccard.dist.fun.76.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.76)
jaccard.dist.fun.76.S4 <- rbind(jaccard.dist.fun.76.S4.1, jaccard.dist.fun.76.S4.2)
length(jaccard.dist.fun.76.S4$Dissimilarity)
jaccard.dist.fun.76.S4$Comparison <- "S4"
mean(jaccard.dist.fun.76.S4$Dissimilarity)

jaccard.dist.fun.76.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S5.76)
jaccard.dist.fun.76.S5 <- rbind(jaccard.dist.fun.76.S5.1, jaccard.dist.fun.76.S5.2)
length(jaccard.dist.fun.76.S5$Dissimilarity)
jaccard.dist.fun.76.S5$Comparison <- "S5"
mean(jaccard.dist.fun.76.S5$Dissimilarity)

jaccard.dist.fun.76.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.76 & Var2 %in% G.141)
jaccard.dist.fun.76.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S6.76)
jaccard.dist.fun.76.S6 <- rbind(jaccard.dist.fun.76.S6.1, jaccard.dist.fun.76.S6.2)
length(jaccard.dist.fun.76.S6$Dissimilarity)
jaccard.dist.fun.76.S6$Comparison <- "S6"
mean(jaccard.dist.fun.76.S6$Dissimilarity)

jaccard.dist.fun.76.S <- rbind(jaccard.dist.fun.76.S1,jaccard.dist.fun.76.S2,jaccard.dist.fun.76.S3,
                               jaccard.dist.fun.76.S4,jaccard.dist.fun.76.S5,jaccard.dist.fun.76.S6)
jaccard.dist.fun.76.S$Days <- 76
###90 days
S1.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S1"))]
S2.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S2"))]
S3.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S3"))]
S4.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S4"))]
S5.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S5"))]
S6.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S6"))]
S7.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S7"))]
S8.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S8"))]
S9.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.90.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.90)
jaccard.dist.fun.90.S1 <- rbind(jaccard.dist.fun.90.S1.1, jaccard.dist.fun.90.S1.2)
length(jaccard.dist.fun.90.S1$Dissimilarity)
jaccard.dist.fun.90.S1$Comparison <- "S1"
mean(jaccard.dist.fun.90.S1$Dissimilarity)

jaccard.dist.fun.90.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.90)
jaccard.dist.fun.90.S2 <- rbind(jaccard.dist.fun.90.S2.1, jaccard.dist.fun.90.S2.2)
length(jaccard.dist.fun.90.S2$Dissimilarity)
jaccard.dist.fun.90.S2$Comparison <- "S2"
mean(jaccard.dist.fun.90.S2$Dissimilarity)

jaccard.dist.fun.90.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.90)
jaccard.dist.fun.90.S3 <- rbind(jaccard.dist.fun.90.S3.1, jaccard.dist.fun.90.S3.2)
length(jaccard.dist.fun.90.S3$Dissimilarity)
jaccard.dist.fun.90.S3$Comparison <- "S3"
mean(jaccard.dist.fun.90.S3$Dissimilarity)

jaccard.dist.fun.90.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.90)
jaccard.dist.fun.90.S4 <- rbind(jaccard.dist.fun.90.S4.1, jaccard.dist.fun.90.S4.2)
length(jaccard.dist.fun.90.S4$Dissimilarity)
jaccard.dist.fun.90.S4$Comparison <- "S4"
mean(jaccard.dist.fun.90.S4$Dissimilarity)

jaccard.dist.fun.90.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S5.90)
jaccard.dist.fun.90.S5 <- rbind(jaccard.dist.fun.90.S5.1, jaccard.dist.fun.90.S5.2)
length(jaccard.dist.fun.90.S5$Dissimilarity)
jaccard.dist.fun.90.S5$Comparison <- "S5"
mean(jaccard.dist.fun.90.S5$Dissimilarity)

jaccard.dist.fun.90.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S6.90)
jaccard.dist.fun.90.S6 <- rbind(jaccard.dist.fun.90.S6.1, jaccard.dist.fun.90.S6.2)
length(jaccard.dist.fun.90.S6$Dissimilarity)
jaccard.dist.fun.90.S6$Comparison <- "S6"
mean(jaccard.dist.fun.90.S6$Dissimilarity)

jaccard.dist.fun.90.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S7.90)
jaccard.dist.fun.90.S7 <- rbind(jaccard.dist.fun.90.S7.1, jaccard.dist.fun.90.S7.2)
length(jaccard.dist.fun.90.S7$Dissimilarity)
jaccard.dist.fun.90.S7$Comparison <- "S7"
mean(jaccard.dist.fun.90.S7$Dissimilarity)

jaccard.dist.fun.90.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S8.90)
jaccard.dist.fun.90.S8 <- rbind(jaccard.dist.fun.90.S8.1, jaccard.dist.fun.90.S8.2)
length(jaccard.dist.fun.90.S8$Dissimilarity)
jaccard.dist.fun.90.S8$Comparison <- "S8"
mean(jaccard.dist.fun.90.S8$Dissimilarity)

jaccard.dist.fun.90.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.90 & Var2 %in% G.141)
jaccard.dist.fun.90.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S9.90)
jaccard.dist.fun.90.S9 <- rbind(jaccard.dist.fun.90.S9.1, jaccard.dist.fun.90.S9.2)
length(jaccard.dist.fun.90.S9$Dissimilarity)
jaccard.dist.fun.90.S9$Comparison <- "S9"
mean(jaccard.dist.fun.90.S9$Dissimilarity)


jaccard.dist.fun.90.S <- rbind(jaccard.dist.fun.90.S1,jaccard.dist.fun.90.S2,jaccard.dist.fun.90.S3,
                               jaccard.dist.fun.90.S4,jaccard.dist.fun.90.S5,jaccard.dist.fun.90.S6,
                               jaccard.dist.fun.90.S7,jaccard.dist.fun.90.S8,jaccard.dist.fun.90.S9)

jaccard.dist.fun.90.S$Days <- 90
###106 days
S1.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S1"))]
S2.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S2"))]
S3.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S3"))]
S4.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S4"))]
S5.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S5"))]
S6.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S6"))]
S7.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S7"))]
S8.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S8"))]
S9.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.106.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.106)
jaccard.dist.fun.106.S1 <- rbind(jaccard.dist.fun.106.S1.1, jaccard.dist.fun.106.S1.2)
length(jaccard.dist.fun.106.S1$Dissimilarity)
jaccard.dist.fun.106.S1$Comparison <- "S1"
mean(jaccard.dist.fun.106.S1$Dissimilarity)

jaccard.dist.fun.106.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.106)
jaccard.dist.fun.106.S2 <- rbind(jaccard.dist.fun.106.S2.1, jaccard.dist.fun.106.S2.2)
length(jaccard.dist.fun.106.S2$Dissimilarity)
jaccard.dist.fun.106.S2$Comparison <- "S2"
mean(jaccard.dist.fun.106.S2$Dissimilarity)

jaccard.dist.fun.106.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.106)
jaccard.dist.fun.106.S3 <- rbind(jaccard.dist.fun.106.S3.1, jaccard.dist.fun.106.S3.2)
length(jaccard.dist.fun.106.S3$Dissimilarity)
jaccard.dist.fun.106.S3$Comparison <- "S3"
mean(jaccard.dist.fun.106.S3$Dissimilarity)

jaccard.dist.fun.106.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.106)
jaccard.dist.fun.106.S4 <- rbind(jaccard.dist.fun.106.S4.1, jaccard.dist.fun.106.S4.2)
length(jaccard.dist.fun.106.S4$Dissimilarity)
jaccard.dist.fun.106.S4$Comparison <- "S4"
mean(jaccard.dist.fun.106.S4$Dissimilarity)

jaccard.dist.fun.106.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S5.106)
jaccard.dist.fun.106.S5 <- rbind(jaccard.dist.fun.106.S5.1, jaccard.dist.fun.106.S5.2)
length(jaccard.dist.fun.106.S5$Dissimilarity)
jaccard.dist.fun.106.S5$Comparison <- "S5"
mean(jaccard.dist.fun.106.S5$Dissimilarity)

jaccard.dist.fun.106.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S6.106)
jaccard.dist.fun.106.S6 <- rbind(jaccard.dist.fun.106.S6.1, jaccard.dist.fun.106.S6.2)
length(jaccard.dist.fun.106.S6$Dissimilarity)
jaccard.dist.fun.106.S6$Comparison <- "S6"
mean(jaccard.dist.fun.106.S6$Dissimilarity)

jaccard.dist.fun.106.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S7.106)
jaccard.dist.fun.106.S7 <- rbind(jaccard.dist.fun.106.S7.1, jaccard.dist.fun.106.S7.2)
length(jaccard.dist.fun.106.S7$Dissimilarity)
jaccard.dist.fun.106.S7$Comparison <- "S7"
mean(jaccard.dist.fun.106.S7$Dissimilarity)

jaccard.dist.fun.106.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S8.106)
jaccard.dist.fun.106.S8 <- rbind(jaccard.dist.fun.106.S8.1, jaccard.dist.fun.106.S8.2)
length(jaccard.dist.fun.106.S8$Dissimilarity)
jaccard.dist.fun.106.S8$Comparison <- "S8"
mean(jaccard.dist.fun.106.S8$Dissimilarity)

jaccard.dist.fun.106.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.106 & Var2 %in% G.141)
jaccard.dist.fun.106.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S9.106)
jaccard.dist.fun.106.S9 <- rbind(jaccard.dist.fun.106.S9.1, jaccard.dist.fun.106.S9.2)
length(jaccard.dist.fun.106.S9$Dissimilarity)
jaccard.dist.fun.106.S9$Comparison <- "S9"
mean(jaccard.dist.fun.106.S9$Dissimilarity)


jaccard.dist.fun.106.S <- rbind(jaccard.dist.fun.106.S1,jaccard.dist.fun.106.S2,jaccard.dist.fun.106.S3,
                                jaccard.dist.fun.106.S4,jaccard.dist.fun.106.S5,jaccard.dist.fun.106.S6,
                                jaccard.dist.fun.106.S7,jaccard.dist.fun.106.S8,jaccard.dist.fun.106.S9)
jaccard.dist.fun.106.S$Days <- 106
###120 days
S1.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S1"))]
S2.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S2"))]
S3.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S3"))]
S4.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S4"))]
S5.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S5"))]
S6.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S6"))]
S7.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S7"))]
S8.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S8"))]
S9.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.120.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.120)
jaccard.dist.fun.120.S1 <- rbind(jaccard.dist.fun.120.S1.1, jaccard.dist.fun.120.S1.2)
length(jaccard.dist.fun.120.S1$Dissimilarity)
jaccard.dist.fun.120.S1$Comparison <- "S1"
mean(jaccard.dist.fun.120.S1$Dissimilarity)

jaccard.dist.fun.120.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.120)
jaccard.dist.fun.120.S2 <- rbind(jaccard.dist.fun.120.S2.1, jaccard.dist.fun.120.S2.2)
length(jaccard.dist.fun.120.S2$Dissimilarity)
jaccard.dist.fun.120.S2$Comparison <- "S2"
mean(jaccard.dist.fun.120.S2$Dissimilarity)

jaccard.dist.fun.120.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.120)
jaccard.dist.fun.120.S3 <- rbind(jaccard.dist.fun.120.S3.1, jaccard.dist.fun.120.S3.2)
length(jaccard.dist.fun.120.S3$Dissimilarity)
jaccard.dist.fun.120.S3$Comparison <- "S3"
mean(jaccard.dist.fun.120.S3$Dissimilarity)

jaccard.dist.fun.120.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.120)
jaccard.dist.fun.120.S4 <- rbind(jaccard.dist.fun.120.S4.1, jaccard.dist.fun.120.S4.2)
length(jaccard.dist.fun.120.S4$Dissimilarity)
jaccard.dist.fun.120.S4$Comparison <- "S4"
mean(jaccard.dist.fun.120.S4$Dissimilarity)

jaccard.dist.fun.120.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S5.120)
jaccard.dist.fun.120.S5 <- rbind(jaccard.dist.fun.120.S5.1, jaccard.dist.fun.120.S5.2)
length(jaccard.dist.fun.120.S5$Dissimilarity)
jaccard.dist.fun.120.S5$Comparison <- "S5"
mean(jaccard.dist.fun.120.S5$Dissimilarity)

jaccard.dist.fun.120.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S6.120)
jaccard.dist.fun.120.S6 <- rbind(jaccard.dist.fun.120.S6.1, jaccard.dist.fun.120.S6.2)
length(jaccard.dist.fun.120.S6$Dissimilarity)
jaccard.dist.fun.120.S6$Comparison <- "S6"
mean(jaccard.dist.fun.120.S6$Dissimilarity)

jaccard.dist.fun.120.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S7.120)
jaccard.dist.fun.120.S7 <- rbind(jaccard.dist.fun.120.S7.1, jaccard.dist.fun.120.S7.2)
length(jaccard.dist.fun.120.S7$Dissimilarity)
jaccard.dist.fun.120.S7$Comparison <- "S7"
mean(jaccard.dist.fun.120.S7$Dissimilarity)

jaccard.dist.fun.120.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S8.120)
jaccard.dist.fun.120.S8 <- rbind(jaccard.dist.fun.120.S8.1, jaccard.dist.fun.120.S8.2)
length(jaccard.dist.fun.120.S8$Dissimilarity)
jaccard.dist.fun.120.S8$Comparison <- "S8"
mean(jaccard.dist.fun.120.S8$Dissimilarity)

jaccard.dist.fun.120.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.120 & Var2 %in% G.141)
jaccard.dist.fun.120.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S9.120)
jaccard.dist.fun.120.S9 <- rbind(jaccard.dist.fun.120.S9.1, jaccard.dist.fun.120.S9.2)
length(jaccard.dist.fun.120.S9$Dissimilarity)
jaccard.dist.fun.120.S9$Comparison <- "S9"
mean(jaccard.dist.fun.120.S9$Dissimilarity)


jaccard.dist.fun.120.S <- rbind(jaccard.dist.fun.120.S1,jaccard.dist.fun.120.S2,jaccard.dist.fun.120.S3,
                                jaccard.dist.fun.120.S4,jaccard.dist.fun.120.S5,jaccard.dist.fun.120.S6,
                                jaccard.dist.fun.120.S7,jaccard.dist.fun.120.S8,jaccard.dist.fun.120.S9)
jaccard.dist.fun.120.S$Days <- 120
###141 days
S1.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S1"))]
S2.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S2"))]
S3.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S3"))]
S4.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S4"))]
S5.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S5"))]
S6.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S6"))]
S7.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S7"))]
S8.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S8"))]
S9.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("S9"))]

jaccard.dist.fun.141.S1.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S1.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S1.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S1.141)
jaccard.dist.fun.141.S1 <- rbind(jaccard.dist.fun.141.S1.1, jaccard.dist.fun.141.S1.2)
length(jaccard.dist.fun.141.S1$Dissimilarity)
jaccard.dist.fun.141.S1$Comparison <- "S1"
mean(jaccard.dist.fun.141.S1$Dissimilarity)

jaccard.dist.fun.141.S2.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S2.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S2.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S2.141)
jaccard.dist.fun.141.S2 <- rbind(jaccard.dist.fun.141.S2.1, jaccard.dist.fun.141.S2.2)
length(jaccard.dist.fun.141.S2$Dissimilarity)
jaccard.dist.fun.141.S2$Comparison <- "S2"
mean(jaccard.dist.fun.141.S2$Dissimilarity)

jaccard.dist.fun.141.S3.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S3.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S3.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S3.141)
jaccard.dist.fun.141.S3 <- rbind(jaccard.dist.fun.141.S3.1, jaccard.dist.fun.141.S3.2)
length(jaccard.dist.fun.141.S3$Dissimilarity)
jaccard.dist.fun.141.S3$Comparison <- "S3"
mean(jaccard.dist.fun.141.S3$Dissimilarity)

jaccard.dist.fun.141.S4.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S4.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S4.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S4.141)
jaccard.dist.fun.141.S4 <- rbind(jaccard.dist.fun.141.S4.1, jaccard.dist.fun.141.S4.2)
length(jaccard.dist.fun.141.S4$Dissimilarity)
jaccard.dist.fun.141.S4$Comparison <- "S4"
mean(jaccard.dist.fun.141.S4$Dissimilarity)

jaccard.dist.fun.141.S5.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S5.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S5.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S5.141)
jaccard.dist.fun.141.S5 <- rbind(jaccard.dist.fun.141.S5.1, jaccard.dist.fun.141.S5.2)
length(jaccard.dist.fun.141.S5$Dissimilarity)
jaccard.dist.fun.141.S5$Comparison <- "S5"
mean(jaccard.dist.fun.141.S5$Dissimilarity)

jaccard.dist.fun.141.S6.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S6.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S6.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S6.141)
jaccard.dist.fun.141.S6 <- rbind(jaccard.dist.fun.141.S6.1, jaccard.dist.fun.141.S6.2)
length(jaccard.dist.fun.141.S6$Dissimilarity)
jaccard.dist.fun.141.S6$Comparison <- "S6"
mean(jaccard.dist.fun.141.S6$Dissimilarity)

jaccard.dist.fun.141.S7.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S7.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S7.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S7.141)
jaccard.dist.fun.141.S7 <- rbind(jaccard.dist.fun.141.S7.1, jaccard.dist.fun.141.S7.2)
length(jaccard.dist.fun.141.S7$Dissimilarity)
jaccard.dist.fun.141.S7$Comparison <- "S7"
mean(jaccard.dist.fun.141.S7$Dissimilarity)

jaccard.dist.fun.141.S8.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S8.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S8.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S8.141)
jaccard.dist.fun.141.S8 <- rbind(jaccard.dist.fun.141.S8.1, jaccard.dist.fun.141.S8.2)
length(jaccard.dist.fun.141.S8$Dissimilarity)
jaccard.dist.fun.141.S8$Comparison <- "S8"
mean(jaccard.dist.fun.141.S8$Dissimilarity)

jaccard.dist.fun.141.S9.1 <- subset(jaccard.dist.fun.melt, Var1 %in% S9.141 & Var2 %in% G.141)
jaccard.dist.fun.141.S9.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% S9.141)
jaccard.dist.fun.141.S9 <- rbind(jaccard.dist.fun.141.S9.1, jaccard.dist.fun.141.S9.2)
length(jaccard.dist.fun.141.S9$Dissimilarity)
jaccard.dist.fun.141.S9$Comparison <- "S9"
mean(jaccard.dist.fun.141.S9$Dissimilarity)


jaccard.dist.fun.141.S <- rbind(jaccard.dist.fun.141.S1,jaccard.dist.fun.141.S2,jaccard.dist.fun.141.S3,
                                jaccard.dist.fun.141.S4,jaccard.dist.fun.141.S5,jaccard.dist.fun.141.S6,
                                jaccard.dist.fun.141.S7,jaccard.dist.fun.141.S8,jaccard.dist.fun.141.S9)

jaccard.dist.fun.141.S$Days <- 141


##Root
R.48<-f.meta.18$SampleID[which(f.meta.18$Days == 48 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.48.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.48 & Var2 %in% G.141)
jaccard.dist.fun.48.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.48)
jaccard.dist.fun.48.R <- rbind(jaccard.dist.fun.48.R.1, jaccard.dist.fun.48.R.2)
length(jaccard.dist.fun.48.R$Dissimilarity)
jaccard.dist.fun.48.R$Comparison <- "R"
jaccard.dist.fun.48.R$Days <- 48
mean(jaccard.dist.fun.48.R$Dissimilarity)

R.62<-f.meta.18$SampleID[which(f.meta.18$Days == 62 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.62.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.62 & Var2 %in% G.141)
jaccard.dist.fun.62.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.62)
jaccard.dist.fun.62.R <- rbind(jaccard.dist.fun.62.R.1, jaccard.dist.fun.62.R.2)
length(jaccard.dist.fun.62.R$Dissimilarity)
jaccard.dist.fun.62.R$Comparison <- "R"
jaccard.dist.fun.62.R$Days <- 62
mean(jaccard.dist.fun.62.R$Dissimilarity)


R.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.76.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.76 & Var2 %in% G.141)
jaccard.dist.fun.76.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.76)
jaccard.dist.fun.76.R <- rbind(jaccard.dist.fun.76.R.1, jaccard.dist.fun.76.R.2)
length(jaccard.dist.fun.76.R$Dissimilarity)
jaccard.dist.fun.76.R$Comparison <- "R"
jaccard.dist.fun.76.R$Days <- 76
mean(jaccard.dist.fun.76.R$Dissimilarity)

R.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.90.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.90 & Var2 %in% G.141)
jaccard.dist.fun.90.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.90)
jaccard.dist.fun.90.R <- rbind(jaccard.dist.fun.90.R.1, jaccard.dist.fun.90.R.2)
length(jaccard.dist.fun.90.R$Dissimilarity)
jaccard.dist.fun.90.R$Comparison <- "R"
jaccard.dist.fun.90.R$Days <- 90
mean(jaccard.dist.fun.90.R$Dissimilarity)


R.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.106.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.106 & Var2 %in% G.141)
jaccard.dist.fun.106.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.106)
jaccard.dist.fun.106.R <- rbind(jaccard.dist.fun.106.R.1, jaccard.dist.fun.106.R.2)
length(jaccard.dist.fun.106.R$Dissimilarity)
jaccard.dist.fun.106.R$Comparison <- "R"
jaccard.dist.fun.106.R$Days <- 106
mean(jaccard.dist.fun.106.R$Dissimilarity)

R.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.120.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.120 & Var2 %in% G.141)
jaccard.dist.fun.120.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.120)
jaccard.dist.fun.120.R <- rbind(jaccard.dist.fun.120.R.1, jaccard.dist.fun.120.R.2)
length(jaccard.dist.fun.120.R$Dissimilarity)
jaccard.dist.fun.120.R$Comparison <- "R"
jaccard.dist.fun.120.R$Days <- 120
mean(jaccard.dist.fun.120.R$Dissimilarity)

R.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("R"))]

jaccard.dist.fun.141.R.1 <- subset(jaccard.dist.fun.melt, Var1 %in% R.141 & Var2 %in% G.141)
jaccard.dist.fun.141.R.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% R.141)
jaccard.dist.fun.141.R <- rbind(jaccard.dist.fun.141.R.1, jaccard.dist.fun.141.R.2)
length(jaccard.dist.fun.141.R$Dissimilarity)
jaccard.dist.fun.141.R$Comparison <- "R"
jaccard.dist.fun.141.R$Days <- 141
mean(jaccard.dist.fun.141.R$Dissimilarity)


jaccard.dist.fun.w.G141 <- rbind(jaccard.dist.fun.48.L, jaccard.dist.fun.62.L,jaccard.dist.fun.76.L,
                               jaccard.dist.fun.90.L,jaccard.dist.fun.106.L,jaccard.dist.fun.120.L,
                               jaccard.dist.fun.141.L,jaccard.dist.fun.48.S,jaccard.dist.fun.62.S,
                               jaccard.dist.fun.76.S,jaccard.dist.fun.90.S,jaccard.dist.fun.106.S,
                               jaccard.dist.fun.120.S,jaccard.dist.fun.141.S,jaccard.dist.fun.48.R,
                               jaccard.dist.fun.62.R,jaccard.dist.fun.76.R,jaccard.dist.fun.90.R,
                               jaccard.dist.fun.106.R,jaccard.dist.fun.120.R,jaccard.dist.fun.141.R)


###Bacteria

G.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("G"))]



## 48
L1.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("L1"))]
L2.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("L2"))]

jaccard.dist.bac.48.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.48 & Var2 %in% G.141)
jaccard.dist.bac.48.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.48)
jaccard.dist.bac.48.L1 <- rbind(jaccard.dist.bac.48.L1.1, jaccard.dist.bac.48.L1.2)
length(jaccard.dist.bac.48.L1$Dissimilarity)
jaccard.dist.bac.48.L1$Comparison <- "L1"
mean(jaccard.dist.bac.48.L1$Dissimilarity)

jaccard.dist.bac.48.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.48 & Var2 %in% G.141)
jaccard.dist.bac.48.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.48)
jaccard.dist.bac.48.L2 <- rbind(jaccard.dist.bac.48.L2.1, jaccard.dist.bac.48.L2.2)
length(jaccard.dist.bac.48.L2$Dissimilarity)
jaccard.dist.bac.48.L2$Comparison <- "L2"
mean(jaccard.dist.bac.48.L2$Dissimilarity)


jaccard.dist.bac.48.L <- rbind(jaccard.dist.bac.48.L1,jaccard.dist.bac.48.L2)
jaccard.dist.bac.48.L$Days <- 48

##62
L1.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L1"))]
L2.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L2"))]
L3.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("L3"))]

jaccard.dist.bac.62.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.62 & Var2 %in% G.141)
jaccard.dist.bac.62.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.62)
jaccard.dist.bac.62.L1 <- rbind(jaccard.dist.bac.62.L1.1, jaccard.dist.bac.62.L1.2)
length(jaccard.dist.bac.62.L1$Dissimilarity)
jaccard.dist.bac.62.L1$Comparison <- "L1"
mean(jaccard.dist.bac.62.L1$Dissimilarity)

jaccard.dist.bac.62.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.62 & Var2 %in% G.141)
jaccard.dist.bac.62.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.62)
jaccard.dist.bac.62.L2 <- rbind(jaccard.dist.bac.62.L2.1, jaccard.dist.bac.62.L2.2)
length(jaccard.dist.bac.62.L2$Dissimilarity)
jaccard.dist.bac.62.L2$Comparison <- "L2"
mean(jaccard.dist.bac.62.L2$Dissimilarity)

jaccard.dist.bac.62.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.62 & Var2 %in% G.141)
jaccard.dist.bac.62.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.62)
jaccard.dist.bac.62.L3 <- rbind(jaccard.dist.bac.62.L3.1, jaccard.dist.bac.62.L3.2)
length(jaccard.dist.bac.62.L3$Dissimilarity)
jaccard.dist.bac.62.L3$Comparison <- "L3"
mean(jaccard.dist.bac.62.L3$Dissimilarity)

jaccard.dist.bac.62.L <- rbind(jaccard.dist.bac.62.L1,jaccard.dist.bac.62.L2,jaccard.dist.bac.62.L3)
jaccard.dist.bac.62.L$Days <- 62


#76
L1.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L1"))]
L2.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L2"))]
L3.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("L3"))]
FL.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.76.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.76 & Var2 %in% G.141)
jaccard.dist.bac.76.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.76)
jaccard.dist.bac.76.L1 <- rbind(jaccard.dist.bac.76.L1.1, jaccard.dist.bac.76.L1.2)
length(jaccard.dist.bac.76.L1$Dissimilarity)
jaccard.dist.bac.76.L1$Comparison <- "L1"
mean(jaccard.dist.bac.76.L1$Dissimilarity)

jaccard.dist.bac.76.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.76 & Var2 %in% G.141)
jaccard.dist.bac.76.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.76)
jaccard.dist.bac.76.L2 <- rbind(jaccard.dist.bac.76.L2.1, jaccard.dist.bac.76.L2.2)
length(jaccard.dist.bac.76.L2$Dissimilarity)
jaccard.dist.bac.76.L2$Comparison <- "L2"
mean(jaccard.dist.bac.76.L2$Dissimilarity)

jaccard.dist.bac.76.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.76 & Var2 %in% G.141)
jaccard.dist.bac.76.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.76)
jaccard.dist.bac.76.L3 <- rbind(jaccard.dist.bac.76.L3.1, jaccard.dist.bac.76.L3.2)
length(jaccard.dist.bac.76.L3$Dissimilarity)
jaccard.dist.bac.76.L3$Comparison <- "L3"
mean(jaccard.dist.bac.76.L3$Dissimilarity)

jaccard.dist.bac.76.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.76 & Var2 %in% G.141)
jaccard.dist.bac.76.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% FL.76)
jaccard.dist.bac.76.FL <- rbind(jaccard.dist.bac.76.FL.1, jaccard.dist.bac.76.FL.2)
length(jaccard.dist.bac.76.FL$Dissimilarity)
jaccard.dist.bac.76.FL$Comparison <- "FL"
mean(jaccard.dist.bac.76.FL$Dissimilarity)

jaccard.dist.bac.76.L <- rbind(jaccard.dist.bac.76.L1,jaccard.dist.bac.76.L2,jaccard.dist.bac.76.L3,
                               jaccard.dist.bac.76.FL)
jaccard.dist.bac.76.L$Days <- 76

###90 days
L1.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L1"))]
L2.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L2"))]
L3.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("L3"))]
FL.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.90.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.90 & Var2 %in% G.141)
jaccard.dist.bac.90.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.90)
jaccard.dist.bac.90.L1 <- rbind(jaccard.dist.bac.90.L1.1, jaccard.dist.bac.90.L1.2)
length(jaccard.dist.bac.90.L1$Dissimilarity)
jaccard.dist.bac.90.L1$Comparison <- "L1"
mean(jaccard.dist.bac.90.L1$Dissimilarity)

jaccard.dist.bac.90.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.90 & Var2 %in% G.141)
jaccard.dist.bac.90.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.90)
jaccard.dist.bac.90.L2 <- rbind(jaccard.dist.bac.90.L2.1, jaccard.dist.bac.90.L2.2)
length(jaccard.dist.bac.90.L2$Dissimilarity)
jaccard.dist.bac.90.L2$Comparison <- "L2"
mean(jaccard.dist.bac.90.L2$Dissimilarity)

jaccard.dist.bac.90.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.90 & Var2 %in% G.141)
jaccard.dist.bac.90.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.90)
jaccard.dist.bac.90.L3 <- rbind(jaccard.dist.bac.90.L3.1, jaccard.dist.bac.90.L3.2)
length(jaccard.dist.bac.90.L3$Dissimilarity)
jaccard.dist.bac.90.L3$Comparison <- "L3"
mean(jaccard.dist.bac.90.L3$Dissimilarity)

jaccard.dist.bac.90.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.90 & Var2 %in% G.141)
jaccard.dist.bac.90.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% FL.90)
jaccard.dist.bac.90.FL <- rbind(jaccard.dist.bac.90.FL.1, jaccard.dist.bac.90.FL.2)
length(jaccard.dist.bac.90.FL$Dissimilarity)
jaccard.dist.bac.90.FL$Comparison <- "FL"
mean(jaccard.dist.bac.90.FL$Dissimilarity)


jaccard.dist.bac.90.L <- rbind(jaccard.dist.bac.90.L1,jaccard.dist.bac.90.L2,jaccard.dist.bac.90.L3,
                               jaccard.dist.bac.90.FL)

jaccard.dist.bac.90.L$Days <- 90

###106 days
L1.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L1"))]
L2.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L2"))]
L3.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("L3"))]
FL.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("FL"))]


jaccard.dist.bac.106.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.106 & Var2 %in% G.141)
jaccard.dist.bac.106.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.106)
jaccard.dist.bac.106.L1 <- rbind(jaccard.dist.bac.106.L1.1, jaccard.dist.bac.106.L1.2)
length(jaccard.dist.bac.106.L1$Dissimilarity)
jaccard.dist.bac.106.L1$Comparison <- "L1"
mean(jaccard.dist.bac.106.L1$Dissimilarity)

jaccard.dist.bac.106.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.106 & Var2 %in% G.141)
jaccard.dist.bac.106.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.106)
jaccard.dist.bac.106.L2 <- rbind(jaccard.dist.bac.106.L2.1, jaccard.dist.bac.106.L2.2)
length(jaccard.dist.bac.106.L2$Dissimilarity)
jaccard.dist.bac.106.L2$Comparison <- "L2"
mean(jaccard.dist.bac.106.L2$Dissimilarity)

jaccard.dist.bac.106.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.106 & Var2 %in% G.141)
jaccard.dist.bac.106.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.106)
jaccard.dist.bac.106.L3 <- rbind(jaccard.dist.bac.106.L3.1, jaccard.dist.bac.106.L3.2)
length(jaccard.dist.bac.106.L3$Dissimilarity)
jaccard.dist.bac.106.L3$Comparison <- "L3"
mean(jaccard.dist.bac.106.L3$Dissimilarity)

jaccard.dist.bac.106.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.106 & Var2 %in% G.141)
jaccard.dist.bac.106.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% FL.106)
jaccard.dist.bac.106.FL <- rbind(jaccard.dist.bac.106.FL.1, jaccard.dist.bac.106.FL.2)
length(jaccard.dist.bac.106.FL$Dissimilarity)
jaccard.dist.bac.106.FL$Comparison <- "FL"
mean(jaccard.dist.bac.106.FL$Dissimilarity)

jaccard.dist.bac.106.L <- rbind(jaccard.dist.bac.106.L1,jaccard.dist.bac.106.L2,jaccard.dist.bac.106.L3,
                                jaccard.dist.bac.106.FL)
jaccard.dist.bac.106.L$Days <- 106
###120 days
L1.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L1"))]
L2.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L2"))]
L3.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("L3"))]
FL.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.120.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.120 & Var2 %in% G.141)
jaccard.dist.bac.120.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.120)
jaccard.dist.bac.120.L1 <- rbind(jaccard.dist.bac.120.L1.1, jaccard.dist.bac.120.L1.2)
length(jaccard.dist.bac.120.L1$Dissimilarity)
jaccard.dist.bac.120.L1$Comparison <- "L1"
mean(jaccard.dist.bac.120.L1$Dissimilarity)

jaccard.dist.bac.120.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.120 & Var2 %in% G.141)
jaccard.dist.bac.120.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.120)
jaccard.dist.bac.120.L2 <- rbind(jaccard.dist.bac.120.L2.1, jaccard.dist.bac.120.L2.2)
length(jaccard.dist.bac.120.L2$Dissimilarity)
jaccard.dist.bac.120.L2$Comparison <- "L2"
mean(jaccard.dist.bac.120.L2$Dissimilarity)

jaccard.dist.bac.120.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.120 & Var2 %in% G.141)
jaccard.dist.bac.120.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.120)
jaccard.dist.bac.120.L3 <- rbind(jaccard.dist.bac.120.L3.1, jaccard.dist.bac.120.L3.2)
length(jaccard.dist.bac.120.L3$Dissimilarity)
jaccard.dist.bac.120.L3$Comparison <- "L3"
mean(jaccard.dist.bac.120.L3$Dissimilarity)

jaccard.dist.bac.120.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.120 & Var2 %in% G.141)
jaccard.dist.bac.120.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% FL.120)
jaccard.dist.bac.120.FL <- rbind(jaccard.dist.bac.120.FL.1, jaccard.dist.bac.120.FL.2)
length(jaccard.dist.bac.120.FL$Dissimilarity)
jaccard.dist.bac.120.FL$Comparison <- "FL"
mean(jaccard.dist.bac.120.FL$Dissimilarity)

jaccard.dist.bac.120.L <- rbind(jaccard.dist.bac.120.L1,jaccard.dist.bac.120.L2,jaccard.dist.bac.120.L3,
                                jaccard.dist.bac.120.FL)

jaccard.dist.bac.120.L$Days <- 120
###141 days
L1.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L1"))]
L2.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L2"))]
L3.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("L3"))]
FL.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("FL"))]

jaccard.dist.bac.141.L1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L1.141 & Var2 %in% G.141)
jaccard.dist.bac.141.L1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L1.141)
jaccard.dist.bac.141.L1 <- rbind(jaccard.dist.bac.141.L1.1, jaccard.dist.bac.141.L1.2)
length(jaccard.dist.bac.141.L1$Dissimilarity)
jaccard.dist.bac.141.L1$Comparison <- "L1"
mean(jaccard.dist.bac.141.L1$Dissimilarity)

jaccard.dist.bac.141.L2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L2.141 & Var2 %in% G.141)
jaccard.dist.bac.141.L2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L2.141)
jaccard.dist.bac.141.L2 <- rbind(jaccard.dist.bac.141.L2.1, jaccard.dist.bac.141.L2.2)
length(jaccard.dist.bac.141.L2$Dissimilarity)
jaccard.dist.bac.141.L2$Comparison <- "L2"
mean(jaccard.dist.bac.141.L2$Dissimilarity)

jaccard.dist.bac.141.L3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% L3.141 & Var2 %in% G.141)
jaccard.dist.bac.141.L3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% L3.141)
jaccard.dist.bac.141.L3 <- rbind(jaccard.dist.bac.141.L3.1, jaccard.dist.bac.141.L3.2)
length(jaccard.dist.bac.141.L3$Dissimilarity)
jaccard.dist.bac.141.L3$Comparison <- "L3"
mean(jaccard.dist.bac.141.L3$Dissimilarity)

jaccard.dist.bac.141.FL.1 <- subset(jaccard.dist.bac.melt, Var1 %in% FL.141 & Var2 %in% G.141)
jaccard.dist.bac.141.FL.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% FL.141)
jaccard.dist.bac.141.FL <- rbind(jaccard.dist.bac.141.FL.1, jaccard.dist.bac.141.FL.2)
length(jaccard.dist.bac.141.FL$Dissimilarity)
jaccard.dist.bac.141.FL$Comparison <- "FL"
mean(jaccard.dist.bac.141.FL$Dissimilarity)

jaccard.dist.bac.141.L <- rbind(jaccard.dist.bac.141.L1,jaccard.dist.bac.141.L2,jaccard.dist.bac.141.L3,
                                jaccard.dist.bac.141.FL)
jaccard.dist.bac.141.L$Days <- 141
##Stem

S1.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S1"))]
S2.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S2"))]
S3.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("S3"))]

jaccard.dist.bac.48.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.48 & Var2 %in% G.141)
jaccard.dist.bac.48.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.48)
jaccard.dist.bac.48.S1 <- rbind(jaccard.dist.bac.48.S1.1, jaccard.dist.bac.48.S1.2)
length(jaccard.dist.bac.48.S1$Dissimilarity)
jaccard.dist.bac.48.S1$Comparison <- "S1"
mean(jaccard.dist.bac.48.S1$Dissimilarity)

jaccard.dist.bac.48.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.48 & Var2 %in% G.141)
jaccard.dist.bac.48.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.48)
jaccard.dist.bac.48.S2 <- rbind(jaccard.dist.bac.48.S2.1, jaccard.dist.bac.48.S2.2)
length(jaccard.dist.bac.48.S2$Dissimilarity)
jaccard.dist.bac.48.S2$Comparison <- "S2"
mean(jaccard.dist.bac.48.S2$Dissimilarity)

jaccard.dist.bac.48.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.48 & Var2 %in% G.141)
jaccard.dist.bac.48.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.48)
jaccard.dist.bac.48.S3 <- rbind(jaccard.dist.bac.48.S3.1, jaccard.dist.bac.48.S3.2)
length(jaccard.dist.bac.48.S3$Dissimilarity)
jaccard.dist.bac.48.S3$Comparison <- "S3"
mean(jaccard.dist.bac.48.S3$Dissimilarity)


jaccard.dist.bac.48.S <- rbind(jaccard.dist.bac.48.S1,jaccard.dist.bac.48.S2,jaccard.dist.bac.48.S3)
jaccard.dist.bac.48.S$Days <- 48


S1.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S1"))]
S2.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S2"))]
S3.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S3"))]
S4.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("S4"))]

jaccard.dist.bac.62.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.62 & Var2 %in% G.141)
jaccard.dist.bac.62.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.62)
jaccard.dist.bac.62.S1 <- rbind(jaccard.dist.bac.62.S1.1, jaccard.dist.bac.62.S1.2)
length(jaccard.dist.bac.62.S1$Dissimilarity)
jaccard.dist.bac.62.S1$Comparison <- "S1"
mean(jaccard.dist.bac.62.S1$Dissimilarity)

jaccard.dist.bac.62.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.62 & Var2 %in% G.141)
jaccard.dist.bac.62.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.62)
jaccard.dist.bac.62.S2 <- rbind(jaccard.dist.bac.62.S2.1, jaccard.dist.bac.62.S2.2)
length(jaccard.dist.bac.62.S2$Dissimilarity)
jaccard.dist.bac.62.S2$Comparison <- "S2"
mean(jaccard.dist.bac.62.S2$Dissimilarity)

jaccard.dist.bac.62.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.62 & Var2 %in% G.141)
jaccard.dist.bac.62.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.62)
jaccard.dist.bac.62.S3 <- rbind(jaccard.dist.bac.62.S3.1, jaccard.dist.bac.62.S3.2)
length(jaccard.dist.bac.62.S3$Dissimilarity)
jaccard.dist.bac.62.S3$Comparison <- "S3"
mean(jaccard.dist.bac.62.S3$Dissimilarity)

jaccard.dist.bac.62.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.62 & Var2 %in% G.141)
jaccard.dist.bac.62.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.62)
jaccard.dist.bac.62.S4 <- rbind(jaccard.dist.bac.62.S4.1, jaccard.dist.bac.62.S4.2)
length(jaccard.dist.bac.62.S4$Dissimilarity)
jaccard.dist.bac.62.S4$Comparison <- "S4"
mean(jaccard.dist.bac.62.S4$Dissimilarity)

jaccard.dist.bac.62.S <- rbind(jaccard.dist.bac.62.S1,jaccard.dist.bac.62.S2,jaccard.dist.bac.62.S3,
                               jaccard.dist.bac.62.S4)
jaccard.dist.bac.62.S$Days <- 62



S1.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S1"))]
S2.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S2"))]
S3.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S3"))]
S4.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S4"))]
S5.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S5"))]
S6.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("S6"))]

jaccard.dist.bac.76.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.76)
jaccard.dist.bac.76.S1 <- rbind(jaccard.dist.bac.76.S1.1, jaccard.dist.bac.76.S1.2)
length(jaccard.dist.bac.76.S1$Dissimilarity)
jaccard.dist.bac.76.S1$Comparison <- "S1"
mean(jaccard.dist.bac.76.S1$Dissimilarity)

jaccard.dist.bac.76.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.76)
jaccard.dist.bac.76.S2 <- rbind(jaccard.dist.bac.76.S2.1, jaccard.dist.bac.76.S2.2)
length(jaccard.dist.bac.76.S2$Dissimilarity)
jaccard.dist.bac.76.S2$Comparison <- "S2"
mean(jaccard.dist.bac.76.S2$Dissimilarity)

jaccard.dist.bac.76.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.76)
jaccard.dist.bac.76.S3 <- rbind(jaccard.dist.bac.76.S3.1, jaccard.dist.bac.76.S3.2)
length(jaccard.dist.bac.76.S3$Dissimilarity)
jaccard.dist.bac.76.S3$Comparison <- "S3"
mean(jaccard.dist.bac.76.S3$Dissimilarity)

jaccard.dist.bac.76.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.76)
jaccard.dist.bac.76.S4 <- rbind(jaccard.dist.bac.76.S4.1, jaccard.dist.bac.76.S4.2)
length(jaccard.dist.bac.76.S4$Dissimilarity)
jaccard.dist.bac.76.S4$Comparison <- "S4"
mean(jaccard.dist.bac.76.S4$Dissimilarity)

jaccard.dist.bac.76.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S5.76)
jaccard.dist.bac.76.S5 <- rbind(jaccard.dist.bac.76.S5.1, jaccard.dist.bac.76.S5.2)
length(jaccard.dist.bac.76.S5$Dissimilarity)
jaccard.dist.bac.76.S5$Comparison <- "S5"
mean(jaccard.dist.bac.76.S5$Dissimilarity)

jaccard.dist.bac.76.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.76 & Var2 %in% G.141)
jaccard.dist.bac.76.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S6.76)
jaccard.dist.bac.76.S6 <- rbind(jaccard.dist.bac.76.S6.1, jaccard.dist.bac.76.S6.2)
length(jaccard.dist.bac.76.S6$Dissimilarity)
jaccard.dist.bac.76.S6$Comparison <- "S6"
mean(jaccard.dist.bac.76.S6$Dissimilarity)

jaccard.dist.bac.76.S <- rbind(jaccard.dist.bac.76.S1,jaccard.dist.bac.76.S2,jaccard.dist.bac.76.S3,
                               jaccard.dist.bac.76.S4,jaccard.dist.bac.76.S5,jaccard.dist.bac.76.S6)
jaccard.dist.bac.76.S$Days <- 76
###90 days
S1.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S1"))]
S2.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S2"))]
S3.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S3"))]
S4.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S4"))]
S5.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S5"))]
S6.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S6"))]
S7.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S7"))]
S8.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S8"))]
S9.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.90.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.90)
jaccard.dist.bac.90.S1 <- rbind(jaccard.dist.bac.90.S1.1, jaccard.dist.bac.90.S1.2)
length(jaccard.dist.bac.90.S1$Dissimilarity)
jaccard.dist.bac.90.S1$Comparison <- "S1"
mean(jaccard.dist.bac.90.S1$Dissimilarity)

jaccard.dist.bac.90.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.90)
jaccard.dist.bac.90.S2 <- rbind(jaccard.dist.bac.90.S2.1, jaccard.dist.bac.90.S2.2)
length(jaccard.dist.bac.90.S2$Dissimilarity)
jaccard.dist.bac.90.S2$Comparison <- "S2"
mean(jaccard.dist.bac.90.S2$Dissimilarity)

jaccard.dist.bac.90.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.90)
jaccard.dist.bac.90.S3 <- rbind(jaccard.dist.bac.90.S3.1, jaccard.dist.bac.90.S3.2)
length(jaccard.dist.bac.90.S3$Dissimilarity)
jaccard.dist.bac.90.S3$Comparison <- "S3"
mean(jaccard.dist.bac.90.S3$Dissimilarity)

jaccard.dist.bac.90.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.90)
jaccard.dist.bac.90.S4 <- rbind(jaccard.dist.bac.90.S4.1, jaccard.dist.bac.90.S4.2)
length(jaccard.dist.bac.90.S4$Dissimilarity)
jaccard.dist.bac.90.S4$Comparison <- "S4"
mean(jaccard.dist.bac.90.S4$Dissimilarity)

jaccard.dist.bac.90.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S5.90)
jaccard.dist.bac.90.S5 <- rbind(jaccard.dist.bac.90.S5.1, jaccard.dist.bac.90.S5.2)
length(jaccard.dist.bac.90.S5$Dissimilarity)
jaccard.dist.bac.90.S5$Comparison <- "S5"
mean(jaccard.dist.bac.90.S5$Dissimilarity)

jaccard.dist.bac.90.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S6.90)
jaccard.dist.bac.90.S6 <- rbind(jaccard.dist.bac.90.S6.1, jaccard.dist.bac.90.S6.2)
length(jaccard.dist.bac.90.S6$Dissimilarity)
jaccard.dist.bac.90.S6$Comparison <- "S6"
mean(jaccard.dist.bac.90.S6$Dissimilarity)

jaccard.dist.bac.90.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S7.90)
jaccard.dist.bac.90.S7 <- rbind(jaccard.dist.bac.90.S7.1, jaccard.dist.bac.90.S7.2)
length(jaccard.dist.bac.90.S7$Dissimilarity)
jaccard.dist.bac.90.S7$Comparison <- "S7"
mean(jaccard.dist.bac.90.S7$Dissimilarity)

jaccard.dist.bac.90.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S8.90)
jaccard.dist.bac.90.S8 <- rbind(jaccard.dist.bac.90.S8.1, jaccard.dist.bac.90.S8.2)
length(jaccard.dist.bac.90.S8$Dissimilarity)
jaccard.dist.bac.90.S8$Comparison <- "S8"
mean(jaccard.dist.bac.90.S8$Dissimilarity)

jaccard.dist.bac.90.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.90 & Var2 %in% G.141)
jaccard.dist.bac.90.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S9.90)
jaccard.dist.bac.90.S9 <- rbind(jaccard.dist.bac.90.S9.1, jaccard.dist.bac.90.S9.2)
length(jaccard.dist.bac.90.S9$Dissimilarity)
jaccard.dist.bac.90.S9$Comparison <- "S9"
mean(jaccard.dist.bac.90.S9$Dissimilarity)


jaccard.dist.bac.90.S <- rbind(jaccard.dist.bac.90.S1,jaccard.dist.bac.90.S2,jaccard.dist.bac.90.S3,
                               jaccard.dist.bac.90.S4,jaccard.dist.bac.90.S5,jaccard.dist.bac.90.S6,
                               jaccard.dist.bac.90.S7,jaccard.dist.bac.90.S8,jaccard.dist.bac.90.S9)

jaccard.dist.bac.90.S$Days <- 90
###106 days
S1.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S1"))]
S2.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S2"))]
S3.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S3"))]
S4.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S4"))]
S5.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S5"))]
S6.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S6"))]
S7.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S7"))]
S8.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S8"))]
S9.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.106.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.106)
jaccard.dist.bac.106.S1 <- rbind(jaccard.dist.bac.106.S1.1, jaccard.dist.bac.106.S1.2)
length(jaccard.dist.bac.106.S1$Dissimilarity)
jaccard.dist.bac.106.S1$Comparison <- "S1"
mean(jaccard.dist.bac.106.S1$Dissimilarity)

jaccard.dist.bac.106.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.106)
jaccard.dist.bac.106.S2 <- rbind(jaccard.dist.bac.106.S2.1, jaccard.dist.bac.106.S2.2)
length(jaccard.dist.bac.106.S2$Dissimilarity)
jaccard.dist.bac.106.S2$Comparison <- "S2"
mean(jaccard.dist.bac.106.S2$Dissimilarity)

jaccard.dist.bac.106.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.106)
jaccard.dist.bac.106.S3 <- rbind(jaccard.dist.bac.106.S3.1, jaccard.dist.bac.106.S3.2)
length(jaccard.dist.bac.106.S3$Dissimilarity)
jaccard.dist.bac.106.S3$Comparison <- "S3"
mean(jaccard.dist.bac.106.S3$Dissimilarity)

jaccard.dist.bac.106.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.106)
jaccard.dist.bac.106.S4 <- rbind(jaccard.dist.bac.106.S4.1, jaccard.dist.bac.106.S4.2)
length(jaccard.dist.bac.106.S4$Dissimilarity)
jaccard.dist.bac.106.S4$Comparison <- "S4"
mean(jaccard.dist.bac.106.S4$Dissimilarity)

jaccard.dist.bac.106.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S5.106)
jaccard.dist.bac.106.S5 <- rbind(jaccard.dist.bac.106.S5.1, jaccard.dist.bac.106.S5.2)
length(jaccard.dist.bac.106.S5$Dissimilarity)
jaccard.dist.bac.106.S5$Comparison <- "S5"
mean(jaccard.dist.bac.106.S5$Dissimilarity)

jaccard.dist.bac.106.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S6.106)
jaccard.dist.bac.106.S6 <- rbind(jaccard.dist.bac.106.S6.1, jaccard.dist.bac.106.S6.2)
length(jaccard.dist.bac.106.S6$Dissimilarity)
jaccard.dist.bac.106.S6$Comparison <- "S6"
mean(jaccard.dist.bac.106.S6$Dissimilarity)

jaccard.dist.bac.106.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S7.106)
jaccard.dist.bac.106.S7 <- rbind(jaccard.dist.bac.106.S7.1, jaccard.dist.bac.106.S7.2)
length(jaccard.dist.bac.106.S7$Dissimilarity)
jaccard.dist.bac.106.S7$Comparison <- "S7"
mean(jaccard.dist.bac.106.S7$Dissimilarity)

jaccard.dist.bac.106.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S8.106)
jaccard.dist.bac.106.S8 <- rbind(jaccard.dist.bac.106.S8.1, jaccard.dist.bac.106.S8.2)
length(jaccard.dist.bac.106.S8$Dissimilarity)
jaccard.dist.bac.106.S8$Comparison <- "S8"
mean(jaccard.dist.bac.106.S8$Dissimilarity)

jaccard.dist.bac.106.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.106 & Var2 %in% G.141)
jaccard.dist.bac.106.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S9.106)
jaccard.dist.bac.106.S9 <- rbind(jaccard.dist.bac.106.S9.1, jaccard.dist.bac.106.S9.2)
length(jaccard.dist.bac.106.S9$Dissimilarity)
jaccard.dist.bac.106.S9$Comparison <- "S9"
mean(jaccard.dist.bac.106.S9$Dissimilarity)


jaccard.dist.bac.106.S <- rbind(jaccard.dist.bac.106.S1,jaccard.dist.bac.106.S2,jaccard.dist.bac.106.S3,
                                jaccard.dist.bac.106.S4,jaccard.dist.bac.106.S5,jaccard.dist.bac.106.S6,
                                jaccard.dist.bac.106.S7,jaccard.dist.bac.106.S8,jaccard.dist.bac.106.S9)
jaccard.dist.bac.106.S$Days <- 106
###120 days
S1.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S1"))]
S2.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S2"))]
S3.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S3"))]
S4.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S4"))]
S5.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S5"))]
S6.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S6"))]
S7.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S7"))]
S8.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S8"))]
S9.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.120.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.120)
jaccard.dist.bac.120.S1 <- rbind(jaccard.dist.bac.120.S1.1, jaccard.dist.bac.120.S1.2)
length(jaccard.dist.bac.120.S1$Dissimilarity)
jaccard.dist.bac.120.S1$Comparison <- "S1"
mean(jaccard.dist.bac.120.S1$Dissimilarity)

jaccard.dist.bac.120.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.120)
jaccard.dist.bac.120.S2 <- rbind(jaccard.dist.bac.120.S2.1, jaccard.dist.bac.120.S2.2)
length(jaccard.dist.bac.120.S2$Dissimilarity)
jaccard.dist.bac.120.S2$Comparison <- "S2"
mean(jaccard.dist.bac.120.S2$Dissimilarity)

jaccard.dist.bac.120.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.120)
jaccard.dist.bac.120.S3 <- rbind(jaccard.dist.bac.120.S3.1, jaccard.dist.bac.120.S3.2)
length(jaccard.dist.bac.120.S3$Dissimilarity)
jaccard.dist.bac.120.S3$Comparison <- "S3"
mean(jaccard.dist.bac.120.S3$Dissimilarity)

jaccard.dist.bac.120.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.120)
jaccard.dist.bac.120.S4 <- rbind(jaccard.dist.bac.120.S4.1, jaccard.dist.bac.120.S4.2)
length(jaccard.dist.bac.120.S4$Dissimilarity)
jaccard.dist.bac.120.S4$Comparison <- "S4"
mean(jaccard.dist.bac.120.S4$Dissimilarity)

jaccard.dist.bac.120.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S5.120)
jaccard.dist.bac.120.S5 <- rbind(jaccard.dist.bac.120.S5.1, jaccard.dist.bac.120.S5.2)
length(jaccard.dist.bac.120.S5$Dissimilarity)
jaccard.dist.bac.120.S5$Comparison <- "S5"
mean(jaccard.dist.bac.120.S5$Dissimilarity)

jaccard.dist.bac.120.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S6.120)
jaccard.dist.bac.120.S6 <- rbind(jaccard.dist.bac.120.S6.1, jaccard.dist.bac.120.S6.2)
length(jaccard.dist.bac.120.S6$Dissimilarity)
jaccard.dist.bac.120.S6$Comparison <- "S6"
mean(jaccard.dist.bac.120.S6$Dissimilarity)

jaccard.dist.bac.120.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S7.120)
jaccard.dist.bac.120.S7 <- rbind(jaccard.dist.bac.120.S7.1, jaccard.dist.bac.120.S7.2)
length(jaccard.dist.bac.120.S7$Dissimilarity)
jaccard.dist.bac.120.S7$Comparison <- "S7"
mean(jaccard.dist.bac.120.S7$Dissimilarity)

jaccard.dist.bac.120.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S8.120)
jaccard.dist.bac.120.S8 <- rbind(jaccard.dist.bac.120.S8.1, jaccard.dist.bac.120.S8.2)
length(jaccard.dist.bac.120.S8$Dissimilarity)
jaccard.dist.bac.120.S8$Comparison <- "S8"
mean(jaccard.dist.bac.120.S8$Dissimilarity)

jaccard.dist.bac.120.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.120 & Var2 %in% G.141)
jaccard.dist.bac.120.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S9.120)
jaccard.dist.bac.120.S9 <- rbind(jaccard.dist.bac.120.S9.1, jaccard.dist.bac.120.S9.2)
length(jaccard.dist.bac.120.S9$Dissimilarity)
jaccard.dist.bac.120.S9$Comparison <- "S9"
mean(jaccard.dist.bac.120.S9$Dissimilarity)


jaccard.dist.bac.120.S <- rbind(jaccard.dist.bac.120.S1,jaccard.dist.bac.120.S2,jaccard.dist.bac.120.S3,
                                jaccard.dist.bac.120.S4,jaccard.dist.bac.120.S5,jaccard.dist.bac.120.S6,
                                jaccard.dist.bac.120.S7,jaccard.dist.bac.120.S8,jaccard.dist.bac.120.S9)
jaccard.dist.bac.120.S$Days <- 120
###141 days
S1.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S1"))]
S2.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S2"))]
S3.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S3"))]
S4.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S4"))]
S5.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S5"))]
S6.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S6"))]
S7.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S7"))]
S8.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S8"))]
S9.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("S9"))]

jaccard.dist.bac.141.S1.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S1.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S1.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S1.141)
jaccard.dist.bac.141.S1 <- rbind(jaccard.dist.bac.141.S1.1, jaccard.dist.bac.141.S1.2)
length(jaccard.dist.bac.141.S1$Dissimilarity)
jaccard.dist.bac.141.S1$Comparison <- "S1"
mean(jaccard.dist.bac.141.S1$Dissimilarity)

jaccard.dist.bac.141.S2.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S2.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S2.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S2.141)
jaccard.dist.bac.141.S2 <- rbind(jaccard.dist.bac.141.S2.1, jaccard.dist.bac.141.S2.2)
length(jaccard.dist.bac.141.S2$Dissimilarity)
jaccard.dist.bac.141.S2$Comparison <- "S2"
mean(jaccard.dist.bac.141.S2$Dissimilarity)

jaccard.dist.bac.141.S3.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S3.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S3.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S3.141)
jaccard.dist.bac.141.S3 <- rbind(jaccard.dist.bac.141.S3.1, jaccard.dist.bac.141.S3.2)
length(jaccard.dist.bac.141.S3$Dissimilarity)
jaccard.dist.bac.141.S3$Comparison <- "S3"
mean(jaccard.dist.bac.141.S3$Dissimilarity)

jaccard.dist.bac.141.S4.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S4.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S4.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S4.141)
jaccard.dist.bac.141.S4 <- rbind(jaccard.dist.bac.141.S4.1, jaccard.dist.bac.141.S4.2)
length(jaccard.dist.bac.141.S4$Dissimilarity)
jaccard.dist.bac.141.S4$Comparison <- "S4"
mean(jaccard.dist.bac.141.S4$Dissimilarity)

jaccard.dist.bac.141.S5.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S5.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S5.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S5.141)
jaccard.dist.bac.141.S5 <- rbind(jaccard.dist.bac.141.S5.1, jaccard.dist.bac.141.S5.2)
length(jaccard.dist.bac.141.S5$Dissimilarity)
jaccard.dist.bac.141.S5$Comparison <- "S5"
mean(jaccard.dist.bac.141.S5$Dissimilarity)

jaccard.dist.bac.141.S6.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S6.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S6.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S6.141)
jaccard.dist.bac.141.S6 <- rbind(jaccard.dist.bac.141.S6.1, jaccard.dist.bac.141.S6.2)
length(jaccard.dist.bac.141.S6$Dissimilarity)
jaccard.dist.bac.141.S6$Comparison <- "S6"
mean(jaccard.dist.bac.141.S6$Dissimilarity)

jaccard.dist.bac.141.S7.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S7.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S7.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S7.141)
jaccard.dist.bac.141.S7 <- rbind(jaccard.dist.bac.141.S7.1, jaccard.dist.bac.141.S7.2)
length(jaccard.dist.bac.141.S7$Dissimilarity)
jaccard.dist.bac.141.S7$Comparison <- "S7"
mean(jaccard.dist.bac.141.S7$Dissimilarity)

jaccard.dist.bac.141.S8.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S8.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S8.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S8.141)
jaccard.dist.bac.141.S8 <- rbind(jaccard.dist.bac.141.S8.1, jaccard.dist.bac.141.S8.2)
length(jaccard.dist.bac.141.S8$Dissimilarity)
jaccard.dist.bac.141.S8$Comparison <- "S8"
mean(jaccard.dist.bac.141.S8$Dissimilarity)

jaccard.dist.bac.141.S9.1 <- subset(jaccard.dist.bac.melt, Var1 %in% S9.141 & Var2 %in% G.141)
jaccard.dist.bac.141.S9.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% S9.141)
jaccard.dist.bac.141.S9 <- rbind(jaccard.dist.bac.141.S9.1, jaccard.dist.bac.141.S9.2)
length(jaccard.dist.bac.141.S9$Dissimilarity)
jaccard.dist.bac.141.S9$Comparison <- "S9"
mean(jaccard.dist.bac.141.S9$Dissimilarity)


jaccard.dist.bac.141.S <- rbind(jaccard.dist.bac.141.S1,jaccard.dist.bac.141.S2,jaccard.dist.bac.141.S3,
                                jaccard.dist.bac.141.S4,jaccard.dist.bac.141.S5,jaccard.dist.bac.141.S6,
                                jaccard.dist.bac.141.S7,jaccard.dist.bac.141.S8,jaccard.dist.bac.141.S9)

jaccard.dist.bac.141.S$Days <- 141


##Root
R.48<-map$SampleID[which(map$Days == 48 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.48.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.48 & Var2 %in% G.141)
jaccard.dist.bac.48.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.48)
jaccard.dist.bac.48.R <- rbind(jaccard.dist.bac.48.R.1, jaccard.dist.bac.48.R.2)
length(jaccard.dist.bac.48.R$Dissimilarity)
jaccard.dist.bac.48.R$Comparison <- "R"
jaccard.dist.bac.48.R$Days <- 48
mean(jaccard.dist.bac.48.R$Dissimilarity)

R.62<-map$SampleID[which(map$Days == 62 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.62.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.62 & Var2 %in% G.141)
jaccard.dist.bac.62.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.62)
jaccard.dist.bac.62.R <- rbind(jaccard.dist.bac.62.R.1, jaccard.dist.bac.62.R.2)
length(jaccard.dist.bac.62.R$Dissimilarity)
jaccard.dist.bac.62.R$Comparison <- "R"
jaccard.dist.bac.62.R$Days <- 62
mean(jaccard.dist.bac.62.R$Dissimilarity)


R.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.76.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.76 & Var2 %in% G.141)
jaccard.dist.bac.76.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.76)
jaccard.dist.bac.76.R <- rbind(jaccard.dist.bac.76.R.1, jaccard.dist.bac.76.R.2)
length(jaccard.dist.bac.76.R$Dissimilarity)
jaccard.dist.bac.76.R$Comparison <- "R"
jaccard.dist.bac.76.R$Days <- 76
mean(jaccard.dist.bac.76.R$Dissimilarity)

R.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.90.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.90 & Var2 %in% G.141)
jaccard.dist.bac.90.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.90)
jaccard.dist.bac.90.R <- rbind(jaccard.dist.bac.90.R.1, jaccard.dist.bac.90.R.2)
length(jaccard.dist.bac.90.R$Dissimilarity)
jaccard.dist.bac.90.R$Comparison <- "R"
jaccard.dist.bac.90.R$Days <- 90
mean(jaccard.dist.bac.90.R$Dissimilarity)


R.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.106.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.106 & Var2 %in% G.141)
jaccard.dist.bac.106.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.106)
jaccard.dist.bac.106.R <- rbind(jaccard.dist.bac.106.R.1, jaccard.dist.bac.106.R.2)
length(jaccard.dist.bac.106.R$Dissimilarity)
jaccard.dist.bac.106.R$Comparison <- "R"
jaccard.dist.bac.106.R$Days <- 106
mean(jaccard.dist.bac.106.R$Dissimilarity)

R.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.120.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.120 & Var2 %in% G.141)
jaccard.dist.bac.120.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.120)
jaccard.dist.bac.120.R <- rbind(jaccard.dist.bac.120.R.1, jaccard.dist.bac.120.R.2)
length(jaccard.dist.bac.120.R$Dissimilarity)
jaccard.dist.bac.120.R$Comparison <- "R"
jaccard.dist.bac.120.R$Days <- 120
mean(jaccard.dist.bac.120.R$Dissimilarity)

R.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("R"))]

jaccard.dist.bac.141.R.1 <- subset(jaccard.dist.bac.melt, Var1 %in% R.141 & Var2 %in% G.141)
jaccard.dist.bac.141.R.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% R.141)
jaccard.dist.bac.141.R <- rbind(jaccard.dist.bac.141.R.1, jaccard.dist.bac.141.R.2)
length(jaccard.dist.bac.141.R$Dissimilarity)
jaccard.dist.bac.141.R$Comparison <- "R"
jaccard.dist.bac.141.R$Days <- 141
mean(jaccard.dist.bac.141.R$Dissimilarity)


jaccard.dist.bac.w.G141 <- rbind(jaccard.dist.bac.48.L, jaccard.dist.bac.62.L,jaccard.dist.bac.76.L,
                               jaccard.dist.bac.90.L,jaccard.dist.bac.106.L,jaccard.dist.bac.120.L,
                               jaccard.dist.bac.141.L,jaccard.dist.bac.48.S,jaccard.dist.bac.62.S,
                               jaccard.dist.bac.76.S,jaccard.dist.bac.90.S,jaccard.dist.bac.106.S,
                               jaccard.dist.bac.120.S,jaccard.dist.bac.141.S,jaccard.dist.bac.48.R,
                               jaccard.dist.bac.62.R,jaccard.dist.bac.76.R,jaccard.dist.bac.90.R,
                               jaccard.dist.bac.106.R,jaccard.dist.bac.120.R,jaccard.dist.bac.141.R)



write.csv(jaccard.dist.bac.w.G141,"jaccard.dist.bac.w.G141.csv")
write.csv(jaccard.dist.fun.w.G141,"jaccard.dist.fun.w.G141.csv")


jaccard.dist.bac.w.G0 <- read.csv("jaccard.dist.bac.w.G0.csv")
jaccard.dist.bac.w.G141 <- read.csv("jaccard.dist.bac.w.G141.csv")


jaccard.dist.fun.w.G0 <- read.csv("jaccard.dist.fun.w.G0.csv")
jaccard.dist.fun.w.G141 <- read.csv("jaccard.dist.fun.w.G141.csv")

jaccard.dist.bac.w.G0$Base <- "G0"
jaccard.dist.bac.w.G141$Base <- "G141"

jaccard.dist.fun.w.G0$Base <- "G0"
jaccard.dist.fun.w.G141$Base <- "G141"


### Merge G0 and G141
# jaccard.dist.bac <- rbind(jaccard.dist.bac.w.G0,jaccard.dist.bac.w.G141)
# jaccard.dist.fun <- rbind(jaccard.dist.fun.w.G0,jaccard.dist.fun.w.G141)


jaccard.dist.bac.w.G0$Kingdom <- "Bacteria"
jaccard.dist.fun.w.G0$Kingdom <- "Fungi"

jaccard.dist.w.G0 <- rbind(jaccard.dist.bac.w.G0,jaccard.dist.fun.w.G0)

### plotting
library(ggpubr)
p <- ggboxplot(data =jaccard.dist.w.G0, x="Days", y="Dissimilarity", fill = "Kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")+
  facet_wrap(~Comparison)
  # +stat_compare_means(aes(group = Base), label = "p.signif", method = "wilcox.test")
p

jaccard.dist.bac.w.G0$Comparison <- factor(jaccard.dist.bac.w.G0$Comparison, levels = c("FL","L3","L2","L1","S9","S8","S7","S6","S5","S4","S3","S2","S1","R","RS","BS"))

jaccard.dist.bac.w.G0$Days <-as.factor(as.character(jaccard.dist.bac.w.G0$Days))
jaccard.dist.bac.w.G0$Days <- factor(jaccard.dist.bac.w.G0$Days, levels = c("48","62","76","90","106","120","141"))

p <- ggplot(data =jaccard.dist.bac.w.G0, aes(x=Days, y=Dissimilarity, fill = Comparison)) + geom_boxplot()+
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")
#+  stat_summary(data =jaccard.dist.bac.w.G0, fun="mean", color="black", shape=23, group = aes(Comparison))
p


jaccard.dist.fun.w.G0$Comparison <- factor(jaccard.dist.fun.w.G0$Comparison, levels = c("FL","L3","L2","L1","S9","S8","S7","S6","S5","S4","S3","S2","S1","R","RS","BS"))

jaccard.dist.fun.w.G0$Days <-as.factor(as.character(jaccard.dist.fun.w.G0$Days))
jaccard.dist.fun.w.G0$Days <- factor(jaccard.dist.fun.w.G0$Days, levels = c("48","62","76","90","106","120","141"))

p <- ggplot(data =jaccard.dist.fun.w.G0, aes(x=Days, y=Dissimilarity, fill = Comparison)) + geom_boxplot()+ 
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")
p



p <- ggboxplot(data = jaccard.dist.fun, x="Days", y="Dissimilarity", fill = "Base") +
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")+
  facet_wrap(~Comparison)+
  stat_compare_means(aes(group = Base), label = "p.signif", method = "wilcox.test")
p



#### Dissimilarity between G0 and developing seeds
G.0<-f.meta.18$SampleID[which(f.meta.18$Days == 0 & f.meta.18$Microhabitat %in% c("G"))]

G.76<-f.meta.18$SampleID[which(f.meta.18$Days == 76 & f.meta.18$Microhabitat %in% c("G"))]
G.90<-f.meta.18$SampleID[which(f.meta.18$Days == 90 & f.meta.18$Microhabitat %in% c("G"))]
G.106<-f.meta.18$SampleID[which(f.meta.18$Days == 106 & f.meta.18$Microhabitat %in% c("G"))]
G.120<-f.meta.18$SampleID[which(f.meta.18$Days == 120 & f.meta.18$Microhabitat %in% c("G"))]
G.141<-f.meta.18$SampleID[which(f.meta.18$Days == 141 & f.meta.18$Microhabitat %in% c("G"))]

jaccard.dist.fun.76.G.1 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.76)
jaccard.dist.fun.76.G.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.76 & Var2 %in% G.0)
jaccard.dist.fun.76.G <- rbind(jaccard.dist.fun.76.G.1, jaccard.dist.fun.76.G.2)
length(jaccard.dist.fun.76.G$Dissimilarity)
jaccard.dist.fun.76.G$Comparison <- "G"
jaccard.dist.fun.76.G$Days <- 76
mean(jaccard.dist.fun.76.G$Dissimilarity)

jaccard.dist.fun.90.G.1 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.90)
jaccard.dist.fun.90.G.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.90 & Var2 %in% G.0)
jaccard.dist.fun.90.G <- rbind(jaccard.dist.fun.90.G.1, jaccard.dist.fun.90.G.2)
length(jaccard.dist.fun.90.G$Dissimilarity)
jaccard.dist.fun.90.G$Comparison <- "G"
jaccard.dist.fun.90.G$Days <- 90
mean(jaccard.dist.fun.90.G$Dissimilarity)

jaccard.dist.fun.106.G.1 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.106)
jaccard.dist.fun.106.G.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.106 & Var2 %in% G.0)
jaccard.dist.fun.106.G <- rbind(jaccard.dist.fun.106.G.1, jaccard.dist.fun.106.G.2)
length(jaccard.dist.fun.106.G$Dissimilarity)
jaccard.dist.fun.106.G$Comparison <- "G"
jaccard.dist.fun.106.G$Days <- 106
mean(jaccard.dist.fun.106.G$Dissimilarity)

jaccard.dist.fun.120.G.1 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.120)
jaccard.dist.fun.120.G.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.120 & Var2 %in% G.0)
jaccard.dist.fun.120.G <- rbind(jaccard.dist.fun.120.G.1, jaccard.dist.fun.120.G.2)
length(jaccard.dist.fun.120.G$Dissimilarity)
jaccard.dist.fun.120.G$Comparison <- "G"
jaccard.dist.fun.120.G$Days <- 120
mean(jaccard.dist.fun.120.G$Dissimilarity)

jaccard.dist.fun.141.G.1 <- subset(jaccard.dist.fun.melt, Var1 %in% G.0 & Var2 %in% G.141)
jaccard.dist.fun.141.G.2 <- subset(jaccard.dist.fun.melt, Var1 %in% G.141 & Var2 %in% G.0)
jaccard.dist.fun.141.G <- rbind(jaccard.dist.fun.141.G.1, jaccard.dist.fun.141.G.2)
length(jaccard.dist.fun.141.G$Dissimilarity)
jaccard.dist.fun.141.G$Comparison <- "G"
jaccard.dist.fun.141.G$Days <- 141
mean(jaccard.dist.fun.141.G$Dissimilarity)




jaccard.dist.fun.G <- rbind(jaccard.dist.fun.76.G, jaccard.dist.fun.90.G,jaccard.dist.fun.106.G,jaccard.dist.fun.120.G,
                                 jaccard.dist.fun.141.G)


##Bacteria
G.0<-map$SampleID[which(map$Days == 0 & map$Microhabitat %in% c("G"))]

G.76<-map$SampleID[which(map$Days == 76 & map$Microhabitat %in% c("G"))]
G.90<-map$SampleID[which(map$Days == 90 & map$Microhabitat %in% c("G"))]
G.106<-map$SampleID[which(map$Days == 106 & map$Microhabitat %in% c("G"))]
G.120<-map$SampleID[which(map$Days == 120 & map$Microhabitat %in% c("G"))]
G.141<-map$SampleID[which(map$Days == 141 & map$Microhabitat %in% c("G"))]

jaccard.dist.bac.76.G.1 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.76)
jaccard.dist.bac.76.G.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.76 & Var2 %in% G.0)
jaccard.dist.bac.76.G <- rbind(jaccard.dist.bac.76.G.1, jaccard.dist.bac.76.G.2)
length(jaccard.dist.bac.76.G$Dissimilarity)
jaccard.dist.bac.76.G$Comparison <- "G"
jaccard.dist.bac.76.G$Days <- 76
mean(jaccard.dist.bac.76.G$Dissimilarity)

jaccard.dist.bac.90.G.1 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.90)
jaccard.dist.bac.90.G.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.90 & Var2 %in% G.0)
jaccard.dist.bac.90.G <- rbind(jaccard.dist.bac.90.G.1, jaccard.dist.bac.90.G.2)
length(jaccard.dist.bac.90.G$Dissimilarity)
jaccard.dist.bac.90.G$Comparison <- "G"
jaccard.dist.bac.90.G$Days <- 90
mean(jaccard.dist.bac.90.G$Dissimilarity)

jaccard.dist.bac.106.G.1 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.106)
jaccard.dist.bac.106.G.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.106 & Var2 %in% G.0)
jaccard.dist.bac.106.G <- rbind(jaccard.dist.bac.106.G.1, jaccard.dist.bac.106.G.2)
length(jaccard.dist.bac.106.G$Dissimilarity)
jaccard.dist.bac.106.G$Comparison <- "G"
jaccard.dist.bac.106.G$Days <- 106
mean(jaccard.dist.bac.106.G$Dissimilarity)

jaccard.dist.bac.120.G.1 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.120)
jaccard.dist.bac.120.G.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.120 & Var2 %in% G.0)
jaccard.dist.bac.120.G <- rbind(jaccard.dist.bac.120.G.1, jaccard.dist.bac.120.G.2)
length(jaccard.dist.bac.120.G$Dissimilarity)
jaccard.dist.bac.120.G$Comparison <- "G"
jaccard.dist.bac.120.G$Days <- 120
mean(jaccard.dist.bac.120.G$Dissimilarity)

jaccard.dist.bac.141.G.1 <- subset(jaccard.dist.bac.melt, Var1 %in% G.0 & Var2 %in% G.141)
jaccard.dist.bac.141.G.2 <- subset(jaccard.dist.bac.melt, Var1 %in% G.141 & Var2 %in% G.0)
jaccard.dist.bac.141.G <- rbind(jaccard.dist.bac.141.G.1, jaccard.dist.bac.141.G.2)
length(jaccard.dist.bac.141.G$Dissimilarity)
jaccard.dist.bac.141.G$Comparison <- "G"
jaccard.dist.bac.141.G$Days <- 141
mean(jaccard.dist.bac.141.G$Dissimilarity)




jaccard.dist.bac.G <- rbind(jaccard.dist.bac.76.G, jaccard.dist.bac.90.G,jaccard.dist.bac.106.G,jaccard.dist.bac.120.G,
                            jaccard.dist.bac.141.G)

jaccard.dist.bac.G$Kingdom <- "Bacteria"
jaccard.dist.fun.G$Kingdom <- "Fungi"

jaccard.dist.G <- rbind(jaccard.dist.bac.G,jaccard.dist.fun.G)

p <- ggboxplot(data = jaccard.dist.G, x="Days", y="Dissimilarity", fill = "Kingdom") +
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")
p



ggplot(jaccard.dist.G, aes(x=Days, y=Dissimilarity, color = Kingdom)) + geom_point(position = 'jitter', alpha = 0.5, aes(shape = Kingdom)) + 
  geom_smooth(method = "lm")+theme_bw() + theme(aspect.ratio=0.5)+ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")


linear.bac<-lm(Days~Dissimilarity, subset(jaccard.dist.G, Kingdom == "Bacteria"))
linear.bac
summary(linear.bac)

linear.fun<-lm(Days~Dissimilarity, subset(jaccard.dist.G, Kingdom == "Fungi"))
linear.fun
summary(linear.fun)



### Table for jaccard distance
jaccard.dist.bac.w.G0
jaccard.dist.fun.w.G0

jaccard.dist.bac.w.G0.summary <- jaccard.dist.bac.w.G0%>% group_by(Comparison, Days) %>% summarise(Mean_dissimilarity = mean(Dissimilarity), STD = sd(Dissimilarity))

jaccard.dist.fun.w.G0.summary <- jaccard.dist.fun.w.G0%>% group_by(Comparison, Days) %>% summarise(Mean_dissimilarity = mean(Dissimilarity), STD = sd(Dissimilarity))


write.csv(jaccard.dist.bac.w.G0.summary,"jaccard.dist.bac.w.G0.summary.csv")
write.csv(jaccard.dist.fun.w.G0.summary,"jaccard.dist.fun.w.G0.summary.csv")
