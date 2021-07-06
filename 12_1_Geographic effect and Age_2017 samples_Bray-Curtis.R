###
bac.clean.log.17 <- subset_samples(bac.clean.log, Year == "year_2017")
bac.clean.log.17<- phyloseq::filter_taxa(bac.clean.log.17, function(x) sum(x) != 0, TRUE)

b.otu.lognorm.17 <- otu_table(bac.clean.log.17)
bray.dist.bac.17<-vegdist(t(b.otu.lognorm.17), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
class(bray.dist.bac.17)

bray.dist.bac.17 <-as.matrix(bray.dist.bac.17)

bray.dist.bac.17.lower<-get_lower_tri(bray.dist.bac.17)
bray.dist.bac.17.melt <- reshape2::melt(as.matrix(bray.dist.bac.17.lower), na.rm = T)
head(bray.dist.bac.17.melt)
names(bray.dist.bac.17.melt)[3] <- "Dissimilarity"
bray.dist.bac.17.melt <- subset(bray.dist.bac.17.melt, Dissimilarity != 0)

### Bulk soil samples
BS.50<-b.meta.17$SampleID[which(b.meta.17$Days == 50 & b.meta.17$Compartment =="Soil")]

bray.dist.bac.17.50.BS <- subset(bray.dist.bac.17.melt, Var1 %in% BS.50 & Var2 %in% BS.50)
bray.dist.bac.17.50.BS$Comparison <- "BS_50"
mean(bray.dist.bac.17.50.BS$Dissimilarity) #0.8396446

BS.80<-b.meta.17$SampleID[which(b.meta.17$Days == 80 & b.meta.17$Compartment =="Soil")]

bray.dist.bac.17.80.BS <- subset(bray.dist.bac.17.melt, Var1 %in% BS.80 & Var2 %in% BS.80)
bray.dist.bac.17.80.BS$Comparison <- "BS_80"
mean(bray.dist.bac.17.80.BS$Dissimilarity) #0.7851222

BS.120<-b.meta.17$SampleID[which(b.meta.17$Days == 120 & b.meta.17$Compartment =="Soil")]

bray.dist.bac.17.120.BS <- subset(bray.dist.bac.17.melt, Var1 %in% BS.120 & Var2 %in% BS.120)
bray.dist.bac.17.120.BS$Comparison <- "BS_120"
mean(bray.dist.bac.17.120.BS$Dissimilarity) # 0.8045414

BS.140<-b.meta.17$SampleID[which(b.meta.17$Days == 140 & b.meta.17$Compartment =="Soil")]

bray.dist.bac.17.140.BS <- subset(bray.dist.bac.17.melt, Var1 %in% BS.140 & Var2 %in% BS.140)
bray.dist.bac.17.140.BS$Comparison <- "BS_140"
mean(bray.dist.bac.17.140.BS$Dissimilarity) #0.7914857

bray.dist.bac.17.BS <- rbind(bray.dist.bac.17.50.BS,bray.dist.bac.17.80.BS,bray.dist.bac.17.120.BS,bray.dist.bac.17.140.BS)

### Root endosphere samples
R.50<-b.meta.17$SampleID[which(b.meta.17$Days == 50 & b.meta.17$Compartment =="Root")]

bray.dist.bac.17.50.R <- subset(bray.dist.bac.17.melt, Var1 %in% R.50 & Var2 %in% R.50)
bray.dist.bac.17.50.R$Comparison <- "R_50"
mean(bray.dist.bac.17.50.R$Dissimilarity) #0.7315499

R.80<-b.meta.17$SampleID[which(b.meta.17$Days == 80 & b.meta.17$Compartment =="Root")]

bray.dist.bac.17.80.R <- subset(bray.dist.bac.17.melt, Var1 %in% R.80 & Var2 %in% R.80)
bray.dist.bac.17.80.R$Comparison <- "R_80"
mean(bray.dist.bac.17.80.R$Dissimilarity) #0.7034279

R.120<-b.meta.17$SampleID[which(b.meta.17$Days == 120 & b.meta.17$Compartment =="Root")]

bray.dist.bac.17.120.R <- subset(bray.dist.bac.17.melt, Var1 %in% R.120 & Var2 %in% R.120)
bray.dist.bac.17.120.R$Comparison <- "R_120"
mean(bray.dist.bac.17.120.R$Dissimilarity) #0.7326756

bray.dist.bac.17.R <- rbind(bray.dist.bac.17.50.R,bray.dist.bac.17.80.R,bray.dist.bac.17.120.R)


### Stem endosphere samples
S.50<-b.meta.17$SampleID[which(b.meta.17$Days == 50 & b.meta.17$Compartment =="Stem")]

bray.dist.bac.17.50.S <- subset(bray.dist.bac.17.melt, Var1 %in% S.50 & Var2 %in% S.50)
bray.dist.bac.17.50.S$Comparison <- "S_50"
mean(bray.dist.bac.17.50.S$Dissimilarity) #0.8548505

S.80<-b.meta.17$SampleID[which(b.meta.17$Days == 80 & b.meta.17$Compartment =="Stem")]

bray.dist.bac.17.80.S <- subset(bray.dist.bac.17.melt, Var1 %in% S.80 & Var2 %in% S.80)
bray.dist.bac.17.80.S$Comparison <- "S_80"
mean(bray.dist.bac.17.80.S$Dissimilarity) # 0.8138796

S.120<-b.meta.17$SampleID[which(b.meta.17$Days == 120 & b.meta.17$Compartment =="Stem")]

bray.dist.bac.17.120.S <- subset(bray.dist.bac.17.melt, Var1 %in% S.120 & Var2 %in% S.120)
bray.dist.bac.17.120.S$Comparison <- "S_120"
mean(bray.dist.bac.17.120.S$Dissimilarity) #0.6080808

bray.dist.bac.17.S <- rbind(bray.dist.bac.17.50.S,bray.dist.bac.17.80.S,bray.dist.bac.17.120.S)


### Leaf endosphere samples
L.50<-b.meta.17$SampleID[which(b.meta.17$Days == 50 & b.meta.17$Compartment =="Leaf")]

bray.dist.bac.17.50.L <- subset(bray.dist.bac.17.melt, Var1 %in% L.50 & Var2 %in% L.50)
bray.dist.bac.17.50.L$Comparison <- "L_50"
mean(bray.dist.bac.17.50.L$Dissimilarity) #0.8663293

L.80<-b.meta.17$SampleID[which(b.meta.17$Days == 80 & b.meta.17$Compartment =="Leaf")]

bray.dist.bac.17.80.L <- subset(bray.dist.bac.17.melt, Var1 %in% L.80 & Var2 %in% L.80)
bray.dist.bac.17.80.L$Comparison <- "L_80"
mean(bray.dist.bac.17.80.L$Dissimilarity) # 0.7758032

L.120<-b.meta.17$SampleID[which(b.meta.17$Days == 120 & b.meta.17$Compartment =="Leaf")]

bray.dist.bac.17.120.L <- subset(bray.dist.bac.17.melt, Var1 %in% L.120 & Var2 %in% L.120)
bray.dist.bac.17.120.L$Comparison <- "L_120"
mean(bray.dist.bac.17.120.L$Dissimilarity) #0.5874034

bray.dist.bac.17.L <- rbind(bray.dist.bac.17.50.L,bray.dist.bac.17.80.L,bray.dist.bac.17.120.L)

### Seed endosphere samples
G.140<-b.meta.17$SampleID[which(b.meta.17$Days == 140 & b.meta.17$Compartment =="Seed")]

bray.dist.bac.17.140.G <- subset(bray.dist.bac.17.melt, Var1 %in% G.140 & Var2 %in% G.140)
bray.dist.bac.17.140.G$Comparison <- "G_140"
mean(bray.dist.bac.17.140.G$Dissimilarity) #0.4767764

G.80<-b.meta.17$SampleID[which(b.meta.17$Days == 80 & b.meta.17$Compartment =="Seed")]

bray.dist.bac.17.80.G <- subset(bray.dist.bac.17.melt, Var1 %in% G.80 & Var2 %in% G.80)
bray.dist.bac.17.80.G$Comparison <- "G_80"
mean(bray.dist.bac.17.80.G$Dissimilarity) # 0.7291162

G.120<-b.meta.17$SampleID[which(b.meta.17$Days == 120 & b.meta.17$Compartment =="Seed")]

bray.dist.bac.17.120.G <- subset(bray.dist.bac.17.melt, Var1 %in% G.120 & Var2 %in% G.120)
bray.dist.bac.17.120.G$Comparison <- "G_120"
mean(bray.dist.bac.17.120.G$Dissimilarity) #0.5200744

bray.dist.bac.17.G <- rbind(bray.dist.bac.17.80.G,bray.dist.bac.17.120.G,bray.dist.bac.17.140.G)

bray.dist.bac.17.L$Compartment <- "Leaf" 
bray.dist.bac.17.S$Compartment <- "Stem" 
bray.dist.bac.17.R$Compartment <- "Root" 
bray.dist.bac.17.G$Compartment <- "Seed" 
bray.dist.bac.17.BS$Compartment <- "Soil" 

bray.dist.bac.17.all <- rbind(bray.dist.bac.17.L,bray.dist.bac.17.S,bray.dist.bac.17.R,bray.dist.bac.17.G,bray.dist.bac.17.BS)
bray.dist.bac.17.all$Comparison <- factor(bray.dist.bac.17.all$Comparison, levels = c("L_50","L_80","L_120", "S_50","S_80","S_120", 
                                                                                      "R_50","R_80","R_120", "G_80","G_120","G_140",
                                                                                      "BS_50","BS_80","BS_120", "BS_140"))

SW <- b.meta.17$SampleID[which(b.meta.17$Location =="Suwon")]
CC1 <- b.meta.17$SampleID[which(b.meta.17$Location =="Chuncheon1")]
CC2 <- b.meta.17$SampleID[which(b.meta.17$Location =="Chuncheon2")]
bray.dist.bac.17.all <- subset(bray.dist.bac.17.all, !(Var1 %in% SW & Var2 %in% SW))
bray.dist.bac.17.all <- subset(bray.dist.bac.17.all, !(Var1 %in% CC1 & Var2 %in% CC1))
bray.dist.bac.17.all <- subset(bray.dist.bac.17.all, !(Var1 %in% CC2 & Var2 %in% CC2))

write.csv(bray.dist.bac.17.all,"bray.dist.bac.17.all.csv")
pp.all <- ggplot(data=bray.dist.bac.17.all, aes(x=Comparison, y=Dissimilarity))
panel_all=pp.all+geom_boxplot(aes(fill=Compartment)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="RA of core OTUs")+theme_bw(base_size=10)+
  scale_fill_manual(values=color_compartment) + theme_new +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_all



#########
fun.clean.log.17 <- subset_samples(fun.clean.log, Year == "year_2017")
fun.clean.log.17<- phyloseq::filter_taxa(fun.clean.log.17, function(x) sum(x) != 0, TRUE)

f.otu.lognorm.17 <- otu_table(fun.clean.log.17)
bray.dist.fun.17<-vegdist(t(f.otu.lognorm.17), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = T)
class(bray.dist.fun.17)

bray.dist.fun.17 <-as.matrix(bray.dist.fun.17)

bray.dist.fun.17.lower<-get_lower_tri(bray.dist.fun.17)
bray.dist.fun.17.melt <- reshape2::melt(as.matrix(bray.dist.fun.17.lower), na.rm = T)
head(bray.dist.fun.17.melt)
names(bray.dist.fun.17.melt)[3] <- "Dissimilarity"
bray.dist.fun.17.melt <- subset(bray.dist.fun.17.melt, Dissimilarity != 0)

### Bulk soil samples
BS.50<-f.meta.17$SampleID[which(f.meta.17$Days == 50 & f.meta.17$Compartment =="Soil")]

bray.dist.fun.17.50.BS <- subset(bray.dist.fun.17.melt, Var1 %in% BS.50 & Var2 %in% BS.50)
bray.dist.fun.17.50.BS$Comparison <- "BS_50"
mean(bray.dist.fun.17.50.BS$Dissimilarity) #0.8396446

BS.80<-f.meta.17$SampleID[which(f.meta.17$Days == 80 & f.meta.17$Compartment =="Soil")]

bray.dist.fun.17.80.BS <- subset(bray.dist.fun.17.melt, Var1 %in% BS.80 & Var2 %in% BS.80)
bray.dist.fun.17.80.BS$Comparison <- "BS_80"
mean(bray.dist.fun.17.80.BS$Dissimilarity) #0.7851222

BS.120<-f.meta.17$SampleID[which(f.meta.17$Days == 120 & f.meta.17$Compartment =="Soil")]

bray.dist.fun.17.120.BS <- subset(bray.dist.fun.17.melt, Var1 %in% BS.120 & Var2 %in% BS.120)
bray.dist.fun.17.120.BS$Comparison <- "BS_120"
mean(bray.dist.fun.17.120.BS$Dissimilarity) # 0.8045414

BS.140<-f.meta.17$SampleID[which(f.meta.17$Days == 140 & f.meta.17$Compartment =="Soil")]

bray.dist.fun.17.140.BS <- subset(bray.dist.fun.17.melt, Var1 %in% BS.140 & Var2 %in% BS.140)
bray.dist.fun.17.140.BS$Comparison <- "BS_140"
mean(bray.dist.fun.17.140.BS$Dissimilarity) #0.7914857

bray.dist.fun.17.BS <- rbind(bray.dist.fun.17.50.BS,bray.dist.fun.17.80.BS,bray.dist.fun.17.120.BS,bray.dist.fun.17.140.BS)

### Root endosphere samples
R.50<-f.meta.17$SampleID[which(f.meta.17$Days == 50 & f.meta.17$Compartment =="Root")]

bray.dist.fun.17.50.R <- subset(bray.dist.fun.17.melt, Var1 %in% R.50 & Var2 %in% R.50)
bray.dist.fun.17.50.R$Comparison <- "R_50"
mean(bray.dist.fun.17.50.R$Dissimilarity) #0.7315499

R.80<-f.meta.17$SampleID[which(f.meta.17$Days == 80 & f.meta.17$Compartment =="Root")]

bray.dist.fun.17.80.R <- subset(bray.dist.fun.17.melt, Var1 %in% R.80 & Var2 %in% R.80)
bray.dist.fun.17.80.R$Comparison <- "R_80"
mean(bray.dist.fun.17.80.R$Dissimilarity) #0.7034279

R.120<-f.meta.17$SampleID[which(f.meta.17$Days == 120 & f.meta.17$Compartment =="Root")]

bray.dist.fun.17.120.R <- subset(bray.dist.fun.17.melt, Var1 %in% R.120 & Var2 %in% R.120)
bray.dist.fun.17.120.R$Comparison <- "R_120"
mean(bray.dist.fun.17.120.R$Dissimilarity) #0.7326756

bray.dist.fun.17.R <- rbind(bray.dist.fun.17.50.R,bray.dist.fun.17.80.R,bray.dist.fun.17.120.R)


### Stem endosphere samples
S.50<-f.meta.17$SampleID[which(f.meta.17$Days == 50 & f.meta.17$Compartment =="Stem")]

bray.dist.fun.17.50.S <- subset(bray.dist.fun.17.melt, Var1 %in% S.50 & Var2 %in% S.50)
bray.dist.fun.17.50.S$Comparison <- "S_50"
mean(bray.dist.fun.17.50.S$Dissimilarity) #0.8548505

S.80<-f.meta.17$SampleID[which(f.meta.17$Days == 80 & f.meta.17$Compartment =="Stem")]

bray.dist.fun.17.80.S <- subset(bray.dist.fun.17.melt, Var1 %in% S.80 & Var2 %in% S.80)
bray.dist.fun.17.80.S$Comparison <- "S_80"
mean(bray.dist.fun.17.80.S$Dissimilarity) # 0.8138796

S.120<-f.meta.17$SampleID[which(f.meta.17$Days == 120 & f.meta.17$Compartment =="Stem")]

bray.dist.fun.17.120.S <- subset(bray.dist.fun.17.melt, Var1 %in% S.120 & Var2 %in% S.120)
bray.dist.fun.17.120.S$Comparison <- "S_120"
mean(bray.dist.fun.17.120.S$Dissimilarity) #0.6080808

bray.dist.fun.17.S <- rbind(bray.dist.fun.17.50.S,bray.dist.fun.17.80.S,bray.dist.fun.17.120.S)


### Leaf endosphere samples
L.50<-f.meta.17$SampleID[which(f.meta.17$Days == 50 & f.meta.17$Compartment =="Leaf")]

bray.dist.fun.17.50.L <- subset(bray.dist.fun.17.melt, Var1 %in% L.50 & Var2 %in% L.50)
bray.dist.fun.17.50.L$Comparison <- "L_50"
mean(bray.dist.fun.17.50.L$Dissimilarity) #0.8663293

L.80<-f.meta.17$SampleID[which(f.meta.17$Days == 80 & f.meta.17$Compartment =="Leaf")]

bray.dist.fun.17.80.L <- subset(bray.dist.fun.17.melt, Var1 %in% L.80 & Var2 %in% L.80)
bray.dist.fun.17.80.L$Comparison <- "L_80"
mean(bray.dist.fun.17.80.L$Dissimilarity) # 0.7758032

L.120<-f.meta.17$SampleID[which(f.meta.17$Days == 120 & f.meta.17$Compartment =="Leaf")]

bray.dist.fun.17.120.L <- subset(bray.dist.fun.17.melt, Var1 %in% L.120 & Var2 %in% L.120)
bray.dist.fun.17.120.L$Comparison <- "L_120"
mean(bray.dist.fun.17.120.L$Dissimilarity) #0.5874034

bray.dist.fun.17.L <- rbind(bray.dist.fun.17.50.L,bray.dist.fun.17.80.L,bray.dist.fun.17.120.L)

### Seed endosphere samples
G.140<-f.meta.17$SampleID[which(f.meta.17$Days == 140 & f.meta.17$Compartment =="Seed")]

bray.dist.fun.17.140.G <- subset(bray.dist.fun.17.melt, Var1 %in% G.140 & Var2 %in% G.140)
bray.dist.fun.17.140.G$Comparison <- "G_140"
mean(bray.dist.fun.17.140.G$Dissimilarity) #0.4767764

G.80<-f.meta.17$SampleID[which(f.meta.17$Days == 80 & f.meta.17$Compartment =="Seed")]

bray.dist.fun.17.80.G <- subset(bray.dist.fun.17.melt, Var1 %in% G.80 & Var2 %in% G.80)
bray.dist.fun.17.80.G$Comparison <- "G_80"
mean(bray.dist.fun.17.80.G$Dissimilarity) # 0.7291162

G.120<-f.meta.17$SampleID[which(f.meta.17$Days == 120 & f.meta.17$Compartment =="Seed")]

bray.dist.fun.17.120.G <- subset(bray.dist.fun.17.melt, Var1 %in% G.120 & Var2 %in% G.120)
bray.dist.fun.17.120.G$Comparison <- "G_120"
mean(bray.dist.fun.17.120.G$Dissimilarity) #0.5200744

bray.dist.fun.17.G <- rbind(bray.dist.fun.17.80.G,bray.dist.fun.17.120.G,bray.dist.fun.17.140.G)

bray.dist.fun.17.L$Compartment <- "Leaf" 
bray.dist.fun.17.S$Compartment <- "Stem" 
bray.dist.fun.17.R$Compartment <- "Root" 
bray.dist.fun.17.G$Compartment <- "Seed" 
bray.dist.fun.17.BS$Compartment <- "Soil" 

bray.dist.fun.17.all <- rbind(bray.dist.fun.17.L,bray.dist.fun.17.S,bray.dist.fun.17.R,bray.dist.fun.17.G,bray.dist.fun.17.BS)
bray.dist.fun.17.all$Comparison <- factor(bray.dist.fun.17.all$Comparison, levels = c("L_50","L_80","L_120", "S_50","S_80","S_120", 
                                                                                      "R_50","R_80","R_120", "G_80","G_120","G_140",
                                                                                      "BS_50","BS_80","BS_120", "BS_140"))

SW <- f.meta.17$SampleID[which(f.meta.17$Location =="Suwon")]
CC1 <- f.meta.17$SampleID[which(f.meta.17$Location =="Chuncheon1")]
CC2 <- f.meta.17$SampleID[which(f.meta.17$Location =="Chuncheon2")]
bray.dist.fun.17.all <- subset(bray.dist.fun.17.all, !(Var1 %in% SW & Var2 %in% SW))
bray.dist.fun.17.all <- subset(bray.dist.fun.17.all, !(Var1 %in% CC1 & Var2 %in% CC1))
bray.dist.fun.17.all <- subset(bray.dist.fun.17.all, !(Var1 %in% CC2 & Var2 %in% CC2))
write.csv(bray.dist.fun.17.all,"bray.dist.fun.17.all.csv")

pp.all <- ggplot(data=bray.dist.fun.17.all, aes(x=Comparison, y=Dissimilarity))
panel_all=pp.all+geom_boxplot(aes(fill=Compartment)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="RA of core OTUs")+theme_bw(base_size=10)+
  scale_fill_manual(values=color_compartment) + theme_new +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_all


###Anova
shapiro.test(bray.dist.fun.17.all$Dissimilarity)
bray.dist.fun.17.all.L <- subset(bray.dist.fun.17.all, Compartment == "Leaf")

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.L)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.L,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn

bray.dist.fun.17.all.S <- subset(bray.dist.fun.17.all, Compartment == "Stem")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.S)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.S,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn


bray.dist.fun.17.all.R <- subset(bray.dist.fun.17.all, Compartment == "Root")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.R)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.R,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn


bray.dist.fun.17.all.G <- subset(bray.dist.fun.17.all, Compartment == "Seed")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.G)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.G,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn

bray.dist.fun.17.all.BS <- subset(bray.dist.fun.17.all, Compartment == "Soil")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.BS)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.fun.17.all.BS,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

##Bateria
shapiro.test(bray.dist.bac.17.all$Dissimilarity)
bray.dist.bac.17.all.L <- subset(bray.dist.bac.17.all, Compartment == "Leaf")

kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.L)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.L,
              method="bh")
PT = DT$res


#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn

bray.dist.bac.17.all.S <- subset(bray.dist.bac.17.all, Compartment == "Stem")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.S)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.S,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn

bray.dist.bac.17.all.R <- subset(bray.dist.bac.17.all, Compartment == "Root")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.R)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.R,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn


bray.dist.bac.17.all.G <- subset(bray.dist.bac.17.all, Compartment == "Seed")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.G)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.G,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn

bray.dist.bac.17.all.BS <- subset(bray.dist.bac.17.all, Compartment == "Soil")
kw<-kruskal.test(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.BS)
kw$p.value
kw
#library(FSA)
DT = dunnTest(Dissimilarity ~ Comparison, data = bray.dist.bac.17.all.BS,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)

dunn
#### Core members
core.bac.seed.17 #56
core.bac.leaf.17 #65
core.bac.stem.17 #146

above.core.17<-Reduce(intersect, list(core.bac.seed.17, core.bac.leaf.17, core.bac.stem.17)) #23

bac.clean.rel.17 <- subset_samples(bac.clean.ss.rel,Year == "year_2017")
bac.clean.rel.17<- phyloseq::filter_taxa(bac.clean.rel.17, function(x) sum(x) != 0, TRUE)


bac.clean.rel.17.core <- prune_taxa(taxa_names(bac.clean.rel.17)%in% above.core.17,bac.clean.rel.17)
above.core.sum.17<-colSums(otu_table(bac.clean.rel.17.core))
above.core.sum.17 <- data.frame(above.core.sum.17)
names(above.core.sum.17)[1] <- "Core_abund"
above.core.sum.17$SampleID <- rownames(above.core.sum.17)


b.meta.17.2 <- b.meta.17
b.meta.17$Comp_age <- paste0(b.meta.17$Compartment,"_",b.meta.17$Days)
b.meta.17 <- data.frame(b.meta.17)
above.core.meta<-merge(b.meta.17,above.core.sum.17, by = "SampleID")

above.core.meta$Comp_age <- factor(above.core.meta$Comp_age, levels = c("Leaf_50","Leaf_80","Leaf_120", "Stem_50","Stem_80","Stem_120", 
                                                                        "Root_50","Root_80","Root_120", "Seed_80","Seed_120","Seed_140",
                                                                        "Soil_50","Soil_80","Soil_120", "Soil_140"))


pp.core <- ggplot(data=above.core.meta, aes(x=Comp_age, y=Core_abund))
panel_core=pp.core+geom_boxplot(aes(fill=Compartment)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="RA of core OTUs")+theme_bw(base_size=10)+
  scale_fill_manual(values=color_compartment) + theme_new +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_core


ggpubr::ggboxplot(data = above.core.meta, x="Location", y="Core_abund", fill = "Comp_age") +
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = Comparison), label = "p.signif", method = "kruskal.test")

###Fungi
core.fun.seed.17 #56
core.fun.leaf.17 #65
core.fun.stem.17 #146

above.core.17<-Reduce(intersect, list(core.fun.seed.17, core.fun.leaf.17, core.fun.stem.17)) #23

fun.clean.rel.17 <- subset_samples(fun.clean.ss.rel,Year == "year_2017")
fun.clean.rel.17<- phyloseq::filter_taxa(fun.clean.rel.17, function(x) sum(x) != 0, TRUE)


fun.clean.rel.17.core <- prune_taxa(taxa_names(fun.clean.rel.17)%in% above.core.17,fun.clean.rel.17)
above.core.sum.17<-colSums(otu_table(fun.clean.rel.17.core))
above.core.sum.17 <- data.frame(above.core.sum.17)
names(above.core.sum.17)[1] <- "Core_abund"
above.core.sum.17$SampleID <- rownames(above.core.sum.17)


f.meta.17.2 <- f.meta.17
f.meta.17$Comp_age <- paste0(f.meta.17$Compartment,"_",f.meta.17$Days)
f.meta.17 <- data.frame(f.meta.17)
above.core.meta<-merge(f.meta.17,above.core.sum.17, by = "SampleID")

above.core.meta$Comp_age <- factor(above.core.meta$Comp_age, levels = c("Leaf_50","Leaf_80","Leaf_120", "Stem_50","Stem_80","Stem_120", 
                                                                        "Root_50","Root_80","Root_120", "Seed_80","Seed_120","Seed_140",
                                                                        "Soil_50","Soil_80","Soil_120", "Soil_140"))


pp.core <- ggplot(data=above.core.meta, aes(x=Comp_age, y=Core_abund))
panel_core=pp.core+geom_boxplot(aes(fill=Compartment)) + geom_point(size=2, position = "jitter",alpha = 0.3, color = "black") +
  scale_y_continuous(name="RA of core OTUs")+theme_bw(base_size=10)+
  scale_fill_manual(values=color_compartment) + theme_new +theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))
panel_core


ggpubr::ggboxplot(data = above.core.meta, x="Location", y="Core_abund", fill = "Comp_age") +
  theme_bw() + theme(aspect.ratio=0.5)+ #scale_fill_manual(values=c("Leaf" = "#336633", "Stem"="#99CC33"))+
  #geom_signif(aes(group = Genus2), step_increase = 0.1, map_signif_level=T)+
  ylab("Jaccard dissimilarity\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")+
  stat_compare_means(aes(group = Comparison), label = "p.signif", method = "kruskal.test")





### Distribution heat map
sample_data(bac.clean.ss.17) <- sample_data(b.meta.17)
bac.clean.ss.17.SW <- subset_samples(bac.clean.ss.17, Location == "Suwon")
bac.clean.ss.17.SW<- phyloseq::filter_taxa(bac.clean.ss.17.SW, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.CC1 <- subset_samples(bac.clean.ss.17, Location == "Chuncheon1")
bac.clean.ss.17.CC1<- phyloseq::filter_taxa(bac.clean.ss.17.CC1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.CC2 <- subset_samples(bac.clean.ss.17, Location == "Chuncheon2")
bac.clean.ss.17.CC2<- phyloseq::filter_taxa(bac.clean.ss.17.CC2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.comp.age <- merge_samples(bac.clean.ss.17.CC2, "Comp_age")
bac.clean.ss.17.comp.age.rel <- microbiome::transform(bac.clean.ss.17.comp.age, "compositional")
otu.bac.clean.ss.17.comp.age.rel<-t(otu_table(bac.clean.ss.17.comp.age.rel))
otu.above.core<-otu.bac.clean.ss.17.comp.age.rel[above.core.17]
otu.above.core <-data.frame(otu.above.core)
otu.above.core$OTU <- rownames(otu.above.core)
bac.list$ordering <- rownames(bac.list)
bac.list_id <- bac.list[c('OTU','OTU_id','ordering')]

above.core.bac<-merge(otu.above.core, bac.list_id, by = "OTU")
rownames(above.core.bac) <- above.core.bac$OTU_id
ncol(above.core.bac)
col.order <-c("Leaf_50","Leaf_80","Leaf_120", "Stem_50","Stem_80","Stem_120", 
              "Root_50","Root_80","Root_120", "Seed_120","Seed_140",
              "Soil_50","Soil_80","Soil_120", "Soil_140","OTU","OTU_id","ordering")
above.core.bac<- above.core.bac[,col.order]

write.xlsx(above.core.bac, "above core bac_SW.xlsx")
write.xlsx(above.core.bac, "above core bac_CC1.xlsx")
write.xlsx(above.core.bac, "above core bac_CC2.xlsx")


###Fungi
sample_data(fun.clean.ss.17) <- sample_data(f.meta.17)
fun.clean.ss.17.SW <- subset_samples(fun.clean.ss.17, Location == "Suwon")
fun.clean.ss.17.SW<- phyloseq::filter_taxa(fun.clean.ss.17.SW, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.CC1 <- subset_samples(fun.clean.ss.17, Location == "Chuncheon1")
fun.clean.ss.17.CC1<- phyloseq::filter_taxa(fun.clean.ss.17.CC1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.CC2 <- subset_samples(fun.clean.ss.17, Location == "Chuncheon2")
fun.clean.ss.17.CC2<- phyloseq::filter_taxa(fun.clean.ss.17.CC2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.comp.age <- merge_samples(fun.clean.ss.17.CC1, "Comp_age")
fun.clean.ss.17.comp.age.rel <- microbiome::transform(fun.clean.ss.17.comp.age, "compositional")
otu.fun.clean.ss.17.comp.age.rel<-t(otu_table(fun.clean.ss.17.comp.age.rel))
otu.above.core<-otu.fun.clean.ss.17.comp.age.rel[above.core.17]
otu.above.core <-data.frame(otu.above.core)
otu.above.core$OTU <- rownames(otu.above.core)
fun.list$ordering <- rownames(fun.list)
fun.list_id <- fun.list[c('OTU','OTU_id','ordering')]

above.core.fun<-merge(otu.above.core, fun.list_id, by = "OTU")
rownames(above.core.fun) <- above.core.fun$OTU_id
ncol(above.core.fun)
col.order <-c("Leaf_50","Leaf_80","Leaf_120", "Stem_50","Stem_80","Stem_120", 
              "Root_50","Root_80","Root_120", "Seed_80","Seed_120","Seed_140",
              "Soil_50","Soil_80","Soil_120", "Soil_140","OTU","OTU_id","ordering")
above.core.fun<- above.core.fun[,col.order]

write.xlsx(above.core.fun, "above core fun_SW.xlsx")
write.xlsx(above.core.fun, "above core fun_CC1.xlsx")
write.xlsx(above.core.fun, "above core fun_CC2.xlsx")
