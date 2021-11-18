#### Vertical transmission
library(tidyverse)
##persistent OTUs
map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)

sample_data(bac.clean.ss) <- sample_data(map) 

(filt.sample <- sample_sums(bac.clean.ss) > 0)
sum(sample_sums(bac.clean.ss) <= 0)  ## 1 sample discarded
bac.clean.ss.f <- prune_samples(filt.sample, bac.clean.ss)
bac.clean.ss.f 


bac.clean.ss.17.G <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
bac.clean.ss.17.G <- phyloseq::filter_taxa(bac.clean.ss.17.G, function(x) sum(x) != 0, TRUE)
bac.clean.ss.18.G <- subset_samples(bac.clean.ss.f, Year == "year_2018" & Compartment == "Seed")
bac.clean.ss.18.G <- phyloseq::filter_taxa(bac.clean.ss.18.G, function(x) sum(x) != 0, TRUE)

sample_data(bac.clean.ss.17.G)
bac.clean.ss.17.G.140 <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment == "Seed" & Location == "Suwon" & Days == 140)
bac.clean.ss.17.G.140 <- phyloseq::filter_taxa(bac.clean.ss.17.G.140, function(x) sum(x) != 0, TRUE)


bac.clean.ss.18.G.0 <-  subset_samples(bac.clean.ss.18.G, Replication == "G_0")
bac.clean.ss.18.G.0 <- phyloseq::filter_taxa(bac.clean.ss.18.G.0, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.G.141 <-  subset_samples(bac.clean.ss.18.G, Replication == "G_141")
bac.clean.ss.18.G.141 <- phyloseq::filter_taxa(bac.clean.ss.18.G.141, function(x) sum(x) != 0, TRUE)


taxa_names(bac.clean.ss.17.G.140)
taxa_names(bac.clean.ss.18.G.0)
taxa_names(bac.clean.ss.18.G.141)

##Only 2018
persistentOTU.bac <- intersect(taxa_names(bac.clean.ss.18.G.0), taxa_names(bac.clean.ss.18.G.141))
tax.persistentOTU.bac <- subset(bac.list, OTU %in% persistentOTU.bac)
#36 OTUs

### persistent OTUs 2017 and 2018
persistent.bac<-Reduce(intersect, list(taxa_names(bac.clean.ss.17.G.140),taxa_names(bac.clean.ss.18.G.0),taxa_names(bac.clean.ss.18.G.141)))
#29 OTUs
tax.persistent.bac <- subset(bac.list, OTU %in% persistent.bac)
#write.csv(tax.persistent.bac,"Bacterial persistent OTUs.csv")

## fungi
(filt.sample <- sample_sums(fun.clean.ss) > 0)
sum(sample_sums(fun.clean.ss) <= 0)  ## 1 sample discarded
fun.clean.ss.f <- prune_samples(filt.sample, fun.clean.ss)
fun.clean.ss.f 

fun.clean.ss.17.G <- subset_samples(fun.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
fun.clean.ss.18.G <- subset_samples(fun.clean.ss.f, Year == "year_2018" & Compartment == "Grain")
fun.clean.ss.18.G <- phyloseq::filter_taxa(fun.clean.ss.18.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.G.140 <- subset_samples(fun.clean.ss.f, Year == "year_2017" & Compartment == "Seed" & Location == "Suwon" & Days == 140)
fun.clean.ss.17.G.140 <- phyloseq::filter_taxa(fun.clean.ss.17.G.140, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.G.0 <-  subset_samples(fun.clean.ss.18.G, Replication == "G_0")
fun.clean.ss.18.G.0 <- phyloseq::filter_taxa(fun.clean.ss.18.G.0, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.G.141 <-  subset_samples(fun.clean.ss.18.G, Replication == "G_141")
fun.clean.ss.18.G.141 <- phyloseq::filter_taxa(fun.clean.ss.18.G.141, function(x) sum(x) != 0, TRUE)


taxa_names(fun.clean.ss.17.G.140)
taxa_names(fun.clean.ss.18.G.0)
taxa_names(fun.clean.ss.18.G.141)


### persistent OTUs 2017 and 2018
persistent.fun<-Reduce(intersect, list(taxa_names(fun.clean.ss.17.G.140),taxa_names(fun.clean.ss.18.G.0),taxa_names(fun.clean.ss.18.G.141)))
#34 OTUs
tax.persistent.fun <- subset(fun.list, OTU %in% persistent.fun)
#write.csv(tax.persistent.fun,"Fungal persistent OTUs.csv")

### Transient OTU
## Bacteria
#2017+2018
transient.bac.1 <- taxa_names(bac.clean.ss.17.G.140)[-which(taxa_names(bac.clean.ss.17.G.140)%in% persistent.bac)]
transient.bac.2 <- taxa_names(bac.clean.ss.18.G.0)[-which(taxa_names(bac.clean.ss.18.G.0)%in% persistent.bac)]
transient.bac.3 <- taxa_names(bac.clean.ss.18.G.141)[-which(taxa_names(bac.clean.ss.18.G.141)%in% persistent.bac)]

transient.bac<-unique(c(transient.bac.1,transient.bac.2,transient.bac.3))
#113 OTUs

## Fungi

#2017+2018
transient.fun.1 <- taxa_names(fun.clean.ss.17.G.140)[-which(taxa_names(fun.clean.ss.17.G.140)%in% persistent.fun)]
transient.fun.2 <- taxa_names(fun.clean.ss.18.G.0)[-which(taxa_names(fun.clean.ss.18.G.0)%in% persistent.fun)]
transient.fun.3 <- taxa_names(fun.clean.ss.18.G.141)[-which(taxa_names(fun.clean.ss.18.G.141)%in% persistent.fun)]

transient.fun<-unique(c(transient.fun.1,transient.fun.2,transient.fun.3))
#189 OTUs


## Proportion of persistent OTUs
length(persistent.bac)/(length(persistent.bac)+length(transient.bac)) #0.2042254
length(persistent.fun)/(length(persistent.fun)+length(transient.fun)) #0.1524664

### Relative abundance of persistent, lost, and gained OTUs in developing seeds
##Bacteria

##2017 + 2018
bac.clean.ss.G <- subset_samples(bac.clean.ss.f, Compartment == "Seed" & Location == "Suwon")
bac.clean.ss.G <- phyloseq::filter_taxa(bac.clean.ss.G, function(x) sum(x) != 0, TRUE)


bac.G.comp.rel <- transform(bac.clean.ss.G, 'compositional')

otu.bac.G.comp.rel <- otu_table(bac.G.comp.rel)
otu.bac.G.comp.rel <- data.frame(otu.bac.G.comp.rel)
otu.bac.G.comp.rel$OTU <- rownames(otu.bac.G.comp.rel)

melt.otu.bac.G.comp.rel <- melt(otu.bac.G.comp.rel)
head(melt.otu.bac.G.comp.rel)
melt.otu.bac.G.comp.rel$Days <- ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "106"), "106", 
                                       ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "120"),"120",
                                              ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "90"), "90",
                                                     ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "76"),"76",
                                                            ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "141"),"141", 
                                                                   ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "80"),"80",
                                                                          ifelse(str_detect(melt.otu.bac.G.comp.rel$variable, "140"),"140","0")))))))

meta.bac.G<-sample_data(bac.clean.ss.G)
meta.bac.G<-data.frame(meta.bac.G)


for (i in as.character(meta.bac.G$SampleID)) {
  if (i %in% as.character(melt.otu.bac.G.comp.rel$variable) == T){
    melt.otu.bac.G.comp.rel$Year[which(melt.otu.bac.G.comp.rel$variable == i)] <- as.character(meta.bac.G$Year[which(meta.bac.G$SampleID == i)])
  }
}

head(melt.otu.bac.G.comp.rel)

otu.bac.G.comp.rel.sum.tab <- melt.otu.bac.G.comp.rel %>% group_by(OTU, Days, Year) %>% summarise(Mean = mean(value))

### Assign origin
persistent.bac
transient.bac

bac.clean.ss.f.soil<- subset_samples(bac.clean.ss.f, Compartment %in% c("Soil","Bulk_soil","Rhizosphere")& Location =="Suwon")
bac.clean.ss.f.soil<- phyloseq::filter_taxa(bac.clean.ss.f.soil, function(x) sum(x) != 0, TRUE)

soil.bac<-taxa_names(bac.clean.ss.f.soil)

transient.bac.soil <- intersect(transient.bac, soil.bac)
transient.bac.nonsoil <- transient.bac[which(!(transient.bac %in% transient.bac.soil))]

otu.bac.G.comp.rel.sum.tab$Class <- ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% persistent.bac, "Persistent",
                                           ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% transient.bac.soil, "Transient_soil",
                                                  ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% transient.bac.nonsoil,"Transient_nonsoil", 
                                                         ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% soil.bac,"Soil","Unknown"))))

otu.bac.G.comp.rel.sum.tab.2 <- otu.bac.G.comp.rel.sum.tab %>% group_by(Days, Class, Year) %>% summarise(Total = sum(Mean))

otu.bac.G.comp.rel.sum.tab.2$Days <- factor(otu.bac.G.comp.rel.sum.tab.2$Days, levels = c("0","76","80","90","106","120","140","141"))
otu.bac.G.comp.rel.sum.tab.2$Class <- factor(otu.bac.G.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))

### plotting
p <- ggplot(otu.bac.G.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Year, scales = 'free')+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p


## Fungi
##2017 + 2018
fun.clean.ss.G <- subset_samples(fun.clean.ss.f, Compartment %in% c("Seed","Grain")  & Location == "Suwon")
fun.clean.ss.G <- phyloseq::filter_taxa(fun.clean.ss.G, function(x) sum(x) != 0, TRUE)


fun.G.comp.rel <- transform(fun.clean.ss.G, 'compositional')

otu.fun.G.comp.rel <- otu_table(fun.G.comp.rel)
otu.fun.G.comp.rel <- data.frame(otu.fun.G.comp.rel)
otu.fun.G.comp.rel$OTU <- rownames(otu.fun.G.comp.rel)

melt.otu.fun.G.comp.rel <- melt(otu.fun.G.comp.rel)
head(melt.otu.fun.G.comp.rel)
melt.otu.fun.G.comp.rel$Days <- ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "106"), "106", 
                                       ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "120"),"120",
                                              ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "90"), "90",
                                                     ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "76"),"76",
                                                            ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "141"),"141", 
                                                                   ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "80"),"80",
                                                                          ifelse(str_detect(melt.otu.fun.G.comp.rel$variable, "140"),"140","0")))))))

meta.fun.G<-sample_data(fun.clean.ss.G)
meta.fun.G<-data.frame(meta.fun.G)


for (i in as.character(meta.fun.G$SampleID)) {
  if (i %in% as.character(melt.otu.fun.G.comp.rel$variable) == T){
    melt.otu.fun.G.comp.rel$Year[which(melt.otu.fun.G.comp.rel$variable == i)] <- as.character(meta.fun.G$Year[which(meta.fun.G$SampleID == i)])
  }
}

head(melt.otu.fun.G.comp.rel)

otu.fun.G.comp.rel.sum.tab <- melt.otu.fun.G.comp.rel %>% group_by(OTU, Days, Year) %>% summarise(Mean = mean(value))

### Assign origin
persistent.fun
transient.fun

fun.clean.ss.f.soil<- subset_samples(fun.clean.ss.f, Compartment %in% c("Soil","Bulk_soil","Rhizosphere")& Location =="Suwon")
fun.clean.ss.f.soil<- phyloseq::filter_taxa(fun.clean.ss.f.soil, function(x) sum(x) != 0, TRUE)

soil.fun<-taxa_names(fun.clean.ss.f.soil)

transient.fun.soil <- intersect(transient.fun, soil.fun)
transient.fun.nonsoil <- transient.fun[which(!(transient.fun %in% transient.fun.soil))]

otu.fun.G.comp.rel.sum.tab$Class <- ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% persistent.fun, "Persistent",
                                           ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% transient.fun.soil, "Transient_soil",
                                                  ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% transient.fun.nonsoil,"Transient_nonsoil", 
                                                         ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% soil.fun,"Soil","Unknown"))))

otu.fun.G.comp.rel.sum.tab.2 <- otu.fun.G.comp.rel.sum.tab %>% group_by(Days, Class, Year) %>% summarise(Total = sum(Mean))

otu.fun.G.comp.rel.sum.tab.2$Days <- factor(otu.fun.G.comp.rel.sum.tab.2$Days, levels = c("0","76","80","90","106","120","140","141"))
otu.fun.G.comp.rel.sum.tab.2$Class <- factor(otu.fun.G.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))

### plotting
p <- ggplot(otu.fun.G.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Year, scales = 'free')+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.fun.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_fungi.csv")
write.csv(otu.bac.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_bacteria.csv")

### Proportion of persistent and transient OTUs

###Distribution of persistent, transient OTUs in other compartment
##Bacteria
bac.clean.ss.wo.G <- subset_samples(bac.clean.ss.f, Year == "year_2018" & Microhabitat != "G")
bac.clean.ss.wo.G <- phyloseq::filter_taxa(bac.clean.ss.wo.G, function(x) sum(x) != 0, TRUE)


bac.wo.G.comp.rel <- transform(bac.clean.ss.wo.G, 'compositional')

otu.bac.wo.G.comp.rel <- otu_table(bac.wo.G.comp.rel)
otu.bac.wo.G.comp.rel <- data.frame(otu.bac.wo.G.comp.rel)
otu.bac.wo.G.comp.rel$OTU <- rownames(otu.bac.wo.G.comp.rel)

melt.otu.bac.wo.G.comp.rel <- melt(otu.bac.wo.G.comp.rel)
head(melt.otu.bac.wo.G.comp.rel)
melt.otu.bac.wo.G.comp.rel$Days <- ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "106"), "106", 
                                       ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "120"),"120",
                                              ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "62"), "62",
                                                     ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "76"),"76",
                                                            ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "141"),"141", 
                                                                   ifelse(str_detect(melt.otu.bac.wo.G.comp.rel$variable, "48"),"48", "90"))))))
                                                                          

meta.bac.wo.G<-sample_data(bac.clean.ss.wo.G)
meta.bac.wo.G<-data.frame(meta.bac.wo.G)


for (i in as.character(meta.bac.wo.G$SampleID)) {
  if (i %in% as.character(melt.otu.bac.wo.G.comp.rel$variable) == T){
    melt.otu.bac.wo.G.comp.rel$Section[which(melt.otu.bac.wo.G.comp.rel$variable == i)] <- as.character(meta.bac.wo.G$Microhabitat[which(meta.bac.wo.G$SampleID == i)])
  }
}

head(melt.otu.bac.wo.G.comp.rel)

otu.bac.wo.G.comp.rel.sum.tab <- melt.otu.bac.wo.G.comp.rel %>% group_by(OTU, Days, Section) %>% summarise(Mean = mean(value))

### Assign origin
otu.bac.wo.G.comp.rel.sum.tab$Class <- ifelse(otu.bac.wo.G.comp.rel.sum.tab$OTU %in% persistent.bac, "Persistent",
                                           ifelse(otu.bac.wo.G.comp.rel.sum.tab$OTU %in% transient.bac.soil, "Transient_soil",
                                                  ifelse(otu.bac.wo.G.comp.rel.sum.tab$OTU %in% transient.bac.nonsoil,"Transient_nonsoil", 
                                                         ifelse(otu.bac.wo.G.comp.rel.sum.tab$OTU %in% soil.bac,"Soil","Unknown"))))

otu.bac.wo.G.comp.rel.sum.tab.2 <- otu.bac.wo.G.comp.rel.sum.tab %>% group_by(Days, Class, Section) %>% summarise(Total = sum(Mean))

otu.bac.wo.G.comp.rel.sum.tab.2$Days <- factor(otu.bac.wo.G.comp.rel.sum.tab.2$Days, levels = c("48","62","76","90","106","120","141"))
otu.bac.wo.G.comp.rel.sum.tab.2$Class <- factor(otu.bac.wo.G.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))
otu.bac.wo.G.comp.rel.sum.tab.2$Section <- factor(otu.bac.wo.G.comp.rel.sum.tab.2$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

### plotting
p <- ggplot(otu.bac.wo.G.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Section, nrow = 2)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p


##Fungi
fun.clean.ss.wo.G <- subset_samples(fun.clean.ss.f, Year == "year_2018" & Microhabitat != "G")
fun.clean.ss.wo.G <- phyloseq::filter_taxa(fun.clean.ss.wo.G, function(x) sum(x) != 0, TRUE)


fun.wo.G.comp.rel <- transform(fun.clean.ss.wo.G, 'compositional')

otu.fun.wo.G.comp.rel <- otu_table(fun.wo.G.comp.rel)
otu.fun.wo.G.comp.rel <- data.frame(otu.fun.wo.G.comp.rel)
otu.fun.wo.G.comp.rel$OTU <- rownames(otu.fun.wo.G.comp.rel)

melt.otu.fun.wo.G.comp.rel <- melt(otu.fun.wo.G.comp.rel)
head(melt.otu.fun.wo.G.comp.rel)
melt.otu.fun.wo.G.comp.rel$Days <- ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "106"), "106", 
                                          ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "120"),"120",
                                                 ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "62"), "62",
                                                        ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "76"),"76",
                                                               ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "141"),"141", 
                                                                      ifelse(str_detect(melt.otu.fun.wo.G.comp.rel$variable, "48"),"48", "90"))))))


meta.fun.wo.G<-sample_data(fun.clean.ss.wo.G)
meta.fun.wo.G<-data.frame(meta.fun.wo.G)


for (i in as.character(meta.fun.wo.G$SampleID)) {
  if (i %in% as.character(melt.otu.fun.wo.G.comp.rel$variable) == T){
    melt.otu.fun.wo.G.comp.rel$Section[which(melt.otu.fun.wo.G.comp.rel$variable == i)] <- as.character(meta.fun.wo.G$Microhabitat[which(meta.fun.wo.G$SampleID == i)])
  }
}

head(melt.otu.fun.wo.G.comp.rel)

otu.fun.wo.G.comp.rel.sum.tab <- melt.otu.fun.wo.G.comp.rel %>% group_by(OTU, Days, Section) %>% summarise(Mean = mean(value))

### Assign origin
otu.fun.wo.G.comp.rel.sum.tab$Class <- ifelse(otu.fun.wo.G.comp.rel.sum.tab$OTU %in% persistent.fun, "Persistent",
                                              ifelse(otu.fun.wo.G.comp.rel.sum.tab$OTU %in% transient.fun.soil, "Transient_soil",
                                                     ifelse(otu.fun.wo.G.comp.rel.sum.tab$OTU %in% transient.fun.nonsoil,"Transient_nonsoil", 
                                                            ifelse(otu.fun.wo.G.comp.rel.sum.tab$OTU %in% soil.fun,"Soil","Unknown"))))

otu.fun.wo.G.comp.rel.sum.tab.2 <- otu.fun.wo.G.comp.rel.sum.tab %>% group_by(Days, Class, Section) %>% summarise(Total = sum(Mean))

otu.fun.wo.G.comp.rel.sum.tab.2$Days <- factor(otu.fun.wo.G.comp.rel.sum.tab.2$Days, levels = c("48","62","76","90","106","120","141"))
otu.fun.wo.G.comp.rel.sum.tab.2$Class <- factor(otu.fun.wo.G.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))
otu.fun.wo.G.comp.rel.sum.tab.2$Section <- factor(otu.fun.wo.G.comp.rel.sum.tab.2$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

### plotting
p <- ggplot(otu.fun.wo.G.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Section, nrow = 2)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.fun.wo.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_fungi_other_sections.csv")
write.csv(otu.bac.wo.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_bacteria_other_sections.csv")


### Proportion of persistent, transient, and other OTUs
#### Proportion of number of OTUs (not relative abundance)
### In seeds

bac.clean.ss.G

otu.bac.G <- otu_table(bac.clean.ss.G)
otu.bac.G <- data.frame(otu.bac.G)
otu.bac.G$OTU <- rownames(otu.bac.G)

melt.otu.bac.G <- melt(otu.bac.G)
head(melt.otu.bac.G)
melt.otu.bac.G$Days <- ifelse(str_detect(melt.otu.bac.G$variable, "106"), "106", 
                                       ifelse(str_detect(melt.otu.bac.G$variable, "120"),"120",
                                              ifelse(str_detect(melt.otu.bac.G$variable, "90"), "90",
                                                     ifelse(str_detect(melt.otu.bac.G$variable, "76"),"76",
                                                            ifelse(str_detect(melt.otu.bac.G$variable, "141"),"141", 
                                                                   ifelse(str_detect(melt.otu.bac.G$variable, "80"),"80",
                                                                          ifelse(str_detect(melt.otu.bac.G$variable, "140"),"140","0")))))))

meta.bac.G<-sample_data(bac.clean.ss.G)
meta.bac.G<-data.frame(meta.bac.G)


for (i in as.character(meta.bac.G$SampleID)) {
  if (i %in% as.character(melt.otu.bac.G$variable) == T){
    melt.otu.bac.G$Year[which(melt.otu.bac.G$variable == i)] <- as.character(meta.bac.G$Year[which(meta.bac.G$SampleID == i)])
  }
}

head(melt.otu.bac.G)


otu.bac.G.comp.rel.sum.tab <- melt.otu.bac.G %>% group_by(OTU, Days, Year) %>% summarise(Mean = mean(value))

otu.bac.G.comp.rel.sum.tab$Class <- ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% persistent.bac, "Persistent",
                                              ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% transient.bac.soil, "Transient_soil",
                                                     ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% transient.bac.nonsoil,"Transient_nonsoil", 
                                                            ifelse(otu.bac.G.comp.rel.sum.tab$OTU %in% soil.bac,"Soil","Unknown"))))



##count
otu.bac.G.comp.rel.sum.tab$Count <- 0
otu.bac.G.comp.rel.sum.tab$Count[which(otu.bac.G.comp.rel.sum.tab$Mean > 0)] <- 1


otu.bac.G.comp.rel.sum.tab.3 <- otu.bac.G.comp.rel.sum.tab %>% group_by(Days, Class, Year) %>% summarise(Total = sum(Count))

otu.bac.G.comp.rel.sum.tab.3.total <- otu.bac.G.comp.rel.sum.tab %>% group_by(Days, Year) %>% summarise(Total_count = sum(Count))


otu.bac.G.comp.rel.sum.tab.3$Total_count <- 0

for (i in as.character(unique(otu.bac.G.comp.rel.sum.tab.3.total$Year))) {
  for (j in as.character(unique(otu.bac.G.comp.rel.sum.tab.3.total$Days))){
    if (i %in% as.character(otu.bac.G.comp.rel.sum.tab.3$Year) == T & j %in% as.character(otu.bac.G.comp.rel.sum.tab.3$Days) == T){
      otu.bac.G.comp.rel.sum.tab.3$Total_count[which(otu.bac.G.comp.rel.sum.tab.3$Year == i & otu.bac.G.comp.rel.sum.tab.3$Days == j)] <- as.character(otu.bac.G.comp.rel.sum.tab.3.total$Total_count[which(otu.bac.G.comp.rel.sum.tab.3.total$Year == i & otu.bac.G.comp.rel.sum.tab.3.total$Days == j)])
    }
  }
}

otu.bac.G.comp.rel.sum.tab.3$Total <- as.numeric(as.character(otu.bac.G.comp.rel.sum.tab.3$Total))
otu.bac.G.comp.rel.sum.tab.3$Total_count <- as.numeric(as.character(otu.bac.G.comp.rel.sum.tab.3$Total_count))

otu.bac.G.comp.rel.sum.tab.3$Proportion <- otu.bac.G.comp.rel.sum.tab.3$Total/otu.bac.G.comp.rel.sum.tab.3$Total_count

otu.bac.G.comp.rel.sum.tab.3$Days <- factor(otu.bac.G.comp.rel.sum.tab.3$Days, levels = c("0","76","80","90","106","120","140","141"))
otu.bac.G.comp.rel.sum.tab.3$Class <- factor(otu.bac.G.comp.rel.sum.tab.3$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))


### plotting
p <- ggplot(otu.bac.G.comp.rel.sum.tab.3, aes(x=Days, y = Proportion, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  facet_wrap(~Year, scales = 'free')+
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

p <- ggplot(otu.bac.G.comp.rel.sum.tab.3, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  facet_wrap(~Year, scales = 'free')+
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Count\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.bac.G.comp.rel.sum.tab.3, "OTU counts of persistent and transient OTUs in seeds_bacteria.csv")


## Other sections
otu.bac.18 <- otu_table(bac.clean.ss.wo.G)
otu.bac.18 <- data.frame(otu.bac.18)
otu.bac.18$OTU <- rownames(otu.bac.18)

melt.otu.bac.18 <- melt(otu.bac.18)
head(melt.otu.bac.18)
melt.otu.bac.18$Days <- ifelse(str_detect(melt.otu.bac.18$variable, "106"), "106", 
                               ifelse(str_detect(melt.otu.bac.18$variable, "120"),"120",
                                      ifelse(str_detect(melt.otu.bac.18$variable, "90"), "90",
                                             ifelse(str_detect(melt.otu.bac.18$variable, "76"),"76",
                                                    ifelse(str_detect(melt.otu.bac.18$variable, "62"),"62",
                                                           ifelse(str_detect(melt.otu.bac.18$variable, "141"),"141",
                                                                  ifelse(str_detect(melt.otu.bac.18$variable, "48"),"48","0")))))))
meta.bac.18<-sample_data(bac.clean.ss.wo.G)
meta.bac.18<-data.frame(meta.bac.18)


for (i in as.character(meta.bac.18$SampleID)) {
  if (i %in% as.character(melt.otu.bac.18$variable) == T){
    melt.otu.bac.18$Section[which(melt.otu.bac.18$variable == i)] <- as.character(meta.bac.18$Microhabitat[which(meta.bac.18$SampleID == i)])
  }
}

head(melt.otu.bac.18)


otu.bac.other.comp.rel.sum.tab <- melt.otu.bac.18 %>% group_by(OTU, Days, Section) %>% summarise(Mean = mean(value))

otu.bac.other.comp.rel.sum.tab$Class <- ifelse(otu.bac.other.comp.rel.sum.tab$OTU %in% persistent.bac, "Persistent",
                                              ifelse(otu.bac.other.comp.rel.sum.tab$OTU %in% transient.bac.soil, "Transient_soil",
                                                     ifelse(otu.bac.other.comp.rel.sum.tab$OTU %in% transient.bac.nonsoil,"Transient_nonsoil", 
                                                            ifelse(otu.bac.other.comp.rel.sum.tab$OTU %in% soil.bac,"Soil","Unknown"))))

##count
otu.bac.other.comp.rel.sum.tab$Count <- 0
otu.bac.other.comp.rel.sum.tab$Count[which(otu.bac.other.comp.rel.sum.tab$Mean > 0)] <- 1


otu.bac.other.comp.rel.sum.tab.3 <- otu.bac.other.comp.rel.sum.tab %>% group_by(Days, Class, Section) %>% summarise(Total = sum(Count))

otu.bac.other.comp.rel.sum.tab.3.total <- otu.bac.other.comp.rel.sum.tab %>% group_by(Days, Section) %>% summarise(Total_count = sum(Count))


otu.bac.other.comp.rel.sum.tab.3$Total_count <- 0

for (i in as.character(unique(otu.bac.other.comp.rel.sum.tab.3.total$Section))) {
  for (j in as.character(unique(otu.bac.other.comp.rel.sum.tab.3.total$Days))){
    if (i %in% as.character(otu.bac.other.comp.rel.sum.tab.3$Section) == T & j %in% as.character(otu.bac.other.comp.rel.sum.tab.3$Days) == T){
      otu.bac.other.comp.rel.sum.tab.3$Total_count[which(otu.bac.other.comp.rel.sum.tab.3$Section == i & otu.bac.other.comp.rel.sum.tab.3$Days == j)] <- as.character(otu.bac.other.comp.rel.sum.tab.3.total$Total_count[which(otu.bac.other.comp.rel.sum.tab.3.total$Section == i & otu.bac.other.comp.rel.sum.tab.3.total$Days == j)])
    }
  }
}

otu.bac.other.comp.rel.sum.tab.3$Total <- as.numeric(as.character(otu.bac.other.comp.rel.sum.tab.3$Total))
otu.bac.other.comp.rel.sum.tab.3$Total_count <- as.numeric(as.character(otu.bac.other.comp.rel.sum.tab.3$Total_count))

otu.bac.other.comp.rel.sum.tab.3$Proportion <- otu.bac.other.comp.rel.sum.tab.3$Total/otu.bac.other.comp.rel.sum.tab.3$Total_count

otu.bac.other.comp.rel.sum.tab.3$Days <- factor(otu.bac.other.comp.rel.sum.tab.3$Days, levels = c("48","62","76","90","106","120","141"))
otu.bac.other.comp.rel.sum.tab.3$Section <- factor(otu.bac.other.comp.rel.sum.tab.3$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))
otu.bac.other.comp.rel.sum.tab.3$Class <- factor(otu.bac.other.comp.rel.sum.tab.3$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))

### plotting
p <- ggplot(otu.bac.other.comp.rel.sum.tab.3, aes(x=Days, y = Proportion, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  facet_wrap(~ Section, nrow =2) +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.bac.other.comp.rel.sum.tab.3, "OTU counts of persistent and transient OTUs in other sections_bacteria.csv")

## Fungi
fun.clean.ss.G

otu.fun.G <- otu_table(fun.clean.ss.G)
otu.fun.G <- data.frame(otu.fun.G)
otu.fun.G$OTU <- rownames(otu.fun.G)

melt.otu.fun.G <- melt(otu.fun.G)
head(melt.otu.fun.G)
melt.otu.fun.G$Days <- ifelse(str_detect(melt.otu.fun.G$variable, "106"), "106", 
                              ifelse(str_detect(melt.otu.fun.G$variable, "120"),"120",
                                     ifelse(str_detect(melt.otu.fun.G$variable, "90"), "90",
                                            ifelse(str_detect(melt.otu.fun.G$variable, "76"),"76",
                                                   ifelse(str_detect(melt.otu.fun.G$variable, "141"),"141", 
                                                          ifelse(str_detect(melt.otu.fun.G$variable, "80"),"80",
                                                                 ifelse(str_detect(melt.otu.fun.G$variable, "140"),"140","0")))))))

meta.fun.G<-sample_data(fun.clean.ss.G)
meta.fun.G<-data.frame(meta.fun.G)


for (i in as.character(meta.fun.G$SampleID)) {
  if (i %in% as.character(melt.otu.fun.G$variable) == T){
    melt.otu.fun.G$Year[which(melt.otu.fun.G$variable == i)] <- as.character(meta.fun.G$Year[which(meta.fun.G$SampleID == i)])
  }
}

head(melt.otu.fun.G)


otu.fun.G.comp.rel.sum.tab <- melt.otu.fun.G %>% group_by(OTU, Days, Year) %>% summarise(Mean = mean(value))

otu.fun.G.comp.rel.sum.tab$Class <- ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% persistent.fun, "Persistent",
                                           ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% transient.fun.soil, "Transient_soil",
                                                  ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% transient.fun.nonsoil,"Transient_nonsoil", 
                                                         ifelse(otu.fun.G.comp.rel.sum.tab$OTU %in% soil.fun,"Soil","Unknown"))))



##count
otu.fun.G.comp.rel.sum.tab$Count <- 0
otu.fun.G.comp.rel.sum.tab$Count[which(otu.fun.G.comp.rel.sum.tab$Mean > 0)] <- 1


otu.fun.G.comp.rel.sum.tab.3 <- otu.fun.G.comp.rel.sum.tab %>% group_by(Days, Class, Year) %>% summarise(Total = sum(Count))

otu.fun.G.comp.rel.sum.tab.3.total <- otu.fun.G.comp.rel.sum.tab %>% group_by(Days, Year) %>% summarise(Total_count = sum(Count))


otu.fun.G.comp.rel.sum.tab.3$Total_count <- 0

for (i in as.character(unique(otu.fun.G.comp.rel.sum.tab.3.total$Year))) {
  for (j in as.character(unique(otu.fun.G.comp.rel.sum.tab.3.total$Days))){
    if (i %in% as.character(otu.fun.G.comp.rel.sum.tab.3$Year) == T & j %in% as.character(otu.fun.G.comp.rel.sum.tab.3$Days) == T){
      otu.fun.G.comp.rel.sum.tab.3$Total_count[which(otu.fun.G.comp.rel.sum.tab.3$Year == i & otu.fun.G.comp.rel.sum.tab.3$Days == j)] <- as.character(otu.fun.G.comp.rel.sum.tab.3.total$Total_count[which(otu.fun.G.comp.rel.sum.tab.3.total$Year == i & otu.fun.G.comp.rel.sum.tab.3.total$Days == j)])
    }
  }
}

otu.fun.G.comp.rel.sum.tab.3$Total <- as.numeric(as.character(otu.fun.G.comp.rel.sum.tab.3$Total))
otu.fun.G.comp.rel.sum.tab.3$Total_count <- as.numeric(as.character(otu.fun.G.comp.rel.sum.tab.3$Total_count))

otu.fun.G.comp.rel.sum.tab.3$Proportion <- otu.fun.G.comp.rel.sum.tab.3$Total/otu.fun.G.comp.rel.sum.tab.3$Total_count

otu.fun.G.comp.rel.sum.tab.3$Days <- factor(otu.fun.G.comp.rel.sum.tab.3$Days, levels = c("0","76","80","90","106","120","140","141"))
otu.fun.G.comp.rel.sum.tab.3$Class <- factor(otu.fun.G.comp.rel.sum.tab.3$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))


### plotting
p <- ggplot(otu.fun.G.comp.rel.sum.tab.3, aes(x=Days, y = Proportion, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  facet_wrap(~Year, scales = 'free')+
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

p <- ggplot(otu.fun.G.comp.rel.sum.tab.3, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  facet_wrap(~Year, scales = 'free')+
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Count\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.fun.G.comp.rel.sum.tab.3, "OTU counts of persistent and transient OTUs in seeds_fungi.csv")


## Other sections
otu.fun.18 <- otu_table(fun.clean.ss.wo.G)
otu.fun.18 <- data.frame(otu.fun.18)
otu.fun.18$OTU <- rownames(otu.fun.18)

melt.otu.fun.18 <- melt(otu.fun.18)
head(melt.otu.fun.18)
melt.otu.fun.18$Days <- ifelse(str_detect(melt.otu.fun.18$variable, "106"), "106", 
                               ifelse(str_detect(melt.otu.fun.18$variable, "120"),"120",
                                      ifelse(str_detect(melt.otu.fun.18$variable, "90"), "90",
                                             ifelse(str_detect(melt.otu.fun.18$variable, "76"),"76",
                                                    ifelse(str_detect(melt.otu.fun.18$variable, "62"),"62",
                                                           ifelse(str_detect(melt.otu.fun.18$variable, "141"),"141",
                                                                  ifelse(str_detect(melt.otu.fun.18$variable, "48"),"48","0")))))))
meta.fun.18<-sample_data(fun.clean.ss.wo.G)
meta.fun.18<-data.frame(meta.fun.18)


for (i in as.character(meta.fun.18$SampleID)) {
  if (i %in% as.character(melt.otu.fun.18$variable) == T){
    melt.otu.fun.18$Section[which(melt.otu.fun.18$variable == i)] <- as.character(meta.fun.18$Microhabitat[which(meta.fun.18$SampleID == i)])
  }
}

head(melt.otu.fun.18)


otu.fun.other.comp.rel.sum.tab <- melt.otu.fun.18 %>% group_by(OTU, Days, Section) %>% summarise(Mean = mean(value))

otu.fun.other.comp.rel.sum.tab$Class <- ifelse(otu.fun.other.comp.rel.sum.tab$OTU %in% persistent.fun, "Persistent",
                                               ifelse(otu.fun.other.comp.rel.sum.tab$OTU %in% transient.fun.soil, "Transient_soil",
                                                      ifelse(otu.fun.other.comp.rel.sum.tab$OTU %in% transient.fun.nonsoil,"Transient_nonsoil", 
                                                             ifelse(otu.fun.other.comp.rel.sum.tab$OTU %in% soil.fun,"Soil","Unknown"))))

##count
otu.fun.other.comp.rel.sum.tab$Count <- 0
otu.fun.other.comp.rel.sum.tab$Count[which(otu.fun.other.comp.rel.sum.tab$Mean > 0)] <- 1


otu.fun.other.comp.rel.sum.tab.3 <- otu.fun.other.comp.rel.sum.tab %>% group_by(Days, Class, Section) %>% summarise(Total = sum(Count))

otu.fun.other.comp.rel.sum.tab.3.total <- otu.fun.other.comp.rel.sum.tab %>% group_by(Days, Section) %>% summarise(Total_count = sum(Count))


otu.fun.other.comp.rel.sum.tab.3$Total_count <- 0

for (i in as.character(unique(otu.fun.other.comp.rel.sum.tab.3.total$Section))) {
  for (j in as.character(unique(otu.fun.other.comp.rel.sum.tab.3.total$Days))){
    if (i %in% as.character(otu.fun.other.comp.rel.sum.tab.3$Section) == T & j %in% as.character(otu.fun.other.comp.rel.sum.tab.3$Days) == T){
      otu.fun.other.comp.rel.sum.tab.3$Total_count[which(otu.fun.other.comp.rel.sum.tab.3$Section == i & otu.fun.other.comp.rel.sum.tab.3$Days == j)] <- as.character(otu.fun.other.comp.rel.sum.tab.3.total$Total_count[which(otu.fun.other.comp.rel.sum.tab.3.total$Section == i & otu.fun.other.comp.rel.sum.tab.3.total$Days == j)])
    }
  }
}

otu.fun.other.comp.rel.sum.tab.3$Total <- as.numeric(as.character(otu.fun.other.comp.rel.sum.tab.3$Total))
otu.fun.other.comp.rel.sum.tab.3$Total_count <- as.numeric(as.character(otu.fun.other.comp.rel.sum.tab.3$Total_count))

otu.fun.other.comp.rel.sum.tab.3$Proportion <- otu.fun.other.comp.rel.sum.tab.3$Total/otu.fun.other.comp.rel.sum.tab.3$Total_count

otu.fun.other.comp.rel.sum.tab.3$Days <- factor(otu.fun.other.comp.rel.sum.tab.3$Days, levels = c("48","62","76","90","106","120","141"))
otu.fun.other.comp.rel.sum.tab.3$Section <- factor(otu.fun.other.comp.rel.sum.tab.3$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))
otu.fun.other.comp.rel.sum.tab.3$Class <- factor(otu.fun.other.comp.rel.sum.tab.3$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))

### plotting
p <- ggplot(otu.fun.other.comp.rel.sum.tab.3, aes(x=Days, y = Proportion, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  facet_wrap(~ Section, nrow =2) +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.fun.other.comp.rel.sum.tab.3, "OTU counts of persistent and transient OTUs in other sections_funteria.csv")


dev.off()

### 2017 samples
###Distribution of persistent, transient OTUs in other compartment
##Bacteria
bac.clean.ss.wo.G.17 <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment != "Seed")
bac.clean.ss.wo.G.17 <- phyloseq::filter_taxa(bac.clean.ss.wo.G.17, function(x) sum(x) != 0, TRUE)


bac.wo.G.17.comp.rel <- transform(bac.clean.ss.wo.G.17, 'compositional')

otu.bac.wo.G.17.comp.rel <- otu_table(bac.wo.G.17.comp.rel)
otu.bac.wo.G.17.comp.rel <- data.frame(otu.bac.wo.G.17.comp.rel)
otu.bac.wo.G.17.comp.rel$OTU <- rownames(otu.bac.wo.G.17.comp.rel)

melt.otu.bac.wo.G.17.comp.rel <- melt(otu.bac.wo.G.17.comp.rel)
head(melt.otu.bac.wo.G.17.comp.rel)
melt.otu.bac.wo.G.17.comp.rel$Days <- ifelse(str_detect(melt.otu.bac.wo.G.17.comp.rel$variable, "50"), "50", 
                                          ifelse(str_detect(melt.otu.bac.wo.G.17.comp.rel$variable, "80"),"80",
                                                 ifelse(str_detect(melt.otu.bac.wo.G.17.comp.rel$variable, "120"), "120","140")))
                                                     

meta.bac.wo.G.17<-sample_data(bac.clean.ss.wo.G.17)
meta.bac.wo.G.17<-data.frame(meta.bac.wo.G.17)


for (i in as.character(meta.bac.wo.G.17$SampleID)) {
  if (i %in% as.character(melt.otu.bac.wo.G.17.comp.rel$variable) == T){
    melt.otu.bac.wo.G.17.comp.rel$Section[which(melt.otu.bac.wo.G.17.comp.rel$variable == i)] <- as.character(meta.bac.wo.G.17$Compartment[which(meta.bac.wo.G.17$SampleID == i)])
    melt.otu.bac.wo.G.17.comp.rel$Location[which(melt.otu.bac.wo.G.17.comp.rel$variable == i)] <- as.character(meta.bac.wo.G.17$Location[which(meta.bac.wo.G.17$SampleID == i)])
  }
}

head(melt.otu.bac.wo.G.17.comp.rel)

otu.bac.wo.G.17.comp.rel.sum.tab <- melt.otu.bac.wo.G.17.comp.rel %>% group_by(OTU, Days, Section, Location) %>% summarise(Mean = mean(value))

### Assign origin
otu.bac.wo.G.17.comp.rel.sum.tab$Class <- ifelse(otu.bac.wo.G.17.comp.rel.sum.tab$OTU %in% persistent.bac, "Persistent",
                                              ifelse(otu.bac.wo.G.17.comp.rel.sum.tab$OTU %in% transient.bac.soil, "Transient_soil",
                                                     ifelse(otu.bac.wo.G.17.comp.rel.sum.tab$OTU %in% transient.bac.nonsoil,"Transient_nonsoil", 
                                                            ifelse(otu.bac.wo.G.17.comp.rel.sum.tab$OTU %in% soil.bac,"Soil","Unknown"))))

otu.bac.wo.G.17.comp.rel.sum.tab.2 <- otu.bac.wo.G.17.comp.rel.sum.tab %>% group_by(Days, Class, Section,Location) %>% summarise(Total = sum(Mean))

otu.bac.wo.G.17.comp.rel.sum.tab.2$Days <- factor(otu.bac.wo.G.17.comp.rel.sum.tab.2$Days, levels = c("50","80","120","140"))
otu.bac.wo.G.17.comp.rel.sum.tab.2$Class <- factor(otu.bac.wo.G.17.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))
otu.bac.wo.G.17.comp.rel.sum.tab.2$Section <- factor(otu.bac.wo.G.17.comp.rel.sum.tab.2$Section, levels = c("Soil","Root","Stem","Leaf"))

### plotting
p <- ggplot(otu.bac.wo.G.17.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Section+Location, nrow = 2)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p


##Fungi
fun.clean.ss.wo.G.17 <- subset_samples(fun.clean.ss.f, Year == "year_2018" & Microhabitat != "wo.G.17")
fun.clean.ss.wo.G.17 <- phyloseq::filter_taxa(fun.clean.ss.wo.G.17, function(x) sum(x) != 0, TRUE)


fun.wo.G.17.comp.rel <- transform(fun.clean.ss.wo.G.17, 'compositional')

otu.fun.wo.G.17.comp.rel <- otu_table(fun.wo.G.17.comp.rel)
otu.fun.wo.G.17.comp.rel <- data.frame(otu.fun.wo.G.17.comp.rel)
otu.fun.wo.G.17.comp.rel$OTU <- rownames(otu.fun.wo.G.17.comp.rel)

melt.otu.fun.wo.G.17.comp.rel <- melt(otu.fun.wo.G.17.comp.rel)
head(melt.otu.fun.wo.G.17.comp.rel)
melt.otu.fun.wo.G.17.comp.rel$Days <- ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "106"), "106", 
                                          ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "120"),"120",
                                                 ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "62"), "62",
                                                        ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "76"),"76",
                                                               ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "141"),"141", 
                                                                      ifelse(str_detect(melt.otu.fun.wo.G.17.comp.rel$variable, "48"),"48", "90"))))))


meta.fun.wo.G.17<-sample_data(fun.clean.ss.wo.G.17)
meta.fun.wo.G.17<-data.frame(meta.fun.wo.G.17)


for (i in as.character(meta.fun.wo.G.17$SampleID)) {
  if (i %in% as.character(melt.otu.fun.wo.G.17.comp.rel$variable) == T){
    melt.otu.fun.wo.G.17.comp.rel$Section[which(melt.otu.fun.wo.G.17.comp.rel$variable == i)] <- as.character(meta.fun.wo.G.17$Microhabitat[which(meta.fun.wo.G.17$SampleID == i)])
  }
}

head(melt.otu.fun.wo.G.17.comp.rel)

otu.fun.wo.G.17.comp.rel.sum.tab <- melt.otu.fun.wo.G.17.comp.rel %>% group_by(OTU, Days, Section) %>% summarise(Mean = mean(value))

### Assign origin
otu.fun.wo.G.17.comp.rel.sum.tab$Class <- ifelse(otu.fun.wo.G.17.comp.rel.sum.tab$OTU %in% persistent.fun, "Persistent",
                                              ifelse(otu.fun.wo.G.17.comp.rel.sum.tab$OTU %in% transient.fun.soil, "Transient_soil",
                                                     ifelse(otu.fun.wo.G.17.comp.rel.sum.tab$OTU %in% transient.fun.nonsoil,"Transient_nonsoil", 
                                                            ifelse(otu.fun.wo.G.17.comp.rel.sum.tab$OTU %in% soil.fun,"Soil","Unknown"))))

otu.fun.wo.G.17.comp.rel.sum.tab.2 <- otu.fun.wo.G.17.comp.rel.sum.tab %>% group_by(Days, Class, Section) %>% summarise(Total = sum(Mean))

otu.fun.wo.G.17.comp.rel.sum.tab.2$Days <- factor(otu.fun.wo.G.17.comp.rel.sum.tab.2$Days, levels = c("48","62","76","90","106","120","141"))
otu.fun.wo.G.17.comp.rel.sum.tab.2$Class <- factor(otu.fun.wo.G.17.comp.rel.sum.tab.2$Class, levels = c("Unknown","Soil","Transient_soil","Transient_nonsoil","Persistent"))
otu.fun.wo.G.17.comp.rel.sum.tab.2$Section <- factor(otu.fun.wo.G.17.comp.rel.sum.tab.2$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

### plotting
p <- ggplot(otu.fun.wo.G.17.comp.rel.sum.tab.2, aes(x=Days, y = Total, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() + 
  facet_wrap(~Section, nrow = 2)+
  scale_fill_manual(values = my_color_collection) +
  theme(aspect.ratio = 1.2)+
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

p

write.csv(otu.fun.wo.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_fungi_other_sections.csv")
write.csv(otu.bac.wo.G.comp.rel.sum.tab.2,"Raw data for RA of persistent and transient_bacteria_other_sections.csv")


#### Taxonomic composition
#Persistent OTUs - Stacked ridge plot or stacked area graph
## Assign taxonomic information
otu.bac.G.comp.rel.sum.tab

bac.list$Genus <- as.character(bac.list$Genus)
bac.list$Genus[is.na(bac.list$Genus)] <- "Unidentified"


##Seeds
otu.bac.G.comp.rel.sum.tab$Genus <- 0
for (i in as.character(bac.list$OTU)) {
  if (i %in% as.character(otu.bac.G.comp.rel.sum.tab$OTU) == T){
    otu.bac.G.comp.rel.sum.tab$Genus[which(otu.bac.G.comp.rel.sum.tab$OTU == i)] <- as.character(bac.list$Genus[which(bac.list$OTU == i)])
  }
}


otu.bac.G.comp.rel.sum.tab.per <- subset(otu.bac.G.comp.rel.sum.tab, Class =="Persistent")
class(otu.bac.G.comp.rel.sum.tab.per$Days)

otu.bac.G.comp.rel.sum.tab.per$Days <- factor(otu.bac.G.comp.rel.sum.tab.per$Days, levels = c("0","76","80","90","106","120","140","141"))

otu.bac.G.comp.rel.sum.tab.per.gen <-otu.bac.G.comp.rel.sum.tab.per%>% group_by(Days, Year, Genus) %>% summarise(Gen_mean = sum(Mean))

otu.bac.G.comp.rel.sum.tab.per.gen$Genus <- factor(otu.bac.G.comp.rel.sum.tab.per.gen$Genus, levels = genus.ordering.new)

write.csv(otu.bac.G.comp.rel.sum.tab.per,"RA of persistent OTUs_in seeds_bacteria.csv")

ggplot(otu.bac.G.comp.rel.sum.tab.per.gen, aes(x=Days, y=Gen_mean, fill=Genus)) + 
  geom_bar(stat = "identity",alpha=0.7, colour="white", position = 'stack') +
  scale_fill_viridis(discrete = T) +
  facet_wrap(~Year, scales ="free")+theme(aspect.ratio = 1.2)+
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())



##Other compartment
otu.bac.wo.G.comp.rel.sum.tab$Genus <- 0
for (i in as.character(bac.list$OTU)) {
  if (i %in% as.character(otu.bac.wo.G.comp.rel.sum.tab$OTU) == T){
    otu.bac.wo.G.comp.rel.sum.tab$Genus[which(otu.bac.wo.G.comp.rel.sum.tab$OTU == i)] <- as.character(bac.list$Genus[which(bac.list$OTU == i)])
  }
}


otu.bac.wo.G.comp.rel.sum.tab.per <- subset(otu.bac.wo.G.comp.rel.sum.tab, Class =="Persistent")
otu.bac.wo.G.comp.rel.sum.tab.per$Days <- as.numeric(as.character(otu.bac.wo.G.comp.rel.sum.tab.per$Days))

otu.bac.wo.G.comp.rel.sum.tab.per$Section <- factor(otu.bac.wo.G.comp.rel.sum.tab.per$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

otu.bac.wo.G.comp.rel.sum.tab.per <- otu.bac.wo.G.comp.rel.sum.tab.per%>% group_by(Days, Section, Genus) %>% summarise(Gen_mean = sum(Mean))

genus.ordering <- otu.bac.wo.G.comp.rel.sum.tab.per%>% group_by(Genus) %>% summarise(Total_value = sum(Gen_mean))%>% arrange(desc(Total_value))
genus.ordering <- genus.ordering$Genus
genus.ordering.wo.unid <- genus.ordering[which(genus.ordering != "Unidentified")]
genus.ordering.new <- c(genus.ordering.wo.unid, "Unidentified")
genus.ordering.new <- rev(genus.ordering.new)
otu.bac.wo.G.comp.rel.sum.tab.per$Genus <- factor(otu.bac.wo.G.comp.rel.sum.tab.per$Genus, levels = genus.ordering.new)

write.csv(otu.bac.wo.G.comp.rel.sum.tab.per,"RA of persistent OTUs_in other sections_bacteria.csv")

# Library
library(viridis)
library(hrbrthemes)

# Plot
ggplot(otu.bac.wo.G.comp.rel.sum.tab.per, aes(x=Days, y=Gen_mean, fill=Genus)) + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  facet_wrap(~Section, nrow =2)+theme(aspect.ratio = 1.2)+
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

###Fungi

fun.list$Genus <- as.character(fun.list$Genus)
fun.list$Genus[is.na(fun.list$Genus)] <- "unidentified"

otu.fun.wo.G.comp.rel.sum.tab$Genus <- 0
for (i in as.character(fun.list$OTU)) {
  if (i %in% as.character(otu.fun.wo.G.comp.rel.sum.tab$OTU) == T){
    otu.fun.wo.G.comp.rel.sum.tab$Genus[which(otu.fun.wo.G.comp.rel.sum.tab$OTU == i)] <- as.character(fun.list$Genus[which(fun.list$OTU == i)])
  }
}

otu.fun.wo.G.comp.rel.sum.tab.per <- subset(otu.fun.wo.G.comp.rel.sum.tab, Class =="Persistent")
otu.fun.wo.G.comp.rel.sum.tab.per$Days <- as.numeric(as.character(otu.fun.wo.G.comp.rel.sum.tab.per$Days))

otu.fun.wo.G.comp.rel.sum.tab.per$Section <- factor(otu.fun.wo.G.comp.rel.sum.tab.per$Section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

otu.fun.wo.G.comp.rel.sum.tab.per <- otu.fun.wo.G.comp.rel.sum.tab.per%>% group_by(Days, Section, Genus) %>% summarise(Gen_mean = sum(Mean))

genus.ordering <- otu.fun.wo.G.comp.rel.sum.tab.per%>% group_by(Genus) %>% summarise(Total_value = sum(Gen_mean))%>% arrange(desc(Total_value))
genus.ordering <- genus.ordering$Genus
genus.ordering.wo.unid <- genus.ordering[which(genus.ordering != "unidentified")]
genus.ordering.new <- c(genus.ordering.wo.unid, "unidentified")
genus.ordering.new <- rev(genus.ordering.new)
otu.fun.wo.G.comp.rel.sum.tab.per$Genus <- factor(otu.fun.wo.G.comp.rel.sum.tab.per$Genus, levels = genus.ordering.new)

write.csv(otu.fun.wo.G.comp.rel.sum.tab.per,"RA of persistent OTUs_in other sections_fungi.csv")

# Plot
ggplot(otu.fun.wo.G.comp.rel.sum.tab.per, aes(x=Days, y=Gen_mean, fill=Genus)) + 
  geom_area(alpha=0.6 , size=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  facet_wrap(~Section, nrow =2)+theme(aspect.ratio = 1.2)+
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())




otu.fun.G.comp.rel.sum.tab$Genus <- 0
for (i in as.character(fun.list$OTU)) {
  if (i %in% as.character(otu.fun.G.comp.rel.sum.tab$OTU) == T){
    otu.fun.G.comp.rel.sum.tab$Genus[which(otu.fun.G.comp.rel.sum.tab$OTU == i)] <- as.character(fun.list$Genus[which(fun.list$OTU == i)])
  }
}

otu.fun.G.comp.rel.sum.tab.per <- subset(otu.fun.G.comp.rel.sum.tab, Class =="Persistent")
class(otu.fun.G.comp.rel.sum.tab.per$Days)

otu.fun.G.comp.rel.sum.tab.per$Days <- factor(otu.fun.G.comp.rel.sum.tab.per$Days, levels = c("0","76","80","90","106","120","140","141"))

otu.fun.G.comp.rel.sum.tab.per.gen <-otu.fun.G.comp.rel.sum.tab.per%>% group_by(Days, Year, Genus) %>% summarise(Gen_mean = sum(Mean))

otu.fun.G.comp.rel.sum.tab.per.gen$Genus <- factor(otu.fun.G.comp.rel.sum.tab.per.gen$Genus, levels = genus.ordering.new)

write.csv(otu.fun.G.comp.rel.sum.tab.per,"RA of persistent OTUs_in seeds_fungi.csv")

ggplot(otu.fun.G.comp.rel.sum.tab.per.gen, aes(x=Days, y=Gen_mean, fill=Genus)) + 
  geom_bar(stat = "identity",alpha=0.7, colour="white", position = 'stack') +
  scale_fill_viridis(discrete = T) +
  facet_wrap(~Year, scales ="free")+theme(aspect.ratio = 1.2)+
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())




###Statistical analysis
### PERMANOVA on Jaccard distance
#Age, Compartment, maternal plant, Genotype

###Subset 2017 and 2018
###Maternal effect in developing seeds
bac.clean.ss.f.17 <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
bac.clean.ss.f.17 <- phyloseq::filter_taxa(bac.clean.ss.f.17, function(x) sum(x) != 0, TRUE)
# Compute Jaccard distance on Presence/absence matrix
dJ.2017.b <- phyloseq::distance(bac.clean.ss.f.17, method = "jaccard", binary = T) 

# Get sample informations
metadata.b <- as(sample_data(bac.clean.ss.f.17), "data.frame") 
head(metadata.b)

metadata.b$motherID <- paste0(metadata.b$Year,"_",metadata.b$Plot,"_",metadata.b$Days)

# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2017.b <- adonis(dJ.2017.b ~ Age + Genotype+motherID, data = metadata.b , 
                      perm = 9999))


fun.clean.ss.f.17 <- subset_samples(fun.clean.ss.f, Year == "year_2017"& Compartment == "Seed")
fun.clean.ss.f.17 <- phyloseq::filter_taxa(fun.clean.ss.f.17, function(x) sum(x) != 0, TRUE)
# Compute Jaccard distance on Presence/absence matrix
dJ.2017.f <- phyloseq::distance(fun.clean.ss.f.17, method = "jaccard", binary = T) 

# Get sample informations
metadata.f <- as(sample_data(fun.clean.ss.f.17), "data.frame") 
metadata.f$motherID <- paste0(metadata.f$Year,"_",metadata.f$Plot,"_",metadata.f$Days)
head(metadata.f)
# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2017.f <- adonis(dJ.2017.f ~ Age + Genotype+motherID, data = metadata.f , 
                      perm = 9999))


##2018 samples
# from 76 days to 141 days

##Bacteria
bac.clean.ss.f.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.f.18.Gs <- subset_samples(bac.clean.ss.f.18, Days %in% c("76","90","106","120","141"))
bac.clean.ss.f.18.Gs <- phyloseq::filter_taxa(bac.clean.ss.f.18.Gs, function(x) sum(x) != 0, TRUE)


# Compute Jaccard distance on Presence/absence matrix
dJ.2018.b <- phyloseq::distance(bac.clean.ss.f.18.Gs, method = "jaccard", binary = T) 

# Get sample informations
metadata.b <- as(sample_data(bac.clean.ss.f.18.Gs), "data.frame") 
metadata.b$motherID <- paste0(metadata.b$Year,"_",metadata.b$Plot,"_",metadata.b$Days)
metadata.b$mother_section <- paste0(metadata.b$Year,"_",metadata.b$Plot,"_",metadata.b$Days,"_",metadata.b$Microhabitat)
head(metadata.b)
# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2018.b <- adonis(dJ.2018.b ~ Age+Compartment+Microhabitat+motherID+mother_section, data = metadata.b , 
                      perm = 9999))

###only seeds
bac.clean.ss.f.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018" & Microhabitat == "G")
bac.clean.ss.f.18.Gs <- subset_samples(bac.clean.ss.f.18, Days %in% c("76","90","106","120","141"))
bac.clean.ss.f.18.Gs <- phyloseq::filter_taxa(bac.clean.ss.f.18.Gs, function(x) sum(x) != 0, TRUE)


# Compute Jaccard distance on Presence/absence matrix
dJ.2018.b <- phyloseq::distance(bac.clean.ss.f.18.Gs, method = "jaccard", binary = T) 

# Get sample informations
metadata.b <- as(sample_data(bac.clean.ss.f.18.Gs), "data.frame") 
metadata.b$motherID <- paste0(metadata.b$Year,"_",metadata.b$Plot,"_",metadata.b$Days)
head(metadata.b)
# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2018.b <- adonis(dJ.2018.b ~ Age+motherID, data = metadata.b , 
                      perm = 9999))



###Fungi
fun.clean.ss.f.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.f.18.Gs <- subset_samples(fun.clean.ss.f.18, Days %in% c("76","90","106","120","141"))
fun.clean.ss.f.18.Gs <- phyloseq::filter_taxa(fun.clean.ss.f.18.Gs, function(x) sum(x) != 0, TRUE)


# Compute Jaccard distance on Presence/absence matrix
dJ.2018.f <- phyloseq::distance(fun.clean.ss.f.18.Gs, method = "jaccard", binary = T) 

# Get sample informations
metadata.f <- as(sample_data(fun.clean.ss.f.18.Gs), "data.frame") 
metadata.f$motherID <- paste0(metadata.f$Year,"_",metadata.f$Plot,"_",metadata.f$Days)
metadata.f$mother_section <- paste0(metadata.f$Year,"_",metadata.f$Plot,"_",metadata.f$Days,"_",metadata.f$Microhabitat)
head(metadata.f)
# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2018.f <- adonis(dJ.2018.f ~ Age+Compartment+Microhabitat+motherID+mother_section, data = metadata.f , 
                      perm = 9999))

###only seeds
fun.clean.ss.f.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018" & Microhabitat == "G")
fun.clean.ss.f.18.Gs <- subset_samples(fun.clean.ss.f.18, Days %in% c("76","90","106","120","141"))
fun.clean.ss.f.18.Gs <- phyloseq::filter_taxa(fun.clean.ss.f.18.Gs, function(x) sum(x) != 0, TRUE)


# Compute Jaccard distance on Presence/absence matrix
dJ.2018.f <- phyloseq::distance(fun.clean.ss.f.18.Gs, method = "jaccard", binary = T) 

# Get sample informations
metadata.f <- as(sample_data(fun.clean.ss.f.18.Gs), "data.frame") 
metadata.f$motherID <- paste0(metadata.f$Year,"_",metadata.f$Plot,"_",metadata.f$Days)
head(metadata.f)
# Run PERMANOVA on the Matrix built using jaccard  distance : Acorn position and site effects on microbiota composition
(aov.2018.f <- adonis(dJ.2018.f ~ Age+motherID, data = metadata.f , 
                      perm = 9999))

