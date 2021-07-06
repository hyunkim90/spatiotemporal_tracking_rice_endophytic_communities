##Possibility of seed to seed movement

##merge files
bac.clean.ss.S.G <- subset_samples(bac.clean.ss.f, Compartment %in% c("Stem", "Seed") & Location == "Suwon")
bac.clean.ss.S.G <- phyloseq::filter_taxa(bac.clean.ss.S.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S.G <- subset_samples(fun.clean.ss2, Compartment %in% c("Stem", "Seed","Grain") & Location == "Suwon")
fun.clean.ss.S.G <- phyloseq::filter_taxa(fun.clean.ss.S.G, function(x) sum(x) != 0, TRUE)
###Find whether core OTUs exist or not
bac.clean.ss.S.G.rep <-  merge_samples(bac.clean.ss.S.G, "Replication")
t(otu_table(bac.clean.ss.S.G.rep))

bac.clean.ss.S.G.rep.rel <- microbiome::transform(bac.clean.ss.S.G.rep, "compositional")

core.bac.90.all <- core_members(bac.clean.ss.S.G.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.all

core.bac.80.all <- core_members(bac.clean.ss.S.G.rep.rel, detection = 0, prevalence = 80/100)
core.bac.80.all

core.bac.80.tax<-subset(bac.list, OTU %in% core.bac.80.all)



### Fungi
fun.clean.ss.S.G.rep <-  merge_samples(fun.clean.ss.S.G, "Replication")
t(otu_table(fun.clean.ss.S.G.rep))
fun.clean.ss.S.G.rep.rel <- microbiome::transform(fun.clean.ss.S.G.rep, "compositional")

core.fun.90.all <- core_members(fun.clean.ss.S.G.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.all

core.fun.80.all <- core_members(fun.clean.ss.S.G.rep.rel, detection = 0, prevalence = 80/100)
core.fun.80.all

core.fun.80.tax<-subset(fun.list, OTU %in% core.fun.80.all)

## Seed and highest stem
# unique(subset(merged.map, Replication %in% c('UF1G80', "UF1G120", 'UF1G140', 'S3_48', 'S4_62','S6_76','S9_90','S9_106','S9_120','S9_141', 'G_0', 'G_76', 'G_90','G_106','G_120','G_141'))$Replication)
# 
# bac.clean.ss.S.G.rep.part <- subset_samples(bac.clean.ss.S.G.rep, Replication %in% c('56', "57", '58', '24', '31','42','55','52','53','54', '1', '2', '3','4','5','6'))
# 
# 
# 
# bac.clean.ss.S.G.rep.part.rel <- microbiome::transform(bac.clean.ss.S.G.rep.part, "compositional")
# t(otu_table(bac.clean.ss.S.G.rep.part.rel))
# 
# core.bac.90 <- core_members(bac.clean.ss.S.G.rep.part.rel, detection = 0, prevalence = 90/100)
# core.bac.90
# 
# core.bac.80 <- core_members(bac.clean.ss.S.G.rep.part.rel, detection = 0, prevalence = 80/100)
# core.bac.80
# 
# tax.table.core.bac <- subset(tax.table.all, rownames(tax.table.all)%in%core.bac.80)




### Stem vs Seed
### Bacteria
###Stem core
bac.clean.ss.S <- subset_samples(bac.clean.ss.S.G, Compartment =="Stem")
bac.clean.ss.S <- phyloseq::filter_taxa(bac.clean.ss.S, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S.rep <- merge_samples(bac.clean.ss.S, "Replication")

bac.clean.ss.S.rep.rel <- microbiome::transform(bac.clean.ss.S.rep, "compositional")

core.bac.90.stem <- core_members(bac.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.stem

core.bac.80.stem <- core_members(bac.clean.ss.S.rep.rel, detection = 0, prevalence = 80/100)
core.bac.80.stem
core.bac.80.stem.tax<-subset(bac.list, OTU %in% core.bac.80.stem)
###Seed core
bac.clean.ss.G <- subset_samples(bac.clean.ss.S.G, Compartment =="Seed"|Compartment =="Grain")
bac.clean.ss.G <- phyloseq::filter_taxa(bac.clean.ss.G, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G.rep <- merge_samples(bac.clean.ss.G, "Replication")

bac.clean.ss.G.rep.rel <- microbiome::transform(bac.clean.ss.G.rep, "compositional")

core.bac.90.grain <- core_members(bac.clean.ss.G.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.grain

core.bac.80.grain <- core_members(bac.clean.ss.G.rep.rel, detection = 0, prevalence = 80/100)
core.bac.80.grain
core.bac.80.grain.tax<-subset(bac.list, OTU %in% core.bac.80.grain)

### Common core
intersect(core.bac.80.stem, core.bac.80.grain)

core.bac.80.commonseedstem.tax<-subset(bac.list, OTU %in% intersect(core.bac.80.stem, core.bac.80.grain))

core.bac.80.grain.tax
core.bac.80.stem.tax


### Fungi
###Stem core
fun.clean.ss.S <- subset_samples(fun.clean.ss.S.G, Compartment =="Stem")
fun.clean.ss.S <- phyloseq::filter_taxa(fun.clean.ss.S, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S.rep <- merge_samples(fun.clean.ss.S, "Replication")

fun.clean.ss.S.rep.rel <- microbiome::transform(fun.clean.ss.S.rep, "compositional")

core.fun.90.stem <- core_members(fun.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.stem

core.fun.80.stem <- core_members(fun.clean.ss.S.rep.rel, detection = 0, prevalence = 80/100)
core.fun.80.stem

core.fun.80.stem.tax<-subset(fun.list, OTU %in% core.fun.80.stem)
###Seed core
fun.clean.ss.G <- subset_samples(fun.clean.ss.S.G, Compartment =="Seed"|Compartment =="Grain")
fun.clean.ss.G <- phyloseq::filter_taxa(fun.clean.ss.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G.rep <- merge_samples(fun.clean.ss.G, "Replication")

fun.clean.ss.G.rep.rel <- microbiome::transform(fun.clean.ss.G.rep, "compositional")

core.fun.90.grain <- core_members(fun.clean.ss.G.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.grain

core.fun.80.grain <- core_members(fun.clean.ss.G.rep.rel, detection = 0, prevalence = 80/100)
core.fun.80.grain
core.fun.80.grain.tax<-subset(fun.list, OTU %in% core.fun.80.grain)

### Common core
intersect(core.fun.80.stem, core.fun.80.grain)

core.fun.80.commonseedstem.tax<-subset(fun.list, OTU %in% intersect(core.fun.80.stem, core.fun.80.grain))

core.fun.80.grain.tax
core.fun.80.stem.tax

### 
bac.clean.ss.G.17 <- subset_samples(bac.clean.ss.G, Replication == "UF1G140")
bac.clean.ss.G.17 <- phyloseq::filter_taxa(bac.clean.ss.G.17, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G.0 <- subset_samples(bac.clean.ss.G, Replication == "G_0")
bac.clean.ss.G.0 <- phyloseq::filter_taxa(bac.clean.ss.G.0, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G.141 <- subset_samples(bac.clean.ss.G, Replication == "G_141")
bac.clean.ss.G.141 <- phyloseq::filter_taxa(bac.clean.ss.G.141, function(x) sum(x) != 0, TRUE)

#### Binary comparison
intersect(taxa_names(bac.clean.ss.G.17), taxa_names(bac.clean.ss.G.0)) #44/76, 44/84
intersect(taxa_names(bac.clean.ss.G.17), taxa_names(bac.clean.ss.G.141)) #37/76, 37/70
intersect(taxa_names(bac.clean.ss.G.0), taxa_names(bac.clean.ss.G.141))#36/84 , 36/70
## 29 OTUs were commonly distributed in 2017 and 2018 seed
common.seedOTUs<-Reduce(intersect, list(taxa_names(bac.clean.ss.G.17),taxa_names(bac.clean.ss.G.0),taxa_names(bac.clean.ss.G.141)))

tax.table.all.common.seed <- subset(bac.list, OTU %in% common.seedOTUs)
write.table(tax.table.all.common.seed, 'Taxonomy of common seed bacterial OTUs.tsv', sep = '\t', quote =F)
### Distribution of each OTU in plant tissue (2018) using ridge plot
get_df_ridge <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% core.bac.90)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  str(df.ridge)
  
  ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Verrucomicrobia","Alphaproteobacteria","Deltaproteobacteria","Gammaproteobacteria","Chloroflexi","Acidobacteria")
  df.ridge$ColPhylum <- 0
  df.ridge$ColPhylum[-which(df.ridge$Phylum %in% ColPhylum)] <- "Other"
  df.ridge$ColPhylum[which(df.ridge$Class =="Alphaproteobacteria")] <- "Alphaproteobacteria"
  df.ridge$ColPhylum[which(df.ridge$Class =="Gammaproteobacteria")] <- "Gammaproteobacteria"
  df.ridge$ColPhylum[which(df.ridge$Class =="Deltaproteobacteria")] <- "Deltaproteobacteria"
  df.ridge$ColPhylum[which(df.ridge$Phylum =="Bacteroidetes")] <- "Bacteroidetes"
  df.ridge$ColPhylum[which(df.ridge$Phylum == "Actinobacteria")] <- "Actinobacteria"
  df.ridge$ColPhylum[which(df.ridge$Phylum == "Verrucomicrobia")] <- "Verrucomicrobia"
  df.ridge$ColPhylum[which(df.ridge$Phylum == "Firmicutes")] <- "Firmicutes"
  df.ridge$ColPhylum[which(df.ridge$Phylum == "Acidobacteria")] <- "Acidobacteria"
  df.ridge$ColPhylum[which(df.ridge$Phylum == "Chloroflexi")] <- "Chloroflexi" 
  
  df.ridge$ColPhylum <- factor(df.ridge$ColPhylum, levels = c("Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",
                                                              "Chloroflexi","Firmicutes","Acidobacteria", "Actinobacteria","Bacteroidetes","Verrucomicrobia","Other"))
  
  
  return(df.ridge)
}


library(ggridges)

compartment = "Leaf"
df.ridge_2.bac.core.leaf<- get_df_ridge(bac.clean.ss.18)
unique(df.ridge_2.bac.core.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.bac.core.leaf$Sample <- factor(df.ridge_2.bac.core.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.bac.core.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +
  theme(legend.text=element_text(size=12)) + scale_fill_manual(values=Palette.phylum)+
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none") +theme(aspect.ratio = 2)



compartment = "Stem"

df.ridge_2.bac.core.stem<- get_df_ridge(bac.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.bac.core.stem$Sample <- factor(df.ridge_2.bac.core.stem$Sample, levels = stem.sample)

ggplot(df.ridge_2.bac.core.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +scale_fill_manual(values=Palette.phylum)+
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none") +theme(aspect.ratio = 2)



compartment = "Root" 
df.ridge_2.bac.core.root<- get_df_ridge(bac.clean.ss.18)

root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')

df.ridge_2.bac.core.root$Sample <- factor(df.ridge_2.bac.core.root$Sample, levels = root.sample)


library(ggridges)

ggplot(df.ridge_2.bac.core.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +scale_fill_manual(values=Palette.phylum)+
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none") +theme(aspect.ratio = 2)


compartment = "Seed" 
df.ridge_2.bac.core.seed<- get_df_ridge(bac.clean.ss.18)

seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')

df.ridge_2.bac.core.seed$Sample <- factor(df.ridge_2.bac.core.seed$Sample, levels = seed.sample)


ggplot(df.ridge_2.bac.core.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +
  theme(legend.text=element_text(size=12)) + scale_fill_manual(values=Palette.phylum)+
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none") + theme(aspect.ratio = 2)



####Fungi
fun.clean.ss.G.17 <- subset_samples(fun.clean.ss.G, Replication == "UF1G140")
fun.clean.ss.G.17 <- phyloseq::filter_taxa(fun.clean.ss.G.17, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G.0 <- subset_samples(fun.clean.ss.G, Replication == "G_0")
fun.clean.ss.G.0 <- phyloseq::filter_taxa(fun.clean.ss.G.0, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G.141 <- subset_samples(fun.clean.ss.G, Replication == "G_141")
fun.clean.ss.G.141 <- phyloseq::filter_taxa(fun.clean.ss.G.141, function(x) sum(x) != 0, TRUE)

#### Binary comparison
intersect(taxa_names(fun.clean.ss.G.17), taxa_names(fun.clean.ss.G.0)) #40/63, 40/121
intersect(taxa_names(fun.clean.ss.G.17), taxa_names(fun.clean.ss.G.141)) #37/63, 37/133
intersect(taxa_names(fun.clean.ss.G.0), taxa_names(fun.clean.ss.G.141))#63/121 , 63/133

tax.table.2018.common.seed.f <- subset(fun.list, OTU %in% intersect(taxa_names(fun.clean.ss.G.0), taxa_names(fun.clean.ss.G.141)))
## 33 OTUs were commonly distributed in 2017 and 2018 seed
common.seedOTUs.fun<-Reduce(intersect, list(taxa_names(fun.clean.ss.G.17),taxa_names(fun.clean.ss.G.0),taxa_names(fun.clean.ss.G.141)))

tax.table.all.common.seed.f <- subset(fun.list, OTU %in% common.seedOTUs.fun)
write.table(tax.table.all.common.seed.f, 'Taxonomy of common seed fungal OTUs.tsv', sep = '\t', quote =F)
### Distribution of each OTU in plant tissue (2018) using ridge plot
get_df_ridge.f <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
 
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% core.fun.90)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  str(df.ridge)
  ColClass<-c("Tremellomycetes","Sordariomycetes","Leotiomycetes","Dothideomycetes","Agaricomycetes","Mortierellomycetes","Ustilaginomycetes","Eurotiomycetes","Saccharomycetes")
  df.ridge$ColClass <- 0
  df.ridge$ColClass[-which(df.ridge$Class %in% ColClass)] <- "Other"
  df.ridge$ColClass[which(df.ridge$Class =="Tremellomycetes")] <- "Tremellomycetes"
  df.ridge$ColClass[which(df.ridge$Class =="Sordariomycetes")] <- "Sordariomycetes"
  df.ridge$ColClass[which(df.ridge$Class =="Leotiomycetes")] <- "Leotiomycetes"
  df.ridge$ColClass[which(df.ridge$Class =="Dothideomycetes")] <- "Dothideomycetes"
  df.ridge$ColClass[which(df.ridge$Class == "Agaricomycetes")] <- "Agaricomycetes"
  df.ridge$ColClass[which(df.ridge$Class == "Mortierellomycetes")] <- "Mortierellomycetes"
  df.ridge$ColClass[which(df.ridge$Class == "Ustilaginomycetes")] <- "Ustilaginomycetes"
  df.ridge$ColClass[which(df.ridge$Class == "Eurotiomycetes")] <- "Eurotiomycetes"
  df.ridge$ColClass[which(df.ridge$Class == "Saccharomycetes")] <- "Saccharomycetes" 
  
  df.ridge$ColClass <- factor(df.ridge$ColClass, levels = c("Sordariomycetes", "Dothideomycetes", "Leotiomycetes",
                                                          "Eurotiomycetes","Saccharomycetes","Tremellomycetes", "Agaricomycetes","Ustilaginomycetes","Mortierellomycetes","Other"))

  return(df.ridge)
}


library(ggridges)

compartment = "Leaf"
df.ridge_2.fun.core.leaf<- get_df_ridge.f(fun.clean.ss.18)
unique(df.ridge_2.fun.core.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.fun.core.leaf$Sample <- factor(df.ridge_2.fun.core.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.fun.core.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +scale_fill_manual(values=Palette.class)+
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")  +theme(aspect.ratio = 2)



compartment = "Stem"

df.ridge_2.fun.core.stem<- get_df_ridge.f(fun.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.fun.core.stem$Sample <- factor(df.ridge_2.fun.core.stem$Sample, levels = stem.sample)

ggplot(df.ridge_2.fun.core.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
  # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
  geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
  #xlab('\n Diagnosis')+
  ylab("Relative abundance (%) \n") +scale_fill_manual(values=Palette.class)+
  theme(legend.text=element_text(size=12)) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(size=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
  theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
  theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
  #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position = "none")  +theme(aspect.ratio = 2)

  
  
  compartment = "Root" 
  df.ridge_2.fun.core.root<- get_df_ridge.f(fun.clean.ss.18)
  
  root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')
  
  df.ridge_2.fun.core.root$Sample <- factor(df.ridge_2.fun.core.root$Sample, levels = root.sample)

  
  ggplot(df.ridge_2.fun.core.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
    #xlab('\n Diagnosis')+
    ylab("Relative abundance (%) \n") +
    theme(legend.text=element_text(size=12)) + scale_fill_manual(values=Palette.class)+
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    theme(legend.position = "none")  +theme(aspect.ratio = 2)
  
  
 compartment = "Grain" 
 df.ridge_2.fun.core.seed<- get_df_ridge.f(fun.clean.ss.18)
 
 seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')
 
 df.ridge_2.fun.core.seed$Sample <- factor(df.ridge_2.fun.core.seed$Sample, levels = seed.sample)
 
 
 ggplot(df.ridge_2.fun.core.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
   # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
   geom_density_ridges2(stat = "identity", scale=10, color='white',size=0.5, alpha = 0.7)+
   #xlab('\n Diagnosis')+
   ylab("Relative abundance (%) \n") +scale_fill_manual(values=Palette.class)+
   theme(legend.text=element_text(size=12)) + 
   # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
   #theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
   theme(legend.title=element_blank()) +
   guides(colour = guide_legend(override.aes = list(size=8)))+
   guides(size=FALSE) +
   theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
   theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
   theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
   theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
   #scale_fill_manual(labels = c('Control'="Control",'RA'= "RA",'Non-differential'='Non-differential'), values = c("Control"= "#6699CC", 'RA'='#CC9900',"Non-differential"= "light grey"))+
   theme(panel.grid.major = element_blank())+
   theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
   theme(legend.position = "none") +theme(aspect.ratio = 2)
 
 
 
