### Origin of non-core OTUs
sample_data(bac.clean.ss.18)
bac.clean.ss.18.L <- subset_samples(bac.clean.ss.18, Compartment == "Leaf")
bac.clean.ss.18.L <- phyloseq::filter_taxa(bac.clean.ss.18.L, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.R <- subset_samples(bac.clean.ss.18, Compartment == "Root")
bac.clean.ss.18.R <- phyloseq::filter_taxa(bac.clean.ss.18.R, function(x) sum(x) != 0, TRUE)

## lost OTU 
non_herited.bac<- taxa_names(bac.clean.ss.G.0)[-which(taxa_names(bac.clean.ss.G.0)%in% intersect(taxa_names(bac.clean.ss.G.0), taxa_names(bac.clean.ss.G.141)))]
tax.non_herited.bac <- subset(bac.list, OTU %in% non_herited.bac)

## Gained OTU
non_herited.bac.2<- taxa_names(bac.clean.ss.G.141)[-which(taxa_names(bac.clean.ss.G.141)%in% intersect(taxa_names(bac.clean.ss.G.0), taxa_names(bac.clean.ss.G.141)))]
tax.non_herited.bac.2 <- subset(bac.list, OTU %in% non_herited.bac.2)


## Distribution of lost and gained OTU

get_df_ridge.lost <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% non_herited.bac)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  
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
  
  
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  str(df.ridge)
  
  return(df.ridge)
}


library(ggridges)

compartment = "Leaf"
df.ridge_2.bac.lost.leaf<- get_df_ridge.lost(bac.clean.ss.18)
unique(df.ridge_2.bac.lost.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.bac.lost.leaf$Sample <- factor(df.ridge_2.bac.lost.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.bac.lost.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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



compartment = "Stem"

df.ridge_2.bac.lost.stem<- get_df_ridge.lost(bac.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.bac.lost.stem$Sample <- factor(df.ridge_2.bac.lost.stem$Sample, levels = stem.sample)

ggplot(df.ridge_2.bac.lost.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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



compartment = "Root" 
df.ridge_2.bac.lost.root<- get_df_ridge.lost(bac.clean.ss.18)

root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')

df.ridge_2.bac.lost.root$Sample <- factor(df.ridge_2.bac.lost.root$Sample, levels = root.sample)


ggplot(df.ridge_2.bac.lost.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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
df.ridge_2.bac.lost.seed<- get_df_ridge.lost(bac.clean.ss.18)

seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')

df.ridge_2.bac.lost.seed$Sample <- factor(df.ridge_2.bac.lost.seed$Sample, levels = seed.sample)


ggplot(df.ridge_2.bac.lost.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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

### Gained
get_df_ridge.gained <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% non_herited.bac.2)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  
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
  
  
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  str(df.ridge)
  
  return(df.ridge)
}


compartment = "Leaf"
df.ridge_2.bac.gained.leaf<- get_df_ridge.gained(bac.clean.ss.18)
unique(df.ridge_2.bac.gained.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.bac.gained.leaf$Sample <- factor(df.ridge_2.bac.gained.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.bac.gained.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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



compartment = "Stem"

df.ridge_2.bac.gained.stem<- get_df_ridge.gained(bac.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.bac.gained.stem$Sample <- factor(df.ridge_2.bac.gained.stem$Sample, levels = stem.sample)

ggplot(df.ridge_2.bac.gained.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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



compartment = "Root" 
df.ridge_2.bac.gained.root<- get_df_ridge.gained(bac.clean.ss.18)

root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')

df.ridge_2.bac.gained.root$Sample <- factor(df.ridge_2.bac.gained.root$Sample, levels = root.sample)


ggplot(df.ridge_2.bac.gained.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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
df.ridge_2.bac.gained.seed<- get_df_ridge.gained(bac.clean.ss.18)

seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')

df.ridge_2.bac.gained.seed$Sample <- factor(df.ridge_2.bac.gained.seed$Sample, levels = seed.sample)


ggplot(df.ridge_2.bac.gained.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColPhylum)) + 
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

### Fungi
## lost OTU 
lost.otu.fun<- taxa_names(fun.clean.ss.G.0)[-which(taxa_names(fun.clean.ss.G.0)%in% intersect(taxa_names(fun.clean.ss.G.0), taxa_names(fun.clean.ss.G.141)))]
tax.lost.otu.fun <- subset(fun.list, OTU %in% lost.otu.fun)

## Gained OTU
gained.otu.fun<- taxa_names(fun.clean.ss.G.141)[-which(taxa_names(fun.clean.ss.G.141)%in% intersect(taxa_names(fun.clean.ss.G.0), taxa_names(fun.clean.ss.G.141)))]
tax.gained.otu.fun <- subset(fun.list, OTU %in% gained.otu.fun)


get_df_ridge.f.lost <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% lost.otu.fun)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  
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

compartment = "Leaf"
df.ridge_2.fun.lost.leaf<- get_df_ridge.f.lost(fun.clean.ss.18)
unique(df.ridge_2.fun.lost.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.fun.lost.leaf$Sample <- factor(df.ridge_2.fun.lost.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.fun.lost.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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



compartment = "Stem"

df.ridge_2.fun.lost.stem<- get_df_ridge.f.lost(fun.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.fun.lost.stem$Sample <- factor(df.ridge_2.fun.lost.stem$Sample, levels = stem.sample)

ggplot(df.ridge_2.fun.lost.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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
  theme(legend.position = "none") +theme(aspect.ratio = 2)



compartment = "Root" 
df.ridge_2.fun.lost.root<- get_df_ridge.f.lost(fun.clean.ss.18)

root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')

df.ridge_2.fun.lost.root$Sample <- factor(df.ridge_2.fun.lost.root$Sample, levels = root.sample)


ggplot(df.ridge_2.fun.lost.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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


compartment = "Grain" 
df.ridge_2.fun.lost.seed<- get_df_ridge.f.lost(fun.clean.ss.18)

seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')

df.ridge_2.fun.lost.seed$Sample <- factor(df.ridge_2.fun.lost.seed$Sample, levels = seed.sample)


ggplot(df.ridge_2.fun.lost.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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


###Gained OTU
get_df_ridge.f.gained <- function(phy.seqs, compartment){
  phy.seqs.tissue <- subset_samples(phy.seqs, Compartment %in% compartment)
  phy.seqs.m <-merge_samples(phy.seqs.tissue, "Replication")
  df.otu <- phy.seqs.m %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% tax.gained.otu.fun.soil_driven$OTU)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  df.ridge <- merge(df.ridge,otu.list, by ='OTU')
  df.ridge$OTU_id <- factor(df.ridge$OTU_id,levels = rev(otu.list$OTU_id))
  
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

compartment = "Leaf"
df.ridge_2.fun.gained.leaf<- get_df_ridge.f.gained(fun.clean.ss.18)
unique(df.ridge_2.fun.gained.leaf$Sample)
leaf.sample <- c('L1_48','L1_62', 'L1_76', 'L1_90','L1_106','L1_120','L1_141', 
                 'L2_48','L2_62', 'L2_76', 'L2_90','L2_106','L2_120','L2_141',
                 'L3_62', 'L3_76', 'L3_90','L3_106','L3_120','L3_141', 
                 'FL_76', 'FL_90','FL_106','FL_120','FL_141')

df.ridge_2.fun.gained.leaf$Sample <- factor(df.ridge_2.fun.gained.leaf$Sample, levels = leaf.sample)

ggplot(df.ridge_2.fun.gained.leaf, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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



compartment = "Stem"

df.ridge_2.fun.gained.stem<- get_df_ridge.f.gained(fun.clean.ss.18)
stem.sample <- c('S1_48','S1_62', 'S1_76', 'S1_90','S1_106','S1_120','S1_141', 
                 'S2_48','S2_62', 'S2_76', 'S2_90','S2_106','S2_120','S2_141',
                 'S3_48','S3_62', 'S3_76', 'S3_90','S3_106','S3_120','S3_141',
                 'S4_62', 'S4_76', 'S4_90','S4_106','S4_120','S4_141',
                 'S5_76', 'S5_90','S5_106','S5_120','S5_141',
                 'S6_76', 'S6_90','S6_106','S6_120','S6_141',
                 'S7_76', 'S7_90','S7_106','S7_120','S7_141',
                 'S8_90','S8_106','S8_120','S8_141',
                 'S9_90','S9_106','S9_120','S9_141')

df.ridge_2.fun.gained.stem$Sample <- factor(df.ridge_2.fun.gained.stem$Sample, levels = stem.sample)
library(ggridges)
ggplot(df.ridge_2.fun.gained.stem, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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
  theme(legend.position = "none") +theme(aspect.ratio = 2)



compartment = "Root" 
df.ridge_2.fun.gained.root<- get_df_ridge.f.gained(fun.clean.ss.18)

root.sample <- c('R_48','R_62', 'R_76', 'R_90','R_106','R_120','R_141')

df.ridge_2.fun.gained.root$Sample <- factor(df.ridge_2.fun.gained.root$Sample, levels = root.sample)


ggplot(df.ridge_2.fun.gained.root, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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


compartment = "Grain" 
df.ridge_2.fun.gained.seed<- get_df_ridge.f.gained(fun.clean.ss.18)

seed.sample<- c('G_0', 'G_76', 'G_90','G_106','G_120','G_141')

df.ridge_2.fun.gained.seed$Sample <- factor(df.ridge_2.fun.gained.seed$Sample, levels = seed.sample)


ggplot(df.ridge_2.fun.gained.seed, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=ColClass)) + 
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





### 
