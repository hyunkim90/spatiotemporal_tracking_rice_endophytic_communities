
###### 2. make barplot with phyloseq
### Bacterial communities
##Sorted by Compartment
bac.clean.ss.17 <- subset_samples(bac.clean.ss.f , Year == "year_2017")
bac.clean.ss.18 <- subset_samples(bac.clean.ss.f , Year == "year_2018")

bac.clean.ss.17.leaf <- subset_samples(bac.clean.ss.17, Compartment == "Leaf")
bac.clean.ss.17.leaf <- phyloseq::filter_taxa(bac.clean.ss.17.leaf, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.stem <- subset_samples(bac.clean.ss.17, Compartment == "Stem")
bac.clean.ss.17.stem <- phyloseq::filter_taxa(bac.clean.ss.17.stem, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.root <- subset_samples(bac.clean.ss.17, Compartment == "Root")
bac.clean.ss.17.root <- phyloseq::filter_taxa(bac.clean.ss.17.root, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.seed <- subset_samples(bac.clean.ss.17, Compartment == "Seed")
bac.clean.ss.17.seed <- phyloseq::filter_taxa(bac.clean.ss.17.seed, function(x) sum(x) != 0, TRUE)

bac.clean.ss.17.soil <- subset_samples(bac.clean.ss.17, Compartment == "Soil")
bac.clean.ss.17.soil <- phyloseq::filter_taxa(bac.clean.ss.17.soil, function(x) sum(x) != 0, TRUE)

#Soil
bac.clean.ss.17.soil.SW <- subset_samples(bac.clean.ss.17.soil, Location == "Suwon")
bac.clean.ss.17.soil.SW <- phyloseq::filter_taxa(bac.clean.ss.17.soil.SW, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.soil.CC1 <- subset_samples(bac.clean.ss.17.soil, Location == "Chuncheon1")
bac.clean.ss.17.soil.CC1 <- phyloseq::filter_taxa(bac.clean.ss.17.soil.CC1, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.soil.CC2 <- subset_samples(bac.clean.ss.17.soil, Location == "Chuncheon2")
bac.clean.ss.17.soil.CC2 <- phyloseq::filter_taxa(bac.clean.ss.17.soil.CC2, function(x) sum(x) != 0, TRUE)

#Leaf
bac.clean.ss.17.leaf.SW <- subset_samples(bac.clean.ss.17.leaf, Location == "Suwon")
bac.clean.ss.17.leaf.SW <- phyloseq::filter_taxa(bac.clean.ss.17.leaf.SW, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.leaf.CC1 <- subset_samples(bac.clean.ss.17.leaf, Location == "Chuncheon1")
bac.clean.ss.17.leaf.CC1 <- phyloseq::filter_taxa(bac.clean.ss.17.leaf.CC1, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.leaf.CC2 <- subset_samples(bac.clean.ss.17.leaf, Location == "Chuncheon2")
bac.clean.ss.17.leaf.CC2 <- phyloseq::filter_taxa(bac.clean.ss.17.leaf.CC2, function(x) sum(x) != 0, TRUE)

#Stem
bac.clean.ss.17.stem.SW <- subset_samples(bac.clean.ss.17.stem, Location == "Suwon")
bac.clean.ss.17.stem.SW <- phyloseq::filter_taxa(bac.clean.ss.17.stem.SW, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.stem.CC1 <- subset_samples(bac.clean.ss.17.stem, Location == "Chuncheon1")
bac.clean.ss.17.stem.CC1 <- phyloseq::filter_taxa(bac.clean.ss.17.stem.CC1, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.stem.CC2 <- subset_samples(bac.clean.ss.17.stem, Location == "Chuncheon2")
bac.clean.ss.17.stem.CC2 <- phyloseq::filter_taxa(bac.clean.ss.17.stem.CC2, function(x) sum(x) != 0, TRUE)

#Root
bac.clean.ss.17.root.SW <- subset_samples(bac.clean.ss.17.root, Location == "Suwon")
bac.clean.ss.17.root.SW <- phyloseq::filter_taxa(bac.clean.ss.17.root.SW, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.root.CC1 <- subset_samples(bac.clean.ss.17.root, Location == "Chuncheon1")
bac.clean.ss.17.root.CC1 <- phyloseq::filter_taxa(bac.clean.ss.17.root.CC1, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.root.CC2 <- subset_samples(bac.clean.ss.17.root, Location == "Chuncheon2")
bac.clean.ss.17.root.CC2 <- phyloseq::filter_taxa(bac.clean.ss.17.root.CC2, function(x) sum(x) != 0, TRUE)

#Seed
bac.clean.ss.17.seed.SW <- subset_samples(bac.clean.ss.17.seed, Location == "Suwon")
bac.clean.ss.17.seed.SW <- phyloseq::filter_taxa(bac.clean.ss.17.seed.SW, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.seed.CC1 <- subset_samples(bac.clean.ss.17.seed, Location == "Chuncheon1")
bac.clean.ss.17.seed.CC1 <- phyloseq::filter_taxa(bac.clean.ss.17.seed.CC1, function(x) sum(x) != 0, TRUE)
bac.clean.ss.17.seed.CC2 <- subset_samples(bac.clean.ss.17.seed, Location == "Chuncheon2")
bac.clean.ss.17.seed.CC2 <- phyloseq::filter_taxa(bac.clean.ss.17.seed.CC2, function(x) sum(x) != 0, TRUE)


#phy.m <-  merge_samples(bac.clean.ss.leaf, "Age")

order.sample.17 <- c("50_1", "50_2", "50_3", "80_1", "80_2", "80_3","80_4","80_5", "80_6", "80_7", "80_8", "80_9","120_1", "120_2", "120_3","120_4","120_5", "120_6", "120_7", "120_8", "120_9","140_1", "140_2", "140_3","140_4","140_5", "140_6", "140_7", "140_8", "140_9")

order.sample.17 <- c("50days", "80days", "120days", "140days")
## let's get phylum sir
plot_RA_bar<-function(phyloseq_form,order.sam,keyword){df.phylum <- phyloseq_form %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

head(df.phylum)
df.phylum$Age <- factor(df.phylum$Age, levels = order.sam)

library(forcats) 
df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))
unique(df.phylum$Phylum2)

levels(df.phylum$Phylum2)
levels(df.phylum$Phylum2) = c(levels(df.phylum$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 5,]$Phylum2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.proteobacteria <- c("Deltaproteobacteria","Gammaproteobacteria","Alphaproteobacteria")
vec.reorder <- append(vec.uniden.Low, vec.order)
vec.reorder <- append(vec.reorder, vec.proteobacteria)

df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Age, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Gammaproteobacteria" = "darkolivegreen3", "Alphaproteobacteria"= "darkolivegreen", 
                               "Actinobacteria"="indianred2","Deltaproteobacteria" ="darkolivegreen1", 
                               "Bacteroidetes"="steelblue1", "Firmicutes" ="tan1", "Acidobacteria"="lightsalmon4", 
                               "Chloroflexi"="gold1", "Verrucomicrobia"="orchid3", "Nitrospirae"="palevioletred2",
                              "Planctomycetes" = "seagreen3", "Cyanobacteria" = "chartreuse4", "Epsilonbacteraeota" = "darkslategray4", 
                             "Spirochaetes" = "bisque4", "Deinococcus-Thermus" = "lightpink3", "Patescibacteria" = "lightblue1", 
                              "Lentisphaerae" = "lightgoldenrod3", "Dependentiae" = "chocolate3",
                               "Fibrobacteres" = "sienna2", "Armatimonadetes" = "coral", "Nitrospinae" = "darkkhaki", "Chlamydiae" = "darkslateblue",
                               "Kiritimatiellaeota" = "hotpink3", "Rokubacteria" = "lightsteelblue3", "Elusimicrobia" = "turquoise3",
                               "Fusobacteria" = "mistyrose3", "BRC1" = "plum2", "Latescibacteria" = "darkorchid2", "FBP" = "burlywood2",
                               "WPS-2" = "cornflowerblue", "Tenericutes" = "paleturquoise3", 
                              "Gemmatimonadetes"= "peachpuff3", "Low abundance" = "light grey", "unidentified" = "black")) +

  xlab('')+ theme(aspect.ratio = 1)+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

write.csv(df.phylum.rel, paste0('Source data of RA plot_',keyword,".csv"))
return(df.phylum.rel.p1)}

#Suwon
plot_RA_bar(bac.clean.ss.17.soil.SW,order.sample.17,"SW_soil")
plot_RA_bar(bac.clean.ss.17.leaf.SW,order.sample.17,"SW_leaf")
plot_RA_bar(bac.clean.ss.17.stem.SW,order.sample.17,"SW_stem")
plot_RA_bar(bac.clean.ss.17.root.SW,order.sample.17,"SW_root")
plot_RA_bar(bac.clean.ss.17.seed.SW,order.sample.17,"SW_seed")

#Chuncheon 1
plot_RA_bar(bac.clean.ss.17.soil.CC1,order.sample.17,"CC1_soil")
plot_RA_bar(bac.clean.ss.17.leaf.CC1,order.sample.17,"CC1_leaf")
plot_RA_bar(bac.clean.ss.17.stem.CC1,order.sample.17,"CC1_stem")
plot_RA_bar(bac.clean.ss.17.root.CC1,order.sample.17,"CC1_root")
plot_RA_bar(bac.clean.ss.17.seed.CC1,order.sample.17,"CC1_seed")

#Chuncheon 2
plot_RA_bar(bac.clean.ss.17.soil.CC2,order.sample.17,"CC2_soil")
plot_RA_bar(bac.clean.ss.17.leaf.CC2,order.sample.17,"CC2_leaf")
plot_RA_bar(bac.clean.ss.17.stem.CC2,order.sample.17,"CC2_stem")
plot_RA_bar(bac.clean.ss.17.root.CC2,order.sample.17,"CC2_root")
plot_RA_bar(bac.clean.ss.17.seed.CC2,order.sample.17,"CC2_seed")

dev.off()


#2018
#Sorted by Compartment
map.18 <- subset(map, Year == "year_2018")
rownames(map.18) <- map.18$SampleID
sample_data(bac.clean.ss.18) <- sample_data(map.18)
bac.clean.ss.18.f.leaf <- subset_samples(bac.clean.ss.18, Compartment == "Leaf")
bac.clean.ss.18.f.leaf <- phyloseq::filter_taxa(bac.clean.ss.18.f.leaf, function(x) sum(x) != 0, TRUE)
sample_data(bac.clean.ss.18.f.leaf)
bac.clean.ss.18.f.stem <- subset_samples(bac.clean.ss.18, Compartment == "Stem")
bac.clean.ss.18.f.stem <- phyloseq::filter_taxa(bac.clean.ss.18.f.stem, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.root <- subset_samples(bac.clean.ss.18, Compartment == "Root")
bac.clean.ss.18.f.root <- phyloseq::filter_taxa(bac.clean.ss.18.f.root, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.seed <- subset_samples(bac.clean.ss.18, Compartment == "Seed")
bac.clean.ss.18.f.seed <- phyloseq::filter_taxa(bac.clean.ss.18.f.seed, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.seed <- phyloseq::filter_taxa(bac.clean.ss.18.f.seed, function(x) sum(x) != 0, TRUE)
(filt.sample <- sample_sums(bac.clean.ss.18.f.seed ) > 0)
sum(sample_sums(bac.clean.ss.18.f.seed ) <= 0)  ## 1 sample discarded
bac.clean.ss.18.f.seed  <- prune_samples(filt.sample,bac.clean.ss.18.f.seed )
bac.clean.ss.18.f.seed   ## 979 samples <- 984 samples

bac.clean.ss.18.f.BS <- subset_samples(bac.clean.ss.18, Compartment == "Bulk_soil")
bac.clean.ss.18.f.BS <- phyloseq::filter_taxa(bac.clean.ss.18.f.BS, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.RS <- subset_samples(bac.clean.ss.18, Compartment == "Rhizosphere")
bac.clean.ss.18.f.RS <- phyloseq::filter_taxa(bac.clean.ss.18.f.RS, function(x) sum(x) != 0, TRUE)

##Ordering for 2018 samples
order.sample.18 <- c("0_1", "0_2", "0_3", "0_4","0_5","0_6","0_7","0_8","0_9",
                  "0_10","0_11","0_12","0_13","0_14","0_15","0_16","0_17","0_18",
                  "0_19","0_20","0_21","0_22","0_23","0_24","0_25","0_26","0_27",
                  "48_1", "48_2", "48_3", "48_4","48_5","48_6","48_7","48_8","48_9",
                  "62_1", "62_2", "62_3", "62_4","62_5","62_6","62_7","62_8","62_9",
                  "76_1", "76_2", "76_3", "76_4","76_5","76_6","76_7","76_8","76_9",
                  "76_10","76_11","76_12","76_13","76_14","76_15","76_16","76_17","76_18",
                  "76_19","76_20","76_21","76_22","76_23","76_24","76_25","76_26","76_27",
                  "90_1", "90_2", "90_3", "90_4","90_5","90_6","90_7","90_8","90_9",
                  "90_10","90_11","90_12","90_13","90_14","90_15","90_16","90_17","90_18",
                  "90_19","90_20","90_21","90_22","90_23","90_24","90_25","90_26","90_27",
                  "106_1", "106_2", "106_3", "106_4","106_5","106_6","106_7","106_8","106_9",
                  "106_10","106_11","106_12","106_13","106_14","106_15","106_16","106_17","106_18",
                  "106_19","106_20","106_21","106_22","106_23","106_24","106_25","106_26","106_27",
                  "120_1", "120_2", "120_3", "120_4","120_5","120_6","120_7","120_8","120_9",
                  "120_10","120_11","120_12","120_13","120_14","120_15","120_16","120_17","120_18",
                  "120_19","120_20","120_21","120_22","120_23","120_24","120_25","120_26","120_27",
                  "141_1", "141_2", "141_3", "141_4","141_5","141_6","141_7","141_8","141_9",
                  "141_10","141_11","141_12","141_13","141_14","141_15","141_16","141_17","141_18",
                  "141_19","141_20","141_21","141_22","141_23","141_24","141_25","141_26","141_27")


### RA plot for each section of leaves and stems
bac.clean.ss.18.f.leaf <- subset_samples(bac.clean.ss.18.f, Compartment == "Leaf")
bac.clean.ss.18.f.leaf <- phyloseq::filter_taxa(bac.clean.ss.18.f.leaf, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.stem <- subset_samples(bac.clean.ss.18.f, Compartment == "Stem")
bac.clean.ss.18.f.stem <- phyloseq::filter_taxa(bac.clean.ss.18.f.stem, function(x) sum(x) != 0, TRUE)

## sectioning leaf compartment
bac.clean.ss.18.f.L1 <- subset_samples(bac.clean.ss.18.f.leaf, Microhabitat == "L1")
bac.clean.ss.18.f.L1 <- phyloseq::filter_taxa(bac.clean.ss.18.f.L1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.L2 <- subset_samples(bac.clean.ss.18.f.leaf, Microhabitat == "L2")
bac.clean.ss.18.f.L2 <- phyloseq::filter_taxa(bac.clean.ss.18.f.L2, function(x) sum(x) != 0, TRUE)
(filt.sample <- sample_sums(bac.clean.ss.18.f.L2 ) > 0)
sum(sample_sums(bac.clean.ss.18.f.L2 ) <= 0)  ## 1 sample discarded
bac.clean.ss.18.f.L2  <- prune_samples(filt.sample,bac.clean.ss.18.f.L2 )
bac.clean.ss.18.f.L2   ## 979 samples <- 984 samples

bac.clean.ss.18.f.L3 <- subset_samples(bac.clean.ss.18.f.leaf, Microhabitat == "L3")
bac.clean.ss.18.f.L3 <- phyloseq::filter_taxa(bac.clean.ss.18.f.L3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.FL <- subset_samples(bac.clean.ss.18.f.leaf, Microhabitat == "FL")
bac.clean.ss.18.f.FL <- phyloseq::filter_taxa(bac.clean.ss.18.f.FL, function(x) sum(x) != 0, TRUE)

## sectioning stem compartment
bac.clean.ss.18.f.S1 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S1")
bac.clean.ss.18.f.S1 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S1 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S1, function(x) sum(x) != 0, TRUE)
(filt.sample <- sample_sums(bac.clean.ss.18.f.S1 ) > 0)
sum(sample_sums(bac.clean.ss.18.f.S1 ) <= 0)  ## 1 sample discarded
bac.clean.ss.18.f.S1  <- prune_samples(filt.sample,bac.clean.ss.18.f.S1 )
bac.clean.ss.18.f.S1   ## 979 samples <- 984 samples

bac.clean.ss.18.f.S2 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S2")
bac.clean.ss.18.f.S2 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S3 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S3")
bac.clean.ss.18.f.S3 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S4 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S4")
bac.clean.ss.18.f.S4 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S4, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S5 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S5")
bac.clean.ss.18.f.S5 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S5, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S6 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S6")
bac.clean.ss.18.f.S6 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S6, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S7 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S7")
bac.clean.ss.18.f.S7 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S7, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S8 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S8")
bac.clean.ss.18.f.S8 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S8, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.f.S9 <- subset_samples(bac.clean.ss.18.f.stem, Microhabitat == "S9")
bac.clean.ss.18.f.S9 <- phyloseq::filter_taxa(bac.clean.ss.18.f.S9, function(x) sum(x) != 0, TRUE)



plot_RA_bar(bac.clean.ss.18.f.BS,order.sample.18,"BS")
plot_RA_bar(bac.clean.ss.18.f.RS,order.sample.18,"RS")
plot_RA_bar(bac.clean.ss.18.f.seed,order.sample.18,"G")
plot_RA_bar(bac.clean.ss.18.f.root,order.sample.18,"R")

plot_RA_bar(bac.clean.ss.18.f.L1,order.sample.18,"L1")
plot_RA_bar(bac.clean.ss.18.f.L2,order.sample.18,"L2")
plot_RA_bar(bac.clean.ss.18.f.L3,order.sample.18,"L3")
plot_RA_bar(bac.clean.ss.18.f.FL,order.sample.18,"FL")

plot_RA_bar(bac.clean.ss.18.f.S1,order.sample.18,"S1")
plot_RA_bar(bac.clean.ss.18.f.S2,order.sample.18,"S2")
plot_RA_bar(bac.clean.ss.18.f.S3,order.sample.18,"S3")
plot_RA_bar(bac.clean.ss.18.f.S4,order.sample.18,"S4")
plot_RA_bar(bac.clean.ss.18.f.S5,order.sample.18,"S5")
plot_RA_bar(bac.clean.ss.18.f.S6,order.sample.18,"S6")
plot_RA_bar(bac.clean.ss.18.f.S7,order.sample.18,"S7")
plot_RA_bar(bac.clean.ss.18.f.S8,order.sample.18,"S8")
plot_RA_bar(bac.clean.ss.18.f.S9,order.sample.18,"S9")

dev.off()


## Merge replicates
plot_RA_bar.18.age<-function(phyloseq_form){df.phylum <- phyloseq_form %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.phylum$Age <- factor(df.phylum$Age, levels = order.age)



library(forcats) 
df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))
unique(df.phylum$Phylum2)

levels(df.phylum$Phylum2)
levels(df.phylum$Phylum2) = c(levels(df.phylum$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 1,]$Phylum2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)


df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Age, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

return(df.phylum.rel.p1)}


plot_RA_bar.18.age(bac.clean.ss.18.f.L1)
plot_RA_bar.18.age(bac.clean.ss.18.f.L2)
plot_RA_bar.18.age(bac.clean.ss.18.f.L3)
plot_RA_bar.18.age(bac.clean.ss.18.f.FL)

plot_RA_bar.18.age(bac.clean.ss.18.f.S1)
plot_RA_bar.18.age(bac.clean.ss.18.f.S2)
plot_RA_bar.18.age(bac.clean.ss.18.f.S3)
plot_RA_bar.18.age(bac.clean.ss.18.f.S4)
plot_RA_bar.18.age(bac.clean.ss.18.f.S5)
plot_RA_bar.18.age(bac.clean.ss.18.f.S6)
plot_RA_bar.18.age(bac.clean.ss.18.f.S7)
plot_RA_bar.18.age(bac.clean.ss.18.f.S8)
plot_RA_bar.18.age(bac.clean.ss.18.f.S9)

plot_RA_bar.18.age(bac.clean.ss.18.f.BS)
plot_RA_bar.18.age(bac.clean.ss.18.f.RS)
plot_RA_bar.18.age(bac.clean.ss.18.f.root)
plot_RA_bar.18.age(bac.clean.ss.18.f.seed)



### Seed genus
plot_RA_bar.genus<-function(phyloseq_form){df.genus <- phyloseq_form %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

order.age <- c('0days', '48days', '62days', '76days', '90days', '106days', '120days', '141days')

df.genus$Age <- factor(df.genus$Age, levels = order.age)



library(forcats) 
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
unique(df.genus$Genus)

levels(df.genus$Genus)
levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 5,]$Genus <- 'Low abundance'
ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)


df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=Age, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_otu2) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 10,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.genus.rel.p1

return(df.genus.rel.p1)}

plot_RA_bar.genus(bac.clean.ss.18.f.seed)
plot_RA_bar.genus(bac.clean.ss.18.f.L1)

dev.off()


### Fungal communities
fun.clean.ss.17 <- subset_samples(fun.clean.ss2, Year == "year_2017")
fun.clean.ss.18 <- subset_samples(fun.clean.ss2, Year == "year_2018")

fun.clean.ss.17.leaf <- subset_samples(fun.clean.ss.17, Compartment == "Leaf")
fun.clean.ss.17.leaf <- phyloseq::filter_taxa(fun.clean.ss.17.leaf, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.stem <- subset_samples(fun.clean.ss.17, Compartment == "Stem")
fun.clean.ss.17.stem <- phyloseq::filter_taxa(fun.clean.ss.17.stem, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.root <- subset_samples(fun.clean.ss.17, Compartment == "Root")
fun.clean.ss.17.root <- phyloseq::filter_taxa(fun.clean.ss.17.root, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.seed <- subset_samples(fun.clean.ss.17, Compartment == "Seed")
fun.clean.ss.17.seed <- phyloseq::filter_taxa(fun.clean.ss.17.seed, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.soil <- subset_samples(fun.clean.ss.17, Compartment == "Soil")
fun.clean.ss.17.soil <- phyloseq::filter_taxa(fun.clean.ss.17.soil, function(x) sum(x) != 0, TRUE)

#Soil
fun.clean.ss.17.soil.SW <- subset_samples(fun.clean.ss.17.soil, Location == "Suwon")
fun.clean.ss.17.soil.SW <- phyloseq::filter_taxa(fun.clean.ss.17.soil.SW, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.soil.CC1 <- subset_samples(fun.clean.ss.17.soil, Location == "Chuncheon1")
fun.clean.ss.17.soil.CC1 <- phyloseq::filter_taxa(fun.clean.ss.17.soil.CC1, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.soil.CC2 <- subset_samples(fun.clean.ss.17.soil, Location == "Chuncheon2")
fun.clean.ss.17.soil.CC2 <- phyloseq::filter_taxa(fun.clean.ss.17.soil.CC2, function(x) sum(x) != 0, TRUE)

#Leaf
fun.clean.ss.17.leaf.SW <- subset_samples(fun.clean.ss.17.leaf, Location == "Suwon")
fun.clean.ss.17.leaf.SW <- phyloseq::filter_taxa(fun.clean.ss.17.leaf.SW, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.leaf.CC1 <- subset_samples(fun.clean.ss.17.leaf, Location == "Chuncheon1")
fun.clean.ss.17.leaf.CC1 <- phyloseq::filter_taxa(fun.clean.ss.17.leaf.CC1, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.leaf.CC2 <- subset_samples(fun.clean.ss.17.leaf, Location == "Chuncheon2")
fun.clean.ss.17.leaf.CC2 <- phyloseq::filter_taxa(fun.clean.ss.17.leaf.CC2, function(x) sum(x) != 0, TRUE)

#Stem
fun.clean.ss.17.stem.SW <- subset_samples(fun.clean.ss.17.stem, Location == "Suwon")
fun.clean.ss.17.stem.SW <- phyloseq::filter_taxa(fun.clean.ss.17.stem.SW, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.stem.CC1 <- subset_samples(fun.clean.ss.17.stem, Location == "Chuncheon1")
fun.clean.ss.17.stem.CC1 <- phyloseq::filter_taxa(fun.clean.ss.17.stem.CC1, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.stem.CC2 <- subset_samples(fun.clean.ss.17.stem, Location == "Chuncheon2")
fun.clean.ss.17.stem.CC2 <- phyloseq::filter_taxa(fun.clean.ss.17.stem.CC2, function(x) sum(x) != 0, TRUE)

#Root
fun.clean.ss.17.root.SW <- subset_samples(fun.clean.ss.17.root, Location == "Suwon")
fun.clean.ss.17.root.SW <- phyloseq::filter_taxa(fun.clean.ss.17.root.SW, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.root.CC1 <- subset_samples(fun.clean.ss.17.root, Location == "Chuncheon1")
fun.clean.ss.17.root.CC1 <- phyloseq::filter_taxa(fun.clean.ss.17.root.CC1, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.root.CC2 <- subset_samples(fun.clean.ss.17.root, Location == "Chuncheon2")
fun.clean.ss.17.root.CC2 <- phyloseq::filter_taxa(fun.clean.ss.17.root.CC2, function(x) sum(x) != 0, TRUE)

#Seed
fun.clean.ss.17.seed.SW <- subset_samples(fun.clean.ss.17.seed, Location == "Suwon")
fun.clean.ss.17.seed.SW <- phyloseq::filter_taxa(fun.clean.ss.17.seed.SW, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.seed.CC1 <- subset_samples(fun.clean.ss.17.seed, Location == "Chuncheon1")
fun.clean.ss.17.seed.CC1 <- phyloseq::filter_taxa(fun.clean.ss.17.seed.CC1, function(x) sum(x) != 0, TRUE)
fun.clean.ss.17.seed.CC2 <- subset_samples(fun.clean.ss.17.seed, Location == "Chuncheon2")
fun.clean.ss.17.seed.CC2 <- phyloseq::filter_taxa(fun.clean.ss.17.seed.CC2, function(x) sum(x) != 0, TRUE)


#phy.m <-  merge_samples(fun.clean.ss.leaf, "Age")

order.sample.17 <- c("50_1", "50_2", "50_3", "80_1", "80_2", "80_3","80_4","80_5", "80_6", "80_7", "80_8", "80_9","120_1", "120_2", "120_3","120_4","120_5", "120_6", "120_7", "120_8", "120_9","140_1", "140_2", "140_3","140_4","140_5", "140_6", "140_7", "140_8", "140_9")


## let's get phylum sir
plot_RA_bar_fun<-function(phyloseq_form,order.sam,keyword){df.phylum <- phyloseq_form %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

head(df.phylum)
df.phylum$Age <- factor(df.phylum$Age, levels = order.sam)

library(forcats) 
df.phylum %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
unique(df.phylum$Class)

levels(df.phylum$Class)
levels(df.phylum$Class) = c(levels(df.phylum$Class), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum %>%  
  group_by(Age) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 5,]$Class <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Class
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)

df.phylum.rel$Class <- factor(df.phylum.rel$Class, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Age, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Mortierellomycetes" = "#4E734E", "Eurotiomycetes"= "#6DA9DC",
                               "Mucoromycetes"="#E4AF2C","Tremellomycetes"="#BE4146", "Microbotryomycetes" ="#DC9A9E",
                               "Dothideomycetes" = "#5195D1", "Agaricomycetes" = "#CC6C71", "Blastocladiomycetes" = "#87AC88",
                               "Sordariomycetes"= "#1E63AF","Leotiomycetes"= "#11335F", "Cystobasidiomycetes" = "#A871AE", "Pezizomycetes" = "#C0DBF3",
                               "unidentified" ="#000000", "Endogonomycetes" ="#CC9900", "Pucciniomycetes" = "#0099FF", "Ustilaginomycetes"="#ffcc33","Arthoniomycetes" = "#cccc99",
                               "Saccharomycetes" = "#666699", "Low abundance" = "light grey", "Malasseziomycetes" = "#99cc99", "Exobasidiomycetes" = "#663366", 
                               "Rhizophydiomycetes"="#ffccff","Wallemiomycetes"="#333300","Agaricostilbomycetes"="#003333","Orbiliomycetes"="#cc9966","Chytridiomycetes" = "#006600",
                               "Spizellomycetes" = "#66cc99","Moniliellomycetes" = "#0099cc","Monoblepharidomycetes" = "#CC0000", "Dacrymycetes" = "#669999",
                               "Spiculogloeomycetes" = "#99cc00")) +

  xlab('')+
  ylab("Relative abundance\n") + theme(aspect.ratio = 1)+
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,200))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

write.csv(df.phylum.rel, paste0("Fungi_",keyword,".csv"))

return(df.phylum.rel.p1)}

#Suwon
plot_RA_bar_fun(fun.clean.ss.17.soil.SW,order.sample.17,"SW_soil")
plot_RA_bar_fun(fun.clean.ss.17.leaf.SW,order.sample.17,"SW_leaf")
plot_RA_bar_fun(fun.clean.ss.17.stem.SW,order.sample.17,"SW_stem")
plot_RA_bar_fun(fun.clean.ss.17.root.SW,order.sample.17,"SW_root")
plot_RA_bar_fun(fun.clean.ss.17.seed.SW,order.sample.17,"SW_seed")

#Chuncheon 1
plot_RA_bar_fun(fun.clean.ss.17.soil.CC1,order.sample.17,"CC1_soil")
plot_RA_bar_fun(fun.clean.ss.17.leaf.CC1,order.sample.17,"CC1_leaf")
plot_RA_bar_fun(fun.clean.ss.17.stem.CC1,order.sample.17,"CC1_stem")
plot_RA_bar_fun(fun.clean.ss.17.root.CC1,order.sample.17,"CC1_root")
plot_RA_bar_fun(fun.clean.ss.17.seed.CC1,order.sample.17,"CC1_seed")

#Chuncheon 2
plot_RA_bar_fun(fun.clean.ss.17.soil.CC2,order.sample.17,"CC2_soil")
plot_RA_bar_fun(fun.clean.ss.17.leaf.CC2,order.sample.17,"CC2_leaf")
plot_RA_bar_fun(fun.clean.ss.17.stem.CC2,order.sample.17,"CC2_stem")
plot_RA_bar_fun(fun.clean.ss.17.root.CC2,order.sample.17,"CC2_root")
plot_RA_bar_fun(fun.clean.ss.17.seed.CC2,order.sample.17,"CC2_seed")

dev.off()


#2018
#Sorted by Compartment
fun.clean.ss.f

fun.clean.ss.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)


fun.clean.ss.18.f.leaf <- subset_samples(fun.clean.ss.18, Compartment == "Leaf")
fun.clean.ss.18.f.leaf <- phyloseq::filter_taxa(fun.clean.ss.18.f.leaf, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.stem <- subset_samples(fun.clean.ss.18, Compartment == "Stem")
fun.clean.ss.18.f.stem <- phyloseq::filter_taxa(fun.clean.ss.18.f.stem, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.root <- subset_samples(fun.clean.ss.18, Compartment == "Root")
fun.clean.ss.18.f.root <- phyloseq::filter_taxa(fun.clean.ss.18.f.root, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.seed <- subset_samples(fun.clean.ss.18, Compartment == "Grain")
fun.clean.ss.18.f.seed <- phyloseq::filter_taxa(fun.clean.ss.18.f.seed, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.BS <- subset_samples(fun.clean.ss.18, Compartment == "Bulk_soil")
fun.clean.ss.18.f.BS <- phyloseq::filter_taxa(fun.clean.ss.18.f.BS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.RS <- subset_samples(fun.clean.ss.18, Compartment == "Rhizosphere")
fun.clean.ss.18.f.RS <- phyloseq::filter_taxa(fun.clean.ss.18.f.RS, function(x) sum(x) != 0, TRUE)

##Ordering for 2018 samples
order.sample.18 <- c("0_1", "0_2", "0_3", "0_4","0_5","0_6","0_7","0_8","0_9",
                     "0_10","0_11","0_12","0_13","0_14","0_15","0_16","0_17","0_18",
                     "0_19","0_20","0_21","0_22","0_23","0_24","0_25","0_26","0_27",
                     "48_1", "48_2", "48_3", "48_4","48_5","48_6","48_7","48_8","48_9",
                     "62_1", "62_2", "62_3", "62_4","62_5","62_6","62_7","62_8","62_9",
                     "76_1", "76_2", "76_3", "76_4","76_5","76_6","76_7","76_8","76_9",
                     "76_10","76_11","76_12","76_13","76_14","76_15","76_16","76_17","76_18",
                     "76_19","76_20","76_21","76_22","76_23","76_24","76_25","76_26","76_27",
                     "90_1", "90_2", "90_3", "90_4","90_5","90_6","90_7","90_8","90_9",
                     "90_10","90_11","90_12","90_13","90_14","90_15","90_16","90_17","90_18",
                     "90_19","90_20","90_21","90_22","90_23","90_24","90_25","90_26","90_27",
                     "106_1", "106_2", "106_3", "106_4","106_5","106_6","106_7","106_8","106_9",
                     "106_10","106_11","106_12","106_13","106_14","106_15","106_16","106_17","106_18",
                     "106_19","106_20","106_21","106_22","106_23","106_24","106_25","106_26","106_27",
                     "120_1", "120_2", "120_3", "120_4","120_5","120_6","120_7","120_8","120_9",
                     "120_10","120_11","120_12","120_13","120_14","120_15","120_16","120_17","120_18",
                     "120_19","120_20","120_21","120_22","120_23","120_24","120_25","120_26","120_27",
                     "141_1", "141_2", "141_3", "141_4","141_5","141_6","141_7","141_8","141_9",
                     "141_10","141_11","141_12","141_13","141_14","141_15","141_16","141_17","141_18",
                     "141_19","141_20","141_21","141_22","141_23","141_24","141_25","141_26","141_27")


## sectioning leaf compartment
fun.clean.ss.18.f.L1 <- subset_samples(fun.clean.ss.18.f.leaf, Microhabitat == "L1")
fun.clean.ss.18.f.L1 <- phyloseq::filter_taxa(fun.clean.ss.18.f.L1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.L2 <- subset_samples(fun.clean.ss.18.f.leaf, Microhabitat == "L2")
fun.clean.ss.18.f.L2 <- phyloseq::filter_taxa(fun.clean.ss.18.f.L2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.L3 <- subset_samples(fun.clean.ss.18.f.leaf, Microhabitat == "L3")
fun.clean.ss.18.f.L3 <- phyloseq::filter_taxa(fun.clean.ss.18.f.L3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.FL <- subset_samples(fun.clean.ss.18.f.leaf, Microhabitat == "FL")
fun.clean.ss.18.f.FL <- phyloseq::filter_taxa(fun.clean.ss.18.f.FL, function(x) sum(x) != 0, TRUE)

## sectioning stem compartment
fun.clean.ss.18.f.S1 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S1")
fun.clean.ss.18.f.S1 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S2 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S2")
fun.clean.ss.18.f.S2 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S3 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S3")
fun.clean.ss.18.f.S3 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S4 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S4")
fun.clean.ss.18.f.S4 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S4, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S5 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S5")
fun.clean.ss.18.f.S5 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S5, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S6 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S6")
fun.clean.ss.18.f.S6 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S6, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S7 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S7")
fun.clean.ss.18.f.S7 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S7, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S8 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S8")
fun.clean.ss.18.f.S8 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S8, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.f.S9 <- subset_samples(fun.clean.ss.18.f.stem, Microhabitat == "S9")
fun.clean.ss.18.f.S9 <- phyloseq::filter_taxa(fun.clean.ss.18.f.S9, function(x) sum(x) != 0, TRUE)



plot_RA_bar_fun(fun.clean.ss.18.f.BS,order.sample.18,"BS")
plot_RA_bar_fun(fun.clean.ss.18.f.RS,order.sample.18,"RS")
plot_RA_bar_fun(fun.clean.ss.18.f.seed,order.sample.18,"G")
plot_RA_bar_fun(fun.clean.ss.18.f.root,order.sample.18,"R")

plot_RA_bar_fun(fun.clean.ss.18.f.L1,order.sample.18,"L1")
plot_RA_bar_fun(fun.clean.ss.18.f.L2,order.sample.18,"L2")
plot_RA_bar_fun(fun.clean.ss.18.f.L3,order.sample.18,"L3")
plot_RA_bar_fun(fun.clean.ss.18.f.FL,order.sample.18,"FL")

plot_RA_bar_fun(fun.clean.ss.18.f.S1,order.sample.18,"S1")
plot_RA_bar_fun(fun.clean.ss.18.f.S2,order.sample.18,"S2")
plot_RA_bar_fun(fun.clean.ss.18.f.S3,order.sample.18,"S3")
plot_RA_bar_fun(fun.clean.ss.18.f.S4,order.sample.18,"S4")
plot_RA_bar_fun(fun.clean.ss.18.f.S5,order.sample.18,"S5")
plot_RA_bar_fun(fun.clean.ss.18.f.S6,order.sample.18,"S6")
plot_RA_bar_fun(fun.clean.ss.18.f.S7,order.sample.18,"S7")
plot_RA_bar_fun(fun.clean.ss.18.f.S8,order.sample.18,"S8")
plot_RA_bar_fun(fun.clean.ss.18.f.S9,order.sample.18,"S9")

dev.off()
