##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
order.compartment <- c("Leaf", "Stem","Root", "Seed", "Soil", "Rhizosphere", "Bulk_soil")
order.microhabitat.leaf <- c('L1', 'L2', "L3", 'FL')
order.microhabitat.stem <- c('S1', 'S2', "S3", 'S4','S5','S6','S7','S8','S9')

map$Compartment <- factor(map$Compartment, levels = order.compartment)



##Unconstrained PCoA
color_compartment <- c("Leaf" = "#336633", "Stem"="#99CC33", "Root" = "#CC9999", "Seed" = "#6699CC", "Soil" = "#663300", "Bulk_soil"="#663300", "Rhizosphere" = "#CC9900")
shape_age <- c("0days" = 18,"50days" = 0, "80days"= 1, "120days"= 2, "140days"=8, "48days" = 15, "62days"= 16, "76days"= 17, "90days"= 5, "106days"= 7, "141days"= 13)

## make Days as numeric
b.meta.17<-sample_data(bac.clean.log.17)
b.meta.17$Days <- as.numeric(as.character(b.meta.17$Days))
b.meta.17$Compartment <- factor(b.meta.17$Compartment, levels = order.compartment)

sample_data(bac.clean.log.17) <- sample_data(b.meta.17)

b.meta.18<-sample_data(bac.clean.log.18)
b.meta.18$Days <- as.numeric(as.character(b.meta.18$Days))
b.meta.18$Compartment <- factor(b.meta.18$Compartment, levels = order.compartment)

sample_data(bac.clean.log.18) <- sample_data(b.meta.18)

bray1.bac <-  ordinate(bac.clean.log.17, 'PCoA', 'bray')
bray2.bac <-  ordinate(bac.clean.log.18, 'PCoA', 'bray')
bray3.bac <-  ordinate(bac.clean.log, 'NMDS', 'bray')

write.csv(bray1.bac$vectors, "Bacteria_PCoA_2017.csv")
write.csv(bray2.bac$vectors, "Bacteria_PCoA_2018.csv")

#Age and compartment
##all samples - NMDS
bray3.bac$stress #0.1475148
plot_ordination(bac.clean.log, bray3.bac, type = "samples", color='Year', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(bac.clean.log, bray3.bac, type = "samples", color='Compartment', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(bac.clean.log, bray3.bac, type = "samples", color='Location', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## 2017 bac
#Compartment
plot_ordination(bac.clean.log.17, bray1.bac, type = "samples", color='Compartment', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(bac.clean.log.17, bray1.bac, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## location

plot_ordination(bac.clean.log.17, bray1.bac, type = "samples", color='Location', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


##2018 bac
#Compartment
plot_ordination(bac.clean.log.18, bray2.bac, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Microhabitat
plot_ordination(bac.clean.log.18, bray2.bac, type = "samples", color='Microhabitat',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  #scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



#Age
plot_ordination(bac.clean.log.18, bray2.bac, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


### Each compartment

microhabitat_shape.L <- c("FL" = 7, "L1"= 0, "L2"= 1, "L3" = 2) 
microhabitat_shape.S <- c('S1' = 0, 'S2' = 1, "S3"= 2, 'S4'= 5,'S5'= 7,'S6'= 8,'S7'= 9,'S8'= 10,'S9'= 13)
## make Days as numeric
b.meta.17<-sample_data(bac.clean.log.17)
b.meta.17$Days <- as.numeric(as.character(b.meta.17$Days))

sample_data(bac.clean.log.17) <- sample_data(b.meta.17)

b.meta.18.L<-sample_data(bac.clean.log.18.f.leaf)
b.meta.18.L$Days <- as.numeric(as.character(b.meta.18.L$Days))

sample_data(bac.clean.log.18.f.leaf) <- sample_data(b.meta.18.L)

b.meta.18.seed<-sample_data(bac.clean.log.18.f.seed)
b.meta.18.seed$Days <- as.numeric(as.character(b.meta.18.seed$Days))

sample_data(bac.clean.log.18.f.seed) <- sample_data(b.meta.18.seed)

b.meta.18.stem<-sample_data(bac.clean.log.18.f.stem)
b.meta.18.stem$Days <- as.numeric(as.character(b.meta.18.stem$Days))

sample_data(bac.clean.log.18.f.stem) <- sample_data(b.meta.18.stem)


## Above and below-ground parts
bac.clean.log.18.above <- subset_samples(bac.clean.log.18, Compartment %in% c("Leaf", "Stem", "Seed"))
bac.clean.log.18.above <- phyloseq::filter_taxa(bac.clean.log.18.above, function(x) sum(x) != 0, TRUE)

bac.clean.log.18.below <- subset_samples(bac.clean.log.18, Compartment %in% c("Root", "Rhizosphere", "Bulk_soil"))
bac.clean.log.18.below <- phyloseq::filter_taxa(bac.clean.log.18.below, function(x) sum(x) != 0, TRUE)

b.meta.18.above<-sample_data(bac.clean.log.18.above)
b.meta.18.above$Days <- as.numeric(as.character(b.meta.18.above$Days))

sample_data(bac.clean.log.18.above) <- sample_data(b.meta.18.above)


b.meta.18.below<-sample_data(bac.clean.log.18.below)
b.meta.18.below$Days <- as.numeric(as.character(b.meta.18.below$Days))

sample_data(bac.clean.log.18.below) <- sample_data(b.meta.18.below)


bray1.bac <-  ordinate(bac.clean.log.17.leaf, 'PCoA', 'bray')
bray2.bac <-  ordinate(bac.clean.log.18.f.leaf, 'PCoA', 'bray')
bray3.bac <-  ordinate(bac.clean.log.18.f.seed, 'PCoA', 'bray')
bray4.bac <-  ordinate(bac.clean.log.18.f.stem, 'PCoA', 'bray')
bray5.bac <-  ordinate(bac.clean.log.18.above, 'PCoA', 'bray')
bray6.bac <-  ordinate(bac.clean.log.18.below, 'PCoA', 'bray')

##2018 bac
plot_ordination(bac.clean.log.18.f.leaf, bray2.bac, type = "samples", color='Days',shape = "Microhabitat", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_shape_manual(values=microhabitat_shape.L)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##2018 bac
plot_ordination(bac.clean.log.18.f.seed, bray3.bac, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(bac.clean.log.18.f.stem, bray4.bac, type = "samples", color='Days',shape = "Microhabitat", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_shape_manual(values=microhabitat_shape.S)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##aboveground part
#Compartment
plot_ordination(bac.clean.log.18.above, bray5.bac, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(bac.clean.log.18.above, bray5.bac, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##belowground part
#Compartment
plot_ordination(bac.clean.log.18.below, bray6.bac, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(bac.clean.log.18.below, bray6.bac, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



##### Fungal community#####
## make Days as numeric
f.meta.17<-sample_data(fun.clean.log.17)
f.meta.17$Days <- as.numeric(as.character(f.meta.17$Days))
f.meta.17$Compartment <- factor(f.meta.17$Compartment, levels = order.compartment)

sample_data(fun.clean.log.17) <- sample_data(f.meta.17)

f.meta.18<-sample_data(fun.clean.log.18)
f.meta.18$Days <- as.numeric(as.character(f.meta.18$Days))
f.meta.18$Compartment <- as.character(f.meta.18$Compartment)
f.meta.18$Compartment[which(f.meta.18$Compartment=="Grain")] <- "Seed"
f.meta.18$Compartment <- factor(f.meta.18$Compartment, levels = order.compartment)

sample_data(fun.clean.log.18) <- sample_data(f.meta.18)

f.meta<-sample_data(fun.clean.log)
f.meta$Days <- as.numeric(as.character(f.meta$Days))
f.meta$Compartment <- as.character(f.meta$Compartment)
f.meta$Compartment[which(f.meta$Compartment=="Grain")] <- "Seed"
f.meta$Compartment <- factor(f.meta$Compartment, levels = order.compartment)

sample_data(fun.clean.log) <- sample_data(f.meta)

bray1.fun <-  ordinate(fun.clean.log.17, 'PCoA', 'bray')
bray2.fun <-  ordinate(fun.clean.log.18, 'PCoA', 'bray')
bray3.fun <-  ordinate(fun.clean.log, 'NMDS', 'bray')

write.csv(bray1.fun$vectors, "Fungi_PCoA_2017.csv")
write.csv(bray2.fun$vectors, "Fungi_PCoA_2018.csv")


#Age and compartment
##all samples - NMDS
bray3.fun$stress #0.1670793
plot_ordination(fun.clean.log, bray3.fun, type = "samples", color='Year', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(fun.clean.log, bray3.fun, type = "samples", color='Compartment', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(fun.clean.log, bray3.fun, type = "samples", color='Location', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


## 2017 fun
#Compartment
plot_ordination(fun.clean.log.17, bray1.fun, type = "samples", color='Compartment', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(fun.clean.log.17, bray1.fun, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

## location

plot_ordination(fun.clean.log.17, bray1.fun, type = "samples", color='Location', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


##2018 fun
#Compartment
plot_ordination(fun.clean.log.18, bray2.fun, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(fun.clean.log.18, bray2.fun, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


### Each compartment

microhabitat_shape.L <- c("FL" = 7, "L1"= 0, "L2"= 1, "L3" = 2) 
microhabitat_shape.S <- c('S1' = 0, 'S2' = 1, "S3"= 2, 'S4'= 5,'S5'= 7,'S6'= 8,'S7'= 9,'S8'= 10,'S9'= 13)
## make Days as numeric
f.meta.17<-sample_data(fun.clean.log.17)
f.meta.17$Days <- as.numeric(as.character(f.meta.17$Days))

sample_data(fun.clean.log.17) <- sample_data(f.meta.17)

f.meta.18.L<-sample_data(fun.clean.log.18.f.leaf)
f.meta.18.L$Days <- as.numeric(as.character(f.meta.18.L$Days))

sample_data(fun.clean.log.18.f.leaf) <- sample_data(f.meta.18.L)

f.meta.18.seed<-sample_data(fun.clean.log.18.f.seed)
f.meta.18.seed$Days <- as.numeric(as.character(f.meta.18.seed$Days))

sample_data(fun.clean.log.18.f.seed) <- sample_data(f.meta.18.seed)

f.meta.18.stem<-sample_data(fun.clean.log.18.f.stem)
f.meta.18.stem$Days <- as.numeric(as.character(f.meta.18.stem$Days))

sample_data(fun.clean.log.18.f.stem) <- sample_data(f.meta.18.stem)


## Above and below-ground parts
fun.clean.log.18.above <- subset_samples(fun.clean.log.18, Compartment %in% c("Leaf", "Stem", "Seed"))
fun.clean.log.18.above <- phyloseq::filter_taxa(fun.clean.log.18.above, function(x) sum(x) != 0, TRUE)

fun.clean.log.18.below <- subset_samples(fun.clean.log.18, Compartment %in% c("Root", "Rhizosphere", "Bulk_soil"))
fun.clean.log.18.below <- phyloseq::filter_taxa(fun.clean.log.18.below, function(x) sum(x) != 0, TRUE)

f.meta.18.above<-sample_data(fun.clean.log.18.above)
f.meta.18.above$Days <- as.numeric(as.character(f.meta.18.above$Days))

sample_data(fun.clean.log.18.above) <- sample_data(f.meta.18.above)


f.meta.18.below<-sample_data(fun.clean.log.18.below)
f.meta.18.below$Days <- as.numeric(as.character(f.meta.18.below$Days))

sample_data(fun.clean.log.18.below) <- sample_data(f.meta.18.below)


bray1.fun <-  ordinate(fun.clean.log.17.leaf, 'PCoA', 'bray')
bray2.fun <-  ordinate(fun.clean.log.18, 'PCoA', 'bray')
bray3.fun <-  ordinate(fun.clean.log.18.f.seed, 'PCoA', 'bray')
bray4.fun <-  ordinate(fun.clean.log.18.f.stem, 'PCoA', 'bray')
bray5.fun <-  ordinate(fun.clean.log.18.above, 'PCoA', 'bray')
bray6.fun <-  ordinate(fun.clean.log.18.below, 'PCoA', 'bray')

##2018 bac
plot_ordination(fun.clean.log.18.f.leaf, bray2.fun, type = "samples", color='Days',shape = "Microhabitat", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_shape_manual(values=microhabitat_shape.L)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##2018 bac
plot_ordination(fun.clean.log.18.f.seed, bray3.fun, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(fun.clean.log.18.f.stem, bray4.fun, type = "samples", color='Days',shape = "Microhabitat", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_shape_manual(values=microhabitat_shape.S)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##aboveground part
#Compartment
plot_ordination(fun.clean.log.18.above, bray5.fun, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(fun.clean.log.18.above, bray5.fun, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

##belowground part
#Compartment
plot_ordination(fun.clean.log.18.below, bray6.fun, type = "samples", color='Compartment',axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_manual(values=color_compartment)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#Age
plot_ordination(fun.clean.log.18.below, bray6.fun, type = "samples", color='Days', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())




#### PCoA indexed by differentially abundant genus
###Pantoea, Methylobacterium, Xanthomonas, Moesziomyces, Alternaria

bac.clean.rel <- transform(bac.clean.ss.f, transform = "compositional")
fun.clean.rel <- transform(fun.clean.ss2, transform = "compositional")

bac.clean.log.pantoea <- subset_taxa(bac.clean.rel, Genus == "Pantoea")

pantoea.rel<-colSums(otu_table(bac.clean.log.pantoea))

pantoea.rel <- data.frame(pantoea.rel)
names(pantoea.rel)[1] <- "Pantoea"


bac.clean.log.methlyo <- subset_taxa(bac.clean.rel, Genus == "Methylobacterium")

methlyo.rel<-colSums(otu_table(bac.clean.log.methlyo))

methlyo.rel <- data.frame(methlyo.rel)
names(methlyo.rel)[1] <- "Methylobacterium"


bac.clean.log.xantho <- subset_taxa(bac.clean.rel, Genus == "Xanthomonas")

xantho.rel<-colSums(otu_table(bac.clean.log.xantho))

xantho.rel <- data.frame(xantho.rel)
names(xantho.rel)[1] <- "Xanthomonas"


bac.clean.log.pseudo <- subset_taxa(bac.clean.rel, Genus == "Pseudomonas")

pseudo.rel<-colSums(otu_table(bac.clean.log.pseudo))

pseudo.rel <- data.frame(pseudo.rel)
names(pseudo.rel)[1] <- "Pseudomonas"


###Enriched in soil

bac.clean.log.methlyo <- subset_taxa(bac.clean.rel, Genus == "Terrabacter")

methlyo.rel<-colSums(otu_table(bac.clean.log.methlyo))

methlyo.rel <- data.frame(methlyo.rel)
names(methlyo.rel)[1] <- "Terrabacter"


bac.clean.log.xantho <- subset_taxa(bac.clean.rel, Genus == "Methylocystis")

xantho.rel<-colSums(otu_table(bac.clean.log.xantho))

xantho.rel <- data.frame(xantho.rel)
names(xantho.rel)[1] <- "Methylocystis"


bac.clean.log.pseudo <- subset_taxa(bac.clean.rel, Genus == "Anaerolinea")

pseudo.rel<-colSums(otu_table(bac.clean.log.pseudo))

pseudo.rel <- data.frame(pseudo.rel)
names(pseudo.rel)[1] <- "Anaerolinea"

bac.clean.log.pantoea <- subset_taxa(bac.clean.rel, Genus == "Variovorax")

pantoea.rel<-colSums(otu_table(bac.clean.log.pantoea))

pantoea.rel <- data.frame(pantoea.rel)
names(pantoea.rel)[1] <- "Variovorax"


###Fungi
fun.clean.log.nigro <- subset_taxa(fun.clean.rel, Genus == "Nigrospora")

nigro.rel<-colSums(otu_table(fun.clean.log.nigro))

nigro.rel <- data.frame(nigro.rel)
names(nigro.rel)[1] <- "Nigrospora"


fun.clean.log.alter <- subset_taxa(fun.clean.rel, Genus == "Alternaria")

alter.rel<-colSums(otu_table(fun.clean.log.alter))

alter.rel <- data.frame(alter.rel)
names(alter.rel)[1] <- "Alternaria"

fun.clean.log.moeszio <- subset_taxa(fun.clean.rel, Genus == "Moesziomyces")

moeszio.rel<-colSums(otu_table(fun.clean.log.moeszio))

moeszio.rel <- data.frame(moeszio.rel)
names(moeszio.rel)[1] <- "Moesziomyces"

fun.clean.log.clado <- subset_taxa(fun.clean.rel, Genus == "Cladosporium")

clado.rel<-colSums(otu_table(fun.clean.log.clado))

clado.rel <- data.frame(clado.rel)
names(clado.rel)[1] <- "Cladosporium"




fun.clean.log.alter <- subset_taxa(fun.clean.rel, Genus == "Guehomyces")

alter.rel<-colSums(otu_table(fun.clean.log.alter))

alter.rel <- data.frame(alter.rel)
names(alter.rel)[1] <- "Guehomyces"

fun.clean.log.clado <- subset_taxa(fun.clean.rel, Genus == "Chaetomium")

clado.rel<-colSums(otu_table(fun.clean.log.clado))

clado.rel <- data.frame(clado.rel)
names(clado.rel)[1] <- "Chaetomium"


meta.18.bac<-sample_data(bac.clean.nolog.18)
meta.18.bac <- data.frame(meta.18.bac)

genus.bac.rel<-merge(pantoea.rel, methlyo.rel,by = 'row.names')
rownames(genus.bac.rel) <-genus.bac.rel$Row.names
genus.bac.rel<-genus.bac.rel[-c(1)]
genus.bac.rel<-merge(genus.bac.rel, xantho.rel,by = 'row.names')
rownames(genus.bac.rel) <-genus.bac.rel$Row.names
genus.bac.rel<-genus.bac.rel[-c(1)]

genus.bac.rel<-merge(genus.bac.rel, pseudo.rel,by = 'row.names')
names(genus.bac.rel)[1] <- "SampleID"

meta.18.bac.rel<-merge(meta.18.bac,genus.bac.rel,by = 'SampleID')
rownames(meta.18.bac.rel) <-meta.18.bac.rel$SampleID

sample_data(bac.clean.log.18) <- sample_data(meta.18.bac.rel)

plot_ordination(bac.clean.log.18, bray2.bac, type = "samples", color='Xanthomonas', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#CCCCCC",
                        high = "#00ffff",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())





meta.18.fun<-sample_data(fun.clean.nolog.18)
meta.18.fun <- data.frame(meta.18.fun)

genus.fun.rel<-merge(nigro.rel, alter.rel,by = 'row.names')
rownames(genus.fun.rel) <-genus.fun.rel$Row.names
genus.fun.rel<-genus.fun.rel[-c(1)]
genus.fun.rel<-merge(genus.fun.rel, moeszio.rel,by = 'row.names')
rownames(genus.fun.rel) <-genus.fun.rel$Row.names
genus.fun.rel<-genus.fun.rel[-c(1)]

genus.fun.rel<-merge(genus.fun.rel, clado.rel,by = 'row.names')
names(genus.fun.rel)[1] <- "SampleID"


meta.18.fun.rel<-merge(meta.18.fun,genus.fun.rel,by = 'SampleID')
rownames(meta.18.fun.rel) <-meta.18.fun.rel$SampleID

sample_data(fun.clean.log.18) <- sample_data(meta.18.fun.rel)

plot_ordination(fun.clean.log.18, bray2.fun, type = "samples", color='Chaetomium', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#CCCCCC",
                        high = "#cc33cc",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())





####
bac.core <- read.xlsx("dis.core.bac_new.xlsx",1)

bac.list.core <- subset(bac.list, OTU_id %in% bac.core$OTU_id)
write.csv(bac.list.core,"taxonomy of core bac.csv")


fun.core <- read.xlsx("dis.core.fun_new.xlsx",1)

fun.list.core <- subset(fun.list, OTU_id %in% fun.core$OTU_id)
write.csv(fun.list.core,"taxonomy of core fun.csv")