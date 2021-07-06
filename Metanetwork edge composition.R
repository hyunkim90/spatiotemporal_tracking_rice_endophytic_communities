#### Module compositions in the metanetwork
##Top 4
module_1 <- read.csv("Metanetwork_module_0_edge.csv")
module_3 <- read.csv("Metanetwork_module_2_edge.csv")
module_6 <- read.csv("Metanetwork_module_5_edge.csv")
module_11 <- read.csv("Metanetwork_module_10_edge.csv")

top4.module <- rbind(module_1,module_3,module_6,module_11)

module.count.compart<-read.csv("Module_compartment_count.csv")

##module 1
for (i in as.character(module.count.compart$Compartment)){
  module.count.compart$Module_1[which(module.count.compart$Compartment == i)]<-length(module_1$compartment[which(module_1$compartment == i)])
}


##module 3
for (i in as.character(module.count.compart$Compartment)){
  module.count.compart$Module_3[which(module.count.compart$Compartment == i)]<-length(module_3$compartment[which(module_3$compartment == i)])
}


##module 6
for (i in as.character(module.count.compart$Compartment)){
  module.count.compart$Module_6[which(module.count.compart$Compartment == i)]<-length(module_6$compartment[which(module_6$compartment == i)])
}


##module 11
for (i in as.character(module.count.compart$Compartment)){
  module.count.compart$Module_11[which(module.count.compart$Compartment == i)]<-length(module_11$compartment[which(module_11$compartment == i)])
}


## Compartment and age
module.count.compart.age<-read.csv("Module_compartment_count_2.csv")

##module 1
for (i in as.character(module.count.compart.age$Compartment)){
  module.count.compart.age$Module_1[which(module.count.compart.age$Compartment == i)]<-length(module_1$net_group[which(module_1$net_group == i)])
}


##module 3
for (i in as.character(module.count.compart.age$Compartment)){
  module.count.compart.age$Module_3[which(module.count.compart.age$Compartment == i)]<-length(module_3$net_group[which(module_3$net_group == i)])
}


##module 6
for (i in as.character(module.count.compart.age$Compartment)){
  module.count.compart.age$Module_6[which(module.count.compart.age$Compartment == i)]<-length(module_6$net_group[which(module_6$net_group == i)])
}


##module 11
for (i in as.character(module.count.compart.age$Compartment)){
  module.count.compart.age$Module_11[which(module.count.compart.age$Compartment == i)]<-length(module_11$net_group[which(module_11$net_group == i)])
}


write.csv(module.count.compart.age, "Distribution of edges in top 4 modules.csv")

#### Graph
###Simplified compartment

module.comp.simple <- read.csv("Edge composition_simple compartment.csv")
melt.compart.module <- melt(module.comp.simple)
melt.compart.module$Compartment <- factor(melt.compart.module$Compartment, levels = as.character(module.comp.simple$Compartment))
melt.compart.module.2<-melt.compart.module %>% group_by(variable) %>% mutate(Prop = value/sum(value))


ggplot(melt.compart.module.2, aes(x=variable, y = Prop, fill = Compartment)) + 
  geom_bar(stat = 'identity',width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())




### Compartment
melt.compart.module <- melt(module.count.compart)
melt.compart.module$Compartment <- factor(melt.compart.module$Compartment, levels = as.character(module.count.compart$Compartment))
melt.compart.module.2<-melt.compart.module %>% group_by(variable) %>% mutate(Prop = value/sum(value))


ggplot(melt.compart.module.2, aes(x=variable, y = Prop, fill = Compartment)) + 
  geom_bar(stat = 'identity',width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_collection) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


#### Compartment and Age
melt.compart.age.module <- melt(module.count.compart.age)
melt.compart.age.module$Compartment <- factor(melt.compart.age.module$Compartment, levels = as.character(module.count.compart.age$Compartment))
melt.compart.age.module.2<-melt.compart.age.module %>% group_by(variable) %>% mutate(Prop = value/sum(value))


ggplot(melt.compart.age.module.2, aes(x=variable, y = Prop, fill = Compartment)) + 
  geom_bar(stat = 'identity',width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  #scale_fill_manual(values = my_color_collection) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 10,reverse = F))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
