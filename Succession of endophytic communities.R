#### Succession of microbiomes
### Seed samples
bac.clean.ss.18.G
fun.clean.ss.18.G

##Bacteria
bac.G.rel <- transform(bac.clean.ss.18.G, "compositional")

melt.bac.G <- psmelt(bac.G.rel)
head(melt.bac.G)
melt.bac.G$Days <- as.numeric(melt.bac.G$Days)
melt.bac.G <- subset(melt.bac.G, Days != 0)



mods <- melt.bac.G %>% group_by(OTU, Replication) %>% mutate(MeanAbund = mean(Abundance)) %>% 
  group_by(OTU) %>% 
  do(mod = summary(lm(MeanAbund ~ Days, .))) %>%
  mutate(
    pval = mod$coefficients[[8]],
    tval = mod$coefficients[[6]],
    group = ifelse(pval < 0.05 & tval < 0, 'Early successional', 'Mid-successional / No trend'),
    group = ifelse(pval < 0.05 & tval > 0, 'Late successional', group),
    group = factor(group, levels = c('Early successional','Mid-successional / No trend','Late successional'))) %>%
  select(OTU, group, pval, tval)

write.csv(mods,"Succession prediction_bacteria.csv")

mods.inherit <- subset(mods, OTU %in% persistent.bac)
length(mods.inherit$OTU[which(mods.inherit$group == "Late successional")])
length(mods.inherit$OTU[which(mods.inherit$group == "Early successional")])

mods.transient <- subset(mods, OTU %in% transient.bac)
length(mods.transient$OTU[which(mods.transient$group == "Late successional")])
length(mods.transient$OTU[which(mods.transient$group == "Early successional")])

###Fungi
fun.G.rel <- transform(fun.clean.ss.18.G, "compositional")
class(melt.fun.G$Days)
melt.fun.G <- psmelt(fun.G.rel)
melt.fun.G$Days <- as.numeric(as.character(melt.fun.G$Days))
melt.fun.G <- subset(melt.fun.G, Days != 0)



mods.f <- melt.fun.G %>% group_by(OTU, Replication) %>% mutate(MeanAbund = mean(Abundance)) %>% 
  group_by(OTU) %>% 
  do(mod = summary(lm(MeanAbund ~ Days, .))) %>%
  mutate(
    pval = mod$coefficients[[8]],
    tval = mod$coefficients[[6]],
    group = ifelse(pval < 0.05 & tval < 0, 'Early successional', 'Mid-successional / No trend'),
    group = ifelse(pval < 0.05 & tval > 0, 'Late successional', group),
    group = factor(group, levels = c('Early successional','Mid-successional / No trend','Late successional'))) %>%
  select(OTU, group, pval, tval)


write.csv(mods.f,"Succession prediction_fungi.csv")

tail(mods.f)
mods.f.inherit <- subset(mods.f, OTU %in% persistent.fun)
length(mods.f.inherit$OTU[which(mods.f.inherit$group == "Late successional")])
length(mods.f.inherit$OTU[which(mods.f.inherit$group == "Early successional")])


mods.f.transient <- subset(mods.f, OTU %in% transient.fun)
length(mods.f.transient$OTU[which(mods.f.transient$group == "Late successional")])
length(mods.f.transient$OTU[which(mods.f.transient$group == "Early successional")])


mods.f.inherit.tax<-merge(fun.list,mods.f.inherit, by = "OTU")
mods.inherit.tax<-merge(bac.list,mods.inherit, by = "OTU")

write.csv(mods.inherit.tax, "tax table and successional group_bacteria.csv")
write.csv(mods.f.inherit.tax, "tax table and successional group_fungi.csv")


### Abundance plot
bac.G.rel <- transform(bac.clean.ss.18.G, "compositional")
melt.bac.G <- psmelt(bac.G.rel)
head(melt.bac.G)
melt.bac.G$Days <- as.numeric(melt.bac.G$Days)

j <- melt.bac.G %>%
  left_join(mods, by = 'OTU') %>%
  group_by(group, Days) %>%
  ungroup()
j$Days <- as.factor(as.character(j$Days))
j$Days <- factor(j$Days, levels = c('76','90','106','120','141', "0"))

j$Inheritance <- ifelse(j$OTU %in% persistent.bac, "Inherited", ifelse(j$OTU %in% transient.bac, "Transient","Non-inherited"))

j2<-j %>%
  group_by(Days, group,OTU,Inheritance) %>% summarise(Mean = mean(Abundance)) %>% group_by(Days, group,Inheritance) %>% summarise(TotalAbund = sum(Mean))

p1<-j2%>%filter(group == 'Early successional') %>% ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Relative abundance', tag = 'a') +
  ggtitle('Early successional') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank())+theme(aspect.ratio = 0.3)

p3 <- j2%>%
  filter(group == 'Mid-successional / No trend') %>%
  ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Relative abundance', tag = 'b') +
  ggtitle('Mid-successional / No trend') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank())+theme(aspect.ratio = 0.3)


p5 <- j2%>%
  filter(group == 'Late successional') %>%
  ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Days after transplanting', y = 'Relative abundance', tag = 'c') +
  ggtitle('Late successional') +
  theme_classic() +
  theme(legend.position = 'none')+theme(aspect.ratio = 0.3)



fun.G.rel <- transform(fun.clean.ss.18.G, "compositional")
class(melt.fun.G$Days)
melt.fun.G <- psmelt(fun.G.rel)
melt.fun.G$Days <- as.numeric(as.character(melt.fun.G$Days))

j3 <- melt.fun.G %>%
  left_join(mods.f, by = 'OTU') %>%
  group_by(group, Days) %>%
  ungroup()
j3$Days <- as.factor(as.character(j3$Days))
j3$Days <- factor(j3$Days, levels = c('76','90','106','120','141',"0"))


j3$Inheritance <- ifelse(j3$OTU %in% persistent.fun, "Inherited", ifelse(j3$OTU %in% transient.fun, "Transient","Non-inherited"))

j5<-j3 %>%
  group_by(Days, group,OTU,Inheritance) %>% summarise(Mean = mean(Abundance)) %>% group_by(Days, group,Inheritance) %>% summarise(TotalAbund = sum(Mean))


p2<-j5%>%filter(group == 'Early successional') %>% ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Relative abundance', tag = 'a') +
  ggtitle('Early successional') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank())+theme(aspect.ratio = 0.3)

p4 <- j5%>%
  filter(group == 'Mid-successional / No trend') %>%
  ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Relative abundance', tag = 'b') +
  ggtitle('Mid-successional / No trend') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank())+theme(aspect.ratio = 0.3)


p6 <- j5%>%
  filter(group == 'Late successional') %>%
  ggplot(aes(x = Days, y = TotalAbund, fill = Inheritance)) +
  geom_bar(stat = 'identity', width = 1) +
  expand_limits(y = 1) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Days after transplanting', y = 'Relative abundance', tag = 'c') +
  ggtitle('Late successional') +
  theme_classic() +
  theme(legend.position = 'none')+theme(aspect.ratio = 0.3)





figTaxOverTime <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, align = 'h', ncol = 2)


write.csv(melt.bac.G, "raw data for succession plot.csv")
write.csv(melt.fun.G, "raw data for succession plot_fungi.csv")

write.csv(j2, "Source data_Fig. 2c_bacteria.csv")
write.csv(j5, "Source data_Fig. 2c_fungi.csv")

### Other compartments
succession.predic<-function(phy.seq,word1,kingd){
  phy.seq.comp <- subset_samples(phy.seq, Microhabitat == word1)
phy.seq.comp <- phyloseq::filter_taxa(phy.seq.comp, function(x) sum(x) != 0, TRUE)
phy.seq.comp <- transform(phy.seq.comp, "compositional")

melt.phy.seq.comp <- psmelt(phy.seq.comp)
head(melt.phy.seq.comp)
melt.phy.seq.comp$Days <- as.numeric(melt.phy.seq.comp$Days)

mods <- melt.phy.seq.comp %>% group_by(OTU, Replication) %>% mutate(MeanAbund = mean(Abundance)) %>% 
  group_by(OTU) %>% 
  do(mod = summary(lm(MeanAbund ~ Days, .))) %>%
  mutate(
    pval = mod$coefficients[[8]],
    tval = mod$coefficients[[6]],
    group = ifelse(pval < 0.05 & tval < 0, 'Early successional', 'Mid-successional / No trend'),
    group = ifelse(pval < 0.05 & tval > 0, 'Late successional', group),
    group = factor(group, levels = c('Early successional','Mid-successional / No trend','Late successional'))) %>%
  select(OTU, group, pval, tval)

write.csv(mods,paste0("Succession prediction_",word1,kingd,".csv"))

return(mods)  
}

word1 <-"L1"
kingd <- "Bacteria"
L1.bac.success <- succession.predic(bac.clean.ss.f, "L1","Bacteria")

word1 <-"L2"
L2.bac.success <- succession.predic(bac.clean.ss.f, "L2","Bacteria")

word1 <-"L3"
L3.bac.success <- succession.predic(bac.clean.ss.f, "L3","Bacteria")

word1 <-"FL"
FL.bac.success <- succession.predic(bac.clean.ss.f, "FL","Bacteria")

word1 <-"S1"
S1.bac.success <- succession.predic(bac.clean.ss.f, "S1","Bacteria")

word1 <-"S2"
S2.bac.success <- succession.predic(bac.clean.ss.f, "S2","Bacteria")

word1 <-"S3"
S3.bac.success <- succession.predic(bac.clean.ss.f, "S3","Bacteria")

word1 <-"S4"
S4.bac.success <- succession.predic(bac.clean.ss.f, "S4","Bacteria")

word1 <-"S5"
S5.bac.success <- succession.predic(bac.clean.ss.f, "S5","Bacteria")

word1 <-"S6"
S6.bac.success <- succession.predic(bac.clean.ss.f, "S6","Bacteria")

word1 <-"S7"
S7.bac.success <- succession.predic(bac.clean.ss.f, "S7","Bacteria")

word1 <-"S8"
S8.bac.success <- succession.predic(bac.clean.ss.f, "S8","Bacteria")

word1 <-"S9"
S9.bac.success <- succession.predic(bac.clean.ss.f, "S9","Bacteria")

word1 <-"L1"
kingd <- "Fungi"
L1.fun.success <- succession.predic(fun.clean.ss.f, "L1","Fungi")

word1 <-"L2"
L2.fun.success <- succession.predic(fun.clean.ss.f, "L2","Fungi")

word1 <-"L3"
L3.fun.success <- succession.predic(fun.clean.ss.f, "L3","Fungi")

word1 <-"FL"
FL.fun.success <- succession.predic(fun.clean.ss.f, "FL","Fungi")

word1 <-"S1"
S1.fun.success <- succession.predic(fun.clean.ss.f, "S1","Fungi")

word1 <-"S2"
S2.fun.success <- succession.predic(fun.clean.ss.f, "S2","Fungi")

word1 <-"S3"
S3.fun.success <- succession.predic(fun.clean.ss.f, "S3","Fungi")

word1 <-"S4"
S4.fun.success <- succession.predic(fun.clean.ss.f, "S4","Fungi")

word1 <-"S5"
S5.fun.success <- succession.predic(fun.clean.ss.f, "S5","Fungi")

word1 <-"S6"
S6.fun.success <- succession.predic(fun.clean.ss.f, "S6","Fungi")

word1 <-"S7"
S7.fun.success <- succession.predic(fun.clean.ss.f, "S7","Fungi")

word1 <-"S8"
S8.fun.success <- succession.predic(fun.clean.ss.f, "S8","Fungi")

word1 <-"S9"
S9.fun.success <- succession.predic(fun.clean.ss.f, "S9","Fungi")


### Inherited OTU
subset(S9.fun.success, OTU %in% persistent.fun)
subset(S8.fun.success, OTU %in% persistent.fun)
subset(S7.fun.success, OTU %in% persistent.fun)
subset(S6.fun.success, OTU %in% persistent.fun)
subset(S5.fun.success, OTU %in% persistent.fun)
subset(S4.fun.success, OTU %in% persistent.fun)
subset(S3.fun.success, OTU %in% persistent.fun)
subset(S2.fun.success, OTU %in% persistent.fun)
subset(S1.fun.success, OTU %in% persistent.fun)


subset(FL.fun.success, OTU %in% persistent.fun)
subset(L3.fun.success, OTU %in% persistent.fun)
subset(L2.fun.success, OTU %in% persistent.fun)
subset(L1.fun.success, OTU %in% persistent.fun)


subset(S9.bac.success, OTU %in% persistent.bac)
subset(S8.bac.success, OTU %in% persistent.bac)
subset(S7.bac.success, OTU %in% persistent.bac)
subset(S6.bac.success, OTU %in% persistent.bac)
subset(S5.bac.success, OTU %in% persistent.bac)
subset(S4.bac.success, OTU %in% persistent.bac)
subset(S3.bac.success, OTU %in% persistent.bac)
subset(S2.bac.success, OTU %in% persistent.bac)
subset(S1.bac.success, OTU %in% persistent.bac)


subset(FL.bac.success, OTU %in% persistent.bac)
subset(L3.bac.success, OTU %in% persistent.bac)
subset(L2.bac.success, OTU %in% persistent.bac)
subset(L1.bac.success, OTU %in% persistent.bac)



###  Supplementary Data file
succession.G.bac <- read.csv('Succession prediction_bacteria.csv')
succession.G.fun <- read.csv('Succession prediction_fungi.csv')


bac.tab.mode <- merge(succession.G.bac, bac.list, by = c('OTU'='OTU'))

bac.tab.mode$Inheritance <- ifelse(bac.tab.mode$OTU %in% persistent.bac, "Inherited", ifelse(bac.tab.mode$OTU %in% transient.bac, "Transient","Non-inherited"))
write.xlsx(bac.tab.mode, "Suppmentary Data 3. succcession mode_bacteria.xlsx")


fun.tab.mode <- merge(succession.G.fun, fun.list, by = c('OTU'='OTU'))

fun.tab.mode$Inheritance <- ifelse(fun.tab.mode$OTU %in% persistent.fun, "Inherited", ifelse(fun.tab.mode$OTU %in% transient.fun, "Transient","Non-inherited"))
write.xlsx(fun.tab.mode, "Suppmentary Data 3. succcession mode_fungi.xlsx")
