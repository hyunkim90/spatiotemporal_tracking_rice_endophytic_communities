### Distance-decay relationship in rice plant
library(betapart)
decay.model()

###Example
## UTM Coordinates (in metres)
data(BCI)
class(BCI)
UTM.EW <- rep(seq(625754, 626654, by=100), each=5)
UTM.NS <- rep(seq(1011569,  1011969, by=100), len=50)

spat.dist<-dist(data.frame(UTM.EW, UTM.NS))

dissim.BCI<-beta.pair.abund(BCI)$beta.bray.bal

plot(spat.dist, dissim.BCI, ylim=c(0,1), xlim=c(0, max(spat.dist)))

BCI.decay.exp<-decay.model(dissim.BCI, spat.dist, y.type="dissim", model.type="exp", perm=100)
BCI.decay.exp<-decay.model(dissim.BCI, spat.dist, y.type="sim", model.type="exp", perm=100)
BCI.decay.pow<-decay.model(dissim.BCI, spat.dist, y.type="dissim", model.type="pow", perm=100)

plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
plot.decay(BCI.decay.pow, col="blue", remove.dots=TRUE, add=TRUE)



### Our data
### Days
b.meta.18.L <- sample_data(bac.clean.log.18.f.leaf)
b.meta.18.L.day <- b.meta.18.L[,c('Days')]
b.meta.18.L.day$Days <- as.numeric(as.character(b.meta.18.L.day$Days))
temp.dist<-dist(b.meta.18.L.day)

otu.bac.log.leaf.18<-otu_table(bac.clean.log.18.f.leaf)
otu.bac.log.leaf.18<-data.frame(otu.bac.log.leaf.18)
dissim.bac.leaf<-beta.pair.abund(t(otu.bac.log.leaf.18))$beta.bray.bal

plot(temp.dist, dissim.bac.leaf, ylim=c(0,1), xlim=c(0, max(temp.dist)))


bac.leaf.decay.exp<-decay.model(dissim.bac.leaf, temp.dist, y.type="sim", model.type="exp", perm=100)


plot.decay(bac.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
theme(aspect.ratio = 1)
dev.off()

##Stem
b.meta.18.S <- sample_data(bac.clean.log.18.f.stem)
b.meta.18.S.day <- b.meta.18.S[,c('Days')]
b.meta.18.S.day$Days <- as.numeric(as.character(b.meta.18.S.day$Days))
temp.dist<-dist(b.meta.18.S.day)

otu.bac.log.stem.18<-otu_table(bac.clean.log.18.f.stem)
otu.bac.log.stem.18<-data.frame(otu.bac.log.stem.18)
dissim.bac.stem<-beta.pair.abund(t(otu.bac.log.stem.18))$beta.bray.bal

plot(temp.dist, dissim.bac.stem, ylim=c(0,1), xlim=c(0, max(temp.dist)))


bac.stem.decay.exp<-decay.model(dissim.bac.stem, temp.dist, y.type="sim", model.type="exp", perm=100)


plot.decay(bac.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
theme(aspect.ratio = 1)
dev.off()


### Topological position
b.meta.18.L <- sample_data(bac.clean.log.18.f.leaf)
b.meta.18.L$Position <- ifelse(b.meta.18.L$Microhabitat == "L1", 20, ifelse(b.meta.18.L$Microhabitat == "L2",40, ifelse(b.meta.18.L$Microhabitat == "L3",60,80)))
b.meta.18.L.position <- b.meta.18.L[,c('Position')]
b.meta.18.L.position$Position <- as.numeric(as.character(b.meta.18.L.position$Position))
position.dist<-dist(b.meta.18.L.position)

otu.bac.log.leaf.18<-otu_table(bac.clean.log.18.f.leaf)
otu.bac.log.leaf.18<-data.frame(otu.bac.log.leaf.18)
dissim.bac.leaf<-beta.pair.abund(t(otu.bac.log.leaf.18))$beta.bray.bal

plot(position.dist, dissim.bac.leaf, ylim=c(0,1), xlim=c(0, max(position.dist)))

bac.leaf.decay.exp<-decay.model(dissim.bac.leaf, position.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(bac.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
bac.leaf.decay.exp$pseudo.r.squared


##Stem
b.meta.18.S <- sample_data(bac.clean.log.18.f.stem)
b.meta.18.S$Position <- ifelse(b.meta.18.S$Microhabitat == "S1", 10, ifelse(b.meta.18.S$Microhabitat == "S2",20,
                                                                            ifelse(b.meta.18.S$Microhabitat == "S3",30,ifelse(b.meta.18.S$Microhabitat == "S4",40,
                                                                                                                              ifelse(b.meta.18.S$Microhabitat == "S5",50,ifelse(b.meta.18.S$Microhabitat == "S6",60,
                                                                                                                                                                                ifelse(b.meta.18.S$Microhabitat == "S7",70,ifelse(b.meta.18.S$Microhabitat == "S8",80,90))))))))
b.meta.18.S.position <- b.meta.18.S[,c('Position')]
b.meta.18.S.position$Position <- as.numeric(as.character(b.meta.18.S.position$Position))
position.dist<-dist(b.meta.18.S.position)

otu.bac.log.stem.18<-otu_table(bac.clean.log.18.f.stem)
otu.bac.log.stem.18<-data.frame(otu.bac.log.stem.18)
dissim.bac.stem<-beta.pair.abund(t(otu.bac.log.stem.18))$beta.bray.bal

plot(position.dist, dissim.bac.stem, ylim=c(0,1), xlim=c(0, max(position.dist)))

bac.stem.decay.exp<-decay.model(dissim.bac.stem, position.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(bac.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
bac.stem.decay.exp$pseudo.r.squared

### Specific age
decay_leaf_age<-function(physeq, dayinteger){b.meta.18.L <- sample_data(physeq)
b.meta.18.L$Position <- ifelse(b.meta.18.L$Microhabitat == "L1", 20, ifelse(b.meta.18.L$Microhabitat == "L2",40, ifelse(b.meta.18.L$Microhabitat == "L3",60,80)))
b.meta.18.L.106 <- subset(b.meta.18.L, Days == dayinteger)
b.meta.18.L.106$Position <- as.numeric(as.character(b.meta.18.L.106$Position))
b.meta.18.L.106 <- b.meta.18.L.106[,c('Position')]
L106.dist<-dist(b.meta.18.L.106)

bac.clean.log.18.leaf.106 <- subset_samples(bac.clean.log.18.f.leaf, Days == dayinteger)
otu.bac.log.leaf.18<-otu_table(bac.clean.log.18.leaf.106)
otu.bac.log.leaf.18<-data.frame(otu.bac.log.leaf.18)
dissim.bac.leaf<-beta.pair.abund(t(otu.bac.log.leaf.18))$beta.bray.bal

plot(L106.dist, dissim.bac.leaf, ylim=c(0,1), xlim=c(0, max(L106.dist)))

bac.leaf.decay.exp<-decay.model(dissim.bac.leaf, L106.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(bac.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
print(bac.leaf.decay.exp)}

dayinteger = 141
decay_leaf_age(bac.clean.log.18.f.leaf, 141)

dayinteger = 120
decay_leaf_age(bac.clean.log.18.f.leaf, 120)

dayinteger = 106
decay_leaf_age(bac.clean.log.18.f.leaf, 106)


dayinteger = 90
decay_leaf_age(bac.clean.log.18.f.leaf, 90)

dayinteger = 76
decay_leaf_age(bac.clean.log.18.f.leaf, 76)

dayinteger = 62
decay_leaf_age(bac.clean.log.18.f.leaf, 62)

dayinteger = 48
decay_leaf_age(bac.clean.log.18.f.leaf, 48)



##Stem
decay_stem_age<-function(physeq, dayinteger){b.meta.18.S <- sample_data(physeq)
b.meta.18.S$Position <- ifelse(b.meta.18.S$Microhabitat == "S1", 10, ifelse(b.meta.18.S$Microhabitat == "S2",20,
                                                                            ifelse(b.meta.18.S$Microhabitat == "S3",30,ifelse(b.meta.18.S$Microhabitat == "S4",40,
                                                                                                                              ifelse(b.meta.18.S$Microhabitat == "S5",50,ifelse(b.meta.18.S$Microhabitat == "S6",60,
                                                                                                                                                                                ifelse(b.meta.18.S$Microhabitat == "S7",70,ifelse(b.meta.18.S$Microhabitat == "S8",80,90))))))))
b.meta.18.S.106 <- subset(b.meta.18.S, Days == dayinteger)
b.meta.18.S.106$Position <- as.numeric(as.character(b.meta.18.S.106$Position))
b.meta.18.S.106 <- b.meta.18.S.106[,c('Position')]
S106.dist<-dist(b.meta.18.S.106)

bac.clean.log.18.stem.106 <- subset_samples(bac.clean.log.18.f.stem, Days == dayinteger)
otu.bac.log.stem.18<-otu_table(bac.clean.log.18.stem.106)
otu.bac.log.stem.18<-data.frame(otu.bac.log.stem.18)
dissim.bac.stem<-beta.pair.abund(t(otu.bac.log.stem.18))$beta.bray.bal

plot(S106.dist, dissim.bac.stem, ylim=c(0,1), xlim=c(0, max(S106.dist)))

bac.stem.decay.exp<-decay.model(dissim.bac.stem, S106.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(bac.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
print(bac.stem.decay.exp)}

dayinteger = 141
decay_stem_age(bac.clean.log.18.f.stem, 141)

dayinteger = 120
decay_stem_age(bac.clean.log.18.f.stem, 120)

dayinteger = 106
decay_stem_age(bac.clean.log.18.f.stem, 106)


dayinteger = 90
decay_stem_age(bac.clean.log.18.f.stem, 90)

dayinteger = 76
decay_stem_age(bac.clean.log.18.f.stem, 76)

dayinteger = 62
decay_stem_age(bac.clean.log.18.f.stem, 62)

dayinteger = 48
decay_stem_age(bac.clean.log.18.f.stem, 48)


####Fungi
### Days
f.meta.18.L <- sample_data(fun.clean.log.18.f.leaf)
f.meta.18.L.day <- f.meta.18.L[,c('Days')]
f.meta.18.L.day$Days <- as.numeric(as.character(f.meta.18.L.day$Days))
temp.dist<-dist(f.meta.18.L.day)

otu.fun.log.leaf.18<-otu_table(fun.clean.log.18.f.leaf)
otu.fun.log.leaf.18<-data.frame(otu.fun.log.leaf.18)
dissim.fun.leaf<-beta.pair.abund(t(otu.fun.log.leaf.18))$beta.bray.bal

plot(temp.dist, dissim.fun.leaf, ylim=c(0,1), xlim=c(0, max(temp.dist)))


fun.leaf.decay.exp<-decay.model(dissim.fun.leaf, temp.dist, y.type="sim", model.type="exp", perm=100)


plot.decay(fun.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)

##Stem
f.meta.18.S <- sample_data(fun.clean.log.18.f.stem)
f.meta.18.S.day <- f.meta.18.S[,c('Days')]
f.meta.18.S.day$Days <- as.numeric(as.character(f.meta.18.S.day$Days))
temp.dist<-dist(f.meta.18.S.day)

otu.fun.log.stem.18<-otu_table(fun.clean.log.18.f.stem)
otu.fun.log.stem.18<-data.frame(otu.fun.log.stem.18)
dissim.fun.stem<-beta.pair.abund(t(otu.fun.log.stem.18))$beta.bray.bal

bray.dist<-vegdist(t(otu.fun.log.stem.18), method="bray", binary=T, diag=FALSE, upper=FALSE, na.rm = T)

plot(temp.dist, bray.dist, ylim=c(0,1), xlim=c(0, max(temp.dist)))


fun.stem.decay.exp<-decay.model(bray.dist, temp.dist, y.type="sim", model.type="exp", perm=100)
fun.stem.decay.pow<-decay.model(bray.dist, temp.dist, y.type="similarities", model.type="pow", perm=100)


plot.decay(fun.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
theme(aspect.ratio = 1)
dev.off()


### Topological position
f.meta.18.L <- sample_data(fun.clean.log.18.f.leaf)
f.meta.18.L$Position <- ifelse(f.meta.18.L$Microhabitat == "L1", 20, ifelse(f.meta.18.L$Microhabitat == "L2",40, ifelse(f.meta.18.L$Microhabitat == "L3",60,80)))
f.meta.18.L.position <- f.meta.18.L[,c('Position')]
f.meta.18.L.position$Position <- as.numeric(as.character(f.meta.18.L.position$Position))
position.dist<-dist(f.meta.18.L.position)

otu.fun.log.leaf.18<-otu_table(fun.clean.log.18.f.leaf)
otu.fun.log.leaf.18<-data.frame(otu.fun.log.leaf.18)
dissim.fun.leaf<-beta.pair.abund(t(otu.fun.log.leaf.18))$beta.bray.bal

plot(position.dist, dissim.fun.leaf, ylim=c(0,1), xlim=c(0, max(position.dist)))

fun.leaf.decay.exp<-decay.model(dissim.fun.leaf, position.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(fun.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
fun.leaf.decay.exp$pseudo.r.squared


##Stem
f.meta.18.S <- sample_data(fun.clean.log.18.f.stem)
f.meta.18.S$Position <- ifelse(f.meta.18.S$Microhabitat == "S1", 10, ifelse(f.meta.18.S$Microhabitat == "S2",20,
                                                                            ifelse(f.meta.18.S$Microhabitat == "S3",30,ifelse(f.meta.18.S$Microhabitat == "S4",40,
                                                                                                                              ifelse(f.meta.18.S$Microhabitat == "S5",50,ifelse(f.meta.18.S$Microhabitat == "S6",60,
                                                                                                                                                                                ifelse(f.meta.18.S$Microhabitat == "S7",70,ifelse(f.meta.18.S$Microhabitat == "S8",80,90))))))))
f.meta.18.S.position <- f.meta.18.S[,c('Position')]
f.meta.18.S.position$Position <- as.numeric(as.character(f.meta.18.S.position$Position))
position.dist<-dist(f.meta.18.S.position)

bray.dist<-vegdist(t(otu.fun.log.stem.18), method="bray", binary=T, diag=FALSE, upper=FALSE, na.rm = T)

otu.fun.log.stem.18<-otu_table(fun.clean.log.18.f.stem)
otu.fun.log.stem.18<-data.frame(otu.fun.log.stem.18)
dissim.fun.stem<-beta.pair.abund(t(otu.fun.log.stem.18))$beta.bray.bal

plot(position.dist, bray.dist, ylim=c(0,1), xlim=c(0, max(position.dist)))

fun.stem.decay.exp<-decay.model(bray.dist, position.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(fun.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
fun.stem.decay.exp$pseudo.r.squared

### Specific age
decay_leaf_age<-function(physeq, dayinteger){f.meta.18.L <- sample_data(physeq)
f.meta.18.L$Position <- ifelse(f.meta.18.L$Microhabitat == "L1", 20, ifelse(f.meta.18.L$Microhabitat == "L2",40, ifelse(f.meta.18.L$Microhabitat == "L3",60,80)))
f.meta.18.L.106 <- subset(f.meta.18.L, Days == dayinteger)
f.meta.18.L.106$Position <- as.numeric(as.character(f.meta.18.L.106$Position))
f.meta.18.L.106 <- f.meta.18.L.106[,c('Position')]
L106.dist<-dist(f.meta.18.L.106)

fun.clean.log.18.leaf.106 <- subset_samples(fun.clean.log.18.f.leaf, Days == dayinteger)
otu.fun.log.leaf.18<-otu_table(fun.clean.log.18.leaf.106)
otu.fun.log.leaf.18<-data.frame(otu.fun.log.leaf.18)
dissim.fun.leaf<-beta.pair.abund(t(otu.fun.log.leaf.18))$beta.bray.bal

plot(L106.dist, dissim.fun.leaf, ylim=c(0,1), xlim=c(0, max(L106.dist)))

fun.leaf.decay.exp<-decay.model(dissim.fun.leaf, L106.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(fun.leaf.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
print(fun.leaf.decay.exp)}

dayinteger = 141
decay_leaf_age(fun.clean.log.18.f.leaf, 141)

dayinteger = 120
decay_leaf_age(fun.clean.log.18.f.leaf, 120)

dayinteger = 106
decay_leaf_age(fun.clean.log.18.f.leaf, 106)


dayinteger = 90
decay_leaf_age(fun.clean.log.18.f.leaf, 90)

dayinteger = 76
decay_leaf_age(fun.clean.log.18.f.leaf, 76)

dayinteger = 62
decay_leaf_age(fun.clean.log.18.f.leaf, 62)

dayinteger = 48
decay_leaf_age(fun.clean.log.18.f.leaf, 48)



##Stem
decay_stem_age<-function(physeq, dayinteger){f.meta.18.S <- sample_data(physeq)
f.meta.18.S$Position <- ifelse(f.meta.18.S$Microhabitat == "S1", 10, ifelse(f.meta.18.S$Microhabitat == "S2",20,
                                                                            ifelse(f.meta.18.S$Microhabitat == "S3",30,ifelse(f.meta.18.S$Microhabitat == "S4",40,
                                                                                                                              ifelse(f.meta.18.S$Microhabitat == "S5",50,ifelse(f.meta.18.S$Microhabitat == "S6",60,
                                                                                                                                                                                ifelse(f.meta.18.S$Microhabitat == "S7",70,ifelse(f.meta.18.S$Microhabitat == "S8",80,90))))))))
f.meta.18.S.106 <- subset(f.meta.18.S, Days == dayinteger)
f.meta.18.S.106$Position <- as.numeric(as.character(f.meta.18.S.106$Position))
f.meta.18.S.106 <- f.meta.18.S.106[,c('Position')]
S106.dist<-dist(f.meta.18.S.106)

fun.clean.log.18.stem.106 <- subset_samples(fun.clean.log.18.f.stem, Days == dayinteger)
otu.fun.log.stem.18<-otu_table(fun.clean.log.18.stem.106)
otu.fun.log.stem.18<-data.frame(otu.fun.log.stem.18)
#dissim.fun.stem<-beta.pair.abund(t(otu.fun.log.stem.18))$beta.bray.bal

bray.dist<-vegdist(t(otu.fun.log.stem.18), method="bray", binary=T, diag=FALSE, upper=FALSE, na.rm = T)


plot(S106.dist, bray.dist, ylim=c(0,1), xlim=c(0, max(S106.dist)))

fun.stem.decay.exp<-decay.model(bray.dist, S106.dist, y.type="sim", model.type="exp", perm=100)
plot.decay(fun.stem.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
print(fun.stem.decay.exp)}

dayinteger = 141
decay_stem_age(fun.clean.log.18.f.stem, 141)

dayinteger = 120
decay_stem_age(fun.clean.log.18.f.stem, 120)

dayinteger = 106
decay_stem_age(fun.clean.log.18.f.stem, 106)


dayinteger = 90
decay_stem_age(fun.clean.log.18.f.stem, 90)

dayinteger = 76
decay_stem_age(fun.clean.log.18.f.stem, 76)

dayinteger = 62
decay_stem_age(fun.clean.log.18.f.stem, 62)

dayinteger = 48
decay_stem_age(fun.clean.log.18.f.stem, 48)




####Aboveground compartment
## Distance decay relationship

bac.clean.log.18.above.woG0

fun.clean.log.18.above


f.meta.18.above <- sample_data(fun.clean.log.18.above)
f.meta.18.above.day <- f.meta.18.above[,c('Days')]
f.meta.18.above.day$Days <- as.numeric(as.character(f.meta.18.above.day$Days))
temp.dist<-dist(f.meta.18.above.day)

otu.fun.log.above.18<-otu_table(fun.clean.log.18.above)
otu.fun.log.above.18<-data.frame(otu.fun.log.above.18)
dissim.fun.above<-beta.pair.abund(t(otu.fun.log.above.18))$beta.bray.bal



fun.above.decay.exp<-decay.model(dissim.fun.above, temp.dist, y.type="similarities", model.type="exp", perm=100)
fun.above.decay.pow<-decay.model(dissim.fun.above, temp.dist, y.type="similarities", model.type="pow", perm=100)

plot(temp.dist, dissim.fun.above, ylim=c(0,1), xlim=c(0, max(temp.dist)))

temp.dist.melt<-melt(as.matrix(temp.dist))
names(temp.dist.melt)[3] <- "Temporal_distance"
dissim.fun.above.melt<-melt(as.matrix(dissim.fun.above))
names(dissim.fun.above.melt)[3] <- "Dissimilarity"

distance.decay.above.fun<-merge(temp.dist.melt,dissim.fun.above.melt, by = c('Var1'="Var1", 'Var2'='Var2'))

ggplot(distance.decay.above.fun, aes(x=Temporal_distance, y=Dissimilarity)) +
  xlab('\n Temporal distance')+
  ylab("Community dissimilarity \n") +
  geom_point(size=2, alpha=0.7) +
  #scale_colour_manual(labels = c('high','medium','low'), values = c("#CC9900", "#0066CC", '#336633'))+
  #scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 0.5)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) +
  theme(plot.title = element_text(size = 20,hjust = 0.5)) +
  #geom_hline(yintercept=quantile(df.norm.degree$y2018_all_deg,prob=1-1/100) , color="maroon4", linetype='dotted')+
  #geom_vline(xintercept=quantile(df.norm.degree$y2017_all_deg,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(0,141,20))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) +
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) +
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+

  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


bray.dist<-vegdist(t(otu.fun.log.above.18), method="bray", binary=T, diag=FALSE, upper=FALSE, na.rm = T)

plot(temp.dist, bray.dist, ylim=c(0,1), xlim=c(0, max(temp.dist)))


fun.above.decay.exp<-decay.model(bray.dist, temp.dist, y.type="dissimilarities", model.type="exp", perm=100)
fun.above.decay.pow<-decay.model(bray.dist, temp.dist, y.type="similarities", model.type="exp", perm=100)

plot.new()
plot.decay(fun.above.decay.exp, col="red", remove.dots=TRUE, add=TRUE)

dev.off()
