### Bacteria

pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR","minpack.lm", "Hmisc", "reltools", "stats4")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

Palette.phylum <-c("Gammaproteobacteria" = "darkolivegreen3", "Alphaproteobacteria"= "darkolivegreen", 
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
                   "Gemmatimonadetes"= "peachpuff3", "Other" = "light grey", "unidentified" = "black")
Palette.shape.1<-c(21, 21, 16)
Palette.shape.2<-c(17,15,19,13,8,7)
Palette.shape.3<-c(16,1)
alpha=0.05
minZ=1.5
maxZ=8.5

Stat<-list()
p<-list()

theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_rect(),
  axis.title.y=element_blank(),
  axis.title.x=element_blank(),
  panel.background = element_rect(fill="white", colour ="black")
)


## Normalized abundance
bac.clean.nolog.18 <- subset_samples(bac.clean.nolog.18, Replication != "G_0")
bac.clean.nolog.18<- phyloseq::filter_taxa(bac.clean.nolog.18, function(x) sum(x) != 0, TRUE)

bac.nolog.habitat <- bac.clean.nolog.18
DT=data.table(otu_table(bac.nolog.habitat), keep.rownames=T, key="rn")
tax <- tax_table(bac.clean.nolog.18)
tax <- data.frame(tax, stringsAsFactors = F)
OTU_B=t(otu_table(bac.clean.nolog.18))
TAXA_B=tax_table(bac.clean.nolog.18)

meta.all<-read.table('Metadata_bac.tsv',sep= '\t', header =T)
meta.18.all <-subset(meta.all, meta.all$SampleID %in% sample_names(bac.clean.nolog.18))

L1.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L1")])
L2.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L2")])
L3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L3")])
FL.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "FL")])
S1.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S1")])
S2.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S2")])
S3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S3")])
S4.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S4")])
S5.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S5")])
S6.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S6")])
S7.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S7")])
S8.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S8")])
S9.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S9")])
R.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "R")])
RS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "RS")])
BS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "BS")])
G.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "G")])

All<-list(L1.series, L2.series, L3.series,FL.series, 
          S1.series, S2.series,S3.series, S4.series, S5.series, S6.series,S7.series, S8.series, S9.series,
          R.series, RS.series, BS.series,G.series)

DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
DT.m=merge(DT,DT.taxa)

ColTaxa<-colnames(DT.taxa)

for (i in All) {
  
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn
  
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_B[i,LIST],taxon=TAXA_B[LIST,])
  
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  
  ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
  DT.sncm$ColPhylum <- ifelse(!DT.sncm$Phylum %in% ColPhylum, "Other", ifelse(DT.sncm$Phylum == "Bacteroidetes", "Bacteroidetes", ifelse(DT.sncm$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.sncm$Phylum == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  
  Richness<-sncm.out$fitstats$Richness
  
  fit.stat<-sncm.out$fitstats
  
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  
  Stat[[i[[1]]]]<-fit.stat
  print("######################################################################")
}

sncm.stat<-rbind(Stat$A1_L2_48_1B, Stat$A1_L1_48_1B, Stat$A1_L1_62_1B, Stat$A1_FL_76_1B, Stat$A1_S1_48_1B, 
                 Stat$A1_S2_48_1B, Stat$A1_S3_48_1B, Stat$A1_S4_62_1B, Stat$A1_S5_76_1B, Stat$A1_S6_76_1B, 
                 Stat$A1_S7_90_1B, Stat$A1_S8_90_1B, Stat$A1_S9_90_1B, Stat$A1_R_48_1B, Stat$A1_RS_48_1B, 
                 Stat$A1_BS_48_1B, Stat$A1_G_76_1B)
rownames(sncm.stat)<-c("L1", "L2", "L3", "FL", "S1", "S2", "S3", "S4", "S5", "S6","S7","S8","S9","R","RS","BS","G")
sncm.stat$comp<-c("L","L","L","L","S","S","S","S","S","S","S","S","S","R","RS","BS", "G")

write.table(sncm.stat, "sncm.stat_all replicates_normalized abundance_bacteria_final_211017.tsv", sep='\t', quote=F)

mid<-mean(sncm.stat$Richness)

pp<-ggplot(sncm.stat, aes(x=Rsqr, y=m, color=Richness, shape=comp))

panel_a=pp+geom_point(size=5)+
  scale_color_gradient2(low="darkcyan", mid="darkorange", high="darkviolet", midpoint=mid)+
  theme_bw() + theme_change + theme(aspect.ratio = 1)
scale_shape_manual(values=Palette.shape.2)

#	pdf("Fig3a.pdf", paper="A4", useDingbats=FALSE)
print(panel_a)
dev.off()


##Grouped by over- and down-represented and fitted OTUs
for (i in All) {
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn	
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_B[i,LIST],taxon=TAXA_B[LIST,])
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Verrucomicrobia","Alphaproteobacteria","Deltaproteobacteria","Gammaproteobacteria","Chloroflexi","Acidobacteria")
  DT.sncm$ColPhylum <- 0
  DT.sncm$ColPhylum[-which(DT.sncm$Phylum %in% ColPhylum)] <- "Other"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Alphaproteobacteria")] <- "Alphaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Gammaproteobacteria")] <- "Gammaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Deltaproteobacteria")] <- "Deltaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum =="Bacteroidetes")] <- "Bacteroidetes"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Actinobacteria")] <- "Actinobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Verrucomicrobia")] <- "Verrucomicrobia"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Firmicutes")] <- "Firmicutes"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Acidobacteria")] <- "Acidobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Chloroflexi")] <- "Chloroflexi" 
  
  DT.sncm$ColPhylum <- factor(DT.sncm$ColPhylum, levels = c("Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",
                                                            "Chloroflexi","Firmicutes","Acidobacteria", "Actinobacteria","Bacteroidetes","Verrucomicrobia","Other"))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  DT.sncm$fit_class_inherit <- ifelse(DT.sncm$rn %in% persistent.bac, "Inherit", "Others")
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  Richness<-sncm.out$fitstats$Richness
  fit.stat<-sncm.out$fitstats
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  Stat[[i[[1]]]]<-fit.stat
  DT.sncm$fit_class <- factor(DT.sncm$fit_class, levels = c("Above prediction", "As predicted", "Below prediction"))
  names(DT.sncm)[1] <- "OTU"
  DT.sncm <- merge(DT.sncm, OTU_id.list, by = "OTU")
  write.table(DT.sncm, paste0("Bacteria_DT.sncm_all reps_211017",i[[1]],".tsv"), sep = '\t', quote = F)
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  p1= pp1 + geom_jitter (aes(color=fit_class, shape = fit_class_inherit),size=2.5) +
    theme_bw() + theme_change + theme(aspect.ratio = 1)+
    scale_color_manual(values=c("#cc6666", "#999999", "#6699cc"))+
    scale_shape_manual(values=Palette.shape.3)+
    geom_line(aes(x=log(p), y=freq.pred), color = '#333333', size = 0.75)+
    geom_line(aes(x=log(p), y=pred.lwr),linetype="longdash",color = '#333333', size = 0.5)+
    geom_line(aes(x=log(p), y=pred.upr),linetype="longdash",color = '#333333', size = 0.5)+
    theme(axis.text.x = element_text(vjust=0.4,size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))
  p[[i[[1]]]]<-p1
  print("######################################################################")
}
p$A1_L2_48_1B
p$A1_L1_48_1B
p$A1_L1_62_1B
p$A1_FL_76_1B

p$A1_S1_48_1B
p$A1_S2_48_1B
p$A1_S3_48_1B
p$A1_S4_62_1B
p$A1_S5_76_1B
p$A1_S6_76_1B
p$A1_S7_90_1B
p$A1_S8_90_1B
p$A1_S9_90_1B

p$A1_R_48_1B
p$A1_RS_48_1B
p$A1_BS_48_1B
p$A1_G_76_1B

dev.off()


###Fungi
Palette.class <-c("Mortierellomycetes" = "#4E734E", "Eurotiomycetes"= "#6DA9DC",
                  "Mucoromycetes"="#E4AF2C","Tremellomycetes"="#BE4146", "Microbotryomycetes" ="#DC9A9E",
                  "Dothideomycetes" = "#5195D1", "Agaricomycetes" = "#CC6C71", "Blastocladiomycetes" = "#87AC88",
                  "Sordariomycetes"= "#1E63AF","Leotiomycetes"= "#11335F", "Cystobasidiomycetes" = "#A871AE", "Pezizomycetes" = "#C0DBF3",
                  "unidentified" ="#000000", "Endogonomycetes" ="#CC9900", "Pucciniomycetes" = "#0099FF", "Ustilaginomycetes"="#ffcc33","Arthoniomycetes" = "#cccc99",
                  "Saccharomycetes" = "#666699", "Other" = "light grey", "Malasseziomycetes" = "#99cc99", "Exobasidiomycetes" = "#663366", 
                  "Rhizophydiomycetes"="#ffccff","Wallemiomycetes"="#333300","Agaricostilbomycetes"="#003333","Orbiliomycetes"="#cc9966","Chytridiomycetes" = "#006600",
                  "Spizellomycetes" = "#66cc99","Moniliellomycetes" = "#0099cc","Monoblepharidomycetes" = "#CC0000", "Dacrymycetes" = "#669999",
                  "Spiculogloeomycetes" = "#99cc00")


fun.clean.nolog.18 <- subset_samples(fun.clean.nolog.18, Replication != "G_0")
fun.clean.nolog.18<- phyloseq::filter_taxa(fun.clean.nolog.18, function(x) sum(x) != 0, TRUE)

fun.nolog.habitat <- fun.clean.nolog.18
DT=data.table(otu_table(fun.nolog.habitat), keep.rownames=T, key="rn")
head(DT)

tax <- tax_table(fun.nolog.habitat)
tax <- data.frame(tax, stringsAsFactors = F)
OTU_F=t(otu_table(fun.nolog.habitat))
TAXA_F=tax_table(fun.nolog.habitat)

meta.18.all<-sample_data(fun.nolog.habitat)
meta.18.all$SampleID <- rownames(meta.18.all)


L1.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L1")])
L2.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L2")])
L3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L3")])
FL.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "FL")])
S1.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S1")])
S2.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S2")])
S3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S3")])
S4.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S4")])
S5.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S5")])
S6.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S6")])
S7.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S7")])
S8.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S8")])
S9.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S9")])
R.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "R")])
RS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "RS")])
BS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "BS")])
G.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "G")])

All<-list(L1.series, L2.series, L3.series,FL.series, 
          S1.series, S2.series,S3.series, S4.series, S5.series, S6.series,S7.series, S8.series, S9.series,
          R.series, RS.series, BS.series,G.series)

DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
DT.m=merge(DT,DT.taxa)

ColTaxa<-colnames(DT.taxa)

for (i in All) {
  
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn
  
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_F[i,LIST],taxon=TAXA_F[LIST,])
  
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  
  ColClass<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
  DT.sncm$ColClass <- ifelse(!DT.sncm$Class %in% ColClass, "Other", ifelse(DT.sncm$Class == "Bacteroidetes", "Bacteroidetes", ifelse(DT.sncm$Class == "Firmicutes", "Firmicutes", ifelse (DT.sncm$Class == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  
  Richness<-sncm.out$fitstats$Richness
  
  fit.stat<-sncm.out$fitstats
  
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  
  Stat[[i[[1]]]]<-fit.stat
  print("######################################################################")
}

sncm.stat.f<-rbind(Stat$A1_L2_48_1F, Stat$A1_L1_48_1F, Stat$A1_L1_106_1F, Stat$A1_FL_106_1F, Stat$A1_S1_106_1F, 
                   Stat$A1_S2_106_1F, Stat$A1_S3_106_1F, Stat$A1_S4_106_1F, Stat$A1_S5_106_1F, Stat$A1_S6_106_1F, 
                   Stat$A1_S7_106_1F, Stat$A1_S8_106_1F, Stat$A1_S9_120_1F, Stat$A1_R_106_1F, Stat$A1_RS_106_1F, 
                   Stat$A1_BS_106_1F, Stat$A1_G_106_1F)
rownames(sncm.stat.f)<-c("L1", "L2", "L3", "FL", "S1", "S2", "S3", "S4", "S5", "S6","S7","S8","S9","R","RS","BS","G")
sncm.stat.f$comp<-c("L","L","L","L","S","S","S","S","S","S","S","S","S","R","RS","BS", "G")

write.table(sncm.stat.f, "sncm.stat_all replicates_normalized abundance_fungi_final.tsv", sep='\t', quote=F)

mid<-mean(sncm.stat.f$Richness)

pp<-ggplot(sncm.stat.f, aes(x=Rsqr, y=m, color=Richness, shape=comp))

panel_a=pp+geom_point(size=5)+
  scale_color_gradient2(low="darkcyan", mid="darkorange", high="darkviolet", midpoint=mid)+
  theme_bw() + theme_change + theme(aspect.ratio = 1)
scale_shape_manual(values=Palette.shape.2)

#	pdf("Fig3a.pdf", paper="A4", useDingbats=FALSE)
print(panel_a)
dev.off()

Stat<-list()
p<-list()

##Grouped by over- and down-represented and fitted OTUs
for (i in All) {
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn	
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_F[i,LIST],taxon=TAXA_F[LIST,])
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  ColClass<-c("Tremellomycetes","Sordariomycetes","Leotiomycetes","Dothideomycetes","Agaricomycetes","Mortierellomycetes","Ustilaginomycetes","Eurotiomycetes","Saccharomycetes")
  DT.sncm$ColClass <- 0
  DT.sncm$ColClass[-which(DT.sncm$Class %in% ColClass)] <- "Other"
  DT.sncm$ColClass[which(DT.sncm$Class =="Tremellomycetes")] <- "Tremellomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Sordariomycetes")] <- "Sordariomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Leotiomycetes")] <- "Leotiomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Dothideomycetes")] <- "Dothideomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Agaricomycetes")] <- "Agaricomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Mortierellomycetes")] <- "Mortierellomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Ustilaginomycetes")] <- "Ustilaginomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Eurotiomycetes")] <- "Eurotiomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Saccharomycetes")] <- "Saccharomycetes" 
  
  DT.sncm$ColClass <- factor(DT.sncm$ColClass, levels = c("Sordariomycetes", "Dothideomycetes", "Leotiomycetes",
                                                          "Eurotiomycetes","Saccharomycetes","Tremellomycetes", "Agaricomycetes","Ustilaginomycetes","Mortierellomycetes","Other"))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  DT.sncm$fit_class_inherit <- ifelse(DT.sncm$rn %in% persistent.bac, "Inherit", "Others")
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  Richness<-sncm.out$fitstats$Richness
  fit.stat<-sncm.out$fitstats
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  Stat[[i[[1]]]]<-fit.stat
  DT.sncm$fit_class <- factor(DT.sncm$fit_class, levels = c("Above prediction", "As predicted", "Below prediction"))
  names(DT.sncm)[1] <- "OTU"
  DT.sncm <- merge(DT.sncm, OTU_id.list, by = "OTU")
  write.table(DT.sncm, paste0("Fungi_DT.sncm_merged_",i[[1]],".tsv"), sep = '\t', quote = F)
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  p1= pp1 + geom_jitter (aes(color=fit_class, shape = fit_class_inherit),size=2.5) +
    theme_bw() + theme_change + theme(aspect.ratio = 1)+
    scale_color_manual(values=c("#cc6666", "#999999", "#6699cc"))+
    scale_shape_manual(values=Palette.shape.3)+
    geom_line(aes(x=log(p), y=freq.pred), color = '#333333', size = 0.75)+
    geom_line(aes(x=log(p), y=pred.lwr),linetype="longdash",color = '#333333', size = 0.5)+
    geom_line(aes(x=log(p), y=pred.upr),linetype="longdash",color = '#333333', size = 0.5)+
    theme(axis.text.x = element_text(vjust=0.4,size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))
  p[[i[[1]]]]<-p1
  print("######################################################################")
}
p$A1_L2_48_1F
p$A1_L1_48_1F
p$A1_L1_106_1F
p$A1_FL_106_1F

p$A1_S1_106_1F
p$A1_S2_106_1F
p$A1_S3_106_1F
p$A1_S4_106_1F
p$A1_S5_106_1F
p$A1_S6_106_1F
p$A1_S7_106_1F
p$A1_S8_106_1F
p$A1_S9_120_1F

p$A1_R_106_1F
p$A1_RS_106_1F
p$A1_BS_106_1F
p$A1_G_106_1F

dev.off()



### Relative abundance

Stat<-list()
p<-list()

bac.clean.ss.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018" &Replication != "G_0")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)

#transformation counts to rel. abund.
bac.ra = transform_sample_counts(bac.clean.ss.18, function(x) x/sum(x))
DT=data.table(otu_table(bac.ra), keep.rownames=T, key="rn")
tax <- tax_table(bac.clean.ss.18)
tax <- data.frame(tax, stringsAsFactors = F)
OTU_B=t(otu_table(bac.clean.ss.18))
TAXA_B=tax_table(bac.clean.ss.18)

meta.18.all<-sample_data(bac.clean.ss.18)
meta.18.all$SampleID <- rownames(meta.18.all)

L1.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L1")])
L2.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L2")])
L3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L3")])
FL.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "FL")])
S1.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S1")])
S2.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S2")])
S3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S3")])
S4.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S4")])
S5.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S5")])
S6.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S6")])
S7.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S7")])
S8.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S8")])
S9.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S9")])
R.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "R")])
RS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "RS")])
BS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "BS")])
G.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "G")])

All<-list(L1.series, L2.series, L3.series,FL.series, 
          S1.series, S2.series,S3.series, S4.series, S5.series, S6.series,S7.series, S8.series, S9.series,
          R.series, RS.series, BS.series,G.series)

DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
DT.m=merge(DT,DT.taxa)
ColTaxa<-colnames(DT.taxa)

for (i in All) {
  
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn
  
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_B[i,LIST],taxon=TAXA_B[LIST,])
  
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  
  ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
  DT.sncm$ColPhylum <- ifelse(!DT.sncm$Phylum %in% ColPhylum, "Other", ifelse(DT.sncm$Phylum == "Bacteroidetes", "Bacteroidetes", ifelse(DT.sncm$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.sncm$Phylum == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  
  Richness<-sncm.out$fitstats$Richness
  
  fit.stat<-sncm.out$fitstats
  
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  
  Stat[[i[[1]]]]<-fit.stat
  print("######################################################################")
}
sncm.stat<-rbind(Stat$A1_L2_48_1B, Stat$A1_L1_48_1B, Stat$A1_L1_106_1B, Stat$A1_FL_106_1B, Stat$A1_S1_106_1B, 
                 Stat$A1_S2_106_1B, Stat$A1_S3_106_1B, Stat$A1_S4_106_1B, Stat$A1_S5_106_1B, Stat$A1_S6_106_1B, 
                 Stat$A1_S7_106_1B, Stat$A1_S8_106_1B, Stat$A1_S9_120_1B, Stat$A1_R_106_1B, Stat$A1_RS_106_1B, 
                 Stat$A1_BS_106_1B, Stat$A1_G_106_1B)
rownames(sncm.stat)<-c("L1", "L2", "L3", "FL", "S1", "S2", "S3", "S4", "S5", "S6","S7","S8","S9","R","RS","BS","G")
sncm.stat$comp<-c("L","L","L","L","S","S","S","S","S","S","S","S","S","R","RS","BS", "G")

write.table(sncm.stat, "sncm.stat_all replicates_relative abundance_bacteria_final_211017.tsv", sep='\t', quote=F)

mid<-mean(sncm.stat$Richness)

pp<-ggplot(sncm.stat, aes(x=Rsqr, y=m, color=Richness, shape=comp))

panel_a=pp+geom_point(size=5)+
  scale_color_gradient2(low="darkcyan", mid="darkorange", high="darkviolet", midpoint=mid)+
  theme_bw() + theme_change + theme(aspect.ratio = 1)
scale_shape_manual(values=Palette.shape.2)

#	pdf("Fig3a.pdf", paper="A4", useDingbats=FALSE)
print(panel_a)
dev.off()


##Grouped by over- and down-represented and fitted OTUs
for (i in All) {
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn	
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_B[i,LIST],taxon=TAXA_B[LIST,])
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Verrucomicrobia","Alphaproteobacteria","Deltaproteobacteria","Gammaproteobacteria","Chloroflexi","Acidobacteria")
  DT.sncm$ColPhylum <- 0
  DT.sncm$ColPhylum[-which(DT.sncm$Phylum %in% ColPhylum)] <- "Other"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Alphaproteobacteria")] <- "Alphaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Gammaproteobacteria")] <- "Gammaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Class =="Deltaproteobacteria")] <- "Deltaproteobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum =="Bacteroidetes")] <- "Bacteroidetes"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Actinobacteria")] <- "Actinobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Verrucomicrobia")] <- "Verrucomicrobia"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Firmicutes")] <- "Firmicutes"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Acidobacteria")] <- "Acidobacteria"
  DT.sncm$ColPhylum[which(DT.sncm$Phylum == "Chloroflexi")] <- "Chloroflexi" 
  
  DT.sncm$ColPhylum <- factor(DT.sncm$ColPhylum, levels = c("Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",
                                                            "Chloroflexi","Firmicutes","Acidobacteria", "Actinobacteria","Bacteroidetes","Verrucomicrobia","Other"))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  DT.sncm$fit_class_inherit <- ifelse(DT.sncm$rn %in% persistent.bac, "Inherit", "Others")
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  Richness<-sncm.out$fitstats$Richness
  fit.stat<-sncm.out$fitstats
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  Stat[[i[[1]]]]<-fit.stat
  DT.sncm$fit_class <- factor(DT.sncm$fit_class, levels = c("Above prediction", "As predicted", "Below prediction"))
  names(DT.sncm)[1] <- "OTU"
  DT.sncm <- merge(DT.sncm, OTU_id.list, by = "OTU")
  write.table(DT.sncm, paste0("Bacteria_DT.sncm_all reps_211017_relabund_",i[[1]],".tsv"), sep = '\t', quote = F)
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  p1= pp1 + geom_jitter (aes(color=fit_class, shape = fit_class_inherit),size=2.5) +
    theme_bw() + theme_change + theme(aspect.ratio = 1)+
    scale_color_manual(values=c("#cc6666", "#999999", "#6699cc"))+
    scale_shape_manual(values=Palette.shape.3)+
    geom_line(aes(x=log(p), y=freq.pred), color = '#333333', size = 0.75)+
    geom_line(aes(x=log(p), y=pred.lwr),linetype="longdash",color = '#333333', size = 0.5)+
    geom_line(aes(x=log(p), y=pred.upr),linetype="longdash",color = '#333333', size = 0.5)+
    theme(axis.text.x = element_text(vjust=0.4,size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))
  p[[i[[1]]]]<-p1
  print("######################################################################")
}
p$A1_L2_48_1B
p$A1_L1_48_1B
p$A1_L1_106_1B
p$A1_FL_106_1B

p$A1_S1_106_1B
p$A1_S2_106_1B
p$A1_S3_106_1B
p$A1_S4_106_1B
p$A1_S5_106_1B
p$A1_S6_106_1B
p$A1_S7_106_1B
p$A1_S8_106_1B
p$A1_S9_120_1B

p$A1_R_106_1B
p$A1_RS_106_1B
p$A1_BS_106_1B
p$A1_G_106_1B

dev.off()


###Fungi
Palette.class <-c("Mortierellomycetes" = "#4E734E", "Eurotiomycetes"= "#6DA9DC",
                  "Mucoromycetes"="#E4AF2C","Tremellomycetes"="#BE4146", "Microbotryomycetes" ="#DC9A9E",
                  "Dothideomycetes" = "#5195D1", "Agaricomycetes" = "#CC6C71", "Blastocladiomycetes" = "#87AC88",
                  "Sordariomycetes"= "#1E63AF","Leotiomycetes"= "#11335F", "Cystobasidiomycetes" = "#A871AE", "Pezizomycetes" = "#C0DBF3",
                  "unidentified" ="#000000", "Endogonomycetes" ="#CC9900", "Pucciniomycetes" = "#0099FF", "Ustilaginomycetes"="#ffcc33","Arthoniomycetes" = "#cccc99",
                  "Saccharomycetes" = "#666699", "Other" = "light grey", "Malasseziomycetes" = "#99cc99", "Exobasidiomycetes" = "#663366", 
                  "Rhizophydiomycetes"="#ffccff","Wallemiomycetes"="#333300","Agaricostilbomycetes"="#003333","Orbiliomycetes"="#cc9966","Chytridiomycetes" = "#006600",
                  "Spizellomycetes" = "#66cc99","Moniliellomycetes" = "#0099cc","Monoblepharidomycetes" = "#CC0000", "Dacrymycetes" = "#669999",
                  "Spiculogloeomycetes" = "#99cc00")


Stat<-list()
p<-list()

fun.clean.ss.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018" &Replication != "G_0")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)

#transformation counts to rel. abund.
fun.ra = transform_sample_counts(fun.clean.ss.18, function(x) x/sum(x))
DT=data.table(otu_table(fun.ra), keep.rownames=T, key="rn")
tax <- tax_table(fun.clean.ss.18)
tax <- data.frame(tax, stringsAsFactors = F)
OTU_F=t(otu_table(fun.clean.ss.18))
TAXA_F=tax_table(fun.clean.ss.18)

meta.18.all<-sample_data(fun.clean.ss.18)
meta.18.all$SampleID <- rownames(meta.18.all)

L1.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L1")])
L2.series <-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L2")])
L3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "L3")])
FL.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "FL")])
S1.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S1")])
S2.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S2")])
S3.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S3")])
S4.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S4")])
S5.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S5")])
S6.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S6")])
S7.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S7")])
S8.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S8")])
S9.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "S9")])
R.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "R")])
RS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "RS")])
BS.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "BS")])
G.series<-as.character(meta.18.all$SampleID[which(meta.18.all$Microhabitat == "G")])

All<-list(L1.series, L2.series, L3.series,FL.series, 
          S1.series, S2.series,S3.series, S4.series, S5.series, S6.series,S7.series, S8.series, S9.series,
          R.series, RS.series, BS.series,G.series)

DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
DT.m=merge(DT,DT.taxa)

ColTaxa<-colnames(DT.taxa)

for (i in All) {
  
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn
  
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_F[i,LIST],taxon=TAXA_F[LIST,])
  
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  
  ColClass<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
  DT.sncm$ColClass <- ifelse(!DT.sncm$Class %in% ColClass, "Other", ifelse(DT.sncm$Class == "Bacteroidetes", "Bacteroidetes", ifelse(DT.sncm$Class == "Firmicutes", "Firmicutes", ifelse (DT.sncm$Class == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  
  Richness<-sncm.out$fitstats$Richness
  
  fit.stat<-sncm.out$fitstats
  
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  
  Stat[[i[[1]]]]<-fit.stat
  print("######################################################################")
}

sncm.stat.f<-rbind(Stat$A1_L2_48_1F, Stat$A1_L1_48_1F, Stat$A1_L1_106_1F, Stat$A1_FL_106_1F, Stat$A1_S1_106_1F, 
                   Stat$A1_S2_106_1F, Stat$A1_S3_106_1F, Stat$A1_S4_106_1F, Stat$A1_S5_106_1F, Stat$A1_S6_106_1F, 
                   Stat$A1_S7_106_1F, Stat$A1_S8_106_1F, Stat$A1_S9_120_1F, Stat$A1_R_106_1F, Stat$A1_RS_106_1F, 
                   Stat$A1_BS_106_1F, Stat$A1_G_106_1F)
rownames(sncm.stat.f)<-c("L1", "L2", "L3", "FL", "S1", "S2", "S3", "S4", "S5", "S6","S7","S8","S9","R","RS","BS","G")
sncm.stat.f$comp<-c("L","L","L","L","S","S","S","S","S","S","S","S","S","R","RS","BS", "G")

write.table(sncm.stat.f, "sncm.stat_all replicates_relative abundance_fungi_final.tsv", sep='\t', quote=F)

mid<-mean(sncm.stat.f$Richness)

pp<-ggplot(sncm.stat.f, aes(x=Rsqr, y=m, color=Richness, shape=comp))

panel_a=pp+geom_point(size=5)+
  scale_color_gradient2(low="darkcyan", mid="darkorange", high="darkviolet", midpoint=mid)+
  theme_bw() + theme_change + theme(aspect.ratio = 1)
scale_shape_manual(values=Palette.shape.2)

#	pdf("Fig3a.pdf", paper="A4", useDingbats=FALSE)
print(panel_a)
dev.off()

Stat<-list()
p<-list()

##Grouped by over- and down-represented and fitted OTUs
for (i in All) {
  #Fiting Sloan Neutral Model
  print(i)
  DT.i<-DT.m[,c("rn",..i)]
  DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
  LIST=DT.i[DT.i$mean!=0]$rn	
  print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
  sncm.out=fit_sncm(spp=OTU_F[i,LIST],taxon=TAXA_F[LIST,])
  DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
  ColClass<-c("Tremellomycetes","Sordariomycetes","Leotiomycetes","Dothideomycetes","Agaricomycetes","Mortierellomycetes","Ustilaginomycetes","Eurotiomycetes","Saccharomycetes")
  DT.sncm$ColClass <- 0
  DT.sncm$ColClass[-which(DT.sncm$Class %in% ColClass)] <- "Other"
  DT.sncm$ColClass[which(DT.sncm$Class =="Tremellomycetes")] <- "Tremellomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Sordariomycetes")] <- "Sordariomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Leotiomycetes")] <- "Leotiomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class =="Dothideomycetes")] <- "Dothideomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Agaricomycetes")] <- "Agaricomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Mortierellomycetes")] <- "Mortierellomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Ustilaginomycetes")] <- "Ustilaginomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Eurotiomycetes")] <- "Eurotiomycetes"
  DT.sncm$ColClass[which(DT.sncm$Class == "Saccharomycetes")] <- "Saccharomycetes" 
  
  DT.sncm$ColClass <- factor(DT.sncm$ColClass, levels = c("Sordariomycetes", "Dothideomycetes", "Leotiomycetes",
                                                          "Eurotiomycetes","Saccharomycetes","Tremellomycetes", "Agaricomycetes","Ustilaginomycetes","Mortierellomycetes","Other"))
  
  DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))
  DT.sncm$fit_class_inherit <- ifelse(DT.sncm$rn %in% persistent.fun, "Inherit", "Others")
  Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
  Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100
  Richness<-sncm.out$fitstats$Richness
  fit.stat<-sncm.out$fitstats
  print(paste0("Percentage OTU Below = ", Blw/Richness))	
  print(paste0("Percentage OTU Above = ", Abv/Richness))
  Stat[[i[[1]]]]<-fit.stat
  DT.sncm$fit_class <- factor(DT.sncm$fit_class, levels = c("Above prediction", "As predicted", "Below prediction"))
  names(DT.sncm)[1] <- "OTU"
  DT.sncm <- merge(DT.sncm, OTU_id.list, by = "OTU")
  write.table(DT.sncm, paste0("Fungi_DT.sncm_all reps_211017_relabund_",i[[1]],".tsv"), sep = '\t', quote = F)
  
  pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))
  p1= pp1 + geom_jitter (aes(color=fit_class, shape = fit_class_inherit),size=2.5) +
    theme_bw() + theme_change + theme(aspect.ratio = 1)+
    scale_color_manual(values=c("#cc6666", "#999999", "#6699cc"))+
    scale_shape_manual(values=Palette.shape.3)+
    geom_line(aes(x=log(p), y=freq.pred), color = '#333333', size = 0.75)+
    geom_line(aes(x=log(p), y=pred.lwr),linetype="longdash",color = '#333333', size = 0.5)+
    geom_line(aes(x=log(p), y=pred.upr),linetype="longdash",color = '#333333', size = 0.5)+
    theme(axis.text.x = element_text(vjust=0.4,size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))
  p[[i[[1]]]]<-p1
  print("######################################################################")
}
p$A1_L2_48_1F
p$A1_L1_48_1F
p$A1_L1_106_1F
p$A1_FL_106_1F

p$A1_S1_106_1F
p$A1_S2_106_1F
p$A1_S3_106_1F
p$A1_S4_106_1F
p$A1_S5_106_1F
p$A1_S6_106_1F
p$A1_S7_106_1F
p$A1_S8_106_1F
p$A1_S9_120_1F

p$A1_R_106_1F
p$A1_RS_106_1F
p$A1_BS_106_1F
p$A1_G_106_1F

dev.off()
