### Number and proportion of shared OTUs between seeds and other compartments
bac.clean.ss.f
fun.clean.ss.f

bac.clean.ss.rel <- transform(bac.clean.ss.f, transform = "compositional")
bac.clean.ss.18 <- subset_samples(bac.clean.ss.rel, Year == "year_2018")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)

melt.bac.clean.rel <-psmelt(bac.clean.ss.18)



fun.clean.ss.rel <- transform(fun.clean.ss.f, transform = "compositional")
fun.clean.ss.18 <- subset_samples(fun.clean.ss.rel, Year == "year_2018")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)

melt.fun.clean.rel <-psmelt(fun.clean.ss.18)

### Number of shared OTUs with G141 and rice compartments
## Bacteria
otu.G.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0)])
otu.L.48 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_48","L2_48") & melt.bac.clean.rel$Abundance > 0)])
otu.S.48 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_48","S2_48","S3_48") & melt.bac.clean.rel$Abundance > 0)])
otu.R.48 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_48") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.48 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_48") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.48 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_48") & melt.bac.clean.rel$Abundance > 0)])

S.48 <-intersect(otu.S.48,otu.G.141)
L.48 <-intersect(otu.L.48,otu.G.141)
R.48 <-intersect(otu.R.48,otu.G.141)
RS.48 <-intersect(otu.RS.48,otu.G.141)
BS.48 <-intersect(otu.BS.48,otu.G.141)


Reduce(intersect, list(otu.L.48,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.48,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.48,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.48,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.48,otu.G.141, persistent.bac))


otu.L.62 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_62","L2_62","L3_62") & melt.bac.clean.rel$Abundance > 0)])
otu.S.62 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_62","S2_62","S3_62","S4_62") & melt.bac.clean.rel$Abundance > 0)])
otu.R.62 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_62") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.62 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_62") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.62 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_62") & melt.bac.clean.rel$Abundance > 0)])

S.62 <-intersect(otu.S.62,otu.G.141)
L.62 <-intersect(otu.L.62,otu.G.141)
R.62 <-intersect(otu.R.62,otu.G.141)
RS.62 <-intersect(otu.RS.62,otu.G.141)
BS.62 <-intersect(otu.BS.62,otu.G.141)

Reduce(intersect, list(otu.L.62,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.62,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.62,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.62,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.62,otu.G.141, persistent.bac))

otu.L.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_76","L2_76","L3_76","FL_76") & melt.bac.clean.rel$Abundance > 0)])
otu.S.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_76","S2_76","S3_76","S4_76","S5_76","S6_76") & melt.bac.clean.rel$Abundance > 0)])
otu.R.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_76") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_76") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_76") & melt.bac.clean.rel$Abundance > 0)])
otu.G.76 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("G_76") & melt.bac.clean.rel$Abundance > 0)])

S.76 <-intersect(otu.S.76,otu.G.141)
L.76 <-intersect(otu.L.76,otu.G.141)
R.76 <-intersect(otu.R.76,otu.G.141)
RS.76 <-intersect(otu.RS.76,otu.G.141)
BS.76 <-intersect(otu.BS.76,otu.G.141)

S.76.2 <-intersect(otu.S.76,otu.G.76)
L.76.2 <-intersect(otu.L.76,otu.G.76)

intersect(S.76.2,L.76.2)

Reduce(intersect, list(otu.L.76,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.76,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.76,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.76,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.76,otu.G.141, persistent.bac))

otu.L.90 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_90","L2_90","L3_90","FL_90") & melt.bac.clean.rel$Abundance > 0)])
otu.S.90 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_90","S2_90","S3_90","S4_90","S5_90","S6_90","S7_90","S8_90","S9_90") & melt.bac.clean.rel$Abundance > 0)])
otu.R.90 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_90") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.90 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_90") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.90 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_90") & melt.bac.clean.rel$Abundance > 0)])

S.90 <-intersect(otu.S.90,otu.G.141)
L.90 <-intersect(otu.L.90,otu.G.141)
R.90 <-intersect(otu.R.90,otu.G.141)
RS.90 <-intersect(otu.RS.90,otu.G.141)
BS.90 <-intersect(otu.BS.90,otu.G.141)

Reduce(intersect, list(otu.L.90,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.90,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.90,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.90,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.90,otu.G.141, persistent.bac))

otu.L.106 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_106","L2_106","L3_106","FL_106") & melt.bac.clean.rel$Abundance > 0)])
otu.S.106 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_106","S2_106","S3_106","S4_106","S5_106","S6_106","S7_106","S8_106","S9_106") & melt.bac.clean.rel$Abundance > 0)])
otu.R.106 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_106") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.106 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_106") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.106 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_106") & melt.bac.clean.rel$Abundance > 0)])

S.106 <-intersect(otu.S.106,otu.G.141)
L.106 <-intersect(otu.L.106,otu.G.141)
R.106 <-intersect(otu.R.106,otu.G.141)
RS.106 <-intersect(otu.RS.106,otu.G.141)
BS.106 <-intersect(otu.BS.106,otu.G.141)

Reduce(intersect, list(otu.L.106,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.106,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.106,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.106,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.106,otu.G.141, persistent.bac))

otu.L.120 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_120","L2_120","L3_120","FL_120") & melt.bac.clean.rel$Abundance > 0)])
otu.S.120 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_120","S2_120","S3_120","S4_120","S5_120","S6_120","S7_120","S8_120","S9_120") & melt.bac.clean.rel$Abundance > 0)])
otu.R.120 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_120") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.120 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_120") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.120 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_120") & melt.bac.clean.rel$Abundance > 0)])

S.120 <-intersect(otu.S.120,otu.G.141)
L.120 <-intersect(otu.L.120,otu.G.141)
R.120 <-intersect(otu.R.120,otu.G.141)
RS.120 <-intersect(otu.RS.120,otu.G.141)
BS.120 <-intersect(otu.BS.120,otu.G.141)

Reduce(intersect, list(otu.L.120,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.120,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.120,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.120,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.120,otu.G.141, persistent.bac))

otu.L.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("L1_141","L2_141","L3_141","FL_141") & melt.bac.clean.rel$Abundance > 0)])
otu.S.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("S1_141","S2_141","S3_141","S4_141","S5_141","S6_141","S7_141","S8_141","S9_141") & melt.bac.clean.rel$Abundance > 0)])
otu.R.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("R_141") & melt.bac.clean.rel$Abundance > 0)])
otu.RS.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("RS_141") & melt.bac.clean.rel$Abundance > 0)])
otu.BS.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication %in% c("BS_141") & melt.bac.clean.rel$Abundance > 0)])

S.141 <-intersect(otu.S.141,otu.G.141)
L.141 <-intersect(otu.L.141,otu.G.141)
R.141 <-intersect(otu.R.141,otu.G.141)
RS.141 <-intersect(otu.RS.141,otu.G.141)
BS.141 <-intersect(otu.BS.141,otu.G.141)

Reduce(intersect, list(otu.L.141,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.S.141,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.R.141,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.RS.141,otu.G.141, persistent.bac))
Reduce(intersect, list(otu.BS.141,otu.G.141, persistent.bac))


### Fungi
otu.G.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0)])
otu.L.48 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_48","L2_48") & melt.fun.clean.rel$Abundance > 0)])
otu.S.48 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_48","S2_48","S3_48") & melt.fun.clean.rel$Abundance > 0)])
otu.R.48 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_48") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.48 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_48") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.48 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_48") & melt.fun.clean.rel$Abundance > 0)])

S.48 <-intersect(otu.S.48,otu.G.141)
L.48 <-intersect(otu.L.48,otu.G.141)
R.48 <-intersect(otu.R.48,otu.G.141)
RS.48 <-intersect(otu.RS.48,otu.G.141)
BS.48 <-intersect(otu.BS.48,otu.G.141)


Reduce(intersect, list(otu.L.48,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.48,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.48,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.48,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.48,otu.G.141, persistent.fun))


otu.L.62 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_62","L2_62","L3_62") & melt.fun.clean.rel$Abundance > 0)])
otu.S.62 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_62","S2_62","S3_62","S4_62") & melt.fun.clean.rel$Abundance > 0)])
otu.R.62 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_62") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.62 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_62") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.62 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_62") & melt.fun.clean.rel$Abundance > 0)])

S.62 <-intersect(otu.S.62,otu.G.141)
L.62 <-intersect(otu.L.62,otu.G.141)
R.62 <-intersect(otu.R.62,otu.G.141)
RS.62 <-intersect(otu.RS.62,otu.G.141)
BS.62 <-intersect(otu.BS.62,otu.G.141)

Reduce(intersect, list(otu.L.62,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.62,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.62,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.62,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.62,otu.G.141, persistent.fun))

otu.L.76 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_76","L2_76","L3_76","FL_76") & melt.fun.clean.rel$Abundance > 0)])
otu.S.76 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_76","S2_76","S3_76","S4_76","S5_76","S6_76") & melt.fun.clean.rel$Abundance > 0)])
otu.R.76 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_76") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.76 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_76") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.76 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_76") & melt.fun.clean.rel$Abundance > 0)])

S.76 <-intersect(otu.S.76,otu.G.141)
L.76 <-intersect(otu.L.76,otu.G.141)
R.76 <-intersect(otu.R.76,otu.G.141)
RS.76 <-intersect(otu.RS.76,otu.G.141)
BS.76 <-intersect(otu.BS.76,otu.G.141)

Reduce(intersect, list(otu.L.76,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.76,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.76,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.76,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.76,otu.G.141, persistent.fun))

otu.L.90 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_90","L2_90","L3_90","FL_90") & melt.fun.clean.rel$Abundance > 0)])
otu.S.90 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_90","S2_90","S3_90","S4_90","S5_90","S6_90","S7_90","S8_90","S9_90") & melt.fun.clean.rel$Abundance > 0)])
otu.R.90 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_90") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.90 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_90") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.90 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_90") & melt.fun.clean.rel$Abundance > 0)])

S.90 <-intersect(otu.S.90,otu.G.141)
L.90 <-intersect(otu.L.90,otu.G.141)
R.90 <-intersect(otu.R.90,otu.G.141)
RS.90 <-intersect(otu.RS.90,otu.G.141)
BS.90 <-intersect(otu.BS.90,otu.G.141)

Reduce(intersect, list(otu.L.90,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.90,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.90,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.90,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.90,otu.G.141, persistent.fun))

otu.L.106 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_106","L2_106","L3_106","FL_106") & melt.fun.clean.rel$Abundance > 0)])
otu.S.106 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_106","S2_106","S3_106","S4_106","S5_106","S6_106","S7_106","S8_106","S9_106") & melt.fun.clean.rel$Abundance > 0)])
otu.R.106 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_106") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.106 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_106") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.106 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_106") & melt.fun.clean.rel$Abundance > 0)])

S.106 <-intersect(otu.S.106,otu.G.141)
L.106 <-intersect(otu.L.106,otu.G.141)
R.106 <-intersect(otu.R.106,otu.G.141)
RS.106 <-intersect(otu.RS.106,otu.G.141)
BS.106 <-intersect(otu.BS.106,otu.G.141)

Reduce(intersect, list(otu.L.106,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.106,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.106,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.106,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.106,otu.G.141, persistent.fun))

otu.L.120 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_120","L2_120","L3_120","FL_120") & melt.fun.clean.rel$Abundance > 0)])
otu.S.120 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_120","S2_120","S3_120","S4_120","S5_120","S6_120","S7_120","S8_120","S9_120") & melt.fun.clean.rel$Abundance > 0)])
otu.R.120 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_120") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.120 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_120") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.120 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_120") & melt.fun.clean.rel$Abundance > 0)])

S.120 <-intersect(otu.S.120,otu.G.141)
L.120 <-intersect(otu.L.120,otu.G.141)
R.120 <-intersect(otu.R.120,otu.G.141)
RS.120 <-intersect(otu.RS.120,otu.G.141)
BS.120 <-intersect(otu.BS.120,otu.G.141)

Reduce(intersect, list(otu.L.120,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.120,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.120,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.120,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.120,otu.G.141, persistent.fun))

otu.L.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("L1_141","L2_141","L3_141","FL_141") & melt.fun.clean.rel$Abundance > 0)])
otu.S.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("S1_141","S2_141","S3_141","S4_141","S5_141","S6_141","S7_141","S8_141","S9_141") & melt.fun.clean.rel$Abundance > 0)])
otu.R.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("R_141") & melt.fun.clean.rel$Abundance > 0)])
otu.RS.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("RS_141") & melt.fun.clean.rel$Abundance > 0)])
otu.BS.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication %in% c("BS_141") & melt.fun.clean.rel$Abundance > 0)])

S.141 <-intersect(otu.S.141,otu.G.141)
L.141 <-intersect(otu.L.141,otu.G.141)
R.141 <-intersect(otu.R.141,otu.G.141)
RS.141 <-intersect(otu.RS.141,otu.G.141)
BS.141 <-intersect(otu.BS.141,otu.G.141)

Reduce(intersect, list(otu.L.141,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.S.141,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.R.141,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.RS.141,otu.G.141, persistent.fun))
Reduce(intersect, list(otu.BS.141,otu.G.141, persistent.fun))


head(melt.bac.clean.rel)
head(melt.fun.clean.rel)
##Cumulative abundance of shared OTUs in plant compartments
melt.bac.clean.rel.mean <- melt.bac.clean.rel %>% group_by(OTU, Replication2) %>% summarise(MeanAbun = mean(Abundance))

melt.fun.clean.rel$Compartment2 <- 0
melt.fun.clean.rel$Compartment2[which(melt.fun.clean.rel$Compartment == "Leaf")] <- "L"
melt.fun.clean.rel$Compartment2[which(melt.fun.clean.rel$Compartment == "Stem")] <- "S"
melt.fun.clean.rel$Compartment2[which(melt.fun.clean.rel$Compartment == "Root")] <- "R"
melt.fun.clean.rel$Compartment2[which(melt.fun.clean.rel$Compartment == "Rhizosphere")] <- "RS"
melt.fun.clean.rel$Compartment2[which(melt.fun.clean.rel$Compartment == "Bulk_soil")] <- "BS"

melt.fun.clean.rel$Replication2 <- paste0(melt.fun.clean.rel$Compartment2,"_",melt.fun.clean.rel$Days)
melt.fun.clean.rel.mean <- melt.fun.clean.rel %>% group_by(OTU, Replication2) %>% summarise(MeanAbun = mean(Abundance))


esti.mean.abun<-function(mean.abun, common.OTU.list, Keyw){
  mean.abun.tab <-subset(mean.abun, OTU %in% common.OTU.list & Replication2 == Keyw)
  sum.abun<-sum(mean.abun.tab$MeanAbun)
  return(print(sum.abun))
}


esti.mean.abun(melt.bac.clean.rel.mean, L.48, "L_48")
esti.mean.abun(melt.bac.clean.rel.mean, S.48, "S_48")
esti.mean.abun(melt.bac.clean.rel.mean, R.48, "R_48")
esti.mean.abun(melt.bac.clean.rel.mean, RS.48, "RS_48")
esti.mean.abun(melt.bac.clean.rel.mean, BS.48, "BS_48")


esti.mean.abun(melt.bac.clean.rel.mean, L.62, "L_62")
esti.mean.abun(melt.bac.clean.rel.mean, S.62, "S_62")
esti.mean.abun(melt.bac.clean.rel.mean, R.62, "R_62")
esti.mean.abun(melt.bac.clean.rel.mean, RS.62, "RS_62")
esti.mean.abun(melt.bac.clean.rel.mean, BS.62, "BS_62")

esti.mean.abun(melt.bac.clean.rel.mean, L.76, "L_76")
esti.mean.abun(melt.bac.clean.rel.mean, S.76, "S_76")
esti.mean.abun(melt.bac.clean.rel.mean, R.76, "R_76")
esti.mean.abun(melt.bac.clean.rel.mean, RS.76, "RS_76")
esti.mean.abun(melt.bac.clean.rel.mean, BS.76, "BS_76")

esti.mean.abun(melt.bac.clean.rel.mean, L.90, "L_90")
esti.mean.abun(melt.bac.clean.rel.mean, S.90, "S_90")
esti.mean.abun(melt.bac.clean.rel.mean, R.90, "R_90")
esti.mean.abun(melt.bac.clean.rel.mean, RS.90, "RS_90")
esti.mean.abun(melt.bac.clean.rel.mean, BS.90, "BS_90")

esti.mean.abun(melt.bac.clean.rel.mean, L.106, "L_106")
esti.mean.abun(melt.bac.clean.rel.mean, S.106, "S_106")
esti.mean.abun(melt.bac.clean.rel.mean, R.106, "R_106")
esti.mean.abun(melt.bac.clean.rel.mean, RS.106, "RS_106")
esti.mean.abun(melt.bac.clean.rel.mean, BS.106, "BS_106")

esti.mean.abun(melt.bac.clean.rel.mean, L.120, "L_120")
esti.mean.abun(melt.bac.clean.rel.mean, S.120, "S_120")
esti.mean.abun(melt.bac.clean.rel.mean, R.120, "R_120")
esti.mean.abun(melt.bac.clean.rel.mean, RS.120, "RS_120")
esti.mean.abun(melt.bac.clean.rel.mean, BS.120, "BS_120")


esti.mean.abun(melt.bac.clean.rel.mean, L.141, "L_141")
esti.mean.abun(melt.bac.clean.rel.mean, S.141, "S_141")
esti.mean.abun(melt.bac.clean.rel.mean, R.141, "R_141")
esti.mean.abun(melt.bac.clean.rel.mean, RS.141, "RS_141")
esti.mean.abun(melt.bac.clean.rel.mean, BS.141, "BS_141")


esti.mean.abun(melt.fun.clean.rel.mean, L.48, "Leaf_48")
esti.mean.abun(melt.fun.clean.rel.mean, S.48, "Stem_48")
esti.mean.abun(melt.fun.clean.rel.mean, R.48, "Root_48")
esti.mean.abun(melt.fun.clean.rel.mean, RS.48, "Rhizosphere_48")
esti.mean.abun(melt.fun.clean.rel.mean, BS.48, "Bulk_soil_48")


esti.mean.abun(melt.fun.clean.rel.mean, L.62, "Leaf_62")
esti.mean.abun(melt.fun.clean.rel.mean, S.62, "Stem_62")
esti.mean.abun(melt.fun.clean.rel.mean, R.62, "Root_62")
esti.mean.abun(melt.fun.clean.rel.mean, RS.62, "Rhizosphere_62")
esti.mean.abun(melt.fun.clean.rel.mean, BS.62, "Bulk_soil_62")

esti.mean.abun(melt.fun.clean.rel.mean, L.76, "Leaf_76")
esti.mean.abun(melt.fun.clean.rel.mean, S.76, "Stem_76")
esti.mean.abun(melt.fun.clean.rel.mean, R.76, "Root_76")
esti.mean.abun(melt.fun.clean.rel.mean, RS.76, "Rhizosphere_76")
esti.mean.abun(melt.fun.clean.rel.mean, BS.76, "Bulk_soil_76")

esti.mean.abun(melt.fun.clean.rel.mean, L.90, "Leaf_90")
esti.mean.abun(melt.fun.clean.rel.mean, S.90, "Stem_90")
esti.mean.abun(melt.fun.clean.rel.mean, R.90, "Root_90")
esti.mean.abun(melt.fun.clean.rel.mean, RS.90, "Rhizosphere_90")
esti.mean.abun(melt.fun.clean.rel.mean, BS.90, "Bulk_soil_90")

esti.mean.abun(melt.fun.clean.rel.mean, L.106, "Leaf_106")
esti.mean.abun(melt.fun.clean.rel.mean, S.106, "Stem_106")
esti.mean.abun(melt.fun.clean.rel.mean, R.106, "Root_106")
esti.mean.abun(melt.fun.clean.rel.mean, RS.106, "Rhizosphere_106")
esti.mean.abun(melt.fun.clean.rel.mean, BS.106, "Bulk_soil_106")

esti.mean.abun(melt.fun.clean.rel.mean, L.120, "Leaf_120")
esti.mean.abun(melt.fun.clean.rel.mean, S.120, "Stem_120")
esti.mean.abun(melt.fun.clean.rel.mean, R.120, "Root_120")
esti.mean.abun(melt.fun.clean.rel.mean, RS.120, "Rhizosphere_120")
esti.mean.abun(melt.fun.clean.rel.mean, BS.120, "Bulk_soil_120")


esti.mean.abun(melt.fun.clean.rel.mean, L.141, "Leaf_141")
esti.mean.abun(melt.fun.clean.rel.mean, S.141, "Stem_141")
esti.mean.abun(melt.fun.clean.rel.mean, R.141, "Root_141")
esti.mean.abun(melt.fun.clean.rel.mean, RS.141, "Rhizosphere_141")
esti.mean.abun(melt.fun.clean.rel.mean, BS.141, "Bulk_soil_141")


##Visualization
tab.summary <- read.xlsx('Number of OTUs shared with G141_scaled_3.xlsx',1)
tab.summary <- tab.summary[c(1:5)]
tab.summary$Days <- as.factor(as.character(tab.summary$Days))
tab.summary$Days <-factor(tab.summary$Days, levels = c("48","62","76","90","106","120","141"))
#tab.summary$Days <-factor(tab.summary$Days, levels = c("141","120","106","90","76","62","48"))
tab.summary$Compartment <-factor(tab.summary$Compartment, levels = c("BS","RS","R","S","L"))

melt.tab.summary <- melt(tab.summary)
melt.tab.summary$value <- as.numeric(as.character(melt.tab.summary$value))

g<-ggplot(subset(melt.tab.summary,Kingdom == "Bacteria") , aes(x= Compartment,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 1)+theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
g
#g+scale_y_continuous(sec.axis = sec_axis(~./52, name="Percent\n"))+coord_flip()

g<-ggplot(subset(melt.tab.summary,Kingdom == "Fungi") , aes(x= Compartment,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 1)+theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
g
#g+scale_y_continuous(sec.axis = sec_axis(~./107, name="Percent\n"))+coord_flip()




### Section
##Bacteria
otu.G.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0)])

section.shared.tab <- data.frame(section = as.character(unique(melt.bac.clean.rel$Replication[which(!(melt.bac.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]), number = integer()))
section.shared.tab$shared.number <- 0
for (i in as.character(unique(melt.bac.clean.rel$Replication[which(!(melt.bac.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication== i & melt.bac.clean.rel$Abundance > 0)])
  shared.otu.name <-intersect(otu.name,otu.G.141)
  section.shared.tab$shared.number[which(section.shared.tab$section == i)]<-length(shared.otu.name)
}

section.shared.tab$persist.number <- 0
for (i in as.character(unique(melt.bac.clean.rel$Replication[which(!(melt.bac.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication== i & melt.bac.clean.rel$Abundance > 0)])
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141, persistent.bac))
  section.shared.tab$persist.number[which(section.shared.tab$section == i)]<-length(persist.otu.name)
}

##Fungi
otu.G.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0)])

section.shared.tab.f <- data.frame(section = as.character(unique(melt.fun.clean.rel$Replication[which(!(melt.fun.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]), number = integer()))
section.shared.tab.f$shared.number <- 0
for (i in as.character(unique(melt.fun.clean.rel$Replication[which(!(melt.fun.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication== i & melt.fun.clean.rel$Abundance > 0)])
  shared.otu.name <-intersect(otu.name,otu.G.141)
  section.shared.tab.f$shared.number[which(section.shared.tab.f$section == i)]<-length(shared.otu.name)
}

section.shared.tab.f$persist.number <- 0
for (i in as.character(unique(melt.fun.clean.rel$Replication[which(!(melt.fun.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication== i & melt.fun.clean.rel$Abundance > 0)])
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141, persistent.fun))
  section.shared.tab.f$persist.number[which(section.shared.tab.f$section == i)]<-length(persist.otu.name)
}


##Cumulative abundance of shared OTUs in plant compartments
melt.bac.clean.rel.mean <- melt.bac.clean.rel %>% group_by(OTU, Replication2) %>% summarise(MeanAbun = mean(Abundance))

melt.fun.clean.rel.mean <- melt.fun.clean.rel %>% group_by(OTU, Replication) %>% summarise(MeanAbun = mean(Abundance))

esti.mean.abun<-function(mean.abun, common.OTU.list, Keyw){
  mean.abun.tab <-subset(mean.abun, OTU %in% common.OTU.list & Replication == Keyw)
  sum.abun<-sum(mean.abun.tab$MeanAbun)
  return(sum.abun)
}



section.shared.tab$shared.abund <- 0
section.shared.tab$persist.abund <- 0
otu.G.141 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0)])

for (i in as.character(unique(melt.bac.clean.rel$Replication[which(!(melt.bac.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication== i & melt.bac.clean.rel$Abundance > 0)])
  shared.otu.name <-intersect(otu.name,otu.G.141)
  shared.abund<-esti.mean.abun(melt.bac.clean.rel.mean,shared.otu.name,i)
  section.shared.tab$shared.abund[which(section.shared.tab$section == i)]<-shared.abund
  
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141, persistent.bac))
  persist.abund<-esti.mean.abun(melt.bac.clean.rel.mean,persist.otu.name,i)
  section.shared.tab$persist.abund[which(section.shared.tab$section == i)]<-persist.abund
  
}


section.shared.tab.f$shared.abund <- 0
section.shared.tab.f$persist.abund <- 0
otu.G.141 <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0)])

for (i in as.character(unique(melt.fun.clean.rel$Replication[which(!(melt.fun.clean.rel$Replication %in% c("G_0","G_76","G_90","G_106","G_120","G_141")))]))){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication== i & melt.fun.clean.rel$Abundance > 0)])
  shared.otu.name <-intersect(otu.name,otu.G.141)
  shared.abund<-esti.mean.abun(melt.fun.clean.rel.mean,shared.otu.name,i)
  section.shared.tab.f$shared.abund[which(section.shared.tab.f$section == i)]<-shared.abund
  
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141, persistent.fun))
  persist.abund<-esti.mean.abun(melt.fun.clean.rel.mean,persist.otu.name,i)
  section.shared.tab.f$persist.abund[which(section.shared.tab.f$section == i)]<-persist.abund
  
}

section.shared.tab$Kingdom <- "Bacteria"
section.shared.tab.f$Kingdom <- "Fungi"

section.shared.tab.merged<-rbind(section.shared.tab,section.shared.tab.f)
write.xlsx(section.shared.tab.merged,"Shared OTUs at the section level.xlsx")


section.shared.tab.merged<- read.xlsx("Shared OTUs at the section level_scaled.xlsx",1)

head(section.shared.tab.merged)
###Visualization
section.shared.tab.merged$Days <- as.factor(as.character(section.shared.tab.merged$Days))
section.shared.tab.merged$Days <-factor(section.shared.tab.merged$Days, levels = c("48","62","76","90","106","120","141"))
#section.shared.tab.merged$Days <-factor(section.shared.tab.merged$Days, levels = c("141","120","106","90","76","62","48"))
section.shared.tab.merged$section <-factor(section.shared.tab.merged$section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))

melt.section.shared.tab.merged <- melt(section.shared.tab.merged)
melt.section.shared.tab.merged$value <- as.numeric(as.character(melt.section.shared.tab.merged$value))


g<-ggplot(subset(melt.section.shared.tab.merged,Kingdom == "Bacteria") , aes(x= section,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 3)+theme(aspect.ratio =0.5)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
g+scale_y_continuous(sec.axis = sec_axis(~./42, name="Percent\n"))+coord_flip()

g<-ggplot(subset(melt.tab.summary,Kingdom == "Fungi") , aes(x= Compartment,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 1)+theme(aspect.ratio =0.4)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
g+scale_y_continuous(sec.axis = sec_axis(~./107, name="Percent\n"))+coord_flip()






section.shared.tab.merged.sub <- section.shared.tab.merged[c(1,2,3,4,5,8)]
section.shared.tab.merged.sub$Days <- as.factor(as.character(section.shared.tab.merged.sub$Days))
section.shared.tab.merged.sub$Days <-factor(section.shared.tab.merged.sub$Days, levels = c("48","62","76","90","106","120","141"))
#section.shared.tab.merged.sub$Days <-factor(section.shared.tab.merged.sub$Days, levels = c("141","120","106","90","76","62","48"))
section.shared.tab.merged.sub$section <-factor(section.shared.tab.merged.sub$section, levels = c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL"))



melt.section.shared.tab.merged.sub <- melt(section.shared.tab.merged.sub)
melt.section.shared.tab.merged.sub$value <- as.numeric(as.character(melt.section.shared.tab.merged.sub$value))



g<-ggplot(subset(melt.section.shared.tab.merged.sub,Kingdom == "Bacteria") , aes(x= section,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 1)+theme(aspect.ratio =0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
#g+scale_y_continuous(sec.axis = sec_axis(~./42, name="Percent\n"))+coord_flip()
g
g<-ggplot(subset(melt.section.shared.tab.merged.sub,Kingdom == "Fungi") , aes(x= section,fill = variable,y= value)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  facet_wrap(~Days, ncol = 1)+theme(aspect.ratio =0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "right")#+ scale_color_manual(values = c("darkred", "steelblue"))
#g+scale_y_continuous(sec.axis = sec_axis(~./107, name="Percent\n"))+coord_flip()
g


#### Considering biological repeats
otu.G.141.A1 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1A1")])
otu.G.141.B2 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1B2")])
otu.G.141.C3 <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication == "G_141" & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1C3")])

##A1
list.A1<-unique(melt.bac.clean.rel$Replication2[which(!(melt.bac.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.bac.clean.rel$Plot == "UF1A1")])

comp.shared.tab.A1 <- data.frame(comp = as.character(as.factor(list.A1)))
comp.shared.tab.A1$shared.number <- 0
comp.shared.tab.A1$persist.number <- 0
comp.shared.tab.A1$shared.abund <- 0
comp.shared.tab.A1$persist.abund <- 0

for (i in as.character(list.A1)){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication2== i & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1A1")])
  shared.otu.name <-intersect(otu.name,otu.G.141.A1)
  comp.shared.tab.A1$shared.number[which(comp.shared.tab.A1$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.A1, persistent.bac))
  comp.shared.tab.A1$persist.number[which(comp.shared.tab.A1$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.bac.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.A1$shared.abund[which(comp.shared.tab.A1$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.bac.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.A1$persist.abund[which(comp.shared.tab.A1$comp == i)]<-persist.abund
  }
comp.shared.tab.A1$Plot <- "A1"


##B2
list.B2<-unique(melt.bac.clean.rel$Replication2[which(!(melt.bac.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.bac.clean.rel$Plot == "UF1B2")])

comp.shared.tab.B2 <- data.frame(comp = as.character(as.factor(list.B2)))
comp.shared.tab.B2$shared.number <- 0
comp.shared.tab.B2$persist.number <- 0
comp.shared.tab.B2$shared.abund <- 0
comp.shared.tab.B2$persist.abund <- 0

for (i in as.character(list.B2)){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication2== i & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1B2")])
  shared.otu.name <-intersect(otu.name,otu.G.141.B2)
  comp.shared.tab.B2$shared.number[which(comp.shared.tab.B2$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.B2, persistent.bac))
  comp.shared.tab.B2$persist.number[which(comp.shared.tab.B2$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.bac.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.B2$shared.abund[which(comp.shared.tab.B2$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.bac.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.B2$persist.abund[which(comp.shared.tab.B2$comp == i)]<-persist.abund
}
comp.shared.tab.B2$Plot <- "B2"


##C3
list.C3<-unique(melt.bac.clean.rel$Replication2[which(!(melt.bac.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.bac.clean.rel$Plot == "UF1C3")])

comp.shared.tab.C3 <- data.frame(comp = as.character(as.factor(list.C3)))
comp.shared.tab.C3$shared.number <- 0
comp.shared.tab.C3$persist.number <- 0
comp.shared.tab.C3$shared.abund <- 0
comp.shared.tab.C3$persist.abund <- 0

for (i in as.character(list.C3)){
  otu.name <- unique(melt.bac.clean.rel$OTU[which(melt.bac.clean.rel$Replication2== i & melt.bac.clean.rel$Abundance > 0 & melt.bac.clean.rel$Plot == "UF1C3")])
  shared.otu.name <-intersect(otu.name,otu.G.141.C3)
  comp.shared.tab.C3$shared.number[which(comp.shared.tab.C3$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.C3, persistent.bac))
  comp.shared.tab.C3$persist.number[which(comp.shared.tab.C3$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.bac.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.C3$shared.abund[which(comp.shared.tab.C3$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.bac.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.C3$persist.abund[which(comp.shared.tab.C3$comp == i)]<-persist.abund
}
comp.shared.tab.C3$Plot <- "C3"

comp.shared.tab.bac<-rbind(comp.shared.tab.A1, comp.shared.tab.B2, comp.shared.tab.C3)

comp.shared.tab.bac$Kingdom <- "Bacteria"



###Fungi
otu.G.141.A1.f <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1A1")])
otu.G.141.B2.f <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1B2")])
otu.G.141.C3.f <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication == "G_141" & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1C3")])

##A1
list.A1<-unique(melt.fun.clean.rel$Replication2[which(!(melt.fun.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.fun.clean.rel$Plot == "UF1A1")])

comp.shared.tab.A1.f <- data.frame(comp = as.character(as.factor(list.A1)))
comp.shared.tab.A1.f$shared.number <- 0
comp.shared.tab.A1.f$persist.number <- 0
comp.shared.tab.A1.f$shared.abund <- 0
comp.shared.tab.A1.f$persist.abund <- 0

for (i in as.character(list.A1)){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication2== i & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1A1")])
  shared.otu.name <-intersect(otu.name,otu.G.141.A1.f)
  comp.shared.tab.A1.f$shared.number[which(comp.shared.tab.A1.f$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.A1.f, persistent.fun))
  comp.shared.tab.A1.f$persist.number[which(comp.shared.tab.A1.f$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.fun.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.A1.f$shared.abund[which(comp.shared.tab.A1.f$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.fun.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.A1.f$persist.abund[which(comp.shared.tab.A1.f$comp == i)]<-persist.abund
}
comp.shared.tab.A1.f$Plot <- "A1"


##B2
list.B2<-unique(melt.fun.clean.rel$Replication2[which(!(melt.fun.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.fun.clean.rel$Plot == "UF1B2")])

comp.shared.tab.B2.f <- data.frame(comp = as.character(as.factor(list.B2)))
comp.shared.tab.B2.f$shared.number <- 0
comp.shared.tab.B2.f$persist.number <- 0
comp.shared.tab.B2.f$shared.abund <- 0
comp.shared.tab.B2.f$persist.abund <- 0

for (i in as.character(list.B2)){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication2== i & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1B2")])
  shared.otu.name <-intersect(otu.name,otu.G.141.B2.f)
  comp.shared.tab.B2.f$shared.number[which(comp.shared.tab.B2.f$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.B2.f, persistent.fun))
  comp.shared.tab.B2.f$persist.number[which(comp.shared.tab.B2.f$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.fun.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.B2.f$shared.abund[which(comp.shared.tab.B2.f$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.fun.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.B2.f$persist.abund[which(comp.shared.tab.B2.f$comp == i)]<-persist.abund
}
comp.shared.tab.B2.f$Plot <- "B2"


##C3
list.C3<-unique(melt.fun.clean.rel$Replication2[which(!(melt.fun.clean.rel$Replication2 %in% c("G_0","G_76","G_90","G_106","G_120","G_141"))& melt.fun.clean.rel$Plot == "UF1C3")])

comp.shared.tab.C3.f <- data.frame(comp = as.character(as.factor(list.C3)))
comp.shared.tab.C3.f$shared.number <- 0
comp.shared.tab.C3.f$persist.number <- 0
comp.shared.tab.C3.f$shared.abund <- 0
comp.shared.tab.C3.f$persist.abund <- 0

for (i in as.character(list.C3)){
  otu.name <- unique(melt.fun.clean.rel$OTU[which(melt.fun.clean.rel$Replication2== i & melt.fun.clean.rel$Abundance > 0 & melt.fun.clean.rel$Plot == "UF1C3")])
  shared.otu.name <-intersect(otu.name,otu.G.141.C3.f)
  comp.shared.tab.C3.f$shared.number[which(comp.shared.tab.C3.f$comp == i)]<-length(shared.otu.name)
  persist.otu.name <-Reduce(intersect,list(otu.name,otu.G.141.C3.f, persistent.fun))
  comp.shared.tab.C3.f$persist.number[which(comp.shared.tab.C3.f$comp == i)]<-length(persist.otu.name)
  shared.abund<-esti.mean.abun(melt.fun.clean.rel.mean,shared.otu.name,i)
  comp.shared.tab.C3.f$shared.abund[which(comp.shared.tab.C3.f$comp == i)]<-shared.abund
  persist.abund<-esti.mean.abun(melt.fun.clean.rel.mean,persist.otu.name,i)
  comp.shared.tab.C3.f$persist.abund[which(comp.shared.tab.C3.f$comp == i)]<-persist.abund
}
comp.shared.tab.C3.f$Plot <- "C3"

comp.shared.tab.fun<-rbind(comp.shared.tab.A1.f, comp.shared.tab.B2.f, comp.shared.tab.C3.f)

comp.shared.tab.fun$Kingdom <- "Fungi"


### Calculate mean and standard deviation values
comp.shared.tab.bac.mean <- comp.shared.tab.bac %>% group_by(comp) %>% summarise(mean.shared = mean(shared.number), mean.persist = mean(persist.number), mean.shared.abund=mean(shared.abund), mean.persist.abund=mean(persist.abund), 
                                                                                 sd.shared = sd(shared.number), sd.persist = sd(persist.number), sd.shared.abund=sd(shared.abund), sd.persist.abund=sd(persist.abund))


comp.shared.tab.fun.mean <- comp.shared.tab.fun %>% group_by(comp) %>% summarise(mean.shared = mean(shared.number), mean.persist = mean(persist.number), mean.shared.abund=mean(shared.abund), mean.persist.abund=mean(persist.abund), 
                                                                                 sd.shared = sd(shared.number), sd.persist = sd(persist.number), sd.shared.abund=sd(shared.abund), sd.persist.abund=sd(persist.abund))




### plotting

write.xlsx(comp.shared.tab.bac.mean,"comp.shared.tab.bac.mean_original.xlsx")
write.xlsx(comp.shared.tab.fun.mean,"comp.shared.tab.fun.mean_original.xlsx")

comp.shared.tab.bac.mean <- read.xlsx("comp.shared.tab.bac.mean_scaled.xlsx",1)
comp.shared.tab.bac.mean$Compartment <- factor(comp.shared.tab.bac.mean$Compartment, levels = c("L","S","R","RS","BS"))
comp.shared.tab.bac.mean$Category <- factor(comp.shared.tab.bac.mean$Category, levels = c("mean.shared","mean.persist","mean.shared.abund","mean.persist.abund"))

g<-ggplot(comp.shared.tab.bac.mean , aes(x= Compartment,fill = Category,y= Mean)) + 
  geom_bar(width = 0.5,position = "dodge2", stat = "identity")+ 
  facet_wrap(~Days, ncol = 1)+ geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),
                                                                      position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")+theme(aspect.ratio = 0.4)#+ scale_color_manual(values = c("darkred", "steelblue"))
g+scale_y_continuous(sec.axis = sec_axis(~./31, name="Cumulative relative abundance\n"))


comp.shared.tab.fun.mean <- read.xlsx("comp.shared.tab.fun.mean_scaled.xlsx",1)
comp.shared.tab.fun.mean$Compartment <- factor(comp.shared.tab.fun.mean$Compartment, levels = c("L","S","R","RS","BS"))
comp.shared.tab.fun.mean$Category <- factor(comp.shared.tab.fun.mean$Category, levels = c("mean.shared","mean.persist","mean.shared.abund","mean.persist.abund"))

g<-ggplot(comp.shared.tab.fun.mean , aes(x= Compartment,fill = Category,y= Mean)) + 
  geom_bar(width = 0.5,position = "dodge2", stat = "identity")+ 
  facet_wrap(~Days, ncol = 1)+ geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),
                                             position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Number of OTUs\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "bottom")+theme(aspect.ratio = 0.4)#+ scale_color_manual(values = c("darkred", "steelblue"))
g+scale_y_continuous(sec.axis = sec_axis(~./56, name="Cumulative relative abundance\n"))


