##load the count tables

## ---- herit0
#BiocManager::install("mixOmics")
library(vegan)
library(MASS)
library(mixOmics)
library(plyr)
source("heritability_functions_optimize.R")

# ##set up loops
# genes=c("ITS", "16S")
# exps=c("ull", "rat", "ada", "ram")
# years=c(2012, 2013)
# 
# ##load the taxonomy
# tb=readRDS("./data/taxonomy_otus_16S_v2.rds")
# tf=readRDS("./data/taxonomy_otus_ITS_v2.rds")
# tax=rbind(tb, tf)


## ---- end-of-herit0

## ---- herit1
##compute the heritabilities of individual OTUs
bac.clean.ss.17.G <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
bac.clean.ss.17.G <- phyloseq::filter_taxa(bac.clean.ss.17.G, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.G <- subset_samples(bac.clean.ss.f, Year == "year_2018" & Compartment == "Seed")
bac.clean.ss.18.G <- phyloseq::filter_taxa(bac.clean.ss.18.G, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G.harv <- subset_samples(bac.clean.ss.f, Developmental_stage == "Harvest" & Compartment == "Seed")
bac.clean.ss.G.harv <- phyloseq::filter_taxa(bac.clean.ss.G.harv, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G <- subset_samples(bac.clean.ss.f, Compartment == "Seed")
bac.clean.ss.G <- phyloseq::filter_taxa(bac.clean.ss.G, function(x) sum(x) != 0, TRUE)



bac.clean.ss.R <- subset_samples(bac.clean.ss.f, Compartment == "Root")
bac.clean.ss.R <- phyloseq::filter_taxa(bac.clean.ss.R, function(x) sum(x) != 0, TRUE)



fun.clean.ss.17.G <- subset_samples(fun.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
fun.clean.ss.17.G <- phyloseq::filter_taxa(fun.clean.ss.17.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.G <- subset_samples(fun.clean.ss.f, Year == "year_2018" & Compartment == "Grain")
fun.clean.ss.18.G <- phyloseq::filter_taxa(fun.clean.ss.18.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G.harv <- subset_samples(fun.clean.ss.f, Developmental_stage == "Harvest" & Compartment %in% c("Seed","Grain"))
fun.clean.ss.G.harv <- phyloseq::filter_taxa(fun.clean.ss.G.harv, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G <- subset_samples(fun.clean.ss.f, Compartment == "Seed")
fun.clean.ss.G <- phyloseq::filter_taxa(fun.clean.ss.G, function(x) sum(x) != 0, TRUE)

otu.list.full <- rbind(bac.list, fun.list)


## Other compartments or sections
### Bacteria
bac.clean.ss.18.L1 <- subset_samples(bac.clean.ss.f, Microhabitat == "L1")
bac.clean.ss.18.L1 <- phyloseq::filter_taxa(bac.clean.ss.18.L1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.L2 <- subset_samples(bac.clean.ss.f, Microhabitat == "L2")
bac.clean.ss.18.L2 <- phyloseq::filter_taxa(bac.clean.ss.18.L2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.L3 <- subset_samples(bac.clean.ss.f, Microhabitat == "L3")
bac.clean.ss.18.L3 <- phyloseq::filter_taxa(bac.clean.ss.18.L3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.FL <- subset_samples(bac.clean.ss.f, Microhabitat == "FL")
bac.clean.ss.18.FL <- phyloseq::filter_taxa(bac.clean.ss.18.FL, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S1 <- subset_samples(bac.clean.ss.f, Microhabitat == "S1")
bac.clean.ss.18.S1 <- phyloseq::filter_taxa(bac.clean.ss.18.S1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S2 <- subset_samples(bac.clean.ss.f, Microhabitat == "S2")
bac.clean.ss.18.S2 <- phyloseq::filter_taxa(bac.clean.ss.18.S2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S3 <- subset_samples(bac.clean.ss.f, Microhabitat == "S3")
bac.clean.ss.18.S3 <- phyloseq::filter_taxa(bac.clean.ss.18.S3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S4 <- subset_samples(bac.clean.ss.f, Microhabitat == "S4")
bac.clean.ss.18.S4 <- phyloseq::filter_taxa(bac.clean.ss.18.S4, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S5 <- subset_samples(bac.clean.ss.f, Microhabitat == "S5")
bac.clean.ss.18.S5 <- phyloseq::filter_taxa(bac.clean.ss.18.S5, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S6 <- subset_samples(bac.clean.ss.f, Microhabitat == "S6")
bac.clean.ss.18.S6 <- phyloseq::filter_taxa(bac.clean.ss.18.S6, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S7 <- subset_samples(bac.clean.ss.f, Microhabitat == "S7")
bac.clean.ss.18.S7 <- phyloseq::filter_taxa(bac.clean.ss.18.S7, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S8 <- subset_samples(bac.clean.ss.f, Microhabitat == "S8")
bac.clean.ss.18.S8 <- phyloseq::filter_taxa(bac.clean.ss.18.S8, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.S9 <- subset_samples(bac.clean.ss.f, Microhabitat == "S9")
bac.clean.ss.18.S9 <- phyloseq::filter_taxa(bac.clean.ss.18.S9, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.R <- subset_samples(bac.clean.ss.f, Microhabitat == "R")
bac.clean.ss.18.R <- phyloseq::filter_taxa(bac.clean.ss.18.R, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.RS <- subset_samples(bac.clean.ss.f, Microhabitat == "RS")
bac.clean.ss.18.RS <- phyloseq::filter_taxa(bac.clean.ss.18.RS, function(x) sum(x) != 0, TRUE)

bac.clean.ss.18.BS <- subset_samples(bac.clean.ss.f, Microhabitat == "BS")
bac.clean.ss.18.BS <- phyloseq::filter_taxa(bac.clean.ss.18.BS, function(x) sum(x) != 0, TRUE)


##Fungi
fun.clean.ss.18.L1 <- subset_samples(fun.clean.ss.f, Microhabitat == "L1")
fun.clean.ss.18.L1 <- phyloseq::filter_taxa(fun.clean.ss.18.L1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.L2 <- subset_samples(fun.clean.ss.f, Microhabitat == "L2")
fun.clean.ss.18.L2 <- phyloseq::filter_taxa(fun.clean.ss.18.L2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.L3 <- subset_samples(fun.clean.ss.f, Microhabitat == "L3")
fun.clean.ss.18.L3 <- phyloseq::filter_taxa(fun.clean.ss.18.L3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.FL <- subset_samples(fun.clean.ss.f, Microhabitat == "FL")
fun.clean.ss.18.FL <- phyloseq::filter_taxa(fun.clean.ss.18.FL, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S1 <- subset_samples(fun.clean.ss.f, Microhabitat == "S1")
fun.clean.ss.18.S1 <- phyloseq::filter_taxa(fun.clean.ss.18.S1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S2 <- subset_samples(fun.clean.ss.f, Microhabitat == "S2")
fun.clean.ss.18.S2 <- phyloseq::filter_taxa(fun.clean.ss.18.S2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S3 <- subset_samples(fun.clean.ss.f, Microhabitat == "S3")
fun.clean.ss.18.S3 <- phyloseq::filter_taxa(fun.clean.ss.18.S3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S4 <- subset_samples(fun.clean.ss.f, Microhabitat == "S4")
fun.clean.ss.18.S4 <- phyloseq::filter_taxa(fun.clean.ss.18.S4, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S5 <- subset_samples(fun.clean.ss.f, Microhabitat == "S5")
fun.clean.ss.18.S5 <- phyloseq::filter_taxa(fun.clean.ss.18.S5, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S6 <- subset_samples(fun.clean.ss.f, Microhabitat == "S6")
fun.clean.ss.18.S6 <- phyloseq::filter_taxa(fun.clean.ss.18.S6, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S7 <- subset_samples(fun.clean.ss.f, Microhabitat == "S7")
fun.clean.ss.18.S7 <- phyloseq::filter_taxa(fun.clean.ss.18.S7, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S8 <- subset_samples(fun.clean.ss.f, Microhabitat == "S8")
fun.clean.ss.18.S8 <- phyloseq::filter_taxa(fun.clean.ss.18.S8, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.S9 <- subset_samples(fun.clean.ss.f, Microhabitat == "S9")
fun.clean.ss.18.S9 <- phyloseq::filter_taxa(fun.clean.ss.18.S9, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.R <- subset_samples(fun.clean.ss.f, Microhabitat == "R")
fun.clean.ss.18.R <- phyloseq::filter_taxa(fun.clean.ss.18.R, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.RS <- subset_samples(fun.clean.ss.f, Microhabitat == "RS")
fun.clean.ss.18.RS <- phyloseq::filter_taxa(fun.clean.ss.18.RS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18.BS <- subset_samples(fun.clean.ss.f, Microhabitat == "BS")
fun.clean.ss.18.BS <- phyloseq::filter_taxa(fun.clean.ss.18.BS, function(x) sum(x) != 0, TRUE)



###Analysis
herit_model.develop <- function(physeqs){Notus=data.frame()
count=1
uotus=c()

ct=t(data.frame(otu_table(physeqs)))  # count table
s=data.frame(sample_data(physeqs)) #Meta data
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
subct=ct
subs=s
tot=sum(subct)
rs=colSums(subct)
subct=as.matrix(subct[, rs>=(0.0001*tot)])
uotus=c(uotus, colnames(subct))
Notus[count,"Notus"]=ncol(subct)
Notus[count,"Nsamples"]=nrow(subct)
r=apply(subct, 1, rank, ties="average")
#Notus[count, "meanRank"]=t(apply(r, 2, mean))
lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
h2=apply(lsubct,2, H2nokin, id=subs$Developmental_stage)
#saveRDS(h2, paste("./res/h2_", gene, "_", e, "_", y,".rds", sep=""))
count=count+1
h2.df <- data.frame(h2)
h2.df$OTU <- rownames(h2.df)
herit.tab<-merge(otu.list.full, h2.df, by = c("OTU"="OTU"))
herit.tab<- herit.tab %>% arrange(desc(h2))
return(herit.tab)
}

bac.G.devel<-herit_model.develop(bac.clean.ss.18.G)
bac.L1.devel<-herit_model.develop(bac.clean.ss.18.L1)
bac.L2.devel<-herit_model.develop(bac.clean.ss.18.L2)
bac.L3.devel<-herit_model.develop(bac.clean.ss.18.L3)
bac.FL.devel<-herit_model.develop(bac.clean.ss.18.FL)

bac.S1.devel<-herit_model.develop(bac.clean.ss.18.S1)
bac.S2.devel<-herit_model.develop(bac.clean.ss.18.S2)
bac.S3.devel<-herit_model.develop(bac.clean.ss.18.S3)
bac.S4.devel<-herit_model.develop(bac.clean.ss.18.S4)
bac.S5.devel<-herit_model.develop(bac.clean.ss.18.S5)
bac.S6.devel<-herit_model.develop(bac.clean.ss.18.S6)
bac.S7.devel<-herit_model.develop(bac.clean.ss.18.S7)
bac.S8.devel<-herit_model.develop(bac.clean.ss.18.S8)
bac.S9.devel<-herit_model.develop(bac.clean.ss.18.S9)

bac.R.devel<-herit_model.develop(bac.clean.ss.18.R)
bac.RS.devel<-herit_model.develop(bac.clean.ss.18.RS)
bac.BS.devel<-herit_model.develop(bac.clean.ss.18.BS)

fun.G.devel<-herit_model.develop(fun.clean.ss.18.G)
fun.L1.devel<-herit_model.develop(fun.clean.ss.18.L1)
fun.L2.devel<-herit_model.develop(fun.clean.ss.18.L2)
fun.L3.devel<-herit_model.develop(fun.clean.ss.18.L3)
fun.FL.devel<-herit_model.develop(fun.clean.ss.18.FL)

fun.S1.devel<-herit_model.develop(fun.clean.ss.18.S1)
fun.S2.devel<-herit_model.develop(fun.clean.ss.18.S2)
fun.S3.devel<-herit_model.develop(fun.clean.ss.18.S3)
fun.S4.devel<-herit_model.develop(fun.clean.ss.18.S4)
fun.S5.devel<-herit_model.develop(fun.clean.ss.18.S5)
fun.S6.devel<-herit_model.develop(fun.clean.ss.18.S6)
fun.S7.devel<-herit_model.develop(fun.clean.ss.18.S7)
fun.S8.devel<-herit_model.develop(fun.clean.ss.18.S8)
fun.S9.devel<-herit_model.develop(fun.clean.ss.18.S9)

fun.R.devel<-herit_model.develop(fun.clean.ss.18.R)
fun.RS.devel<-herit_model.develop(fun.clean.ss.18.RS)
fun.BS.devel<-herit_model.develop(fun.clean.ss.18.BS)

write.csv(bac.G.devel, "Bacteria_fixed effect_developmental stage_2018 developing seeds.csv")
write.csv(fun.G.devel, "Fungi_fixed effect_developmental stage_2018 developing seeds.csv")


### Proportion of affected OTUs
length(bac.BS.devel$OTU[which(bac.BS.devel$h2 >= 0.15)])
length(bac.BS.devel$OTU)

length(fun.G.devel$OTU[which(fun.G.devel$h2 >= 0.15)])
length(fun.G.devel$OTU)

microhabitat.list <- c("BS","RS","R","S1","S2","S3","S4","S5","S6","S7","S8","S9","L1","L2","L3","FL","G")


####bootstrap model #100
herit_model.develop.boot <- function(physeqs){Notus=data.frame()
count=1
uotus=c()

ct=t(data.frame(otu_table(physeqs)))  # count table
s=data.frame(sample_data(physeqs)) #Meta data
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
subct=ct
subs=s
tot=sum(subct)
rs=colSums(subct)
subct=as.matrix(subct[, rs>=(0.0001*tot)])
uotus=c(uotus, colnames(subct))
Notus[count,"Notus"]=ncol(subct)
Notus[count,"Nsamples"]=nrow(subct)
r=apply(subct, 1, rank, ties="average")
#Notus[count, "meanRank"]=t(apply(r, 2, mean))
lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
h2=apply(lsubct,2, H2bootNokin, id=subs$Age, rep = subs$Replication)
#saveRDS(h2, paste("./res/h2_", gene, "_", e, "_", y,".rds", sep=""))
count=count+1
# h2.df <- data.frame(h2)
# h2.df$OTU <- rownames(h2.df)
# herit.tab<-merge(otu.list.full, h2.df, by = c("OTU"="OTU"))
#herit.tab<- herit.tab %>% dplyr::arrange(desc(h2))
return(h2)
}

bac.G.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.G)
bac.L1.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.L1)
bac.L2.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.L2)
bac.L3.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.L3)
bac.FL.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.FL)

bac.S1.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S1)
bac.S2.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S2)
bac.S3.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S3)
bac.S4.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S4)
bac.S5.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S5)
bac.S6.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S6)
bac.S7.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S7)
bac.S8.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S8)
bac.S9.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.S9)

bac.R.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.R)
bac.RS.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.RS)
bac.BS.devel.boot<-herit_model.develop.boot(bac.clean.ss.18.BS)

write.csv(bac.G.devel.boot, "G_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.FL.devel.boot, "FL_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.L3.devel.boot, "L3_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.L2.devel.boot, "L2_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.L1.devel.boot, "L1_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S9.devel.boot, "S9_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S8.devel.boot, "S8_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S7.devel.boot, "S7_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S6.devel.boot, "S6_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S5.devel.boot, "S5_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S4.devel.boot, "S4_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S3.devel.boot, "S3_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S2.devel.boot, "S2_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.S1.devel.boot, "S1_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.R.devel.boot, "R_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.RS.devel.boot, "RS_Bacteria_h2_bootstrap_age.csv")
write.csv(bac.BS.devel.boot, "BS_Bacteria_h2_bootstrap_age.csv")

# get_h2_tab<-function(h2.boot.tab){
# h2.boot.tab<-data.frame(h2.boot.tab)
# otu.mean<-colMeans(h2.boot.tab)
# otu.mean <- data.frame(otu.mean)
# names(otu.mean)[1] <- "mean_h2"
# otu.mean$h2_sd <- matrixStats::colSds(as.matrix(h2.boot.tab[sapply(h2.boot.tab, is.numeric)]))
# otu.mean$OTU <- rownames(otu.mean)
# h2.tab<-merge(otu.list.full, otu.mean, by = c("OTU" = "OTU"))
# return(h2.tab)}
# 
# h2.tab.G.bac<-get_h2_tab(bac.G.devel.boot)
# h2.tab.FL.bac<-get_h2_tab(bac.FL.devel.boot)
# h2.tab.L3.bac<-get_h2_tab(bac.L3.devel.boot)
# h2.tab.L2.bac<-get_h2_tab(bac.L2.devel.boot)
# h2.tab.L1.bac<-get_h2_tab(bac.L1.devel.boot)


fun.G.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.G)
fun.L1.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.L1)
fun.L2.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.L2)
fun.L3.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.L3)
fun.FL.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.FL)

fun.S1.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S1)
fun.S2.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S2)
fun.S3.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S3)
fun.S4.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S4)
fun.S5.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S5)
fun.S6.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S6)
fun.S7.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S7)
fun.S8.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S8)
fun.S9.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.S9)

fun.R.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.R)
fun.RS.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.RS)
fun.BS.devel.boot<-herit_model.develop.boot(fun.clean.ss.18.BS)


write.csv(fun.G.devel.boot, "G_Fungi_h2_bootstrap_age.csv")
write.csv(fun.FL.devel.boot, "FL_Fungi_h2_bootstrap_age.csv")
write.csv(fun.L3.devel.boot, "L3_Fungi_h2_bootstrap_age.csv")
write.csv(fun.L2.devel.boot, "L2_Fungi_h2_bootstrap_age.csv")
write.csv(fun.L1.devel.boot, "L1_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S9.devel.boot, "S9_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S8.devel.boot, "S8_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S7.devel.boot, "S7_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S6.devel.boot, "S6_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S5.devel.boot, "S5_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S4.devel.boot, "S4_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S3.devel.boot, "S3_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S2.devel.boot, "S2_Fungi_h2_bootstrap_age.csv")
write.csv(fun.S1.devel.boot, "S1_Fungi_h2_bootstrap_age.csv")
write.csv(fun.R.devel.boot, "R_Fungi_h2_bootstrap_age.csv")
write.csv(fun.RS.devel.boot, "RS_Fungi_h2_bootstrap_age.csv")
write.csv(fun.BS.devel.boot, "BS_Fungi_h2_bootstrap_age.csv")


herit_model.maternal <- function(physeqs){Notus=data.frame()
count=1
uotus=c()

ct=t(data.frame(otu_table(physeqs)))  # count table
s=data.frame(sample_data(physeqs)) #Meta data
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
subct=ct
subs=s
tot=sum(subct)
rs=colSums(subct)
subct=as.matrix(subct[, rs>=(0.0001*tot)])
uotus=c(uotus, colnames(subct))
Notus[count,"Notus"]=ncol(subct)
Notus[count,"Nsamples"]=nrow(subct)
r=apply(subct, 1, rank, ties="average")
#Notus[count, "meanRank"]=t(apply(r, 2, mean))
lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
h2=apply(lsubct,2, H2nokin, id=subs$id, rep = subs$Replication)
#saveRDS(h2, paste("./res/h2_", gene, "_", e, "_", y,".rds", sep=""))
count=count+1
h2.df <- data.frame(h2)
h2.df$OTU <- rownames(h2.df)
herit.tab<-merge(otu.list.full, h2.df, by = c("OTU"="OTU"))
herit.tab<- herit.tab %>% arrange(desc(h2))
return(herit.tab)
}

fun.G.maternal<-herit_model.maternal(fun.clean.ss.18.G)
bac.G.maternal<-herit_model.maternal(bac.clean.ss.18.G)

##Genotype and developmental stage
herit_model.develop.geno <- function(physeqs){Notus=data.frame()
count=1
uotus=c()

ct=t(data.frame(otu_table(physeqs)))  # count table
s=data.frame(sample_data(physeqs)) #Meta data
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
subct=ct
subs=s
tot=sum(subct)
rs=colSums(subct)
subct=as.matrix(subct[, rs>=(0.0001*tot)])
uotus=c(uotus, colnames(subct))
Notus[count,"Notus"]=ncol(subct)
Notus[count,"Nsamples"]=nrow(subct)
r=apply(subct, 1, rank, ties="average")
#Notus[count, "meanRank"]=t(apply(r, 2, mean))
lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
h2=apply(lsubct,2, H2nokin.multi, id1=subs$Developmental_stage,id2 =subs$Genotype, rep = subs$Replication)
#saveRDS(h2, paste("./res/h2_", gene, "_", e, "_", y,".rds", sep=""))
count=count+1
h2.df <- data.frame(h2)
h2.df$OTU <- rownames(h2.df)
herit.tab<-merge(otu.list.full, h2.df, by = c("OTU"="OTU"))
herit.tab<- herit.tab %>% arrange(desc(h2))
return(herit.tab)
}
fun.G.devel.gen<-herit_model.develop.geno(fun.clean.ss.G)
bac.G.devel.gen<-herit_model.develop.geno(bac.clean.ss.G)

##2017+2018

##set up loops
years=c("year_2017", "year_2018")

##compute the heritabilities of individual OTUs

Notus=data.frame()
count=1
uotus=c()
    ct=t(data.frame(otu_table(bac.clean.ss.R)))
    s=data.frame(sample_data(bac.clean.ss.R))
    s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
    s$Year=as.factor(s$Year)
    #s$plate=as.factor(s$plate)
        for(y in years){
            subct=ct[s$Year==y,]
            subs=droplevels(s[s$Year==y,])
            tot=sum(subct)
            rs=colSums(subct)
            subct=as.matrix(subct[, rs>=(0.0001*tot)])
            uotus=c(uotus, colnames(subct))
            #Notus[count,"exp"]=e
            Notus[count,"year"]=y
            Notus[count,"gene"]="16S"
            Notus[count,"Notus"]=ncol(subct)
            Notus[count,"Nsamples"]=nrow(subct)
            r=apply(subct, 1, rank, ties="average")
            #Notus[count, "meanRank"]=t(apply(r, 2, mean))
            lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
            h2=apply(lsubct,2, H2nokin, id=subs$Developmental_stage) #H2bootNokin H2nokin
            #saveRDS(h2, paste("./res/h2_", "16S", "_", e, "_", y,".rds", sep=""))
            count=count+1
        }


h2.df <- data.frame(h2)
h2.df$OTU <- rownames(h2.df)
bac.herit<-merge(bac.list, h2.df, by = c("OTU"="OTU"))
fun.herit<-merge(fun.list, h2.df, by = c("OTU"="OTU"))

bac.herit%>% arrange(desc(bac.herit$h2))
write.csv(bac.herit,"Heritability of seed bacteria_2017 and 2018_SW_nonboot.csv")
write.csv(h2.df,"Heritability of seed bacteria_2017 and 2018_SW_boot.csv")

fun.herit%>% arrange(desc(fun.herit$h2))
write.csv(fun.herit,"Heritability of seed fungi_2017 and 2018_SW_nonboot.csv")
write.csv(h2.df,"Heritability of seed bacteria_2017 and 2018_SW_boot.csv")



### 2017 Genotype
bac.clean.ss.17.G <- subset_samples(bac.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
bac.clean.ss.17.G <- phyloseq::filter_taxa(bac.clean.ss.17.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.17.G <- subset_samples(fun.clean.ss.f, Year == "year_2017" & Compartment == "Seed")
fun.clean.ss.17.G <- phyloseq::filter_taxa(fun.clean.ss.17.G, function(x) sum(x) != 0, TRUE)


exps=c("Akibare", "Gohyangchalbyeo", "Daeanbyeo")

Notus=data.frame()
count=1
uotus=c()
    ct=t(data.frame(otu_table(fun.clean.ss.17.G)))
    s=data.frame(sample_data(fun.clean.ss.17.G))
    s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
    #s$Year=as.factor(s$Year)
    for(e in exps){
        for(y in years){
            subct=ct[s$Genotype==e,]
            subs=droplevels(s[s$Genotype==e,])
            tot=sum(subct)
            rs=colSums(subct)
            subct=as.matrix(subct[, rs>=(0.0001*tot)])
            uotus=c(uotus, colnames(subct))
            Notus[count,"Genotype"]=e
            #Notus[count,"Year"]=y
            Notus[count,"gene"]="16S"
            Notus[count,"Notus"]=ncol(subct)
            Notus[count,"Nsamples"]=nrow(subct)
            r=apply(subct, 1, rank, ties="average")
            #Notus[count, "meanRank"]=t(apply(r, 2, mean))
            lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
            h2=apply(lsubct,2, H2nokin, id=subs$id)
            #saveRDS(h2, paste("./res/h2_", "16S", "_", e, "_", y,".rds", sep=""))
            count=count+1
        }
    }

    
    h2.df <- data.frame(h2)
    h2.df$OTU <- rownames(h2.df)
    bac.herit<-merge(bac.list, h2.df, by = c("OTU"="OTU"))
    fun.herit<-merge(fun.list, h2.df, by = c("OTU"="OTU"))
    
    bac.herit%>% arrange(desc(bac.herit$h2))
    write.csv(bac.herit,"Heritability of seed bacteria_2017_Genotype_nonboot.csv")
    
    fun.herit%>% arrange(desc(fun.herit$h2))
    write.csv(fun.herit,"Heritability of seed fungi_2017_Genotype_nonboot.csv")
    


saveRDS(Notus, "./res/Notus_Nsamples.rds")

##### Heritability of community phenotypes
##Bray Curtis distance based PCoA
bac.clean.ss.18.G

b.meta.18<-sample_data(bac.clean.ss.18)
b.meta.18$Days <- as.numeric(as.character(b.meta.18$Days))
sample_data(bac.clean.log.18) <- sample_data(b.meta.18)

bac.clean.log.18.G <- subset_samples(bac.clean.log.18, Compartment == "Seed")
bac.clean.log.18.G <- phyloseq::filter_taxa(bac.clean.log.18.G, function(x) sum(x) != 0, TRUE)

bray.bac.18.G <-  ordinate(bac.clean.log.18.G , 'PCoA', 'bray')
bray.Pcos.bac <- bray.bac.18.G$vectors


f.meta.18<-sample_data(fun.clean.ss.18)
f.meta.18$Days <- as.numeric(as.character(f.meta.18$Days))
sample_data(fun.clean.log.18) <- sample_data(f.meta.18)

fun.clean.log.18.G <- subset_samples(fun.clean.log.18, Compartment == "Grain")
fun.clean.log.18.G <- phyloseq::filter_taxa(fun.clean.log.18.G, function(x) sum(x) != 0, TRUE)

bray.fun.18.G <-  ordinate(fun.clean.log.18.G , 'PCoA', 'bray')
bray.Pcos.fun <- bray.fun.18.G$vectors



bray2.bac <-  ordinate(bac.clean.log.18, 'PCoA', 'bray')

### Jaccard distance based PCoA
jaccard.bac.18.G <-  ordinate(bac.clean.log.18.G , 'PCoA', 'jaccard')
jaccard.Pcos.bac <- jaccard.bac.18.G$vectors


jaccard.fun.18.G <-  ordinate(fun.clean.log.18.G , 'PCoA', 'jaccard')
jaccard.Pcos.fun <- jaccard.fun.18.G$vectors




Notus=data.frame()
count=1
uotus=c()
ct=data.frame(bray.Pcos.fun)
s=data.frame(sample_data(fun.clean.ss.18.G))
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
s$Year=as.factor(s$Year)
#s$plate=as.factor(s$plate)
    subct=ct
    subs=s
    tot=sum(subct)
    rs=colSums(subct)
    subct=as.matrix(subct[, rs>=(0.0001*tot)])
    uotus=c(uotus, colnames(subct))
    #Notus[count,"exp"]=e
    #Notus[count,"year"]=y
    Notus[count,"gene"]="16S"
    Notus[count,"Notus"]=ncol(subct)
    Notus[count,"Nsamples"]=nrow(subct)
    r=apply(subct, 1, rank, ties="average")
    #Notus[count, "meanRank"]=t(apply(r, 2, mean))
    lsubct=logratio.transfo(subct+1, logratio="CLR", offset=0)
    h2=apply(lsubct,2, H2nokin, id=subs$Developmental_stage, rep = subs$Replication) 
   #H2bootNokin H2nokin
    #saveRDS(h2, paste("./res/h2_", "16S", "_", e, "_", y,".rds", sep=""))
    count=count+1

## ---- end-of-herit1

## ---- datasize1

Notus=readRDS("./res/Notus_Nsamples.rds")
kable(Notus, caption="Number of OTUs and number of samples for each location/year combination")

## ---- end-of-datasize1

## ---- datasize2
#print(paste("Nsamples on average=", mean(Notus$Nsamples), " / from= ", min(Notus$Nsamples), "/ to= ", max(Notus$Nsamples), sep=""))
kable(ddply(Notus, "gene", function(x){c(Av_Notus=mean(x$Notus), from=min(x$Notus), to=max(x$Notus))}), caption="Number of OTUs and samples in each experiment used on which heritability analysis were performed.")
## ---- end-of-datasize2

## ---- distherit

prop_herit=c()
pdf("./figures/otu_heritability_distribution.pdf", paper="special", height=6, width=4, pointsize=8)
## m=c(1,1,2,3,4,4,5,6,rep(1, 4), rep(4, 4))
## mall=rep(m, 4)
## add=sort(rep(seq(0, (6*3), 6), 16))
## mat=matrix(mall+add, ncol=8, byrow=T)
## layout(mat, width=c(1, 1, 3, 2, 1, 1, 3, 2), heights=rep(c(3, 1), 4))
m=c(1,1,2,3,3,4,rep(1, 3), rep(3, 3))
mall=rep(m, 4)
add=sort(rep(seq(0, (4*3), 4), 12))
mat=matrix(mall+add, ncol=6, byrow=T)
layout(mat, width=c(1, 1, 2, 1, 1, 2), heights=rep(c(3, 1), 4))
#layout.show(16)
count=1
for(e in exps){
    e2=c("SU", "SR", "NM", "NA")[match(e,c("ull", "rat", "ram", "ada"))]
    for(y in years){
        par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.6, 0))
        h2_1=readRDS(paste("./res/h2_ITS_", e, "_", y,".rds", sep=""))
        h2_2=readRDS(paste("./res/h2_16S_", e, "_", y,".rds", sep=""))
        h2boot_1=readRDS(paste("./res/h2boot_ITS_", e, "_", y,".rds", sep=""))
        h2boot_2=readRDS(paste("./res/h2boot_16S_", e, "_", y,".rds", sep=""))
        ##reconcile OTUs that are in one but not in the other
        h2_1=h2_1[names(h2_1)%in%colnames(h2boot_1)]
        h2boot_1=h2boot_1[, colnames(h2boot_1)%in%names(h2_1)]
        h2_2=h2_2[names(h2_2)%in%colnames(h2boot_2)]
        h2boot_2=h2boot_2[, colnames(h2boot_2)%in%names(h2_2)]
        h2boot_1=h2boot_1[, names(h2_1)]
        h2boot_2=h2boot_2[, names(h2_2)]
        h2=c(h2_1, h2_2)
        h2boot=cbind(h2boot_1, h2boot_2)
        x=data.frame(otus=names(h2),exp=e, year=y, h2=unlist(h2), medboot=unlist(apply(h2boot,2, median)), CIlow=unlist(apply(h2boot,2,  quantile, 0.025)), CIhigh=unlist(apply(h2boot,2, quantile, 0.975)))
        ##remove OTUs that have heritabilities outside of 95%confidence interval
        x=x[x$h2>=x$CIlow & x$h2<=x$CIhigh,]
        x$org=substring(x$otus, 1, 1)
        prop_herit[paste(e, y, sep="_")]=(100*(nrow(x[x$CIlow>0.001,])/nrow(x)))
        col=c("firebrick2", "gold","darkblue",  "dodgerblue")[match(e,c("ull", "rat", "ada", "ram"))]
        hist(x$h2, probability=T, xlim=c(0, 0.35), ylim=c(0, 70), main=paste(e2, y), breaks=seq(0, 0.35, by=0.01), xlab="heritability", ylab="frequency (bins of 1%)", col=col)
        text(0.22, 2, paste("% of OTUs with low CI\nabove 0.1%: ", round(prop_herit[paste(e, y, sep="_")], 2), "%", sep=""), adj=c(0, 0), cex=0.8, bty="o", col="firebrick3")
        box()
        par(mar=c(3,3, 3.5, 1), mgp=c(1.8, 0.6, 0))
        ## ed=ecdf(x$h2)
        ## plot(1-ed(x$h2), x$h2, main="", cex=0.5, pch=16, col="firebrick", cex.axis=0.8, las=2, xlab="quantiles", ylab="heritability", cex.lab=0.8)
        ## par(mar=c(3, 1.5, 3.5, 1.5), mgp=c(1.5, 0.6, 0))
        ## ##add example on ULL 2012.
        ## if(e=="ull" & y==2012){
        ##     hex=0.05
        ##     z=1-ed(hex) ##function of heritability
        ##     segments(z, -0.01, z, hex, lwd=0.8, lty=1, col="grey60")
        ##     segments(-0.01, hex,z,hex, lwd=0.8, lty=1, col="grey60")
        ## }
        #par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.6, 0))
        boxplot(h2~org, data=x, cex=0.5, pch="", lwd=0.5, cex.axis=0.8, las=2, xlab="")
        stripchart(h2~org, data=x, jitter=0.2, method="jitter", vertical=T, add=T, pch=16, col="Dodgerblue", cex=0.3)
        if(count==1){X=x}else{X=rbind(X, x)}
        count=count+1
    }
}
dev.off()
saveRDS(X, "./res/h2_ind_otus.rds")
saveRDS(prop_herit, "./res/prop_herit.rds")

## ---- end-of-distherit




physeqs <- fun.clean.ss.18.G

ct=t(data.frame(otu_table(physeqs)))  # count table
s=data.frame(sample_data(physeqs)) #Meta data
s$id <- paste0(s$Plot,"_",s$Year,"_",s$Days)
subct=ct
subs=s
tot=sum(subct)
rs=colSums(subct)
subct=as.matrix(subct[, rs>=(0.0001*tot)])
uotus=c(uotus, colnames(subct))

