### Core and inherited OTU
#### Core analysis
map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
rownames(map) <- map$SampleID

map <- subset(map, map$Replication != "Negative")
nrow(map)

sample_data(bac.clean.ss) <- sample_data(map) 

(filt.sample <- sample_sums(bac.clean.ss) > 0)
sum(sample_sums(bac.clean.ss) <= 0)  ## 1 sample discarded
bac.clean.ss.f <- prune_samples(filt.sample, bac.clean.ss)
bac.clean.ss.f 


## fungi
(filt.sample <- sample_sums(fun.clean.ss) > 0)
sum(sample_sums(fun.clean.ss) <= 0)  ## 1 sample discarded
fun.clean.ss.f <- prune_samples(filt.sample, fun.clean.ss)
fun.clean.ss.f 


bac.clean.ss.L1 <- subset_samples(bac.clean.ss.f, Microhabitat == "L1")
bac.clean.ss.L1 <- phyloseq::filter_taxa(bac.clean.ss.L1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.L2 <- subset_samples(bac.clean.ss.f, Microhabitat == "L2")
bac.clean.ss.L2 <- phyloseq::filter_taxa(bac.clean.ss.L2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.L3 <- subset_samples(bac.clean.ss.f, Microhabitat == "L3")
bac.clean.ss.L3 <- phyloseq::filter_taxa(bac.clean.ss.L3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.FL <- subset_samples(bac.clean.ss.f, Microhabitat == "FL")
bac.clean.ss.FL <- phyloseq::filter_taxa(bac.clean.ss.FL, function(x) sum(x) != 0, TRUE)


bac.clean.ss.S1 <- subset_samples(bac.clean.ss.f, Microhabitat == "S1")
bac.clean.ss.S1 <- phyloseq::filter_taxa(bac.clean.ss.S1, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S2 <- subset_samples(bac.clean.ss.f, Microhabitat == "S2")
bac.clean.ss.S2 <- phyloseq::filter_taxa(bac.clean.ss.S2, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S3 <- subset_samples(bac.clean.ss.f, Microhabitat == "S3")
bac.clean.ss.S3 <- phyloseq::filter_taxa(bac.clean.ss.S3, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S4 <- subset_samples(bac.clean.ss.f, Microhabitat == "S4")
bac.clean.ss.S4 <- phyloseq::filter_taxa(bac.clean.ss.S4, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S5 <- subset_samples(bac.clean.ss.f, Microhabitat == "S5")
bac.clean.ss.S5 <- phyloseq::filter_taxa(bac.clean.ss.S5, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S6 <- subset_samples(bac.clean.ss.f, Microhabitat == "S6")
bac.clean.ss.S6 <- phyloseq::filter_taxa(bac.clean.ss.S6, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S7 <- subset_samples(bac.clean.ss.f, Microhabitat == "S7")
bac.clean.ss.S7 <- phyloseq::filter_taxa(bac.clean.ss.S7, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S8 <- subset_samples(bac.clean.ss.f, Microhabitat == "S8")
bac.clean.ss.S8 <- phyloseq::filter_taxa(bac.clean.ss.S8, function(x) sum(x) != 0, TRUE)

bac.clean.ss.S9 <- subset_samples(bac.clean.ss.f, Microhabitat == "S9")
bac.clean.ss.S9 <- phyloseq::filter_taxa(bac.clean.ss.S9, function(x) sum(x) != 0, TRUE)

bac.clean.ss.R <- subset_samples(bac.clean.ss.f, Microhabitat == "R")
bac.clean.ss.R <- phyloseq::filter_taxa(bac.clean.ss.R, function(x) sum(x) != 0, TRUE)

bac.clean.ss.RS <- subset_samples(bac.clean.ss.f, Microhabitat == "RS")
bac.clean.ss.RS <- phyloseq::filter_taxa(bac.clean.ss.RS, function(x) sum(x) != 0, TRUE)

bac.clean.ss.BS <- subset_samples(bac.clean.ss.f, Microhabitat == "BS")
bac.clean.ss.BS <- phyloseq::filter_taxa(bac.clean.ss.BS, function(x) sum(x) != 0, TRUE)

bac.clean.ss.G <- subset_samples(bac.clean.ss, Microhabitat == "G")
bac.clean.ss.G <- phyloseq::filter_taxa(bac.clean.ss.G, function(x) sum(x) != 0, TRUE)

bac.clean.ss.L1.rep <-  merge_samples(bac.clean.ss.L1, "Replication")
bac.clean.ss.L2.rep <-  merge_samples(bac.clean.ss.L2, "Replication")
bac.clean.ss.L3.rep <-  merge_samples(bac.clean.ss.L3, "Replication")
bac.clean.ss.FL.rep <-  merge_samples(bac.clean.ss.FL, "Replication")

bac.clean.ss.S1.rep <-  merge_samples(bac.clean.ss.S1, "Replication")
bac.clean.ss.S2.rep <-  merge_samples(bac.clean.ss.S2, "Replication")
bac.clean.ss.S3.rep <-  merge_samples(bac.clean.ss.S3, "Replication")
bac.clean.ss.S4.rep <-  merge_samples(bac.clean.ss.S4, "Replication")
bac.clean.ss.S5.rep <-  merge_samples(bac.clean.ss.S5, "Replication")
bac.clean.ss.S6.rep <-  merge_samples(bac.clean.ss.S6, "Replication")
bac.clean.ss.S7.rep <-  merge_samples(bac.clean.ss.S7, "Replication")
bac.clean.ss.S8.rep <-  merge_samples(bac.clean.ss.S8, "Replication")
bac.clean.ss.S9.rep <-  merge_samples(bac.clean.ss.S9, "Replication")

bac.clean.ss.R.rep <-  merge_samples(bac.clean.ss.R, "Replication")
bac.clean.ss.RS.rep <-  merge_samples(bac.clean.ss.RS, "Replication")
bac.clean.ss.BS.rep <-  merge_samples(bac.clean.ss.BS, "Replication")

bac.clean.ss.G.rep <-  merge_samples(bac.clean.ss.G, "Replication")


bac.clean.ss.L1.rep.rel <- microbiome::transform(bac.clean.ss.L1.rep, "compositional")
bac.clean.ss.L2.rep.rel <- microbiome::transform(bac.clean.ss.L2.rep, "compositional")
bac.clean.ss.L3.rep.rel <- microbiome::transform(bac.clean.ss.L3.rep, "compositional")
bac.clean.ss.FL.rep.rel <- microbiome::transform(bac.clean.ss.FL.rep, "compositional")

bac.clean.ss.S1.rep.rel <- microbiome::transform(bac.clean.ss.S1.rep, "compositional")
bac.clean.ss.S2.rep.rel <- microbiome::transform(bac.clean.ss.S2.rep, "compositional")
bac.clean.ss.S3.rep.rel <- microbiome::transform(bac.clean.ss.S3.rep, "compositional")
bac.clean.ss.S4.rep.rel <- microbiome::transform(bac.clean.ss.S4.rep, "compositional")
bac.clean.ss.S5.rep.rel <- microbiome::transform(bac.clean.ss.S5.rep, "compositional")
bac.clean.ss.S6.rep.rel <- microbiome::transform(bac.clean.ss.S6.rep, "compositional")
bac.clean.ss.S7.rep.rel <- microbiome::transform(bac.clean.ss.S7.rep, "compositional")
bac.clean.ss.S8.rep.rel <- microbiome::transform(bac.clean.ss.S8.rep, "compositional")
bac.clean.ss.S9.rep.rel <- microbiome::transform(bac.clean.ss.S9.rep, "compositional")

bac.clean.ss.R.rep.rel <- microbiome::transform(bac.clean.ss.R.rep, "compositional")
bac.clean.ss.RS.rep.rel <- microbiome::transform(bac.clean.ss.RS.rep, "compositional")
bac.clean.ss.BS.rep.rel <- microbiome::transform(bac.clean.ss.BS.rep, "compositional")
bac.clean.ss.G.rep.rel <- microbiome::transform(bac.clean.ss.G.rep, "compositional")

##Fungi
fun.clean.ss.L1 <- subset_samples(fun.clean.ss.f, Microhabitat == "L1")
fun.clean.ss.L1 <- phyloseq::filter_taxa(fun.clean.ss.L1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.L2 <- subset_samples(fun.clean.ss.f, Microhabitat == "L2")
fun.clean.ss.L2 <- phyloseq::filter_taxa(fun.clean.ss.L2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.L3 <- subset_samples(fun.clean.ss.f, Microhabitat == "L3")
fun.clean.ss.L3 <- phyloseq::filter_taxa(fun.clean.ss.L3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.FL <- subset_samples(fun.clean.ss.f, Microhabitat == "FL")
fun.clean.ss.FL <- phyloseq::filter_taxa(fun.clean.ss.FL, function(x) sum(x) != 0, TRUE)


fun.clean.ss.S1 <- subset_samples(fun.clean.ss.f, Microhabitat == "S1")
fun.clean.ss.S1 <- phyloseq::filter_taxa(fun.clean.ss.S1, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S2 <- subset_samples(fun.clean.ss.f, Microhabitat == "S2")
fun.clean.ss.S2 <- phyloseq::filter_taxa(fun.clean.ss.S2, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S3 <- subset_samples(fun.clean.ss.f, Microhabitat == "S3")
fun.clean.ss.S3 <- phyloseq::filter_taxa(fun.clean.ss.S3, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S4 <- subset_samples(fun.clean.ss.f, Microhabitat == "S4")
fun.clean.ss.S4 <- phyloseq::filter_taxa(fun.clean.ss.S4, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S5 <- subset_samples(fun.clean.ss.f, Microhabitat == "S5")
fun.clean.ss.S5 <- phyloseq::filter_taxa(fun.clean.ss.S5, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S6 <- subset_samples(fun.clean.ss.f, Microhabitat == "S6")
fun.clean.ss.S6 <- phyloseq::filter_taxa(fun.clean.ss.S6, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S7 <- subset_samples(fun.clean.ss.f, Microhabitat == "S7")
fun.clean.ss.S7 <- phyloseq::filter_taxa(fun.clean.ss.S7, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S8 <- subset_samples(fun.clean.ss.f, Microhabitat == "S8")
fun.clean.ss.S8 <- phyloseq::filter_taxa(fun.clean.ss.S8, function(x) sum(x) != 0, TRUE)

fun.clean.ss.S9 <- subset_samples(fun.clean.ss.f, Microhabitat == "S9")
fun.clean.ss.S9 <- phyloseq::filter_taxa(fun.clean.ss.S9, function(x) sum(x) != 0, TRUE)

fun.clean.ss.R <- subset_samples(fun.clean.ss.f, Microhabitat == "R")
fun.clean.ss.R <- phyloseq::filter_taxa(fun.clean.ss.R, function(x) sum(x) != 0, TRUE)

fun.clean.ss.RS <- subset_samples(fun.clean.ss.f, Microhabitat == "RS")
fun.clean.ss.RS <- phyloseq::filter_taxa(fun.clean.ss.RS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.BS <- subset_samples(fun.clean.ss.f, Microhabitat == "BS")
fun.clean.ss.BS <- phyloseq::filter_taxa(fun.clean.ss.BS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.G <- subset_samples(fun.clean.ss, Microhabitat == "G")
fun.clean.ss.G <- phyloseq::filter_taxa(fun.clean.ss.G, function(x) sum(x) != 0, TRUE)

fun.clean.ss.L1.rep <-  merge_samples(fun.clean.ss.L1, "Replication")
fun.clean.ss.L2.rep <-  merge_samples(fun.clean.ss.L2, "Replication")
fun.clean.ss.L3.rep <-  merge_samples(fun.clean.ss.L3, "Replication")
fun.clean.ss.FL.rep <-  merge_samples(fun.clean.ss.FL, "Replication")

fun.clean.ss.S1.rep <-  merge_samples(fun.clean.ss.S1, "Replication")
fun.clean.ss.S2.rep <-  merge_samples(fun.clean.ss.S2, "Replication")
fun.clean.ss.S3.rep <-  merge_samples(fun.clean.ss.S3, "Replication")
fun.clean.ss.S4.rep <-  merge_samples(fun.clean.ss.S4, "Replication")
fun.clean.ss.S5.rep <-  merge_samples(fun.clean.ss.S5, "Replication")
fun.clean.ss.S6.rep <-  merge_samples(fun.clean.ss.S6, "Replication")
fun.clean.ss.S7.rep <-  merge_samples(fun.clean.ss.S7, "Replication")
fun.clean.ss.S8.rep <-  merge_samples(fun.clean.ss.S8, "Replication")
fun.clean.ss.S9.rep <-  merge_samples(fun.clean.ss.S9, "Replication")

fun.clean.ss.R.rep <-  merge_samples(fun.clean.ss.R, "Replication")
fun.clean.ss.RS.rep <-  merge_samples(fun.clean.ss.RS, "Replication")
fun.clean.ss.BS.rep <-  merge_samples(fun.clean.ss.BS, "Replication")

fun.clean.ss.G.rep <-  merge_samples(fun.clean.ss.G, "Replication")


fun.clean.ss.L1.rep.rel <- microbiome::transform(fun.clean.ss.L1.rep, "compositional")
fun.clean.ss.L2.rep.rel <- microbiome::transform(fun.clean.ss.L2.rep, "compositional")
fun.clean.ss.L3.rep.rel <- microbiome::transform(fun.clean.ss.L3.rep, "compositional")
fun.clean.ss.FL.rep.rel <- microbiome::transform(fun.clean.ss.FL.rep, "compositional")

fun.clean.ss.S1.rep.rel <- microbiome::transform(fun.clean.ss.S1.rep, "compositional")
fun.clean.ss.S2.rep.rel <- microbiome::transform(fun.clean.ss.S2.rep, "compositional")
fun.clean.ss.S3.rep.rel <- microbiome::transform(fun.clean.ss.S3.rep, "compositional")
fun.clean.ss.S4.rep.rel <- microbiome::transform(fun.clean.ss.S4.rep, "compositional")
fun.clean.ss.S5.rep.rel <- microbiome::transform(fun.clean.ss.S5.rep, "compositional")
fun.clean.ss.S6.rep.rel <- microbiome::transform(fun.clean.ss.S6.rep, "compositional")
fun.clean.ss.S7.rep.rel <- microbiome::transform(fun.clean.ss.S7.rep, "compositional")
fun.clean.ss.S8.rep.rel <- microbiome::transform(fun.clean.ss.S8.rep, "compositional")
fun.clean.ss.S9.rep.rel <- microbiome::transform(fun.clean.ss.S9.rep, "compositional")

fun.clean.ss.R.rep.rel <- microbiome::transform(fun.clean.ss.R.rep, "compositional")
fun.clean.ss.RS.rep.rel <- microbiome::transform(fun.clean.ss.RS.rep, "compositional")
fun.clean.ss.BS.rep.rel <- microbiome::transform(fun.clean.ss.BS.rep, "compositional")
fun.clean.ss.G.rep.rel <- microbiome::transform(fun.clean.ss.G.rep, "compositional")


#Core microbiota analysis
#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds.

### whole body core
### Whole core regardless of location and year
bac.clean.ss.f.rel <- microbiome::transform(bac.clean.ss.f, "compositional")
bac.clean.rel.rep <-  merge_samples(bac.clean.ss.f.rel, "Replication")

core.bac.whole<- core_members(bac.clean.rel.rep, detection = 0, prevalence = 80/100)
core.bac.whole

###Endosphere
bac.clean.ss.endo <- subset_samples(bac.clean.ss.f, Compartment %in% c("Seed","Leaf","Stem","Root"))

bac.clean.ss.endo <- phyloseq::filter_taxa(bac.clean.ss.endo, function(x) sum(x) != 0, TRUE)
bac.clean.ss.endo.rep <-  merge_samples(bac.clean.ss.endo, "Replication")
bac.clean.ss.endo.rep.rel <- microbiome::transform(bac.clean.ss.endo.rep, "compositional")

core.bac.endo <- core_members(bac.clean.ss.endo.rep.rel, detection = 0, prevalence = 80/100)
core.bac.endo

###Above
bac.clean.ss.above <- subset_samples(bac.clean.ss.f, Compartment %in% c("Seed","Leaf","Stem"))

bac.clean.ss.above <- phyloseq::filter_taxa(bac.clean.ss.above, function(x) sum(x) != 0, TRUE)
bac.clean.ss.above.rep <-  merge_samples(bac.clean.ss.above, "Replication2")
bac.clean.ss.above.rep.rel <- microbiome::transform(bac.clean.ss.above.rep, "compositional")

core.bac.above <- core_members(bac.clean.ss.above.rep.rel, detection = 0, prevalence = 90/100)
core.bac.above

####Below
bac.clean.ss.below <- subset_samples(bac.clean.ss.f, Compartment %in% c("Root","Soil","Bulk_soil","Rhizosphere"))

bac.clean.ss.below <- phyloseq::filter_taxa(bac.clean.ss.below, function(x) sum(x) != 0, TRUE)
bac.clean.ss.below.rep <-  merge_samples(bac.clean.ss.below, "Replication2")
bac.clean.ss.below.rep.rel <- microbiome::transform(bac.clean.ss.below.rep, "compositional")

core.bac.below <- core_members(bac.clean.ss.below.rep.rel, detection = 0, prevalence = 90/100)
core.bac.below



###Without Chuncheon samples
### Whole core regardless of location and year
bac.clean.ss.f.rel <- microbiome::transform(bac.clean.ss.f, "compositional")
bac.clean.rel.rep <-  merge_samples(bac.clean.ss.f.rel, "Replication")

core.bac.whole<- core_members(bac.clean.rel.rep, detection = 0, prevalence = 80/100)
core.bac.whole

###Endosphere
bac.clean.ss.endo <- subset_samples(bac.clean.ss.f, Location == "Suwon"&Compartment %in% c("Seed","Leaf","Stem","Root"))

bac.clean.ss.endo <- phyloseq::filter_taxa(bac.clean.ss.endo, function(x) sum(x) != 0, TRUE)
bac.clean.ss.endo.rep <-  merge_samples(bac.clean.ss.endo, "Replication")
bac.clean.ss.endo.rep.rel <- microbiome::transform(bac.clean.ss.endo.rep, "compositional")

core.bac.endo <- core_members(bac.clean.ss.endo.rep.rel, detection = 0, prevalence = 90/100)
core.bac.endo

###Above
bac.clean.ss.above <- subset_samples(bac.clean.ss.f, Location == "Suwon"&Compartment %in% c("Seed","Leaf","Stem"))

bac.clean.ss.above <- phyloseq::filter_taxa(bac.clean.ss.above, function(x) sum(x) != 0, TRUE)
bac.clean.ss.above.rep <-  merge_samples(bac.clean.ss.above, "Replication2")
bac.clean.ss.above.rep.rel <- microbiome::transform(bac.clean.ss.above.rep, "compositional")

core.bac.above <- core_members(bac.clean.ss.above.rep.rel, detection = 0, prevalence = 95/100)
core.bac.above

####Below
bac.clean.ss.below <- subset_samples(bac.clean.ss.f, Location == "Suwon"&Compartment %in% c("Root","Soil","Bulk_soil","Rhizosphere"))

bac.clean.ss.below <- phyloseq::filter_taxa(bac.clean.ss.below, function(x) sum(x) != 0, TRUE)
bac.clean.ss.below.rep <-  merge_samples(bac.clean.ss.below, "Replication2")
bac.clean.ss.below.rep.rel <- microbiome::transform(bac.clean.ss.below.rep, "compositional")

core.bac.below <- core_members(bac.clean.ss.below.rep.rel, detection = 0, prevalence = 90/100)
core.bac.below



bac.clean.ss.endo <- subset_samples(bac.clean.ss.f, Year == "year_2018")

bac.clean.ss.endo <- phyloseq::filter_taxa(bac.clean.ss.endo, function(x) sum(x) != 0, TRUE)
bac.clean.ss.endo.rep <-  merge_samples(bac.clean.ss.endo, "Replication")
bac.clean.ss.endo.rep.rel <- microbiome::transform(bac.clean.ss.endo.rep, "compositional")

core.bac.endo <- core_members(bac.clean.ss.endo.rep.rel, detection = 0, prevalence = 75/100)
core.bac.endo

fun.clean.ss.endo <- subset_samples(fun.clean.ss.f, Year == "year_2018")

fun.clean.ss.endo <- phyloseq::filter_taxa(fun.clean.ss.endo, function(x) sum(x) != 0, TRUE)
fun.clean.ss.endo.rep <-  merge_samples(fun.clean.ss.endo, "Replication")
fun.clean.ss.endo.rep.rel <- microbiome::transform(fun.clean.ss.endo.rep, "compositional")

core.fun.endo <- core_members(fun.clean.ss.endo.rep.rel, detection = 0, prevalence = 75/100)
core.fun.endo

bac.list.core.endo <- subset(bac.list, OTU %in% core.bac.endo)
write.csv(bac.list.core.endo,"whole core_bac.csv")

fun.list.core.endo <- subset(fun.list, OTU %in% core.fun.endo)
write.csv(fun.list.core.endo,"whole core_fun.csv")

###Endophyte-specific core
bac.clean.ss.endo <- subset_samples(bac.clean.ss.f, Compartment %in% c("Leaf","Stem","Root","Seed"))
bac.clean.ss.endo <- subset_samples(bac.clean.ss.endo, Year == "year_2018")

bac.clean.ss.endo <- phyloseq::filter_taxa(bac.clean.ss.endo, function(x) sum(x) != 0, TRUE)
bac.clean.ss.endo.rep <-  merge_samples(bac.clean.ss.endo, "Replication")
bac.clean.ss.endo.rep.rel <- microbiome::transform(bac.clean.ss.endo.rep, "compositional")

core.bac.endo <- core_members(bac.clean.ss.endo.rep.rel, detection = 0, prevalence = 75/100)
core.bac.endo

fun.clean.ss.endo <- subset_samples(fun.clean.ss.f, Compartment %in% c("Leaf","Stem","Root","Grain"))
fun.clean.ss.endo <- subset_samples(fun.clean.ss.endo, Year == "year_2018")

fun.clean.ss.endo <- phyloseq::filter_taxa(fun.clean.ss.endo, function(x) sum(x) != 0, TRUE)
fun.clean.ss.endo.rep <-  merge_samples(fun.clean.ss.endo, "Replication")
fun.clean.ss.endo.rep.rel <- microbiome::transform(fun.clean.ss.endo.rep, "compositional")

core.fun.endo <- core_members(fun.clean.ss.endo.rep.rel, detection = 0, prevalence = 75/100)
core.fun.endo

bac.list.core.endo <- subset(bac.list, OTU %in% core.bac.endo)
write.csv(bac.list.core.endo,"Endophytic core_bac.csv")

fun.list.core.endo <- subset(fun.list, OTU %in% core.fun.endo)
write.csv(fun.list.core.endo,"Endophytic core_fun.csv")


### Compartment-specific core
core.bac.90.L1 <- core_members(bac.clean.ss.L1.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.L1

core.bac.90.L2 <- core_members(bac.clean.ss.L2.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.L2

core.bac.90.L3 <- core_members(bac.clean.ss.L3.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.L3

core.bac.90.FL <- core_members(bac.clean.ss.FL.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.FL

core.bac.90.S1 <- core_members(bac.clean.ss.S1.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S1

core.bac.90.S2 <- core_members(bac.clean.ss.S2.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S2

core.bac.90.S3 <- core_members(bac.clean.ss.S3.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S3

core.bac.90.S4 <- core_members(bac.clean.ss.S4.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S4

core.bac.90.S5 <- core_members(bac.clean.ss.S5.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S5

core.bac.90.S6 <- core_members(bac.clean.ss.S6.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S6

core.bac.90.S7 <- core_members(bac.clean.ss.S7.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S7

core.bac.90.S8 <- core_members(bac.clean.ss.S8.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S8

core.bac.90.S9 <- core_members(bac.clean.ss.S9.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.S9

core.bac.90.R <- core_members(bac.clean.ss.R.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.R

core.bac.90.RS <- core_members(bac.clean.ss.RS.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.RS

core.bac.90.BS <- core_members(bac.clean.ss.BS.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.BS

core.bac.90.G <- core_members(bac.clean.ss.G.rep.rel, detection = 0, prevalence = 90/100)
core.bac.90.G

##fungal core
core.fun.90.L1 <- core_members(fun.clean.ss.L1.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.L1

core.fun.90.L2 <- core_members(fun.clean.ss.L2.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.L2

core.fun.90.L3 <- core_members(fun.clean.ss.L3.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.L3

core.fun.90.FL <- core_members(fun.clean.ss.FL.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.FL

core.fun.90.S1 <- core_members(fun.clean.ss.S1.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S1

core.fun.90.S2 <- core_members(fun.clean.ss.S2.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S2

core.fun.90.S3 <- core_members(fun.clean.ss.S3.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S3

core.fun.90.S4 <- core_members(fun.clean.ss.S4.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S4

core.fun.90.S5 <- core_members(fun.clean.ss.S5.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S5

core.fun.90.S6 <- core_members(fun.clean.ss.S6.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S6

core.fun.90.S7 <- core_members(fun.clean.ss.S7.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S7

core.fun.90.S8 <- core_members(fun.clean.ss.S8.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S8

core.fun.90.S9 <- core_members(fun.clean.ss.S9.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.S9

core.fun.90.R <- core_members(fun.clean.ss.R.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.R

core.fun.90.RS <- core_members(fun.clean.ss.RS.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.RS

core.fun.90.BS <- core_members(fun.clean.ss.BS.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.BS

core.fun.90.G <- core_members(fun.clean.ss.G.rep.rel, detection = 0, prevalence = 90/100)
core.fun.90.G

### leaf and stem core
##Bacteria
bac.clean.ss.L <- subset_samples(bac.clean.ss.f, Microhabitat %in%c("L1","L2","L3","FL"))
bac.clean.ss.L <- phyloseq::filter_taxa(bac.clean.ss.L, function(x) sum(x) != 0, TRUE)
bac.clean.ss.L.rep <-  merge_samples(bac.clean.ss.L, "Replication2")
bac.clean.ss.L.rep.rel <- microbiome::transform(bac.clean.ss.L.rep, "compositional")

core.bac.L <- core_members(bac.clean.ss.L.rep.rel, detection = 0, prevalence = 90/100)
core.bac.L

intersect(core.bac.L, persistent.bac)


bac.clean.ss.S <- subset_samples(bac.clean.ss.f, Microhabitat %in%c("S1","S2","S3","S4","S5","S6","S7","S8","S9"))
bac.clean.ss.S <- phyloseq::filter_taxa(bac.clean.ss.S, function(x) sum(x) != 0, TRUE)
bac.clean.ss.S.rep <-  merge_samples(bac.clean.ss.S, "Replication2")
bac.clean.ss.S.rep.rel <- microbiome::transform(bac.clean.ss.S.rep, "compositional")

core.bac.S <- core_members(bac.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
core.bac.S

intersect(core.bac.S, persistent.bac)

stem.core.inherit<-Reduce(intersect, list(core.bac.S,persistent.bac))
leaf.core.inherit<-Reduce(intersect, list(core.bac.L,persistent.bac))
leaf.stem.core.inherit <-Reduce(intersect, list(core.bac.S,core.bac.L,persistent.bac))
subset(bac.list, OTU%in%leaf.stem.core.inherit)
##Fungi
# meta.fun<-sample_data(fun.clean.ss.f)
# write.csv(meta.fun,"meta.fun_needtoedit.csv")
meta.fun <- read.csv('meta.fun_needtoedit.csv')
rownames(meta.fun) <- meta.fun$SampleID
sample_data(fun.clean.ss.f) <- sample_data(meta.fun)

fun.clean.ss.L <- subset_samples(fun.clean.ss.f, Microhabitat %in%c("L1","L2","L3","FL"))
fun.clean.ss.L <- phyloseq::filter_taxa(fun.clean.ss.L, function(x) sum(x) != 0, TRUE)
fun.clean.ss.L.rep <-  merge_samples(fun.clean.ss.L, "Replication2")
fun.clean.ss.L.rep.rel <- microbiome::transform(fun.clean.ss.L.rep, "compositional")

core.fun.L <- core_members(fun.clean.ss.L.rep.rel, detection = 0, prevalence = 90/100)
core.fun.L

intersect(core.fun.L, persistent.fun)
intersect(core.fun.L, core.fun.90.G)

fun.clean.ss.S <- subset_samples(fun.clean.ss.f, Microhabitat %in%c("S1","S2","S3","S4","S5","S6","S7","S8","S9"))
fun.clean.ss.S <- phyloseq::filter_taxa(fun.clean.ss.S, function(x) sum(x) != 0, TRUE)
fun.clean.ss.S.rep <-  merge_samples(fun.clean.ss.S, "Replication2")
fun.clean.ss.S.rep.rel <- microbiome::transform(fun.clean.ss.S.rep, "compositional")

core.fun.S <- core_members(fun.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
core.fun.S

stem.core.inherit.f<-Reduce(intersect, list(core.fun.S,persistent.fun))
leaf.core.inherit.f<-Reduce(intersect, list(core.fun.L,persistent.fun))
leaf.stem.core.inherit.f <-Reduce(intersect, list(core.fun.S,core.fun.L,persistent.fun,core.fun.90.G))
subset(fun.list, OTU%in%leaf.stem.core.inherit.f)

## Table
write.csv(core.fun.L, "Fungi_leaf core.csv")
write.csv(core.fun.S, "Fungi_stem core.csv")
write.csv(core.fun.90.G, "Fungi_seed core.csv")
write.csv(persistent.fun, "Fungi_inhertied OTU.csv")

write.csv(core.bac.L, "Bacteria_leaf core.csv")
write.csv(core.bac.S, "Bacteria_stem core.csv")
write.csv(core.bac.90.G, "Bacteria_seed core.csv")
write.csv(persistent.bac, "Bacteria_inhertied OTU.csv")


### Assigning groups
tax.tab.inhertied.bac <- subset(bac.list, OTU%in%persistent.bac)
tax.tab.inhertied.bac$Group <- 0
tax.tab.inhertied.bac$Group[tax.tab.inhertied.bac$OTU%in%Reduce(intersect, list(core.bac.S,core.bac.L,persistent.bac,core.bac.90.G))] <- "Inherited above core"
tax.tab.inhertied.bac$Group[tax.tab.inhertied.bac$OTU == "EF182718.1.1422"] <- "Inherited seed and leaf core"
tax.tab.inhertied.bac$Group[tax.tab.inhertied.bac$OTU %in% c("KJ184878.1.1455","EU239123.1.1223","MLFS01000124.50.1600", "EU705754.1.1215","KU725926.1.1331")] <- "Inherited seed core"
tax.tab.inhertied.bac$Group[tax.tab.inhertied.bac$OTU %in% c("KY393355.1.1390", "KJ606801.1.1322","EF098942.1.1345","JF224904.1.1305","HM274363.1.1305","JQ977010.1.1397")] <- "Inherited stem core"
tax.tab.inhertied.bac$Group[tax.tab.inhertied.bac$Group == "0"] <- "Inherited noncore"


write.csv(tax.tab.inhertied.bac,"tax.tab.inhertied.bac.csv")

tax.tab.inhertied.fun <- subset(fun.list, OTU%in%persistent.fun)
tax.tab.inhertied.fun$Group <- 0
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU%in%Reduce(intersect, list(core.fun.S,core.fun.L,persistent.fun,core.fun.90.G))] <- "Inherited above core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU == "SH194976.07FU_AY015439_refs"] <- "Inherited leaf core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU %in% c("SH640759.07FU_KX515816_reps","SH523895.07FU_KT581852_reps","SH177636.07FU_HM535391_reps")] <- "Inherited seed core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU %in% c("SH180705.07FU_FJ372389_reps", "aa8b76ce061ef5e33efbaa07228b3c74","SH182913.07FU_KC916678_reps","SH182985.07FU_AF191549_reps")] <- "Inherited stem core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU %in% c( "SH198996.07FU_JF740262_refs","SH180757.07FU_KU574712_reps")] <- "Inherited leaf and stem core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$OTU %in% c( "SH187563.07FU_JN192379_refs","b5747d266e5e82ef6c7308be1bd3e5ce")] <- "Inherited seed and stem core"
tax.tab.inhertied.fun$Group[tax.tab.inhertied.fun$Group == "0"] <- "Inherited noncore"


write.csv(tax.tab.inhertied.fun,"tax.tab.inhertied.fun.csv")


### Inherited OTU
persistent.bac
persistent.fun

L1.common<-intersect(core.bac.90.L1,persistent.bac)
L2.common<-intersect(core.bac.90.L2,persistent.bac)
L3.common<-intersect(core.bac.90.L3,persistent.bac)
FL.common<-intersect(core.bac.90.FL,persistent.bac)

S1.common<-intersect(core.bac.90.S1,persistent.bac)
S2.common<-intersect(core.bac.90.S2,persistent.bac)
S3.common<-intersect(core.bac.90.S3,persistent.bac)
S4.common<-intersect(core.bac.90.S4,persistent.bac)
S5.common<-intersect(core.bac.90.S5,persistent.bac)
S6.common<-intersect(core.bac.90.S6,persistent.bac)
S7.common<-intersect(core.bac.90.S7,persistent.bac)
S8.common<-intersect(core.bac.90.S8,persistent.bac)
S9.common<-intersect(core.bac.90.S9,persistent.bac)

R.common<-intersect(core.bac.90.R,persistent.bac)
RS.common<-intersect(core.bac.90.RS,persistent.bac)
BS.common<-intersect(core.bac.90.BS,persistent.bac)


##Fungi
L1.common.f<-intersect(core.fun.90.L1,persistent.fun)
L2.common.f<-intersect(core.fun.90.L2,persistent.fun)
L3.common.f<-intersect(core.fun.90.L3,persistent.fun)
FL.common.f<-intersect(core.fun.90.FL,persistent.fun)

S1.common.f<-intersect(core.fun.90.S1,persistent.fun)
S2.common.f<-intersect(core.fun.90.S2,persistent.fun)
S3.common.f<-intersect(core.fun.90.S3,persistent.fun)
S4.common.f<-intersect(core.fun.90.S4,persistent.fun)
S5.common.f<-intersect(core.fun.90.S5,persistent.fun)
S6.common.f<-intersect(core.fun.90.S6,persistent.fun)
S7.common.f<-intersect(core.fun.90.S7,persistent.fun)
S8.common.f<-intersect(core.fun.90.S8,persistent.fun)
S9.common.f<-intersect(core.fun.90.S9,persistent.fun)

R.common.f<-intersect(core.fun.90.R,persistent.fun)
RS.common.f<-intersect(core.fun.90.RS,persistent.fun)
BS.common.f<-intersect(core.fun.90.BS,persistent.fun)

### G core
L1.G<-intersect(core.bac.90.L1,core.bac.90.G)
L2.G<-intersect(core.bac.90.L2,core.bac.90.G)
L3.G<-intersect(core.bac.90.L3,core.bac.90.G)
FL.G<-intersect(core.bac.90.FL,core.bac.90.G)

S1.G<-intersect(core.bac.90.S1,core.bac.90.G)
S2.G<-intersect(core.bac.90.S2,core.bac.90.G)
S3.G<-intersect(core.bac.90.S3,core.bac.90.G)
S4.G<-intersect(core.bac.90.S4,core.bac.90.G)
S5.G<-intersect(core.bac.90.S5,core.bac.90.G)
S6.G<-intersect(core.bac.90.S6,core.bac.90.G)
S7.G<-intersect(core.bac.90.S7,core.bac.90.G)
S8.G<-intersect(core.bac.90.S8,core.bac.90.G)
S9.G<-intersect(core.bac.90.S9,core.bac.90.G)

R.G<-intersect(core.bac.90.R,core.bac.90.G)
RS.G<-intersect(core.bac.90.RS,core.bac.90.G)
BS.G<-intersect(core.bac.90.BS,core.bac.90.G)

intersect(persistent.bac, core.bac.90.G)
intersect(persistent.fun, core.fun.90.G)


##Fungi
L1.G.f<-intersect(core.fun.90.L1,core.fun.90.G)
L2.G.f<-intersect(core.fun.90.L2,core.fun.90.G)
L3.G.f<-intersect(core.fun.90.L3,core.fun.90.G)
FL.G.f<-intersect(core.fun.90.FL,core.fun.90.G)

S1.G.f<-intersect(core.fun.90.S1,core.fun.90.G)
S2.G.f<-intersect(core.fun.90.S2,core.fun.90.G)
S3.G.f<-intersect(core.fun.90.S3,core.fun.90.G)
S4.G.f<-intersect(core.fun.90.S4,core.fun.90.G)
S5.G.f<-intersect(core.fun.90.S5,core.fun.90.G)
S6.G.f<-intersect(core.fun.90.S6,core.fun.90.G)
S7.G.f<-intersect(core.fun.90.S7,core.fun.90.G)
S8.G.f<-intersect(core.fun.90.S8,core.fun.90.G)
S9.G.f<-intersect(core.fun.90.S9,core.fun.90.G)

R.G.f<-intersect(core.fun.90.R,core.fun.90.G)
RS.G.f<-intersect(core.fun.90.RS,core.fun.90.G)
BS.G.f<-intersect(core.fun.90.BS,core.fun.90.G)


intersect(S9.G.f, S9.common.f)


Reduce(intersect, list(S1.common, S2.common, S3.common, S4.common, S5.common, S6.common,S7.common,S8.common,S9.common))



### Relative abundance of inherited core OTUs (heat map)

fun.clean.ss.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)
fun.clean.ss.18.rep <-  merge_samples(fun.clean.ss.18, "Replication2")
fun.clean.ss.18.rep.rel <- microbiome::transform(fun.clean.ss.18.rep, "compositional")

bac.clean.ss.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)
bac.clean.ss.18.rep <-  merge_samples(bac.clean.ss.18, "Replication2")
bac.clean.ss.18.rep.rel <- microbiome::transform(bac.clean.ss.18.rep, "compositional")


otu.bac.18<-otu_table(bac.clean.ss.18.rep.rel)
otu.fun.18<-otu_table(fun.clean.ss.18.rep.rel)
otu.fun.18<-t(otu.fun.18)
otu.bac.18<-t(otu.bac.18)

df.otu.bac.18 <- data.frame(otu.bac.18)
df.otu.fun.18 <- data.frame(otu.fun.18)

df.otu.bac.18.inherit<-subset(df.otu.bac.18, rownames(df.otu.bac.18)%in%persistent.bac)
df.otu.fun.18.inherit<-subset(df.otu.fun.18, rownames(df.otu.fun.18)%in%persistent.fun)

write.xlsx(df.otu.bac.18.inherit,"df.otu.bac.18.inherit.xlsx")
write.xlsx(df.otu.fun.18.inherit,"df.otu.fun.18.inherit.xlsx")


otu.tab.fun <- read.xlsx("df.otu.fun.18.inherit.xlsx",1)
otu.tab.fun<-merge(otu.tab.fun, tax.tab.inhertied.fun, by = c('OTU' = "OTU"))
write.xlsx(otu.tab.fun,"df.otu.fun.18.inherit.xlsx")



otu.tab.bac <- read.xlsx("df.otu.bac.18.inherit.xlsx",1)
otu.tab.bac<-merge(otu.tab.bac, tax.tab.inhertied.bac, by = c('OTU' = "OTU"))
write.xlsx(otu.tab.bac,"df.otu.bac.18.inherit.xlsx")

### Ridge plot
##
bac.ridge.input <- read.xlsx("df.otu.bac.18.inherit_input.xlsx",1)

melt.bac.ridge<-melt(bac.ridge.input)
melt.bac.ridge<-subset(melt.bac.ridge, variable != 'total')
melt.bac.ridge$AbbrOTU <- factor(melt.bac.ridge$AbbrOTU, levels = rev(bac.ridge.input$AbbrOTU))

fun.ridge.input <- read.xlsx("df.otu.fun.18.inherit.xlsx",1)
melt.fun.ridge<-melt(fun.ridge.input)
melt.fun.ridge<-subset(melt.fun.ridge, variable != 'total')
melt.fun.ridge$AbbrOTU <- factor(melt.fun.ridge$AbbrOTU, levels = rev(fun.ridge.input$AbbrOTU))


library(ggridges)
  ggplot(melt.fun.ridge, aes(x=variable, y=AbbrOTU, height=value, group = AbbrOTU, fill=Group)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=5, color='white',size=0.5)+
    xlab('\n Compart')+
    ylab("Relative abundance \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="bottom", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size="none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    #scale_fill_manual(labels = c('Non-fer','Fer','ns'), values = c("Non-fer"= "Black", 'Fer'='Red','ns'='gray50'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
    
    #geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)
  

### transient OTU and cores
  ### leaf and stem core
  ##Bacteria
  bac.clean.ss.L <- subset_samples(bac.clean.ss.f, Microhabitat %in%c("L1","L2","L3","FL"))
  bac.clean.ss.L <- phyloseq::filter_taxa(bac.clean.ss.L, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.L.rep <-  merge_samples(bac.clean.ss.L, "Replication2")
  bac.clean.ss.L.rep.rel <- microbiome::transform(bac.clean.ss.L.rep, "compositional")
  
  core.bac.L <- core_members(bac.clean.ss.L.rep.rel, detection = 0, prevalence = 90/100)
  core.bac.L
  
  intersect(core.bac.L, transient.bac)
  
  
  bac.clean.ss.S <- subset_samples(bac.clean.ss.f, Microhabitat %in%c("S1","S2","S3","S4","S5","S6","S7","S8","S9"))
  bac.clean.ss.S <- phyloseq::filter_taxa(bac.clean.ss.S, function(x) sum(x) != 0, TRUE)
  bac.clean.ss.S.rep <-  merge_samples(bac.clean.ss.S, "Replication2")
  bac.clean.ss.S.rep.rel <- microbiome::transform(bac.clean.ss.S.rep, "compositional")
  
  core.bac.S <- core_members(bac.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
  core.bac.S
  
  intersect(core.bac.S, transient.bac)
  
  stem.core.transient<-Reduce(intersect, list(core.bac.S,transient.bac))
  leaf.core.transient<-Reduce(intersect, list(core.bac.L,transient.bac))
  leaf.stem.core.transient <-Reduce(intersect, list(core.bac.S,core.bac.L,transient.bac))
  subset(bac.list, OTU%in%leaf.stem.core.transient)
  ##Fungi
  # meta.fun<-sample_data(fun.clean.ss.f)
  # write.csv(meta.fun,"meta.fun_needtoedit.csv")
  meta.fun <- read.csv('meta.fun_needtoedit.csv')
  rownames(meta.fun) <- meta.fun$SampleID
  sample_data(fun.clean.ss.f) <- sample_data(meta.fun)
  
  fun.clean.ss.L <- subset_samples(fun.clean.ss.f, Microhabitat %in%c("L1","L2","L3","FL"))
  fun.clean.ss.L <- phyloseq::filter_taxa(fun.clean.ss.L, function(x) sum(x) != 0, TRUE)
  fun.clean.ss.L.rep <-  merge_samples(fun.clean.ss.L, "Replication2")
  fun.clean.ss.L.rep.rel <- microbiome::transform(fun.clean.ss.L.rep, "compositional")
  
  core.fun.L <- core_members(fun.clean.ss.L.rep.rel, detection = 0, prevalence = 90/100)
  core.fun.L
  
  intersect(core.fun.L, transient.fun)

  fun.clean.ss.S <- subset_samples(fun.clean.ss.f, Microhabitat %in%c("S1","S2","S3","S4","S5","S6","S7","S8","S9"))
  fun.clean.ss.S <- phyloseq::filter_taxa(fun.clean.ss.S, function(x) sum(x) != 0, TRUE)
  fun.clean.ss.S.rep <-  merge_samples(fun.clean.ss.S, "Replication2")
  fun.clean.ss.S.rep.rel <- microbiome::transform(fun.clean.ss.S.rep, "compositional")
  
  core.fun.S <- core_members(fun.clean.ss.S.rep.rel, detection = 0, prevalence = 90/100)
  core.fun.S
  
  stem.core.transient.f<-Reduce(intersect, list(core.fun.S,transient.fun))
  leaf.core.transient.f<-Reduce(intersect, list(core.fun.L,transient.fun))
  leaf.stem.core.transient.f <-Reduce(intersect, list(core.fun.S,core.fun.L,transient.fun,core.fun.90.G))
  subset(fun.list, OTU%in%leaf.stem.core.transient.f)
  
  ## Table
  write.csv(core.fun.L, "Fungi_leaf core.csv")
  write.csv(core.fun.S, "Fungi_stem core.csv")
  write.csv(core.fun.90.G, "Fungi_seed core.csv")
  write.csv(transient.fun, "Fungi_transient OTU.csv")
  
  write.csv(core.bac.L, "Bacteria_leaf core.csv")
  write.csv(core.bac.S, "Bacteria_stem core.csv")
  write.csv(core.bac.90.G, "Bacteria_seed core.csv")
  write.csv(transient.bac, "Bacteria_transient OTU.csv")
  
  
  ### Assigning groups
  tax.tab.transient.bac <- subset(bac.list, OTU%in%transient.bac)
  tax.tab.transient.bac$Group <- 0
  tax.tab.transient.bac$Group[tax.tab.transient.bac$OTU %in% c("KU513667.1.1232","GASZ01001652.185.1608")] <- "transient seed and stem core"
  tax.tab.transient.bac$Group[tax.tab.transient.bac$OTU %in% c("KT803965.1.1472")] <- "transient seed core"
  tax.tab.transient.bac$Group[tax.tab.transient.bac$OTU %in% c("KJ200413.1.1428", "JQ311912.1.1466","KJ601761.1.1362","KM253192.1.1336")] <- "transient stem core"
  tax.tab.transient.bac$Group[tax.tab.transient.bac$Group == "0"] <- "transient noncore"
  
  
  write.csv(tax.tab.transient.bac,"tax.tab.transient.bac.csv")
  
  tax.tab.transient.fun <- subset(fun.list, OTU%in%transient.fun)
  tax.tab.transient.fun$Group <- 0
  tax.tab.transient.fun$Group[tax.tab.transient.fun$OTU%in%Reduce(intersect, list(core.fun.S,core.fun.L,transient.fun,core.fun.90.G))] <- "transient above core"
  tax.tab.transient.fun$Group[tax.tab.transient.fun$OTU %in% c("SH220700.07FU_AB586992_refs")] <- "transient seed core"
  tax.tab.transient.fun$Group[tax.tab.transient.fun$OTU %in% c("SH188231.07FU_KM246160_reps","SH212071.07FU_KM266105_reps","SH191159.07FU_DQ249199_refs",
                                                               "c1116d32a27f3ba9f45f45ae8389aee9","SH214135.07FU_JQ666320_reps","b1d5037a16626858eaf9775237f972a3",
                                                               "210067274460ee2a076d9283e14aa737","SH206475.07FU_GU721466_reps","SH641103.07FU_KU534864_reps",
                                                               "0cc51c048cc5c73cdf8acba7586e88a5","SH212843.07FU_GU721983_reps","SH376403.07FU_KT693730_refs",
                                                               "SH209530.07FU_AF294699_refs", "SH198028.07FU_FJ265959_reps","SH344881.07FU_AY373880_refs",
                                                               "SH194047.07FU_KM032316_reps", "SH527185.07FU_KJ156314_reps")] <- "transient stem core"
  tax.tab.transient.fun$Group[tax.tab.transient.fun$OTU %in% c( "SH199891.07FU_HQ331085_reps","SH640247.07FU_KY102347_reps",
                                                                "SH194039.07FU_KT959323_reps", "SH215392.07FU_GU721438_reps",
                                                                "SH221540.07FU_AF145586_refs")] <- "transient leaf and stem core"
  tax.tab.transient.fun$Group[tax.tab.transient.fun$OTU %in% c( "SH085061.07FU_FJ196778_refs")] <- "transient seed and stem core"
  tax.tab.transient.fun$Group[tax.tab.transient.fun$Group == "0"] <- "transient noncore"
  
  
  write.csv(tax.tab.transient.fun,"tax.tab.transient.fun.csv")
  
  
  
  df.otu.bac.18.transient<-subset(df.otu.bac.18, rownames(df.otu.bac.18)%in%transient.bac)
  df.otu.fun.18.transient<-subset(df.otu.fun.18, rownames(df.otu.fun.18)%in%transient.fun)
  df.otu.bac.18.transient$OTU <- rownames(df.otu.bac.18.transient)
  df.otu.fun.18.transient$OTU <- rownames(df.otu.fun.18.transient)
  
  
  otu.tab.fun<-merge( df.otu.fun.18.transient, tax.tab.transient.fun, by = c('OTU' = "OTU"))
  write.xlsx(otu.tab.fun,"df.otu.fun.18.transient.xlsx")
  
  otu.tab.bac<-merge( df.otu.bac.18.transient, tax.tab.transient.bac, by = c('OTU' = "OTU"))
  write.xlsx(otu.tab.bac,"df.otu.bac.18.transient.xlsx")

  ### Ridge plot
  ##
  bac.ridge.input <- read.xlsx("df.otu.bac.18.transient_input.xlsx",1)
  
  melt.bac.ridge<-melt(bac.ridge.input)

  melt.bac.ridge$AbbrOTU <- factor(melt.bac.ridge$AbbrOTU, levels = rev(bac.ridge.input$AbbrOTU))
  
  fun.ridge.input <- read.xlsx("df.otu.fun.18.transient_input.xlsx",1)
  melt.fun.ridge<-melt(fun.ridge.input)

  melt.fun.ridge$AbbrOTU <- factor(melt.fun.ridge$AbbrOTU, levels = rev(fun.ridge.input$AbbrOTU))
  
  

  ggplot(melt.fun.ridge, aes(x=variable, y=AbbrOTU, height=value, group = AbbrOTU, fill=Group)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=5, color='white',size=0.5)+
    xlab('\n Compart')+
    ylab("Relative abundance \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="bottom", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size="none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    #scale_fill_manual(labels = c('Non-fer','Fer','ns'), values = c("Non-fer"= "Black", 'Fer'='Red','ns'='gray50'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
  
  #geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)
  
  dev.off()
  
  
  
  ### Supplementary Data 5 & Source data for Venn diagram
  core.bac.90.BS
  core.bac.90.RS
  core.bac.90.R
  core.bac.S
  core.bac.L
  
  core.fun.90.BS
  core.fun.90.RS
  core.fun.90.R
  core.fun.S
  core.fun.L
  
  core.list.bac.BS<-subset(bac.list, OTU%in%core.bac.90.BS)
  core.list.bac.BS$Compartment <- "BS"
  core.list.bac.RS<-subset(bac.list, OTU%in%core.bac.90.RS)
  core.list.bac.RS$Compartment <- "RS"
  core.list.bac.R<-subset(bac.list, OTU%in%core.bac.90.R)
  core.list.bac.R$Compartment <- "R"
  core.list.bac.S<-subset(bac.list, OTU%in%core.bac.S)
  core.list.bac.S$Compartment <- "S"
  core.list.bac.L<-subset(bac.list, OTU%in%core.bac.L)
  core.list.bac.L$Compartment <- "L"
  core.list.bac.G<-subset(bac.list, OTU%in%core.bac.90.G)
  core.list.bac.G$Compartment <- "G"
  
  core.list.bac <- rbind(core.list.bac.BS,core.list.bac.RS,core.list.bac.R,
                         core.list.bac.S,core.list.bac.L,core.list.bac.G)
  write.csv(core.list.bac,"core.list.bac.csv")

    
  core.list.fun.BS<-subset(fun.list, OTU%in%core.fun.90.BS)
  core.list.fun.BS$Compartment <- "BS"
  core.list.fun.RS<-subset(fun.list, OTU%in%core.fun.90.RS)
  core.list.fun.RS$Compartment <- "RS"
  core.list.fun.R<-subset(fun.list, OTU%in%core.fun.90.R)
  core.list.fun.R$Compartment <- "R"
  core.list.fun.S<-subset(fun.list, OTU%in%core.fun.S)
  core.list.fun.S$Compartment <- "S"
  core.list.fun.L<-subset(fun.list, OTU%in%core.fun.L)
  core.list.fun.L$Compartment <- "L"
  core.list.fun.G<-subset(fun.list, OTU%in%core.fun.90.G)
  core.list.fun.G$Compartment <- "G"
  
  core.list.fun <- rbind(core.list.fun.BS,core.list.fun.RS,core.list.fun.R,
                         core.list.fun.S,core.list.fun.L,core.list.fun.G)
  write.csv(core.list.fun,"core.list.fun.csv")
  

  df.example.1<-data.frame(t(otu_table(fun.clean.ss.L.rep)) )
  df.example.2<-data.frame(t(otu_table(fun.clean.ss.L.rep.rel)) )

  write.csv(df.example.1,"df.example.1.csv")  
  write.csv(df.example.2,"df.example.2.csv")  
  