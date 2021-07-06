###PERMANOVA
#2018 all samples
bac.clean.log.18.woG0 <- subset_samples(bac.clean.log.18, Replication != "G_0")
bac.clean.log.18.woG0 <- phyloseq::filter_taxa(bac.clean.log.18.woG0, function(x) sum(x) != 0, TRUE)

b.otu <- otu_table(bac.clean.log.18.woG0)
b.meta <- sample_data(bac.clean.log.18.woG0)
b.meta <- data.frame(b.meta)

b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment+Microhabitat), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Age), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis



fun.clean.log.18.woG0 <- subset_samples(fun.clean.log.18, Replication != "G_0")
fun.clean.log.18.woG0 <- phyloseq::filter_taxa(fun.clean.log.18.woG0, function(x) sum(x) != 0, TRUE)

f.otu <- otu_table(fun.clean.log.18.woG0)
f.meta <- sample_data(fun.clean.log.18.woG0)
f.meta <- data.frame(f.meta)

f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment+Microhabitat), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Age), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

#2018 belowground
b.otu <- otu_table(bac.clean.log.18.below)
b.meta <- sample_data(bac.clean.log.18.below)
b.meta <- data.frame(b.meta)
b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Age), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis


f.otu <- otu_table(fun.clean.log.18.below)
f.meta <- sample_data(fun.clean.log.18.below)
f.meta <- data.frame(f.meta)
f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Age), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis


#2018 aboveground
bac.clean.log.18.above.woG0 <- subset_samples(bac.clean.log.18.above, Replication != "G_0")
bac.clean.log.18.above.woG0 <- phyloseq::filter_taxa(bac.clean.log.18.above.woG0, function(x) sum(x) != 0, TRUE)

b.otu <- otu_table(bac.clean.log.18.above.woG0)
b.meta <- sample_data(bac.clean.log.18.above.woG0)
b.meta <- data.frame(b.meta)
b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment+Microhabitat), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis
unique(b.meta$Developmental_stage)
abund_table.adonis <- adonis(formula = t(b.otu) ~ (Age), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis


fun.clean.log.18.above.woG0 <- subset_samples(fun.clean.log.18.above, Replication != "G_0")
fun.clean.log.18.above.woG0 <- phyloseq::filter_taxa(fun.clean.log.18.above.woG0, function(x) sum(x) != 0, TRUE)

f.otu <- otu_table(fun.clean.log.18.above.woG0)
f.meta <- sample_data(fun.clean.log.18.above.woG0)
f.meta <- data.frame(f.meta)
f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Air_temp+Humidity+Soil_condition+Developmental_stage+Compartment+Microhabitat), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Age), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis



#2017 all samples
b.otu <- otu_table(bac.clean.log.17)
b.meta <- sample_data(bac.clean.log.17)
b.meta <- data.frame(b.meta)

b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Location+Soil_condition+Air_temp+Humidity+Developmental_stage+Compartment), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Age), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

f.otu <- otu_table(fun.clean.log.17)
f.meta <- sample_data(fun.clean.log.17)
f.meta <- data.frame(f.meta)

f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Location+Soil_condition+Air_temp+Humidity+Developmental_stage+Compartment), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Age), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Location+Age), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

#2017 belowground
bac.clean.log.17.below <- subset_samples(bac.clean.log.17, Compartment %in% c("Root", "Soil"))
bac.clean.log.17.below <- phyloseq::filter_taxa(bac.clean.log.17.below, function(x) sum(x) != 0, TRUE)

b.otu <- otu_table(bac.clean.log.17.below)
b.meta <- sample_data(bac.clean.log.17.below)
b.meta <- data.frame(b.meta)
b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Location+Soil_condition+Age+Air_temp+Humidity+Developmental_stage+Compartment), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis


fun.clean.log.17.below <- subset_samples(fun.clean.log.17, Compartment %in% c("Root", "Soil"))
fun.clean.log.17.below <- phyloseq::filter_taxa(fun.clean.log.17.below, function(x) sum(x) != 0, TRUE)

f.otu <- otu_table(fun.clean.log.17.below)
f.meta <- sample_data(fun.clean.log.17.below)
f.meta <- data.frame(f.meta)
f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Location+Soil_condition+Age+Air_temp+Humidity+Developmental_stage+Compartment), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis





# aboveground
bac.clean.log.above <- subset_samples(bac.clean.log, Compartment %in% c("Leaf", "Stem", "Seed","Grain"))
bac.clean.log.above <- phyloseq::filter_taxa(bac.clean.log.above, function(x) sum(x) != 0, TRUE)

bac.clean.log.above <- subset_samples(bac.clean.log.above, Replication != "G_0")
bac.clean.log.above <- phyloseq::filter_taxa(bac.clean.log.above, function(x) sum(x) != 0, TRUE)


b.otu <- otu_table(bac.clean.log.above)
b.meta <- sample_data(bac.clean.log.above)
b.meta <- data.frame(b.meta)
b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))
b.meta$Compartment[which(b.meta$Compartment == "Grain")] <- "Seed"
unique(b.meta$Compartment)


abund_table.adonis <- adonis(formula = t(b.otu) ~ (Location+Age+Developmental_stage+Compartment), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Soil_condition+Air_temp+Humidity), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis


# belowground
bac.clean.log.below <- subset_samples(bac.clean.log, Compartment %in% c("Soil", "Bulk_soil", "Rhizosphere","Root"))
bac.clean.log.below <- phyloseq::filter_taxa(bac.clean.log.below, function(x) sum(x) != 0, TRUE)


b.otu <- otu_table(bac.clean.log.below)
b.meta <- sample_data(bac.clean.log.below)
b.meta <- data.frame(b.meta)
b.meta$Air_temp <- as.numeric(as.character(b.meta$Air_temp))
b.meta$Humidity <- as.numeric(as.character(b.meta$Humidity))
b.meta$Compartment[which(b.meta$Compartment == "Soil")] <- "Bulk_soil"
unique(b.meta$Compartment)

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Location+Age+Developmental_stage+Compartment), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(b.otu) ~ (Soil_condition+Air_temp+Humidity), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis


# aboveground
fun.clean.log.above <- subset_samples(fun.clean.log, Compartment %in% c("Leaf", "Stem", "Seed","Grain"))
fun.clean.log.above <- phyloseq::filter_taxa(fun.clean.log.above, function(x) sum(x) != 0, TRUE)

fun.clean.log.above <- subset_samples(fun.clean.log.above, Replication != "G_0")
fun.clean.log.above <- phyloseq::filter_taxa(fun.clean.log.above, function(x) sum(x) != 0, TRUE)


f.otu <- otu_table(fun.clean.log.above)
f.meta <- sample_data(fun.clean.log.above)
f.meta <- data.frame(f.meta)
f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))
f.meta$Compartment[which(f.meta$Compartment == "Grain")] <- "Seed"
unique(f.meta$Compartment)
abund_table.adonis <- adonis(formula = t(f.otu) ~ (Location+Age+Developmental_stage+Compartment), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Soil_condition+Air_temp+Humidity), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis


# belowground
fun.clean.log.below <- subset_samples(fun.clean.log, Compartment %in% c("Soil", "Bulk_soil", "Rhizosphere","Root"))
fun.clean.log.below <- phyloseq::filter_taxa(fun.clean.log.below, function(x) sum(x) != 0, TRUE)


f.otu <- otu_table(fun.clean.log.below)
f.meta <- sample_data(fun.clean.log.below)
f.meta <- data.frame(f.meta)
f.meta$Air_temp <- as.numeric(as.character(f.meta$Air_temp))
f.meta$Humidity <- as.numeric(as.character(f.meta$Humidity))
f.meta$Compartment[which(f.meta$Compartment == "Soil")] <- "Bulk_soil"
unique(f.meta$Compartment)

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Location+Age+Developmental_stage+Compartment), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

abund_table.adonis <- adonis(formula = t(f.otu) ~ (Soil_condition+Air_temp+Humidity), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis


