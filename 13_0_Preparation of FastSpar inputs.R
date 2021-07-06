#### Seed network (temporal changes in role of node)
bac.clean.nolog.18.f.seed
sample_names(bac.clean.nolog.18.f.seed)
fun.clean.nolog.18.f.seed


fun.clean.nolog.18.seed <- subset_samples(fun.clean.nolog, Compartment == "Grain" & Year == "year_2018")
fun.clean.nolog.18.seed<- phyloseq::filter_taxa(fun.clean.nolog.18.seed, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18.seed <- subset_samples(bac.clean.nolog, Compartment == "Seed" & Year == "year_2018")
bac.clean.nolog.18.seed<- phyloseq::filter_taxa(bac.clean.nolog.18.seed, function(x) sum(x) != 0, TRUE)


### G0
bac.clean.nolog.18.seed.G0 <- subset_samples(bac.clean.nolog, Replication == "G_0" & Year == "year_2018")
bac.clean.nolog.18.seed.G0<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G0, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18.seed.G76 <- subset_samples(bac.clean.nolog, Replication == "G_76" & Year == "year_2018")
bac.clean.nolog.18.seed.G76<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G76, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18.seed.G90 <- subset_samples(bac.clean.nolog, Replication == "G_90" & Year == "year_2018")
bac.clean.nolog.18.seed.G90<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G90, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18.seed.G106 <- subset_samples(bac.clean.nolog, Replication == "G_106" & Year == "year_2018")
bac.clean.nolog.18.seed.G106<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G106, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18.seed.G120 <- subset_samples(bac.clean.nolog, Replication == "G_120" & Year == "year_2018")
bac.clean.nolog.18.seed.G120<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G120, function(x) sum(x) != 0, TRUE)


bac.clean.nolog.18.seed.G141 <- subset_samples(bac.clean.nolog, Replication == "G_141" & Year == "year_2018")
bac.clean.nolog.18.seed.G141<- phyloseq::filter_taxa(bac.clean.nolog.18.seed.G141, function(x) sum(x) != 0, TRUE)


fun.clean.nolog.18.seed.G0 <- subset_samples(fun.clean.nolog, Replication == "G_0" & Year == "year_2018")
fun.clean.nolog.18.seed.G0<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G0, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18.seed.G76 <- subset_samples(fun.clean.nolog, Replication == "G_76" & Year == "year_2018")
fun.clean.nolog.18.seed.G76<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G76, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18.seed.G90 <- subset_samples(fun.clean.nolog, Replication == "G_90" & Year == "year_2018")
fun.clean.nolog.18.seed.G90<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G90, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18.seed.G106 <- subset_samples(fun.clean.nolog, Replication == "G_106" & Year == "year_2018")
fun.clean.nolog.18.seed.G106<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G106, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18.seed.G120 <- subset_samples(fun.clean.nolog, Replication == "G_120" & Year == "year_2018")
fun.clean.nolog.18.seed.G120<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G120, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18.seed.G141 <- subset_samples(fun.clean.nolog, Replication == "G_141" & Year == "year_2018")
fun.clean.nolog.18.seed.G141<- phyloseq::filter_taxa(fun.clean.nolog.18.seed.G141, function(x) sum(x) != 0, TRUE)


### Extract OTU table (TSV format)
extract_otu_table_bac<-function(physeqfile, keyword){
  otu.table <- otu_table(physeqfile)
  otu.table <- data.frame(otu.table)
  otu.table$OTU <- rownames(otu.table)
  bac.list.id <- bac.list[c("OTU","OTU_id")]
  otu.table.bac <- merge(otu.table, bac.list.id, by = "OTU")
  rownames(otu.table.bac) <- otu.table.bac$OTU_id
  n <- ncol(otu.table.bac)
  otu.table.bac <- otu.table.bac[-c(1,n)]
  write.table(otu.table.bac,paste0('otu_norm_nolog_bac_', keyword,'.tsv'),quote=F, sep = '\t')}

extract_otu_table_fun<-function(physeqfile, keyword){
  otu.table <- otu_table(physeqfile)
  otu.table <- data.frame(otu.table)
  otu.table$OTU <- rownames(otu.table)
  fun.list.id <- fun.list[c("OTU","OTU_id")]
  otu.table.fun <- merge(otu.table, fun.list.id, by = "OTU")
  rownames(otu.table.fun) <- otu.table.fun$OTU_id
  n <- ncol(otu.table.fun)
  otu.table.fun <- otu.table.fun[-c(1,n)]
  write.table(otu.table.fun,paste0('otu_norm_nolog_fun_', keyword,'.tsv'),quote=F, sep = '\t')}

extract_otu_table_bac(bac.clean.nolog.18.seed.G0,"G0")
extract_otu_table_bac(bac.clean.nolog.18.seed.G76,"G76")
extract_otu_table_bac(bac.clean.nolog.18.seed.G90,"G90")
extract_otu_table_bac(bac.clean.nolog.18.seed.G106,"G106")
extract_otu_table_bac(bac.clean.nolog.18.seed.G120,"G120")
extract_otu_table_bac(bac.clean.nolog.18.seed.G141,"G141")

extract_otu_table_fun(fun.clean.nolog.18.seed.G0,"G0")
extract_otu_table_fun(fun.clean.nolog.18.seed.G76,"G76")
extract_otu_table_fun(fun.clean.nolog.18.seed.G90,"G90")
extract_otu_table_fun(fun.clean.nolog.18.seed.G106,"G106")
extract_otu_table_fun(fun.clean.nolog.18.seed.G120,"G120")
extract_otu_table_fun(fun.clean.nolog.18.seed.G141,"G141")

