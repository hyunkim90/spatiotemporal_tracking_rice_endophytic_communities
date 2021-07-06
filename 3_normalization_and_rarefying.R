### Normalization

## Let's do normalization with CSS
## phyloseq to metagenomeSeq


#phy.clean? or phy.clean --> let's start from phy.clean

bac.clean.ss

arch.clean.ss

fun.clean.ss2


## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(bac.clean.ss) == 0)
taxa_sums(bac.clean.ss)

bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(bac.clean.ss) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.clean.ss))


## only keep samples over 0 read
(filt.sample <- sample_sums(bac.clean.ss) > 0)
sum(sample_sums(bac.clean.ss) <= 0)  ## 1 sample discarded
bac.clean.ss.f <- prune_samples(filt.sample, bac.clean.ss)
bac.clean.ss.f  ## 979 samples <- 984 samples


sort(sample_sums(bac.clean.ss.f))

# #Archaea
# sum(taxa_sums(arch.clean.ss) == 0)
# taxa_sums(arch.clean.ss)
# 
# arch.clean.ss <- phyloseq::filter_taxa(arch.clean.ss, function(x) sum(x) != 0, TRUE)
# sum(taxa_sums(arch.clean.ss) == 0)
# 
# sort(sample_sums(arch.clean.ss))
# 
# 
# ## only keep samples over 0 read
# (filt.sample <- sample_sums(arch.clean.ss) > 0)
# sum(sample_sums(arch.clean.ss) <= 0)  
# arch.clean.ss.f <- prune_samples(filt.sample, arch.clean.ss)
# arch.clean.ss.f  ## 326 samples <- 1164 samples
# sort(sample_sums(arch.clean.ss.f))


## CODE for CSS normalization using preloaded data
bac.clean.filt <- bac.clean.ss.f
bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog <- bac.clean.filt
otu_table(bac.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log <- bac.clean.filt
otu_table(bac.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


#2017
bac.clean.filt <- bac.clean.ss.17
bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17 <- bac.clean.filt
otu_table(bac.clean.nolog.17) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17 <- bac.clean.filt
otu_table(bac.clean.log.17) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


bac.clean.filt <- bac.clean.ss.18
bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18 <- bac.clean.filt
otu_table(bac.clean.nolog.18) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18 <- bac.clean.filt
otu_table(bac.clean.log.18) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


### Normalization sorted by compartment
## Leaf
bac.clean.filt <- bac.clean.ss.17.leaf

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17.leaf <- bac.clean.filt
otu_table(bac.clean.nolog.17.leaf) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17.leaf <- bac.clean.filt
otu_table(bac.clean.log.17.leaf) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Stem
bac.clean.filt <- bac.clean.ss.17.stem

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17.stem <- bac.clean.filt
otu_table(bac.clean.nolog.17.stem) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17.stem <- bac.clean.filt
otu_table(bac.clean.log.17.stem) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Root
bac.clean.filt <- bac.clean.ss.17.root

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17.root <- bac.clean.filt
otu_table(bac.clean.nolog.17.root) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17.root <- bac.clean.filt
otu_table(bac.clean.log.17.root) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Seed
bac.clean.filt <- bac.clean.ss.17.seed

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17.seed <- bac.clean.filt
otu_table(bac.clean.nolog.17.seed) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17.seed <- bac.clean.filt
otu_table(bac.clean.log.17.seed) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Soil
bac.clean.filt <- bac.clean.ss.17.soil

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.17.soil <- bac.clean.filt
otu_table(bac.clean.nolog.17.soil) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.17.soil <- bac.clean.filt
otu_table(bac.clean.log.17.soil) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

### 2018 Samples

## Leaf
bac.clean.filt <- bac.clean.ss.18.f.leaf

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.leaf <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.leaf) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.leaf <- bac.clean.filt
otu_table(bac.clean.log.18.f.leaf) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Stem
bac.clean.filt <- bac.clean.ss.18.f.stem

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.stem <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.stem) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.stem <- bac.clean.filt
otu_table(bac.clean.log.18.f.stem) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Root
bac.clean.filt <- bac.clean.ss.18.f.root

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.root <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.root) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.root <- bac.clean.filt
otu_table(bac.clean.log.18.f.root) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Seed
bac.clean.filt <- bac.clean.ss.18.f.seed

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.seed <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.seed) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.seed <- bac.clean.filt
otu_table(bac.clean.log.18.f.seed) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Bulk soil
bac.clean.filt <- bac.clean.ss.18.f.BS

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.BS <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.BS) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.BS <- bac.clean.filt
otu_table(bac.clean.log.18.f.BS) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



## Rhizosphere
bac.clean.filt <- bac.clean.ss.18.f.RS

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.18.f.RS <- bac.clean.filt
otu_table(bac.clean.nolog.18.f.RS) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.18.f.RS <- bac.clean.filt
otu_table(bac.clean.log.18.f.RS) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



### fungal community ###
## only keep samples over 0 read
(filt.sample <- sample_sums(fun.clean.ss2) > 10)
sum(sample_sums(fun.clean.ss2) <= 10)  ## 1 sample discarded
fun.clean.ss2.f <- prune_samples(filt.sample, fun.clean.ss2)
fun.clean.ss2.f  ## 979 samples <- 984 samples


sort(sample_sums(fun.clean.ss2.f))

## CODE for CSS normalization using preloaded data
fun.clean.filt <- fun.clean.ss2.f
fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog <- fun.clean.filt
otu_table(fun.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log <- fun.clean.filt
otu_table(fun.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


#2018
fun.clean.filt <- fun.clean.ss2.18.f

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18 <- fun.clean.filt
otu_table(fun.clean.nolog.18) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18 <- fun.clean.filt
otu_table(fun.clean.log.18) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


### Normalization sorted by compartment
## Leaf
fun.clean.filt <- fun.clean.ss.17.leaf

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17.leaf <- fun.clean.filt
otu_table(fun.clean.nolog.17.leaf) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17.leaf <- fun.clean.filt
otu_table(fun.clean.log.17.leaf) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Stem
fun.clean.filt <- fun.clean.ss.17.stem

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17.stem <- fun.clean.filt
otu_table(fun.clean.nolog.17.stem) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17.stem <- fun.clean.filt
otu_table(fun.clean.log.17.stem) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Root
fun.clean.filt <- fun.clean.ss.17.root

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17.root <- fun.clean.filt
otu_table(fun.clean.nolog.17.root) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17.root <- fun.clean.filt
otu_table(fun.clean.log.17.root) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Seed
fun.clean.filt <- fun.clean.ss.17.seed

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17.seed <- fun.clean.filt
otu_table(fun.clean.nolog.17.seed) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17.seed <- fun.clean.filt
otu_table(fun.clean.log.17.seed) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Soil
fun.clean.filt <- fun.clean.ss.17.soil

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17.soil <- fun.clean.filt
otu_table(fun.clean.nolog.17.soil) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17.soil <- fun.clean.filt
otu_table(fun.clean.log.17.soil) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

### 2018 Samples

## Leaf
fun.clean.filt <- fun.clean.ss.18.f.leaf

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.leaf <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.leaf) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.leaf <- fun.clean.filt
otu_table(fun.clean.log.18.f.leaf) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

## Stem
fun.clean.filt <- fun.clean.ss.18.f.stem

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.stem <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.stem) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.stem <- fun.clean.filt
otu_table(fun.clean.log.18.f.stem) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Root
fun.clean.filt <- fun.clean.ss.18.f.root

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.root <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.root) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.root <- fun.clean.filt
otu_table(fun.clean.log.18.f.root) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Seed
fun.clean.filt <- fun.clean.ss.18.f.seed

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.seed <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.seed) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.seed <- fun.clean.filt
otu_table(fun.clean.log.18.f.seed) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


## Bulk soil
fun.clean.filt <- fun.clean.ss.18.f.BS

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.BS <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.BS) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.BS <- fun.clean.filt
otu_table(fun.clean.log.18.f.BS) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



## Rhizosphere
fun.clean.filt <- fun.clean.ss.18.f.RS

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18.f.RS <- fun.clean.filt
otu_table(fun.clean.nolog.18.f.RS) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18.f.RS <- fun.clean.filt
otu_table(fun.clean.log.18.f.RS) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



#### 
fun.clean.filt <- fun.clean.ss.17
fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.17 <- fun.clean.filt
otu_table(fun.clean.nolog.17) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.17 <- fun.clean.filt
otu_table(fun.clean.log.17) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


fun.clean.filt <- fun.clean.ss.18
(filt.sample <- sample_sums(fun.clean.filt) > 10)
sum(sample_sums(fun.clean.filt) <= 10)  ## 1 sample discarded
fun.clean.filt <- prune_samples(filt.sample, fun.clean.filt)
fun.clean.filt ## 979 samples <- 984 samples

fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.18 <- fun.clean.filt
otu_table(fun.clean.nolog.18) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.18 <- fun.clean.filt
otu_table(fun.clean.log.18) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


##### merged samples
fun.clean.filt <- fun.clean.ss.18.microhabitat
(filt.sample <- sample_sums(fun.clean.filt) > 10)
sum(sample_sums(fun.clean.filt) <= 10)  ## 1 sample discarded
fun.clean.filt <- prune_samples(filt.sample, fun.clean.filt)
fun.clean.filt ## 979 samples <- 984 samples

fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog.habitat <- fun.clean.filt
otu_table(fun.clean.nolog.habitat) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.habitat <- fun.clean.filt
otu_table(fun.clean.log.habitat) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)





bac.clean.filt <- bac.clean.ss.18.microhabitat
(filt.sample <- sample_sums(bac.clean.filt) > 10)
sum(sample_sums(bac.clean.filt) <= 10)  ## 1 sample discarded
bac.clean.filt <- prune_samples(filt.sample, bac.clean.filt)
bac.clean.filt ## 979 samples <- 984 samples

bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.habitat <- bac.clean.filt
otu_table(bac.clean.nolog.habitat) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.habitat <- bac.clean.filt
otu_table(bac.clean.log.habitat) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)
