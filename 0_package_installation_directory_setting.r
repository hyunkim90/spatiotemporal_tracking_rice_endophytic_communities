## import libraries : phyloseq and microbiome
library(dplyr)
library(forcats) 
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(scales)
library(grid)
library(reshape2)
library(seqtime)
library(agricolae)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
library(multifunc)
library(FSA)
library(rcompanion)
library(seqinr)

# Set plotting theme
theme_set(theme_bw())

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("OTU_table_final_bac.biom")

### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
map <- read.table(file = 'Metadata_bac.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID
rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("rooted_tree_bac.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)
phy   ## 17133 otus

## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy  ## 17133 OTUs



## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("otu_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = 'sample_metadata.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("tree.nwk") #rooted tree
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 6159 otus

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun
summarize_phyloseq(fun)
sort(colSums(otu_table(fun)))
