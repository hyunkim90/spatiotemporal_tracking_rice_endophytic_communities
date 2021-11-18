
### Source tracking using FEAST
## Loading required packages
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
lapply(Packages, library, character.only = TRUE)

library(FEAST)

### Count table
bac.clean.ss.18 <- subset_samples(bac.clean.ss.f, Year == "year_2018")
bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)

fun.clean.ss.18 <- subset_samples(fun.clean.ss.f, Year == "year_2018")
fun.clean.ss.18 <- phyloseq::filter_taxa(fun.clean.ss.18, function(x) sum(x) != 0, TRUE)

otu.bac.18 <- otu_table(bac.clean.ss.18)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.ss.18)
otu.fun.18 <- data.frame(otu.fun.18)


###Designating source-sink
## grouped by id
bac.meta.tab  <- sample_data(bac.clean.ss.18)
bac.meta.tab  <- data.frame(bac.meta.tab)
# 
fun.meta.tab  <- sample_data(fun.clean.ss.18)
fun.meta.tab  <- data.frame(fun.meta.tab)
# 
#  
# write.csv(bac.meta.tab,"bac.meta.tab_for source tracking.csv")
# write.csv(fun.meta.tab,"fun.meta.tab_for source tracking.csv")

bac.meta.tab<-read.csv("bac.meta.tab_for source tracking.csv")
fun.meta.tab<-read.csv("fun.meta.tab_for source tracking.csv")

rownames(bac.meta.tab) <- bac.meta.tab$SampleID
rownames(fun.meta.tab) <- fun.meta.tab$SampleID

sample_data(bac.clean.ss.18) <- sample_data(bac.meta.tab)
sample_data(fun.clean.ss.18) <- sample_data(fun.meta.tab)

### Merge replicates
bac.clean.rep<-merge_samples(bac.clean.ss.18, "SampleID_2", fun = sum)
fun.clean.rep<-merge_samples(fun.clean.ss.18, "SampleID_2", fun = sum)

# bac.meta.tab  <- sample_data(bac.clean.rep)
# bac.meta.tab  <- data.frame(bac.meta.tab)
# # 
# fun.meta.tab  <- sample_data(fun.clean.rep)
# fun.meta.tab  <- data.frame(fun.meta.tab)
# # 
# #  
# write.csv(bac.meta.tab,"bac.meta.tab_for source tracking_2.csv")
# write.csv(fun.meta.tab,"fun.meta.tab_for source tracking_2.csv")
# 
### load edited meta file
bac.meta.tab<-read.csv("bac.meta.tab_for source tracking_2.csv")
fun.meta.tab<-read.csv("fun.meta.tab_for source tracking_2.csv")
rownames(bac.meta.tab) <- bac.meta.tab$SampleID
rownames(fun.meta.tab) <- fun.meta.tab$SampleID


sample_data(bac.clean.rep) <- sample_data(bac.meta.tab)
sample_data(fun.clean.rep) <- sample_data(fun.meta.tab)

### OTU table
otu.bac.18 <- otu_table(bac.clean.rep)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep)
otu.fun.18 <- data.frame(otu.fun.18)

source("./FEAST/FEAST_opt.R")
source("./FEAST/utils.R")
source("./FEAST/utils_plot.R")
source("./FEAST/Infer_LatentVariables.R")
source("./FEAST/RcppExports.R")

### Take all samples into account

FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = bac.meta.tab, different_sources_flag = 1, dir_path = "./",
                              outfile="Bac_test1")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = fun.meta.tab, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_test1")

### Source contribution without developing seeds
### OTU table
bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(Env%in% c('G_76','G_90','G_106','G_120')))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)

fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(Env%in% c('G_76','G_90','G_106','G_120')))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)

otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)

meta.bac.18 <- sample_data(bac.clean.rep.sub)
meta.bac.18 <- data.frame(meta.bac.18)

meta.fun.18 <- sample_data(fun.clean.rep.sub)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="bac_test3")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_test3")


### Source contribution without developing seeds and compartments in 141 days
### OTU table
bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(Env%in% c('G_0','G_76','G_90','G_106','G_120')))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)

fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(Env%in% c('G_0','G_76','G_90','G_106','G_120')))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)

remove.141<-sample_data(subset_samples(bac.clean.rep.sub, Days == '141' &  SourceSink == "Source"))$SampleID
bac.clean.rep.sub <- subset_samples(bac.clean.rep.sub, !(SampleID %in% remove.141))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)


fun.clean.rep.sub <- subset_samples(fun.clean.rep.sub, !(SampleID %in% remove.141))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)

otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)

meta.bac.18 <- sample_data(bac.clean.rep.sub)
meta.bac.18 <- data.frame(meta.bac.18)

meta.fun.18 <- sample_data(fun.clean.rep.sub)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="bac_wo141_woG_")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_wo141_woG_")


bac.result <- read.csv('bac_wo141_woG__source_contributions_matrix_for_plot.csv')
head(bac.result)
bac.result$Sink1 <- as.numeric(as.character(bac.result$Sink1))
bac.result$Sink2 <- as.numeric(as.character(bac.result$Sink2))
bac.result$Sink3 <- as.numeric(as.character(bac.result$Sink3))
bac.result.melt <- melt(bac.result, na.rm = T)
head(bac.result.melt)

bac.result.compart <- bac.result.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))

bac.result.mean.compart <- bac.result.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


fun.result <- read.csv('fun_wo141_woG__source_contributions_matrix_for_plot.csv')
head(fun.result)
fun.result$Sink1 <- as.numeric(as.character(fun.result$Sink1))
fun.result$Sink2 <- as.numeric(as.character(fun.result$Sink2))
fun.result$Sink3 <- as.numeric(as.character(fun.result$Sink3))
fun.result.melt <- melt(fun.result, na.rm = T)
head(fun.result.melt)

fun.result.compart <- fun.result.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))

fun.result.mean.compart <- fun.result.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))




##Compartment level
bac.result.mean.compart$Days <- factor(bac.result.mean.compart$Days, levels = c("48","62","76","90","106","120","Unknown"))
bac.result.mean.compart$Compartment <- factor(bac.result.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))

fun.result.mean.compart$Days <- factor(fun.result.mean.compart$Days, levels = c("48","62","76","90","106","120","Unknown"))
fun.result.mean.compart$Compartment <- factor(fun.result.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))


g<-ggplot(fun.result.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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






### Source contribution of rice compartments
bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(Env%in% c('G_0','G_76','G_90','G_106','G_120')))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)

fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(Env%in% c('G_0','G_76','G_90','G_106','G_120')))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)

otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)

meta.bac.18 <- sample_data(bac.clean.rep.sub)
meta.bac.18 <- data.frame(meta.bac.18)

meta.fun.18 <- sample_data(fun.clean.rep.sub)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="bac_test4")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_test4")



#### Source contribution to each developmental stage (from 76 days to 141 days)
### load metadata
bac.meta.tab<-read.csv("bac.meta.tab_for source tracking_3.csv")
fun.meta.tab<-read.csv("fun.meta.tab_for source tracking_3.csv")
rownames(bac.meta.tab) <- bac.meta.tab$SampleID
rownames(fun.meta.tab) <- fun.meta.tab$SampleID
### OTU table

bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(Days%in% c('0','48','62')))
bac.clean.rep.sub <- subset_samples(bac.clean.rep.sub, !(Plot == "UF1C3" & Days == "76"))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)


fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(Days%in% c('0','48','62')))
fun.clean.rep.sub <- subset_samples(fun.clean.rep.sub, !(Plot == "UF1C3" & Days == "76"))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)


otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = bac.meta.tab, different_sources_flag = 1, dir_path = "./",
                              outfile="Bac_test2")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = fun.meta.tab, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_test2")



## Plotting
## Import raw table
bac.result1 <- read.csv('Bac_test1_source_contributions_matrix_for_plot.csv')
head(bac.result1)
bac.result1$Sink_1 <- as.numeric(as.character(bac.result1$Sink_1))
bac.result1$Sink_2 <- as.numeric(as.character(bac.result1$Sink_2))
bac.result1$Sink_3 <- as.numeric(as.character(bac.result1$Sink_3))
bac.result1.melt <- melt(bac.result1, na.rm = T)
head(bac.result1.melt)

bac.result1.section <- bac.result1.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
bac.result1.compart <- bac.result1.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

bac.result1.section <- bac.result1.melt %>% group_by(variable, Section,Days) %>% summarise(SumCont = sum(value))
bac.result1.compart <- bac.result1.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))

bac.result1.mean.section <- bac.result1.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
bac.result1.mean.compart <- bac.result1.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))

bac.result1.mean.compart$Kingdom <- "Bacteria"


fun.result1 <- read.csv('Fun_test1_source_contributions_matrix_for_plot.csv')
head(fun.result1)
fun.result1$Sink_1 <- as.numeric(as.character(fun.result1$Sink_1))
fun.result1$Sink_2 <- as.numeric(as.character(fun.result1$Sink_2))
fun.result1$Sink_3 <- as.numeric(as.character(fun.result1$Sink_3))
fun.result1.melt <- melt(fun.result1, na.rm = T)
head(fun.result1.melt)

fun.result1.section <- fun.result1.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
fun.result1.compart <- fun.result1.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

fun.result1.section <- fun.result1.melt %>% group_by(variable, Section,Days) %>% summarise(SumCont = sum(value))
fun.result1.compart <- fun.result1.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))


fun.result1.mean.section <- fun.result1.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
fun.result1.mean.compart <- fun.result1.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))

fun.result1.mean.compart$Kingdom <- "Fungi"


bac.fun.result1.mean.compart <- rbind(bac.result1.mean.compart,fun.result1.mean.compart)

##Bar plot
bac.fun.result1.mean.compart$Compartment <- factor(bac.fun.result1.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
g<-ggplot(bac.fun.result1.mean.compart , aes(x= Kingdom,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


g<-ggplot(bac.result1.compart , aes(x= variable,fill = Compartment,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


g<-ggplot(fun.result1.section , aes(x= variable,fill = Section,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


g<-ggplot(fun.result1.compart , aes(x= variable,fill = Compartment,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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






#### Each time
bac.result2 <- read.csv("Bac_test2_source_contributions_matrix_for_plot.csv")
bac.result2$Days <- as.factor(as.character(bac.result2$Days))
bac.result2.melt <- melt(bac.result2, na.rm = T)
head(bac.result2.melt)

bac.result2.section <- bac.result2.melt %>% group_by(variable, Microhabitat, Days) %>% summarise(SumCont = sum(value))
bac.result2.compart <- bac.result2.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))


g<-ggplot(bac.result2.section , aes(x= variable,fill = Microhabitat,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


g<-ggplot(bac.result2.compart , aes(x= variable,fill = Compartment,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


fun.result2 <- read.csv("fun_test2_source_contributions_matrix_for_plot.csv")
fun.result2$Days <- as.factor(as.character(fun.result2$Days))
fun.result2.melt <- melt(fun.result2, na.rm = T)
head(fun.result2.melt)

fun.result2.section <- fun.result2.melt %>% group_by(variable, Microhabitat, Days) %>% summarise(SumCont = sum(value))
fun.result2.compart <- fun.result2.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))


g<-ggplot(fun.result2.section , aes(x= variable,fill = Microhabitat,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


g<-ggplot(fun.result2.compart , aes(x= variable,fill = Compartment,y= SumCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


##Without G
bac.result3 <- read.csv('bac_test3_source_contributions_matrix_for_plot.csv')
head(bac.result3)
bac.result3$Sink1 <- as.numeric(as.character(bac.result3$Sink1))
bac.result3$Sink2 <- as.numeric(as.character(bac.result3$Sink2))
bac.result3$Sink3 <- as.numeric(as.character(bac.result3$Sink3))
bac.result3.melt <- melt(bac.result3, na.rm = T)
head(bac.result3.melt)

bac.result3.section <- bac.result3.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
bac.result3.compart <- bac.result3.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))



fun.result3 <- read.csv('fun_test3_source_contributions_matrix_for_plot.csv')
head(fun.result3)
fun.result3$Sink1 <- as.numeric(as.character(fun.result3$Sink1))
fun.result3$Sink2 <- as.numeric(as.character(fun.result3$Sink2))
fun.result3$Sink3 <- as.numeric(as.character(fun.result3$Sink3))
fun.result3.melt <- melt(fun.result3, na.rm = T)
head(fun.result3.melt)

fun.result3.section <- fun.result3.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
fun.result3.compart <- fun.result3.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))


##Without all G
bac.result4 <- read.csv('bac_test4_source_contributions_matrix_for_plot.csv')
head(bac.result4)
bac.result4$Sink1 <- as.numeric(as.character(bac.result4$Sink1))
bac.result4$Sink2 <- as.numeric(as.character(bac.result4$Sink2))
bac.result4$Sink3 <- as.numeric(as.character(bac.result4$Sink3))
bac.result4.melt <- melt(bac.result4, na.rm = T)
head(bac.result4.melt)

bac.result4.section <- bac.result4.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
bac.result4.compart <- bac.result4.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))



fun.result4 <- read.csv('fun_test4_source_contributions_matrix_for_plot.csv')
head(fun.result4)
fun.result4$Sink1 <- as.numeric(as.character(fun.result4$Sink1))
fun.result4$Sink2 <- as.numeric(as.character(fun.result4$Sink2))
fun.result4$Sink3 <- as.numeric(as.character(fun.result4$Sink3))
fun.result4.melt <- melt(fun.result4, na.rm = T)
head(fun.result4.melt)

fun.result4.section <- fun.result4.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
fun.result4.compart <- fun.result4.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))





### Mean contribution
bac.result1.melt
head(bac.result1.melt)

bac.result1.mean.section <- bac.result1.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result1.mean.compart <- bac.result1.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result1.melt
head(fun.result1.melt)

fun.result1.mean.section <- fun.result1.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result1.mean.compart <- fun.result1.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.result2.melt
head(bac.result2.melt)

bac.result2.mean.section <- bac.result2.section %>% group_by(Microhabitat,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result2.mean.compart <- bac.result2.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result2.melt
head(fun.result2.melt)

fun.result2.mean.section <- fun.result2.section %>% group_by(Microhabitat,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result2.mean.compart <- fun.result2.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.result3.mean.section <- bac.result3.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result3.mean.compart <- bac.result3.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result3.mean.section <- fun.result3.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result3.mean.compart <- fun.result3.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.result4.mean.section <- bac.result4.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result4.mean.compart <- bac.result4.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result4.mean.section <- fun.result4.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result4.mean.compart <- fun.result4.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

### Include days
bac.result4.section <- bac.result4.melt %>% group_by(variable, Section, Days) %>% summarise(SumCont = sum(value))
bac.result4.compart <- bac.result4.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))

fun.result4.section <- fun.result4.melt %>% group_by(variable, Section, Days) %>% summarise(SumCont = sum(value))
fun.result4.compart <- fun.result4.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))


bac.result4.mean.section <- bac.result4.section %>% group_by(Section, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result4.mean.compart <- bac.result4.compart %>% group_by(Compartment, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result4.mean.section <- fun.result4.section %>% group_by(Section, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result4.mean.compart <- fun.result4.compart %>% group_by(Compartment, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.result1.section <- bac.result1.melt %>% group_by(variable, Section, Days) %>% summarise(SumCont = sum(value))
bac.result1.compart <- bac.result1.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))

fun.result1.section <- fun.result1.melt %>% group_by(variable, Section, Days) %>% summarise(SumCont = sum(value))
fun.result1.compart <- fun.result1.melt %>% group_by(variable, Compartment, Days) %>% summarise(SumCont = sum(value))


bac.result1.mean.section <- bac.result1.section %>% group_by(Section, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
bac.result1.mean.compart <- bac.result1.compart %>% group_by(Compartment, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

fun.result1.mean.section <- fun.result1.section %>% group_by(Section, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))
fun.result1.mean.compart <- fun.result1.compart %>% group_by(Compartment, Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))




##Plotting
##Compartment level
bac.result2.mean.compart$Days <- factor(bac.result2.mean.compart$Days, levels = c("76","90","106","120","141"))

bac.result2.mean.compart$Compartment <- factor(bac.result2.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))

g<-ggplot(bac.result2.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



fun.result2.mean.compart$Days <- factor(fun.result2.mean.compart$Days, levels = c("76","90","106","120","141"))

fun.result2.mean.compart$Compartment <- factor(fun.result2.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))

g<-ggplot(fun.result2.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



bac.result3.mean.compart$Compartment <- factor(bac.result3.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))

g<-ggplot(bac.result3.mean.compart, aes(x= Compartment,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



fun.result3.mean.compart$Compartment <- factor(fun.result3.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))

g<-ggplot(fun.result3.mean.compart, aes(x= Compartment,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


bac.result1.mean.compart$Compartment <- factor(bac.result1.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
bac.result1.mean.compart$Days <- factor(bac.result1.mean.compart$Days, levels = c("0","48","62","76","90","106","120","141", "Unknown"))

g<-ggplot(bac.result1.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



fun.result1.mean.compart$Compartment <- factor(fun.result1.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
fun.result1.mean.compart$Days <- factor(fun.result1.mean.compart$Days, levels = c("0","48","62","76","90","106","120","141", "Unknown"))

g<-ggplot(fun.result1.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


bac.result4.mean.compart$Compartment <- factor(bac.result4.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))
bac.result4.mean.compart$Days <- factor(bac.result4.mean.compart$Days, levels = c("48","62","76","90","106","120","141", "Unknown"))


g<-ggplot(bac.result4.mean.compart, aes(x= Compartment,fill = Days,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



fun.result4.mean.compart$Compartment <- factor(fun.result4.mean.compart$Compartment, levels = c("L","S","R","RS","BS","Unknown"))
fun.result4.mean.compart$Days <- factor(fun.result4.mean.compart$Days, levels = c("48","62","76","90","106","120","141", "Unknown"))

g<-ggplot(fun.result4.mean.compart, aes(x= Compartment,fill = Days,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


##Section level
bac.result2.mean.section$Days <- factor(bac.result2.mean.section$Days, levels = c("76","90","106","120","141"))

bac.result2.mean.section$Microhabitat <- factor(bac.result2.mean.section$Microhabitat, levels = c("FL","L3","L2","L1","S9","S8","S7","S6","S5","S4","S3","S2","S1","R","RS","BS","Unknown"))

g<-ggplot(bac.result2.mean.section, aes(x= Days,fill = Microhabitat,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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



fun.result2.mean.section$Days <- factor(fun.result2.mean.section$Days, levels = c("76","90","106","120","141"))

fun.result2.mean.section$Microhabitat <- factor(fun.result2.mean.section$Microhabitat, levels = c("FL","L3","L2","L1","S9","S8","S7","S6","S5","S4","S3","S2","S1","R","RS","BS","Unknown"))

g<-ggplot(fun.result2.mean.section, aes(x= Days,fill = Microhabitat,y= MeanCont)) + 
  geom_bar(width = 0.5,position = position_dodge(preserve = 'single'), stat = "identity")+
  theme(aspect.ratio =0.3)+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                         position=position_dodge(0.5), width = 0.2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


### Without 141 days samples
remove.141<-sample_data(subset_samples(bac.clean.rep, Days == '141' &  SourceSink == "Source"))$SampleID
bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(SampleID %in% remove.141))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)


fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(SampleID %in% remove.141))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)


otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)


meta.bac.18 <- sample_data(bac.clean.rep.sub)
meta.bac.18 <- data.frame(meta.bac.18)

meta.fun.18 <- sample_data(fun.clean.rep.sub)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="Bac_wo141")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_wo141")


## plot
bac.wo141 <- read.csv('Bac_wo141_source_contributions_matrix_for_plot.csv')
head(bac.wo141)
bac.wo141$Sink1 <- as.numeric(as.character(bac.wo141$Sink1))
bac.wo141$Sink2 <- as.numeric(as.character(bac.wo141$Sink2))
bac.wo141$Sink3 <- as.numeric(as.character(bac.wo141$Sink3))
bac.wo141.melt <- melt(bac.wo141, na.rm = T)
head(bac.wo141.melt)

bac.wo141.section <- bac.wo141.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
bac.wo141.compart <- bac.wo141.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

bac.wo141.mean.section <- bac.wo141.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
bac.wo141.mean.compart <- bac.wo141.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))

bac.wo141.mean.compart$Kingdom <- "Bacteria"


fun.wo141 <- read.csv('Fun_wo141_source_contributions_matrix_for_plot.csv')
head(fun.wo141)
fun.wo141$Sink_1 <- as.numeric(as.character(fun.wo141$Sink_1))
fun.wo141$Sink_2 <- as.numeric(as.character(fun.wo141$Sink_2))
fun.wo141$Sink_3 <- as.numeric(as.character(fun.wo141$Sink_3))
fun.wo141.melt <- melt(fun.wo141, na.rm = T)
head(fun.wo141.melt)

fun.wo141.section <- fun.wo141.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
fun.wo141.compart <- fun.wo141.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

fun.wo141.mean.section <- fun.wo141.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
fun.wo141.mean.compart <- fun.wo141.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))
fun.wo141.mean.compart$Kingdom <- "Fungi"

bac.fun.wo141.mean.compart <- rbind(bac.wo141.mean.compart,fun.wo141.mean.compart)
bac.fun.wo141.mean.compart$Compartment <- factor(bac.fun.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
g<-ggplot(bac.fun.result1.mean.compart , aes(x= Kingdom,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


### Contribution of each time
##Sorted by Days
bac.wo141.section <- bac.wo141.melt %>% group_by(variable, Section,Days) %>% summarise(SumCont = sum(value))
bac.wo141.compart <- bac.wo141.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))

bac.wo141.mean.compart <- bac.wo141.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))

write.csv(bac.wo141.mean.compart)

fun.wo141.section <- fun.wo141.melt %>% group_by(variable, Section,Days) %>% summarise(SumCont = sum(value))
fun.wo141.compart <- fun.wo141.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))

fun.wo141.mean.compart <- fun.wo141.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.wo141.mean.compart$Compartment <- factor(bac.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
bac.wo141.mean.compart$Days <- factor(bac.wo141.mean.compart$Days, levels = c("0","48","62","76","90","106","120","Unknown"))


fun.wo141.mean.compart$Compartment <- factor(fun.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
fun.wo141.mean.compart$Days <- factor(fun.wo141.mean.compart$Days, levels = c("0","48","62","76","90","106","120","Unknown"))

g<-ggplot(fun.wo141.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity", position = position_dodge())+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                                         position=position_dodge(0.5), width = 0.2)+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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

#### Considering only persistent (inherited) OTUs
persistent.bac
persistent.fun

## all samples
otu.bac.18 <- otu_table(bac.clean.rep)
otu.bac.18 <- data.frame(otu.bac.18)

otu.fun.18 <- otu_table(fun.clean.rep)
otu.fun.18 <- data.frame(otu.fun.18)

otu.bac.18.persist <- subset(t(otu.bac.18), colnames(otu.bac.18) %in% persistent.bac)
otu.bac.18.persist<-t(otu.bac.18.persist)

otu.bac.18.persist<-subset(otu.bac.18.persist, rowSums(otu.bac.18.persist) > 0)

rowSums(otu.bac.18.persist)

meta.bac.18 <- sample_data(bac.clean.rep)
meta.bac.18 <- data.frame(meta.bac.18)


otu.fun.18.persist <- subset(t(otu.fun.18), colnames(otu.fun.18) %in% persistent.fun)
otu.fun.18.persist<-t(otu.fun.18.persist)

otu.fun.18.persist<-subset(otu.fun.18.persist, rowSums(otu.fun.18.persist) > 0)

rowSums(otu.fun.18.persist)

meta.fun.18 <- sample_data(fun.clean.rep)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18.persist, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="Bac_persist1")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18.persist, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="Fun_persist1")


### wo 141 samples
remove.141<-sample_data(subset_samples(bac.clean.rep, Days == '141' &  SourceSink == "Source"))$SampleID
bac.clean.rep.sub <- subset_samples(bac.clean.rep, !(SampleID %in% remove.141))
bac.clean.rep.sub  <- phyloseq::filter_taxa(bac.clean.rep.sub , function(x) sum(x) != 0, TRUE)


fun.clean.rep.sub <- subset_samples(fun.clean.rep, !(SampleID %in% remove.141))
fun.clean.rep.sub  <- phyloseq::filter_taxa(fun.clean.rep.sub , function(x) sum(x) != 0, TRUE)


otu.bac.18 <- otu_table(bac.clean.rep.sub)
otu.bac.18 <- data.frame(otu.bac.18)

otu.bac.18.persist <- subset(t(otu.bac.18), colnames(otu.bac.18) %in% persistent.bac)
otu.bac.18.persist<-t(otu.bac.18.persist)

otu.bac.18.persist<-subset(otu.bac.18.persist, rowSums(otu.bac.18.persist) > 0)



otu.fun.18 <- otu_table(fun.clean.rep.sub)
otu.fun.18 <- data.frame(otu.fun.18)

otu.fun.18.persist <- subset(t(otu.fun.18), colnames(otu.fun.18) %in% persistent.fun)
otu.fun.18.persist<-t(otu.fun.18.persist)

otu.fun.18.persist<-subset(otu.fun.18.persist, rowSums(otu.fun.18.persist) > 0)



meta.bac.18 <- sample_data(bac.clean.rep.sub)
meta.bac.18 <- data.frame(meta.bac.18)

meta.fun.18 <- sample_data(fun.clean.rep.sub)
meta.fun.18 <- data.frame(meta.fun.18)


FEAST_output.bac <- FEAST_opt(C = otu.bac.18.persist, metadata = meta.bac.18, different_sources_flag = 1, dir_path = "./",
                              outfile="Bac_wo141_persist")
FEAST_output.fun <- FEAST_opt(C = otu.fun.18.persist, metadata = meta.fun.18, different_sources_flag = 1, dir_path = "./",
                              outfile="fun_wo141_persist")

##plot
bac.wo141 <- read.csv('Bac_wo141_persist_source_contributions_matrix_for_plot.csv')
head(bac.wo141)
bac.wo141$Sink1 <- as.numeric(as.character(bac.wo141$Sink1))
bac.wo141$Sink2 <- as.numeric(as.character(bac.wo141$Sink2))
bac.wo141$Sink3 <- as.numeric(as.character(bac.wo141$Sink3))
bac.wo141.melt <- melt(bac.wo141, na.rm = T)
head(bac.wo141.melt)

bac.wo141.section <- bac.wo141.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
bac.wo141.compart <- bac.wo141.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

bac.wo141.mean.section <- bac.wo141.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
bac.wo141.mean.compart <- bac.wo141.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))

bac.wo141.mean.compart$Kingdom <- "Bacteria"

fun.wo141 <- read.csv('Fun_wo141_persist_source_contributions_matrix_for_plot.csv')
head(fun.wo141)
fun.wo141$Sink1 <- as.numeric(as.character(fun.wo141$Sink1))
fun.wo141$Sink2 <- as.numeric(as.character(fun.wo141$Sink2))
fun.wo141$Sink3 <- as.numeric(as.character(fun.wo141$Sink3))
fun.wo141.melt <- melt(fun.wo141, na.rm = T)
head(fun.wo141.melt)

fun.wo141.section <- fun.wo141.melt %>% group_by(variable, Section) %>% summarise(SumCont = sum(value))
fun.wo141.compart <- fun.wo141.melt %>% group_by(variable, Compartment) %>% summarise(SumCont = sum(value))

fun.wo141.mean.section <- fun.wo141.section %>% group_by(Section) %>% summarise(MeanCont = mean(SumCont))
fun.wo141.mean.compart <- fun.wo141.compart %>% group_by(Compartment) %>% summarise(MeanCont = mean(SumCont))
fun.wo141.mean.compart$Kingdom <- "Fungi"

bac.fun.wo141.mean.compart <- rbind(bac.wo141.mean.compart,fun.wo141.mean.compart)
bac.fun.wo141.mean.compart$Compartment <- factor(bac.fun.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
g<-ggplot(bac.fun.wo141.mean.compart , aes(x= Kingdom,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity")+
  theme(aspect.ratio =2)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


### Sorted by days
bac.wo141.compart <- bac.wo141.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))
bac.wo141.mean.compart <- bac.wo141.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


bac.wo141.mean.compart$Compartment <- factor(bac.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
bac.wo141.mean.compart$Days <- factor(bac.wo141.mean.compart$Days, levels = c("0","48","62","76","90","106","120"))


g<-ggplot(bac.wo141.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity", position = position_dodge())+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                                        position=position_dodge(0.5), width = 0.2)+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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


fun.wo141.compart <- fun.wo141.melt %>% group_by(variable, Compartment,Days) %>% summarise(SumCont = sum(value))
fun.wo141.mean.compart <- fun.wo141.compart %>% group_by(Compartment,Days) %>% summarise(MeanCont = mean(SumCont), SDCont = sd(SumCont))


fun.wo141.mean.compart$Compartment <- factor(fun.wo141.mean.compart$Compartment, levels = c("G","L","S","R","RS","BS","Unknown"))
fun.wo141.mean.compart$Days <- factor(fun.wo141.mean.compart$Days, levels = c("0","48","62","76","90","106","120"))


g<-ggplot(fun.wo141.mean.compart, aes(x= Days,fill = Compartment,y= MeanCont)) + 
  geom_bar(width = 0.5,stat = "identity", position = position_dodge())+geom_errorbar(aes(ymin=MeanCont-SDCont,ymax=MeanCont+SDCont),
                                                                                     position=position_dodge(0.5), width = 0.2)+
  theme(aspect.ratio =0.3)+
  #scale_linetype_manual(values = c("BS"="solid","RS"="solid","R"="solid","S1"="twodash","S2"="twodash","S3"="twodash","S4"="longdash","S5"="longdash","S6"="longdash","S7"="longdash","S8"="longdash","S9"="longdash","L1"="dashed","L2"="dashed","L3"="dashed","FL"="dashed","G"="dotted"))+
  ylab("Source contribution\n") + theme(plot.title = element_text(size = 5,hjust = 0.5, face='bold')) + 
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

