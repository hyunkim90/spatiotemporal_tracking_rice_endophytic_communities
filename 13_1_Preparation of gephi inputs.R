
### sparcc let's go!!
### making gephi input!!
### and also making edge and node list
library(reshape2)
library(tidyr)
## import edge list
## 

### Gephi input
get_edge_node_withp <- function(cor_mat,p_mat,edgename,nodename,threshold){
  head(cor_mat)
  dim(cor_mat)
  
  any(is.na(cor_mat))
  cor_mat[upper.tri(cor_mat)] <- NA
  df_cor <- reshape2::melt(cor_mat, varnames = c('Source', 'Target'), na.rm = TRUE)
  
  dim(df_cor) ## 258840
  head(df_cor)
  unique(df_cor$variable)
  tail(df_cor)
  any(is.na(df_cor))
  
  # threshold <-  0.3
  df_cor_sig <- df_cor %>% filter(abs(value)>=threshold)  ## 0.3 is 22 edges   ## 0.2 is 124 edges  ## 0.25 is 48
  print(dim(df_cor_sig))  ## total 94 edges
  # colnames(df_cor_sig) <- c('Source','Target','Cor')
  
  head(p_mat)
  dim(p_mat)
  
  any(is.na(p_mat))
  p_mat[upper.tri(p_mat)] <- NA
  df_p <- reshape2::melt(p_mat, varnames = c('Source', 'Target'), na.rm = TRUE)
  
  dim(df_p) ## 4186
  head(df_p)
  tail(df_p)
  any(is.na(df_p))
  
  df_cor_p <- left_join(df_cor, df_p, by=c('OTU_id'='OTU_id','variable'='variable'))
  head(df_cor_p)
  
  colnames(df_cor_p) <- c('Source','Target','Cor','pseudo_p')
  df_cor_sig <- df_cor_p %>% filter(pseudo_p < 0.05) %>% filter(abs(Cor)>threshold)  ## 15
  dim(df_cor_sig)   ## 94
  head(df_cor_sig)
  df_cor_sig$positive <- df_cor_sig$Cor > 0
  
  write.table(df_cor_sig, file=edgename, quote=FALSE, sep='\t', row.names = F)
  v.source <- df_cor_sig$Source
  v.target <- df_cor_sig$Target
  v.id <- union(v.source,v.target)
  length(v.id)
  
  df.node <- data.frame(id=v.id)
  df.node$label <- df.node$id
  #df.node$sep <- df.node$id
  
  df.node$Bacteria <- ifelse(grepl("^B",df.node$id),"Bacteria","Fungi")
 
  write.table(df.node, file=nodename, quote=FALSE, sep='\t', row.names = F)
}


get_gephi_input_withp <- function(keyword,front_number, threshold){
  # keyword <- 'mer.rel.ss.5'
  coreword <- paste0(keyword,'_median_correlation','.csv')
  coreword
  cor_mat <- read.csv(file = coreword)
  keyword
  pword <- paste0(keyword,'_pvalues','.csv')
  p_mat <- read.csv(file = pword)
  edgeword <- paste0(as.character(front_number),'_',as.character(threshold),'_','edge_',keyword,'.tsv')
  # threshold <- 0.2
  # front_number <- 1
  edgeword
  nodeword <- paste0(as.character(front_number),'_',as.character(threshold),'_','node_',keyword,'.tsv')
  nodeword
  
  get_edge_node_withp(cor_mat,p_mat,edgeword,nodeword,threshold)
}

###Thresholds were decided via Random matrix theory
##G0 0.11
##G76 0.15
##G90 0.1
##G106 0.116
##G120 0.1
##G141 0.1
#G0
get_gephi_input_withp('G0',0,0.11)
get_gephi_input_withp('G0',0,0.2)
get_gephi_input_withp('G0',1,0.3)
get_gephi_input_withp('G0',2,0.4)
get_gephi_input_withp('G0',3,0.5)

#G76
get_gephi_input_withp('G76',0,0.15)
get_gephi_input_withp('G76',0,0.2)
get_gephi_input_withp('G76',1,0.3)
get_gephi_input_withp('G76',2,0.4)
get_gephi_input_withp('G76',3,0.5)


#G90
get_gephi_input_withp('G90',0,0.1)
get_gephi_input_withp('G90',0,0.2)
get_gephi_input_withp('G90',1,0.3)
get_gephi_input_withp('G90',2,0.4)
get_gephi_input_withp('G90',3,0.5)

#G106
get_gephi_input_withp('G106',0,0.116)
get_gephi_input_withp('G106',0,0.2)
get_gephi_input_withp('G106',1,0.3)
get_gephi_input_withp('G106',2,0.4)
get_gephi_input_withp('G106',3,0.5)

#G120
get_gephi_input_withp('G120',0,0.1)
get_gephi_input_withp('G120',0,0.2)
get_gephi_input_withp('G120',1,0.3)
get_gephi_input_withp('G120',2,0.4)
get_gephi_input_withp('G120',3,0.5)

#G141
get_gephi_input_withp('G141',0,0.1)
get_gephi_input_withp('G141',0,0.2)
get_gephi_input_withp('G141',1,0.3)
get_gephi_input_withp('G141',2,0.4)
get_gephi_input_withp('G141',3,0.5)



#L1_48
get_gephi_input_withp('L1_48',1,0.5)

#L1_62
get_gephi_input_withp('L1_62',1,0.5)

#L1_76
get_gephi_input_withp('L1_76',1,0.5)

#L1_90
get_gephi_input_withp('L1_90',1,0.5)

#L1_106
get_gephi_input_withp('L1_106',1,0.5)

#L1_120
get_gephi_input_withp('L1_120',1,0.5)

#L1_141
get_gephi_input_withp('L1_141',1,0.5)




#L2_48
get_gephi_input_withp('L2_48',1,0.5)

#L2_62
get_gephi_input_withp('L2_62',1,0.5)

#L2_76
get_gephi_input_withp('L2_76',1,0.5)

#L2_90
get_gephi_input_withp('L2_90',1,0.5)

#L2_106
get_gephi_input_withp('L2_106',1,0.5)

#L2_120
get_gephi_input_withp('L2_120',1,0.5)

#L2_141
get_gephi_input_withp('L2_141',1,0.5)



#L3_62
get_gephi_input_withp('L3_62',1,0.5)

#L3_76
get_gephi_input_withp('L3_76',1,0.5)

#L3_90
get_gephi_input_withp('L3_90',1,0.5)

#L3_106
get_gephi_input_withp('L3_106',1,0.5)

#L3_120
get_gephi_input_withp('L3_120',1,0.5)

#L3_141
get_gephi_input_withp('L3_141',1,0.5)



#FL_76
get_gephi_input_withp('FL_76',1,0.5)

#FL_90
get_gephi_input_withp('FL_90',1,0.3)

#FL_106
get_gephi_input_withp('FL_106',1,0.5)

#FL_120
get_gephi_input_withp('FL_120',1,0.5)

#FL_141
get_gephi_input_withp('FL_141',1,0.5)






#S1_48
get_gephi_input_withp('S1_48',1,0.5)

#S1_62
get_gephi_input_withp('S1_62',1,0.5)

#S1_76
get_gephi_input_withp('S1_76',1,0.5)

#S1_90
get_gephi_input_withp('S1_90',1,0.5)

#S1_106
get_gephi_input_withp('S1_106',1,0.5)

#S1_120
get_gephi_input_withp('S1_120',1,0.5)

#S1_141
get_gephi_input_withp('S1_141',1,0.5)



#S2_48
get_gephi_input_withp('S2_48',1,0.5)

#S2_62
get_gephi_input_withp('S2_62',1,0.5)

#S2_76
get_gephi_input_withp('S2_76',1,0.5)

#S2_90
get_gephi_input_withp('S2_90',1,0.5)

#S2_106
get_gephi_input_withp('S2_106',1,0.5)

#S2_120
get_gephi_input_withp('S2_120',1,0.5)


#S2_141
get_gephi_input_withp('S2_141',1,0.5)


#S3_48
get_gephi_input_withp('S3_48',1,0.5)

#S3_62
get_gephi_input_withp('S3_62',1,0.5)

#S3_76
get_gephi_input_withp('S3_76',1,0.5)

#S3_90
get_gephi_input_withp('S3_90',1,0.5)

#S3_106
get_gephi_input_withp('S3_106',1,0.5)

#S3_120
get_gephi_input_withp('S3_120',1,0.5)

#S3_141
get_gephi_input_withp('S3_141',1,0.5)


#S4_62
get_gephi_input_withp('S4_62',1,0.5)

#S4_76
get_gephi_input_withp('S4_76',1,0.5)

#S4_90
get_gephi_input_withp('S4_90',1,0.5)

#S4_106
get_gephi_input_withp('S4_106',1,0.5)

#S4_120
get_gephi_input_withp('S4_120',1,0.5)

#S4_141
get_gephi_input_withp('S4_141',1,0.5)



#S5_76
get_gephi_input_withp('S5_76',1,0.5)

#S5_90
get_gephi_input_withp('S5_90',1,0.5)

#S5_106
get_gephi_input_withp('S5_106',1,0.5)

#S5_120
get_gephi_input_withp('S5_120',1,0.5)

#S5_141
get_gephi_input_withp('S5_141',1,0.5)



#S6_76
get_gephi_input_withp('S6_76',1,0.5)

#S6_90
get_gephi_input_withp('S6_90',1,0.5)

#S6_106
get_gephi_input_withp('S6_106',1,0.5)

#S6_120
get_gephi_input_withp('S6_120',1,0.5)

#S6_141
get_gephi_input_withp('S6_141',1,0.5)



#S7_90
get_gephi_input_withp('S7_90',1,0.5)

#S7_106
get_gephi_input_withp('S7_106',1,0.5)

#S7_120
get_gephi_input_withp('S7_120',1,0.5)

#S7_141
get_gephi_input_withp('S7_141',1,0.5)


#S8_90
get_gephi_input_withp('S8_90',1,0.5)

#S8_106
get_gephi_input_withp('S8_106',1,0.5)

#S8_120
get_gephi_input_withp('S8_120',1,0.5)

#S8_141
get_gephi_input_withp('S8_141',1,0.5)



#S9_90
get_gephi_input_withp('S9_90',1,0.5)

#S9_106
get_gephi_input_withp('S9_106',1,0.5)

#S9_120
get_gephi_input_withp('S9_120',1,0.5)

#S9_141
get_gephi_input_withp('S9_141',1,0.5)



#R_48
get_gephi_input_withp('R_48',1,0.5)

#R_62
get_gephi_input_withp('R_62',1,0.5)

#R_76
get_gephi_input_withp('R_76',1,0.5)

#R_90
get_gephi_input_withp('R_90',1,0.5)

#R_106
get_gephi_input_withp('R_106',1,0.5)

#R_120
get_gephi_input_withp('R_120',1,0.5)

#R_141
get_gephi_input_withp('R_141',1,0.5)



#RS_48
get_gephi_input_withp('RS_48',1,0.5)

#RS_62
get_gephi_input_withp('RS_62',1,0.5)

#RS_76
get_gephi_input_withp('RS_76',1,0.5)

#RS_90
get_gephi_input_withp('RS_90',1,0.5)

#RS_106
get_gephi_input_withp('RS_106',1,0.5)

#RS_120
get_gephi_input_withp('RS_120',1,0.5)

#RS_141
get_gephi_input_withp('RS_141',1,0.5)


#BS_48
get_gephi_input_withp('BS_48',1,0.5)

#BS_62
get_gephi_input_withp('BS_62',1,0.5)

#BS_76
get_gephi_input_withp('BS_76',1,0.5)

#BS_90
get_gephi_input_withp('BS_90',1,0.5)

#BS_106
get_gephi_input_withp('BS_106',1,0.5)

#BS_120
get_gephi_input_withp('BS_120',1,0.5)

#BS_141
get_gephi_input_withp('BS_141',1,0.5)