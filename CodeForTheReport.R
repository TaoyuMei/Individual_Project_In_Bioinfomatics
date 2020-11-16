##### formal code that may be submitted as supplementary material for the final report ########
{
library(GEOmetadb)
library(tidyverse)
library(stringr)
library(crossmeta)
library(hash)
library(clusterProfiler)
library(pathview)
library(limma)
library(pheatmap)
library(msigdbr)
library(exploBATCH)
library(sva)
setwd("D:/datasets and methods for all projects/Microarray_for_mRNA")
}


##### mapping DEGs from previously published meta-analysis to the KEGG mitophagy pathway

MapMitoMeta <- function(region){  # better use EntrezID from the original supplementary table
  GeneSymbol_FoldChange <- read.table(paste0("./MetaAnal/", region, ".rnk"),
                                      row.names = 1)  # the original table only include those with p-adj < 0.05
  GeneSymbol_FoldChange.M <- as.matrix(GeneSymbol_FoldChange)

  pv.out.mitophagy <- pathview(gene.data = GeneSymbol_FoldChange.M[, 1],
                               gene.idtype = "SYMBOL", 
                               kegg.dir = "./pathview",  # provide the xml file by myself in this dir
                               pathway.id = "04137", species = "hsa",
                               same.layer = F, kegg.native = T, 
                               out.suffix = paste0(region, ".mitophagy.meta"))
}

MapMitoMeta("Frontal Lobe")
MapMitoMeta("Cerebellum")
MapMitoMeta("Parietal Lobe")
MapMitoMeta("Temporal Lobe")


MapMitoAlzData <- function(region){
  tbltmp <- read_tsv(paste0("./AlzData/", region, "_all_adj.txt"))
  tbltmp <- dplyr::filter(tbltmp, `P-adj` < 0.05)  # remove those with P-adj ≥ 0.05
  GeneSymbol_FoldChange <- tbltmp$log2FC
  names(GeneSymbol_FoldChange) <- tbltmp$GENE
  
  pv.out.mitophagy <- pathview(gene.data = GeneSymbol_FoldChange,
                               gene.idtype = "SYMBOL", 
                               kegg.dir = "./pathview",  # provide the xml file by myself in this dir
                               pathway.id = "04137", species = "hsa",
                               same.layer = F, kegg.native = T, 
                               out.suffix = paste0(region, ".mitophagy.AlzData"))
}

MapMitoAlzData("EC")
MapMitoAlzData("FC")
MapMitoAlzData("HP")
MapMitoAlzData("TC")


##### GSEA via clusterProfiler and MSigDb

# ALL.msigdb <- read.gmt("MSigDb/msigdb.v7.0.symbols.gmt")
ALL.msigdb <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::select(gs_name, entrez_gene)

region <- "Frontal Lobe"

# create a named vector of log2FC, with gene symbols as names
geneList <- read.table(paste0("./MetaAnal/", region, "_entrez.rnk"), row.names = 1)  
# the original table only include genes with p-adj < 0.05
geneList <- as.matrix(geneList)[, 1]
geneList <- sort(geneList, decreasing = TRUE)
geneList <- 2 ^ geneList  # reverse log2 transform, useless


# create a vector of gene symbols
gene <- names(geneList)

# GSEA and over-representation analysis (default p-adj < 0.05)
msigdb <- GSEA(geneList, TERM2GENE = ALL.msigdb)

msigdb.down <- GSEA(geneList[geneList < 0], TERM2GENE = ALL.msigdb)
msigdb.up <- GSEA(geneList[geneList > 0], TERM2GENE = ALL.msigdb)



##### choosing, collecting and re-analysing GEO data


# read in all GSE numbers that have been manually collected from GEO and PubMed
gse_numbers <- as.character(read.table("GSE.txt")[[1]])

##### collect information to summarise all datasets and double check by inclusion criteria 

# prepare the SQL database of metadata
{
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
db = src_sqlite('GEOmetadb.sqlite')
GSE = tbl(db,'gse')
GSE_GSM <- tbl(db, 'gse_gsm')
GSM <- tbl(db,'gsm')
GSE_GPL <- tbl(db, 'gse_gpl')
GPL <- tbl(db,'gpl')
}

# collect the original text describing source tissue of each GSM and remove redundancy

CollectTissue <- function(gse_numbers){
  # gse_numbers can also be only one GSE number
  
  tissue_names <- c()
  for(gse_num in gse_numbers){  
    gse_certain_gsm <- filter(GSE_GSM, gse == gse_num) 
    gse_certain_gsm <- as_tibble(gse_certain_gsm)
    gsm_numbers <- gse_certain_gsm$gsm
  
    for(gsm_num in gsm_numbers){
      gsm_certain <- filter(GSM, gsm == gsm_num)
      gsm_certain <- as_tibble(gsm_certain)
      if((length(gsm_certain$source_name_ch1) != 0) && !(gsm_certain$source_name_ch1 %in% tissue_names)){
        tissue_names <- c(tissue_names, gsm_certain$source_name_ch1)
      }
    }
  }
  return(tissue_names)
}

tissue_names <- CollectTissue(gse_numbers)

write(tissue_names, "tissue_names.xls")

# check the text file "tissue_names.xls" and further remove redundancy using the function below

ConvertTissueName <- function(name){
  # convert the original tissue's name as are shown in source_name_ch1 to hash key
  # hash key can also act as a relatively 'clean' tissue name for the table in the final report
  
  if(str_detect(name, "^brain,")){  # after "brain,", before the next ","
    tmp <- gsub("^brain,\\s*([^,]+),?.*$", "\\1", name)
  }
  else if(str_detect(name, "of")){
    tmp <- gsub("(^.+)\\sof.+", "\\1", name)
  }
  else if(str_detect(name, ".+_\\S*male_ad")){
    tmp <- gsub("(.+)_\\S*male_ad", "\\1", name)
  }
  else if(str_detect(name, "_[A-Z][A-Z]_")){
    tmp <- gsub(".+_([A-Z][A-Z])_.+", "\\1", name)
  }
  else if(str_detect(name, "_")){
    tmp <- gsub("(^[^_]+)_.+", "\\1", name)
  }
  else if(str_detect(name, "frozen")){
    tmp <- gsub(".+frozen\\s(.+)", "\\1", name)
  }
  else if(str_detect(name, ",")){
    tmp <- gsub("(^[^,]+),.+", "\\1", name)
  }
  else {
    tmp <- name
  }
  tmp <- tolower(tmp)
  return(tmp)
}

# read in the file "tissue_names.xls", 
# and create a file "tissue_name_less.txt" with even less redundancy
tissue_names <- read_tsv("tissue_names.xls", col_names = FALSE)[[1]]
tissue_name_less <- c()
tmp <- ""
for(name in tissue_names){
  tmp <- ConvertTissueName(name)  # the function is defined above
  tissue_name_less <- c(tissue_name_less, tmp)
}
tissue_name_less <- unique.default(tissue_name_less)
write(tissue_name_less, "tissue_name_less.txt")

# manually edit the file 'tissue_name_less.txt' to match each tissue name
# to 1 or 2 brain regions of interest that the tissue belongs to.
# add column names as 'Original_tissue_name',	'Belongs_to_version1' and	'Belongs_to_version2'
# save as a new file "tissue_name_mapping.txt"

# read in the file 'tissue_name_mapping.txt' and make a hash table 
# with the first and second column as keys and values respectively
old_new <- read_tsv("tissue_names_mapping.txt")
old.tissue_new.tissue <- hash(old_new$Original_tissue_name, old_new$Belongs_to_version1)


# extract information for each GSE from the tables in the database
datasetInfo <- tibble()
for(gse_num in gse_numbers){
  gse_certain <- as_tibble(filter(GSE, gse == gse_num))
  if(dim(gse_certain)[1] == 0){  # this GSE has not been recorded
    gse_certain <- as_tibble(t(as_tibble(rep(NA, 20))))
    colnames(gse_certain) <- c(colnames(GSE), "original_tissue", "mapping_tissue")
    gse_certain$gse <- gse_num
  }
  else {  # add tissue information
    tissue_names <- CollectTissue(gse_num)  # function defined above
    tissue_names <- unique.default(tissue_names)
    original_tissue <- c()
    mapping_tissue <- c()
    for(tissue in tissue_names){
      clean_tissue <- ConvertTissueName(tissue)  # function defined above
      original_tissue <- c(original_tissue, clean_tissue)
      mapping_tissue <- c(mapping_tissue, old.tissue_new.tissue[[clean_tissue]])
    }
    original_tissue <- paste0(unique.default(original_tissue), collapse = "; ") 
    mapping_tissue <- paste0(unique.default(mapping_tissue), collapse = "; ")
    twoNewColumns <- tibble(original_tissue = original_tissue, 
                           mapping_tissue = mapping_tissue)
    
    gse_certain <- bind_cols(gse_certain, twoNewColumns)
  }
  
  gpl_certain <- as_tibble(filter(GSE_GPL, gse == gse_num))
  if(dim(gpl_certain)[1] == 0){  # no GPL info found for this GSE
    gpl_certain_detail <- as_tibble(t(as_tibble(rep(NA, 20))))
    colnames(gpl_certain_detail) <- colnames(GPL)
    gse_certain <- bind_cols(gse_certain, gpl_certain_detail)
  }
  else{
    gpl_certain_detail <- as_tibble(filter(GPL, gpl == !!gpl_certain$gpl[1]))
    gse_certain <- bind_cols(gse_certain, gpl_certain_detail)
    if(dim(gpl_certain)[1] > 1){  # this gse has more than one gpl
      for(n in 2:dim(gpl_certain)[1]){
        tmp <- as_tibble(t(as_tibble(rep(NA, 20))))
        colnames(tmp) <- c(colnames(GSE), "original_tissue", "mapping_tissue")
        gpl_certain_detail <- as_tibble(filter(GPL, gpl == !!gpl_certain$gpl[n]))
        tmp <- bind_cols(tmp, gpl_certain_detail)
        gse_certain <- bind_rows(gse_certain, tmp)
      }
    }
  }
  
  datasetInfo <- bind_rows(datasetInfo, gse_certain)
}

datasetInfoFinal <- dplyr::select(datasetInfo, title, gse, submission_date, last_update_date,
                           pubmed_id, original_tissue, mapping_tissue, summary, type, 
                           overall_design, gpl, platform = title1,
                           manufacturer, data_row_count, bioc_package)

# write the table to a file and double check each GSE for inclusion criteria
write_tsv(datasetInfoFinal, "dataset_info.xls")
# GSE118553 and GSE132903 is too new to be included by the GEOmetadb.sqlite
# GSE131617, GSE122063 and GSE125583 have no tissue information included
# the missing information need to be added manually

# edit the file "dataset_info.xls" to create a table for the report's main text


##### collecting, fixing and annotating data from selected GSE

crossmeta::get_raw(gse_numbers, getwd())
# no supplementary files for GSE15222, GSE4226, GSE4227, GSE4229 (-6, -17, -18, -19). 
# need to use GEOquery::getGEO() to create the eset from matrix

# the following GSE are generated by illumina, the column name of their table 
# was edited using NotePad++'s regular expression

# GSE29378
# find what: (\S+)\tDetection Pval
# replace with: AVG_Signal-\1\tDetection-\1

# GSE26927
# find what: (SAMPLE \d+)\tDetection Pval
# replace with: AVG_Signal-\1\tDetection-\1

# GSE118553
# find what: (\S+)\t\1_Detection_Pval
# replace with: AVG_Signal-\1\tDetection-\1

# GSE63060
# find what: \t(\S+)-AVG_Signal
# replace with: \tAVG_Signal-\1
# no dection p value present

# GSE63061
# find what: \t(\d{10}_[A-Z])
# replace with: \tAVG_Signal-\1
# no dection p value present

# loading and annotating data by crossmeta
esets1 <- crossmeta::load_raw(gse_numbers[c(-6, -17, -18, -19, -16)])

# some packages cannot be installed automatically: 
# pd.hg.u133.plus.2, illuminaHumanv3.db, hgu133plus2.db, 
# huex10sttranscriptcluster.db, hgu133a.db, illuminaHumanv4.db, 
# hgu133b.db, pd.hg.u133b, pd.hg.u133a, pd.hugene.1.0.st.v1, 
# pd.hugene.1.1.st.v1, pd.hugene.2.0.st, pd.hugene.2.1.st, 
# hugene10sttranscriptcluster.db, hugene20sttranscriptcluster.db
# install them manually

# GSE45596 is a dual-channel array which is not currently supported by crossmeta (-16)

# platform type not supported be crossmeta: 
# GSE122063, GSE125583, GSE44771, GSE44768, GSE44770, GSE33000
# GSE125583 is a high-throughput sequencing dataset (-23)
# (-22, -23, -31, -32, -33, -7)

# Affymatrix data that cannot be loaded by crossmeta: GSE110226 (-26)

# illumina data that cannot be loaded by crossmeta: 
# GSE26927, GSE118553 (-2, -13), 
# GSE63060, GSE63061, GSE132903 are laking detection p value (-14, -15, -20)
# Agilent data that cannot be loaded by crossmeta: GSE85426 (-30)

# the eset object of each GSE is stored as .rds file in the folder of corresponding GSE
# for each GSE that is successfully loaded and annotated. 
# the .rds file will be re-use for next loading

# re-load datasets that are actually loadable by crossmeta
gse_numbers_crossmeta <- gse_numbers[c(-6, -17, -18, -19, -16, -22, -23, -31, -32, -33, 
                                       -7, -26, -2, -13, -14, -15, -20, -30)]
esets1 <- crossmeta::load_raw(gse_numbers_crossmeta)


# for each microaeeay dataset that cannot be processed by crossmeta,
# get and save the ExpressionSet directly from GEO using GEOquery
gse_numbers_not_crossmeta <- gse_numbers[c(6, 17, 18, 19, 16, 22, 31, 32, 33, 
                                           7, 26, 2, 13, 14, 15, 20, 30)]

GetExpressionSets <- function(gse_numbers_not_crossmeta){
  # a function to get, annotate and save ExpressionSet for given GSEs
  
  for(gse_num in gse_numbers_not_crossmeta){
    dir.create(gse_num)
    eset <- GEOquery::getGEO(gse_num, AnnotGPL = TRUE, destdir = gse_num)
    names(eset) <- gse_num
    fData(eset[[1]]) <- dplyr::select(fData(eset[[1]]), ENTREZID = `Gene ID`, 
                                      PROBE = ID, SYMBOL = `Gene symbol`)
    saveRDS(eset, file = paste0("./", gse_num, "/", gse_num, "_eset.rds"))
  }
}

# GetExpressionSets(gse_numbers_not_crossmeta)  # error

# GSE122063 (-22) and GSE85426 (-30) are lacking fData and/or pData
# skip the 2 and try the rest

gse_numbers_not_crossmeta <- gse_numbers[c(31, 32, 33, 7, 26, 2, 13, 14, 15, 20)]
GetExpressionSets(gse_numbers_not_crossmeta)

# all GSEs that are loaded by GEOquery instead of crossmeta
gse_numbers_not_crossmeta <- gse_numbers[c(6, 17, 18, 19, 16, 31, 32, 33, 
                                           7, 26, 2, 13, 14, 15, 20)]
esets2 <- crossmeta::load_raw(gse_numbers_not_crossmeta)  # load the above prepared esets


##### summary of datasets loaded

# esets1 includes 15 GSEs processed by crossmeta;
# (GSE131617 has 2 GPL, GSE84422 has 3, so 18 elements)
# esets2 include 15 GSEs with processed expression data directly get from GEO;
# 3 GSEs are discarded: GSE122063, GSE85426 and GSE125583
# 4 GSEs have no tissues information in GEOmetadb.sqlite, 
# but such info can be collected from esets: 
# GSE33000, GSE118553, GSE132903, GSE131617

# check the 4 datasets lacking tissue info in GEOmetadb.sqlite 
# and modify the ConvertTissueName() function if necessary
{
View(pData(esets1$GSE131617.GPL5175))
View(pData(esets1$GSE131617.GPL5188))
# the 2 have tissue info in characterstics_ch1.2, as "tissue region: Frontal cortex (FC)"/ 
# "tissue region: Entorhinal cortex (EC)" / "tissue region: Temporal cortex (TC)"

View(pData(esets2$GSE33000))
# 2 channels, source_name_ch1 is as "HBTRC_PF_Pool_1" while source_name_ch2 is as PFC_1

View(pData(esets2$GSE118553))
# source_name_ch1 contains "Temporal_Cortex" / "Cerebellum" 
# / "Frontal_Cortex" / "Entorhinal_Cortex"

View(pData(esets2$GSE132903))
# source_name_ch1 contain middle_temporal_gyrus

}


##### adding group information according to AD-control and tissues


# show all GSEs that include:

# hippocampus 
datasetInfoFinal %>%
  filter(str_detect(mapping_tissue, "hippocampus")) %>%
  dplyr::select(gse)

# GSE5281  
# GSE48350 
# GSE36980 
# GSE28146 
# GSE29378 
# GSE1297  
# GSE110298
# GSE84422 
# none of them are lacking tissue info in the GEOmetadb.sqlite
# all of them are in the esets1 (loaded by crossmeta)

# cingulate cortex
{
datasetInfoFinal %>%
  filter(str_detect(mapping_tissue, "cingulate cortex")) %>%
  dplyr::select(gse)

# GSE5281 
# GSE26927
# GSE84422
}
# entorhinal cortex

# extract the AD-Ctrl info and add the group column one by one for each GSE
{
  View(pData(esets1$GSE5281))
  group <- gsub(".+control.+", "ctrl", pData(esets1$GSE5281)$title)
  group <- gsub(".+affected.+", "AD", group)
  group <- make.names(group)
  pData(esets1$GSE5281)$group <- group
  
  for(i in 1:length(group)){  # remove those under 65 in the ctrl group
    tmp <- gsub(".ge: \\D*(\\S+) [a-z]+$", "\\1", esets1$GSE5281$characteristics_ch1.11[i])
    tmp <- as.numeric(tmp)
    if(tmp < 65){
      esets1$GSE5281$group[i] <- "young"
    }
  }
  
  
  View(pData(esets1$GSE48350))
  # AD patients are all late-onset, healthy control include both young and old individuals
  group <- gsub(".+_AD$", "AD", pData(esets1$GSE48350)$source_name_ch1)
  group <- gsub(".+,\\s*(\\d+)\\s*years$", "\\1", group)
  
  for(i in 1:length(group)){  # only set healthy elderly as control
    if(str_detect(group[i], "\\d+")){
      t <- as.numeric(group[i])
      if(t > 65){
        group[i] <- "ctrl"
      }
      else{group[i] <- "young"}
    }
  }
  group <- make.names(group)
  pData(esets1$GSE48350)$group <- group
  
  
  View(pData(esets1$GSE36980))
  group <- gsub("^AD_.+$", "AD", pData(esets1$GSE36980)$title)
  group <- gsub("^non-AD_.+$", "ctrl", group)
  group <- make.names(group)
  pData(esets1$GSE36980)$group <- group
  
  for(i in 1:length(group)){  # remove those under 65 in the ctrl group
    tmp <- gsub("age: (\\S+)$", "\\1", esets1$GSE36980$characteristics_ch1.1[i])
    tmp <- as.numeric(tmp)
    if(tmp < 65){
      esets1$GSE36980$group[i] <- "young"
    }
  }
  
  
  View(pData(esets1$GSE28146))
  # no one under 65, but only use moderate and severe patients as cases
  group <- gsub("^Control.+$", "ctrl", pData(esets1$GSE28146)$title)
  group <- gsub("^Moderate.+$", "AD", group)
  group <- gsub("^Severe.+$", "AD", group)
  group <- gsub("^Incipient.+$", "out", group)
  group <- make.names(group)
  pData(esets1$GSE28146)$group <- group
  
  
  View(pData(esets1$GSE29378))
  
  group <- gsub(".*Control.*", "ctrl", pData(esets1$GSE29378)$title)
  group <- gsub(".*AD.*", "AD", group)
  group <- make.names(group)
  pData(esets1$GSE29378)$group <- group
  for(i in 1:length(group)){  # remove those under 65 from AD or control group
    if(as.numeric(esets1$GSE29378$`age(years):ch1`[i]) < 65){
      esets1$GSE29378$group[i] <- "young"
    }
  }
  
  
  View(pData(esets1$GSE1297))
  # no one under 65, but only use moderate and severe patients as cases
  group <- gsub("^Control.+$", "ctrl", pData(esets1$GSE1297)$title)
  group <- gsub("^Moderate.+$", "AD", group)
  group <- gsub("^Severe.+$", "AD", group)
  group <- gsub("^Incipient.+$", "out", group)
  group <- make.names(group)
  pData(esets1$GSE1297)$group <- group
  
  
  View(pData(esets1$GSE110298))
  # cannot figure out AD and control
  
  View(pData(esets1$GSE84422.GPL570))
  View(pData(esets1$GSE84422.GPL96))
  View(pData(esets1$GSE84422.GPL97))  # only this include hippocampus sample
  
  group <- pData(esets1$GSE84422.GPL97)$`neuropathological category:ch1`
  group <- gsub("definite AD", "AD", group)
  group <- gsub("Probable AD", "AD", group)
  group <- gsub("Normal", "ctrl", group)
  group <- gsub("Possible AD", "out", group)  # only include definite and probable AD as cases
  
  for(i in 1:length(group)){
    if(as.numeric(esets1$GSE84422.GPL97$`age:ch1`[i]) < 65){
      group[i] <- "young"
    }
  }
  
  group <- make.names(group)
  pData(esets1$GSE84422.GPL97)$group <- group
  
}

esets_hippocampus <- c(esets1$GSE5281, esets1$GSE48350, esets1$GSE36980, esets1$GSE28146, 
                       esets1$GSE29378, esets1$GSE1297, esets1$GSE84422.GPL97)
names(esets_hippocampus) <- c("GSE5281", "GSE48350", "GSE36980", "GSE28146", 
                              "GSE29378", "GSE1297", "GSE84422")


# extract tissue info and append it to the group column for all GSE
AddTissue <- function(old.tissue_new.tissue, ConvertTissueName, esets){
  # input the tissue conversion hash table, the conversion function 
  #and a list of ExpressionSet objects

  for(g in 1:length(esets)){  # if use gse in esets, the group column won't change
    keys <- ConvertTissueName(pData(esets[[g]])$source_name_ch1)
    
    tissue <- c()
    for(n in 1:length(keys)){
      if(length(old.tissue_new.tissue[[keys[n]]]) == 0){
        keys[n] <- gsub("(.+)_\\S*male_ad", "\\1", keys[n])
      }
      tissue[n] <- old.tissue_new.tissue[[keys[n]]]
    }
    group_tissue <- paste(pData(esets[[g]])$group, tissue)
    group_tissue <- make.names(group_tissue)
    pData(esets[[g]])$group <- group_tissue
  }
  return(esets)
}

esets_hippocampus <- AddTissue(old.tissue_new.tissue, ConvertTissueName, esets_hippocampus)
# saveRDS(esets_hippocampus, file = "esets_hippocampus_before_QC.rds")  # for re-use






##### try quality control (QC), including PCA, heatmap, hierarchical clustering tree and sva comBat 

esets_hippocampus <- readRDS("esets_hippocampus_before_QC.rds")

table(str_detect(esets_hippocampus$GSE5281$group, "AD.hippocampus"))  # check the size of AD group
length(unique.default(na.omit(fData(esets_hippocampus$GSE5281)$ENTREZID)))  #  count the number of genes 
# many probes cannot be annotated, many genes have more than one probe

# subset an eset
esets_test <- esets_hippocampus$GSE5281[, (esets_hippocampus$GSE5281$group == "AD.hippocampus" 
                                           | esets_hippocampus$GSE5281$group == "ctrl.hippocampus")]

GroupDiffVis <- function(esets_test, gse){
  # check the difference between AD and control group
  esets_test <- esets_test[, (esets_test$group == "AD.hippocampus" 
                              | esets_test$group == "ctrl.hippocampus")]
  
  # hierarchical clustering tree
  pdf(file = paste0(gse, "_tree.pdf"))
  expr.matr <- exprs(esets_test)
  expr.matr.t <- t(expr.matr)
  rownames(expr.matr.t) <- paste0(esets_test$group, ".", esets_test$geo_accession)
  eu.dist <- dist(expr.matr.t, method = "euclidean")
  tree <- hclust(eu.dist, method = "average")
  plot(tree, main = paste0("Hierachical Clustering Tree for ", gse), xlab = "Sample", sub = "")
  dev.off()
  
  # PCA
  expr.matr.t.ti <- as_tibble(expr.matr.t)
  expr.pca <- prcomp(expr.matr.t.ti, center = TRUE, scale. = FALSE)
  
  expr.pca.plot <- expr.pca$x %>% 
    as_tibble() %>%
    bind_cols(group = esets_test$group, sample = esets_test$geo_accession)  
    # prepare to plot the sameple names and experimental groups
  
  ggplot(expr.pca.plot, aes(x = PC1, y = PC2)) + 
    geom_point(mapping = aes(colour = group)) +
    geom_text(mapping = aes(label = sample), size = 3) +
    labs(title = paste0("PCA Plot for ", gse), 
         x = "PC1", y =  "PC2") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(filename = paste0(gse, "_pca.pdf"))
  
  # heatmap (something wrong, extremely time-consuming and computational resource consuming)
  # pheatmap(expr.matr, show_rownames = FALSE, fontsize_row = 5, main = "Heatmap")
  
}

dir.create("pca_and_tree_plots")
setwd("pca_and_tree_plots")
{
GroupDiffVis(esets_hippocampus$GSE5281, "GSE5281")
GroupDiffVis(esets_hippocampus$GSE48350, "GSE48350")
GroupDiffVis(esets_hippocampus$GSE36980, "GSE36980")
GroupDiffVis(esets_hippocampus$GSE28146, "GSE28146")
GroupDiffVis(esets_hippocampus$GSE29378, "GSE29378")
GroupDiffVis(esets_hippocampus$GSE1297, "GSE1297")
GroupDiffVis(esets_hippocampus$GSE84422, "GSE84422")
}
setwd("..")

# remove batch effect!!

# try gender & age on GSE48350
esets_test <- esets_hippocampus$GSE48350
esets_test <- esets_test[, (esets_test$group == "AD.hippocampus" 
                            | esets_test$group == "ctrl.hippocampus")]


# A. limma's removeBatchEffect

exprs(esets_test) <- removeBatchEffect(x = exprs(esets_test), 
                                       batch = esets_test$`gender:ch1`, 
                                       covariates = as.numeric(esets_test$`age (yrs):ch1`), 
                                       design = matrix(data = as.numeric(gsub(pattern = "ctrl.hippocampus",
                                                                              replacement = 2,
                                                                              x = gsub("AD.hippocampus", 1,
                                                                                       esets_test$group)), 
                                                                         nrow = ncol(exprs(esets_test)), 
                                                                         ncol = 1)))

                  
# B. sva's comBat (does not work)

group <- factor(esets_test$group)
group <- relevel(group, "ctrl.hippocampus")
design <- model.matrix(~group)

exprs(esets_test) <- sva::ComBat(dat = exprs(esets_test),
                                 batch = esets_test$`gender:ch1`, 
                                 # covariates = as.numeric(esets_test$`age (yrs):ch1`), 
                                 mod = design)



# try exploBATCH on GSE48350

# test the package
# require(exploBATCHbreast)
# data(Breast)  # column - gene, row - sample
data(batchBreast)
expBATCH(D=Breast,batchCL=batchBreast,
         Conf=NA,mindim=2,maxdim=9,method="ppcca",SDselect=2)

require(exploBATCH)
require(exploBATCHcolon)   
data(Colon)                       
data(batchColon)            
data(bioCL)

expBATCH(
  D=Colon,                 
  batchCL=batchColon,     
  Conf=bioCL,              # this is the variable defining the biological effect of interest
  mindim=2,                
  maxdim=11,                  
  method="ppcca",          
  SDselect=1.2             # set SD at 2 to reduce computational time.
)


# test expoBATCH on GSE48350

group <- factor(esets_test$group)
levels(group) <- c(1, 2)
group <- as.numeric(group)
exprMatr <- t(exprs(esets_test))
exploBATCH::expBATCH(
  D = exprMatr,                 
  batchCL = rep(1, length(group)),     
  Conf = group,           
  mindim = 2,                
  maxdim = 3,                  
  method = "ppcca",          
  SDselect = 1
)


# try downloading the processed data ready for DE analysis from GEO directly
# GSE48350
esets_test <- GEOquery::getGEO("GSE48350", AnnotGPL = TRUE)
names(esets_test) <- "GSE48350"
fData(esets_test[[1]]) <- dplyr::select(fData(esets_test[[1]]), ENTREZID = `Gene ID`, 
                                        PROBE = ID, SYMBOL = `Gene symbol`)

group <- gsub(".+_AD$", "AD", pData(esets_test$GSE48350)$source_name_ch1)
group <- gsub(".+,\\s*(\\d+)\\s*years$", "\\1", group)

for(i in 1:length(group)){  # only set healthy elderly as control
  if(str_detect(group[i], "\\d+")){
    t <- as.numeric(group[i])
    if(t > 65){
      group[i] <- "ctrl"
    }
    else{group[i] <- "young"}
  }
}
group <- make.names(group)
pData(esets_test$GSE48350)$group <- group
esets_test <- AddTissue(old.tissue_new.tissue, ConvertTissueName, esets_test)


# GSE84422
esets_test <- GEOquery::getGEO("GSE84422", AnnotGPL = TRUE)
group <- pData(esets_test$`GSE84422-GPL97_series_matrix.txt.gz`)$`neuropathological category:ch1`
group <- gsub("definite AD", "AD", group)
group <- gsub("Probable AD", "AD", group)
group <- gsub("Normal", "ctrl", group)
group <- gsub("Possible AD", "out", group)  # only include definite and probable AD as cases

for(i in 1:length(group)){
  if(as.numeric(esets_test$`GSE84422-GPL97_series_matrix.txt.gz`$`age:ch1`[i]) < 65){
    group[i] <- "young"
  }
}

group <- make.names(group)
pData(esets_test$`GSE84422-GPL97_series_matrix.txt.gz`)$group <- group
esets_test <- AddTissue(old.tissue_new.tissue, ConvertTissueName, esets_test[3])
names(esets_test) <- "GSE84422"
fData(esets_test[[1]]) <- dplyr::select(fData(esets_test[[1]]), ENTREZID = `Gene ID`, 
                                        PROBE = ID, SYMBOL = `Gene symbol`)



# try sva's sva() on GSE84422 & visualisation
esets_test <- esets_hippocampus$GSE84422[, (esets_hippocampus$GSE84422$group == "AD.hippocampus" 
                                           | esets_hippocampus$GSE84422$group == "ctrl.hippocampus")]
esets_test <- list("GSE84422" = esets_test)

pheno <- pData(esets_test)
edata <- exprs(esets_test)
mod <- model.matrix(~as.factor(group), data = pheno)
mod0 = model.matrix(~1, data = pheno)

n.sv = num.sv(edata, mod, 
              method = "leek")  # if use default "method = c("be", "leek")", then n.sv = 6
svobj = sva(edata, mod, mod0, n.sv = n.sv)

n.pc <- n.sv
edata <- sva_network(edata, 6)  # does not work, cannot be used to replace 
exprs(esets_test$GSE84422) <- edata

svaBatchCor <- function(dat, mmi, mm0, n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  library(sva)
  if(is.null(n.sv)){n.sv <- num.sv(dat, mmi, method = "leek")}
  o <- sva(dat, mmi, mm0, n.sv = n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
}

lnj.corr <- svaBatchCor(edata, mod, mod0, n.sv = 1)
co <- lnj.corr$corrected
exprs(esets_test$GSE84422) <- co

# check the result

GroupDiffVis(esets_test$GSE48350, "GSE48350")
GroupDiffVis(esets_test$GSE84422, "GSE84422")

DEanalysis(esets_test[1], "hippocampus")


##### differential expression analysis using crossmeta

DEanalysis <- function(esets, tissue){
  # must input a list of ExpressionSet objects, not an ExpressionSet
  # only for one tissue
  tissue <- make.names(tissue)
  contrasts = paste0("AD.", tissue, "-", "ctrl.", tissue)
  #anals <- c()
  for(i in length(esets)){
    sel <- crossmeta::setup_prev(esets[i], contrasts)
    anal <- crossmeta::diff_expr(esets[i], getwd(), prev_anals = sel, annot = "ENTREZID")
    # anals <- c(anals, anal)
  }
  # return(anals)
}

# anals_hippocampus <- DEanalysis(esets_hippocampus, "hippocampus")
{
DEanalysis(esets_hippocampus[1], "hippocampus")
DEanalysis(esets_hippocampus[2], "hippocampus")
DEanalysis(esets_hippocampus[3], "hippocampus")
DEanalysis(esets_hippocampus[4], "hippocampus")
DEanalysis(esets_hippocampus[5], "hippocampus")
DEanalysis(esets_hippocampus[6], "hippocampus")
DEanalysis(esets_hippocampus[7], "hippocampus")
}

# GSE5281_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 8443 
# GSE48350_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 0 
# GSE36980_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 2 
# GSE28146_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 0 
# GSE29378_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 451 
# GSE1297_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 0 
# GSE84422_AD.hippocampus-ctrl.hippocampus (# p < 0.05): 0 

gse_hippocampus <- c("GSE5281", "GSE48350", "GSE36980", "GSE28146",
                     "GSE29378", "GSE1297", "GSE84422")
# re-use for the next time if not overwritten by DE analysis for other tissues
anals_hippocampus <- crossmeta::load_diff(gse_hippocampus, getwd(), annot = "ENTREZID")


##### meta-analysis using crossmeta

es_res_hippocampus <- crossmeta::es_meta(anals_hippocampus)

PrepareIntput <- function(es_res){
  
  # preparing the inputs for clusterProfiler, pathview and so on
  es_res$all$filt %>% 
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    dplyr::select(gene, mu, var, z, pval, fdr) -> InPutForCP

  # filter the DEGs by p-value and false discover rate
  InPutForCP <- filter(InPutForCP, pval < 0.05)  #  && fdr < 0.25

  # prepare a named vector of effect size (z-score) for GSEA and pathview
  # names are corresponding gene ID
  GeneNamedVector <- InPutForCP$z
  names(GeneNamedVector) <- InPutForCP$gene
  GeneNamedVector <- sort(GeneNamedVector, decreasing=TRUE)
  return(list( one = InPutForCP, two = GeneNamedVector))
}

tmp <- PrepareIntput(es_res_hippocampus)
InPutForCP_hippo <- tmp$one
GeneNamedVector_hippo <- tmp$two

View(InPutForCP_hippo)  
# 865 genes have a combined p-value above 0.05
# but all fdr is 1

# may need to try other DE analysis R packages or use RNA-Seq



##### overrepresentation analysis using clusterProfiler (need to be devided into up and down!)

GO <- clusterProfiler::enrichGO(InPutForCP_hippo$gene, org.Hs.eg.db)
KEGGpathway <- clusterProfiler::enrichKEGG(InPutForCP_hippo$gene, keyType = "ncbi-geneid")
KEGGmodule <- clusterProfiler::enrichMKEGG(InPutForCP_hippo$gene, keyType = "ncbi-geneid")


##### GSEA using clusterProfiler (need to be devided into up and down!)

GO.gsea <- clusterProfiler::gseGO(geneList = GeneNamedVector_hippo, 
                                  OrgDb = org.Hs.eg.db)
KEGGpathway.gsea <- clusterProfiler::gseKEGG(GeneNamedVector_hippo, 
                                                keyType = "ncbi-geneid")
KEGGmodule.gsea <- clusterProfiler::gseMKEGG(GeneNamedVector_hippo, 
                                                keyType = "ncbi-geneid")


##### using pathview to visualise enriched pathways that are related to the 3 keywords

pv.out.mitophagy <- pathview(gene.data = GeneNamedVector_hippo,
                             gene.idtype = "entrez", 
                             kegg.dir = "./pathview",  # provide the xml file by myself in this dir
                             pathway.id = "04137", species = "hsa",
                             same.layer = F, kegg.native = T, 
                             out.suffix = "hippocampus")


