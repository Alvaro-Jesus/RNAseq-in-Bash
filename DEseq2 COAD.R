## Loading packages
pacman::p_load(TCGAbiolinks,tidyverse,maftools,pheatmap,SummarizedExperiment)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('CPTAC-2')

summary(gdcprojects)

# building a query
query_TCGA <- GDCquery(project = 'CPTAC-2',
                       data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-COAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')

getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)


# prepare data
tcga_coad_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
coad_matrix <- assay(tcga_coad_data, 'unstranded')

write.table(coad_matrix,file="TCGA_COAD_counts.txt", quote=FALSE, sep="\t")
TCGA_COAD_counts <- read.table("TCGA_COAD_counts.txt", sep="\t", header=TRUE, quote = "", check.names = FALSE)


getwd()

samples_tumor <- TCGAquery_SampleTypes(barcode = colnames(coad_matrix), typesample = "TP")
samples_normal <- TCGAquery_SampleTypes(barcode = colnames(coad_matrix), typesample = "NT")

#Dataframe con los diferentes tipos de muestra Normal y tumor
df <- data.frame(
  barcode = colnames(coad_matrix),
  type = case_when(
    substr(colnames(coad_matrix), 14, 15) == "01" ~ "Tumor",
    substr(colnames(coad_matrix), 14, 15) == "02" ~ "Tumor", #Recurrent Tumor
    substr(colnames(coad_matrix), 14, 15) == "06" ~ "Tumor", #Metastasis
    substr(colnames(coad_matrix), 14, 15) == "11" ~ "Normal",
    substr(colnames(coad_matrix), 14, 15) == "12" ~ "Blood Normal",
    TRUE ~ "Other"
  )
)


# Saving colDATA

write.table(df,file="TCGA_COAD_phenotype.txt", quote=FALSE, sep="\t")
TCGA_COAD_phenotype <- read.table("TCGA_COAD_phenotype.txt", sep="\t", header=TRUE, quote = "", check.names = FALSE)


unique(df$type)
table(df$type)

rownames(df) <- df[,1]

#------------------#
#-------DEG--------#
#------------------#


# load libraries
pacman::p_load(DESeq2,tidyverse,airway)
# Step 1: preparing count data ----------------

# read in counts data
TCGA_COAD_counts <- read.table("TCGA_COAD_counts.txt", sep="\t", header=TRUE, quote = "", check.names = FALSE)
head(TCGA_COAD_counts)


# read in sample info
TCGA_COAD_phenotype <- read.table("TCGA_COAD_phenotype.txt", sep="\t", header=TRUE, quote = "", check.names = FALSE)
unique(TCGA_COAD_phenotype$type)
table(TCGA_COAD_phenotype$type)

# making sure the row names in colData matches to column names in counts_data
all(colnames(TCGA_COAD_counts) %in% rownames(TCGA_COAD_phenotype))

# are they in the same order?
all(colnames(TCGA_COAD_counts) == rownames(TCGA_COAD_phenotype))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = TCGA_COAD_counts,
                              colData = TCGA_COAD_phenotype,
                              design = ~ type)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$type <- relevel(dds$type, ref = "Normal")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)



# MA plot
plotMA(res)


#------------------#
#-------Volcano Plot--------#
#------------------#


pacman::p_load(EnhancedVolcano)
library(org.Hs.eg.db)
rownames(res.df) <- gsub("\\..*", "", rownames(res.df))

#Pasando ENSEMBL a SYMBOL
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(org.Hs.eg.db, keys=rownames(res.df), keytype= "ENSEMBL", column="SYMBOL")


setwd("C:/Users/user/OneDrive - Kagoshima University (1)/CRC/CRC/TCGA_DEG")

write.table(res.df,file="TCGA_COAD_DEG_tumorvsnormal.txt", quote=FALSE, sep="\t")

TCGA_COAD_DEG_tumorvsnormal <- read.table("TCGA_COAD_DEG_tumorvsnormal.txt", sep="\t", header=TRUE, quote = "", check.names = TRUE)
res.df <- TCGA_COAD_DEG_tumorvsnormal
##
str(res.df)
EnhancedVolcano(res.df, x= "log2FoldChange", y="padj", lab= res.df$symbol)

############### CPTAC ########## 

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-COAD')

summary(gdcprojects)

# building a query
query_TCGA <- GDCquery(project = 'TCGA-COAD',
                       data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-COAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')

getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)


# prepare data
tcga_coad_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
coad_matrix <- assay(tcga_coad_data, 'unstranded')




