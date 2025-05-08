#Author: Dr. Clara Apicella
#Script to perform deconvolution of bulk RNA-seq data from murine hearts with MuSic
#https://xuranw.github.io/MuSiC/articles/MuSiC.html

#load libraries 
#library(Seurat)
#library(SeuratDisk)
library(MuSiC)
library(openxlsx)
library(Biobase)
library(dplyr)
library(SingleCellExperiment)

#Set seed for reproducibility (important for any random operation later)
set.seed(1234)

##----  Read bulk dataset and reference datasets

#set working directory
#wdir='W:/01_GenEpi_Projects/Schulte_2025/6.Deconvolution/Music'
#setwd(wdir)

#Load reference data from snRNA-seq from human study 'Tucker et al 2020 - Transcriptional and Cellular diversity of the human heart'
#Load data in h5ad format (from scanpy)
# Convert("W:/01_GenEpi_Projects/Schulte_2025/6.Deconvolution/2020_TranscriptionalandCellularDiversityOfTheHumanHeart/data_from_broadInstitute/healthy_human_4chamber_map_unnormalized_V4.h5ad", dest="h5Seurat", overwrite = FALSE)
#load h5Seurat file into a Seurat object
# Tucker2020_Seuratv4 <-LoadH5Seurat("W:/01_GenEpi_Projects/Schulte_2025/6.Deconvolution/2020_TranscriptionalandCellularDiversityOfTheHumanHeart/data_from_broadInstitute/healthy_human_4chamber_map_unnormalized_V4.h5seurat")
# str(Tucker2020_Seuratv4)

# Tucker2020_Seuratv4@assays$RNA$data


#--- Create the ExpressionSet object
#https://www.bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf 


#Load bulk RNA raw data and covariates
bulk_groups <-read.table("Groups.txt", header=TRUE, row.names = 1)
bulk_counts <- read.table("CountMatrix.txt", header=TRUE, row.names = 1)
#Check that samples are in the same order 
all(rownames(bulk_groups)==colnames(bulk_counts))

library(Biobase)
#create matrix from count data.frame
bulk_counts_mtx <- as.matrix(bulk_counts)

#create AnnotatedDataFrame object from covariates data.frame for phenoData
#Determine colnames covariates dataframe
#These will be the rownames of the metadata dataframe
names(bulk_groups)
#Create metadata dataframe to describe the covariates, with a single column "labelDescription"
metadata <- data.frame(labelDescription=
                         c("Sex+Genotype group combination",
                           "Sample Sex",
                           "Genotype group"),
                       row.names = c("Condition", "Sex","Group"))
#create the phenoData annotated dataframe
phenoData <- new("AnnotatedDataFrame", data = bulk_groups, varMetadata= metadata)
head(phenoData)
#to extract the original phenotypic data use the function pData()
head(pData(phenoData))
#to extract the variables metadata varMetadata()
head(varMetadata(phenoData))



#Generate test object
bulk_obj_test <- ExpressionSet(assayData = bulk_counts_mtx[,1:2], phenoData = phenoData[1:2,])

#Generate the bulk ExpressionSet object
bulk_obj <-ExpressionSet(assayData = bulk_counts_mtx, phenoData=phenoData)
class(bulk_obj@assayData)
str(bulk_obj)
class(exprs(bulk_obj))

#check dimensions of the epression matrix, should return "genes, samples"
dim(exprs(bulk_obj))
dim(exprs(bulk_obj_test))

#---- Create the SingeCellExperiment Object
#https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html


#Load data from the Azimuth pansci dataset
#https://azimuth.hubmapconsortium.org/references/#Mouse%20-%20Pansci

pansci_Seurat <-readRDS("pansci_filtered.rds")
str(pansci_Seurat)
head(pansci_Seurat@assays$RNA@counts, 5)

#Subset to keep only wild type data
pansci_Seurat_WT <- subset(pansci_Seurat, subset=Genotype=="WT")
#saveRDS(object = pansci_Seurat_WT, file="pansci_Seurat_WT.rds")
#pansci_Seurat_WT <- readRDS("pansci_Seurat_WT.rds")


pansci_Seurat@meta.data %>%
  group_by(Genotype) %>%
  summarise(no_rows = length(Genotype))
pansci_Seurat_WT@meta.data %>%
  group_by(Genotype) %>%
  summarise(no_rows = length(Genotype))


#Subset to keep only cells relative to the Heart
pansci_Seurat_WT_Heart<-subset(x=pansci_Seurat_WT, subset=Organ_name=="Heart")

#Remove Seurat objects not needed to free up space
rm(pansci_Seurat)
rm(pansci_Seurat_WT)

#Check that only desired cells have been included and which cell types
pansci_Seurat_WT_Heart@meta.data %>%
  group_by(Main_cell_type) %>%
  summarise(no_rows = length(Organ_Name))

PanSci_Heart_Cell_Types <- pansci_Seurat_WT_Heart@meta.data %>%
  group_by(Main_cell_type) %>%
  summarise(no_rows = length(Main_cell_type))

write.xlsx(PanSci_Heart_Cell_Types, "PanSci_Heart_Cell_Types.xlsx")                        

#Create singleCellExperiment object
library(SingleCellExperiment)
#sce_obj <- SingleCellExperiment(pansci_Seurat@assays$RNA@counts)
#even better we can just convert the Seurat object to a sce object from the Seurat package
library(Seurat)
sce_obj <- as.SingleCellExperiment(pansci_Seurat_WT_Heart, assay = "RNA")
sce_obj
head(sce_obj@colData)

#List the colnames of the colData in the sce_obj 
#--> will need to select columns to assign to the musci_prop(clusters = , samples =) arguments
colnames(colData(sce_obj))


#---- Perform the deconvolution with MuSic


#1.Check interect of genes between the two datasets
length(intersect(rownames(exprs(bulk_obj_test)), rownames(sce_obj)))
rownames(exprs(bulk_obj_test))
rownames(sce_obj)
#single cell data ENSEMBLE IDs have specified ENSEMBL version with .# suffix
#strip the ENSEMBLE version from the sce_obj data
rownames(sce_obj) <- sub("\\..*", "", rownames(sce_obj))
length(intersect(rownames(exprs(bulk_obj_test)), rownames(sce_obj)))


#Estimate cell type proportions
#The buylk.mtx=bulk_obj returned an issue with rowMeans(bulk.mtx) not not having two dimensions
#https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet I then used the exprs() function to only access the expression matrix

#Test a minimal dataset and subset of cell types to see if it works
#Est.prop.SchulteBulk_test = music_prop(bulk.mtx=exprs(bulk_obj_test), #the ExpressionSet object containing the bulk RNA-seq data
sc.sce=sce_obj, #the SingleCellExperiment object containing the singleCell data for the reference dataset
clusters="Main_cell_type" , #character, one of the column names of the colData(sce_obj)
samples="Cell_name", #character, one of the column names of the colData(sce_obj) data used as samples
verbose=TRUE,
select.ct = c("Atrial cardiomyocytes","Ventricular cardiomyocytes"))


#Run on the full dataset
Est.prop.SchulteBulk = music_prop(bulk.mtx=exprs(bulk_obj), #the ExpressionSet object containing the bulk RNA-seq data
                                  sc.sce=sce_obj, #the SingleCellExperiment object containing the singleCell data for the reference dataset
                                  clusters="Main_cell_type" , #character, one of the column names of the colData(sce_obj)
                                  samples="Cell_name", #character, one of the column names of the colData(sce_obj) data used as samples
                                  verbose=TRUE)



head(Est.prop.SchulteBulk$Est.prop.weighted, 5)

#Write table
write.xlsx(Est.prop.SchulteBulk$Est.prop.weighted, "CellProportions_MuSic.xlsx")





