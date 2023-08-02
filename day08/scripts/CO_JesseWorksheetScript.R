###Seurat Worksheet###

install.packages("devtools")
library(devtools)

devtools::install_github("satijalab/seurat") ###if this throws an error, use the line of code below instead (after the hash)
#install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/", repos = NULL, type="source") 
library(Seurat)

workdir <- '/Users/christopherozeroff/Desktop/workshop/cellranger_outputs_R_objects/Seurat_4.3.0.tar.gz' ###Change the directory path to your own

obj.4dpi.10x <- Read10X(paste(workdir, "dpi4_filtered_feature_bc_matrix/", sep = ""))
obj.4dpi.seurat <- CreateSeuratObject(counts = obj.4dpi.10x,
                                      projects = "4dpi",
                                      min.cells = 3,
                                      min.features = 200)

obj.7dpi.10x <- Read10X(paste(workdir, "dpi7_filtered_feature_bc_matrix/", sep = ""))
obj.7dpi.seurat <- CreateSeuratObject(counts = obj.7dpi.10x,
                                      projects = "7dpi",
                                      min.cells = 3,
                                      min.features = 200)

obj.7dpi.seurat$cellid <- 'dpi7' 
obj.4dpi.seurat$cellid <- 'dpi4' 


##Quality control##

obj.7dpi.seurat[["percent.mt"]] <- PercentageFeatureSet(obj.7dpi.seurat, pattern = "mt-")

#visualize precleaning
VlnPlot(obj.7dpi.seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

obj.7dpi.seurat <- subset(obj.7dpi.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#visualize post cleaning
VlnPlot(obj.7dpi.seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

##repeat for 4dpi##

obj.4dpi.seurat[["percent.mt"]] <- PercentageFeatureSet(obj.4dpi.seurat, pattern = "mt-")

#visualize precleaning
VlnPlot(obj.4dpi.seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)

obj.4dpi.seurat <- subset(obj.4dpi.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#visualize post cleaning
VlnPlot(obj.4dpi.seurat, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)


##Merging/integration##
merged <- merge(obj.4dpi.seurat, obj.7dpi.seurat)

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(object = merged))
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

list <- SplitObject(merged, split.by = "cellid")

anchors.rpca <- FindIntegrationAnchors(list, dims = 1:20, anchor.features = 2000, reduction = 'rpca')


##Clustering and visualization##
load(paste0(workdir, 'final_myogenic_subset_seurat.RData'))

#set parameters
dims = 1:15
min.dist = 0.05
n.neighbors = 200
npcs = 30
resolution = 0.2
algorithm = 1

#Scaling
DefaultAssay(myogenic.subset.umap.final) <- "integrated"
myogenic.subset.umap.final <- ScaleData(myogenic.subset.umap.final, verbose = FALSE)

#PCA

myogenic.subset.umap.final <- RunPCA(myogenic.subset.umap.final, npcs = npcs, verbose = FALSE)

#UMAP
myogenic.subset.umap.final <- RunUMAP(myogenic.subset.umap.final, reduction = "pca",
                                      n.neighbors = n.neighbors,
                                      min.dist = min.dist,
                                      dims = dims)

#identify clusters
myogenic.subset.umap.final <- FindNeighbors(myogenic.subset.umap.final)
myogenic.subset.umap.final <<- FindClusters(myogenic.subset.umap.final, resolution = resolution,
                                           algorithm = algorithm)

DefaultAssay(myogenic.subset.umap.final) <- "RNA"
DimPlot(myogenic.subset.umap.final, pt.size = 1, label = T)
DimPlot(myogenic.subset.umap.final, pt.size = 1, label = T, group.by = 'cell.id')
DimPlot(myogenic.subset.umap.final, pt.size = 1, label = T, group.by = 'injury.dpi')

