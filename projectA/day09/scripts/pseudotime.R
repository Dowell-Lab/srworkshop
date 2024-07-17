#########INSTALLS############

install.packages("R.utils")
library(R.utils)
devtools::install_github('satijalab/seurat-wrappers')

########Library loads############
library(monocle3)
library(Seurat)
library(SeuratWrappers)

#load("/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/t21_official_umap_clust_06.25.24.Rdata")
load("/tmp/t21_official_umap_clust_06.25.24.Rdata")

#"fig1b_fbm_scaled_gex_updated_dr_20210104"       "fig5a_downs_fbm_scaled_gex_updated_dr_20210119"


ls()

t21@meta.data
t21@meta.data$broad_extfig7A_cell.labels
unique(t21@meta.data$broad_extfig7A_cell.labels)
unique(t21@meta.data$orig.ident)


table(t21@meta.data$cell.labels)

# pull only the red blood cells
Idents(object = t21) <- "broad_extfig7A_cell.labels"
cell_values <- c("erythroid")
red_T21 <- subset(t21, idents = cell_values, invert = FALSE)

table(red_T21@meta.data$cell.labels)

#change into a monocle cell type
cds <- SeuratWrappers::as.cell_data_set(red_T21)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(red_T21[["RNA"]])

#what is in my data?
monocle3::pData(cds)

#cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- monocle3::preprocess_cds(cds, method="PCA", num_dim = 3)

colnames(monocle3::pData(cds))

monocle3::plot_pc_variance_explained(cds)


cds <- monocle3::reduce_dimension(cds, reduction_method="UMAP",  preprocess_method = 'PCA')
#what are my options
#https://rdrr.io/github/cole-trapnell-lab/monocle3/
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds)


#what does the clustered data look like?
monocle3::plot_cells(cds)

#what can I color by?
colnames(monocle3::pData(cds))

#what is the data clustering by?
monocle3::plot_cells(cds, color_cells_by="broad_extfig7A_cell.labels")
monocle3::plot_cells(cds, color_cells_by="cell.labels")
monocle3::plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
monocle3::plot_cells(cds, color_cells_by="age", label_cell_groups=FALSE)


#The data looks like it clusters by age---- but there are only 4 people in this study. 
#So how many people are what age?
agedf <- as.data.frame(monocle3::pData(cds)['age'])
persondf <- as.data.frame(monocle3::pData(cds)['orig.ident'])
agepersondf <- as.data.frame(monocle3::pData(cds)[c("age",'orig.ident')])

table(agedf$age)
table(persondf$orig.ident)
table(agepersondf)

agepersonlanedf <- as.data.frame(monocle3::pData(cds)[c("lanes", "age",'orig.ident')])

table(agepersonlanedf)



### color by my favorite gene

rownames(cds)[grepl(pattern = "RUNX", x = rownames(cds))]

my_favorite_genes <- c("RUNX3", "RUNX2", "RUNX1")



monocle3::plot_cells(cds,
           genes=my_favorite_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

### optional code

cds_subset <- monocle3::choose_cells(cds)

cds <- monocle3::order_cells(cds)

monocle3::plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


isthisreallytime <- monocle3::pseudotime(cds, reduction_method = "UMAP")
