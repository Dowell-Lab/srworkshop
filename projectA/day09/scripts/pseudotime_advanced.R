#########INSTALLS############

#install.packages("R.utils")
library(R.utils)
#devtools::install_github('satijalab/seurat-wrappers')


########Library loads############
library(monocle3)
library(Seurat)
library(SeuratWrappers)


#load("/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/t21_official_umap_clust_06.25.24.Rdata")
load("/tmp/t21_official_umap_clust_06.25.24.Rdata")

ls()

#Lets take a look at the data before we use it
t21@meta.data
t21@meta.data$broad_extfig7A_cell.labels
unique(t21@meta.data$broad_extfig7A_cell.labels)
unique(t21@meta.data$orig.ident)


#how many uniuqe cell.lables are there?
table(t21@meta.data$cell.labels)

#how many uniuqe broad cell types are there?

# How does the cell cycle look in these cells?
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#cds <- CellCycleScoring(t21, s.features = s.genes, g2m.features = g2m.genes)

# We have too many cells for this to go quick. Lets pull only the red blood cells.
Idents(object = t21) <- "broad_extfig7A_cell.labels"
cell_values <- c("erythroid")
red_T21 <- subset(t21, idents = cell_values, invert = FALSE)

table(red_T21@meta.data$cell.labels)

#change into a monocle cell type
cds <- SeuratWrappers::as.cell_data_set(red_T21)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(red_T21[["RNA"]])

#what is in my data?
monocle3::pData(cds)
colnames(monocle3::pData(cds))

#what are my options
#https://rdrr.io/github/cole-trapnell-lab/monocle3/



#the next few lines are all about setting up the monocole obeject so we can run pseudotime and are required:
#required preprocess_cds, reduce_dimension, cluster_cells, learn_graph

#cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- monocle3::preprocess_cds(cds, method="PCA", num_dim = 3)

#how much varience explained by pca
monocle3::plot_pc_variance_explained(cds)


cds <- monocle3::reduce_dimension(cds, reduction_method="UMAP",  preprocess_method = 'PCA')

cds <- monocle3::cluster_cells(cds)

#Do you want the program to consider all the cells as part of the same cell type progression for psudotime? 
#If yes, then use_partition = FALSE
#otherwise use the defaults for learn_graphs

#cds <- monocle3::learn_graph(cds)
cds <- monocle3::learn_graph(cds, use_partition = FALSE)

#this part is not required to run pseudotime 
#but if you don't know what the data looks like you won't know where to set your start and end points


#what does the clustered data look like?
monocle3::plot_cells(cds)

#what can I color by?
colnames(monocle3::pData(cds))

#what is the data clustering by?
monocle3::plot_cells(cds, color_cells_by="broad_extfig7A_cell.labels")
monocle3::plot_cells(cds, color_cells_by="cell.labels")
monocle3::plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
monocle3::plot_cells(cds, color_cells_by="age", label_cell_groups=FALSE)

#I want pseudotime acorss early mid and late red blood cells. 
#Which of the graphs above shows me where to start and end the pseudo time?




### color by my favorite gene
#How do I find out what my favorite gene is named in this object?
rownames(cds)[grepl(pattern = "RUNX", x = rownames(cds))]

#This is a list of genes I want to see where they are in this plot
my_favorite_genes <- c("RUNX3", "RUNX2", "RUNX1", "HBB", "PCNA")



monocle3::plot_cells(cds,
           genes=my_favorite_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


#if you only want to pick start sites use order_cells

cds <- monocle3::order_cells(cds, reduction_method = "UMAP")

monocle3::plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
deg <- subset(cds_pr_test_res, q_value < 0.0005) 
deg <- subset(deg, morans_I > 0.5) 
deg <- deg[order(deg$q_value),]
deg
pr_deg_ids <- row.names(deg)


genes_look_interesting <-c("HMGN2", "RPL22", 'ATPIF1', 'HBB', "BLVRB", "ATP5E")

monocle3::plot_cells(cds,
                     genes=genes_look_interesting,
                     label_cell_groups=FALSE,
                     show_trajectory_graph=FALSE)

plot_genes_by_group(cds,                    # our CDS object
                    genes_look_interesting,                 # a list of gene names to show in the plot
                    group_cells_by="cell.labels",         # how to group cells when labeling
                    ordering_type="maximal_on_diag")  # how to order the genes / groups on the dot plot



#MOST of these are ribosomes, maybe I don't care about ribosomes
pr_deg_ids_ribosomal_large <- pr_deg_ids[grepl(x=pr_deg_ids, pattern="RPL")]
pr_deg_ids_ribosomal_small <- pr_deg_ids[grepl(x=pr_deg_ids, pattern="RPS")]
pr_deg_ids_other <- pr_deg_ids[is.na(match(pr_deg_ids, pr_deg_ids_ribosomal_small))]
pr_deg_ids_other <- pr_deg_ids_other[is.na(match(pr_deg_ids_other, pr_deg_ids_ribosomal_large))]



plot_genes_by_group(cds,                    # our CDS object
                    pr_deg_ids_other,                 # a list of gene names to show in the plot
                    group_cells_by="cell.labels",         # how to group cells when labeling
                    ordering_type="maximal_on_diag")  # how to order the genes / groups on the dot plot


sort(pr_deg_ids_other)
SLC_genes_that_change <- pr_deg_ids_other[grepl(pattern = "SLC", x = pr_deg_ids_other)]


plot_genes_by_group(cds,                    # our CDS object
                    SLC_genes_that_change,                 # a list of gene names to show in the plot
                    group_cells_by="cell.labels",         # how to group cells when labeling
                    ordering_type="maximal_on_diag")  # how to order the genes / groups on the dot plot



subsetgenes_lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes_look_interesting,
                       colData(cds)$cell.labels %in% c("early erythroid", "mid erythroid")]


monocle3::plot_genes_in_pseudotime(subsetgenes_lineage_cds, vertical_jitter = TRUE, cell_size = 0.1)

df<- as.data.frame(subsetgenes_lineage_cds$data)


#if you  want to pick start sites and end sites use choose_graph_segments

cds_sub <- choose_graph_segments(cds)

monocle3::plot_cells(cds_sub,
                     color_cells_by = "pseudotime",
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE,
                     graph_label_size=1.5)




### optional code

cds_subset <- monocle3::choose_cells(cds)

cds_subset <- monocle3::reduce_dimension(cds_subset, reduction_method="UMAP",  preprocess_method = 'PCA')

cds_subset <- monocle3::cluster_cells(cds_subset)

#what would you do after clustering cells?


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

#How does cell cycle look in these cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cds <- CellCycleScoring(cds, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



#This is another way too look at genes acorss the clusters

#~ in front of a column name is what model formula wants
gene_fits <- fit_models(cds, model_formula_str = "~cell.labels")
fit_coefs <- coefficient_table(gene_fits)

