# Annotating Cell Types
Author: Hope Townsend June 24, 2024

## How do I use this worksheet?
* This worksheet corresponds to the script in ../scripts/Cell_Type_Annotations.Rmd. The information (not including photos) is both here and in that script so feel free to just open the script in RStudio and not look at this worksheet if that's easier.

## Overview
We will go over two ways we can annotate cell types automatically:
* **Known atlas**: In this case, we will use a cell type atlas with cell types that we expect to be in our dataset to annotate ours. We do this by mapping each of our cells to the cell types identified in the annotated set. 
The Seurat based Tutorial can be found [here](https://satijalab.org/seurat/articles/integration_mapping.html)
![alt text](<Screenshot 2024-06-25 at 9.01.55 PM.png>)
* **(Un)Known marker genes**: SC-Type uses marker genes that are either expected or not expected in certain cell types. You can provide your own set, but it has a database for the following tissues: Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen Stomach, Thymus.
![alt text](<Screenshot 2024-06-25 at 9.08.28 PM.png>)

## A) Known Atlas
Some great resources to find already defined atlases (good for references) include https://data.humancellatlas.org/, https://www.synapse.org/, and most importantly PubMed/google scholar. Most single-cell papers have a section with the finished data in the "Data Availability" or "Code Availability" section. These atlases will include the celltypes in the metadata. 
* *Note:* Sometimes, they are provided as AnnData objects rather than Seurat objects (.h5ad). These can be converted to Seurat objects via the function h5ad2seurat() from the package [schard](https://github.com/cellgeni/schard).


### Example
**0. Get the data**

For demonstration purposes, we are using an object which contains half of the T21 dataset as the reference atlas to label the cell types for the other half of the T21 cells. I first will split up randomly where each cell type is equally split between the "reference" and "query"

**1. Graph the original UMAP to see that the datasets are in fact split**
    
Let's first check that there is an equal split between the cells by looking at the UMAP (stored under reduction Xumap_)

This UMAP will look different from what you did yesterday because we are using only cells belonging to individuals with T21 for the variable genes, clustering and UMAP.

**2. Follow Seurat instructions to transfer information**

We are treating the t21_ref as our "reference" since we expect the same cell types to exist. In reality, we would not have the annotations for our data and therefore would need to use a different dataset. Examples of where to get these are found in the previous descriptions. The full Seurat based Tutorial followed below can be found [here](https://satijalab.org/seurat/articles/integration_mapping.html)

   A. Need to find anchors (FindTransferAnchors): pairs of cells from the reference and query datasets that are mutual nearest neighbors. This means that when you look at the dimension reductions of both, the cells are very similar to one another. Details are found [here](https://satijalab.org/seurat/reference/findtransferanchors).
    
   B. Transfer the reference data labels of interest to the query data (TransferData): This uses the anchors to help compare individual cells between one another annd transfer labels accordingly. Details are found [here].(https://satijalab.org/seurat/reference/transferdata)
    
   C. Add the determined labels to the metadata of the query (AddMetaData).
    
**3. See how well this did**
Since we have the truth set, we can get an idea of how many of the cell types we predicted were accurate for the individual cells.
    
### Do on your own
**4. Repeat the analysis on your own but with the smaller scale cell labels.**
* Hint: the smaller scale cell types are under the metadata column "cell.labels" 

#### Reflection Questions
    
1. How well did we do with the large scale cell types compared to the small scale cell states?
2. Therefore, how much does clustering resolution and pattern play a role in cell type annotation?
3. Which parts of the code would you change to look at the seurat clusters (in this case the query clusters used would be the $seurat_clusters instead of broad_extfig7A_cell.labels -- check that you know how to do this)?
    
    
## B) (Un)Known Marker Genes
SC-Type uses marker genes that are either expected or not expected in certain cell types. The full instructions on how to use it can be found at their [Github README](https://github.com/IanevskiAleksandr/sc-type). They also have an interactive [web implementation](https://sctype.app/) for it, but this can be more difficult to debug or use if you have really large files. Therefore, we will be showing the use in R.

### Example
**1. Load Functions and Database for SC-Type**
* SC-Type allows you to load the functions directly which is what we're doing with "source()"
* We then need a database that has the list of genes that are considered "marker" and "non-marker" genes for cell types of interest. In this case we use the default one they already made for the "Immune system", however you can also use your own by uploading a XLSX file with the following headings (one line is provided as an example):

|tissueType|cellName|geneSymbolmore1|geneSymbolmore2|
|:---|:---|:---|:---|
|Immune system|Platelets|CD41,CD42b|SELL,CCR7|

**2. Prep the SC data for analysis**

We want to get the scaled gene data matrix. SC-Type provides suggested code for getting this according to which version of Seurat you're using.

**3. Run ScType**

sctype_score() gets the scores for each cell for beinging within all the different possible cell types. This will be used downstream to predict the clusters.
* scaled=TRUE is used since the data we're using is scaled
* gs = gs_positive -- the gene sets we expect to be in the cell types considering
* gs2 = gs_negative -- the gene sets we don't expect to be in teh cell types considering

**4. Based on single cell scores, assign the clusters with cell type annotations**
* Since we know the cell types, instead of using Seurat based clusters as the clusters to assign cell types too, we can use broad_extfig7A_cell.labels as the cluster annotations we need to assign. Usually, we'd use seurat_clusters to replace broad_extfig7A_cell.labels.
* We then use the es.max object that has the individual cell scores to predict which of cell type oof the clusters.
* Any clusters that have low confidence scores regarding the linked cell type are labeled as such.

**5. Look at results**
* Let's look at the annotated cell types (cluster = original cluster name (if using seurat clusters, these would be numbers), type = the predicted cell type, scores= the scores corresponding to the predicted cell type)

* **Interpretation**: Stroma was mislabeled as macrophages with a high score still, although macrophages and stromal cells have lots of cross-talk which might help explain the result. These are not the same cell type although the two are similarly related. Therefore, 174 cells were mislabeled compared to the 315 mislabeled with the Seurat integration mapping. Therefore, looking at multiple sources of information and marker genes is still important, regardless of the confidence of the tool sometimes.

* Now graph the UMAP to see how well the cell types compare visually. I made the colors of corresponding cell states in the graphs the same by using scale_color_manual.

### Do on your own 
**6. Your turn with the smaller scale clusters**
Use the code previously shown and the outline provided to repeform similar analysis but using the smaller scale clusters (cell.labels) as the clusters of interest.


#### Reflection Questions
    
1. Was the default database used by this tool sensitive enough to distinguish between the more finely tuned cell states? (*Hint*: do we see a lot of the same cell type for different clusters?)
2. If yes, how might you be able to still use this tool but have a more sensitive database? (*Hint*: look at step 1 of this section)
