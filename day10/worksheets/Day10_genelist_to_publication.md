# Suggestions on getting from Day 9 to a publication. 
<img width="1120" height="628" alt="image" src="https://github.com/user-attachments/assets/ba16c93d-6b7a-41c5-9b90-4a3124fda52a" />

## Sequencing experiments are just like any other experiment...
They are best with a clear hypothesis and good controls. You should be able to succinctly articulate what question you'll ask of your data and *have a plan to process the data to answer it*. Otherwise, your experiment probably needs more planning. 

Sometimes people sequence because they don't know what else to do. The problem? Sequencing generates *massive amounts of data* that now you need to interpret. With nowhere to start. That's a recipe for a lot of time and money spent for little payout. 

**TL;DR**  The best way to get from a gene list to a publication is to 
1. Have a clear question
2. Know your field well and why your question is important
3. Learn to make good figures

Read on to learn more about making good figures! 


### First, is a gene list good for anything? YES!!
A gene list is good for lots of things! Here are the top things to use it for: 
**1. Check your positive and negative controls.** Every experiment should have control genes you expect to change in a predictable way. Make sure your controls are acting as you expect. You don't want to analyze a dataset if your treatment didn't even work. 
**2. Are the top changing genes interesting to people in your field?** Sometimes the most differentially expressed genes are a direct hint as to what sort of changes are going on in your experiment. 


# Steps to go beyond a gene list

1. Know your field
2. Pathway analysis
3. Mine public databases
4. Visualize your nextworks and pathways (aka, put your differentially expressed genes into their biological context)
5. Combine with additional datasets, timepoints, etc
6. **Know your field** and find a story.

# 1. Know your field 
Read every day. 

Genes are just genes. Mostly just three letter acronyms. A bioinformatics pipeline can't tell you why they're important. Even ChatGPT isn't very good at that yet. As scientists, we have to keep struggling until we figure out why the genes matter in the context we're studying. 

Trust me, this is harder than learning how to code (sorry). 

# 2. Use pathway analysis 

For detailed instructions, please refer to [Day10_GSEA.md](Day10_GSEA.md) and [Day10_GO_analysis_walkthrough](Day10_GO_analysis_walkthrough). 

**What's the point of pathway analysis?** 
We want to put our differentially expressed genes into the context of cellular functions that are easier to understand. <ins> In my experience, people are overly optimistic that pathway analysis will tell them THE answer. </ins> It's an extremely helpful tool, but can be confusing. 

The most common challenges with pathway analysis are: 
* The pathways are too general, and don't lead to specific hypotheses about what's going on.
* The pathways are related to a different cell type.
* There are so many pathways! Where does anyone start figuring out which ones are important?

When this happens, look into which pathways you're testing. For **GSEA**, there are multitudinous different [pathway collections](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp) you can test. Test sets of pathways specifically instead of testing them all at once. **GO** relies on different databases, and similarly you can test specific databases. For example, [KEGG](https://www.genome.jp/kegg/), [Reactome](https://reactome.org/), [WikiPathways](https://www.wikipathways.org/). 

Ultimately, you're hoping to get a figure that compares multiple pathways and shows only the ones significantly enriched. Like this one below: 

<figure>
<img width="414" height="236" alt="image" src="https://github.com/user-attachments/assets/70f519ab-a54f-46d1-b553-7bed1a21f161" />
<figurecaption>https://bioinformatics.sdstate.edu/go/</figurecaption>
</figure>

### Similar tools to check out: 
* Transcription Factor Enrichment Analysis ([TFEA](https://github.com/Dowell-Lab/TFEA))
* [TFEA-ChIP](https://pubmed.ncbi.nlm.nih.gov/31347689/)
* ChIP Enrichment Analysis ([ChEA](https://pmc.ncbi.nlm.nih.gov/articles/PMC6602523/))
* Qiagen's [IPA](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/) (Very $$$ but actually pretty good)


# 3. Mine public databases

Public databases are incredible ways to search for what other proteins might be interacting with your differentially expressed genes, what TF's are known to regulate it, and if it's associated with any known disease or SNPs. 

**Interaction databases** \
Use these databases to figure out what other proteins your specific gene-of-interest could be interacting with. This is a *hypothesis generating* exercise, and can be very helpful in figuring out what you want your next experiment to be. 

[Reactome](https://reactome.org) \
[String](https://string-db.org) 

<figure>
<img width="264" height="265" alt="image" src="https://github.com/user-attachments/assets/4e929b9c-dd45-418b-bdd8-5438f1f6e57f" />
<figurecaption> STAT1 interaction network, generated from String</figurecaption>
</figure>



**Known mutatios and clinical relevance** \
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) \
Online Mendelian Inheritance in Man [OMIM](https://www.omim.org) 

Below, I've searched STAT1 in ClinVar and am looking at a table of known variants that are pathogenic. 

<img width="631" height="280" alt="image" src="https://github.com/user-attachments/assets/0bea8080-eb2c-4ea9-8fc7-73cf297fb408" />


# 4. Visualize networks and pathways 

Differential expression data can be visualized in the context of networks and pathways using R packages like [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) or GUI interfaces like [Cytoscape](https://cytoscape.org/). 

This is a [good tutorial](https://cytoscape.org/cytoscape-tutorials/protocols/basic-data-visualization/#/) for getting started with Cytoscape. My [<ins>favorite</ins> way to use Cytoscape is with Wikipathways](https://cytoscape.org/cytoscape-tutorials/protocols/wikipathways-app/#/title), so check out this tutorial. It allows you to make pathways like the one below, where differential expression is mapped onto each node (blue is decreased expression, red is increased expression). 
<figure>
<img width="674" height="515" alt="image" src="https://github.com/user-attachments/assets/9b87be23-08d4-4433-8c66-09c561ecf400" />
<figurecaption>https://cytoscape.org/cytoscape-tutorials/protocols/wikipathways-app/#/3/5</figurecaption>
</figure>


Cytoscape is a little clunky to learn, but extremely useful once you know how to use it. And you can export your figures as .svg's and take them into Adobe Illustrator or Inkscape to finalize your figures. 

# 5. Integrate your data with other type of omics 

<img width="408" height="230" alt="image" src="https://github.com/user-attachments/assets/0c5c2b7f-b0a1-499b-830c-d53126ddd867" />

#### Working with published data: 
When working with previously published data, you'll first need to process the data as appropriate. Then, you'll be able to bring the processed data into R or Python and work with it similar to how we did on Day 9 in the [R graphing worksheet](https://github.com/Dowell-Lab/srworkshop/blob/main/projectB/day09/worksheets/01_r_graphing_worksheet.md). 

#### Performing your own multiomics experiment 
If you're performing your own multiomics experiments, there are most advanced tools you can use to integrate data. However, usually these tools come with the caveat that all the omics where collected from the same population of cells at the same time! So be careful wtih your experimental design. A package I've liked using in the past is [mixOmics](https://mixomics.org/). 


# 6. Pick a story to tell! 

<img width="600" height="337" alt="image" src="https://github.com/user-attachments/assets/b396dbb2-08f9-46cb-bdd5-218954d2663b" />

Sequencing datasets are very rich, with massive amounts of data. If you're not sure what the story is, either there is one, or there are many different stories you could tell! A PhD or postdoc is only so long. Eventually you have to pick the story that is the most interesting to you, your PI, and the field and tell that story. You can always come back and use the dataset again in another paper. 
