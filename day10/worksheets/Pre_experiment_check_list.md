# Pre-Experiment Checklist (Data Science / Bioinformatics Version)

A structured, falsification-first workflow for validating a hypothesis before committing to full analysis. Inspired by Strong Inference (Platt, 1964) combined with bioinformatics-specific data QA and tool validation practices.

## 1. State the Hypothesis

Write a specific, falsifiable claim before touching any data.

- What is the exact hypothesis?
- What would disprove it?

## 2. What graph(s) would show the hypothesis was true?

Sketch it on paper first — before analysis — to avoid confirmation bias in later interpretation. Draw several versions if more than one plot type could show the effect (e.g., PCA plot, volcano plot, boxplot, heatmap).

- Graph 1:
- Graph 2 (alternate view):

## 3. How else could I get those graphs and my hypothesis still be false?

List every plausible confound or artifact that could produce the same visual pattern without the hypothesis being true.

- Batch effects
- Sample size imbalance
- Confounding covariates
- Technical/platform artifacts
- Other:

## 4. If the answer above was anything other than "none" — what graph would still show the hypothesis is true?

Design a graph that explicitly controls for the confound(s) identified above (e.g., PCA colored by batch, covariate-adjusted residual plot).

- Adjusted graph:

## 5. What data do I need to make that graph?

- Data source(s):
- Sample size / coverage needed:
- Metadata / covariates required:

## 6. What issues will I have to account for in that data?

- Normalization
- Scaling
- Units
- Missing values
- Known mistakes / mislabeling
- Batch/platform differences

## 7. How am I going to clean that data?

- Cleaning steps:
- QC thresholds:
- Outlier handling:

## 8. How am I going to combine the multiple datasets?

- Unique identifiers per "sample"
- Join/merge keys:
- Harmonization strategy (units, IDs, formats):
- Batch correction method (if needed):

## 9. What tools am I going to use?

- Tool(s):
- Language/environment:

## 10. Are they installed, or do I have to install them?

- [ ] Installed and verified
- [ ] Needs installation — installation plan:

## 11. Is the tool made for my type of data?

- Intended data type(s) per documentation:
- Match to my data type: Yes / No / Partial — explain:

## 12. What artifacts could the tool cause?

- What are the assumptions of the tool? Does my data violate those assumptions?
- What are the optional parameters for that tool? Should I change any of the default parameters?

## 13. Should I test the tool with simulated data?

- [ ] Yes — simulate data with known ground truth first
- Simulation plan:
- Expected result if tool works correctly:

## 14. Start the Experiment

Once all steps above are resolved, proceed to run the analysis on real data.

