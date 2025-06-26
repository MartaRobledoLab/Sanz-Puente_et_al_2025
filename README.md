## Github Repository for
# Seed-mediated vertical transmission of Pantoea core endophytes
#### by Irene Sanz-Puente, Santiago Redondo-Salvo, Gloria Torres-Cortés, María de Toro, Susana Fernandes, Andreas Börner, Óscar Lorenzo, Fernando de la Cruz, and Marta Robledo.

Repository associated with the analysis used in Sanz-Puente et al., 2025, focus on the metabarcoding analysis of wheat bacterial endophytic community (doi: https://doi.org/10.1101/2025.01.06.628327). 
<i>This work is under review.</i> 

## Repository structure
- `data/`: Files containing the metadata, ASV sequences table and taxonomy table data are provided here 
- `scripts/`: bioinformatics and statistical analysis  
  - `qiime.sh` – read quality control and trimming  
  - `Metabarcoding_analysis.R` – taxonomic profiling, statistical tests and visualizations   
  - `qPCR.R` – data processing, statistical tests and visualizations  
- `results/`: figures and tables generated from the analysis

## Data analysis
### Running the pipeline
Raw read sequences of this study have been deposited in the NCBI Sequence Read Archive (RSA) under the BioProject accession number PRJNA1282304. Other data used in the analysis can be found in the Supplementary material of the paper or as data frame in scripts.
```bash
bash scripts/preprocess.sh data/raw/*.fastq -o data/processed/
````
### Requirements
- R (>= 4.2) with the following packages:  
```r
biomformat      # Handling BIOM format files
ComplexUpset    # Advanced UpSet plots for set intersections
FSA             # Fisheries Stock Analysis, used for Dunn’s test
ggplot2         # Core plotting library
ggsci           # Scientific journal themes for ggplot2
ggsignif        # Significance annotations for ggplot2
officer         # Exporting Word documents from R
pairwiseAdonis  # Pairwise comparisons for PERMANOVA
phyloseq        # Microbiome data analysis
pheatmap        # Pretty heatmaps
reshape2        # Data reshaping
rvg             # Vector graphics export (e.g., for Word/PowerPoint)
tidyverse       # Data wrangling and visualization (includes dplyr, ggplot2, etc.)
vegan           # Community ecology and multivariate statistics
````

#### Citation
If you use this repository or its contents, please cite:
Sanz-Puente et al. (2025). Seed-mediated vertical transmission of Pantoea core endophytes. bioRxiv. https://doi.org/10.1101/2025.01.06.628327v1
