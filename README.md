# Marker Gene Extraction Tool

A Python script that takes a phylogenetic structure as csv file and a large number of gene alignments.
Herafter it identifies a subselection of genes (=<10) capable of recovering the phylogenetic structure used as input.

## About this Project
I started this project as part of my Bachelor Thesis in order to find new marker genes for a novel phylogenetic structure based off of 2000+ conserved genes. 
It focuses on ease of use and applicability for other unresolved or revised taxonomies. 

### Usage:
  #Bash
*  conda create env -f environment.yml
*  conda activate markertool
*  python markerfinder.py config.yaml

### Input
*  Directory containing Multiple Sequence Alignment (MSA) files
*  .csv file indicating the phylogenetic structure (clade - species - isolate_ID; OR: species - isolate_ID)

### Output
Directory containing:
*  .csv file containing a diagnostic positions per gene per taxon matrix
*  .csv file containing a gene information matrix
*  .csv file containing a per taxon descriptive statistics matrix
*  .csv file containing an advanced gene information matrix
*  .txt file containing the selection of marker genes
*  .png file containing a heatmap of the marker genes diagnostic sites per taxon

### Important Functions:
*   **Per-site alignment analysis** 
*   **Diagnostic site identification** 
*   **Automated scoring framework**
*   **Automated marker gene selection**

## Technologies used:
*   **Python** 
*   **Pandas, Numpy, Counter, Seaborn, SeqIO** 

