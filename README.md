# Marker Gene Extraction Tool

A Python script that takes a phylogenetic structure as csv file and a large number of gene alignments.
Herafter it identifies a subselection of genes (=<10) capable of recovering the phylogenetic structure used as input.

## About this Project
I started this project as part of my Bachelor Thesis in order to find new marker genes for a novel phylogenetic structure based off of 2000+ conserved genes. 
It focuses on ease of use and applicability for other unresolved or revised taxonomies. 

### Usage
  #Bash
conda create env -f environment.yml
conda activate markertool
python markerfinder.py config.yml

### Important Functions:
*   **Per-site alignment analysis** 
*   **Diagnostic site identification** 
*   **Automated scoring framework**
*   **Automated marker gene selection**

## Technologies used:
*   **Python (Currently Jupyter Notebook)** 
*   **Pandas, Numpy, Counter, Seaborn, SeqIO** 

