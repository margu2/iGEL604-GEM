# iGEL601
This repository contains Matlab scripts and supporting data for the manuscript "_In silico_ analysis of the fast-growing thermophile _Geobacillus_ sp. LC300 using a novel genome-scale metabolic model" by Emil Ljungqvist and Martin Gustavsson. 

The included scripts are dependent on the [RAVEN](https://github.com/SysBioChalmers/RAVEN) and [COBRA](https://opencobra.github.io/cobratoolbox/stable/) Matlab toolboxes for model reconstruction and analysis. To recreate the procedure described in the manuscript, start from the script file "master.m" to generate the model iGEL604.

# File structure

Code: Folder containing all scripts used in the project.

Data: Folder containing input data needed to generate iGEL601, and data obtained from flux sampling

Models: Folder containing the final model, iGEL601, in different formats for compatibility.

Memote: Folder containing the Memote report of iGEL601


#### Description

{{ fill in a short description or the paper abstract }}


#### Citation

### Keywords

**Utilisation:** Experimental data reconstruction; _in silico_ strain design; model template 
**Field:** Metabolic-network reconstruction  
**Type of model:** Genome-scale  
**Taxonomic name:** _Geobacillus_ sp. LC300  
**Taxonomy ID:** [taxonomy:1519377](https://identifiers.org/taxonomy:1519377)  
**Genome ID:** [GCA_001191625.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001191625.1)
**Metabolic system:** General metabolism
**Bioreactor:** Belach Bioteknik GRETA system
**Condition:** Aerobic; Defined media; 

### Usage

The model can be imported into either MATLAB or Cobrapy, using either of the model formats found in the model folder


### Contributors
Scripts and modelling have been written and performed by Emil Ljungqvist and Martin Gustavsson
