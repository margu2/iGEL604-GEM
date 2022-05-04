# iGEL601
This repository contains Matlab scripts and supporting data for the manuscript "In silicoanalysis of the fast-growing thermophile *Geobacillus* sp. LC300 using a novel genome-scale metabolic model" by Emil Ljungqvist and Martin Gustavsson. 

The included scripts are dependent on the RAVEN (https://github.com/SysBioChalmers/RAVEN) and COBRA (https://opencobra.github.io/cobratoolbox/stable/) toolboxes for model reconstruction and analysis. To recreate the procedure described in the manuscript, start from the script file "master.m" to generate the model iGEL601.

# File structure

Scripts: Folder containing all scripts used in the project.

Data: Folder containing input data needed to generate iGEL601

Models: Folder containing the final model, iGEL601, in different formats for compatibility.

SamplingData: Folder containing flux distributions for the different scenarios analysed by flux sampling

## {{repository name}}: {{repository description}}

[![Version](https://badge.fury.io/gh/{{organization or username}}%2F{{repository name}}.svg)](https://badge.fury.io/gh/sysbiochalmers/yeast-gem)  
[![Zenodo](https://zenodo.org/badge/{{Zenodo ID}}.svg)](https://zenodo.org/badge/latestdoi/{{Zenodo ID}})  
[![Gitter chat](https://badges.gitter.im/{{organization or username}}/{{repository name}}.svg)](https://gitter.im/{{organization or username}}/{{repository name}})


#### Description

{{ fill in a short description or the paper abstract }}


#### Citation


#### Keywords

> Keywords are be separated by semicolons.
> The `Model source` field contains the source(s) of the current model, eg existing GEMs. If possible, use the Markdown format to add the URL with the DOI. The (NCBI) taxonomy ID should be provided in the [format from identifiers.org](https://registry.identifiers.org/registry/taxonomy). For the genome identifier, please provide the ENA/GenBank/RefSeq identifier via *identifiers.org*, or from other sources such as PATRIC or KBase.  

**Utilisation:** {{ experimental data reconstruction; multi-omics integrative analysis;, _in silico_ strain design; model template }}  
**Field:** {{ metabolic-network reconstruction }}  
**Type of model:** {{ reconstruction; curated }}  
**Model source:** {{ [YeastMetabolicNetwork](http://doi.org/10.1038/nbt1492) }}  
**Omic source:** {{ genomics; metabolomics }}  
**Taxonomic name:** {{ _Saccharomyces cerevisiae_ }}  
**Taxonomy ID:** {{ [taxonomy:559292](https://identifiers.org/taxonomy:559292) }}  
**Genome ID:** {{ [insdc.gca:GCA_000146045.2](https://identifiers.org/insdc.gca:GCA_000146045.2)  }}  
**Metabolic system:** {{ general metabolism }}  
**Tissue:**  
**Bioreactor:**    
**Cell type:**  
**Cell line:**  
**Strain:** {{ S288C }}  
**Condition:** {{ aerobic; glucose-limited; defined media }}  


### Installation

{{ Be mindful of users who do not have a typical background - provide a clear overview of the required software. Also, there might be different requirements for users and collaborators. }}


### Usage

The model can be imported into either MATLAB or Cobrapy, using either of the model formats found in the model folder


### Contributors
Scripts and modelling have been written and performed by Emil Ljungqvist and Martin Gustavsson
