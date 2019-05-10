# PhosMap
---
## A Comprehensive R Package For Analyzing Quantitative Phosphoproteomics Data
---

### Functions summary

*1. Extracting the confidence probability of phosphorylation sites at peptide level from identification results searched by Mascot.*<br> 
*2. Generating the quality control file of phosphorylation sites based on score of sites from Mascot.*<br> 
*3. Pre-processing phosphoproteomic data.*<br> 
*4. Kinase substrate enrichment analysis.*<br> 
*5. Motif enrichment analysis.*<br> 
*6. Data visualization*<br> 


### Imports
`graphics`, `grDevices`, `stats`, `utils`, `stringr`, `ggseqlogo`, `samr`, `limma`, `e1071`, `ClueR`, `Rtsne`, `glmnet`, `yaml`, `impute`
<br> 


### Installation
```R
install.packages('devtools')
require(devtools)
install_github('ecnuzdd/PhosMap')
```

### Source location of the built-in reference library
[PhosMap_datasets](https://github.com/ecnuzdd/PhosMap_datasets)  [https://github.com/ecnuzdd/PhosMap_datasets](https://github.com/ecnuzdd/PhosMap_datasets)  <br> 
* [data](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/data) (Demo data: BRAFi.RData)
* [fasta_library](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/fasta_libarary)
  * Refseq (Human, Mouse, Rattus)
  * Uniprot (Human, Mouse, Rattus)
* [id_coversion_table](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/id_coversion_table)
  * Human, Mouse, Rattus
* [kinase_substrate_regulation_relationship_table](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/kinase_substrate_regulation_relationship_table)
  * Human, Mouse, Rattus
* [motif_library](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/motif_library)
  * Refseq (Human, Mouse, Rattus)
  * Uniprot (Human, Mouse, Rattus)
