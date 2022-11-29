# PhosMap
---
## A Comprehensive R Package For Analyzing Quantitative Phosphoproteomics Data
---

### Functions summary

*1. Extracting the confidence probability of phosphorylation sites at peptide level from identification results searched by Mascot.*<br> 
*2. Generating the quality control file of phosphorylation sites based on score of sites from Mascot.*<br> 
*3. Pre-processing phosphoproteomic data.*<br> 
*4. Kinase activity prediction.*<br> 
*5. Motif enrichment analysis.*<br> 
*6. Data visualization*<br> 


### Imports
`graphics`, `grDevices`, `stats`, `utils`, `stringr`, `ggseqlogo`, `samr`, `limma`, `e1071`, `ClueR`, `Rtsne`, `glmnet`, `yaml`, `impute`
<br> 

### External dependencies
* `ksea` (https://github.com/evocellnet/ksea)
```R
install.packages('devtools')
require(devtools)
install_github('evocellnet/ksea')
```
* `rmotifx` (https://github.com/omarwagih/rmotifx)
```R
install.packages('devtools')
require(devtools)
install_github('omarwagih/rmotifx')
```

### Installation of PhosMap
```R
install.packages('devtools')
require(devtools)
install_github('ecnuzdd/PhosMap')
```

### User guide of PhosMap
* [Tutorial of PhosMap](https://github.com/ecnuzdd/PhosMap_datasets/blob/master/Tutorial_of_PhosMap.pdf) (Tutorial_of_PhosMap.pdf)
* [Functions manual of PhosMap](https://github.com/ecnuzdd/PhosMap_datasets/blob/master/Manual_of_PhosMap.pdf) (Manual_of_PhosMap.pdf)
* [The help for installation](https://github.com/ecnuzdd/PhosMap_datasets/blob/master/lnstallation_help.pdf) (lnstallation_help.pdf)

### Case study
You can go to this directory to find the files needed for the case study.

``` ftp://111.198.139.72:4000/PhosMap/data/ ```

Or download them directly through the following links.

#### Complete workflow
This .zip contains PhosMap's demo code and raw datasets. After unzip, start the R program with setting the new folder as your working directory, and run the R script in it.

``` ftp://111.198.139.72:4000/PhosMap/data/PhosMap_Demo_With_BRAFi_Data.zip ```

#### Partial data
The .RData file stores the input objects required to invoke functions in the PhosMap script demo, which can be downloaded from the following link.

``` ftp://111.198.139.72:4000/PhosMap/data/BRAFi.RData ```

After downloading, put it into your working directory. Use the following command to introduce these objects into the global environment,  then you can smoothly execute the statements in the script.

```R
load(BRAFi.RData)
```

### Source location of the built-in reference library
[PhosMap_datasets](https://github.com/ecnuzdd/PhosMap_datasets)  [https://github.com/ecnuzdd/PhosMap_datasets](https://github.com/ecnuzdd/PhosMap_datasets)  <br> 
* [data](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/data) (Demo data: BRAFi.RData)
* [fasta_library](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/fasta_library)
  * Refseq (Human, Mouse, Rattus)
  * Uniprot (Human, Mouse, Rattus)
* [id_coversion_table](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/id_coversion_table)
  * Human, Mouse, Rattus
* [kinase_substrate_regulation_relationship_table](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/kinase_substrate_regulation_relationship_table)
  * Human, Mouse, Rattus
* [motif_library](https://github.com/ecnuzdd/PhosMap_datasets/tree/master/motif_library)
  * Refseq (Human, Mouse, Rattus)
  * Uniprot (Human, Mouse, Rattus)

### Complete demo data processed by [Firmiana](http://www.firmiana.org)
#### Raw data from deposited in the ProteomeXchange Consortium
[PXD007740](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD007740)

#### Processed data by Firmiana
```
ftp://111.198.139.72:4000/PhosMap/PhosMap_Demo_With_BRAFi_Data.zip
```

#### *Reference*
*1. Ressa, A., et al. (2018) A System-wide Approach to Monitor Responses to Syner-gistic BRAF and EGFR Inhibition in Colorectal Cancer Cells, Molecular & cellular proteomics : MCP, 17, 1892-1908.* <br>
*2. Feng, J., et al. (2017) Firmiana: towards a one-stop proteomic cloud platform for data processing and analysis, Nature biotechnology, 35, 409-412.*











