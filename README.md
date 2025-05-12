## ETBF infection alters stromal cell gene expression in implant fibrosis
This repository contains processed data and code supporting the publication [Murine gut microbiota dysbiosis via enteric infection modulates the foreign body response to a distal biomaterial implant.](https://www.pnas.org/doi/10.1073/pnas.2422169122)

### Overview

This repository provides the bash and R scripts used to process and analyze the bulk RNA sequencing data used in the study, which investigates the role of ETBF infection on the transcriptomic profile of stromal cells surrounding a biomaterial (polycaprolactone, PCL) implant.

The raw data and aligned counts files have been deposited in NCBI's Gene Expression Omnibus (Edgar et al., 2002) and are accessible through GEO Series accession number [GSE293564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293564), containing RNA sequencing data from mice 6 weeks after a volumetric muscle loss injury and implantation of a synthetic biomaterial, with and without ETBF infection (n = 4, each group). On Day -11, both groups were treated with antibiotics. On Day -7, antibiotics were removed, the control group was inoculated with DPBS and the experimental group was inoculated with ETBF. On day 0, both groups underwent volumetric muscle loss surgery with PCL implantation in the defect space. Data was collected 6 week after implantation. Stromal (singlet, live, CD31-/CD45-, CD29+) cells were sorted from the harvested quadriceps and subsequently underwent bulk RNA sequencing.

### Repository Contents

The aligned counts files are available in the data directory of this repository, and the results of this analysis are included in the results directory. All scripts required for the analysis are included in the scripts directory. This project used the `renv` package to ensure a reproducible environment, which can be used as:

```r
install.packages("renv")
renv::restore()
```

The analysis scripts have been numbered in the appropriate order for use:

```r
source("scripts/utils.R")
source("scripts/1_processing.R")
source("scripts/2_dge.R")
source("scripts/3_gsea.R")
```

### Results and Figures

These scripts generate data used in the manuscript figures:

 - Figure 3B is a volcano plot created from 2_ETBF_pos_v_neg_dge.csv
 - Figure 3C contains manually selected genes based on function, subsetted from 2_ETBF_pos_v_neg_dge.csv
 - Figure 3D reflects the results of 3_reactome_gsea.csv
 - Figure S10 contains a miscellaneous category of significantly differentially expressed genes from 2_ETBF_pos_v_neg_dge.csv

Figures were created using these exported results.

### Citation

If you use this data or code, please cite our publication:

Yang, B. et al. Murine gut microbiota dysbiosis via enteric infection modulates the foreign body response to a distal biomaterial implant. Proceedings of the National Academy of Sciences 122, e2422169122 (2025).


### Contact

For questions, comments, or concerns, please contact Kavita Krishnan at [kkrishnan@jhmi.edu](mailto:kkrishnan@jhmi.edu) or open an issue in this repository.
