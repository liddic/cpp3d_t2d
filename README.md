# cpp3d_t2d
Examining overlaps and hypothetical contributions of altered functional capacities in degraded ecosystem soil microbiomes to metabolic anomalies encoded in type 2 diabetes (T2D) gut microbiomes

**Contains code for the article:**

*Degraded ecosystem soil and type 2 diabetes gut microbiomes share altered potential metabolism for sugars, lignin and branched-chain fatty acids: a blind spot for global health?*

by Craig Liddicoat, Bart A. Eijkelkamp, Timothy R. Cavagnaro, Jake M. Robinson, Kiri Joy Wallace, Andrew D. Barnes, Garth Harmsworth, Damien J. Keating, Robert A. Edwards and Martin F. Breed

Preprint server: [BIORXIV/2025/642605](https://doi.org/10.1101/2025.03.11.642605)

**About this repository:** This repository documents code used in the above article, and is not intended for ongoing code development. Sub-folders contain separate workflow documentation and scripts relevant to the case study datasets described below.

**Data sources:** comprise gut metagenome data from type 2 diabetes versus normal healthy subjects from [Forslund et al 2015](https://www.nature.com/articles/nature15766), based on earlier cohort data published by [Karlsson et al 2013](https://www.nature.com/articles/nature12198) from Sweden and [Qin et al 2012](https://www.nature.com/articles/nature11450) from China; together with soil metagenome data from a post-mining ecosystem restoration case study data published by [Sun and Badgley 2019](https://www.sciencedirect.com/science/article/abs/pii/S0038071719301385?via%3Dihub). These data are available from the NCBI [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) accessions PRJEB1786 and PRJNA422434, and [MG-RAST](https://www.mg-rast.org/) project mgp16379.


**Method overview:** This study uses microbial metagenome functional profiling at the resolution of individual compounds, as outlined in [cpp3d](https://github.com/liddic/cpp3d) and the [manuscript](https://doi.org/10.1101/2025.03.11.642605). Metagenomics data were processed in three steps: (i) raw sequences were accessed/downloaded, (ii) QA/QC was performed to obtain good sequences, then (iii) functional profiles were derived using [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS), all performed on Flinders University [DeepThought HPC](https://deepthoughtdocs.flinders.edu.au/en/latest/) linux high performance computer. SUPER-FOCUS results were downloaded to a local machine for further visualisation and analysis via [R code](cpp3d-t2d-R-code-June2025.R). SUPER-FOCUS functional profile data were previously generated for the [restoration gradient (soils)](https://github.com/liddic/compound_potential/tree/main/sunbad-resto/sunbad-resto_3_superfocus_fxns) and the [Swedish T2D cohort (gut)](https://github.com/liddic/compound_potential/tree/main/forslund-t2d/ft2d_3_superfocus_fxns) case studies. A further [Chinese T2D cohort (gut)](forslund-t2d-chn/3_fxn_superfocus) case study has been added to confirm and consolidate overlapping potential metabolism trends. Note that folder/filepath structures used will need to be adjusted to run on other HPCs. Software versions used were: Python (v3.8.5), FastQC (v0.11.9), Snakemake (v5.22.0), Fastp (v0.23.2), Diamond (v0.9.19), SUPER-FOCUS (v0.0.0), R (v4.2.2).
