# DISC-seq
Cell heterogeneity is a fundamental feature of biological systems, driving diverse responses to stimuli and stressors, including developmental cues, diseases, and drug treatments. While single-cell RNA sequencing (scRNA-seq) has revolutionized our ability to characterize this diversity by profiling gene expression at single cell levels, a critical gap remains in that it cannot directly link transcriptional profiles to functional cellular outcomes, such as stress-induced damage. To bridge this gap, we developed DISC-seq (Damage Identification in Single-Cell RNA sequencing), a method compatible with standard scRNA-seq workflows that simultaneously quantifies transcriptome-wide gene expression and evaluates the extent of cell damage at single-cell resolution. Applied to both cancer cell lines and clinical peripheral blood mononuclear cells (PBMCs) from pediatric hematology patients, DISC-seq uncovered key molecular pathways and gene expression determinants that govern heterogeneous treatment and stress responses. Our approach enables the systematic discovery of regulatory mechanisms underlying heterogeneous cellular stress sensitivity within and across cell types, providing a powerful tool for dissecting the molecular basis of cell heterogeneity.

This repo contains essential scripts used for implementation of DISC-seq as well as codes to reproduce analysis results in the referred publication. We used MGI DNBelab C-series High-throughput Single-cell RNA Library Preparation Kit V3.0 for scRNA-seq library preparation and sequencing. With modification, these scripts can be applied for other conventional scRNA-seq libraries and sequencing data.

* Script *1.annx-r2-bgilist-1id2CB.sh* is used for demultiplexing of DISC-seq barcodes from raw sequencing data.
* Optionally, if CASB sample pooling is used, script *2.CASB-r2-bgilist-1id2cB.sh* is used for demultiplexing of CASB barcodes from raw sequencing data.
* Script *endo-merged-chip4.ipynb* is used for integrating DISC-seq and CASB information with endogenous expression matrix.
* Supplement scripts *whitelist.py* and *whitelist_methods.py* are used in *1.annx-r2-bgilist-1id2CB.sh* and *2.CASB-r2-bgilist-1id2cB.sh* for getting barcode and UMI information of CASB and DISC reads.
* Supplement file *anx5-barcodes.txt* is DISC-seq barcode information.
* Supplement file *CASB-barcode-example.txt* is CASB barcodes information. NNNNNNNN should be filled by corresponding sample barcodes provided in our paper.
* Supplement file *BGI_droplet_scRNA_readStructureV2_cDNA_T1-2.sort.zip* is the BGI cell barcode whitelist file used for filtering out incorrect cell barcode. Uncompress the file before use.

Raw and processed data can be retrieved with the following links:
* DISC-seq and CASB raw sequencing files: https://doi.org/10.17632/v9kdgw7nrh.2
* Processed matrices of scRNA-seq data: https://doi.org/10.17632/yd3jfm2g3j.1    
* Raw sequencing data of scRNA-seq: https://ngdc.cncb.ac.cn/gsa-human/, accession number HRA014683 and HRA016180

Reference
* Hu Q, Wang Y, Zhao Y, Kong L, Tang X, Lin Q, Zhou Y, Wang Y, Wang H, Jiang H, Luo X, He J, Liu S, Hu Y. DISC-seq: deciphering cell stress heterogeneity through joint mapping of cellular damage and transcriptomic landscapes in scRNA-seq. *In revision*
* Fang L, Li G, Sun Z, Zhu Q, Cui H, Li Y, Zhang J, Liang W, Wei W, Hu Y, Chen W. CASB: a concanavalin A-based sample barcoding strategy for single-cell sequencing. Mol Syst Biol. 2021 Apr;17(4):e10060. doi: 10.15252/msb.202010060IF: 7.7 Q1 . PMID: 33821571; PMCID: PMC8022202.
