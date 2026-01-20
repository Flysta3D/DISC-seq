DISC-seq Pipeline
========================
Cell heterogeneity is a fundamental feature of biological systems, driving diverse responses to stimuli and stressors, including developmental cues, diseases, and drug treatments. While single-cell RNA sequencing (scRNA-seq) has revolutionized our ability to characterize this diversity by profiling gene expression at single cell levels, a critical gap remains in that it cannot directly link transcriptional profiles to functional cellular outcomes, such as stress-induced damage. To bridge this gap, we developed DISC-seq (Damage Identification in Single-Cell RNA sequencing), a method compatible with standard scRNA-seq workflows that simultaneously quantifies transcriptome-wide gene expression and evaluates the extent of cell damage at single-cell resolution. Applied to both cancer cell lines and clinical peripheral blood mononuclear cells (PBMCs) from pediatric hematology patients, DISC-seq uncovered key molecular pathways and gene expression determinants that govern heterogeneous treatment and stress responses. Our approach enables the systematic discovery of regulatory mechanisms underlying heterogeneous cellular stress sensitivity within and across cell types, providing a powerful tool for dissecting the molecular basis of cell heterogeneity.

This repo contains essential scripts used for implementation of DISC-seq as well as codes to reproduce analysis results in the referred publication. We used **MGI DNBelab C-series High-throughput Single-cell RNA Library Preparation Kit V3.0** for scRNA-seq library preparation and sequencing. With modification, these scripts can be applied for other conventional scRNA-seq libraries and sequencing data.

## DISC-seq Overview:

<img width="3086" height="766" alt="image" src="https://github.com/user-attachments/assets/587e90c6-c9cf-448f-9038-39e137a6e961" />

## This repo includes:
* **Figure_Code** including codes for generating figures in the referred publication.
* **Pipeline_Code** including **DISC-Seq Data Demultiplexing and Integration Pipeline**.
* **Supplement** including example files:
```
Supplement/
├── ANXV barcode reference: anx5-barcodes.txt
├── CASB barcode reference: CASB-barcode-example.txt
├── config file example: config_template.txt
├── barcode translate file example: barcodeTranslate-example.txt
└── cellbarcode whitelist reference: BGI_droplet_scRNA_readStructureV2_cDNA_T1-2.sort.zip (unzip before use).
```

## Pipeline requires these data:
* DISC-seq library​ for ANXV label counting.
* (Optional) CASB library​ for sample demultiplexing.
* Endogenous RNA expression matrix​ for gene expression analysis.

The pipeline performs quality control, barcode extraction, sample demultiplexing, and integrates all data into a Seurat object for downstream analysis.

# Quick Start
## 1. Clone Repository
```Bash
git clone git@github.com:Flysta3D/DISC-seq.git
cd DISC-seq/Pipeline_Code
```
## 2. Install Dependencies
### System Requirements
* Linux or macOS
* Python 3.7+
* R 4.0+
### Required Tools
Tool|Version|Installation
|--|--|--|
|fastp|1.0.1|https://github.com/OpenGene/fastp|
|fastq-multx|1.4.3|https://github.com/brwnj/fastq-multx|

### Python Packages
```Bash
pip install numpy matplotlib scipy regex umi_tools
```
### R Packages
```Bash
install.packages(c("Seurat", "ggplot2", "data.table", "dplyr", 
                   "reshape2", "optparse", "stringr"))
```
## 3. Configure the Pipeline
```Bash
# Copy configuration template
cp ../Supplement/config_template.txt my_config.txt

# Edit configuration
nano my_config.txt
```
## 4. Run the Pipeline
```Bash
# Test run (dry run)
bash ./DISC_pipeline.sh -c my_config.txt --dry-run

# Full analysis
bash ./DISC_pipeline.sh -c my_config.txt
```
#### Debug Mode
Run with debug mode for detailed output:
```Bash
bash ./DISC_pipeline.sh -c config.txt --debug
```

# Input/Output
## Input Structure
```
data/
├── ANXV_library/                    # ANXV library FASTQ files
│   ├── [chip_name]
│   │   ├── ANXV_R1.fastq.gz
│   │   └── ANXV_R2.fastq.gz
│   └── ...     # ANXV library for other chips
├── CASB_library/                    # CASB library FASTQ files
│   ├── [chip_name]
│   │   ├── CASB_R1.fastq.gz
│   │   └── CASB_R2.fastq.gz
│   └── ...     # CASB library for other chips
├── rna_matrix/                      # [chip_name] RNA matrix
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
├── ...     # rna_matrix for other chips
└── library reference/               # Library reference files
    ├── Cellbarcode_whitelist.txt    # Cell barcode whitelist
    ├── barcodeTranslate.txt         # Barcode translation table
    ├── anx5-barcodes.txt            # ANXV barcode list
    └── CASB-barcodes.txt            # CASB barcode list
```

## Output Structure
```
results_YYYYMMDD_HHMMSS/
├── ANXV_lib/                                             # Processed ANXV libraries
│   └── [chip_name]/
│       ├── CellBarcodes/                                 # Barcode lists and cutoff curve picture
│       │       ├── comm.ANX5-labeled_whitelist.txt       # After barcode filtering
│       │       ├── ANX5-labeled_whitelist.txt            # Before barcode filtering
│       │       └── ANX5-labeled_expect_whitelist_cell_barcode_count_cutoff.png
│       ├── [chip_name]_ANXV_stats.txt                    # Summary statistics
│       └── tmp.*/                                        # Temporary files (Will be removed after running)
├── CASB_lib/                                             # Processed CASB libraries (Similar as ANXV libraries)
│   └── [chip_name]/
│       ├── CellBarcodes/
│       │       ├── comm.[sample_name]_whitelist.txt
│       │       ├── [sample_name]_whitelist.txt
│       │       ├── [sample_name]_expect_whitelist_cell_barcode_count_cutoff.png
│       │       └── ...
│       ├── [chip_name]_CASB_stats.txt
│       └── tmp.*/              
├── merged/                                               # Integrated analysis results
│   ├── final_seurat_object.rds                           # Seurat object ***
│   ├── cell_metadata.csv                                 # Cell metadata ***
│   ├── cell_statistics.csv                               # Summary statistics
│   ├── sample_statistics.csv                             # Sample counts
│   ├── analysis_parameters.txt                           # Analysis parameters
│   ├── r_analysis_report.txt                             # R analysis report
│   ├── analysis_report.md                                # Summary report
│   └── qc_plots/                                         # Quality control plots
│       ├── sample_distribution.png
│       ├── hto_heatmap.png
│       ├── qc_correlation.png
│       ├── anxv_distribution.png
│       └── qc_violin.png
├── qc_reports/                                           # QC reports
│   ├── [chip_name].html
│   └── [chip_name].json
├── logs/                                                 # Log files
│   ├── process.log
│   ├── error.log
│   ├── fastp_[chip_name]_ANXV/CASB.log
│   ├── fastq_multx_[chip_name]_ANXV/CASB.log
│   ├── whitelist_[chip_name]_ANXV/CASB.log
│   ├── umi_tools_[chip_name]_ANXV/CASB.log
│   └── R_analysis.log
└── pipeline_config.txt                                   # Pipeline configuration
```

## Key Output Files
|File|Description|
|--|--|
|final_seurat_object.rds|Seurat object with integrated data|
|cell_metadata.csv|Complete cell metadata including HTO classification|
|analysis_report.md|Summary of analysis results|
|hto_heatmap.png|Heatmap of HTO counts|

### Seurat Object Metadata
The Seurat object contains the following metadata columns:
|Column|Description|
|--|--|
|HTO_classification|Sample assignment (Sample1/Sample2...)|
|HTO_classification.global|Singlet assignment (Only singlet kept)|
|HTO_maxID|Maximum HTO identity|
|HTO_maxCount|Maximum HTO count|
|HTO_secondID|Second highest HTO identity|
|HTO_secondCount|Second highest HTO count|
|HTO_magrin|The difference between signals for hash.maxID and hash.secondID|
|hash.ID|Classification result where doublet IDs are collapsed|
|ANXV|ANX5 labeling counts|
|log2ANXV|log2 transformed ANX5 counts|
|log10ANXV|log10 transformed ANX5 counts|
|percent.mt|Percentage of mitochondrial genes|

# Loading Results in R
```R
# Load the Seurat object
library(Seurat)
seurat_obj <- readRDS("path/to/final_seurat_object.rds")

# View metadata
head(seurat_obj@meta.data)

# Visualize results
DimPlot(seurat_obj, group.by = "HTO_classification")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"))

# Downstream analysis
seurat_obj <- NormalizeData(seurat_obj)

# ...
```

# Data Availability
Raw and processed data in the referred publication can be retrieved with the following links:
* DISC-seq and CASB raw sequencing files: https://doi.org/10.17632/v9kdgw7nrh.2
* Processed matrices of scRNA-seq data: https://doi.org/10.17632/yd3jfm2g3j.1    
* Raw sequencing data of scRNA-seq: https://ngdc.cncb.ac.cn/gsa-human/, accession number HRA014683 and HRA016180

# Reference
* Hu Q, Wang Y, Zhao Y, Kong L, Tang X, Lin Q, Zhou Y, Wang Y, Wang H, Jiang H, Luo X, He J, Liu S, Hu Y. DISC-seq: deciphering cell stress heterogeneity through joint mapping of cellular damage and transcriptomic landscapes in scRNA-seq. *In revision*
* Fang L, Li G, Sun Z, Zhu Q, Cui H, Li Y, Zhang J, Liang W, Wei W, Hu Y, Chen W. CASB: a concanavalin A-based sample barcoding strategy for single-cell sequencing. Mol Syst Biol. 2021 Apr;17(4):e10060. doi: 10.15252/msb.202010060. PMCID: PMC8022202.
