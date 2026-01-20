#!/usr/bin/env Rscript

# DISC-Seq data merging and demultiplexing script
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(optparse)
  library(stringr)
})

# Define command-line options
option_list <- list(
  make_option(c("--ANXV_dir"), type="character", default=NULL, 
              help="ANXV library directory", metavar="DIR"),
  make_option(c("--CASB_dir"), type="character", default=NULL,
              help="CASB library directory", metavar="DIR"),
  make_option(c("--rna_matrix"), type="character", default=NULL,
              help="RNA matrix directory", metavar="DIR"),
  make_option(c("--barcode_translate"), type="character", default=NULL,
              help="Cell barcode translation file", metavar="FILE"),
  make_option(c("--output_dir"), type="character", default="./merged_results",
              help="Output directory", metavar="DIR"),
  make_option(c("--quantile"), type="numeric", default=0.99,
              help="HTODemux quantile threshold", metavar="NUM"),
  make_option(c("--output_rds"), type="character", default="final_seurat_object.rds",
              help="Output RDS filename", metavar="FILE"),
  make_option(c("--install_packages"), type="logical", default=FALSE,
              help="R package installation flag (true/false)", metavar="BOOL")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$ANXV_dir) || is.null(opt$CASB_dir) || 
    is.null(opt$rna_matrix) || is.null(opt$barcode_translate)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}


# Create output directory
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Record analysis parameters
writeLines(c(
  paste("Analysis time:", Sys.time()),
  paste("ANXV library directory:", opt$ANXV_dir),
  paste("CASB library directory:", opt$CASB_dir),
  paste("RNA matrix directory:", opt$rna_matrix),
  paste("Barcode translation file:", opt$barcode_translate),
  paste("Output directory:", opt$output_dir),
  paste("HTODemux quantile threshold:", opt$quantile),
  paste("Output RDS:", opt$output_rds),
  ""
), file.path(opt$output_dir, "analysis_parameters.txt"))

# Load ANXV library data
load_ANXV_library <- function(ANXV_dir) {
  cat("Load ANXV library data...\n")
  
  # Find all ANXV whitelist files
  ANXV_files <- list.files(ANXV_dir, pattern = "comm\\.ANX5-labeled_whitelist\\.txt", 
                           recursive = TRUE, full.names = TRUE)
  
  if (length(ANXV_files) == 0) {
    stop("Could not find ANXV library result files")
  }
  
  # Read all ANXV library results
  ANXV_data <- data.frame()
  for (f in ANXV_files) {
    chip_name <- basename(dirname(dirname(f)))
    temp <- fread(f, header = FALSE, sep = "\t")
    if (ncol(temp) >= 3) {
      temp <- temp[, c(1, 3)]  # Extract CellBarcodes and Count
      colnames(temp) <- c("CellBarcodes", "Count")
      temp$Sample <- chip_name
      ANXV_data <- rbind(ANXV_data, temp)
    }
  }
  
  cat("Find ", nrow(ANXV_data), " ANXV barcodes\n")
  return(ANXV_data)
}

# Load CASB library data
load_CASB_library <- function(CASB_dir) {
  cat("Load CASB library data...\n")
  
  # Find all CASB whitelist files
  CASB_files <- list.files( CASB_dir, pattern = "^comm\\..*_whitelist\\.txt$", 
                            recursive = TRUE, full.names = TRUE)
  
  if (length(CASB_files) == 0) {
    stop("Could not find CASB library result files")
  }
  
  # Read all CASB library results
  CASB_data <- data.frame()
  for (f in CASB_files) {
    # Extract sample name from filename
    fname <- basename(f)
    sample_name <- gsub("^comm\\.(.+)_whitelist\\.txt$", "\\1", fname)
    
    temp <- fread(f, header = FALSE, sep = "\t")
    if (ncol(temp) >= 3) {
      temp <- temp[, c(1, 3)]  #  Extract CellBarcodes and Count
      colnames(temp) <- c("CellBarcodes", "Count")
      temp$Sample <- sample_name
      CASB_data <- rbind(CASB_data, temp)
    }
  }
  
  cat("Find ", nrow(CASB_data), " CASB barcodes\n")
  return(CASB_data)
}

# Main analysis function
main_analysis <- function(ANXV_data, CASB_data, rna_matrix, barcode_translate, quantile, outdir) {
  cat("Start main analysis...\n")
  
  # Read RNA matrix
  cat("Read RNA matrix...\n")
  rna_data <- tryCatch({
    Read10X(data.dir = rna_matrix, gene.column = 1, unique.features=TRUE)
  }, error = function(e) {
    stop("Failed to read RNA matrix: ", e$message)
  })
  if (ncol(rna_data) == 0 || nrow(rna_data) == 0) {
    stop("RNA matrix is empty")
  }
  cat("RNA matrix dimension:", dim(rna_data), "\n")
  cat("Cell number:", ncol(rna_data), "\n")
  cat("Gene number:", nrow(rna_data), "\n")
  
  # Read barcode translation table
  cat("Read barcode translation table...\n")
  cell_barcodes <- fread(barcode_translate, header = FALSE, sep = "\t")
  colnames(cell_barcodes) <- c("CellBarcodes", "Cell")
  
  # Filter barcodes present in RNA data
  cell_barcodes <- cell_barcodes[Cell %in% colnames(rna_data), ]
  cat("Number of annotated cell barcodes in endogeneous library:", nrow(cell_barcodes), "\n")
  cat("Number of detected cell barcodes in CASB library:", length(unique(CASB_data$CellBarcodes)), "\n")
  
  # Merge barcode information
  cat("Merging barcode information...\n")
  merged_barcodes <- merge(CASB_data, cell_barcodes, by = "CellBarcodes")
  cat("Number of intersected cell barcodes in both endogeneous and CASB library:", nrow(merged_barcodes), "\n")
  
  # Get unique intersected barcodes
  joint_cbs <- unique(merged_barcodes$CellBarcodes)
  cat("Number of unique intersected CBs: ", length(joint_cbs), "\n")
  
  # Filter RNA data to keep only intersected barcodes
  cell_barcodes <- cell_barcodes[CellBarcodes %in% joint_cbs, ]
  rna_data <- rna_data[, colnames(rna_data) %in% cell_barcodes$Cell]
  
  # Create HTO matrix
  cat("Creating HTO matrix...\n")
  hto_matrix <- dcast(merged_barcodes[, c("Sample", "CellBarcodes", "Count")], 
                      Sample ~ CellBarcodes, value.var = "Count", fill = 0)
  rownames(hto_matrix) <- hto_matrix$Sample
  hto_matrix$Sample <- NULL
  
  
  # Reorder HTO matrix columns to match cell barcodes
  common_cells <- intersect(cell_barcodes$CellBarcodes, colnames(hto_matrix))
  hto_matrix <- hto_matrix[, common_cells, drop = FALSE]
  cell_barcodes <- cell_barcodes[CellBarcodes %in% common_cells, ]

  # Create Seurat object
  cat("Creating Seurat object...\n")
  idx = match(cell_barcodes$Cell, colnames(rna_data))
  rna_data.umis <- rna_data[, idx, drop = FALSE]
  colnames(rna_data.umis) = cell_barcodes$CellBarcodes

  seurat_obj <- CreateSeuratObject(
    counts =  rna_data.umis,
    project = "DISC_Merge"
  )
  
  # Add cell barcode metadata
  seurat_obj$CellBarcodes <- cell_barcodes$CellBarcodes
  seurat_obj$CellNames <- colnames(rna_data)[idx]

  # Add HTO assay
  cat("Adding HTO assay...\n")
  # Ensure HTO matrix columns match Seurat object cell order
  idx = match(cell_barcodes$CellBarcodes, colnames(hto_matrix))
  hto_matrix = hto_matrix[, idx, drop = FALSE]
  colnames(hto_matrix) =  cell_barcodes$CellBarcodes

  seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto_matrix)
  
  # Create pseudoHTO assay for demultiplexing (adding 1 to avoid zeros)
  cat("Adding pseudoHTO assay...\n")
  seurat_obj[["pseudoHTO"]] <- CreateAssayObject(
    counts = as.matrix(hto_matrix) + 1
  )
  
  # HTODemux analysis
  cat("Running HTODemux (quantile threshold:", quantile,")...\n")
  seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
  seurat_obj <- NormalizeData(seurat_obj, assay = "pseudoHTO", normalization.method = "CLR")
  
  seurat_obj <- HTODemux(
    seurat_obj,
    assay = "pseudoHTO",
    positive.quantile = quantile,
    nstarts = 100
  )
  
  # Rename HTO metadata columns
  colnames(seurat_obj@meta.data) <- gsub("pseudoHTO_", "HTO_", 
                                         colnames(seurat_obj@meta.data))

  # HTO heatmap
  cat("Generating HTO heatmap...\n")
  hto_heatmap <- HTOHeatmap(seurat_obj, assay = "HTO") +
    labs(title = paste("HTO Heatmap (Quantile =", quantile, ")"))
    
  dir.create(file.path(outdir, "qc_plots"), showWarnings = FALSE)
  ggsave(file.path(outdir, "qc_plots/hto_heatmap.png"), hto_heatmap, 
         width = 12, height = 8)

  # Save cell statistics before singlet extraction
  hto_stats <- data.frame(
    Singlets = sum(seurat_obj$HTO_classification.global == "Singlet"),
    Doublets = sum(seurat_obj$HTO_classification.global == "Doublet"),
    Negative = sum(seurat_obj$HTO_classification.global == "Negative")
  )
  
  write.csv(hto_stats, 
            file.path(outdir, "hto_statistics.csv"), 
            row.names = FALSE)

  # Extract singlets
  cat("Extracting singlet cells...\n")
  singlet_cells <- WhichCells(seurat_obj, expression = HTO_classification.global == "Singlet")
  seurat_obj.singlet <- subset(seurat_obj, cells = singlet_cells)

  seurat_obj.singlet = seurat_obj.singlet[, !duplicated(colnames(seurat_obj.singlet))]
  seurat_obj.singlet = seurat_obj.singlet[!duplicated(rownames(seurat_obj.singlet)), ]
  metadf = seurat_obj.singlet@meta.data
  metadf %>% filter(HTO_classification.global == 'Singlet') %>%
    group_by(CellNames) %>%
    filter(nCount_HTO == max(nCount_HTO)) %>%
    ungroup() %>%
    filter(!duplicated(CellNames)) %>% as.data.frame() -> singlet.metadf
  seurat_obj.singlet = seurat_obj.singlet[, colnames(seurat_obj.singlet) %in% singlet.metadf$CellBarcodes]
  rownames(singlet.metadf) = singlet.metadf$CellNames
  colnames(seurat_obj.singlet) = singlet.metadf$CellNames
  Idents(seurat_obj.singlet) <- "HTO_classification.global"

  
  # Add ANXV data
  cat("Adding ANXV data...\n")
  ANXV_merged <- merge(ANXV_data, singlet.metadf, by = "CellBarcodes", all.y = T, sort=F)
  ANXV_wide <- dcast(ANXV_merged, Sample ~ CellBarcodes, 
                       value.var = "Count", fill = 0)
  ANXV_wide$Sample <- NULL
  
  # Reorder ANXV matrix columns to match Seurat object cell order
  idx = match( seurat_obj.singlet$CellBarcodes, ANXV_wide$CellBarcodes)
  seurat_obj.singlet$ANXV <- as.numeric(ANXV_wide[1, ])
  seurat_obj.singlet$log2ANXV <- log2(seurat_obj.singlet$ANXV + 1)
  seurat_obj.singlet$log10ANXV <- log10(seurat_obj.singlet$ANXV + 1)
  
  # Calculate percentage of mitochondrial genes
  cat("Calculating percentage of mitochondrial genes...\n")
  seurat_obj.singlet[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj.singlet, 
    pattern = "^MT-"
  )

  
  cat("Cell number after processing: ", ncol(seurat_obj.singlet), "\n")
  cat("Gene number after processing: ", nrow(seurat_obj.singlet), "\n")
  
  return(seurat_obj.singlet)
}

# Generate QC plots
generate_qc_plots <- function(seurat_obj, output_dir) {
  cat("Generating QC plots...\n")
  
  # Create QC plots directory
  qc_dir <- file.path(output_dir, "qc_plots")
  dir.create(qc_dir, showWarnings = FALSE)
  
  # Sample distribution
  sample_plot <- ggplot(seurat_obj@meta.data, 
                        aes(x = HTO_classification, fill = HTO_classification)) +
    geom_bar() +
    geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, color = "white")
    theme_minimal() +
    labs(title = "Sample Distribution", 
         x = "Sample", 
         y = "Number of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(qc_dir, "sample_distribution.png"), sample_plot, 
         width = 8, height = 6)
  
  # QC correlation plot
  qc_cor <- ggplot(seurat_obj@meta.data, 
                   aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(title = "QC Metrics Correlation", 
         x = "UMI Counts", 
         y = "Gene Counts", 
         color = "MT%")
  
  ggsave(file.path(qc_dir, "qc_correlation.png"), qc_cor, 
         width = 8, height = 6)
  
  # ANXV distribution plot
  anxv_plot <- ggplot(seurat_obj@meta.data, aes(x = log2ANXV)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    labs(title = "ANXV Expression Distribution", 
         x = "log2ANXV", 
         y = "Number of Cells")

  ggsave(file.path(qc_dir, "anxv_distribution.png"), anxv_plot, 
         width = 8, height = 6)
  
  # QC violin plots
  qc_violin <- VlnPlot(seurat_obj, 
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3, pt.size = 0)

  ggsave(file.path(qc_dir, "qc_violin.png"), qc_violin, 
         width = 12, height = 5)
}

# Save results
save_results <- function(seurat_obj, output_dir, output_rds) {
  cat("Saving results...\n")
  
  # Save Seurat object
  saveRDS(seurat_obj, file.path(output_dir, output_rds))
  cat("Saving Seurat object to:", file.path(output_dir, output_rds), "\n")
  
  # Save cell metadata
  write.csv(seurat_obj@meta.data, 
            file.path(output_dir, "cell_metadata.csv"), 
            row.names = TRUE)
  
  # Save cell statistics
  cell_stats <- data.frame(
    Total_Cells = ncol(seurat_obj),
    Total_Genes = nrow(seurat_obj),
    Mean_UMI = mean(seurat_obj$nCount_RNA),
    Mean_Genes = mean(seurat_obj$nFeature_RNA),
    Mean_MT = mean(seurat_obj$percent.mt)
  )
  
  write.csv(cell_stats, 
            file.path(output_dir, "cell_statistics.csv"), 
            row.names = FALSE )
  
  # Save sample statistics
  sample_stats <- as.data.frame(table(seurat_obj$HTO_classification))
  colnames(sample_stats) <- c("Sample", "Cell_Count")
  write.csv(sample_stats, 
            file.path(output_dir, "sample_statistics.csv"), 
            row.names = FALSE)
  
  # Generate analysis report
  generate_report(seurat_obj, output_dir)
}

# Generate analysis report
generate_report <- function(seurat_obj, output_dir) {
  cat("Generating analysis report...\n")
  
  report_file <- file.path(output_dir, "r_analysis_report.txt")
  hto_counts <- read.csv(file.path(output_dir, "hto_statistics.csv"))
  
  report_content <- c(
    "= DISC-Seq Analysis Report =",
    paste("Generate time:", Sys.time()),
    paste("Seurat package version:", packageVersion("Seurat")),
    "",
    "1. Data overview",
    paste("   - Total cell number:", ncol(seurat_obj)),
    paste("   - Total gene number:", nrow(seurat_obj)),
    paste("   - Sample number:", length(unique(seurat_obj$HTO_classification))),
    "",
    "2. Cell quality statistics",
    paste("   - Average UMI:", round(mean(seurat_obj$nCount_RNA), 2)),
    paste("   - Average Gene number:", round(mean(seurat_obj$nFeature_RNA), 2)),
    paste("   - Average percentage of mitochondrial genes:", round(mean(seurat_obj$percent.mt), 2), "%"),
    "",
    "3. HTO classification statistics",
    paste("   - Singlets:", hto_counts$Singlets),
    paste("   - Doublets:", hto_counts$Doublets),
    paste("   - Negative:", hto_counts$Negative),
    "",
    "4. ANXV statistics",
    paste("   - Average ANXV counts:", round(mean(seurat_obj$ANXV), 2)),
    paste("   - Cells with ANXV labeled:", sum(seurat_obj$ANXV > 0)),
    "",
    "5. Files generated:",
    paste("   - Seurat object:", file.path(output_dir, basename(opt$output_rds))),
    paste("   - Cell metadata: cell_metadata.csv"),
    paste("   - Cell statistics: cell_statistics.csv"),
    paste("   - HTO statistics: hto_statistics.csv"),
    paste("   - Sample statistics: sample_statistics.csv"),
    paste("   - QC plots in: qc_plots/")
  )
  
  writeLines(report_content, report_file)
  cat("Report file is saved in:", report_file, "\n")
}

# Main function
main <- function() {
  # Install required packages if needed
  if (opt$install_packages) {
    cat("Installing required R packages...\n")
    packages <- c("Seurat", "ggplot2", "data.table", "dplyr", 
                  "reshape2", "optparse", "stringr")
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }
  
  # Load library data
  ANXV_data <- load_ANXV_library(opt$ANXV_dir)
  CASB_data <- load_CASB_library(opt$CASB_dir)
  
  # Run main analysis
  seurat_obj <- main_analysis(
    ANXV_data = ANXV_data,
    CASB_data = CASB_data,
    rna_matrix = opt$rna_matrix,
    barcode_translate = opt$barcode_translate,
    quantile = opt$quantile,
    outdir = opt$output_dir
  )
  
  # Generate QC plots
  generate_qc_plots(seurat_obj, opt$output_dir)
  
  # Save results
  save_results(seurat_obj, opt$output_dir, basename(opt$output_rds))
  
  cat("\nAnalysis completed successfully!\n")
  cat("Output files saved in:", opt$output_dir, "\n")
}

# Run main function
main()