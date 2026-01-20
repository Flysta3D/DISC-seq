#!/bin/bash

# DISC-Seq Data Demultiplexing and Integration Pipeline
# Version：1.0.0

set -euo pipefail

# Default parameters
CONFIG_FILE=""
DEBUG=false
SCRIPT_DIR="./"

OUTPUT_DIR=""
LOG_DIR=""
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

# Logging functions
log_info() { echo -e "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_warn() { echo -e "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_error() { echo -e "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_debug() {
    [[ "${DEBUG:-false}" == "true" ]] || return 0
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "${timestamp} - $*"
}

# Show help message
show_help() {
    cat << EOF
DISC-Seq Data Demultiplexing and Integration Pipeline
Version: 1.0.0

Usage: $0 [Options]

Options:
  -c, --config FILE      Configure file (required)
  -o, --output DIR       Output dir (overrides config file)
  -q, --quantile NUM     HTODemux quantile value (overrides config file)
  -t, --threads NUM      Threads (overrides config file)
  --debug               Enable debug mode
  -h, --help            Show this help message

examples:
  $0 -c config.txt
  $0 -c config.txt -o ./results -q 0.95
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config) CONFIG_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -q|--quantile) HTODEMUX_QUANTILE="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --debug) DEBUG=true; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) log_error "Unknown arg: $1"; show_help; exit 1 ;;
    esac
done

# Check config file
if [ -z "$CONFIG_FILE" ]; then
    log_error "Must provide a configuration file!"
    show_help
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    log_error "Config file not found:$CONFIG_FILE"
    exit 1
fi

# Load configuration
log_info "Load config file: $CONFIG_FILE"
source "$CONFIG_FILE"
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")
LOG_DIR=$(realpath -m "$LOG_DIR")
ANXV_LIB_DIR=$(realpath -m "$ANXV_LIB_DIR")
CASB_LIB_DIR=$(realpath -m "$CASB_LIB_DIR")
RNA_MATRIX=$(realpath -m "$RNA_MATRIX")
CB_WHITELIST=$(realpath -m "$CB_WHITELIST")
BARCODE_TRANSLATE=$(realpath -m "$BARCODE_TRANSLATE")
ANXV_BARCODES=$(realpath -m "$ANXV_BARCODES")
CASB_BARCODES=$(realpath -m "$CASB_BARCODES")
SCRIPT_DIR=$(realpath -m "$SCRIPT_DIR")

# Override config with command line args
[ -n "${OUTPUT_DIR:-}" ] && OUTPUT_DIR="$OUTPUT_DIR"
[ -n "${HTODEMUX_QUANTILE:-}" ] && HTODEMUX_QUANTILE="$HTODEMUX_QUANTILE"
[ -n "${THREADS:-}" ] && THREADS="$THREADS"

# Set default values if not set
OUTPUT_DIR="${OUTPUT_DIR:-./results_$TIMESTAMP}"
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")
LOG_DIR="${LOG_DIR:-$OUTPUT_DIR/logs}"
THREADS="${THREADS:-8}"
HTODEMUX_QUANTILE="${HTODEMUX_QUANTILE:-0.99}"
KEEP_TEMP="${KEEP_TEMP:-false}"

# Check required tools
check_tools() {
    local tools=("fastp" "fastq-multx" "python3" "umi_tools" "Rscript")
    local missing=()
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing+=("$tool")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Cannot find the tool: ${missing[*]}"
        return 1
    fi
    
    log_info "All required tools are available."
    return 0
}

# Check R packages
check_r_packages() {
    log_info "Checking R packages..."
    
    local required_packages=("Seurat" "ggplot2" "data.table" "dplyr" "reshape2" "optparse" "stringr")
    local missing_packages=()
    
    for pkg in "${required_packages[@]}"; do
        if ! Rscript -e "if(!require('$pkg', quietly=TRUE)) q(status=1)" &>/dev/null; then
            missing_packages+=("$pkg")
        fi
    done
    
    if [ ${#missing_packages[@]} -gt 0 ]; then
        log_warn "Missing R packages: ${missing_packages[*]}"
        log_info "Attempting to install missing R packages..."
            
        Rscript -e "
        packages <- c('Seurat', 'ggplot2', 'data.table', 'dplyr', 'reshape2', 'optparse', 'stringr')
        for (pkg in packages) {
            if (!require(pkg, character.only = TRUE)) {
                    install.packages(pkg, repos='https://cloud.r-project.org')
            }
        }
        "
        if [ $? -eq 0 ]; then
            log_info "R packages installed successfully"
        else
            log_error "Failed to install R packages"
            return 1
        fi
    else
        log_info "All required R packages are available"
    fi
    return 0
}

# Check input files
validate_inputs() {
    log_info "Check input files..."
    
    local errors=0
    
    # Check input directories
    for dir in "$ANXV_LIB_DIR" "$CASB_LIB_DIR"; do
        if [ ! -d "$dir" ]; then
            log_error "Directory not found: $dir"
            ((errors++))
        fi
    done
    
    # Check required files
    for file in "$CB_WHITELIST" "$BARCODE_TRANSLATE" "$ANXV_BARCODES" "$CASB_BARCODES"; do
        if [ ! -f "$file" ]; then
            log_error "File not found: $file"
            ((errors++))
        fi
    done
    
    # Check RNA matrix files
    if [ ! -d "$RNA_MATRIX" ]; then
        log_error "RNA matrix directory not found: $RNA_MATRIX"
        ((errors++))
    else
        local required_files=("matrix.mtx" "features.tsv" "barcodes.tsv")
        for f in "${required_files[@]}"; do
            if [ ! -f "$RNA_MATRIX/$f" ] && [ ! -f "$RNA_MATRIX/${f}.gz" ]; then
                log_warn "RNA matrix lack file: $f"
            fi
        done
    fi
    
    # Check for FASTQ files in input directories
    local ANXV_files=$(find "$ANXV_LIB_DIR" -name "*R1*.fastq*" -o -name "*R1*.fq*" 2>/dev/null | head -5)
    local CASB_files=$(find "$CASB_LIB_DIR" -name "*R1*.fastq*" -o -name "*R1*.fq*" 2>/dev/null | head -5)
    
    if [ -z "$ANXV_files" ]; then
        log_error "Cannot find R1 FASTQ in ANXV directory: $ANXV_LIB_DIR"
        ((errors++))
    fi
    
    if [ -z "$CASB_files" ]; then
        log_error "Cannot find R1 FASTQ in CASB directory: $CASB_LIB_DIR"
        ((errors++))
    fi
    
    if [ $errors -gt 0 ]; then
        log_error "Find $errors input errors. Please check!"
        return 1
    fi
    
    log_info "Input files are valid."
    return 0
}

# Create output directories
create_directories() {
    log_info "Creating output directories..."
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$LOG_DIR"
    mkdir -p "$OUTPUT_DIR/ANXV_lib"
    mkdir -p "$OUTPUT_DIR/CASB_lib"
    mkdir -p "$OUTPUT_DIR/merged"
    mkdir -p "$OUTPUT_DIR/qc_reports"
    
    # Configure logging to files
    exec 1> >(tee -a "$LOG_DIR/process.log")
    exec 2> >(tee -a "$LOG_DIR/error.log")
    
    log_info "Output dir: $OUTPUT_DIR"
    log_info "Log dir: $LOG_DIR"
}

# Record configuration
log_configuration() {
    cat > "$OUTPUT_DIR/pipeline_config.txt" << EOF
# DISC-Seq Data Demultiplexing and Integration Pipeline
# Creating time: $(date)
# Command line: $0 $@

Input parameters:
- ANXV dir: $ANXV_LIB_DIR
- CASB dir: $CASB_LIB_DIR
- RNA matrix dir: $RNA_MATRIX
- Cell barcode whitelist: $CB_WHITELIST
- Barcode translation file: $BARCODE_TRANSLATE
- ANXV barcode file: $ANXV_BARCODES
- CASB barcode file: $CASB_BARCODES

Output parameters:
- Output dir: $OUTPUT_DIR
- Log dir: $LOG_DIR
- HTOdemux quantile value: $HTODEMUX_QUANTILE
- Threads number: $THREADS
- Keep temp file: $KEEP_TEMP

Pipeline parameters:
- ANXV_fix sequence pattern: $ANXV_FIX_SEQ_PATTERN
- CASB_fix sequence pattern: $CASB_FIX_SEQ_PATTERN
- Tools versions:
  $(fastp --version 2>/dev/null || echo "fastp: not found")
  $(umi_tools --version 2>/dev/null || echo "umi_tools: not found")
EOF
}


# Process ANXV library
process_ANXV_library() {
    local chip_name="$1"
    local fq1="$2"
    local fq2="$3"
    
    log_info "Process ANXV library: $chip_name"
    
    local sample_dir="$OUTPUT_DIR/ANXV_lib/$chip_name"
    mkdir -p "$sample_dir/CellBarcodes"
    
    cd "$sample_dir"
    
    # Create temporary directory
    local tmp_dir
    if ! tmp_dir=$(mktemp -d -p . "tmp.${chip_name}.XXXXXX"); then
        log_error "Failed to create temporary directory"
        return 1
    fi
    
    log_debug "Temporary directory created: $tmp_dir"    
    
    local fq1_tmp="$tmp_dir/${chip_name}_R1.fastq"
    local fq2_tmp="$tmp_dir/${chip_name}_R2.fastq"
    
    # Define cleanup function
    cleanup_temp() {
        if [ -n "${tmp_dir:-}" ] && [ -d "$tmp_dir" ]; then
            log_info "Cleaning up temporary directory: $tmp_dir"
            rm -rf "$tmp_dir" 2>/dev/null || true
        fi
    }
    trap cleanup_temp EXIT
  
    
    # Uncompress FASTQ files to temporary directory
    if [[ "$fq1" == *.gz ]]; then
        log_debug "Uncompress FASTQ files to temporary directory: $fq1"
        gunzip -c "$fq1" > "$fq1_tmp"
        gunzip -c "$fq2" > "$fq2_tmp"
    else
        # Create hard links in temporary directory
        log_debug "Create hard links in temporary directory: $fq1"
        ln "$fq1" "$fq1_tmp"
        ln "$fq2" "$fq2_tmp"
    fi
    
    # Run fastp quality control
    log_debug "Run fastp quality control"

    fastp -n 10 -q 35 -w "$THREADS" --detect_adapter_for_pe \
          -i "$fq1_tmp" -I "$fq2_tmp" \
        -j "$OUTPUT_DIR/qc_reports/${chip_name}.json" -h "$OUTPUT_DIR/qc_reports/${chip_name}.html" -R "$chip_name" \
          2>&1 | tee "$LOG_DIR/fastp_${chip_name}_ANXV.log"
    
    # Run fastq-multx with ANX5 barcodes
    log_debug "Using ANX5 barcode to split reads"
    fastq-multx -B "$ANXV_BARCODES" -m 4 -b "$fq1_tmp" "$fq2_tmp" \
                -o $tmp_dir/%.R1.fastq -o $tmp_dir/%.R2.fastq \
                 2>&1 | tee "$LOG_DIR/fastq_multx_${chip_name}_ANXV.log"
    
    # Process ANX5-labeled reads
    process_barcoded_reads "ANX5-labeled" "$chip_name" "ANXV" "$ANXV_FIX_SEQ_PATTERN" "$tmp_dir"

    mv *whitelist.txt *.png *.txt CellBarcodes/ 2>/dev/null || true

    # Generate statistics
    generate_statistics "$chip_name" "$sample_dir" "ANXV" "$tmp_dir"

    log_info "Completed ANXV library process: $chip_name"
}

# Process CASB library
process_CASB_library() {
    local chip_name="$1"
    local fq1="$2"
    local fq2="$3"
    
    log_info "Processing CASB library: $chip_name"
    
    local sample_dir="$OUTPUT_DIR/CASB_lib/$chip_name"
    mkdir -p "$sample_dir/CellBarcodes"
    
    cd "$sample_dir"
    
    # Create temporary directory
    local tmp_dir
    if ! tmp_dir=$(mktemp -d -p . "tmp.${chip_name}.XXXXXX"); then
        log_error "Failed to create temporary directory"
        return 1
    fi
    
    local fq1_tmp="$tmp_dir/${chip_name}_R1.fastq"
    local fq2_tmp="$tmp_dir/${chip_name}_R2.fastq"
    
    # Define cleanup function
    cleanup_temp() {
        if [ -n "${tmp_dir:-}" ] && [ -d "$tmp_dir" ]; then
            log_info "Cleaning up temporary directory: $tmp_dir"
            rm -rf "$tmp_dir" 2>/dev/null || true
        fi
    }
    trap cleanup_temp EXIT
    
    
    # Uncompress FASTQ files to temporary directory
    if [[ "$fq1" == *.gz ]]; then
        log_debug "Uncompress FASTQ files to temporary directory: $fq1"
        gunzip -c "$fq1" > "$fq1_tmp"
        gunzip -c "$fq2" > "$fq2_tmp"
    else
        # Create hard links in temporary directory
        log_debug "Create hard links in temporary directory: $fq1"
        ln "$fq1" "$fq1_tmp"
        ln "$fq2" "$fq2_tmp"
    fi
    
    # Run fastp quality control
    log_debug "Run fastp quality control"

    fastp -n 10 -q 35 -w "$THREADS" --detect_adapter_for_pe \
          -i "$fq1_tmp" -I "$fq2_tmp" \
          -j "$OUTPUT_DIR/qc_reports/${chip_name}.json" -h "$OUTPUT_DIR/qc_reports/${chip_name}.html" -R "$chip_name" \
           2>&1 | tee "$LOG_DIR/fastp_${chip_name}_CASB.log"
    
    # Run fastq-multx with CASB barcodes
    log_debug "Using CASB barcode to split reads"
    fastq-multx -B "$CASB_BARCODES" -m 4 -b "$fq1_tmp" "$fq2_tmp" \
                -o $tmp_dir/%.R1.fastq -o $tmp_dir/%.R2.fastq  \
                2>&1 | tee "$LOG_DIR/fastq_multx_${chip_name}_CASB.log"
    
    # Get sample names from CASB barcode file
    local samples=$(cut -f1 "$CASB_BARCODES")
    
    # Process each sample
    for pre in $samples; do
        if [ -f "$tmp_dir/${pre}.R1.fastq" ]; then
            process_barcoded_reads "$pre" "$chip_name" "CASB" "$CASB_FIX_SEQ_PATTERN" "$tmp_dir"
        fi
    done
    
    # Merge sample results
    merge_sample_results "$chip_name"
    
    # Generate statistics
    generate_statistics "$chip_name" "$sample_dir" "CASB" "$tmp_dir"
    
    log_info "Completed CASB library process: $chip_name"
}

# Process barcoded reads
process_barcoded_reads() {
    local prefix="$1"
    local chip_name="$2"
    local lib_type="$3"
    local FIX_SEQ_PATTERN="$4"
    local tmp_dir="$5"
    
    log_debug "Process barcoded reads: $prefix"
    
    local nfq1="$tmp_dir/${prefix}.R1.fastq"
    local nfq2="$tmp_dir/${prefix}.R2.fastq"
    
    [ ! -f "$nfq1" ] && {
        log_warn "File not found: $nfq1"
        return
    }
    
    # Whitelist script to generate whitelist
    python3 "$SCRIPT_DIR/whitelist.py" \
            --stdin "$nfq1" \
            --read2-in="$nfq2" \
            --error-correct-threshold=0 \
            --method=umis \
            --extract-method=regex \
            --bc-pattern2=$FIX_SEQ_PATTERN \
            --plot-prefix="${prefix}_expect_whitelist" \
            --log2stderr \
            --expect-count=0 \
            --knee-method=None \
            > "${prefix}_whitelist.txt" 2>> "$LOG_DIR/whitelist_${chip_name}_${lib_type}.log"

    # Check whitelist file
    if [ ! -f "${prefix}_whitelist.txt" ]; then
        log_error "Failed to generate whitelist for $prefix"
        return 1
    fi
    
    # Filter whitelist
    join "${prefix}_whitelist.txt" "$CB_WHITELIST" 2>/dev/null | \
    awk '{print $1"\t\t"$2"\t"}' > "comm.${prefix}_whitelist.txt"
    
    # Extract true cell reads using umi_tools
    umi_tools extract \
            -I "$nfq1" \
            --read2-in "$nfq2" \
            --extract-method regex \
            --bc-pattern2 $FIX_SEQ_PATTERN \
            --stdout "$tmp_dir/truecell.${prefix}.1.fastq" \
            --read2-out "$tmp_dir/truecell.${prefix}.2.fastq" \
            --filtered-out "$tmp_dir/Nocell.${prefix}.1.fastq" \
            --filtered-out2 "$tmp_dir/Nocell.${prefix}.2.fastq" \
            --whitelist "comm.${prefix}_whitelist.txt" \
            2>> "$LOG_DIR/umi_tools_${chip_name}_${lib_type}.log"    
}

# Merge sample results
merge_sample_results() {
    local chip_name="$1"
    
    log_debug "Merge sample results: $chip_name"
    
    [ ! -f "CBlist_${chip_name}.txt" ] && {
        {
            echo -e "CellBarcodes\tCount\tSample"
            for whitelist in comm.*_whitelist.txt; do
                [ -f "$whitelist" ] || continue
                sample_name=$(basename "$whitelist" | sed 's/comm\.\(.*\)_whitelist\.txt/\1/')
                awk -v sample="$sample_name" 'BEGIN {OFS="\t"} {print $0, sample}' "$whitelist"
            done
        } > "CBlist_${chip_name}.txt"
        
        mv *whitelist.txt *.png *.txt CellBarcodes/ 2>/dev/null || true
    }
}

# Generate statistics
generate_statistics() {
    local chip_name="$1"
    local sample_dir="$2"
    local lib_type="$3"
    local tmp_dir="$4"
    
    log_debug "Generate statistics: $chip_name ($lib_type)"
        
    local outfile="${sample_dir}/${chip_name}_${lib_type}_stats.txt"
    

    # Initialize statistics file
    echo "===== Statistics: $chip_name ($lib_type) =====" > "$outfile"
    echo "Process time: $(date)" >> "$outfile"
    echo "==========================================" >> "$outfile"
        
    # Count raw reads
    local raw_reads=0
    if [ -f "$tmp_dir/${chip_name}_R1.fastq" ]; then
        raw_reads=$(wc -l < "$tmp_dir/${chip_name}_R1.fastq")
        raw_reads=$((raw_reads/4))
    fi
        
    echo "RawReads: $raw_reads" >> "$outfile"
    echo "" >> "$outfile"
        
    # Process statistics based on library type
    if [ "$lib_type" = "ANXV" ]; then
        # ANXV - handle ANX5-labeled
        process_anxv_stats "$chip_name" "$outfile" "$tmp_dir"
    elif [ "$lib_type" = "CASB" ]; then
        # CASB - handle CASB-labeled
        process_casb_stats "$chip_name" "$outfile" "$tmp_dir"
    fi
        
    # Add footer
    echo "==========================================" >> "$outfile"
    echo "Statistics file: $outfile" >> "$outfile"
        
    log_info "Statistics saved in: $outfile"
}
# Process ANXV statistics
process_anxv_stats() {
    local chip_name="$1"
    local outfile="$2"
    local tmp_dir="$3"
    
    local pre="ANX5-labeled"
    local b=0  # labeled reads
    local c=0  # cell reads
    local n_raw=0  # raw barcode types
    local n_filtered=0  # filtered barcode types
    

    echo "=== ANXV Statistics ===" >> "$outfile"

    # Check file existence
    if [ -f "$tmp_dir/${pre}.R1.fastq" ]; then
        # Count labeled reads
        b=$(wc -l < "$tmp_dir/${pre}.R1.fastq")
        b=$((b/4))
        
        # Count raw whitelist barcode types
        if [ -f "CellBarcodes/${pre}_whitelist.txt" ]; then
            n_raw=$(wc -l < "CellBarcodes/${pre}_whitelist.txt")
        fi
        
        # Count whitelist barcode types after filtering
        if [ -f "CellBarcodes/comm.${pre}_whitelist.txt" ]; then
            n_filtered=$(wc -l < "CellBarcodes/comm.${pre}_whitelist.txt")
        fi
        
        # Count cell reads
        if [ -f "$tmp_dir/truecell.${pre}.1.fastq" ]; then
            c=$(wc -l < "$tmp_dir/truecell.${pre}.1.fastq")
            c=$((c/4))
        fi
        
        # Write statistics
        echo "${pre}.labeled: $b" >> "$outfile"
        echo "${pre}.CB: $c" >> "$outfile"
        echo "raw.${pre}.CBtypes: $n_raw" >> "$outfile"
        echo "filtered.${pre}.CBtypes: $n_filtered" >> "$outfile"
        echo "" >> "$outfile"
        
        # Cumulative totals
        local d="$c"  # withCBTotal
        local e="$b"  # withSLTotal
        
        echo "withCBTotal: $d" >> "$outfile"
        echo "withSLTotal: $e" >> "$outfile"
    else
        echo "Warning: files related with ${pre} not found" >> "$outfile"
    fi
}

# Process CASB statistics
process_casb_stats() {
    local chip_name="$1"
    local outfile="$2"
    local tmp_dir="$3"

    echo "=== CASB Statistics ===" >> "$outfile"

    # Get sample names from CASB barcode file
    local samples=()
    if [ -f "$CASB_BARCODES" ]; then
        mapfile -t samples < <(cut -f1 "$CASB_BARCODES")
    fi
    
    # Init cumulative totals
    local d=0  # withCBTotal
    local e=0  # withSLTotal
    
    # Process each sample
    for pre in "${samples[@]}"; do
        local b=0  # labeled reads
        local c=0  # cell reads
        local n_raw=0  # raw barcode types
        local n_filtered=0  # filtered barcode types
        
        # Check file existence
        if [ -f "$tmp_dir/${pre}.R1.fastq" ]; then
            # Count labeled reads
            b=$(wc -l < "$tmp_dir/${pre}.R1.fastq")
            b=$((b/4))
            
            # Count raw whitelist barcode types
            if [ -f "CellBarcodes/${pre}_whitelist.txt" ]; then
                n_raw=$(wc -l < "CellBarcodes/${pre}_whitelist.txt")
            fi
            
            # Count whitelist barcode types after filtering
            if [ -f "CellBarcodes/comm.${pre}_whitelist.txt" ]; then
                n_filtered=$(wc -l < "CellBarcodes/comm.${pre}_whitelist.txt")
            fi
            
            # Count cell reads
            if [ -f "$tmp_dir/truecell.${pre}.1.fastq" ]; then
                c=$(wc -l < "$tmp_dir/truecell.${pre}.1.fastq")
                c=$((c/4))
            fi
            
            # Write statistics
            echo "${pre}.labeled: $b" >> "$outfile"
            echo "${pre}.CB: $c" >> "$outfile"
            echo "raw.${pre}.CBtypes: $n_raw" >> "$outfile"
            echo "${pre}.CBtypes: $n_filtered" >> "$outfile"
            echo "" >> "$outfile"
            
            # Update cumulative totals
            d=$((d + c))
            e=$((e + b))
            
        else
            echo "Warning: files related with  ${pre} not found" >> "$outfile"
        fi
    done
    
    # Write cumulative totals
    echo "withCBTotal: $d" >> "$outfile"
    echo "withSLTotal: $e" >> "$outfile"
}

# R analysis for data integration
run_r_analysis() {
    log_info "Run R analysis for data integration..."
    
    local rds_output="${OUTPUT_DIR}/merged/final_seurat_object.rds"
    
     Rscript "$SCRIPT_DIR/DISC_merge.R" \
            --ANXV_dir "$OUTPUT_DIR/ANXV_lib" \
            --CASB_dir "$OUTPUT_DIR/CASB_lib" \
            --rna_matrix "$RNA_MATRIX" \
            --barcode_translate "$BARCODE_TRANSLATE" \
            --output_dir "$OUTPUT_DIR/merged" \
            --quantile "$HTODEMUX_QUANTILE" \
            --output_rds "$rds_output" \
            2>&1 | tee "$LOG_DIR/R_analysis.log"
        
    # Check if RDS file is generated
    if [ -f "$rds_output" ]; then
            log_info "RDS file is generated: $rds_output"
        else
            log_error "RDS file generation failed!"
            return 1
        fi           
    
    return 0
}

# Generate final report
generate_final_report() {
    log_info "Generate final report..."
    
    local report_file="$OUTPUT_DIR/merged/analysis_report.md"
    
    {
        cat > "$report_file" << EOF
# DISC-Seq RNA-seq Data Demultiplexing and Integration Report

## Overview
- Generate time: $(date)
- Output dir: $OUTPUT_DIR
- HTOdemux quantile: $HTODEMUX_QUANTILE
- Threads number: $THREADS

## Input Files
1. ANXV library dir: $ANXV_LIB_DIR
2. CASB library dir: $CASB_LIB_DIR
3. RNA matrix dir: $RNA_MATRIX
4. Cell barcode whitelist file: $CB_WHITELIST
5. Barcode translation file: $BARCODE_TRANSLATE
6. ANXV barcode file: $ANXV_BARCODES
7. CASB barcode file: $CASB_BARCODES

## Output Structure
\`\`\`
$OUTPUT_DIR/
├── ANXV_lib/          # ANXV library processing results
│   ├── chipname1/     
│   │   ├── CellBarcodes/
│   │   │   ├── comm.ANX5-labeled_whitelist.txt
│   │   │   ├── ANX5-labeled_whitelist.txt
│   │   │   └── ANX5-labeled_expect_whitelist_cell_barcode_count_cutoff.png
│   │   ├── chipname1_ANXV_stats.txt   
│   └── ...           
├── CASB_lib/         # CASB library processing results
│   ├── chipname1/    
│   │   ├── CellBarcodes/
│   │   │   ├── comm.sample1_whitelist.txt
│   │   │   ├── sample1_whitelist.txt
│   │   │   ├── sample1_expect_whitelist_cell_barcode_count_cutoff.png
│   │   │   └── ...  
│   │   ├── chipname1_CASB_stats.txt
│   └── ...
├── merged/            # Merged analysis results
│   ├── final_seurat_object.rds
│   ├── cell_metadata.csv
│   ├── cell_statistics.csv
│   ├── hto_statistics.csv
│   ├── sample_statistics.csv
│   ├── qc_plots/
│   ├── analysis_report.md
│   ├── analysis_parameters.txt
│   └── r_analysis_report.txt
├── logs/              # Log files
└── qc_reports/        # QC reports
\`\`\`

## Analysis Steps
1. Quality control (fastp)
2. Barcode split (fastq-multx)
3. Cell barcode extraction (whitelist.py)
4. Data merging (R/Seurat)
5. Sample split (HTODemux)
6. Statistics generation

## Seurat Object Metadata Description
- HTO_classification: Cell classification (Singlet/Doublet/Negative)
- HTO_classification.global: Cell classification across samples
- HTO_maxID: Maximum HTO identity
- HTO_maxCount: Maximum HTO count
- HTO_secondID: Second highest HTO identity
- HTO_secondCount: Second highest HTO count
- HTO_ratio: Ratio of max to second HTO counts
- nCount_HTO: HTO UMI count
- nFeature_HTO: HTO feature count
- ANXV: ANX5 label
- log2ANXV: log2 transformed ANX5 count
- log10ANXV: log10 transformed ANX5 count
- percent.mt: percentage of mitochondrial genes
- nCount_RNA: UMI count
- nFeature_RNA: Feature count

## Manual for Loading Seurat Object
Load the Seurat object in R using the following command:
\`\`\`r
library(Seurat)
seurat_obj <- readRDS("$OUTPUT_DIR/merged/final_seurat_object.rds")
\`\`\`

## Log Information
Detailed log in: $LOG_DIR/
EOF
    }
    
    log_info "Log report in: $report_file"
}

# Main function
main() {
    log_info "DISC-Seq Data Demultiplexing and Integration Pipeline"
    log_info "========================================="
    
    # Check required tools
    check_tools || exit 1
    
    # Validate inputs
    validate_inputs || exit 1
    
    # Create output directories
    create_directories
    
    # Record configuration
    log_configuration "$@"
    
    # Get FASTQ files
    mapfile -t ANXV_files < <(find "$ANXV_LIB_DIR" -name "*R1*.fastq*" -o -name "*R1*.fq*" 2>/dev/null | sort)
    mapfile -t CASB_files < <(find "$CASB_LIB_DIR" -name "*R1*.fastq*" -o -name "*R1*.fq*" 2>/dev/null | sort)
    
    if [ ${#ANXV_files[@]} -eq 0 ] || [ ${#CASB_files[@]} -eq 0 ]; then
        log_error "Cannot find FASTQ files for processing!"
        exit 1
    fi
    
    log_info "Found ${#ANXV_files[@]} ANXV library samples"
    log_info "Found ${#CASB_files[@]} CASB library samples"
    
    # Handle each ANXV library
    for ((i=0; i<${#ANXV_files[@]}; i++)); do
        ANXV_fq1="${ANXV_files[$i]}"
        ANXV_fq2="${ANXV_fq1/R1/R2}"
        
        # Extract sample name
        local base_name=$(basename "$ANXV_fq1")
    
        # Delete file extensions
        local chip_name=${base_name%.fastq.gz}
        chip_name=${chip_name%.fq.gz}
        chip_name=${chip_name%.fastq}
        chip_name=${chip_name%.fq}
    
        # Remove trailing _R1 or .R1 or similar patterns
        chip_name=$(echo "$chip_name" | sed -E '
            s/_R[12](_[0-9]+)?$//;  # Remove _R1, _R2, _R1_1, etc.
            s/\.[Rr][12]$//;         # Remove .R1 or .r2 at the end
        ')
    
        # If there is a trailing _number, remove it
        chip_name=$(echo "$chip_name" | sed 's/_[0-9]\+$//')
            
        log_info "ANXV library: $chip_name"
        
        # Process ANXV library
        process_ANXV_library "$chip_name" "$ANXV_fq1" "$ANXV_fq2"
    done
    
    for ((i=0; i<${#CASB_files[@]}; i++)); do
        CASB_fq1="${CASB_files[$i]}"
        CASB_fq2="${CASB_fq1/R1/R2}"
        
        # Extract sample name
        local base_name=$(basename "$CASB_fq1")
    
        # Delete file extensions
        local chip_name=${base_name%.fastq.gz}
        chip_name=${chip_name%.fq.gz}
        chip_name=${chip_name%.fastq}
        chip_name=${chip_name%.fq}
    
        # Remove trailing _R1 or .R1 or similar patterns
        chip_name=$(echo "$chip_name" | sed -E '
            s/_R[12](_[0-9]+)?$//;  # Remove _R1, _R2, _R1_1, etc.
            s/\.[Rr][12]$//;         # Remove .R1 or .r2 at the end
        ')
    
        # If there is a trailing _number, remove it
        chip_name=$(echo "$chip_name" | sed 's/_[0-9]\+$//')  

        log_info "CASB library: $chip_name"
        
        # Process CASB library
        process_CASB_library "$chip_name" "$CASB_fq1" "$CASB_fq2"
    done
    
    # R analysis for data integration
    run_r_analysis || {
        log_error "R analysis failed!"
        exit 1
    }

    # Cleanup remaining temporary directories
    cleanup_remaining_tmp() {
    log_info "Cleaning up remaining temporary directories..."
    find "$OUTPUT_DIR" -name "tmp.*" -type d 2>/dev/null | while read dir; do
        if [ -d "$dir" ]; then
            log_info "Cleaning: $dir"
            rm -rf "$dir" 2>/dev/null || true
        fi
    done
    }

    cleanup_remaining_tmp

    # Generate final report
    generate_final_report
    
    log_info "========================================="
    log_info "Completed all processes!"
    log_info "Output dir: $OUTPUT_DIR"
    log_info "RDS file: $OUTPUT_DIR/merged/final_seurat_object.rds"
    log_info "Log dir: $LOG_DIR"
    
    # Preview RDS metadata
    if [ -f "$OUTPUT_DIR/merged/final_seurat_object.rds" ]; then
        log_info "Preview RDS metadata..."
        Rscript -e "
            library(Seurat)
            obj <- readRDS('$OUTPUT_DIR/merged/final_seurat_object.rds')
            cat('\nSeurat object metadata preview:\n')
            print(head(obj@meta.data))
            cat('\nCell number:', ncol(obj), '\n')
            cat('Gene number:', nrow(obj), '\n')
            cat('Sample distribution:\n')
            print(table(obj\$HTO_classification))
        " 2>/dev/null || log_warn "Cannot read RDS file to preview metadata."
    fi
}

# Run main function
main "$@"