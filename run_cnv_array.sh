#!/bin/bash
#SBATCH --job-name=cnv_parallel
#SBATCH --output=logs/cnv_parallel_%A_%a.out
#SBATCH --error=logs/cnv_parallel_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --array=1-50

# ==============================================================================
# 1. ENVIRONMENT CONFIGURATION
# ==============================================================================
module load StdEnv/2023 gcc/12.3 r-bundle-bioconductor/3.18
module load r/4.4.0 samtools/1.19

BAM_DIR="./results/5_markdup"
BED_FILE="./assets/C_parapsilosis_CDC317.bed"
CNV_DIR="./results/11_cnv"
GFF_FILE="./assets/C_parapsilosis_CDC317_current_features.gff"
GENES_CSV="./assets/Genes-of-interests.csv"
CNV_RSCRIPT="./scripts/process_cnv.R"

mkdir -p logs $CNV_DIR

# ==============================================================================
# 2. FILE INTERROGATION & ARRAY MAPPING
# ==============================================================================
FILES=("$BAM_DIR"/*.marked.bam)
TOTAL_FILES=${#FILES[@]}

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "❌ No BAM files found in $BAM_DIR"
    exit 1
fi

BAM_FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}
SAMPLE_ID=$(basename "$BAM_FILE" .marked.bam)
DEPTH_FILE="$CNV_DIR/${SAMPLE_ID}.depth"

echo "📊 Processing sample: $SAMPLE_ID"
echo "   BAM: $BAM_FILE"
echo "   Depth output: $DEPTH_FILE"

# ==============================================================================
# 3. STEP 1: Calculate depth (if needed)
# ==============================================================================
if [ ! -f "$DEPTH_FILE" ]; then
    echo "🔢 Calculating read depth..."
    samtools depth -b "$BED_FILE" "$BAM_FILE" > "$DEPTH_FILE"
    
    if [ $? -ne 0 ]; then
        echo "❌ Depth calculation failed for $SAMPLE_ID"
        exit 1
    fi
    echo "✅ Depth calculation complete"
else
    echo "⏭️ Depth file exists, skipping depth calculation"
fi

# ==============================================================================
# 4. STEP 2: Process CNV with R (if needed)
# ==============================================================================
if [ -f "${DEPTH_FILE}.csv" ]; then
    echo "⏭️ Skipping CNV analysis — output already exists"
    exit 0
fi

echo "🧬 Running CNV analysis..."
Rscript "$CNV_RSCRIPT" "$DEPTH_FILE" "$GFF_FILE" "$GENES_CSV"

if [ $? -eq 0 ]; then
    echo "✅ CNV analysis complete: $SAMPLE_ID"
else
    echo "❌ CNV analysis failed: $SAMPLE_ID"
    exit 1
fi