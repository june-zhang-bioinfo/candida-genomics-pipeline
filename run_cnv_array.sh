#!/bin/bash
#SBATCH --job-name=cnv_parallel
#SBATCH --output=logs/cnv_parallel_%A_%a.out
#SBATCH --error=logs/cnv_parallel_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00          # Adjusted to a standard benchmark runtime string
#SBATCH --array=1-50             # Adjust upper boundary to match your sample size

# ==============================================================================
# 1. ENVIRONMENT CONFIGURATION
# ==============================================================================
# Load system-level cluster dependencies (Adjust to your local HPC profile)
module load StdEnv/2023 gcc/12.3 r-bundle-bioconductor/3.18
module load r/4.4.0

# Define your file path infrastructures
CNV_DIR="./11_cnv"                     # Location containing generated .depth files
GFF_FILE="./assets/C_parapsilosis_CDC317_current_features.gff"
GENES_CSV="./assets/Genes-of-interests-251115.csv"
CNV_RSCRIPT="./scripts/process_cnv.R"  # Path to the refactored R analysis worker

# Make a dedicated logs folder for cleaner workspace housekeeping
mkdir -p logs

# ==============================================================================
# 2. FILE INTERROGATION & ARRAY MAPPING
# ==============================================================================
# Collect target depth tracking grids into an index array
FILES=("$CNV_DIR"/*.depth)
TOTAL_FILES=${#FILES[@]}

# Safety catch if array indexing overflows or no datasets match boundaries
if [ "$TOTAL_FILES" -eq 0 ] || [ ! -e "${FILES[0]}" ]; then
    echo "❌ Execution Error: No valid '.depth' file matrix targets detected in $CNV_DIR"
    exit 1
fi

# Map current Slurm task thread block to specific target file path
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

echo "📊 HPC Batch Processing Stats:"
echo "   - Collective Matching Metrics Pool Size: $TOTAL_FILES"
echo "   - Current Active Slurm Task ID: $SLURM_ARRAY_TASK_ID"
echo "   - Target Processing Record Locus: $FILE"

# ==============================================================================
# 3. ANALYSIS EXECUTION LOOP
# ==============================================================================
# Skip logic check to protect existing calculations
if [ -f "${FILE}.csv" ]; then
    echo "⏭️ Skipping $FILE — target analysis profile already computed."
    exit 0
fi

# Route arguments downstream directly into our flexible R framework
Rscript "$CNV_RSCRIPT" "$FILE" "$GFF_FILE" "$GENES_CSV"