#!/bin/bash
#SBATCH --job-name=variant_array
#SBATCH --output=logs/variant_array_%A_%a.out
#SBATCH --error=logs/variant_array_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8        
#SBATCH --mem=16G                
#SBATCH --time=12:00:00
#SBATCH --array=1-50             

# ==============================================================================
# 1. ENVIRONMENT SETTINGS
# ==============================================================================
module load StdEnv/2023 gcc/12.3 bwa/0.7.17 samtools/1.19 gatk/4.5.0.0 seqkit/2.8.0

INPUT_RAW_DIR="/path/to/your/raw_fastq_directory"
REF_FASTA="/path/to/your/reference/C_parapsilosis_reference.fasta"
SNPEFF_JAR="/path/to/your/apps/snpEff/snpEff.jar"
GENES_CSV="./assets/Genes-of-interests.csv"
OUTPUT_SCRATCH="./results"
SCRIPT="./scripts/run_pipeline.py"

mkdir -p logs

# ==============================================================================
# 2. RESOLVE SAMPLE CORRESPONDING TO TASK INDEX
# ==============================================================================
SAMPLES=($(find "$INPUT_RAW_DIR" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sort))
TOTAL_SAMPLES=${#SAMPLES[@]}

if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_SAMPLES" ]; then
    echo "❌ Array Error: Slurm Task index overflows available sample pool."
    exit 1
fi

ACTIVE_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# ==============================================================================
# 3. RUN COHORT VARIANT PIPELINE NODE
# ==============================================================================
python "$SCRIPT" \
  --input-dir "$INPUT_RAW_DIR" \
  --ref-fasta "$REF_FASTA" \
  --snpeff-jar "$SNPEFF_JAR" \
  --genes-csv "$GENES_CSV" \
  --output-dir "$OUTPUT_SCRATCH" \
  --sample-id "$ACTIVE_SAMPLE" \
  --threads 8