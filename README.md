# Candida parapsilosis Genomics Pipeline

An automated, high-throughput bioinformatics suite designed to analyze next-generation sequencing (NGS) data from *Candida parapsilosis*. This pipeline combines end-to-end variant calling (SNVs/Indels) with parallelized Copy Number Variation (CNV) detection to profile genomic mechanisms of antifungal resistance.

---

## 📌 Introduction & Background

*Candida parapsilosis* is a major opportunistic fungal pathogen, with drug-resistant clinical isolates rising globally. A primary driver of resistance to azole antifungals (like fluconazole) is mutations and copy number expansions in the targeted *ERG11* gene. 

This repository provides a reliable, reproducible framework to process raw sequencing reads into clean, actionable variant matrices and copy number profiles. This exact framework was utilized to process clinical datasets and identify genomic resistance mechanisms for our peer-reviewed research:

> **Zhang Z**, Wang Y, et al. *Persistence and spread of fluconazole-resistant Candida parapsilosis clinical isolates associated with increased ERG11 copies in Qatar.* **Microbial Genomics**, 2026;12(2):001653.

---

## 🧬 Pipeline Architecture

To maximize compute efficiency on High-Performance Computing (HPC) clusters, the suite is divided into two decoupled modules. This modular design prevents large memory bottlenecks and isolates processing steps.

1. **Module 1: Short-Variant Pipeline (`run_variant_pipeline.py`)**
   A standalone Python utility that automates secondary NGS analysis. It streams data sequentially per-sample from raw FASTQ reads all the way to functional variant annotation and candidate gene filtering.
   
2. **Module 2: Scalable CNV Suite (`run_cnv_array.sh`)**
   A decoupled module optimized for Slurm workload managers. Because mapping per-base coverage over large cohorts is computationally intensive, this step utilizes Slurm job arrays to run deep-depth profiling in parallel across dozens of samples concurrently.

---

## 🛠️ Prerequisites

### System Requirements
The wrappers assume a standard Linux environment with the following binaries available in your `$PATH` or loaded via environment modules:
* **Alignment & QC:** BWA, Samtools, SeqKit
* **Trimming & Core GATK:** Trimmomatic, GATK4
* **Functional Annotation:** snpEff (configured with the `Candida_parapsilosis_cdc317` database)

### Package Dependencies
* **Python:** `pandas`, `numpy`, `PyVCF3`
* **R:** `dplyr`, `data.table`, `rtracklayer`

---

## 🚀 Quick Start Guide

### 1. Installation
Clone this repository and install the Python dependencies:
```bash
git clone https://github.com/june-zhang-bioinfo/candida-genomics-pipeline.git
cd candida-genomics-pipeline
pip install -r requirements.txt
```

### 2. Run Variant Calling (Module 1)
Point the master Python script to your raw data directory. It will automatically resolve fastq structures, handle multi-lane merging, and execute alignment through variant annotation:
```bash
python run_variant_pipeline.py   -i /path/to/raw_fastq_folders   -r /path/to/C_parapsilosis_reference.fasta   --snpeff-jar /path/to/snpEff.jar   --genes-csv assets/Genes-of-interests-251115.csv   -o ./results   --threads 8
```

### 3. Run Parallel CNV Profiling (Module 2)
Once Module 1 completes, use the coordinate-sorted, deduplicated BAM files to calculate regional read depth and estimate copy states:

```bash
# Step A: Compute raw targeted depth metrics across genomic intervals
python scripts/1_depth.py

# Step B: Submit the parallel R normalization engine to the Slurm cluster
sbatch run_cnv_array.sh
```

---

## 📁 Data Provenance & Output Structure

All execution outputs land cleanly in segregated, intuitive step-directories under your designated output folder:

```text
results/
├── 1_trimming/            # Quality-filtered and adapter-clipped reads
├── 2_alignment/           # Raw sequence alignment streams (.sam)
├── 3_conversion/          # Compressed binary format tracks (.bam)
├── 4_sorting/             # Coordinate-sorted alignment maps (.sorted.bam)
├── 5_markdup/             # Deduplicated tracking metrics and indexed BAM files (.marked.bam)
├── 6_variant_calling/     # Raw multi-ploidy germline variant sheets (.vcf)
├── 7_variant_filtering/   # GATK hard-filtered variant call sets (.filtered.vcf)
├── 8_update_chr/          # Chromosome accessions regularized to clean numeric strings
├── 9.1_annotation/        # HGVS-compliant variant effects and HTML report summaries
├── 9.2_annotation/        # Legacy EFF format multi-effect classifications
├── 10_vcf_extraction/     # Final high-confidence (QUAL >= 500) parsed variant matrices (.csv)
└── 11_cnv/                # Locus-specific normalized copy-number state summaries (.csv)
```

---

## 🔮 Future Roadmap
* **Workflow Orchestration:** Port the current infrastructure to native Nextflow/Snakemake DSL formats for superior step-caching and container portability.
* **Vector Optimization:** Replace the row-wise matching loops in the legacy CNV script with vectorized `GenomicRanges` structures or fast C++ implementations using `samtools bedcov`.
