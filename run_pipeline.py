#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse
import vcf
import numpy as np
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="End-to-end variant calling and annotation pipeline for Candida parapsilosis.")
    
    # Required Arguments
    parser.add_argument("-i", "--input-dir", required=True, help="Directory containing raw sample folders (FASTQ files)")
    parser.add_argument("-r", "--ref-fasta", required=True, help="Path to the reference FASTA file")
    parser.add_argument("--snpeff-jar", required=True, help="Path to the snpEff.jar file")
    parser.add_argument("--genes-csv", required=True, help="Path to the Genes of interest CSV file")
    
    # Optional Arguments
    parser.add_argument("-o", "--output-dir", default=os.getcwd(), help="Directory to save all pipeline outputs (default: current directory)")
    parser.add_argument("-e", "--exclude-list", help="Path to a text file containing sample IDs to skip (one per line)")
    parser.add_argument("-t", "--threads", default=8, type=int, help="Number of threads to use for alignment (default: 8)")
    
    # Tool-specific Overrides (Useful if not running on an HPC with Lmod)
    parser.add_argument("--trimmomatic-jar", help="Path to trimmomatic jar (overrides $EBROOTTRIMMOMATIC)")
    parser.add_argument("--adapter-file", help="Path to Illumina adapter FASTA for trimming")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # ==========================================
    # 1. SETUP AND CONFIGURATION
    # ==========================================
    os.makedirs(args.output_dir, exist_ok=True)
    os.chdir(args.output_dir)

    # Load Exclude List
    exclude_samples = set()
    if args.exclude_list and os.path.exists(args.exclude_list):
        with open(args.exclude_list, 'r') as f:
            exclude_samples = {line.strip() for line in f if line.strip()}
        print(f"ℹ️ Loaded {len(exclude_samples)} samples to exclude.")

    # Load Gene List
    if not os.path.exists(args.genes_csv):
        raise FileNotFoundError(f"❌ Cannot find {args.genes_csv}.")
    csv_genes = pd.read_csv(args.genes_csv)
    gene_list = list(csv_genes['Cparap ID'])

    # Determine Trimmomatic paths
    trimmomatic_jar = args.trimmomatic_jar
    adapter_file = args.adapter_file

    if not trimmomatic_jar or not adapter_file:
        trimmomatic_root = os.environ.get("EBROOTTRIMMOMATIC")
        if trimmomatic_root:
            trimmomatic_jar = trimmomatic_jar or os.path.join(trimmomatic_root, "trimmomatic-0.39.jar")
            adapter_file = adapter_file or os.path.join(trimmomatic_root, "adapters", "NexteraPE-PE.fa")
        else:
            raise EnvironmentError("❌ Trimmomatic paths not provided and $EBROOTTRIMMOMATIC is not set. "
                                   "Please load the module or use --trimmomatic-jar and --adapter-file.")

    # Check for seqkit
    try:
        subprocess.run(["seqkit", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        raise EnvironmentError("❌ seqkit is not available in PATH. Please install or load it.")

    # Setup Output Directories
    step_dirs = {
        "trimmed": "1_trimming",
        "aligned": "2_alignment",
        "converted": "3_conversion",
        "sorted": "4_sorting",
        "marked": "5_markdup",
        "called": "6_variant_calling",
        "filtered": "7_variant_filtering",
        "updated_chr": "8_update_chr",
        "annotation_91": "9.1_annotation",
        "annotation_92": "9.2_annotation",
        "extraction": "10_vcf_extraction"
    }

    for folder in step_dirs.values():
        os.makedirs(folder, exist_ok=True)

    # ==========================================
    # 2. SAMPLE PROCESSING LOOP
    # ==========================================
    for sample_folder in sorted(os.listdir(args.input_dir)):
        sample_path = os.path.join(args.input_dir, sample_folder)
        if not os.path.isdir(sample_path):
            continue

        sample_id = sample_folder
        if sample_id in exclude_samples:
            print(f"⏭️ Skipping excluded sample: {sample_id}")
            continue

        try:
            # --- R1 & R2 FastQ Gathering ---
            def sort_by_lane(file_list):
                def extract_lane(filename):
                    for part in filename.split("_"):
                        if part.startswith("L") and part[1:].isdigit():
                            return int(part[1:])
                    return 999
                return sorted(file_list, key=extract_lane)

            r1_list = sort_by_lane(glob.glob(os.path.join(sample_path, "*_1.fq.gz")))
            r2_list = sort_by_lane(glob.glob(os.path.join(sample_path, "*_2.fq.gz")))

            if not r1_list or not r2_list:
                print(f"⚠️ Missing R1 or R2 FASTQ for {sample_id}. Skipping.")
                continue

            print(f"\n🔄 Processing sample: {sample_id}")

            # === Merge & Harmonize FASTQs ===
            if len(r1_list) > 1:
                merged_r1 = os.path.join(step_dirs["trimmed"], f"{sample_id}_merged_R1.fq.gz")
                if not os.path.exists(merged_r1):
                    with open(merged_r1, "wb") as combined_file:
                        for r1_file in r1_list:
                            with open(r1_file, "rb") as f:
                                combined_file.write(f.read())
                harmonized_r1 = os.path.join(step_dirs["trimmed"], f"{sample_id}_harmonized_R1.fq.gz")
                if not os.path.exists(harmonized_r1):
                    subprocess.run(f"seqkit rename {merged_r1} -o {harmonized_r1}", shell=True, check=True)
                r1 = harmonized_r1
            else:
                r1 = r1_list[0]

            if len(r2_list) > 1:
                merged_r2 = os.path.join(step_dirs["trimmed"], f"{sample_id}_merged_R2.fq.gz")
                if not os.path.exists(merged_r2):
                    with open(merged_r2, "wb") as combined_file:
                        for r2_file in r2_list:
                            with open(r2_file, "rb") as f:
                                combined_file.write(f.read())
                harmonized_r2 = os.path.join(step_dirs["trimmed"], f"{sample_id}_harmonized_R2.fq.gz")
                if not os.path.exists(harmonized_r2):
                    subprocess.run(f"seqkit rename {merged_r2} -o {harmonized_r2}", shell=True, check=True)
                r2 = harmonized_r2
            else:
                r2 = r2_list[0]

            # === Step 1: Trimmomatic ===
            trimmed_r1_paired = os.path.join(step_dirs["trimmed"], f"{sample_id}_R1_paired.fq.gz")
            trimmed_r1_unpaired = os.path.join(step_dirs["trimmed"], f"{sample_id}_R1_unpaired.fq.gz")
            trimmed_r2_paired = os.path.join(step_dirs["trimmed"], f"{sample_id}_R2_paired.fq.gz")
            trimmed_r2_unpaired = os.path.join(step_dirs["trimmed"], f"{sample_id}_R2_unpaired.fq.gz")
            if not os.path.exists(trimmed_r1_paired):
                print("🚿 Trimming reads with Trimmomatic...")
                trim_cmd = f"""java -jar {trimmomatic_jar} PE {r1} {r2} \
                {trimmed_r1_paired} {trimmed_r1_unpaired} {trimmed_r2_paired} {trimmed_r2_unpaired} \
                ILLUMINACLIP:{adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"""
                subprocess.run(trim_cmd, shell=True, check=True)

            # === Step 2: BWA MEM Alignment ===
            sam_file = os.path.join(step_dirs["aligned"], f"{sample_id}.sam")
            if not os.path.exists(sam_file):
                print("🧬 Aligning reads with BWA MEM...")
                bwa_cmd = f"""bwa mem -R '@RG\\tID:{sample_id}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{sample_id}' \
                -t {args.threads} {args.ref_fasta} {trimmed_r1_paired} {trimmed_r2_paired} > {sam_file}"""
                subprocess.run(bwa_cmd, shell=True, check=True)

            # === Step 3: Convert SAM to BAM ===
            bam_file = os.path.join(step_dirs["converted"], f"{sample_id}.bam")
            if not os.path.exists(bam_file):
                print("📥 Converting SAM to BAM...")
                subprocess.run(f"samtools view -bS {sam_file} -o {bam_file}", shell=True, check=True)

            # === Step 4: Sort BAM ===
            sorted_bam = os.path.join(step_dirs["sorted"], f"{sample_id}.sorted.bam")
            if not os.path.exists(sorted_bam):
                print("📊 Sorting BAM file...")
                subprocess.run(f"samtools sort -@ {args.threads} {bam_file} -o {sorted_bam}", shell=True, check=True)

            # === Step 5: MarkDuplicates ===
            marked_bam = os.path.join(step_dirs["marked"], f"{sample_id}.marked.bam")
            metrics_file = os.path.join(step_dirs["marked"], f"{sample_id}.metrics.txt")
            if not os.path.exists(marked_bam):
                print("🧪 Marking duplicates with GATK...")
                subprocess.run(f"gatk MarkDuplicates -I {sorted_bam} -O {marked_bam} -M {metrics_file}", shell=True, check=True)
                subprocess.run(f"samtools index {marked_bam}", shell=True, check=True)

            # === Step 6: HaplotypeCaller ===
            vcf_file = os.path.join(step_dirs["called"], f"{sample_id}.vcf")
            if not os.path.exists(vcf_file):
                print("🔍 Calling variants with GATK HaplotypeCaller...")
                subprocess.run(f"gatk HaplotypeCaller -R {args.ref_fasta} -I {marked_bam} -O {vcf_file} -ploidy 2 -mbq 20", shell=True, check=True)

            # === Step 7: Variant Filtration ===
            filtered_vcf = os.path.join(step_dirs["filtered"], f"{sample_id}.filtered.vcf")
            if not os.path.exists(filtered_vcf):
                print("🧹 Filtering variants...")
                filter_cmd = f"""gatk VariantFiltration -R {args.ref_fasta} -V {vcf_file} -O {filtered_vcf} \
                --genotype-filter-expression 'GQ < 20' --genotype-filter-name GQFilter \
                --genotype-filter-expression 'DP < 10' --genotype-filter-name DPFilter"""
                subprocess.run(filter_cmd, shell=True, check=True)

            # === Step 8: Update Chromosomes ===
            updated_vcf = os.path.join(step_dirs["updated_chr"], f"{sample_id}.updated_chr.vcf")
            if not os.path.exists(updated_vcf):
                print("📝 Updating chromosome names...")
                update_cmd = (f"cat {filtered_vcf} | "
                              f"sed 's/^Contig005504_C_parapsilosis_CDC317/5504/' | "
                              f"sed 's/^Contig005569_C_parapsilosis_CDC317/5569/' | "
                              f"sed 's/^Contig005806_C_parapsilosis_CDC317/5806/' | "
                              f"sed 's/^Contig005807_C_parapsilosis_CDC317/5807/' | "
                              f"sed 's/^Contig005809_C_parapsilosis_CDC317/5809/' | "
                              f"sed 's/^Contig006110_C_parapsilosis_CDC317/6110/' | "
                              f"sed 's/^Contig006139_C_parapsilosis_CDC317/6139/' | "
                              f"sed 's/^Contig006372_C_parapsilosis_CDC317/6372/' | "
                              f"sed 's/^mito_C_parapsilosis_CDC317/Mt/' > {updated_vcf}")
                subprocess.run(update_cmd, shell=True, check=True)

            # === Step 9.1 & 9.2: snpEff Annotation ===
            ann_vcf_91 = os.path.join(step_dirs["annotation_91"], f"{sample_id}.ann.vcf")
            ann_html_91 = os.path.join(step_dirs["annotation_91"], f"{sample_id}.html")
            if not os.path.exists(ann_vcf_91):
                print("🧬 Annotating variants (ANN format)...")
                subprocess.run(f"java -Xmx8g -jar {args.snpeff_jar} -v -s {ann_html_91} Candida_parapsilosis_cdc317 {updated_vcf} > {ann_vcf_91}", shell=True, check=True)

            ann_vcf_92 = os.path.join(step_dirs["annotation_92"], f"{sample_id}.ann.vcf")
            ann_html_92 = os.path.join(step_dirs["annotation_92"], f"{sample_id}.html")
            if not os.path.exists(ann_vcf_92):
                print("🧬 Annotating variants (EFF format)...")
                subprocess.run(f"java -Xmx8g -jar {args.snpeff_jar} -formatEff -v -s {ann_html_92} Candida_parapsilosis_cdc317 {updated_vcf} > {ann_vcf_92}", shell=True, check=True)

            # === Step 10: VCF Extraction and Merge ===
            final_csv = os.path.join(step_dirs["extraction"], f"{sample_id}.csv")
            if not os.path.exists(final_csv):
                print("📊 Extracting VCF data and filtering genes...")
                
                # Arrays for 9.1 Data
                EFF, REF, ALT, QUAL, GENE, HET, POS, CHROM = [], [], [], [], [], [], [], []
                vcf_reader1 = vcf.Reader(open(ann_vcf_91, 'r'))
                for record in vcf_reader1:
                    ref_val, alt_val, qual_val = record.REF, record.ALT, record.QUAL
                    pos, chrom = record.POS, record.CHROM
                    het = record.heterozygosity
                    if 'ANN' in record.INFO:
                        for ann in record.INFO['ANN']:
                            parts = ann.split('|')
                            if len(parts) > 3:
                                effect, gene_id = parts[1], parts[3]
                                for gene in gene_list:
                                    if gene_id in gene:
                                        EFF.append(effect)
                                        GENE.append(gene_id)
                                        REF.append(ref_val)
                                        ALT.append(alt_val[0] if alt_val else None)
                                        QUAL.append(qual_val)
                                        HET.append(het)
                                        POS.append(pos)
                                        CHROM.append(chrom)

                # Arrays for 9.2 Data
                AA, GENE2 = [], []
                vcf_reader2 = vcf.Reader(open(ann_vcf_92, 'r'))
                for record in vcf_reader2:
                    if 'EFF' in record.INFO:
                        for eff in record.INFO['EFF']:
                            parts = eff.split('|')
                            if len(parts) > 5:
                                gene_id, aa = parts[5], parts[3]
                                if gene_id:
                                    for gene in gene_list:
                                        if gene_id in gene:
                                            GENE2.append(gene_id)
                                            AA.append(aa)

                # Create Pandas DataFrame
                df = pd.DataFrame(GENE, columns=['Gene'])
                df['Reference'] = REF
                df['Alternative'] = ALT
                df['Effect'] = EFF
                df['Chromosome'] = CHROM
                df['Position'] = POS
                df['Gene_id'] = GENE
                df['Quality Score'] = QUAL
                
                if len(df) == len(AA):
                    df['Amino Acid'] = AA
                else:
                    df['Amino Acid'] = AA[:len(df)] + ["Unknown"] * max(0, len(df) - len(AA))
                    
                df['Heterozygosity'] = HET
                df = df[df['Quality Score'] >= 500]
                
                outcome = df.merge(csv_genes, how='inner', left_on='Gene_id', right_on='Cparap ID')
                outcome.to_csv(final_csv, index=False)
            
            print(f"✅ Finished sample: {sample_id}")

        except subprocess.CalledProcessError as e:
            print(f"❌ Error running shell command for sample {sample_id}:\n{e}\nSkipping to next sample.\n")
        except Exception as e:
            print(f"❌ Unexpected error with sample {sample_id}: {e}\nSkipping to next sample.\n")

    print("\n🎉 All samples processed.")

if __name__ == "__main__":
    main()