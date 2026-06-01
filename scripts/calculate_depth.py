#!/usr/bin/env python3
import os
import glob
import argparse
import subprocess

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate locus-specific read depth across BED coordinates.")
    
    # Define customizable command-line interfaces
    parser.add_argument("-i", "--input-dir", required=True, 
                        help="Directory containing coordinate-sorted duplicate-marked BAM files (Step 5 output)")
    parser.add_argument("-b", "--bed-file", required=True, 
                        help="Path to the target BED file specifying feature genomic intervals")
    parser.add_argument("-o", "--output-dir", required=True, 
                        help="Landing directory to save generated per-base .depth matrices")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Create landing directory safely
    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve sample IDs directly from the target input directory
    input_files = glob.glob(os.path.join(args.input_dir, "*.marked.bam"))
    if not input_files:
        print(f"⚠️ Warning: No valid '.marked.bam' allocation patterns found in {args.input_dir}")
        return

    sample_ids = [os.path.basename(f).replace(".marked.bam", "") for f in input_files]

    # Dynamically look for existing profiles in the custom output folder to prevent redundant work
    search_path = os.path.join(args.output_dir, "*.depth")
    existing_depth_files = {os.path.basename(f).replace('.depth', '') for f in glob.glob(search_path)}

    for sample_id in sorted(sample_ids):
        if sample_id in existing_depth_files:
            print(f"⏭️ Skipping {sample_id} — depth matrix already exists")
            continue

        bam_path = os.path.join(args.input_dir, f"{sample_id}.marked.bam")
        output_file = os.path.join(args.output_dir, f"{sample_id}.depth")

        print(f"📊 Extrapolating base depth for sample: {sample_id}...")
        
        # Production-grade process invocation via subprocess stream pipes
        cmd = f"samtools depth -b {args.bed_file} {bam_path} > {output_file}"
        result = subprocess.run(cmd, shell=True)

        if result.returncode != 0:
            print(f"❌ Process Error: Failure executing sample {sample_id} (Code: {result.returncode})")
        else:
            print(f"✅ Depth calculation complete: {sample_id}")

    print("\n🎉 Bulk depth calculation sequence complete.")

if __name__ == "__main__":
    main()