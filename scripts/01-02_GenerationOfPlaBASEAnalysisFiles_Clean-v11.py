import os
import subprocess
import pandas as pd
from Bio import SeqIO

# --- CONFIGURATION ---
INPUT_DIR = "PLABASE_gbff_Input"  # Directory containing .gbff files
FASTA_OUTPUT_DIR = "DIAMOND_FASTA_Input"  # Central FASTA storage
DIAMOND_DB = "~/mount_point/PlaBase/PGPT_BASE_nr_Aug2021n_ul_1.dmnd"  # DIAMOND database
DIAMOND_OUTPUT_DIR = "diamond_results"  # DIAMOND output directory
RENAMED_OUTPUT_DIR = "renamed_diamond_results"  # Final renamed files
PROCESSED_DIR = "processed_diamond_results"  # Processed DIAMOND results
FILTERED_DIR = "filtered_diamond_results"  # Filtered results
FINAL_DIR = "final_processed_results"  # Final results
MERGED_DIR = "merged_results"  # Merged results
REFINED_DIR = "refined_results"  # Refined results with only one best hit per PGPT family
ONTOLOGY_FILE = "PLABASE-Ontology-finaln6911.txt"  # Ontology file
LOG_FILE = "extraction.log"  # Log file for CDS extraction

# Prompt for DIAMOND sensitivity mode
sensitivity_mode = input("Enter DIAMOND sensitivity mode (faster/fast/mid-sensitive/sensitive/more-sensitive/very-sensitive/ultra-sensitive) [default]:").strip()
sensitivity_flag = f"--{sensitivity_mode}" if sensitivity_mode in ["faster", "fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"] else ""

diamond_cmd_template = (
    "diamond blastp -q {input_fasta} -d {db} -o {output_file} "
    "--outfmt 6 qtitle qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
    "--evalue 1e-5 --threads 8 {sensitivity_flag}"
)

# --- STEP 1: Extract Sequences from GBFF ---
def extract_sequences(skip_fasta=False):
    global FASTA_OUTPUT_DIR
    if skip_fasta:
        FASTA_OUTPUT_DIR = input("Enter existing FASTA input directory: ").strip()
        if not os.path.exists(FASTA_OUTPUT_DIR):
            print(f" Error: Directory '{FASTA_OUTPUT_DIR}' does not exist. Exiting.")
            exit(1)
        print(f" Using existing FASTA directory: {FASTA_OUTPUT_DIR}")
        return
        print(" Skipping FASTA extraction step...")
        return

    os.makedirs(FASTA_OUTPUT_DIR, exist_ok=True)
    genome_count = 0

    log = open(LOG_FILE, "w")  # Open log file outside the loop to keep it open
    try:
        for file in os.listdir(INPUT_DIR):
            if file.endswith(".gbff"):
                genome_count += 1
                cds_count = 0  # Track the number of CDSs extracted
                print(f" Processing genome {genome_count}: {file}")

                fasta_path = os.path.join(FASTA_OUTPUT_DIR, f"{file}.fasta")
                with open(fasta_path, "w") as fasta_out:
                    for record in SeqIO.parse(os.path.join(INPUT_DIR, file), "genbank"):
                        for feature in record.features:
                            if feature.type == "CDS":
                                cds_count += 1
                                locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                                product = feature.qualifiers.get("product", ["unknown"])[0]
                                sequence = feature.qualifiers.get("translation", [""])[0]
                                fasta_out.write(f">{record.id}_{locus_tag}_{gene_name} | Source: {file} | Product: {product}\n{sequence}\n")

# Log file name and number of extracted CDSs
                log.write(f"{file}: {cds_count} CDSs extracted\n")
                log.flush()  # Ensure data is written to the log file

    finally:
        log.close()

# --- STEP 2: Run DIAMOND ---
def run_diamond():
    os.makedirs(DIAMOND_OUTPUT_DIR, exist_ok=True)
    os.makedirs(RENAMED_OUTPUT_DIR, exist_ok=True)
    fasta_files = [f for f in os.listdir(FASTA_OUTPUT_DIR) if f.endswith(".fasta")]
    total_files = len(fasta_files)
    
    for i, fasta_file in enumerate(fasta_files, start=1):
        print(f" Running DIAMOND ({i}/{total_files}): {fasta_file}")
        input_fasta = os.path.join(FASTA_OUTPUT_DIR, fasta_file)
        output_file = os.path.join(DIAMOND_OUTPUT_DIR, fasta_file.replace(".fasta", ".tsv"))
        cmd = diamond_cmd_template.format(input_fasta=input_fasta, db=DIAMOND_DB, output_file=output_file, sensitivity_flag=sensitivity_flag)
        subprocess.run(cmd, shell=True, check=True)
        os.rename(output_file, os.path.join(RENAMED_OUTPUT_DIR, fasta_file.replace(".fasta", "_diamond_output.tsv")))

# --- STEP 3: Process DIAMOND Results ---
def process_results():
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    tsv_files = [f for f in os.listdir(RENAMED_OUTPUT_DIR) if f.endswith(".tsv")]
    for i, file in enumerate(tsv_files, start=1):
        print(f" Processing DIAMOND results ({i}/{len(tsv_files)}): {file}")
        df = pd.read_csv(os.path.join(RENAMED_OUTPUT_DIR, file), sep="\t", header=None, dtype=str)
        df_split = df[0].str.extract(r"^([^_]+)_([^_]+)_(.+) \| (.+) \| (.+)$")
        df_split.columns = ["Genome_Accession", "CDS", "Gene_Name", "Source", "Product"]
        df_split["Source"] = df_split["Source"].str.replace(r"^Source: ", "", regex=True)
        df_split["Source"] = df_split["Source"].str.replace(r".gbff", "", regex=True)
        df_final = pd.concat([df_split, df.iloc[:, 1:]], axis=1)
        df_final.to_csv(os.path.join(PROCESSED_DIR, file), sep="\t", index=False, header=False)

# --- STEP 4: Filter DIAMOND Results ---
def filter_results():
    os.makedirs(FILTERED_DIR, exist_ok=True)

    def parse_e_value(e_value):
        """Convert e-value to float safely, handling potential errors."""
        try:
            return float(e_value)
        except ValueError:
            return float('inf')  # Assign a high value if parsing fails

    for file in os.listdir(PROCESSED_DIR):
        if file.endswith(".tsv"):
            file_path = os.path.join(PROCESSED_DIR, file)
            df = pd.read_csv(file_path, sep="\t", header=None)

            # Define column names
            df.columns = [
                "Accession", "Locus", "Tag", "Source", "Product", "CDS", "PGPT_Hit", "Identity",
                "Alignment_Length", "Mismatches", "Gap_Openings", "Query_Start", "Query_End",
                "Subject_Start", "Subject_End", "E_Value", "Bit_Score"
            ]

            # Extract PGPT family (before first underscore)
            df["PGPT_Family"] = df["PGPT_Hit"].astype(str).apply(lambda x: x.split("_")[0] if "_" in x else x)

            # Convert e-values safely
            df["E_Value_Float"] = df["E_Value"].apply(parse_e_value)

            #  Remove problematic rows where E-value conversion failed
            df = df[df["E_Value_Float"] != float('inf')]

            #  Primary filter: Keep only hits with E-value ≤ 1E-50
            df = df[df["E_Value_Float"] <= 1E-50]

            #  Secondary filter: Keep **best E-value per PGPT family within each CDS**
            df = df.loc[df.groupby(["CDS", "PGPT_Family"])["E_Value_Float"].idxmin()]

            #  Final filter: Keep **best overall hit per CDS**
            df = df.loc[df.groupby("CDS")["E_Value_Float"].idxmin()]

            # Remove helper columns
            df = df.drop(columns=["PGPT_Family", "E_Value_Float"])

            # Save filtered results
            df.to_csv(os.path.join(FILTERED_DIR, file), sep="\t", index=False)

# --- STEP 5: Process Final Results ---
def process_final_results():
    os.makedirs(FINAL_DIR, exist_ok=True)
    tsv_files = [f for f in os.listdir(FILTERED_DIR) if f.endswith(".tsv")]
    for i, file in enumerate(tsv_files, start=1):
        print(f" Final processing ({i}/{len(tsv_files)}): {file}")
        df = pd.read_csv(os.path.join(FILTERED_DIR, file), sep="\t")
        df[["PGPT_Family", "PGPT_Number"]] = df["PGPT_Hit"].str.split("_", expand=True)
        df.to_csv(os.path.join(FINAL_DIR, file), sep="\t", index=False)


# --- STEP 6: Merge Results with Ontology ---
def merge_with_ontology():
    os.makedirs(MERGED_DIR, exist_ok=True)
    ontology_df = pd.read_csv(ONTOLOGY_FILE, sep="\t", dtype=str)

    summary_path = os.path.join(MERGED_DIR, "PGPT_summary.txt")

    with open(summary_path, "a") as summary_file:  # Append mode

        for file in os.listdir(FINAL_DIR):
            if file.endswith(".tsv"):
                df = pd.read_csv(os.path.join(FINAL_DIR, file), sep="\t", dtype=str)
                df = df.merge(ontology_df[["ID", "PGPT_ID", "PATHS"]], left_on="PGPT_Family", right_on="ID", how="left").drop(columns=["ID"]) # Merge with ontology
                df.drop_duplicates(subset=["CDS", "PGPT_Family"], keep="first", inplace=True)  # Remove duplicate PGPT families
                unique_pgpt_count = df["PGPT_Family"].nunique()  # Get unique PGPT families count
                df.to_csv(os.path.join(MERGED_DIR, file), sep="\t", index=False)
                summary_file.write(f"{file}: {unique_pgpt_count} unique PGPT families\n")

# --- STEP 7: Refine PGPT Hits ---
def refine_pgpt_hits():
    for file in os.listdir(MERGED_DIR):
        if file.endswith(".tsv"):
            file_path = os.path.join(MERGED_DIR, file)
            df = pd.read_csv(file_path, sep="\t")
            
            # Keeping all "best" CDS hits per PGPT_Family
            df.to_csv(os.path.join(MERGED_DIR, file), sep="\t", index=False)
            
            # Keeping only the best overall hit per PGPT_Family
            best_hits = df.sort_values(by="Bit_Score", ascending=False).drop_duplicates(subset=["PGPT_Family"])
            best_hits.to_csv(os.path.join(MERGED_DIR, f"best_{file}"), sep="\t", index=False)
    print(" PGPT hits refined and saved.")

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    skip_fasta = input("Skip FASTA extraction? (yes/no): ").strip().lower() == "yes"
    output_folder = input("Enter output folder name: ").strip()
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
        FASTA_OUTPUT_DIR = os.path.join(output_folder, FASTA_OUTPUT_DIR)
        DIAMOND_OUTPUT_DIR = os.path.join(output_folder, DIAMOND_OUTPUT_DIR)
        RENAMED_OUTPUT_DIR = os.path.join(output_folder, RENAMED_OUTPUT_DIR)
        PROCESSED_DIR = os.path.join(output_folder, PROCESSED_DIR)
        FILTERED_DIR = os.path.join(output_folder, FILTERED_DIR)
        FINAL_DIR = os.path.join(output_folder, FINAL_DIR)
        MERGED_DIR = os.path.join(output_folder, MERGED_DIR)
        REFINED_DIR = os.path.join(output_folder, REFINED_DIR)
        LOG_FILE = os.path.join(output_folder, LOG_FILE)
        os.makedirs(FASTA_OUTPUT_DIR, exist_ok=True)
        os.makedirs(DIAMOND_OUTPUT_DIR, exist_ok=True)
        os.makedirs(RENAMED_OUTPUT_DIR, exist_ok=True)
        os.makedirs(PROCESSED_DIR, exist_ok=True)
        os.makedirs(FILTERED_DIR, exist_ok=True)
        os.makedirs(FINAL_DIR, exist_ok=True)
        os.makedirs(MERGED_DIR, exist_ok=True)
    
    extract_sequences(skip_fasta)
    run_diamond()
    process_results()
    filter_results()
    process_final_results()
    merge_with_ontology()
    refine_pgpt_hits()
    print(" Processing complete!") 
