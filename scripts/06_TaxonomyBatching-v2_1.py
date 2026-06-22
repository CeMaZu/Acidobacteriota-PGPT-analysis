import os
import re
import argparse
import pandas as pd


# =========================
# Helper functions
# =========================

def sanitize_filename(value):
    """
    Make safe filenames from taxonomy names.
    """
    value = str(value).strip()
    value = value.replace("/", "_")
    value = value.replace("\\", "_")
    value = re.sub(r"[^\w\-.]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_")


def is_unknown(value):
    """
    Identify missing or unknown taxonomy labels.
    """
    if pd.isna(value):
        return True

    value = str(value).strip()

    if value == "":
        return True

    if value.lower() in ["unknown", "unclassified", "nan", "none", "na"]:
        return True

    return False


def load_tsv_files(input_dir, deduplicated_mode):
    """
    Load TSV files from input directory.

    deduplicated_mode:
    - best: only files starting with best_
    - raw: only files not starting with best_
    - all: all TSV files
    """
    dataframes = []
    loaded_files = []

    files = sorted(
        f for f in os.listdir(input_dir)
        if f.endswith(".tsv")
    )

    for filename in files:

        if deduplicated_mode == "best" and not filename.startswith("best_"):
            continue

        if deduplicated_mode == "raw" and filename.startswith("best_"):
            continue

        filepath = os.path.join(input_dir, filename)

        df = pd.read_csv(filepath, sep="\t", dtype=str)
        df["Input_File"] = filename

        dataframes.append(df)
        loaded_files.append(filename)

    return dataframes, loaded_files


def get_ranks_to_process(full_df, rank_argument):
    """
    Determine which taxonomy columns should be used for batching.
    """
    default_ranks = [
        "Custom_Group",
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
    ]

    available_ranks = [
        rank for rank in default_ranks
        if rank in full_df.columns
    ]

    if rank_argument.lower() == "all":
        return available_ranks

    requested = [
        r.strip()
        for r in rank_argument.split(",")
        if r.strip()
    ]

    missing = [
        r for r in requested
        if r not in full_df.columns
    ]

    if missing:
        raise ValueError(
            f"Requested rank(s) not found in input files: {missing}\n"
            f"Available columns are: {list(full_df.columns)}"
        )

    return requested


def count_genomes(df):
    """
    Count unique genomes in a batch.

    Prefer Matched_user_genome if available,
    otherwise use Source.
    """
    if "Matched_user_genome" in df.columns:
        valid = df["Matched_user_genome"].dropna()
        valid = valid[valid.astype(str).str.lower() != "unknown"]

        if len(valid) > 0:
            return valid.nunique()

    if "Source" in df.columns:
        return df["Source"].nunique()

    return "NA"


def create_batches_for_rank(full_df, rank, output_dir):
    """
    Create batch files for one taxonomic rank.
    """
    rank_output_dir = os.path.join(output_dir, rank)
    os.makedirs(rank_output_dir, exist_ok=True)

    report_dir = os.path.join(output_dir, "Reports")
    os.makedirs(report_dir, exist_ok=True)

    summary_records = []
    unknown_records = []

    groups = full_df[rank].dropna().unique()

    for group in sorted(groups):

        if is_unknown(group):
            continue

        batch_df = full_df[full_df[rank] == group].copy()

        safe_group = sanitize_filename(group)
        batch_filename = f"{rank}_{safe_group}_batch.tsv"
        batch_path = os.path.join(rank_output_dir, batch_filename)

        batch_df.to_csv(batch_path, sep="\t", index=False)

        genome_count = count_genomes(batch_df)

        summary_records.append({
            "Rank": rank,
            "Group": group,
            "Batch_File": batch_filename,
            "Rows": len(batch_df),
            "Genome_Count": genome_count
        })

        print(
            f"  Created {batch_filename}: "
            f"{len(batch_df)} rows, {genome_count} genomes"
        )

    # Unknown / unclassified entries
    unknown_df = full_df[
        full_df[rank].apply(is_unknown)
    ].copy()

    if not unknown_df.empty:

        unknown_filename = f"{rank}_Unknown_unclassified_batch.tsv"
        unknown_path = os.path.join(rank_output_dir, unknown_filename)

        unknown_df.to_csv(unknown_path, sep="\t", index=False)

        genome_count = count_genomes(unknown_df)

        summary_records.append({
            "Rank": rank,
            "Group": "Unknown_or_unclassified",
            "Batch_File": unknown_filename,
            "Rows": len(unknown_df),
            "Genome_Count": genome_count
        })

        print(
            f"  Created {unknown_filename}: "
            f"{len(unknown_df)} rows, {genome_count} genomes"
        )

        # More detailed unknown report
        if "Source" in unknown_df.columns:
            unknown_sources = (
                unknown_df[["Source"]]
                .drop_duplicates()
                .copy()
            )

            if "Matched_user_genome" in unknown_df.columns:
                matched_info = (
                    unknown_df[["Source", "Matched_user_genome", "Match_Type"]]
                    .drop_duplicates()
                    .copy()
                )
                unknown_sources = matched_info

            unknown_sources.insert(0, "Rank", rank)
            unknown_sources.insert(1, "Problem", f"Unknown {rank}")

            unknown_records.append(unknown_sources)

    summary_df = pd.DataFrame(summary_records)

    summary_path = os.path.join(
        report_dir,
        f"{rank}_batch_summary.tsv"
    )

    summary_df.to_csv(summary_path, sep="\t", index=False)

    unknown_report_df = pd.DataFrame()

    if unknown_records:
        unknown_report_df = pd.concat(
            unknown_records,
            ignore_index=True
        )

        unknown_report_path = os.path.join(
            report_dir,
            f"{rank}_unknown_report.tsv"
        )

        unknown_report_df.to_csv(
            unknown_report_path,
            sep="\t",
            index=False
        )

    return summary_df, unknown_report_df


def create_taxonomic_batches(input_dir, output_dir, ranks, deduplicated_mode):
    """
    Main processing function.
    """
    os.makedirs(output_dir, exist_ok=True)

    report_dir = os.path.join(output_dir, "Reports")
    os.makedirs(report_dir, exist_ok=True)

    dataframes, loaded_files = load_tsv_files(
        input_dir,
        deduplicated_mode
    )

    if not dataframes:
        print("No matching TSV files found.")
        return

    print(f"Loaded {len(loaded_files)} files:")
    for f in loaded_files:
        print(f"  {f}")

    full_df = pd.concat(dataframes, ignore_index=True)

    if "Source" not in full_df.columns:
        raise ValueError(
            "Input files must contain a 'Source' column."
        )

    ranks_to_process = get_ranks_to_process(full_df, ranks)

    if not ranks_to_process:
        raise ValueError(
            "No valid taxonomy ranks found. "
            "Make sure TaxonomyMerge-v2 was run before this script."
        )

    print("\nRanks to process:")
    for rank in ranks_to_process:
        print(f"  {rank}")

    all_summaries = []
    all_unknowns = []

    for rank in ranks_to_process:

        print(f"\nProcessing rank: {rank}")

        summary_df, unknown_df = create_batches_for_rank(
            full_df,
            rank,
            output_dir
        )

        if not summary_df.empty:
            all_summaries.append(summary_df)

        if not unknown_df.empty:
            all_unknowns.append(unknown_df)

    # Combined summary report
    if all_summaries:
        combined_summary = pd.concat(
            all_summaries,
            ignore_index=True
        )

        combined_summary_path = os.path.join(
            report_dir,
            "ALL_taxonomic_batch_summary.tsv"
        )

        combined_summary.to_csv(
            combined_summary_path,
            sep="\t",
            index=False
        )

        print(f"\nCombined summary saved to:")
        print(combined_summary_path)

    # Combined unknown report
    if all_unknowns:
        combined_unknown = pd.concat(
            all_unknowns,
            ignore_index=True
        )

        combined_unknown_path = os.path.join(
            report_dir,
            "ALL_unknown_taxonomy_report.tsv"
        )

        combined_unknown.to_csv(
            combined_unknown_path,
            sep="\t",
            index=False
        )

        print(f"\nCombined unknown report saved to:")
        print(combined_unknown_path)

    print("\nDone.")


# =========================
# Main
# =========================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Create taxonomic batch files from taxonomy-annotated PGPT TSV files. "
            "Batches can be generated for all taxonomy ranks or selected ranks."
        )
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing taxonomy-annotated TSV files from TaxonomyMerge-v2"
    )

    parser.add_argument(
        "output_dir",
        help="Directory where taxonomic batch files will be saved"
    )

    parser.add_argument(
        "--ranks",
        default="all",
        help=(
            "Taxonomic ranks to batch by. "
            "Use 'all' for all available ranks, or comma-separated values, "
            "e.g. Family,Custom_Group,Order"
        )
    )

    parser.add_argument(
        "--deduplicated",
        choices=["best", "raw", "all"],
        default="best",
        help=(
            "Which TSV files to process: "
            "'best' = only files starting with best_, "
            "'raw' = only files not starting with best_, "
            "'all' = all TSV files. Default: best"
        )
    )

    args = parser.parse_args()

    create_taxonomic_batches(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        ranks=args.ranks,
        deduplicated_mode=args.deduplicated
    )
