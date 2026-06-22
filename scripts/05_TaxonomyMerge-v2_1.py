import os
import re
import argparse
import pandas as pd


# =========================
# Helper functions
# =========================

def clean_taxonomy_value(value):
    """
    Remove GTDB prefixes such as d__, p__, c__, o__, f__, g__, s__.
    """
    if pd.isna(value):
        return "Unknown"

    value = str(value).strip()

    if value == "":
        return "Unknown"

    value = re.sub(r"^[a-z]__", "", value)

    return value if value else "Unknown"


def normalize_id(value):
    """
    Normalize genome/sample identifiers for fuzzy matching.
    Removes punctuation and makes lowercase.

    Example:
    P.alni_DSM44104_PHUJ01 -> palnidsm44104phuj01
    PAlni_DSM44104-PHUJ01.1 -> palnidsm44104phuj011
    """
    value = str(value).lower()
    value = re.sub(r"[^a-z0-9]", "", value)
    return value


def extract_possible_ids(source_value):
    """
    Extract possible genome/sample identifiers from Source values.

    Supports:
    - GCA workflow:
      GCA_002747255.1_ASM274725v1_genomic.gbff

    - FASTA/custom workflow:
      HAG010016_QUXAI209AT
      FH1421_QUXAI207AD
      P.alni_DSM44104_PHUJ01

    Returns possible IDs from most specific to least specific.
    """

    source_value = str(source_value).strip()
    possible_ids = []

    # Case 1: GCA accession
    gca_match = re.search(r"GCA_\d+\.\d+", source_value)
    if gca_match:
        possible_ids.append(gca_match.group())

    # Case 2: remove common wrappers/prefixes/suffixes
    cleaned = source_value

    cleaned = cleaned.replace("best_", "")
    cleaned = cleaned.replace("_diamond_output.tsv", "")
    cleaned = cleaned.replace("_diamond_output", "")
    cleaned = cleaned.replace(".gbff", "")
    cleaned = cleaned.replace(".gb", "")
    cleaned = cleaned.replace(".faa", "")
    cleaned = cleaned.replace(".fasta", "")
    cleaned = cleaned.strip()

    if cleaned and cleaned not in possible_ids:
        possible_ids.append(cleaned)

    # Case 3: shortened ID before underscore
    # Example: HAG010016_QUXAI209AT -> HAG010016
    if "_" in cleaned:
        short_id = cleaned.split("_")[0]
        if short_id and short_id not in possible_ids:
            possible_ids.append(short_id)

    # Remove duplicates while preserving order
    unique_ids = []
    for item in possible_ids:
        if item not in unique_ids:
            unique_ids.append(item)

    return unique_ids


def find_taxonomy(possible_ids, taxonomy_df):
    """
    Match possible IDs against taxonomy table.

    Matching order:
    1. Exact match
    2. Substring match
    3. Normalized fuzzy match
    """

    # Exact match
    for query in possible_ids:
        query = str(query).strip()

        hit = taxonomy_df[taxonomy_df["user_genome"] == query]

        if not hit.empty:
            return hit.iloc[0], query, "exact"

    # Substring match
    for query in possible_ids:
        query = str(query).strip()

        hit = taxonomy_df[
            taxonomy_df["user_genome"].str.contains(
                query,
                regex=False,
                na=False
            )
        ]

        if not hit.empty:
            return hit.iloc[0], query, "substring"

    # Normalized fuzzy match
    for query in possible_ids:
        query_norm = normalize_id(query)

        if query_norm == "":
            continue

        hit = taxonomy_df[
            taxonomy_df["user_genome_normalized"].apply(
                lambda x: query_norm in x or x in query_norm
            )
        ]

        if not hit.empty:
            return hit.iloc[0], query, "normalized"

    return None, " / ".join(possible_ids), "unmatched"


def prepare_taxonomy_table(taxonomy_file):
    """
    Load and clean taxonomy table.
    """

    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", dtype=str)

    # Remove accidental whitespace in column names, e.g. "Order "
    taxonomy_df.columns = taxonomy_df.columns.str.strip()

    if "user_genome" not in taxonomy_df.columns:
        raise ValueError("Taxonomy file must contain a 'user_genome' column.")

    taxonomy_df["user_genome"] = (
        taxonomy_df["user_genome"]
        .astype(str)
        .str.strip()
    )

    # Treat 'classification' as Domain if Domain does not exist
    if "Domain" not in taxonomy_df.columns:
        if "classification" in taxonomy_df.columns:
            taxonomy_df["Domain"] = taxonomy_df["classification"]
        else:
            taxonomy_df["Domain"] = "Unknown"

    # Ensure all expected taxonomy columns exist
    expected_columns = [
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
    ]

    if "Custom_Group" not in taxonomy_df.columns:
    	taxonomy_df["Custom_Group"] = "NA"
    else:
    	taxonomy_df["Custom_Group"] = taxonomy_df["Custom_Group"].fillna("NA").astype(str).str.strip()

    for col in expected_columns:
        if col not in taxonomy_df.columns:
            taxonomy_df[col] = "Unknown"

        taxonomy_df[col] = taxonomy_df[col].apply(clean_taxonomy_value)

    taxonomy_df["user_genome_normalized"] = taxonomy_df["user_genome"].apply(
        normalize_id
    )

    return taxonomy_df


def insert_taxonomy_columns(df, taxonomy_columns):
    """
    Insert taxonomy columns directly after Source.
    If Product exists, this places taxonomy between Source and Product.
    """

    original_cols = list(df.columns)

    if "Source" not in original_cols:
        return df

    # Remove taxonomy columns from current position if present
    remaining_cols = [
        col for col in original_cols
        if col not in taxonomy_columns
    ]

    source_index = remaining_cols.index("Source")

    new_order = (
        remaining_cols[:source_index + 1]
        + taxonomy_columns
        + remaining_cols[source_index + 1:]
    )

    return df[new_order]


def process_file(input_file, output_file, taxonomy_df, custom_group_label):
    """
    Merge one DIAMOND/PGPT TSV file with taxonomy data.
    """

    df = pd.read_csv(input_file, sep="\t", dtype=str)

    if "Source" not in df.columns:
        print(f"Skipping {os.path.basename(input_file)}: no Source column")
        return None

    df["Source"] = df["Source"].astype(str).str.strip()

    taxonomy_records = []
    report_records = []

    for source in df["Source"]:

        possible_ids = extract_possible_ids(source)

        tax_row, matched_query, match_type = find_taxonomy(
            possible_ids,
            taxonomy_df
        )

        if tax_row is None:
            record = {
                "Custom_Group": custom_group_label,
                "Domain": "Unknown",
                "Phylum": "Unknown",
                "Class": "Unknown",
                "Order": "Unknown",
                "Family": "Unknown",
                "Genus": "Unknown",
                "Species": "Unknown",
                "Matched_user_genome": "Unknown",
                "Match_Type": "unmatched"
            }
        else:
            record = {
                "Custom_Group": tax_row["Custom_Group"],
                "Domain": tax_row["Domain"],
                "Phylum": tax_row["Phylum"],
                "Class": tax_row["Class"],
                "Order": tax_row["Order"],
                "Family": tax_row["Family"],
                "Genus": tax_row["Genus"],
                "Species": tax_row["Species"],
                "Matched_user_genome": tax_row["user_genome"],
                "Match_Type": match_type
            }

        taxonomy_records.append(record)

        report_records.append({
            "Source": source,
            "Possible_IDs": " / ".join(possible_ids),
            "Matched_Query": matched_query,
            "Matched_user_genome": record["Matched_user_genome"],
            "Match_Type": record["Match_Type"],
            "Custom_Group": record["Custom_Group"],
            "Domain": record["Domain"],
            "Phylum": record["Phylum"],
            "Class": record["Class"],
            "Order": record["Order"],
            "Family": record["Family"],
            "Genus": record["Genus"],
            "Species": record["Species"]
        })

    taxonomy_added = pd.DataFrame(taxonomy_records)

    merged_df = pd.concat(
        [
            df.reset_index(drop=True),
            taxonomy_added.reset_index(drop=True)
        ],
        axis=1
    )

    taxonomy_columns = [
        "Custom_Group",
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Matched_user_genome",
        "Match_Type"
    ]

    merged_df = insert_taxonomy_columns(
        merged_df,
        taxonomy_columns
    )

    merged_df.to_csv(output_file, sep="\t", index=False)

    report_df = pd.DataFrame(report_records)

    return report_df


def process_directory(input_dir, taxonomy_file, output_dir, custom_group_label):
    """
    Process all TSV files in input_dir.
    """

    taxonomy_df = prepare_taxonomy_table(taxonomy_file)

    os.makedirs(output_dir, exist_ok=True)

    files = sorted(
        file for file in os.listdir(input_dir)
        if file.endswith(".tsv")
    )

    print(f"Found {len(files)} TSV files.")

    all_reports = []

    for filename in files:

        input_file = os.path.join(input_dir, filename)
        output_file = os.path.join(output_dir, filename)

        print(f"\nProcessing {filename}")

        report_df = process_file(
            input_file=input_file,
            output_file=output_file,
            taxonomy_df=taxonomy_df,
            custom_group_label=custom_group_label
        )

        if report_df is not None:
            matched = (report_df["Match_Type"] != "unmatched").sum()
            total = len(report_df)

            print(f"  matched rows: {matched}/{total}")
            print(f"  saved to {output_file}")

            report_df.insert(0, "File", filename)
            all_reports.append(report_df)

    if all_reports:
        final_report = pd.concat(all_reports, ignore_index=True)
        report_file = os.path.join(output_dir, "TaxonomyMerge_match_report.tsv")
        final_report.to_csv(report_file, sep="\t", index=False)

        print(f"\nMatch report saved to: {report_file}")

    print("\nDone.")


# =========================
# Main
# =========================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Merge GTDB taxonomy information into DIAMOND/PGPT TSV files. "
            "Supports GCA-based and FASTA/custom genome identifiers."
        )
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing DIAMOND/PGPT .tsv files"
    )

    parser.add_argument(
        "taxonomy_file",
        help="GTDB taxonomy table"
    )

    parser.add_argument(
        "output_dir",
        help="Directory for taxonomy-annotated output files"
    )

    parser.add_argument(
        "--custom_group",
        default="NA",
        help=(
            "User-defined label added as Custom_Group column "
            "(e.g. GlcNAc, Acidobacteria, TestSet)"
        )
    )

    args = parser.parse_args()

    process_directory(
        input_dir=args.input_dir,
        taxonomy_file=args.taxonomy_file,
        output_dir=args.output_dir,
        custom_group_label=args.custom_group
    )
