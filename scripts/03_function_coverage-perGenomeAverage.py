import os
import pandas as pd
from collections import defaultdict

def parse_paths_column(paths_column):
    """Extracts hierarchical PGPT levels from the PATHS column."""
    parsed_data = []
    for entry in paths_column.dropna():
        for path in entry.split('#'):
            levels = path.strip().split(';')[1:]  # Skip root 'PGPT'
            parsed_data.append(levels)
    return parsed_data

def compute_function_coverage_per_genome(batch_df, ontology_df):
    print(" Starting PGPT function coverage computation...")
    
    if batch_df.empty or ontology_df.empty:
        print(" Input data is empty!")
        return pd.DataFrame()

    unique_genomes = batch_df['Genome_ID'].unique()
    genome_count = len(unique_genomes)
    print(f" Found {genome_count} unique genomes.")

    # Precompute total PGPTs at each level from ontology
    total_pgpt_counts = defaultdict(int)
    for path in ontology_df['PATHS'].dropna():
        for subpath in path.split('#'):
            levels = subpath.strip().split(';')[1:]
            for i, level in enumerate(levels):
                total_pgpt_counts[(i, level)] += 1

    print(f" Total PGPTs computed across levels: {len(total_pgpt_counts)}")

    # Collect data per genome
    coverage_records = []

    for genome_id in unique_genomes:
        genome_df = batch_df[batch_df['Genome_ID'] == genome_id]
        parsed_paths = parse_paths_column(genome_df['PATHS'])

        observed_counts = defaultdict(int)
        for levels in parsed_paths:
            for i, level in enumerate(levels):
                observed_counts[(i, level)] += 1

        for (level_idx, level_name), total_count in total_pgpt_counts.items():
            observed = observed_counts.get((level_idx, level_name), 0)
            coverage_pct = (observed / total_count) * 100
            coverage_records.append([
                genome_id, level_idx, level_name,
                observed, total_count, coverage_pct
            ])

    coverage_df = pd.DataFrame(coverage_records, columns=[
        'Genome_ID', 'Level', 'Function', 'Observed_PGPTs',
        'Total_PGPTs', 'Coverage_Percentage'
    ])

    #  FIXED: mean and std are calculated properly
    grouped = coverage_df.groupby(['Level', 'Function']).agg(
        Sum_Observed_PGPTs=('Observed_PGPTs', 'sum'),
        Total_PGPTs=('Total_PGPTs', 'first'),
        Std_Coverage_Percentage=('Coverage_Percentage', 'std')
    ).reset_index()

    grouped['Genome_Count'] = genome_count
    grouped['Mean_Coverage_Percentage'] = (
        grouped['Sum_Observed_PGPTs'] / (grouped['Total_PGPTs'] * genome_count)
    ) * 100

    # Final tidy format
    final_df = grouped[[
        'Genome_Count', 'Level', 'Function',
        'Mean_Coverage_Percentage', 'Std_Coverage_Percentage', 'Sum_Observed_PGPTs'
    ]]

    print(" Final summary:")
    print(final_df.head())
    return final_df

def generate_excel_outputs(coverage_df, output_file):
    """Write coverage data into separate sheets per PGPT level."""
    with pd.ExcelWriter(output_file) as writer:
        for level_idx in sorted(coverage_df['Level'].unique()):
            subset = coverage_df[coverage_df['Level'] == level_idx]
            subset.to_excel(writer, sheet_name=f'Level_{level_idx}', index=False)
    print(f" Excel file saved: {output_file}")

def process_batches(batch_dir, ontology_file, output_dir):
    print(" Processing directory:", batch_dir)

    ontology_df = pd.read_csv(ontology_file, sep='\t')
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(batch_dir):
        if filename.endswith('.tsv'):
            batch_path = os.path.join(batch_dir, filename)
            print(f" Reading: {filename}")

            batch_df = pd.read_csv(batch_path, sep='\t')
            batch_df['Genome_ID'] = batch_df['Source']

            coverage_df = compute_function_coverage_per_genome(batch_df, ontology_df)

            output_file = os.path.join(output_dir, f"{filename}_function_coverage.xlsx")
            generate_excel_outputs(coverage_df, output_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="PGPT function coverage per genome batch.")
    parser.add_argument("batch_dir", help="Directory with .tsv genome batch files")
    parser.add_argument("ontology_file", help="Ontology reference .tsv")
    parser.add_argument("output_dir", help="Where to save Excel summaries")
    args = parser.parse_args()

    process_batches(args.batch_dir, args.ontology_file, args.output_dir)

