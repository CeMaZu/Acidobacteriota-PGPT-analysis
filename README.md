# Acidobacteriota-PGPT-analysis
Scripts and workflows for the genome-wide identification and comparative analysis of plant growth-promoting traits (PGPTs) in Acidobacteriota using DIAMOND and PLABASE.

## Reference data

The `reference_data/` directory contains ontology and mapping files required for PGPT annotation and downstream analyses. In particular, it includes the PLaBAse ontology used to assign protein hits to plant growth-promoting trait (PGPT) categories.

PLaBAse is a curated web resource for analysing the plant growth-promoting potential of plant-associated bacteria.

Reference:
Patz S, Gautam A, Becker M, Ruppel S, Rodríguez-Palenzuela P, Huson DH. PLaBAse: A comprehensive web resource for analyzing the plant growth-promoting potential of plant-associated bacteria. bioRxiv, 2021.

ReadMe PGPT Pipeline

--------------------------------------------------
PIPELINE ORDER (IMPORTANT)
--------------------------------------------------

Correct order of execution:

1. Run main PGPT pipeline (01-02)
   → produces: merged_results/

2. Run taxonomy merging script (05)
   → produces: Merged_Tax/

3. Optional Run taxonomic batch script (06) (if you want to analyse taxonomic groups)
   → produces: Merged_Tax_Batches/

4. Run Function Coverage to determine the number of PGPT detected with the theoretical maximum (03)
   → could also be run on 01-02 output but then there will be no info on taxonomic affiliation. 

If you skip step 2, the batch script will fail or produce empty outputs.


# PLABASE / PGPT DIAMOND Processing Pipeline - 01-02_GenerationOfPlaBASEAnalysisFiles_Clean-v11.py
====================

	This script extracts protein-coding sequences from GenBank (.gbff) files, runs DIAMOND blastp searches against the PLABASE PGPT protein database, filters the hits, merges results with a PLABASE ontology file, and produces final PGPT family summaries.

	--------------------------------------------------
	REQUIRED SOFTWARE
	--------------------------------------------------

	1. Python
	   Version: Python >= 3.8

	2. Python packages
	   Install with:
	   pip install pandas biopython

	   Required libraries:
	   - pandas
	   - biopython

	3. DIAMOND
	   Must be installed and available in PATH.

	   Check:
	   diamond --version

	   Install with conda:
	   conda install -c bioconda diamond

	--------------------------------------------------
	REQUIRED INPUT FILES AND DATABASES
	--------------------------------------------------

	1. GenBank input folder

	   Default:
	   PLABASE_gbff_Input/

	   Must contain:
	   *.gbff files

	   Each file must include CDS features with:
	   - locus_tag
	   - gene
	   - product
	   - translation

	2. DIAMOND database

	   Default path needs to be changed in the script

	   This must already exist (not created by the script).

	   If needed, create with:
	   diamond makedb --in input_proteins.faa -d PGPT_BASE_nr_Aug2021n_ul_1

	3. PLABASE ontology file

	   File:
	   PLABASE-Ontology-finaln6911.txt

	   Must be tab-separated and include columns:
	   - ID
	   - PGPT_ID
	   - PATHS

	Provided with the script 

	--------------------------------------------------
	EXPECTED DIRECTORY STRUCTURE
	--------------------------------------------------

	project/
	├── GenerationOfPlaBASEAnalysisFiles_Clean-v11.py
	├── PLABASE-Ontology-finaln6911.txt
	├── PLABASE_gbff_Input/
	│   ├── genome_1.gbff
	│   ├── genome_2.gbff
	│   └── genome_3.gbff
	└── PGPT_BASE_nr_Aug2021n_ul_1.dmnd

	--------------------------------------------------
	OUTPUT FOLDERS
	--------------------------------------------------

	If an output folder is provided, all results are created inside it.

	1. DIAMOND_FASTA_Input/
	   Extracted protein FASTA files from .gbff files

	2. diamond_results/
	   Raw DIAMOND output

	3. renamed_diamond_results/
	   Renamed DIAMOND output files

	4. processed_diamond_results/
	   Parsed DIAMOND results with metadata

	5. filtered_diamond_results/
	   Filtered hits (E-value ≤ 1e-50)

	6. final_processed_results/
	   Adds PGPT_Family and PGPT_Number

	7. merged_results/
	   Final merged output with ontology

	   Contains:
	   - merged .tsv files
	   - best_*.tsv (best hit per PGPT family)
	   - PGPT_summary.txt

	--------------------------------------------------
	LOG FILE - ONLY Generated if .gbff files are used
	--------------------------------------------------

	extraction.log

	Contains number of CDS extracted per genome.

	--------------------------------------------------
	HOW TO RUN
	--------------------------------------------------

	Run:
	GenerationOfPlaBASEAnalysisFiles_Clean-v11.py

	Prompts:

	1. DIAMOND sensitivity mode:
	   Options:
	   faster, fast, mid-sensitive, sensitive,
	   more-sensitive, very-sensitive, ultra-sensitive

	   Press Enter for default.

	2. Skip FASTA extraction:
	   yes / no

	3. Output folder name:
	   Example:
	   run_01

			--------------------------------------------------
			USING EXISTING FASTA FILES (SKIP EXTRACTION)
			--------------------------------------------------

			The pipeline allows you to skip the GenBank (.gbff) → FASTA extraction step and use pre-existing FASTA files.

			This is useful if:
			- FASTA files were already generated previously
			- You want to rerun DIAMOND with different parameters
			- You are working with externally generated protein FASTA files

			--------------------------------------------------
			REQUIREMENTS FOR FASTA FILES
			--------------------------------------------------

			FASTA files must:

			1. Be located in a single directory
			2. Have the extension:
			   .fasta

			3. Contain protein sequences (NOT nucleotide sequences)

			4. Use the exact header format expected by the pipeline:

			   >GenomeID_LocusTag_GeneName | Source: filename.gbff | c

			Example:

			   >NC_000001_tag123_geneX | Source: genome1.gbff | Product: hypothetical protein

			--------------------------------------------------
			IMPORTANT: HEADER FORMAT
			--------------------------------------------------

			The script parses FASTA headers using a strict pattern.

			If headers do NOT follow this structure:
			- Downstream parsing will fail
			- Metadata columns (Genome, CDS, Gene, Source, Product) may be incorrect or missing

			If your FASTA files were generated outside this pipeline, you may need to reformat headers.

			--------------------------------------------------
			HOW TO USE EXISTING FASTA FILES
			--------------------------------------------------

			Run the script:

			python script.py

			When prompted:

			Skip FASTA extraction? (yes/no):
			→ Enter: yes

			Then provide the FASTA directory:

			Enter existing FASTA input directory:
			→ Example:
			DIAMOND_FASTA_Input/

			The script will:
			- Skip sequence extraction
			- Use the provided FASTA files directly for DIAMOND

			--------------------------------------------------
			EXPECTED DIRECTORY STRUCTURE (FASTA MODE)
			--------------------------------------------------

			project/
			├── script.py
			├── PLABASE-Ontology-finaln6911.txt
			├── DIAMOND_FASTA_Input/
			│   ├── genome_1.gbff.fasta
			│   ├── genome_2.gbff.fasta
			│   └── genome_3.gbff.fasta
			└── PGPT_BASE_nr_Aug2021n_ul_1.dmnd

			--------------------------------------------------
			COMMON ERRORS
			--------------------------------------------------

			1. Directory does not exist

			   Error:
			   Directory 'X' does not exist

			   Fix:
			   Check the path you entered

			2. No FASTA files detected

			   Ensure files end with:
			   .fasta

			3. Incorrect header format

			   Symptoms:
			   - Missing columns in output
			   - Parsing errors in process_results step

			   Fix:
			   Reformat FASTA headers to match expected pattern

			--------------------------------------------------
			TIP
			--------------------------------------------------

			If you are unsure whether your FASTA files are compatible, run one file through the full pipeline first and inspect:

			processed_diamond_results/

			If columns like Genome_Accession, CDS, or Gene_Name look broken, your FASTA headers might be the problem.
			

	--------------------------------------------------
	DIAMOND PARAMETERS
	--------------------------------------------------

	Command:
	diamond blastp

	Main settings:
	- --outfmt 6
	- --evalue 1e-5
	- --threads 8

	Final filtering:
	E-value ≤ 1e-50


	--------------------------------------------------
	FINAL OUTPUT
	--------------------------------------------------

	Most relevant results:

	merged_results/*.tsv
	merged_results/best_*.tsv
	merged_results/PGPT_summary.txt
	
For some questions taxonomic batching or custom might be required for comparative analysis. To do so, two subsequent scripts are used. 

05-TaxonomyMerge-v2_1.py
================

Purpose
-------
TaxonomyMerge-v2 merges GTDB taxonomy information into DIAMOND/PGPT result files.

The script supports both:

1. NCBI genome workflows using GCA accession numbers
2. Custom FASTA workflows where no GCA accession is available

Taxonomy information is inserted directly after the Source column to simplify
downstream analyses and taxonomic grouping.


Input Files
-----------

1. DIAMOND/PGPT Results

Input directory containing one or more TSV files.

Required column:

Source

Example:

HAG010016_QUXAI209AT
FH1421_QUXAI207AD
GCA_018268895.1_ASM1826889v1_genomic.gbff


2. GTDB Taxonomy Table

Required columns:

user_genome
classification
Phylum
Class
Order
Family
Genus
Species

Optional column:

Custom_Group

Example:

user_genome    Family               Custom_Group
HAG010016      Streptomycetaceae    Group2
FH1421         Streptomycetaceae    Group4


Matching Strategy
-----------------

The script uses several matching approaches:

1. Exact Match

Source:
HAG010016

GTDB:
HAG010016


2. Prefix Match

Source:
HAG010016_QUXAI209AT

GTDB:
HAG010016


3. Normalized Match

Source:
P.alni_DSM44104_PHUJ01

GTDB:
PAlni_DSM44104-PHUJ01.1

Special characters such as:
.
_
-

are removed before comparison.


Added Columns
-------------

The following columns are inserted immediately after Source:

Custom_Group
Domain
Phylum
Class
Order
Family
Genus
Species
Matched_user_genome
Match_Type

The GTDB column "classification" is interpreted as:

Domain

Example:

d__Bacteria

becomes:

Bacteria

Taxonomic prefixes are removed:

d__
p__
c__
o__
f__
g__
s__


Output
------

For each input TSV file, a taxonomy-annotated TSV file is generated.

Example:

Input:
Sample1.tsv

Output:
Sample1.tsv

with additional taxonomy columns inserted after Source.


Match Report
------------

The script automatically generates:

TaxonomyMerge_match_report.tsv

The report contains:

File
Source
Possible_IDs
Matched_Query
Matched_user_genome
Match_Type
Custom_Group
Domain
Phylum
Class
Order
Family
Genus
Species

This report can be used to identify unmatched genomes and verify taxonomy assignments.


Usage
-----

Basic:

python3 TaxonomyMerge-v2.py INPUT_DIRECTORY GTDB_RESULTS.txt OUTPUT_DIRECTORY

Example:

python3 TaxonomyMerge-v2.py \
DiamondResults \
GTDB_GlCNacCorrect.txt \
Merged_Taxonomy


Optional custom label:

python3 TaxonomyMerge-v2.py \
DiamondResults \
GTDB_GlCNacCorrect.txt \
Merged_Taxonomy \
--custom_group GlcNAc

The command-line label is only used if no Custom_Group column exists in the taxonomy file.


Notes
-----

- Supports both GCA-based and FASTA-based workflows.
- Taxonomy information is inserted directly after Source.
- Supports grouping by:
  Domain
  Phylum
  Class
  Order
  Family
  Genus
  Species
  Custom_Group
- Automatically generates a match report for quality control.
- Recommended for all downstream PGPT coverage and taxonomic analyses.	

06-TaxonomicBatching-v2_1.py
====================

Purpose
-------
TaxonomicBatching-v2 creates taxonomic batch files from taxonomy-annotated PGPT/DIAMOND TSV files.

The script is designed to run after TaxonomyMerge-v2. It uses the taxonomy columns added to each TSV file to generate batch files by taxonomic rank or custom grouping.

The resulting batch files can be used as input for the PGPT function coverage analysis script.


Input
-----

Input directory containing taxonomy-annotated TSV files.

These files should already contain taxonomy columns added by TaxonomyMerge-v2.

Required columns:

Source
PATHS

Recommended taxonomy columns:

Custom_Group
Domain
Phylum
Class
Order
Family
Genus
Species
Matched_user_genome
Match_Type


Output
------

The script creates one output folder per selected taxonomic rank.

Example:

Taxonomic_Batches/
  Family/
    Family_Streptomycetaceae_batch.tsv
    Family_Pseudonocardiaceae_batch.tsv

  Genus/
    Genus_Streptomyces_batch.tsv
    Genus_Pseudonocardia_batch.tsv

  Custom_Group/
    Custom_Group_Group1_batch.tsv
    Custom_Group_Group2_batch.tsv

  Reports/
    Family_batch_summary.tsv
    Genus_batch_summary.tsv
    Custom_Group_batch_summary.tsv
    ALL_taxonomic_batch_summary.tsv
    ALL_unknown_taxonomy_report.tsv


Batch Files
-----------

Each batch file contains all rows belonging to one taxonomic group.

Example:

Family_Streptomycetaceae_batch.tsv

contains all PGPT/DIAMOND hits assigned to the family Streptomycetaceae.

Batch files retain the original columns, including:

Source
taxonomy columns
PGPT annotation columns
PATHS

These files are suitable for downstream PGPT function coverage analysis.


Reports
-------

All report files are stored in:

Reports/

The script generates summary reports for each rank and combined reports across all ranks.

Summary reports include:

Rank
Group
Batch_File
Rows
Genome_Count

The Genome_Count column reports the number of unique genomes assigned to each batch.

Unknown reports list genomes or rows where the selected taxonomic rank was missing, Unknown, or unclassified.


Taxonomic Ranks
---------------

By default, the script can generate batches for all available ranks:

Custom_Group
Domain
Phylum
Class
Order
Family
Genus
Species

You can also select specific ranks.

Examples:

Family only
Order only
Family and Custom_Group
Genus and Species


Usage
-----

Basic usage:

python3 TaxonomicBatching-v2.py INPUT_DIRECTORY OUTPUT_DIRECTORY

This uses the default settings:

--ranks all
--deduplicated best


Example:

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches


Select specific ranks:

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches \
--ranks Family,Custom_Group

Batch only by Order:

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches \
--ranks Order

Process all TSV files:

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches \
--deduplicated all


Process only non-deduplicated/raw files:

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches \
--deduplicated raw


Deduplication Mode
------------------

The script supports three file-selection modes:

best
  Only processes files starting with best_

raw
  Only processes files not starting with best_

all
  Processes all TSV files

Default:

best

This is recommended for PGPT coverage analysis because deduplicated best-hit files prevent inflated counts caused by redundant hits.


Unknown or Unclassified Entries
-------------------------------

Entries are considered unknown if the selected rank contains:

Unknown
unclassified
NA
NaN
empty value

These entries are written to unknown report files and, where applicable, to Unknown/unclassified batch files.


Important Notes
---------------

- This script does not calculate PGPT coverage.
- It only creates batch files for downstream analysis.
- It should be run after TaxonomyMerge-v2.
- The output batch files can be used as input for the PGPT function coverage script.
- Report files are saved separately in the Reports folder so that downstream scripts do not accidentally process them as batch files.
- For publication workflows, use deduplicated best_ files unless intentionally comparing raw outputs.


Example Workflow
----------------

python3 TaxonomyMerge-v2.py \
DIAMOND_Outputs \
GTDB_results.txt \
Merged_Taxonomy \
--custom_group GlcNAc

python3 TaxonomicBatching-v2.py \
Merged_Taxonomy \
Taxonomic_Batches \
--ranks Family,Custom_Group \
--deduplicated best

python3 03_functionCoverageScript.py \
Taxonomic_Batches/Family \
PLABASE_Ontology.txt \
Function_Coverage_Family
	

# PGPT Function Coverage Calculation

## Overview
This script calculates PGPT (Plant Growth-Promoting Traits) functional coverage for one or more genomes based on DIAMOND search results against the PLABASE database.

It works on merged DIAMOND output files merged_results, but also on Output files that have been merged with Taxonomy data or grouped using TaxonomyBatch_v2.py

Coverage is computed per genome and subsequently aggregated across genomes to obtain:
- mean PGPT coverage
- standard deviation
- total observed PGPT counts

The output is written as an Excel file with separate sheets for each PGPT hierarchy level.

---

## Input Files

### 1. Batch input files (.tsv)
Directory containing one or more tab-separated files with PGPT hit annotations.

Required columns:
- Source → genome identifier (used as Genome_ID)
- PATHS → hierarchical PGPT annotation (PLABASE format)

Example PATHS entry:
PGPT;DIRECT_EFFECTS;BIO-FERTILIZATION#PGPT;INDIRECT_EFFECTS;STRESS_CONTROL

---

### 2. Ontology file (.tsv)
PLABASE ontology file containing all PGPT categories. - Provided with scripts

Required column:
- PATHS → defines hierarchical PGPT structure

This file is used to determine the total number of PGPT families per functional category, which serves as the denominator for coverage calculations.

---

## Output

For each input .tsv file, an Excel file is generated:
<filename>_function_coverage.xlsx

Each Excel file contains:
- One sheet per PGPT level (Level_0, Level_1, Level_2, ...)
- Columns:
  - Genome_Count
  - Level
  - Function
  - Mean_Coverage_Percentage
  - Std_Coverage_Percentage
  - Sum_Observed_PGPTs

---

## Method Summary

### 1. Parsing PGPT annotations
The PATHS column is parsed to extract hierarchical PGPT levels.

### 2. Per-genome coverage calculation
For each genome:
- PGPT occurrences are counted per functional category
- Coverage is calculated as:

Coverage (%) = (Observed PGPTs / Total PGPTs in ontology) × 100

### 3. Aggregation across genomes
For multiple genomes:
- Mean coverage is calculated as:

Mean Coverage (%) =
   Sum of observed PGPTs across genomes
   -------------------------------------
   (Total PGPTs × Number of genomes)

- Standard deviation is calculated from per-genome coverage values

Important:
Coverage is calculated per genome first, then aggregated.
This avoids inflation caused by cumulative counting across genomes.

---

## Usage

python compute_pgpt_coverage_per_genome.py <batch_dir> <ontology_file> <output_dir>

Arguments:
- batch_dir → directory containing .tsv input files
- ontology_file → PLABASE ontology file
- output_dir → directory for Excel output files

---

## Example

python compute_pgpt_coverage_per_genome.py \
    Tax_Batch/ \
    PLABASE_Ontology.tsv \
    Output_Coverage/


## Notes

- Input files should contain non-redundant PGPT hits (e.g. best_ DIAMOND output)
- The script supports:
  - single genomes
  - genome batches
  - taxonomic groupings

