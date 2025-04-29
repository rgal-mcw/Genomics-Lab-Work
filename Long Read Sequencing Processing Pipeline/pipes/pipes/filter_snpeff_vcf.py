#############
### Title: Filtering snpeff .vcfs
### Author: Ryan Gallagher
### Date: 04-09-2024
#############

## This script is the python version of format_snpeff_vcf.R. R does not handle large dataframes well,
## so I am using python to handle the large files.

## INFO Fields: ##

# The .vcf annotations are found at the top of the .vcf file
# Our .vcf files are annotated using 
#   - ClinVar GRCh38: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240325.vcf.gz
#     -> With fields names: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/README_VCF.txt
#
#   - dbSNP: https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/#file-update
#     -> With field names: https://www.ncbi.nlm.nih.gov/variation/docs/oldglossary_dbSNP1_vcf/
#
#
## USAGE: 
#   `python format_snpeff_vcf.py <file_path>`

# The INFO field on our *.snpEFF.vcf have annotations added to the INFO field, and the interpretation
# comes from understanding these fields. 


import os 
import re
import time
import sys
import pandas as pd
import numpy as np
import pyranges as pr
# Change to the directory with the .vcf files
#os.chdir("~/Desktop/Project/SNP_Anno/data/annotations/")
#Use absolute file path, or let user input file path

# Check for correct number of script arguments and file type
if len(sys.argv) != 2:
    print("ERROR: Insufficient Arguments \nUsage: python format_snpeff_vcf.py <file_path to .vcf>")
    sys.exit(1)

file_path = os.path.abspath(sys.argv[1])

if not os.path.exists(file_path):
    print(f"ERROR: File {file_path} does not exist.")
    sys.exit(1)

if not file_path.endswith(".vcf"):
    print(f"ERROR: File {file_path} is not a .vcf file.")
    sys.exit(1)

# Create output file file
output_path = file_path.replace(".vcf", ".filtered.vcf")
output_dir = os.path.dirname(output_path)
# Need to figure out how to allow custom filters in arguments
print("Current filters set to:\n1.Remove 'CLNSIG=benign' OR 'CLNSIG=NULL', unless HIGH Impact AND GENEINFO not blank.\n2.Remove MC=Intronic")

name = input("What will you name your .tsv output?: ")
start_time = time.time()

# Read in lines from .vcf
with open(file_path, "r") as file:
    lines = [line.strip() for line in file.readlines()]


# Separate lines that start with "##" and those that do not
hash_lines = [line for line in lines if re.match("^##", line)]
data_lines = [line for line in lines if not re.match("^##", line)]

print("Finished reading vcf lines.")

## Headers ##

# Extract ids, types, and descriptions from each line in hash_lines
ids = [re.sub(r".*ID=([^,]+),.*", r"\1", line) for line in hash_lines]
types = [re.sub(r".*Type=([^,]+),.*", r"\1", line) for line in hash_lines]
descriptions = [re.sub(r".*Description=\"([^\"]+)\".*", r"\1", line) for line in hash_lines]


# Remove the first elements
ids = ids[1:]
types = types[1:]
descriptions = descriptions[1:]

# Create pandas dataframe 
hash_df = pd.DataFrame({"ID": ids, "Type": types, "Description": descriptions, "from": 'VCF4.2'})
from_val = "VCF4.2"

# Loop through hash_df and update `from` column
for i, row in hash_df.iterrows():
    id_val = row['ID']

    if id_val.startswith("##"):
        if 'SnpEff' in id_val:
            from_val = "SnpEff"
        if 'SnpSift' in id_val:
            if '00-All.vcf.gz' in id_val:
                from_val = "dbSNP"
            if 'clinvar' in id_val:
                from_val = "ClinVar"

    hash_df.at[i, 'from'] = from_val


# Remove ## rows 
hash_df = hash_df[~hash_df['Description'].str.startswith("##")]

# Update duplicate column names
hash_df.loc[hash_df['ID'] == 'REF', 'ID'] = 'REFA'

# Filter out VCF4.2 rows
headers = hash_df[hash_df['from'] != 'VCF4.2']
print("Headers formatted")


## Put the data in a dataframe

# Create a pandas dataframe from data_lines, then rearrange columns
data = data_lines[1:]
data_df = pd.DataFrame([line.split("\t") for line in data])
data_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
data_df = data_df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT", "SAMPLE", "INFO"]]

# Extract ANN from INFO column
def extract_ann(info):
    ann_match = re.search(r"ANN=([^;]+)", info)
    if ann_match:
        ann = ann_match.group(1)  # Extract only the content after 'ANN='
        return ann, info.replace(ann_match.group(0), "")
    else:
        return None, info

data_df["ANN"], data_df['Else'] = zip(*data_df['INFO'].apply(extract_ann))

# Drop INFO column and rename Else to INFO
data_df.drop(columns=["INFO"], inplace=True)
data_df.rename(columns={"Else": "INFO"}, inplace=True)



# Filter ANN from headers
headers_df = headers[headers['ID'] != 'ANN']

# Iterate over rows and add columns to data_df
for i, row in headers_df.iterrows():
    data_df[row['ID']] = pd.NA

print("Beginning column fill")
# Fill columns with data from INFO field
def vcfCols(input_df, header, headers_df):
    if not isinstance(header, str):
        raise ValueError("header must be a string")
    if not isinstance(headers_df, pd.DataFrame):
        raise ValueError("headers_df must be a pandas DataFrame")

    df = headers_df[headers_df['ID'] == header]

    if len(df) > 1:
        raise ValueError("ERROR: Duplicate IDs in dataframe")

    # Flag Conditions
    if df['Type'].values[0] == 'Flag':
        input_df[header] = input_df['INFO'].apply(lambda x: 1 if header in x else 0)

    # String, Integer, Float Conditions

    if df['Type'].values[0] in ['String', 'Integer', 'Float']:
        pattern = rf'{header}=([^;]*);.*'

        def extract_value(info_str):
            match = re.search(pattern, info_str)
            return match.group(1) if match else None

        input_df[header] = input_df['INFO'].apply(lambda x: extract_value(x))

    return input_df

unique_ids = headers_df['ID'].unique()

for id in unique_ids:
    data_df = vcfCols(data_df, id, headers_df)

print("building ANN_details columns")
## Add ANN stuff
def count_annotation_details(ann):
    annotations = [ann_entry for ann_entry in str(ann).split(',') if ann_entry]

    mc_counts = {}
    impact_counts = {}
    gene_counts = {}
    
    for ann_entry in annotations:
        fields = ann_entry.split('|')
        # SnpEff ANN fields should have at least 4 fields: Allele | Annotation | Impact | Gene_Name
        if len(fields) >= 4:
            mc = fields[1].strip()  # Annotation
            impact = fields[2].strip()  # Annotation_Impact
            gene_name = fields[3].strip()  # Gene_Name
            
            if mc:
                mc_counts[mc] = mc_counts.get(mc, 0) + 1
                
            if impact:
                impact_counts[impact] = impact_counts.get(impact, 0) + 1
                
            if gene_name:
                gene_counts[gene_name] = gene_counts.get(gene_name, 0) + 1
        else:
            # Log or handle annotations with insufficient fields
            print(f"Warning: Incomplete ANN entry skipped: {ann_entry}")
    
    # Format counts as 'VALUE(COUNT)' strings
    mc_summary = ', '.join([f"{key}({value})" for key, value in mc_counts.items()]) if mc_counts else None
    impact_summary = ', '.join([f"{key}({value})" for key, value in impact_counts.items()]) if impact_counts else None
    gene_summary = ', '.join([f"{key}({value})" for key, value in gene_counts.items()]) if gene_counts else None
        
    return mc_summary, impact_summary, gene_summary

# Apply the function to each row of the ANN column
data_df[['MC_details', 'IMPACT_details', 'GENE_details']] = data_df['ANN'].apply(
    lambda x: pd.Series(count_annotation_details(x))
)

print("finished building ANN_details columns")

## Add OMIM annotations
print("Beginning OMIM annotation attachment")
omim_file = "/data/ref/OMIM/OMIM_Annotations.tsv"
omim = pd.read_csv(omim_file, sep='\t')

# Make gene_list from GENE_details
def extract_genes(gene_details):
    if isinstance(gene_details, str):
        return re.findall(r'(\w+)\(\d+\)', gene_details)
    else:
        return []
# Apply the extraction to create a new column 'gene_list'
data_df['gene_list'] = data_df['GENE_details'].apply(extract_genes)
data_df['first_gene'] = data_df['gene_list'].apply(lambda genes: genes[0] if genes else np.nan)

if omim['Approved_Gene_Symbol'].duplicated().any():
    print("WARNING: Duplicate gene symbols found in omim. Only the first occurrence will be used.")
    omim = omim.drop_duplicates(subset='Approved_Gene_Symbol')

omim_mapping = omim.set_index('Approved_Gene_Symbol')['OMIM_INFO']

print("Applying OMIM_INFO")
data_df['OMIM_INFO'] = data_df['first_gene'].map(omim_mapping)

##
data_df.drop(columns=['gene_list', 'first_gene'], inplace=True)
len_data_df = len(data_df)


print("Finished column fill \nApplying Custom Filters")

### Applying Custom Filters

def filter_vcf(row):
    """
    Filter rows based on the following criteria:
    
    1. Filter out if ANN is empty.
    2. If IMPACT_details is empty, filter out unless ANN contains specific strings.
    3. Filter out if IMPACT_details only has LOW or MODIFIER or a combination of the two.
    4. Filter out if MC_details contains 'intron'.
    
    Parameters:
    row (pd.Series): A row from the DataFrame.
    
    Returns:
    bool: True if the row should be kept, False otherwise.
    """
    
    # 1. Filter out if ANN is empty
    if pd.isna(row['ANN']) or row['ANN'].strip() == '':
        return False
    
    # 2. If IMPACT_details is empty, filter out unless ANN contains specific strings
    if pd.isna(row['IMPACT_details']) or row['IMPACT_details'].strip() == '':
        # Define the specific substrings to look for in ANN
        required_substrings = ['|HIGH|', 'missense_variant', 'nonsense', 'stop_lost', 'stop_gained']
        ann = row['ANN']
        # Check if any of the required substrings are present in ANN
        if not any(substring in ann for substring in required_substrings):
            return False
    
    # 3. Filter out if IMPACT_details only has LOW or MODIFIER or a combination
    if pd.notna(row['IMPACT_details']) and row['IMPACT_details'].strip() != '':
        # Extract individual impact levels from IMPACT_details
        # Example format: "HIGH(1), MODERATE(2)"
        impacts = [item.split('(')[0].strip().upper() for item in row['IMPACT_details'].split(',')]
        unique_impacts = set(impacts)
        # Define the set of impacts that should lead to filtering out
        filter_out_impacts = {'LOW', 'MODIFIER'}
        # If all unique impacts are in the filter_out_impacts set, filter out the row
        if unique_impacts and unique_impacts.issubset(filter_out_impacts):
            return False
    
    # 4. Filter out if MC_details contains 'intron'
    if pd.notna(row['MC']) and 'intron' in row['MC'].lower():
        return False
    
    # If none of the filter conditions matched, keep the row
    return True

# Apply the filter function to each row in the dataframe
data_df_filtered = data_df[data_df.apply(filter_vcf, axis=1)]
len_data_df_filtered = len(data_df_filtered)

# Split FORMAT column (DEBUGGING)

temp_df = data_df_filtered['SAMPLE'].str.split(":", expand=True)

# Function to split 'SAMPLE' column based on 'FORMAT' column
def split_sample(row):
    format_columns = row['FORMAT'].split(':')
    sample_values = row['SAMPLE'].split(':')
    return pd.Series(sample_values, index=format_columns)

# Apply the function to each row
split_samples = data_df_filtered.apply(split_sample, axis=1)

# Concatenate the original DataFrame with the new split columns & Drop Format/Sample
data_df_filtered = pd.concat([data_df_filtered, split_samples], axis=1)
data_df_filtered.drop(columns=['FORMAT', 'SAMPLE'], inplace=True)


#### ADD CUSTOM FIELDS HERE ####

## Determine homozygous or heterozygous
def determine_homo_hetero(gt):
    alleles = gt.replace('|', '/').split('/')
    if alleles[0] == alleles[1]:
        return 'Homozygous'
    else:
        return 'Heterozygous'

data_df_filtered['Homo_Hetero'] = data_df_filtered['GT'].apply(determine_homo_hetero)


##### BEGIN CONVERSION TO TSV #########

# Print to .tsv
data_df_filtered.to_csv(f"{output_dir}/{name}_vcf_calls.tsv", sep="\t", index=False)

### Filter original .vcf for matching CHROM + POS in data_df_filtered
chrom_pos_set = set(zip(data_df_filtered['CHROM'], data_df_filtered['POS'].astype(str)))


print("Filtering Original .vcf & Writing filtered .vcf file")
filtered_vcf_lines = []

with open(file_path, 'r') as vcf_file:
    for line in vcf_file:
        # Skip header lines
        if line.startswith("#"):
            filtered_vcf_lines.append(line)
            continue

        # Split the lines and extract CHROM and POS
        line_parts = line.split("\t")
        chrom = line_parts[0]
        pos = line_parts[1]

        # Check if CHROM + POS is in chrom_pos_set
        if (chrom, pos) in chrom_pos_set:
            filtered_vcf_lines.append(line)

# Write filtered_vcf_lines to a new file
with open(output_path, 'w') as filtered_vcf_file:
    for line in filtered_vcf_lines:
        filtered_vcf_file.write(line)


end_time = time.time()
duration = round(end_time - start_time, 2)
print(f"{len_data_df_filtered} rows filtered from {len_data_df} rows")
print(f"Script Duration: {duration} seconds")
