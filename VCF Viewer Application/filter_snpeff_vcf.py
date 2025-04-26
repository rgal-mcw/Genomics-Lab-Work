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
import csv

#############
# Functions #
#############

## FILTERING (should be VAGUE) ##
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

    # 2. Filter out if MC_details has only one entry and it's 'intron_variant'
    mc_details = row['MC_details'].split(', ')
    if len(mc_details) == 1 and mc_details[0].strip().startswith('intron_variant'):
        return False

    # 3. Filter out if the only entry in IMPACT_details is 'MODIFIER(*)'
    # Check if IMPACT_details matches "MODIFIER(number)" exactly
    if re.fullmatch(r"MODIFIER\(\d+\)", row['IMPACT_details']):
        return False
    
    # If none of the filter conditions matched, keep the row
    return True


## Extract ANN from INFO field
def extract_ann(info):
    ann_match = re.search(r"ANN=([^;]+)", info)
    if ann_match:
        ann = ann_match.group(1)  # Extract only the content after 'ANN='
        # Remove 'ANN=...' from info
        info_without_ann = info[:ann_match.start()] + info[ann_match.end():]
        if info_without_ann.startswith(';'):
            info_without_ann = info_without_ann[1:]
        return ann, info_without_ann
    else:
        return None, info


## Make VCF into dataframe
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
        pattern = rf'{header}=([^;]*)'

        def extract_value(info_str):
            match = re.search(pattern, info_str)
            return match.group(1) if match else None

        input_df[header] = input_df['INFO'].apply(lambda x: extract_value(x))

    return input_df


## Parse the ANN field for information
def parse_ann_field(ann_str):
    if pd.isnull(ann_str) or ann_str.strip() == '':
        return None, None, None

    # Use StringIO and csv.reader to handle complex splitting
    ann_list = next(csv.reader([ann_str], delimiter=',', quotechar='"'))

    records = []
    for ann in ann_list:
        fields = ann.split('|')
        if len(fields) >= 4:
            record = {
                'Annotation': fields[1].strip(),
                'Impact': fields[2].strip(),
                'Gene_Name': fields[3].strip()
            }
            records.append(record)

    if not records:
        return None, None, None

    df_records = pd.DataFrame(records)
    mc_counts = df_records['Annotation'].value_counts().to_dict()
    impact_counts = df_records['Impact'].value_counts().to_dict()
    gene_counts = df_records['Gene_Name'].value_counts().to_dict()
    mc_summary = ', '.join([f"{key}({value})" for key, value in mc_counts.items()])
    impact_summary = ', '.join([f"{key}({value})" for key, value in impact_counts.items()])
    gene_summary = ', '.join([f"{key}({value})" for key, value in gene_counts.items()])
    return mc_summary, impact_summary, gene_summary


## Make gene_list from GENE_details
def extract_genes(gene_details):
    if isinstance(gene_details, str):
        return re.findall(r'(\w+)\(\d+\)', gene_details)
    else:
        return []

## For attaching OMIM information based on all genes in call
def get_omim_info(genes):
    # Fetch OMIM info for each gene in the list, if available
    omim_info_list = [omim_dict.get(gene, '') for gene in genes]
    # Concatenate all OMIM info, filtering out empty results, and add the separator " | "
    return " =+=+ ".join(filter(None, omim_info_list))


## Determine homozygous or heterozygous
def determine_homo_hetero(gt):
    alleles = gt.replace('|', '/').split('/')
    if alleles[0] == alleles[1]:
        return 'Homozygous'
    else:
        return 'Heterozygous'

# Function to split 'SAMPLE' column based on 'FORMAT' column
def split_sample(row):
    format_columns = row['FORMAT'].split(':')
    sample_values = row['Sample'].split(':')
    return pd.Series(sample_values, index=format_columns)


############
## Begin Executable Code ##
############
## Inputs and Checks ##

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
print("Current filters set to:\n1.Remove Empty Ann\n2.Remove MC=Intronic ONLY\n3.Remove if Impact is ONLY Modifier")

name = input("What will you name your .tsv output?: ")
start_time = time.time()

## Read Files ##

# Initialize lists to store header lines and data lines
hash_lines = []
data_lines = []
headers = None

with open(file_path, 'r') as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('##'):
            hash_lines.append(line)
        elif line.startswith('#CHROM'):
            headers = line.lstrip('#').strip().split('\t')
        else:
            data_lines.append(line)

# Check if headers were found
if not headers:
    print("ERROR: Header line starting with #CHROM not found in VCF file.")
    sys.exit(1)

# Extract ids, types, and descriptions from hash_lines
ids = []
types = []
descriptions = []
from_vals = []
from_val = "VCF4.2"

for line in hash_lines:
    if line.startswith('##INFO='):
        id_match = re.search(r'ID=([^,>]+)', line)
        type_match = re.search(r'Type=([^,>]+)', line)
        desc_match = re.search(r'Description="([^"]+)"', line)
        if id_match and type_match and desc_match:
            ids.append(id_match.group(1))
            types.append(type_match.group(1))
            descriptions.append(desc_match.group(1))
            from_vals.append(from_val)
    else:
        if 'SnpEff' in line:
            from_val = "SnpEff"
        elif 'SnpSift' in line:
            if '00-All.vcf.gz' in line:
                from_val = "dbSNP"
            if 'clinvar' in line:
                from_val = "ClinVar"

# Create pandas dataframe 
hash_df = pd.DataFrame({"ID": ids, "Type": types, "Description": descriptions, "from": from_vals})

# Remove duplicates
hash_df.loc[hash_df['ID'] == 'REF', 'ID'] = 'REF_ANN'

# Remove 'ANN' from headers (we'll handle it separately)
headers_df = hash_df[hash_df['ID'] != 'ANN']



print("Obtained VCF Headers")

# Create a pandas dataframe from data_lines, using the extracted headers
data_df = pd.DataFrame([line.split('\t') for line in data_lines], columns=headers)

# Extract ANN from INFO column
data_df["ANN"], data_df['INFO'] = zip(*data_df['INFO'].apply(extract_ann))

print("Building VCF Dataframe")
# Fill columns with data from INFO field
unique_ids = headers_df['ID'].unique()

## Apply vcfCols function ##
for id in unique_ids:
    data_df = vcfCols(data_df, id, headers_df)

print("Parsing ANN Field - Building Details Columns (eta ~= 50mins for full .vcf")

## Add ANN stuff ## 
# Apply the function to the 'ANN' column and create new columns
data_df[['MC_details', 'IMPACT_details', 'GENE_details']] = data_df['ANN'].apply(
    lambda x: pd.Series(parse_ann_field(x))
)

###  Add OMIM annotations ###
print("Beginning OMIM annotation attachment\nAnnotations grabbed from ./db/OMIM_Annotations.tsv - see OMIM.py.")
omim_file = "./db/OMIM_Annotations.tsv"
omim = pd.read_csv(omim_file, sep='\t')

# Apply the extraction to create a new column 'gene_list'
data_df['gene_list'] = data_df['GENE_details'].apply(extract_genes)
omim_dict = omim.set_index('Approved_Gene_Symbol')['OMIM_INFO'].to_dict()
print("Applying OMIM Annotation")
data_df['OMIM_INFO'] = data_df['gene_list'].apply(get_omim_info)
data_df.drop(columns=['gene_list'], inplace=True)
len_data_df = len(data_df)



print("Finished Creating VCF Dataframe \nApplying Filters")

### Applying Custom Filters ###

# Apply the filter function to each row in the dataframe
data_df_filtered = data_df[data_df.apply(filter_vcf, axis=1)]
len_data_df_filtered = len(data_df_filtered)

# Split FORMAT and SAMPLE columns
split_samples = data_df_filtered.apply(split_sample, axis=1)

# Concatenate the original DataFrame with the new split columns & Drop Format/Sample
data_df_filtered = pd.concat([data_df_filtered, split_samples], axis=1)
data_df_filtered.drop(columns=['FORMAT', 'Sample'], inplace=True)

#### ADD CUSTOM FIELDS HERE ####

#data_df_filtered['Homo_Hetero'] = data_df_filtered['GT'].apply(determine_homo_hetero)

##### BEGIN CONVERSION TO TSV ####

# Print to .tsv
data_df_filtered.to_csv(f"{output_dir}/{name}_vcf_calls.tsv", sep="\t", index=False)

### Filter original .vcf for matching CHROM + POS in data_df_filtered
chrom_pos_set = set(zip(data_df_filtered['CHROM'], data_df_filtered['POS']))

print("Parsing Calls in OG Dataframe. Writing Filtered VCF")
filtered_vcf_lines = []

with open(file_path, 'r') as vcf_file:
    for line in vcf_file:
        # Write header lines
        if line.startswith("#"):
            filtered_vcf_lines.append(line)
            continue

        # Split the line and extract CHROM and POS
        line_parts = line.strip().split("\t")
        if len(line_parts) > 1:
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
