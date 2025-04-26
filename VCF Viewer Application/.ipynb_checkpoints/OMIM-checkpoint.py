import os
import re
import sys
from io import StringIO
import pandas as pd
from pathlib import Path

Print("Hi! This file is currently written to work on my MacBook (ry28926). This script currently:\n1. Finds and verifies OMIM files (downloaded from OMIM.org).\n2. Converts genemap2.txt to a dataframe.\n3. Deletes some columns from genemap2_df.\n4. Writes OMIM_ANN df and saves to OMIM_Annotatons.tsv. This is used to attach OMIM annotations to our annotated .vcfs.\n\n")

## Read in OMIM Files
dir = Path("/Users/ry28926/Desktop/Project/maple/OMIM")
print(f"OMIM directory set to: {dir}\n")

genemap2 = dir / "genemap2.txt"
#mim2gene = dir / "mim2gene.txt"
#mimTitles = dir / "mimTitles.txt"
#morbidmap = dir / "morbidmap.txt"

# Check existence
print("** Verifying Existence of OMIM files **")
#files = [genemap2, mim2gene, mimTitles, morbidmap]
file = [genemap2]

for file_name in files:
    if  file_name.exists():
        print(f"{file_name} found")
    else: 
        print(f"{file_name} does not exists.")
        sys.exit()
print(f"** All files found **")

# Format Headers & Write to pandas DF
for file_name in files:
    with open(file_name) as f:
        lines = f.readlines()
    
    header = []
    data_lines = []
    
    for i in range(len(lines)):
        # Keep updating until the last "#" line, so saves as last (headers)
        if lines[i].startswith("#"):
            header = lines[i].lstrip("#").strip() + "\n"
        else:
            data_lines = lines[i:] #Remaining lines are data
            break

    # Combine headers and data lines
    cleaned_data = header + "".join(data_lines)

    # Write to DF
    df = pd.read_csv(StringIO(cleaned_data), sep='\t')
    del cleaned_data
    df_name = f"{os.path.splitext(os.path.basename(file_name))[0]}_df"
    globals()[df_name] = df

# Format Genemap2 
genemap2_df = genemap2_df.drop(columns=['Genomic Position Start', 'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'Mouse Gene Symbol/ID'])
# Create OMIM_ANN from Formatted Genemap2
OMIM_ANN = genemap2_df[['Chromosome', 'Approved Gene Symbol']].assign(
    OMIM_INFO=genemap2_df.drop(columns=['Chromosome', 'Approved Gene Symbol']).astype(str).agg(' | '.join, axis=1))

OMIM_ANN = OMIM_ANN.rename(columns={"Approved Gene Symbol": "Approved_Gene_Symbol"})

# Write OMIM ANN
out_path = dir / "OMIM_Annotations.tsv"
OMIM_ANN.to_csv(out_path, sep='\t', index=False)

print("DONT FORGET TO COPY THIS OUTPUT TO toucan/db/ FOR IT TO BE USED IN THE ANNOTATION SCRIPT!!!!!!!")