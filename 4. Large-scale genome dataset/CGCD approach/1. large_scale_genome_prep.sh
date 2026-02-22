#!/bin/bash

### Large-scale genome preparation
# All available A.baumannii genomes were first downloaded from the NCBI database as compressed archive files.

# Unzip all NCBI downloads
mkdir -p unzipped_all

for zip in ncbi_*.zip; do
  species_dir="${zip%.zip}"
  unzip -o "$zip" -d "unzipped_all/$species_dir"
done

# Collect all .fna files
mkdir -p fasta
find unzipped_all/ -name "*.fna" -exec cp {} fasta/ \;

#Rename genomes using FASTA headers
input_dir="fasta"
output_dir="fasta/renamed"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all .fna files in the current directory
for file in "$input_dir"/*.fna; do
  filename=$(basename "$file")

  # Read the first header line
  header=$(grep -m 1 "^>" "$file")

  # Extract species (first word after 'Acinetobacter')
  species=$(echo "$header" | sed -n 's/.*Acinetobacter \([a-zA-Z]\+\).*/\1/p')

  # Try to extract strain if the word "strain" exists
  if echo "$header" | grep -qi "strain"; then
    strain=$(echo "$header" | sed -n 's/.*strain \([^,>]*\).*/\1/p')
  else
    # No "strain" keyword; take word(s) after species
    strain=$(echo "$header" | sed -n "s/.*Acinetobacter $species \([^,>]*\).*/\1/p")
  fi

  # Clean strain: remove spaces and special characters, strip "chromosome"
  strain_clean=$(echo "$strain" | tr -d ' ' | tr -cd 'A-Za-z0-9-' | sed 's/[Cc]hromosome//g')

  # Construct new file name
  if [ -n "$species" ] && [ -n "$strain_clean" ]; then
    newname="${species}${strain_clean}.fna"
    dest="$output_dir/$newname"
    if [ ! -f "$dest" ]; then
      cp "$file" "$dest"
      echo "Copied: $filename → $newname"
    else
      echo "Skipping: $filename → $newname already exists"
    fi
  else
    echo "Could not parse header from: $filename"
  fi
done

