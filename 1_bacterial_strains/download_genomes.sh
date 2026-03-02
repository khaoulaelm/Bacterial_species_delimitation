#!/bin/bash

#Downloads genomes from NCBI Assembly
#Input: strains.txt (47 VUB strain names)
#Output: <strain_name>.fna

# Directory to save the downloaded genomes
mkdir -p genomes

# Read strains from the comma-separated list in strains.txt
IFS=',' read -r -a strains < strains.txt

# Loop through each strain
for strain in "${strains[@]}"; do
    # Remove leading and trailing spaces from the strain name
    strain=$(echo "$strain" | xargs)
    
    # Skip empty strains
    if [[ -z "$strain" ]]; then
        continue
    fi
    
    # Print the current strain being processed
    echo "Processing strain: $strain"
    
    # Query NCBI for the assembly ID of the strain
    assembly_id=$(esearch -db assembly -query "$strain" | esummary | xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)
    
    # If assembly ID is empty, skip this strain
    if [[ -z "$assembly_id" ]]; then
        echo "No assembly found for $strain"
        continue
    fi
    
    # Get the FTP path of the genome from NCBI
    ftp_path=$(esummary -db assembly -id "$assembly_id" | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
    
    # If no FTP path is found, try GenBank FTP path
    if [[ -z "$ftp_path" ]]; then
        ftp_path=$(esummary -db assembly -id "$assembly_id" | xtract -pattern DocumentSummary -element FtpPath_GenBank)
    fi
    
    # If no FTP path is found, skip this strain
    if [[ -z "$ftp_path" ]]; then
        echo "No FTP path found for $strain"
        continue
    fi
    
    # Convert FTP to HTTPS
    https_path=${ftp_path//ftp:/https:}
    
    # Construct the URL for the genomic file (assumed .fna.gz format)
    genomic_url="${https_path}/$(basename "$ftp_path")_genomic.fna.gz"
    
    # Debug: Show the URL
    echo "Downloading from: $genomic_url"
    
    # Download the genome to the genomes directory
    wget -O "genomes/$strain.fasta.gz" "$genomic_url"
    
    # If the download was successful, unzip the file
    if [[ $? -eq 0 ]]; then
        # Remove existing .fasta file if it exists
        if [[ -f "genomes/$strain.fasta" ]]; then
            rm "genomes/$strain.fasta"
        fi
        
        # Attempt to unzip the downloaded .fasta.gz file
        gunzip "genomes/$strain.fasta.gz"
        
        # Check if gunzip was successful
        if [[ $? -eq 0 ]]; then
            echo "Downloaded and unzipped genome for $strain"
        else
            echo "Error unzipping $strain.fasta.gz"
        fi
    else
        echo "Failed to download genome for $strain"
    fi
done

echo "All genome downloads are complete!"
