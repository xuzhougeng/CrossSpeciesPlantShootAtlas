#!/bin/bash

# Directory containing protein files
prot_dir="processed"
n_threads=40  # Define the number of threads here

# Create an array to hold the base names of the protein files
declare -a names

# Read all pep files from the processed directory
index=0
for file in ${prot_dir}/*.pep; do
  base_name=$(basename $file .pep)  # Extract the base name without extension
  names[$index]=$base_name
  let "index++"
  makeblastdb -in $file -dbtype prot  # Create a BLAST database for each file
done

# Function to run blastp for two given sequences
run_blastp() {
  local query=$1
  local db=$2
  mkdir -p maps/${query}${db}
  blastp -query ${prot_dir}/${query}.pep -db ${prot_dir}/${db}.pep -outfmt 6 -out maps/${query}${db}/${query}_to_${db}.txt -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6 &
  blastp -query ${prot_dir}/${db}.pep -db ${prot_dir}/${query}.pep -outfmt 6 -out maps/${query}${db}/${db}_to_${query}.txt -num_threads ${n_threads} -max_hsps 1 -evalue 1e-6 &
}

# Loop through all pairs of names for the two-way BLASTp
for (( i=0; i<${#names[@]}; i++ )); do
  for (( j=i+1; j<${#names[@]}; j++ )); do
    run_blastp ${names[i]} ${names[j]}
  done
done

# Wait for all BLAST processes to finish
wait