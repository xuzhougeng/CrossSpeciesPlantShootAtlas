#!/bin/bash
#
# DeepTFactor prediction script
# Usage: Run this script from the deeptfactor directory
#   ./TF_prediction.sh input.fasta [gpu|cpu]
# Source: https://bitbucket.org/kaistsystemsbiology/deeptfactor/src/master/
#

# Check if running from deeptfactor directory
if [ ! -f "tf_running.py" ]; then
  echo "Error: This script must be run from the deeptfactor directory"
  echo "Please cd to the deeptfactor directory and try again"
  exit 1
fi

# Check for input file parameter
if [ -z "$1" ]; then
  echo "Error: Missing input fasta file"
  echo "Usage: $0 input.fasta [gpu|cpu]"
  exit 1
fi

# Set GPU/CPU option (default to CPU if not specified)
GPU_OPT="cpu"
if [ "$2" = "gpu" ]; then
  GPU_OPT="cuda:0"
fi

# Run prediction with specified parameters
python tf_running.py -i "$1" -o ./result -g $GPU_OPT