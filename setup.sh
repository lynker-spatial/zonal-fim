#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Name of the Conda environment
ENV_NAME="coastal_fim_vis"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Conda and try again."
    exit 1
fi

echo "Setting up the Conda environment..."

# Create or update the Conda environment from the environment.yml file
if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Updating it..."
    conda env update --name "$ENV_NAME" --file environment.yml --prune
else
    echo "Creating new environment '$ENV_NAME'..."
    conda env create --file environment.yml
fi

# Activate the Conda environment
echo "Activating the environment..."
conda activate "$ENV_NAME"

# Install the package
echo "Installing the package..."
pip install .

echo "Setup completed successfully!"
