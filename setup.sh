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

# Ensure Conda is initialized
if ! grep -q "conda initialize" ~/.bashrc ~/.zshrc ~/.profile &> /dev/null; then
    echo "Initializing Conda for your shell..."
    conda init
    echo "Conda initialized. Please restart your shell or source your shell configuration file, then re-run this script."
    exit 1
fi

# Source the Conda profile script to enable Conda commands
source "$(conda info --base)/etc/profile.d/conda.sh"

# Check if the environment exists
if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Updating it..."
    conda env update --name "$ENV_NAME" --file environment.yml --prune
else
    echo "Creating new environment '$ENV_NAME'..."
    conda env create --file environment.yml
fi

# Activate the Conda environment
echo "Activating the Conda environment..."
conda activate "$ENV_NAME"

# Install the package in the activated environment
echo "Installing the package..."
pip install .

echo "Setup completed successfully!"