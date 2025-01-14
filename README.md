# Coastal-FIM

A Python package for implementing barycentric interpolation using DuckDB, specifically designed to compute SCHISM-derived depths on 30m grids. This tool provides an efficient, scalable solution for geospatial computations in large coastal domains.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Preprocessing Workflow](#preprocessing-workflow)
- [Testing](#testing)
- [Current Progress](#current-progress)
- [Report](#report)

---

## Features

1. **Efficient Barycentric Interpolation**: Leverages DuckDB for handling large-scale geospatial data efficiently.
2. **Leverages Preprocessed Data**: Preprocessed data from DEM, zonal coverage fraction, barycentric weight, and geosptial masks allows fast and effective interpolations.
3. **Efficient Storage**: Data is converted and stored in DuckDB that significantly reduces storage volume.
4. **Customizable Pipelines**: Modular structure allows easy adaptation to different datasets and use cases.

---

## Installation
This software requires a working anaconda/miniconda installation, please visit [miniconda page](https://docs.anaconda.com/miniconda/install/)

### Steps

1. Clone the repository:
   ```bash
   cd your_custom_path
   git clone https://github.com/owp-spatial/zonal-fim.git
   cd zonal-fim
   ```

2. Then the package can be installed via executing
    ```shell
    ./setup.sh 
    ```

3. In the case there are permission issues execute
    ```shell
    chmod u+x setup.sh
    ```

### Environment

To activate pre-configured environment execute
```shell
conda activate coastal_fim_vis
```

---

## Usage

1. **Preprocessing and Prepare Input Data**:
   - visit [preprocessing folder](preprocesing/README.md) for instructions

2. **Run the Barycentric Computation**:
   Check out notebooks

3. **Output**:
   - Barycentric weights and associated data saved in DuckDB tables and parquet files.
   - Can write WSE interpolation and depth values as .tif

---

## Preprocessing Workflow

1. **Setup Databases**:
   - Convert files to DuckDB 

2. **Data Validation**:
   - Ensure CRS consistency across input datasets.
   - Validate node and element data for missing values.

3. **Pre-Processing**:
   - Use the provided tools to recompute any pt pre-processing steps as needed

4. **Post-Processing**:
   - Perform barycentric interpolation

---


## Testing

### Test Cases

1. **Tampa Region**:
   - Executed entire process: **Pass**
2. **Atlantic and Gulf Domain**:
   - Executed all steps except coverage fraction interpolation: **Pass**
3. **Comparison with Linear Interpolation**:
   - Validated results against linear interpolation: **Pass**

---

## Current Progress

1. **Implemented**:
   - Barycentric computations.
   - Batch processing for DEM and zonal data.
2. **Next Steps**:
   - Re-index cell IDs globally to match the original DEM.
   - Write comprehensive tests for the package.

---

## Report

Detailed documentation and implementation notes are available in the report:
[Report Link](https://docs.google.com/document/d/1DoPeE0IRVHkjqabqTUaX5aWCnPZn9Mdv/edit?usp=sharing&ouid=110666552849114372265&rtpof=true&sd=true)

---