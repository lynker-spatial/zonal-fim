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
   The script zonal_fim.py performs execution of pipeline for generating barycentric interpolation, generating masks, and preprocessing pipeline
1. **Preprocessing and Prepare Input Data**:
   1. The first step is to generate an overall mask using a set of different masks and the following procedure
   
      This can be achieved by executing zonal_fim.py and activating generate_mask flag while other flags are set to `False` also one can alter the name of masks shapefile names in case they have changed by these flags:

      --water_table_name <br>
      --nwm_table_name <br>
      --levee_table_name <br>
      --state_table_name <br>
      --schisim_table_name <br>

      Several paths must be specified including:

      -k '/path/shape_file_folder' that is the folder path where the schisim element shape file lives <br>
      -l '/path/output_elements_folder' that is the folder path where the schisim element parquet and geopackage will be saved to <br>
      -a '/path/masks.duckdb' path to the mask databse that will store the generated overall mask <br>

      and there is an optional item to dissolve all geometries into a single multipolygon feature that can be set but recommend the default `True` value

      --dissolve

      ```shell
      python zonal_fim.py --generate_mask True --preprocess False --generate_wse False --generate_depth False --zarr_format False  --execute False  --dissolve True -k '/path/shape_file_folder' -l '/path/output_elements_folder' -a '/path/masks.duckdb' -p 'exterior_mask_schism_boundary_atlantic_buffer_atlgulf' -n 'exterior_mask_state_boundaries_conus_atlgulf' -x 'interior_mask_levee_protected_area_conus_atlgulf' -t 'interior_mask_nwm_lakes_conus_atlgulf' -f 'interior_mask_water_polygon_conus_atlgulf'
      ```
     
   2. Next step is to generate coverage fractions from zonal in R
      - visit [preprocessing folder](preprocesing/README.md) for instructions


   3. Final step in preprocessing is to generate barycentric weights for all none-masked schisim node. This can be done by calling  the script with the following configuration.

      Setting all flags to `False` except for <br>
      --preprocess True <br>
     
      ```shell
      python zonal_fim.py --generate_mask False --preprocess True --generate_wse False --generate_depth False --zarr_format False  --execute False -u '/path/ElementPolygons.parquet' -o '/path/TBDEM_AtlanticGulf_Mosaic_NWM_3_Revised_v4_COG.tif' -i '/path/agGridfile.gr3' -c '/path/schisim_database.duckdb' -w '/path/coverage_fraction.parquet'
      ```
   
2. **Run the Barycentric Computation**:
    General pipeline for executing barycentric interpolation is given that a preprocessing has been done once, we can pass a new .gr3 file and specify the output path. There is an option `--zarr_format` to produce the outputs as zarr instead of a .tif file one does not need to change .tif to .zarr in dpeth and wse inputs this conversion will be done automatically.    
    
    `Note:` make sure that the DEM is in 4326 if not performing preprocessing!

    ```shell
    python zonal_fim.py --generate_mask False --preprocess False --generate_wse False --generate_depth True --zarr_format False  --execute True  --dissolve False -i '/path/agGridfile.gr3' -c '/path/schisim_database.duckdb' -m '/path/depth_raser_v1.tif' -q '/path/wse_raser_v1.tif'
    ```

3. **Output**:
   - Barycentric interpolation is saved as depth table in the DuckDB database.
   - Can write WSE interpolation and depth values as .tif and .zarr file if specified. 

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

