# This subdirectory contains the preproccing steps file for the project.

## Boundary Masking

See this rebuild_mask.ipynb notebook for [masking pipeline]('preprocesing/rebuild_mask.ipynb') 

**Inputs**: Path to the masks, path to the duckdb <br>
**Process**: Merge into single mask in duckdb <br>
**Output**: generates step_5 mask that is the merged mask and masks the SCHISIM triangle elements all within duckdb <br>

## Elevation Extractions

See this extract_point_elevation.ipynb notebook for [storing elevation data]('preprocesing/extract_point_elevation.ipynb') 

**Inputs**: Path to the DEM file, path to the duckdb <br>
**Process**: Generates a table of elevation values in duckdb <br>
**Output**: generates a table of elevation values in duckdb 

## Barycentric Coordinates / Weights


## Zonal Weights

