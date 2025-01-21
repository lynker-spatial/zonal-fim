# This subdirectory contains the preproccing steps file for the project.

This step requires:
1. the ElementPolygons.gpkg file generated from preprocessing step 1
2. A TB DEM grid

## 1. Write the mask file to .gpkg 
By executing the `generate_zonal_req.ipynb` we add the mask layer to the `ElementPolygons.gpkg`

## 2. Generate coverage-fraction and elevation values for all cells 
The `weights_mask_elev_table.R` file then does the following

1. Takes the triangles and TB dataset to compute the relative coverage fraction of each cell and polygon ID. 
2. For each cell, the cell the area is determined as mask or not, and an elevation is assigned

The result is a table that contains the following columns:

```r
> head(fin)
    pg_id      cell coverage_fraction elevation  mask
    <num>     <num>             <num>     <num> <num>
1: 122881 720766867       0.023299048  5.497267    NA
2: 122881 720801591       0.552611947  5.367274    NA
3: 122881 720801592       0.040903941  5.823075    NA
4: 122881 720836314       0.006548808  4.383194    NA
5: 122881 720836315       0.748685956  4.613204    NA
6: 122881 720836316       0.334290981  4.148118    NA
```