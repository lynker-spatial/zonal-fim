# bary_interpolation.py

import duckdb
import ibis
from ibis import _
import rasterio
import numpy as np
import zarr

def interpolate(database_path: str) -> None:
    """
    The `interpolate` function processes spatial data to generate a raster file representing 
    water surface elevation (WSE) using barycentric interpolation. It produces 
    a GeoTIFF file as output.

    Input:
        - database_path (str): 
            The file path to the DuckDB database containing required spatial data tables:
              * `triangle_barycentric`: Includes `pg_id` and `wse_weighted_average`.
              * `z_w`: Maps cells to polygons with `pg_id`, `cell`, and `coverage_fraction`.

        - s3_path (str): 
            The path to a reference GeoTIFF file on the local file system or an S3-compatible 
            storage. This file provides raster dimensions and metadata.

    Output:
        - None: 
            Generates a GeoTIFF raster file (`wse_barycentric_interpolation.tif`) containing 
            interpolated WSE values.

    Example:
        1. Ensure the following data exists:
           - DuckDB database (`my_database.duckdb`) with `triangle_barycentric` and `z_w` tables.
           - Reference raster file (`DEM_masked_4326.tif`).

        2. Call the function:
            interpolate("my_database.duckdb", "data/DEM_masked_4326.tif")

        3. Result:
            A GeoTIFF file named `wse_barycentric_interpolation.tif` is created in the `data/` 
            directory, containing interpolated WSE values.

    Notes:
        - The output file uses `float32` data type for raster values.

    """
    data_conn = ibis.duckdb.connect(database_path)
    # out_data_conn = ibis.duckdb.connect(output_database_path)
    for conn in [data_conn]:
        try:
            conn.raw_sql('LOAD spatial')
        except:
            conn.raw_sql('INSTALL spatial')
            conn.raw_sql('LOAD spatial')

    # out_data_conn.raw_sql(f"ATTACH '{database_path}' AS compute_db;")

    merged_polys = data_conn.table("triangle_barycentric").select(['pg_id', 'wse_weighted_average'])
    z_w = data_conn.table("coverage_fraction").select(['pg_id', 'cell', 'coverage_fraction', 'elevation'])      
    
    # left join on the 'pg_id' column
    z_w_merged = z_w.join(merged_polys, z_w.pg_id == merged_polys.pg_id, how='left')
    z_w_merged = z_w_merged.drop(z_w_merged.pg_id_right)
    grouped = z_w_merged.group_by("cell")

    # Calculate 'wse_cell_weighted_average'
    numerator = grouped.aggregate(
        weighted_sum = (z_w_merged.wse_weighted_average * z_w_merged.coverage_fraction).sum()
    )
    denominator = grouped.aggregate(
        coverage_sum = z_w_merged.coverage_fraction.sum()
    )
    z_w_merged = z_w_merged.left_join(numerator, z_w_merged.cell == numerator.cell)
    z_w_merged = z_w_merged.drop(z_w_merged.cell_right)
    z_w_merged = z_w_merged.left_join(denominator, z_w_merged.cell == denominator.cell)
    z_w_merged = z_w_merged.drop(z_w_merged.cell_right)
    z_w_merged = z_w_merged.mutate(
        wse_cell_weighted_average=z_w_merged.weighted_sum / z_w_merged.coverage_sum
    )

    z_w_merged = z_w_merged.group_by('cell').aggregate(
        wse_cell_weighted_average=z_w_merged.wse_cell_weighted_average.first()
    )

    z_w_cell = (
        z_w
        .group_by("cell")  
        .aggregate(*[z_w[col].first().name(col) for col in z_w.columns if col != "cell"])  # Aggregate all columns except 'cell'
    )

    z_w_merged = z_w_merged.join(z_w_cell.select(['cell', 'elevation']), z_w_merged.cell == z_w_cell.cell, how='left')
    z_w_merged = z_w_merged.drop(z_w_merged.cell_right)
    z_w_merged = z_w_merged.filter(~z_w_merged["elevation"].isnan())

    z_w_with_depth = z_w_merged.mutate(depth=z_w_merged["wse_cell_weighted_average"] - z_w_merged["elevation"])

    # Filter the rows where 'depth' is greater than or equal to 0
    filtered_z_w = z_w_with_depth.filter(z_w_with_depth["depth"] >= 0)
    data_conn.create_table("depth", filtered_z_w,  overwrite=True)

    data_conn.con.close()
    # out_data_conn.con.close()
    return

def make_wse_depth_rasters(database_path: str, s3_path: str, generate_wse: bool=False, generate_depth: bool=True,
                           output_depth_path: str='data/output_depth.tif',
                            output_wse_path: str =  "data/output_wse_barycentric_interpolation.tif", 
                             zarr_format: bool = True, mask_negative: bool = True) -> None:
    """
    The `make_depth_raster` function calculates a depth raster by computing the difference between 
    a water surface elevation (WSE) raster and a digital elevation model (DEM) raster. The result 
    is saved as a new GeoTIFF file.

    Input:
        - dem_path (str): 
            File path to the DEM raster file. This file represents the elevation of the terrain.

        - wse_path (str, optional): 
            File path to the WSE raster file. The default value is 
            "data/wse_barycentric_interpolation.tif". This raster represents water surface elevation.

        - mask_negative (bool, optional): 
            If True, negative values in the resulting depth raster are masked (set to `NaN`). 
            Default is True.

    Output:
        - None: 
            Generates a GeoTIFF raster file (`depth_barycentric_interpolation.tif`) in the `data/` 
            directory, representing the depth values.

    Example:
        1. Ensure the following input files are available:
           - A DEM raster file, e.g., `data/DEM_masked_4326.tif`.
           - A WSE raster file, e.g., `data/wse_barycentric_interpolation.tif`.

        2. Call the function:
            make_depth_raster("data/DEM_masked_4326.tif", "data/wse_barycentric_interpolation.tif")

        3. Result:
            A GeoTIFF file named `depth_barycentric_interpolation.tif` is created in the `data/` 
            directory, containing the depth values.

    Notes:
        - Both input rasters must have the same Coordinate Reference System (CRS). 
          If they differ, the function raises a `ValueError`.
        - NoData values in either raster are preserved in the output raster as `NaN`.
    """
    data_conn = ibis.duckdb.connect(database_path)
    # out_data_conn = ibis.duckdb.connect(output_database_path)
    for conn in [data_conn]:
        try:
            conn.raw_sql('LOAD spatial')
        except:
            conn.raw_sql('INSTALL spatial')
            conn.raw_sql('LOAD spatial')
    flags = [generate_depth, generate_wse]
    # flag_count = sum(flags)
    for i, flag in enumerate(flags):
        if not flag:
            continue
        if i == 0:
            target = 'depth'
            output_raster_path = output_depth_path
        if i == 1:
            target = 'wse_cell_weighted_average'
            output_raster_path = output_wse_path

        # Load DEM and depth
        df = data_conn.table('depth').execute()
        # Extract raster metadata
        with rasterio.open(s3_path) as src:
            raster_meta = src.meta
            width = raster_meta['width']
            height = raster_meta['height']
        
        # Create an empty raster array filled with NaN
        nodata_value = np.nan
        raster_array = np.full((height, width), nodata_value, dtype=np.float32)

        # Map cell IDs to raster indices and fill the raster array
        for _, row in df.iterrows():
            cell_id = int(row['cell'])
            depth_value = row[target]
            
            row_idx = int((cell_id - 1) // width) 
            col_idx = int((cell_id - 1) % width)  

            raster_array[row_idx, col_idx] = depth_value

        # Update raster metadata for nodata value
        raster_meta.update(dtype='float32', nodata=nodata_value)
        if zarr_format:
            output_raster_path = output_raster_path.replace('.tif', '.zarr')
            z = zarr.open(output_raster_path, mode='w', shape=raster_array.shape, 
                        dtype='float32', chunks=(1000, 1000), fill_value=nodata_value)
            z[:, :] = raster_array
        else:
            # Write the raster array to a GeoTIFF file
            with rasterio.open(output_raster_path, 'w', **raster_meta) as dst:
                dst.write(raster_array, 1) 
        print(f"Saved file to: {output_raster_path}") 

    data_conn.con.close()
    # out_data_conn.con.close()
    return