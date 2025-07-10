# bary_interpolation.py

import duckdb
import ibis
from ibis import _
import rasterio
from rasterio.crs import CRS
import numpy as np
import zarr
from affine import Affine
import json
from numcodecs import Blosc, Zlib
# from osgeo import gdal, osr

def interpolate(database_path: str) -> None:
    """
    The `interpolate` function processes spatial data to generate a raster file representing 
    water surface elevation (WSE) using barycentric interpolation. It produces 
    a GeoTIFF file as output.

    Input:
        - database_path (str): 
            The file path to the DuckDB database containing required spatial data tables:
              * `triangle_barycentric`: Includes `pg_id` and `wse_weighted_average`.
              * `z_w`: Maps cells to polygons with `pg_id`, `cell`, and `masked_coverage_fraction`.

        - dem_path (str): 
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
    z_w = data_conn.table("masked_coverage_fraction").select(['pg_id', 'cell', 'coverage_fraction', 'elevation'])      
    
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

def make_wse_depth_rasters(database_path: str, generate_wse: bool=False, generate_depth: bool=True,
                           output_depth_path: str='data/output_depth.tif',
                            output_wse_path: str =  "data/output_wse_barycentric_interpolation.tif", 
                             zarr_format: bool = True, depth_threshold: float = None) -> None:
    """
    The `make_wse_depth_rasters` function generates depth and/or water surface elevation (WSE) rasters 
    from a spatial database and saves them in GeoTIFF or Zarr format.

    Input:
        - database_path (str): 
            File path to the spatial database containing DEM, WSE, and depth data.
        
        - generate_wse (bool, optional):
            Whether to generate a WSE raster. Default is False.
        
        - generate_depth (bool, optional):
            Whether to generate a depth raster. Default is True.
        
        - output_depth_path (str, optional):
            File path for the output depth raster. Default is 'data/output_depth.tif'.
        
        - output_wse_path (str, optional):
            File path for the output WSE raster. Default is 'data/output_wse_barycentric_interpolation.tif'.
        
        - zarr_format (bool, optional):
            Whether to save the output in Zarr format instead of GeoTIFF. Default is True.
        
        - depth_threshold (int, optional):
            Performs a depth mask based on specific threshold i.e., 0.5 feet or 0.1524 m. Default is None.
    
    Output:
        - None:
            Generates one or both of the specified raster files, storing depth or WSE values.

    Example:
        Call the function:
        make_wse_depth_rasters("data/spatial_database.duckdb", generate_depth=True, generate_wse=True)
        
        Result:
        - A depth raster saved as 'data/output_depth.tif' or '.zarr'.
        - A WSE raster saved as 'data/output_wse_barycentric_interpolation.tif' or '.zarr'.
    
    Notes:
        - The function extracts raster metadata (CRS, height, width, etc.) from the database.
        - Outputs are stored in either GeoTIFF (compressed) or Zarr format.
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
        table = data_conn.table('dem_metadata')
        query = table.filter(table['dem_metadata'] == 'DEM')
        dem_meta = query.execute()
        # Reconstruct metadata
        crs_str = dem_meta['crs'][0]
        height = int(dem_meta['height'][0])
        width = int(dem_meta['width'][0])
        nodata = float(dem_meta['nodata'][0])
        transform_str = dem_meta['transform'][0]
        transform_values = json.loads(transform_str)
        transform = Affine(*transform_values)
        raster_meta = {
            'driver': 'GTiff',  
            'dtype': 'float32',  
            'nodata': nodata, 
            'width': width,
            'height': height,
            'count': 1,       
            'crs': CRS.from_string(crs_str),   
            'transform': transform  
        }
    
        # Create an empty raster array filled with NaN
        nodata_value = np.nan
        raster_array = np.full((height, width), nodata_value, dtype=np.float32)

        # Map cell IDs to raster indices and fill the raster array
        df['row_idx'] = ((df['cell'] - 1) // width).astype(int)
        df['col_idx'] = ((df['cell'] - 1) % width).astype(int)
        raster_array[df['row_idx'].values, df['col_idx'].values] = df[target].values
        # Mask depth raster
        if i == 0 and depth_threshold is not None:
            print(f"Applying depth threshold: filtering for depth >= {depth_threshold}")
            raster_array[raster_array < depth_threshold] = np.nan # 0.5 foot or 0.1524 m 
        # gdal_wirte = True
        # if gdal_wirte:
        #     driver = gdal.GetDriverByName("GTiff")
        #     dataset = driver.Create(
        #         output_raster_path,
        #         raster_meta.get('width'),
        #         raster_meta.get('height'),
        #         raster_meta.get('count'),  
        #         raster_meta.get('dtype'),
        #         options=["COMPRESS=DEFLATE"]  # Enable Deflate compression
        #     )
        #     dataset.SetGeoTransform(raster_meta.get('transform'))
        #     # Set the projection (CRS)
        #     srs = osr.SpatialReference()
        #     srs.ImportFromEPSG(int(crs_str.split(":")[1])) 
        #     dataset.SetProjection(srs.ExportToWkt())
        #     # Write the array to the first band
        #     band = dataset.GetRasterBand(1)
        #     band.WriteArray(raster_array)

        #     # Set no data value
        #     band.SetNoDataValue(raster_meta.get('nodata'))
        #     band.FlushCache()
        #     dataset = None
        # else:
        # Update raster metadata for nodata value
        raster_meta.update(dtype='float32', nodata=nodata_value, compress='deflate')
        if zarr_format:
            compressor = Zlib(level=9) # For space-saving, use: compressor = Blosc(cname='zstd', clevel=9)
            
            output_raster_path = output_raster_path.replace('.tif', '.zarr')
            z = zarr.open(output_raster_path, mode='w', shape=raster_array.shape, compressor=compressor,
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