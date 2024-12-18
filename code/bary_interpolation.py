# bary_interpolation.py

import duckdb
import ibis
from ibis import _
import rasterio
import pandas as pd
import numpy as np

def interpolate(database_path: str, s3_path: str) -> None:
    data_conn = ibis.duckdb.connect(database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
    except: 
        data_conn.raw_sql('INSTALL spatial')
        data_conn.raw_sql('LOAD spatial')

    merged_polys = data_conn.table("triangle_barycentric").select(['pg_id', 'wse_weighted_average'])
    z_w = data_conn.table("z_w").select(['pg_id', 'cell', 'coverage_fraction'])      
    
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
    z_w_merged.execute()

    # Load DEM 
    with rasterio.open("data/DEM_masked_4326.tif") as src:
        raster_meta = src.meta
        total_pixels = raster_meta['width'] * raster_meta['height']
        width = raster_meta['width']
        height = raster_meta['height']
    
    df = pd.DataFrame({'cell': range(1, total_pixels + 1), 'wse_cell_weighted_average': 0})
    data_conn.register(df, 'cell_range_table') 
    cell_range_table = data_conn.table('cell_range_table')

    # Assuming z_w_merged is already an Ibis table
    z_w_merged_with_missing = cell_range_table.left_join(
        z_w_merged,
        cell_range_table['cell'] == z_w_merged['cell']
    ).mutate(
        wse_cell_weighted_average=ibis.coalesce(
            z_w_merged['wse_cell_weighted_average'],
            cell_range_table['wse_cell_weighted_average']
        )
    ).select(
        'cell',  
        'wse_cell_weighted_average'  
    )
    z_w_merged_with_missing = z_w_merged_with_missing.order_by('cell')

    # Write as a raster file
    df_complete = z_w_merged_with_missing.execute()
    # Extract the column with cell indices and the values
    cell_indices = df_complete.index.to_numpy()  
    values = df_complete['wse_cell_weighted_average'].to_numpy()

    row_indices = ((cell_indices - 1) // width).astype(int) 
    col_indices = ((cell_indices - 1) % width).astype(int)   

    raster_array = np.zeros((height, width))
    raster_array[row_indices, col_indices] = values

    # Read metadata from the existing GeoTIFF file (DEM file)
    with rasterio.open(s3_path) as src:
        raster_meta = src.meta.copy()

    raster_meta.update({
        'dtype': 'float32',  
        'count': 1,          
        'compress': 'deflate',   
    })

    output_path = "data/barycentric_interpolation.tif"
    with rasterio.open(output_path, 'w', **raster_meta) as dst:
        dst.write(raster_array.astype('float32'), 1)
    data_conn.con.close()
    return