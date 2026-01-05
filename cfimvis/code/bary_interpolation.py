# bary_interpolation.py
from typing import Dict, Optional
import ibis
from ibis import _
import rasterio
from rasterio.crs import CRS
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import zarr
from affine import Affine
import json
import gc
import os
from numcodecs import Zlib
from osgeo import gdal, osr
import uuid

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS")
gdal.SetCacheMax(1024 * 1024 * 512)  # 512 MB cache

def sparse_tile_worker(args):
    """
    The `sparse_tile_worker` function processes a DataFrame chunk for a specific tile, 
    creating a dense NumPy array that represents that tile's data.

    Input:
        - args (tuple): 
            A tuple containing the parameters for tile processing. It includes:
            - x_offset (int): The starting column offset of the tile.
            - y_offset (int): The starting row offset of the tile.
            - block_size (int): The height and width of the tile.
            - nodata_value (float): The value to fill in empty cells.
            - block_df (pd.DataFrame): DataFrame containing sparse data ('row', 'col', target_column) for the tile.
            - target_column (str): The name of the column with the raster values.
    
    Output:
        - tuple or None:
            - On success, returns a tuple containing (x_offset, y_offset, sparse_tile_array).
            - Returns None if an exception occurs during tile creation.
    """
    x_offset, y_offset, block_size, nodata_value, block_df, target_column = args
    try:
        sparse_tile_array = np.full((block_size, block_size), nodata_value, dtype=np.float32)
        relative_rows = (block_df['row'] - y_offset).values
        relative_cols = (block_df['col'] - x_offset).values
        values = block_df[target_column].values.astype(np.float32)

        mask_in_bounds = (relative_rows >= 0) & (relative_rows < block_size) & (relative_cols >= 0) & (relative_cols < block_size)
        if not np.any(mask_in_bounds):
            return (x_offset, y_offset, sparse_tile_array)

        sparse_tile_array[relative_rows[mask_in_bounds], relative_cols[mask_in_bounds]] = values[mask_in_bounds]
        return (x_offset, y_offset, sparse_tile_array)
    except Exception as e:
        print(f"Error preparing tile at ({x_offset}, {y_offset}): {e}")
        return None
    
def make_wse_depth_cogs(data_table: ibis.expr.types.Table, dem_meta: pd.DataFrame, target_column: str,
                        output_cog_path: str, depth_threshold: float = None) -> None:
    """
    The `make_wse_depth_cogs` function creates a Cloud Optimized GeoTIFF (COG) from an Ibis table expression 
    and saves it to a specified path.

    Input:
        - data_table (ibis.expr.types.Table): 
            An Ibis table expression containing 'cell' and the target data column (e.g., 'depth' or 'wse').
        
        - dem_meta (pd.DataFrame):
            A DataFrame containing raster metadata such as 'height', 'width', 'transform', and 'crs'.
        
        - target_column (str):
            The name of the column in `data_table` to extract raster values from.
        
        - output_cog_path (str):
            File path for the output COG raster.
        
        - depth_threshold (float, optional):
            Performs a mask on the 'depth' column, setting values below this threshold to NoData. Default is None.
    
    Output:
        - None:
            Generates a Cloud Optimized GeoTIFF (COG) file, storing the specified raster values.

    Example:
        Call the function:
        make_wse_depth_cogs(
            data_table=my_ibis_table, 
            dem_meta=my_dem_metadata, 
            target_column='depth', 
            output_cog_path='data/output_depth.tif', 
            depth_threshold=0.5
        )
        
        Result:
        - A COG raster saved as 'data/output_depth.tif', where depth values below 0.5 are masked.
    """
    height, width = int(dem_meta['height'][0]), int(dem_meta['width'][0])
    transform_list = json.loads(dem_meta['transform'][0])
    gdal_transform = (transform_list[2], transform_list[0], transform_list[1], 
                      transform_list[5], transform_list[3], transform_list[4])
    crs_string_from_db = dem_meta['crs'][0]

    nodata_value = -9999.0 #np.nan
    block_size = 512
    temp_tif_path = output_cog_path.replace(".tif", "_temp.tif")

    driver = gdal.GetDriverByName("GTiff")
    creation_options = [
        "TILED=YES", f"BLOCKXSIZE={block_size}", f"BLOCKYSIZE={block_size}",
        "COMPRESS=DEFLATE", "BIGTIFF=YES", "SPARSE_OK=TRUE"
    ]

    dst_ds = driver.Create(temp_tif_path, width, height, 1, gdal.GDT_Float32, creation_options)
    dst_ds.SetGeoTransform(gdal_transform)
    
    srs = osr.SpatialReference()
    srs.SetFromUserInput(crs_string_from_db)
    dst_ds.SetProjection(srs.ExportToWkt())
    
    band = dst_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata_value)

    max_workers = os.cpu_count() or 1
    
    data_table = data_table.order_by('cell')
    total_rows = data_table.count().execute()
    chunk_size = 5000000

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for offset in range(0, total_rows, chunk_size):
            print(f"  Processing rows {offset:,} to {min(offset + chunk_size, total_rows):,} / {total_rows:,}")
            chunk_df = data_table.select('cell', target_column).limit(chunk_size, offset=offset).execute()
            if chunk_df.empty:
                continue

            if target_column == 'depth' and depth_threshold is not None:
                chunk_df.loc[chunk_df[target_column] < depth_threshold, target_column] = nodata_value

            chunk_df['row'] = ((chunk_df['cell'] - 1) // width).astype(np.int64)
            chunk_df['col'] = ((chunk_df['cell'] - 1) % width).astype(np.int64)
            chunk_df['block_row_off'] = (chunk_df['row'] // block_size) * block_size
            chunk_df['block_col_off'] = (chunk_df['col'] // block_size) * block_size

            grouped = chunk_df.groupby(['block_row_off', 'block_col_off'])
            tasks = [
                (int(col_off), int(row_off), block_size, nodata_value, block_df_grouped[['row','col',target_column]].copy(), target_column)
                for (row_off, col_off), block_df_grouped in grouped
            ]
            del chunk_df, grouped
            gc.collect()

            for result in executor.map(sparse_tile_worker, tasks):
                if not result:
                    continue
                x_offset, y_offset, sparse_tile = result 
                existing_tile = band.ReadAsArray(x_offset, y_offset, block_size, block_size)
                if existing_tile is None:
                    continue
                existing_tile = np.array(existing_tile, dtype=np.float32)

                h_existing, w_existing = existing_tile.shape
                h_sparse, w_sparse = sparse_tile.shape
                h = min(h_existing, h_sparse)
                w = min(w_existing, w_sparse)

                sparse_crop = sparse_tile[:h, :w]
                existing_crop = existing_tile[:h, :w]

                valid_mask = sparse_crop != nodata_value
                existing_crop[valid_mask] = sparse_crop[valid_mask]
                existing_tile[:h, :w] = existing_crop
                band.WriteArray(existing_tile[:h, :w], xoff=x_offset, yoff=y_offset)

    band.FlushCache()
    dst_ds = None

    print("Creating COG file.")
    gdal.Translate(
        output_cog_path,
        temp_tif_path,
        format="COG",
        creationOptions=[
            "COMPRESS=DEFLATE",
            "PREDICTOR=2",
            "OVERVIEWS=AUTO",
            "RESAMPLING=AVERAGE",
            "NUM_THREADS=ALL_CPUS"
        ]
    )

    try:
        os.remove(temp_tif_path)
    except OSError as e:
        print(f"Error removing temporary file: {e}")

    print(f"Output saved to: {output_cog_path}")
    return

def compute_barycentric_values(database_path: str, generate_depth: bool,
                               generate_wse: bool, depth_threshold: float = 0.0) -> Optional[str]:
    """
    Computes weighted average WSE and Depth from the database and saves the result 
    to a temporary Parquet file.
    
    Returns:
        str: Path to the temporary Parquet file containing columns ['cell', 'depth', 'wse_cell_weighted_average'].
        None: If neither flag is set or no data is generated.
    """
    if not generate_depth and not generate_wse:
        print("Both generate_depth and generate_wse are False. No output will be generated.")
        return None

    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql('LOAD spatial')
    temp_parquet_path = f"temp_results_{uuid.uuid4().hex}.parquet"

    try:
        merged_polys = data_conn.table("triangle_barycentric").select('pg_id', 'wse_weighted_average')
        z_w = data_conn.table("masked_coverage_fraction").select('pg_id', 'cell', 'coverage_fraction', 'elevation')

        print("Calculating cell wise weighted average WSE.")
        z_w_merged_expr = z_w.join(merged_polys, z_w.pg_id == merged_polys.pg_id, how='left')
        wse_avg_expr = (
            z_w_merged_expr.group_by('cell')
            .aggregate(
                wse_cell_weighted_average=(
                    (z_w_merged_expr.wse_weighted_average * z_w_merged_expr.coverage_fraction).sum() /
                    z_w_merged_expr.coverage_fraction.sum()
                )
            )
        )
        data_conn.create_table("temp_cell_wse_avg", wse_avg_expr, overwrite=True)

        print("Extracting unique elevation per cell.")
        cell_elevation_agg = (z_w.group_by("cell").aggregate(elevation=z_w.elevation.first()))
        cell_elevation_expr = cell_elevation_agg.filter(~cell_elevation_agg.elevation.isnan())
        data_conn.create_table("temp_cell_elevation", cell_elevation_expr, overwrite=True)

        wse_table = data_conn.table("temp_cell_wse_avg")
        elev_table = data_conn.table("temp_cell_elevation")

        final_table_expr = wse_table.join(
            elev_table, wse_table.cell == elev_table.cell, how='inner'
        ).mutate(depth=(ibis._.wse_cell_weighted_average - ibis._.elevation))

        if depth_threshold is None:
            depth_threshold = 0.0
        final_table_expr = final_table_expr.filter(final_table_expr.depth >= depth_threshold)

        columns_to_select = ['cell']
        if generate_wse:
            columns_to_select.append('wse_cell_weighted_average')
        if generate_depth:
            columns_to_select.append('depth')
        
        final_table_expr = final_table_expr.select(*columns_to_select)
        final_table_expr.to_parquet(temp_parquet_path)
        
        check_conn = ibis.duckdb.connect()
        count = check_conn.read_parquet(temp_parquet_path).count().execute()
        check_conn.con.close()
        
        print(f"Materialized {count:,} valid cells with depth >= {depth_threshold}.")
        return temp_parquet_path

    except Exception as e:
        print(f"Error during computation: {e}")
        if os.path.exists(temp_parquet_path):
            os.remove(temp_parquet_path)
        raise e
    finally:
        try:
            data_conn.drop_table("temp_cell_wse_avg", force=True)
            data_conn.drop_table("temp_cell_elevation", force=True)
        except Exception as e:
            print(f"Warning: Could not drop temporary tables. {e}")
        if data_conn and data_conn.con:
             data_conn.con.close()

def rasterize_barycentric_values(parquet_path: str, database_path: str, source_file_path: str,
                                 output_format: str, generate_depth: bool, generate_wse: bool,
                                 depth_path: str, wse_path: str) -> Optional[Dict[str, gdal.Dataset]]:
    """
    Reads the intermediate Parquet file and produces rasters (COG, ZARR, or IN_MEMORY).
    Uses specific depth_path and wse_path for output locations.
    """
    if not parquet_path or not os.path.exists(parquet_path):
        print("No intermediate Parquet file found. Skipping rasterization.")
        return None

    output_format = output_format.upper()
    valid_formats = ['COG', 'ZARR', 'IN_MEMORY']
    if output_format not in valid_formats:
        raise ValueError(f"output_format must be one of {valid_formats}, got '{output_format}'")

    data_conn = ibis.duckdb.connect(database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
        dem_metadata_table = data_conn.table('dem_metadata')
        dem_meta = dem_metadata_table.filter(dem_metadata_table['dem_metadata'] == 'DEM').execute()
        
        data_table_from_parquet = data_conn.read_parquet(parquet_path)
        base_name = os.path.splitext(os.path.basename(source_file_path))[0]
        gdal_datasets = {}

        if output_format == 'COG':
            if generate_wse:
                print(f"Generating WSE COG at {wse_path}")
                make_wse_depth_cogs(data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                                    dem_meta=dem_meta,
                                    target_column='wse_cell_weighted_average',
                                    output_cog_path=wse_path)
            if generate_depth:
                print(f"Generating Depth COG at {depth_path}")
                make_wse_depth_cogs(data_table=data_table_from_parquet.select("cell", "depth"),
                                    dem_meta=dem_meta,
                                    target_column='depth',
                                    output_cog_path=depth_path)
            return None

        elif output_format == 'ZARR':
            zarr_store_path = depth_path if depth_path.endswith('.zarr') else os.path.join(os.path.dirname(depth_path), "atlgulf_fim.zarr")
            
            if generate_wse:
                write_ibis_to_zarr(data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                                   dem_meta=dem_meta,
                                   target_column='wse_cell_weighted_average',
                                   output_zarr_path=zarr_store_path,
                                   array_name=f"{base_name}_wse")
            if generate_depth:
                write_ibis_to_zarr(data_table=data_table_from_parquet.select("cell", "depth"),
                                   dem_meta=dem_meta,
                                   target_column='depth',
                                   output_zarr_path=zarr_store_path,
                                   array_name=f"{base_name}_depth")
            return None

        elif output_format == 'IN_MEMORY':
            if generate_wse:
                gdal_datasets['wse'] = create_in_memory_gdal_array(
                    data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                    dem_meta=dem_meta, target_column='wse_cell_weighted_average'
                )
            if generate_depth:
                gdal_datasets['depth'] = create_in_memory_gdal_array(
                    data_table=data_table_from_parquet.select("cell", "depth"),
                    dem_meta=dem_meta, target_column='depth'
                )
            return gdal_datasets

    finally:
        if data_conn and data_conn.con:
            data_conn.con.close()
            
def interpolate_and_rasterize(database_path: str, source_file_path: str, output_format: str, generate_depth: bool,
                                generate_wse: bool, output_path: str, depth_threshold: float = 0.0) -> Optional[Dict[str, gdal.Dataset]]:
    """
    The `interpolate_and_rasterize` function orchestrates the generation of depth and/or WSE 
    raster data by performing barycentric interpolation. It can produce outputs as Cloud 
    Optimized GeoTIFFs (COG), Zarr stores, or in-memory GDAL Dataset objects.

    Input:
        - database_path (str): 
            The file path to the DuckDB database containing required spatial data tables:
              * `triangle_barycentric`: Includes `pg_id` and `wse_weighted_average`.
              * `masked_coverage_fraction`: Maps cells to polygons with `pg_id`, `cell`, etc.
              * `dem_metadata`: Contains the georeferencing information for the output raster.

        - source_file_path (str): 
            The path to a source data file (e.g., a NetCDF file) which is used to derive a 
            unique base name for the output arrays or files.

        - output_format (str):
            The desired output format, which must be one of 'COG', 'ZARR', or 'IN_MEMORY'.

        - generate_depth (bool):
            A flag to indicate whether to generate a depth raster.

        - generate_wse (bool):
            A flag to indicate whether to generate a water surface elevation (WSE) raster.

        - output_path (str):
            The destination path for the output. Its meaning varies with `output_format`:
              * For 'COG', this is the full file path for the depth raster.
              * For 'ZARR', this path is used to determine the parent directory for the Zarr store.

        - depth_threshold (float, optional):
            The minimum depth value to include in the final raster data. Defaults to 0.0.

    Output:
        - Optional[Dict[str, gdal.Dataset]]: 
            - For 'COG' and 'ZARR' formats, the function generates files on disk and returns `None`.
            - For the 'IN_MEMORY' format, it returns a dictionary mapping 'depth' and/or 'wse' 
              to their corresponding in-memory `gdal.Dataset` objects.

    Example:
        1. Ensure the following data exists:
           - DuckDB database (`spatial_data.duckdb`) with all required tables.
           - A source file for naming (`nwm.t00z.analysis.nc`).

        2. Call the function to generate a Zarr store:
            interpolate_and_rasterize(
                database_path="data/spatial_data.duckdb",
                source_file_path="data/nwm.t00z.analysis.nc",
                output_format="ZARR",
                generate_depth=True,
                generate_wse=True,
                output_path="output/rasters/"
            )

        3. Result:
            A Zarr store named `atlgulf_fim.zarr` is created in the `output/rasters/`
            directory. It will contain two new arrays, `nwm.t00z.analysis_wse` and
            `nwm.t00z.analysis_depth`, with the interpolated raster values.

    """
    output_format = output_format.upper()
    valid_formats = ['COG', 'ZARR', 'IN_MEMORY']
    if output_format not in valid_formats:
        raise ValueError(f"output_format must be one of {valid_formats}, got '{output_format}'")

    if not generate_depth and not generate_wse:
        print("Both generate_depth and generate_wse are False. No output will be generated.")
        return None

    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql('LOAD spatial')
    temp_parquet_path = f"temp_results_{uuid.uuid4().hex}.parquet"
    gdal_datasets = {}

    try:
        merged_polys = data_conn.table("triangle_barycentric").select('pg_id', 'wse_weighted_average')
        z_w = data_conn.table("masked_coverage_fraction").select('pg_id', 'cell', 'coverage_fraction', 'elevation')

        print("Calculating cell wise weighted average WSE.")
        z_w_merged_expr = z_w.join(merged_polys, z_w.pg_id == merged_polys.pg_id, how='left')
        wse_avg_expr = (
            z_w_merged_expr.group_by('cell')
            .aggregate(
                wse_cell_weighted_average=(
                    (z_w_merged_expr.wse_weighted_average * z_w_merged_expr.coverage_fraction).sum() /
                    z_w_merged_expr.coverage_fraction.sum()
                )
            )
        )
        data_conn.create_table("temp_cell_wse_avg", wse_avg_expr, overwrite=True)

        print("Extracting unique elevation per cell.")
        cell_elevation_agg = (
            z_w.group_by("cell")
            .aggregate(elevation=z_w.elevation.first())
        )
        cell_elevation_expr = cell_elevation_agg.filter(~cell_elevation_agg.elevation.isnan())
        data_conn.create_table("temp_cell_elevation", cell_elevation_expr, overwrite=True)

        wse_table = data_conn.table("temp_cell_wse_avg")
        elev_table = data_conn.table("temp_cell_elevation")

        final_table_expr = wse_table.join(
            elev_table, wse_table.cell == elev_table.cell, how='inner'
        ).mutate(
            depth=(ibis._.wse_cell_weighted_average - ibis._.elevation)
        )

        if depth_threshold is None:
            depth_threshold = 0.0
        final_table_expr = final_table_expr.filter(final_table_expr.depth >= depth_threshold)

        columns_to_select = ['cell']
        if generate_wse:
            columns_to_select.append('wse_cell_weighted_average')
        if generate_depth:
            columns_to_select.append('depth')
        
        final_table_expr = final_table_expr.select(*columns_to_select)
        final_table_expr.to_parquet(temp_parquet_path)

        data_table_from_parquet = data_conn.read_parquet(temp_parquet_path)
        
        result_count = data_table_from_parquet.count().execute()
        print(f"Materialized {result_count:,} valid data cells with depth >= {depth_threshold}.")
        if result_count == 0:
            print("No data remains after filtering. Output will be empty but valid.")

        dem_metadata_table = data_conn.table('dem_metadata')
        dem_meta = dem_metadata_table.filter(dem_metadata_table['dem_metadata'] == 'DEM').execute()
        base_name = os.path.splitext(os.path.basename(source_file_path))[0]

        if output_format == 'COG':
            if generate_wse:
                make_wse_depth_cogs(data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                                    dem_meta=dem_meta,
                                    target_column='wse_cell_weighted_average',
                                    output_cog_path=output_path.replace('depth', 'wse'))
            if generate_depth:
                make_wse_depth_cogs(data_table=data_table_from_parquet.select("cell", "depth"),
                                    dem_meta=dem_meta,
                                    target_column='depth',
                                    output_cog_path=output_path)
            return None

        elif output_format == 'ZARR':
            output_dir = os.path.dirname(output_path)
            zarr_store_path = os.path.join(output_dir, "atlgulf_fim.zarr")
            if generate_wse:
                write_ibis_to_zarr(data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                                   dem_meta=dem_meta,
                                   target_column='wse_cell_weighted_average',
                                   output_zarr_path=zarr_store_path,
                                   array_name=f"{base_name}_wse")
            if generate_depth:
                write_ibis_to_zarr(data_table=data_table_from_parquet.select("cell", "depth"),
                                   dem_meta=dem_meta,
                                   target_column='depth',
                                   output_zarr_path=zarr_store_path,
                                   array_name=f"{base_name}_depth")
            return None

        elif output_format == 'IN_MEMORY':
            if generate_wse:
                gdal_datasets['wse'] = create_in_memory_gdal_array(
                    data_table=data_table_from_parquet.select("cell", "wse_cell_weighted_average"),
                    dem_meta=dem_meta, target_column='wse_cell_weighted_average'
                )
            if generate_depth:
                gdal_datasets['depth'] = create_in_memory_gdal_array(
                    data_table=data_table_from_parquet.select("cell", "depth"),
                    dem_meta=dem_meta, target_column='depth'
                )
            return gdal_datasets

    finally:
        try:
            data_conn.drop_table("temp_cell_wse_avg", force=True)
            data_conn.drop_table("temp_cell_elevation", force=True)
        except Exception as e:
            print(f"Could not drop temporary tables. Manual cleanup may be needed. Error: {e}")
        
        if os.path.exists(temp_parquet_path):
            try:
                os.remove(temp_parquet_path)
            except OSError as e:
                print(f"Error removing temporary parquet file '{temp_parquet_path}': {e}")
        
        if data_conn and data_conn.con:
             data_conn.con.close()
        
        print("Interpolation process complete.")

def write_ibis_to_zarr(data_table: ibis.expr.types.Table, dem_meta: pd.DataFrame, target_column: str,
                    output_zarr_path: str, array_name: str, chunk_size: int = 512, db_chunk_size: int = 5000000) -> None:
    """
    The `write_ibis_to_zarr` function processes tabular data to generate a raster array representing 
    a specified variable (e.g., depth or WSE). It produces a Zarr store as output,
    adding the new raster array to it.

    Input:
        - data_table (ibis.expr.types.Table): 
            An Ibis table expression containing the source data for the raster. It must include:
              * `cell`: The unique identifier for each raster cell.
              * A target column (e.g., `depth`) with the values to be rasterized.

        - dem_meta (pd.DataFrame): 
            A Pandas DataFrame providing the raster's metadata, including `height`, `width`,
            `transform`, and `crs`.

        - target_column (str):
            The name of the column in `data_table` that contains the raster cell values.

        - output_zarr_path (str):
            The file path to the root Zarr store directory. If the store does not exist,
            it will be created.

        - array_name (str):
            The name of the new raster dataset to be created within the Zarr store.

        - chunk_size (int, optional):
            The tile size (height and width) for the Zarr array's chunks. Default is 512.
            
        - db_chunk_size (int, optional):
            The number of rows to process from the `data_table` in each batch. Default is 5000000.

    Output:
        - None: 
            Creates or updates a Zarr store at the specified path, adding a new array containing the rasterized values.

    """
    print(f"Preparing to write array '{array_name}' to Zarr store at '{output_zarr_path}'.")
    height, width = int(dem_meta['height'][0]), int(dem_meta['width'][0])
    nodata_value = -9999.0

    root_group = zarr.open_group(output_zarr_path, mode='a')
    
    z_array = root_group.create_dataset(
        name=array_name,
        shape=(height, width),
        chunks=(chunk_size, chunk_size),
        dtype='float32',
        compressor=None,  # Use Zarr's optimized default (Blosc)
        fill_value=nodata_value,
        overwrite=True
    )
    z_array.attrs['geotransform'] = json.loads(dem_meta['transform'][0])
    z_array.attrs['crs'] = dem_meta['crs'][0]
    z_array.attrs['nodata'] = nodata_value
    z_array.attrs['_ARRAY_DIMENSIONS'] = ['lat', 'lon']

    total_rows = data_table.count().execute()
    if total_rows == 0:
        print(f"No data rows found for '{array_name}'. An empty array will be created.")
        zarr.consolidate_metadata(output_zarr_path)
        return

    for offset in range(0, total_rows, db_chunk_size):
        print(f"    Processing database rows {offset:,} to {min(offset + db_chunk_size, total_rows):,}.")
        chunk_df = data_table.select('cell', target_column).limit(db_chunk_size, offset=offset).execute()
        
        if chunk_df.empty:
            continue

        row_indices = ((chunk_df['cell'] - 1) // width).to_numpy(dtype=np.int64)
        col_indices = ((chunk_df['cell'] - 1) % width).to_numpy(dtype=np.int64)
        values = chunk_df[target_column].to_numpy(dtype=np.float32)
        
        chunk_row_indices = (row_indices // chunk_size) * chunk_size
        chunk_col_indices = (col_indices // chunk_size) * chunk_size
        
        df_for_grouping = pd.DataFrame({
            'row': row_indices, 'col': col_indices, 'value': values,
            'chunk_row': chunk_row_indices, 'chunk_col': chunk_col_indices
        })
        
        grouped = df_for_grouping.groupby(['chunk_row', 'chunk_col'])

        for (y_off, x_off), block_df in grouped:
            tile = np.full((chunk_size, chunk_size), nodata_value, dtype=np.float32)
            relative_rows = block_df['row'].values - y_off
            relative_cols = block_df['col'].values - x_off
            tile[relative_rows, relative_cols] = block_df['value'].values
            
            y_slice = slice(y_off, y_off + chunk_size)
            x_slice = slice(x_off, x_off + chunk_size)
            z_array[y_slice, x_slice] = tile
        
        del chunk_df, df_for_grouping, grouped
        gc.collect()

    zarr.consolidate_metadata(output_zarr_path)
    return

def create_in_memory_gdal_array(data_table: ibis.expr.types.Table, dem_meta: pd.DataFrame, target_column: str,
                                db_chunk_size: int = 5000000, block_size: int = 512) -> gdal.Dataset:                       
    """
    The `create_in_memory_gdal_array` function processes tabular data to generate a GDAL Dataset 
    object representing a raster. The entire process occurs in memory using GDAL's 
    virtual file system (`/vsimem/`).

    Input:
        - data_table (ibis.expr.types.Table): 
            An Ibis table expression containing the source data for the raster. It must include:
              * `cell`: The unique identifier for each raster cell.
              * A target column (e.g., `depth` or `wse`) with the values to be rasterized.

        - dem_meta (pd.DataFrame): 
            A Pandas DataFrame providing the raster's metadata, including `height`, `width`, 
            `transform`, and `crs`.

        - target_column (str):
            The name of the column in `data_table` that contains the raster cell values.
            
        - db_chunk_size (int, optional):
            The number of rows to process from the `data_table` in each batch. Default is 5000000.
            
        - block_size (int, optional):
            The internal tile size (height and width) for processing the raster in chunks. Default is 512.

    Output:
        - gdal.Dataset: 
            Returns a GDAL Dataset object that exists entirely in memory. This object can be
            passed directly to other GDAL functions for further processing or file creation.

    """
    print(f"Creating GDAL array for '{target_column}' using VFS.")
    height, width = int(dem_meta['height'][0]), int(dem_meta['width'][0])
    transform_list = json.loads(dem_meta['transform'][0])
    gdal_transform = (transform_list[2], transform_list[0], transform_list[1],
                      transform_list[5], transform_list[3], transform_list[4])
    crs_wkt = dem_meta['crs'][0]
    nodata_value = -9999.0
    temp_vsimem_path = f'/vsimem/{uuid.uuid4().hex}.tif'

    # Using the GTiff driver, but point it to the virtual path (TILED=YES and SPARSE_OK=TRUE are critical for memory efficiency).
    driver = gdal.GetDriverByName('GTiff')
    creation_options = [
        "TILED=YES",
        f"BLOCKXSIZE={block_size}",
        f"BLOCKYSIZE={block_size}",
        "SPARSE_OK=TRUE",
        "COMPRESS=NONE"  
    ]
    gdal_ds = driver.Create(temp_vsimem_path, width, height, 1, gdal.GDT_Float32, creation_options)
    gdal_ds.SetGeoTransform(gdal_transform)
    gdal_ds.SetProjection(crs_wkt)
    band = gdal_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata_value)

    total_rows = data_table.count().execute()
    if total_rows == 0:
        print("No data to write for this array.")
        band.FlushCache()
        return gdal_ds

    for offset in range(0, total_rows, db_chunk_size):
        print(f"    Processing database rows {offset:,} to {min(offset + db_chunk_size, total_rows):,}...")
        chunk_df = data_table.select('cell', target_column).limit(db_chunk_size, offset=offset).execute()
        if chunk_df.empty:
            continue

        row_indices = ((chunk_df['cell'] - 1) // width).to_numpy(dtype=np.int64)
        col_indices = ((chunk_df['cell'] - 1) % width).to_numpy(dtype=np.int64)
        values = chunk_df[target_column].to_numpy(dtype=np.float32)

        df_for_grouping = pd.DataFrame({
            'row': row_indices, 'col': col_indices, 'value': values,
            'block_row': (row_indices // block_size) * block_size,
            'block_col': (col_indices // block_size) * block_size
        })

        grouped = df_for_grouping.groupby(['block_row', 'block_col'])

        # Write each tile to the virtual raster
        for (y_off, x_off), block_df in grouped:
            tile = np.full((block_size, block_size), nodata_value, dtype=np.float32)
            relative_rows = block_df['row'].values - y_off
            relative_cols = block_df['col'].values - x_off
            tile[relative_rows, relative_cols] = block_df['value'].values
            band.WriteArray(tile, xoff=x_off, yoff=y_off)

        del chunk_df, df_for_grouping, grouped
        # gc.collect()

    band.FlushCache()
    print("  In-memory VFS array created successfully.")
    return gdal_ds

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
    return