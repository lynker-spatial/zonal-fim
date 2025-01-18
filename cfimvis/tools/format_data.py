# format_data.py

import duckdb
import ibis
from ibis import _
import json
import geopandas as gpd
import rasterio
from pathlib import Path
import rasterio
from rasterio.crs import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling



def convert_elements_file(shape_file_folder_path: str, output_folder_path: str) -> None:
    """
    Converts a shapefile containing elements to  parquet, and gpkg files.
    
    Args:
        shape_file_folder_path (str): Path to the shapefile folder.
        output_folder_path (str): Path to save the output folder.
    """
    gdf_polys = gpd.read_file(shape_file_folder_path)

    # Write the GeoDataFrame as a GeoParquet file
    gdf_polys.to_parquet(output_folder_path+"/ElementPolygons.parquet", engine="pyarrow", index=False)
    # Write the GeoDataFrame as a Geopackage file
    gdf_polys.to_file(output_folder_path+"/ElementPolygons.gpkg", driver="GPKG", layer='ElementPolygons')
    return

def exctract_mask(mask_database_path: str, output_folder_path: str) -> None:
    """
    Extracts spatial data from a DuckDB mask database and saves it as a GeoPackage.

    Args:
        mask_database_path (str): The file path to the DuckDB database containing spatial data.
        output_folder_path (str): The folder path where the GeoPackage will be saved.

    Functionality:
        - Saves the GeoDataFrame to a GeoPackage file named "ElementPolygons.gpkg" with layer name mask in the specified output folder.

    Returns:
        None
    """
    mask_conn = ibis.duckdb.connect(mask_database_path)
    try:
        mask_conn.raw_sql('LOAD spatial')
    except:
        mask_conn.raw_sql('INSTALL spatial')
        mask_conn.raw_sql('LOAD spatial')

    geopackage_path = output_folder_path+'/ElementPolygons.gpkg'  
    mask_table = mask_conn.table("step_5")
    gdf = mask_table.execute()
    gdf=gdf[['geometry']]
    gdf = gdf.set_crs("EPSG:4326")
    # Write the GeoDataFrame to GeoPackage
    gdf.to_file(geopackage_path, layer='mask', driver="GPKG")
    mask_conn.con.close()
    return

def store_metadata(dem_path: str, database_path: str, table_name: str) -> None:
    data_conn = ibis.duckdb.connect(database_path)
    with rasterio.open(dem_path) as src:
        raster_meta = src.meta
        width = raster_meta['width']
        height = raster_meta['height']
        nodata = raster_meta['nodata']
        crs = src.crs.to_string() if isinstance(src.crs, CRS) else str(src.crs)
        transform = raster_meta['transform']  
        transform_values = list(transform)
        transform_string = json.dumps(transform_values)
        meta = 'DEM'
        table_name = "dem_metadata"
        # Empty table creation 
        if table_name not in data_conn.list_tables():
                data_conn.raw_sql(
                    f"""
                    CREATE TABLE IF NOT EXISTS {table_name} (
                        dem_metadata STRING,
                        crs STRING,
                        height INTEGER,
                        width INTEGER,
                        nodata FLOAT,
                        transform STRING
                    )
                    """
                )

        # Data insertion
        data_conn.raw_sql(
            f"""
            INSERT INTO {table_name} (dem_metadata, crs, height, width, nodata, transform)
            VALUES ('{meta}', '{crs}', {height}, {width}, {nodata}, '{transform_string}');
            """
        )
    data_conn.con.close()
    return

def reproject_dem(dem_path:str, output_dem_path:str) -> None:
    # Define the target CRS
    target_crs = "EPSG:4326"
    with rasterio.open(dem_path) as src:
        raster_crs = src.crs
        # Check if the CRS matches the target CRS
        if raster_crs == target_crs:
            print("CRS matches the target CRS.")
            return
        else:
            # Calculate the transform and dimensions for the target CRS
            transform, width, height = calculate_default_transform(
                src.crs, target_crs, src.width, src.height, *src.bounds
            )
            
            # Update the metadata to reflect the new CRS
            profile = src.profile.copy()
            profile.update({
                'crs': target_crs,
                'transform': transform,
                'width': width,
                'height': height
            })
            
            # Reproject and save the output raster
            with rasterio.open(output_dem_path, 'w', **profile) as dst:
                reproject(
                    source=rasterio.band(src, 1),  # Input raster's first band
                    destination=rasterio.band(dst, 1),  # Output raster's first band
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest
                )

        print(f"Reprojected raster saved at: {output_dem_path}")
    return
