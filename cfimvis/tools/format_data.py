# format_data.py

import duckdb
import ibis
from ibis import _
import geopandas as gpd
import rasterio
from pathlib import Path


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
