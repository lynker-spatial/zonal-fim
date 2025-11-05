# format_data.py
import ibis
from ibis import _
import json
import geopandas as gpd
import pandas as pd 
import rasterio
import rasterio
from rasterio.crs import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
import zipfile
import os
import tempfile


def setup_databse(database_path: str) -> None:
    """
    Sets up a DuckDB database with spatial extension loaded.

    Args:
        database_path (str): Path to the DuckDB database file.

    The function establishes a connection to the DuckDB database and loads the spatial extension.
    """
    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql('INSTALL spatial')
    data_conn.raw_sql('LOAD spatial')
    data_conn.con.close()
    return

def convert_elements_file(zip_file_path: str, output_folder_path: str, save_parqeut: bool=True) -> None:
    """
    Unzips a shapefile zip file, processes its contents, and saves the data as a GeoPackage and optionally as a Parquet file.

    Args:
        zip_file_path (str): The file path to the zip file containing the shapefile.
        output_folder_path (str): The folder path where the shapefile contents will be saved.
        save_parqeut (bool): Flag to determine whether to save the data as a Parquet file (default is True).

    Functionality:
        - Extracts the contents of the zip file to the specified output folder.

    Returns:
        None
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(zip_file_path, 'r') as z:
            z.extractall(temp_dir)
            
            shapefile_path = None
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith('.shp'):
                        shapefile_path = os.path.join(root, file)
                        break
                if shapefile_path:
                    break
            
            if not shapefile_path:
                raise FileNotFoundError("No .shp file found within the extracted contents of the zip archive.")
            gdf = gpd.read_file(shapefile_path)

    gdf = gdf.rename(columns={'elementID': 'pg_id'})
    if 'pg_id' in gdf.columns:
        gdf = gdf[['pg_id', 'geometry']]
    else:
        print("'elementID' column not found. The schema may be different than expected.")
    
    os.makedirs(output_folder_path, exist_ok=True)
    gdf.to_file(os.path.join(output_folder_path, 'ElementPolygons.gpkg'), layer="ElementPolygons", driver='GPKG')
    
    if save_parqeut:
        if not gdf.empty:
            gdf.to_parquet(os.path.join(output_folder_path, 'agElementPolygons.parquet'), index=False)
        else:
            print("GeoDataFrame is empty. Skipping parquet file creation.")
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
    """
    Extracts metadata from a DEM (Digital Elevation Model) file and stores it in a specified DuckDB table.

    Args:
        dem_path (str): Path to the DEM file (usually a raster file such as GeoTIFF).
        database_path (str): Path to the DuckDB database where metadata will be stored.
        table_name (str): The name of the table where the metadata will be stored.

    Returns:
        None
    """
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
    """
    Reprojects a DEM (Digital Elevation Model) raster file to the target CRS (EPSG:4326).

    Args:
        dem_path (str): Path to the input DEM file (typically in a raster format like GeoTIFF).
        output_dem_path (str): Path where the reprojected DEM file will be saved.

    Returns:
        None
    """
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

def setup_crosswalk_table(database_path: str, node_id_path: str) -> None:
    """
    Sets up a crosswalk table in the database by processing a CSV file with node IDs and storing it as a DuckDB table.

    Args:
        database_path (str): The file path to the DuckDB database where the crosswalk table will be created.
        node_id_path (str): The file path to the CSV file containing node IDs to be used for the crosswalk.

    Returns:
        None
    """
    data_conn = ibis.duckdb.connect(database_path)
    
    cross_walk = pd.read_csv(node_id_path, header=None, names=['node_id_gr3'])
    cross_walk.reset_index(inplace=True)
    cross_walk.rename(columns={'index': 'node_id_nc'}, inplace=True)
    cross_walk['node_id_nc'] += 1
    temp_file = None

    directory = 'temp'
    if not os.path.exists(directory):
        os.makedirs(directory)
    temp_file = os.path.join(directory, 'temp.parquet')
    cross_walk.to_parquet(temp_file, index=False)
    del cross_walk

    # Create crosswalk table
    data_conn.raw_sql(
        f"""
        CREATE OR REPLACE TABLE node_cross_walk AS
        SELECT * FROM '{temp_file}'
        """
    )

    # Clean up temporary file and directory
    if temp_file:
        os.remove(temp_file)
        os.rmdir(directory)

    data_conn.con.close()
    return
