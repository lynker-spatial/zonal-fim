# tools/make_masks_duckdb.py
import os
import zipfile
import ibis
import shutil
from pathlib import Path

def make_masks_duckdb(zip_folder_path: str, duckdb_path: str) -> None:
    """
    Extracts shapefiles from zipped archives in a folder and saves them as tables in a DuckDB database.
    Utilizes the DuckDB Ibis spatial extension for handling geospatial data.

    Args:
        zip_folder_path (str): 
            The path to the folder containing `.zip` files with shapefiles (`.shp`).
        duckdb_path (str): 
            The path to the DuckDB database file where the extracted shapefiles will be saved.

    Returns:
        None: 
            The function does not return any value. It saves the extracted shapefiles as tables in the DuckDB database.

    Example:
        zip_folder = "/path/to/zipped/shapefiles"
        duckdb_file = "/path/to/database/masks.duckdb"
        make_masks_duckdb(zip_folder, duckdb_file)

    Notes:
        - The function assumes valid `.zip` archives and `.shp` files.
        - Errors during processing are logged, allowing the function to continue with other files.
    """
        
    # Spinup to the DuckDB database
    mask_conn = ibis.duckdb.connect(duckdb_path)
    mask_conn.raw_sql('INSTALL spatial')
    mask_conn.raw_sql('LOAD spatial')

    # Temporary extraction folder
    extraction_path = os.path.join(zip_folder_path, "temp_extracted")

    # Cleanup the temporary directory
    if os.path.exists(extraction_path):
        shutil.rmtree(extraction_path)
    os.makedirs(extraction_path, exist_ok=True)

    # Iterate through all files in the folder
    for file_name in os.listdir(zip_folder_path):
        if file_name.endswith('.zip'):
            zip_path = os.path.join(zip_folder_path, file_name)
            try:
                # Extract zip contents to the temporary folder
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extraction_path)
                    # Add to database
                    for extracted_file in zip_ref.namelist():
                        if extracted_file.endswith('.shp'):
                            shapefile_path = os.path.join(extraction_path, extracted_file)
                            table_name = Path(file_name).stem + '_atlgulf'
                            try:
                                mask_conn.raw_sql(f"""
                                                    CREATE OR REPLACE TABLE {table_name} AS 
                                                    SELECT * FROM ST_Read('{shapefile_path}')
                                                    """)

                                print(f"Saved {table_name} to DuckDB.")
                            except Exception as e:
                                print(f"Error saving {table_name} to DuckDB: {e}")
            except Exception as e:
                print(f"Error processing {zip_path}: {e}")
    
    # Save crs info
    mask_conn.raw_sql("""
    CREATE TABLE IF NOT EXISTS metadata (
        table_name STRING,
        crs STRING
        )
    """)
    mask_conn.raw_sql("""
        INSERT INTO metadata (table_name, crs)
        VALUES ('masks', 'EPSG:4326')
    """)

    # Cleanup the temporary folder
    shutil.rmtree(extraction_path, ignore_errors=True)
    mask_conn.close()
    print("All shapefiles have been successfully saved to the DuckDB database.")
    return
