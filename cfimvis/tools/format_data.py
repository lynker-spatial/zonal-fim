# format_data.py

import duckdb
import ibis
from ibis import _
import rasterio
from rasterio.features import shapes
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape
from rasterio.session import AWSSession
from rio_tiler.errors import PointOutsideBounds
from rio_tiler.io import COGReader
import rasterio
import boto3
import numpy as np
from tqdm import tqdm
import os

def convert_elements_file(shape_file_path: str, output_folder_path: str) -> None:
    """
    Converts a shapefile containing elements to  parquet, and gpkg files.
    
    Args:
        shape_file_path (str): Path to the shapefile.
        output_folder_path (str): Path to save the output folder.
    """
    gdf_polys = gpd.read_file(shape_file_path)

    # Write the GeoDataFrame as a GeoParquet file
    gdf_polys.to_parquet(output_folder_path+"/agElementPolygons.parquet", engine="pyarrow", index=False)
    # Write the GeoDataFrame as a Geopackage file
    gdf_polys.to_file(output_folder_path+"/agElementPolygons.gpkg", driver="GPKG")
    return

