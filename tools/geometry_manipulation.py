# geometry_manipulation.py
import duckdb
import ibis
from ibis import _
import rasterio
from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape

def add_point_geo(database_path: str, lat_col_nam: str, long_col_name: str):
    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql('LOAD spatial')
    data_conn.list_tables()
    data_conn.raw_sql(f""" 
                  ALTER TABLE nodes ADD COLUMN geometry GEOMETRY; 
                  UPDATE nodes SET geometry = ST_Point({long_col_name}, {lat_col_nam}); 
                  """)
    data_conn.close()
    return

def get_none_overlapping(s3_path: str, point_gdf: gpd.GeoDataFrame([])) -> gpd.GeoDataFrame([]):

    # Vectorize raster
    with rasterio.open(s3_path) as src:
        raster_crs = src.crs
        mask = src.dataset_mask()
        # Extract shapes from the raster
        raster_polygons = [
            shape(geom) for geom, val in shapes(mask, transform=src.transform) if val > 0
        ]
    # Create a GeoDataFrame from raster polygons
    raster_gdf = gpd.GeoDataFrame(geometry=raster_polygons, crs=raster_crs)

    # Match crs
    if raster_gdf.crs != point_gdf.crs:
        points_gdf = points_gdf.to_crs(raster_gdf.crs)
    print(f"Computing under crs: {raster_gdf.crs}")

    points_outside_raster = points_gdf[~points_gdf.within(raster_gdf.geometry.iloc[0])]
    return points_outside_raster
