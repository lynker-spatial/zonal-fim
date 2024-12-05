# barrycentric.py

import duckdb
import ibis
from ibis import _
import geopandas as gpd
from rio_tiler.io import COGReader
from rio_tiler.errors import PointOutsideBounds
from tqdm import tqdm
import numpy as np
import os


def extract_elevation(raster_path: str, points_df: gpd.GeoDataFrame()) -> gpd.GeoDataFrame():
    with COGReader(raster_path) as cog:
        # Fix crs mismatches
        raster_crs = cog.dataset.crs
        points_crs = points_df.crs
        if points_crs != raster_crs:
            print("CRS mismatch detected. Transforming points to match raster CRS.")
            points_df = points_df.to_crs(raster_crs)
            print(f"new CRS: {points_df.crs}")
        else:
            print("CRS match. No transformation required.")

        # Extract raster values for the transformed points
        coords = [(geom.x, geom.y) for geom in points_df["geometry"]]
        values = []
        for x, y in tqdm(coords):
            try:
                value = cog.point(x, y, coord_crs=raster_crs).array[0]
                values.append(float(value))
            except PointOutsideBounds:
                # Assign NaN if the point is outside the raster bounds
                values.append(np.nan)


    # Add the extracted raster values to the Ibis DuckDB table
    points_df["elevation"] = values
    # transform back to the 4326
    points_df.to_crs('EPSG:4326')
    return points_df

def calculate_slope(vertex, other_vertex1, other_vertex2):
    # Vectors from the vertex to the other two vertices
    vector1 = other_vertex1 - vertex
    vector2 = other_vertex2 - vertex
    
    # Calculate the change in elevation over the distance
    dz1 = vector1[2]  
    dz2 = vector2[2]  
    
    # Calculate the distances between the vertices in the horizontal plane 
    distance1 = np.linalg.norm(vector1[:2])  
    distance2 = np.linalg.norm(vector2[:2])  
    
    # Calculate the slope (rate of change in elevation) for each vector
    slope1 = abs(dz1 / distance1) if distance1 != 0 else 0
    slope2 = abs(dz2 / distance2) if distance2 != 0 else 0
    
    # The final slope is the average of the two slopes
    return (slope1 + slope2) / 2

def calculate_barycentric_weights(triangle_points):
    # Calculate the slopes for each vertex
    A_point = triangle_points[0]
    B_point = triangle_points[1]
    C_point = triangle_points[2]
    slope_A = calculate_slope(A_point, B_point, C_point)
    slope_B = calculate_slope(B_point, A_point, C_point)
    slope_C = calculate_slope(C_point, A_point, B_point)

    # Normalize the slopes so that they sum to 1
    total_slope = slope_A + slope_B + slope_C
    # Set equal weights if it is a flat region
    if total_slope == 0: 
        w_A = 0.33333333333
        w_B = 0.33333333333
        w_C = 0.33333333333
    else:
        w_A = slope_A / total_slope
        w_B = slope_B / total_slope
        w_C = slope_C / total_slope

    # Compute the weighted centroid based on the slope weights
    centroid = (w_A * A_point + w_B * B_point + w_C * C_point)
    
    # Calculate vectors from each vertex to the centroid
    vec_A = centroid - A_point
    vec_B = centroid - B_point
    vec_C = centroid - C_point

    # Calculate the normal vector to the triangle 
    normal = np.cross(B_point - A_point, C_point - A_point)
    normal_length = np.linalg.norm(normal)

    # Calculate areas for weights as the magnitude of the cross products
    area_A = np.linalg.norm(np.cross(vec_B, vec_C)) / normal_length
    area_B = np.linalg.norm(np.cross(vec_C, vec_A)) / normal_length
    area_C = np.linalg.norm(np.cross(vec_A, vec_B)) / normal_length

    # Total area for normalization
    total_area = area_A + area_B + area_C

    # Calculate normalized barycentric weights
    weight_A = area_A / total_area
    weight_B = area_B / total_area
    weight_C = area_C / total_area

    weights_3d = np.array([weight_A, weight_B, weight_C])
    return weights_3d

def compute_3d_barycentric(database_path: str, raster_path: str, table_name: str):
    data_conn = duckdb.connect(database_path)
    data_conn.raw_sql('LOAD spatial')
    point_gdf = data_conn.table(table_name).execute()
    point_gdf = point_gdf.set_crs(epsg='4326')

    # Extract elevation
    point_gdf = extract_elevation(raster_path, point_gdf)
    
