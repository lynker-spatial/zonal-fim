# barrycentric.py

import warnings
import duckdb
import ibis
from ibis import _
import geopandas as gpd
import pandas as pd
from rio_tiler.io import COGReader
from rio_tiler.errors import PointOutsideBounds
from tqdm import tqdm
import numpy as np
import os
import shutil



# def extract_elevation(raster_path: str, points_df: gpd.GeoDataFrame()) -> gpd.GeoDataFrame():
#     with COGReader(raster_path) as cog:
#         # Fix crs mismatches
#         raster_crs = cog.dataset.crs
#         points_crs = points_df.crs
#         if points_crs != raster_crs:
#             print("CRS mismatch detected. Transforming points to match raster CRS.")
#             points_df = points_df.to_crs(raster_crs)
#             print(f"new CRS: {points_df.crs}")
#         else:
#             print("CRS match. No transformation required.")

#         # Extract raster values for the transformed points
#         coords = [(geom.x, geom.y) for geom in points_df["geometry"]]
#         values = []
#         for x, y in tqdm(coords):
#             try:
#                 value = cog.point(x, y, coord_crs=raster_crs).array[0]
#                 values.append(float(value))
#             except PointOutsideBounds:
#                 # Assign NaN if the point is outside the raster bounds
#                 values.append(np.nan)


#     # Add the extracted raster values to the Ibis DuckDB table
#     points_df["elevation"] = values
#     # Transform back to the 4326
#     points_df.to_crs('EPSG:4326')
#     return points_df

def calculate_slope(vertex: np.ndarray, other_vertex1: np.ndarray, other_vertex2: np.ndarray) -> float:
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

def calculate_barycentric_weights(triangle_points: np.ndarray) -> np.ndarray:
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

def compute_3d_barycentric(database_path: str, mask_database_path: str, node_table_name: str, 
                           element_table_name: str) -> None:
    
    data_conn = ibis.duckdb.connect(database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
    except: 
        data_conn.raw_sql('INSTALL spatial')
        data_conn.raw_sql('LOAD spatial')
    nodes_df = data_conn.table(node_table_name).execute()
    nodes_df = nodes_df.set_crs(epsg='4326')
    points_dem = nodes_df[['long', 'lat', 'elevation']].to_numpy()  

    output_folder = 'temp'  
    os.makedirs(output_folder, exist_ok=True)
    
    # Mapping polygons to nodes 
    data_conn.raw_sql(
        f"""
        -- Add idx mapping to nodes and save it
        CREATE OR REPLACE TABLE node_element_crosswalk AS
        SELECT 
            node_id,
            ROW_NUMBER() OVER () - 1 AS idx
        FROM {node_table_name};

        -- Map elements (triangles) to node indices using idx
        CREATE OR REPLACE TABLE indexed_triangles AS
        SELECT
            el.pg_id,
            n1.idx AS node_id_1,
            n2.idx AS node_id_2,
            n3.idx AS node_id_3
        FROM {element_table_name} el
        LEFT JOIN node_element_crosswalk n1 ON el.node_id_1 = n1.node_id
        LEFT JOIN node_element_crosswalk n2 ON el.node_id_2 = n2.node_id
        LEFT JOIN node_element_crosswalk n3 ON el.node_id_3 = n3.node_id
        WHERE n1.idx IS NOT NULL AND n2.idx IS NOT NULL AND n3.idx IS NOT NULL;

        -- Restore original node IDs for indexed triangles
        CREATE OR REPLACE TABLE original_triangles AS
        SELECT
            it.pg_id,
            n1.node_id AS node_id_1,
            n2.node_id AS node_id_2,
            n3.node_id AS node_id_3
        FROM indexed_triangles it
        LEFT JOIN node_element_crosswalk n1 ON it.node_id_1 = n1.idx
        LEFT JOIN node_element_crosswalk n2 ON it.node_id_2 = n2.idx
        LEFT JOIN node_element_crosswalk n3 ON it.node_id_3 = n3.idx;
        """
    )
    triangles_df = data_conn.table('indexed_triangles').execute()
    triangles = triangles_df[['node_id_1', 'node_id_2', 'node_id_3']].to_numpy()
    # Map node IDs to indices for the triangulation
    # id_columns = ['pg_id', 'node_id_1', 'node_id_2', 'node_id_3']
    # triangles_df = data_conn.table(element_table_name).select(id_columns).execute()
    # id_columns.remove('pg_id')
    # node_id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_df['node_id'])}
    # triangles_df[id_columns] = triangles_df[id_columns].applymap(node_id_to_index.get)
    # triangles_df.dropna(inplace=True)
    # triangles_df.reset_index(inplace=True, drop=True)
    # triangles = triangles_df.drop(columns=['pg_id']).to_numpy()
    
    
    # Process the triangle data in batches
    batch_size = 100000
    num_batches = len(triangles) // batch_size + (1 if len(triangles) % batch_size != 0 else 0)
    

    for batch_num in tqdm(range(num_batches), desc="Processing batches"):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(triangles))
        
        # Get the current batch of triangles
        batch_triangles = triangles[start_idx:end_idx]
        
        # Lists to hold results for the current batch
        weights_raw = []
        
        for tri_indices in batch_triangles:
            # Get triangle vertices using the point indices
            triangle_points_dem = points_dem[tri_indices.astype(int)]

            # Calculate Barycentric weights
            weights = calculate_barycentric_weights(triangle_points_dem)
            weights_raw.append(weights)
            
        
        # Save raw weights as well
        weights_raw_path = os.path.join(output_folder, f'batch_{batch_num}_weights.csv')
        pd.DataFrame(weights_raw).to_csv(weights_raw_path, index=False, header=False)

    print("Processing completed. Outputs saved to:", output_folder)
    del weights_raw, batch_triangles
    
    # Stich patches to a single parquet file
    weights_raw_list = []

    # Read and concatenate weighted centers and raw weights from all batches
    num_batches = len(os.listdir(output_folder)) // 2 

    for batch_num in tqdm(range(num_batches), desc="Processing batches"):
        # Read raw weights
        weights_raw_path = os.path.join(output_folder, f'batch_{batch_num}_weights.csv')
        weights_raw_df = pd.read_csv(weights_raw_path, header=None) 
        weights_raw_list.append(weights_raw_df)

    # Concatenate all batches into single DataFrames
    final_weights_raw = pd.concat(weights_raw_list, ignore_index=True)


    new_column_names = ['node1_weight', 'node2_weight', 'node3_weight']

    # Assign the new column names to the DataFrame
    final_weights_raw.columns = new_column_names

    elements_df_combined = pd.concat([triangles_df, final_weights_raw], axis=1)
    
    # Merge with triangles
    output_path = os.path.join(output_folder, 'bary_weights.parquet')
    elements_df_combined.to_parquet(output_path, index=False)
    data_conn.raw_sql(
        f"""
        CREATE OR REPLACE TABLE bary_weights AS
        SELECT * FROM '{output_path}';
        """
    )
    
    data_conn.raw_sql(
        """
        CREATE OR REPLACE TABLE triangle_weights AS
        SELECT
                el.*,
                bw.* EXCLUDE (pg_id)
        FROM 
            elements AS el
        LEFT JOIN
            bary_weights AS bw
        ON
            el.pg_id = bw.pg_id;
        """
    )

    # Save crs info
    data_conn.raw_sql(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            table_name STRING,
            crs STRING
        );
        """
    )
    data_conn.raw_sql(
        """
        INSERT INTO metadata (table_name, crs)
        VALUES 
            ('bary_weights', 'EPSG:4326'),
            ('triangle_weights', 'EPSG:4326');             
        """
    )

    # Add geometry back
    data_conn.raw_sql(f"ATTACH '{mask_database_path}' AS mask_db;")

    # Create A new table for masked elements <---- we can keep geometry form the start to avoid this something to fix later on
    data_conn.raw_sql(
        """ 
        ALTER TABLE triangle_weights ADD COLUMN geometry GEOMETRY;
        UPDATE triangle_weights
        SET geometry = mask_db.triangles.geometry
        FROM mask_db.triangles
        WHERE triangle_weights.pg_id = mask_db.triangles.pg_id;
        """
    )

    # Look for any problems
    triangle_elements = data_conn.table('triangle_weights')
    nan_counts = triangle_elements.aggregate(
        **{
            column: triangle_elements[column].isnull().sum().name(f"{column}_nan_count")
            for column in triangle_elements.columns
        }
    )
    nan_counts_result = nan_counts.execute()
    nan_flag = nan_counts_result.sum().sum()
    if nan_flag == 0:
        print("Found nan in triangles check previous steps")
        print(nan_counts_result)
    
    # Save for R
    triangle_elements = data_conn.table('triangle_weights').execute()
    triangle_elements = triangle_elements.set_crs(epsg=4326)
    print("Saving to gpkg for R ...")
    triangle_elements.to_file("data/bary_triangles.gpkg", layer="triangles", driver="GPKG")

    # Clean up
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    else:
        print(f"The folder '{output_folder}' does not exist.")
    data_conn.con.close()
    return
