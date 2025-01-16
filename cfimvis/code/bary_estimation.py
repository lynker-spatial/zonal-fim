# bary_estimation.py

import duckdb
import ibis
from ibis import _
import rasterio
import pandas as pd
import numpy as np

def estimate(database_path: str, output_database_path: str) -> None:
    """
    The `estimate` function processes spatial data stored in a DuckDB database to compute 
    weighted centroids and water surface elevation (WSE) for triangular elements. It leverages 
    the spatial extension of DuckDB and performs SQL-based operations for data preprocessing 
    and computation. The results are stored in a table named `triangle_barycentric`.

    Input:
        - database_path (str): The file path to the DuckDB database containing the spatial 
          data tables. These tables are expected to include:
            * nodes: Contains `node_id`, `long`, `lat`, and `wse` (water surface elevation).
            * elements: Contains `pg_id` (triangle ID) and the three nodes defining each triangle 
              with their barycentric weights (`node1_weight`, `node2_weight`, `node3_weight`).

    Output:
        - None: The function modifies the database by creating or updating the table 
          `triangle_barycentric` with the following computed columns:
            * `centroid_long`: Longitude of the weighted centroid.
            * `centroid_lat`: Latitude of the weighted centroid.
            * `wse_weighted_average`: Weighted average of the water surface elevation for the triangle.

    Example:
        1. Ensure `masked_nodes` and `elements` tables exist in the DuckDB database:
            * `nodes` table includes `node_id`, `long`, `lat`, `wse`.
            * `elements` table includes `pg_id`, `node_id_1`, `node_id_2`, `node_id_3`,
              and barycentric weights for the three nodes.

        2. Call the function:
            estimate("path/to/database.duckdb")

        3. Result:
            A new table `triangle_barycentric` will be created with the original triangle data, 
            computed centroids, and WSE averages.

    """

    data_conn = ibis.duckdb.connect(database_path)
    out_data_conn = ibis.duckdb.connect(output_database_path)
    try:
        out_data_conn.raw_sql('LOAD spatial')
    except: 
        out_data_conn.raw_sql('INSTALL spatial')
        out_data_conn.raw_sql('LOAD spatial')    
    out_data_conn.raw_sql(f"ATTACH '{database_path}' AS compute_db;")
    # Compute WSE based on barycentric weights for all elements
    out_data_conn.raw_sql(
        """
        -- Map node indices to their coordinates and elevation (wse)
        CREATE OR REPLACE TEMP TABLE nodes_data AS 
        SELECT 
            node_id,
            long,
            lat,
            wse
        FROM compute_db.masked_nodes;

        -- Join triangle table with node data to get vertex coordinates
        CREATE OR REPLACE TEMP TABLE triangle_vertices AS 
        SELECT 
            t.pg_id,
            t.node_id_1, t.node_id_2, t.node_id_3,
            t.node1_weight, t.node2_weight, t.node3_weight,
            n1.long AS node1_long, n1.lat AS node1_lat, n1.wse AS node1_wse,
            n2.long AS node2_long, n2.lat AS node2_lat, n2.wse AS node2_wse,
            n3.long AS node3_long, n3.lat AS node3_lat, n3.wse AS node3_wse
        FROM compute_db.triangle_weights t
        JOIN nodes_data n1 ON t.node_id_1 = n1.node_id 
        JOIN nodes_data n2 ON t.node_id_2 = n2.node_id 
        JOIN nodes_data n3 ON t.node_id_3 = n3.node_id; 

        -- Compute weighted centroids
        CREATE OR REPLACE TEMP TABLE weighted_centroids AS 
        SELECT
            pg_id,
            (node1_weight * node1_long + node2_weight * node2_long + node3_weight * node3_long) AS centroid_long,
            (node1_weight * node1_lat + node2_weight * node2_lat + node3_weight * node3_lat) AS centroid_lat,
            (node1_weight * node1_wse + node2_weight * node2_wse + node3_weight * node3_wse) AS wse_weighted_average
        FROM triangle_vertices;

        -- Combine original triangle data with computed centroids
        CREATE OR REPLACE TABLE triangle_barycentric AS 
        SELECT
            t.*,
            wc.centroid_long,
            wc.centroid_lat,
            wc.wse_weighted_average
        FROM compute_db.triangle_weights t
        LEFT JOIN weighted_centroids wc ON t.pg_id = wc.pg_id;
        """
    )

    data_conn.con.close()
    out_data_conn.con.close()
    return 



