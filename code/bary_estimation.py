# bary_estimation.py

import duckdb
import ibis
from ibis import _
import rasterio
import pandas as pd
import numpy as np

def estimate(database_path: str) -> None:
    data_conn = ibis.duckdb.connect(database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
    except: 
        data_conn.raw_sql('INSTALL spatial')
        data_conn.raw_sql('LOAD spatial')    
    
    # Mapping polygons to nodes ----------- << this could be part of preprocessing and table stored beforehand
    data_conn.raw_sql(
        """
        -- Map node indices to their elements
        CREATE OR REPLACE TABLE index_table AS
        WITH node_mapping AS (
            SELECT
                node_id,
                ROW_NUMBER() OVER () - 1 AS idx
            FROM nodes
        ),
        mapped_triangles AS (
            SELECT
                el.pg_id,
                nm1.idx AS node_id_1,
                nm2.idx AS node_id_2,
                nm3.idx AS node_id_3
            FROM elements el
            LEFT JOIN node_mapping nm1 ON el.node_id_1 = nm1.node_id
            LEFT JOIN node_mapping nm2 ON el.node_id_2 = nm2.node_id
            LEFT JOIN node_mapping nm3 ON el.node_id_3 = nm3.node_id
            WHERE nm1.idx IS NOT NULL AND nm2.idx IS NOT NULL AND nm3.idx IS NOT NULL
        )
        SELECT 
            pg_id,
            node_id_1,
            node_id_2,
            node_id_3
        FROM mapped_triangles;
        """
    )

    # Compute WSE based on barycentric weights for all elements
    data_conn.raw_sql(
        """
        -- Map node indices to their coordinates and elevation (wse)
        WITH nodes_data AS (
            SELECT 
                node_id,
                long,
                lat,
                wse
            FROM nodes_elevation
        ),

        -- Join triangle table with node data to get vertex coordinates
        triangle_vertices AS (
            SELECT 
                t.pg_id,
                t.node_id_1, t.node_id_2, t.node_id_3,
                t.node1_weight, t.node2_weight, t.node3_weight,
                n1.long AS node1_long, n1.lat AS node1_lat, n1.wse AS node1_wse,
                n2.long AS node2_long, n2.lat AS node2_lat, n2.wse AS node2_wse,
                n3.long AS node3_long, n3.lat AS node3_lat, n3.wse AS node3_wse
            FROM triangle t
            JOIN nodes_data n1 ON t.node_id_1 = n1.node_id 
            JOIN nodes_data n2 ON t.node_id_2 = n2.node_id 
            JOIN nodes_data n3 ON t.node_id_3 = n3.node_id 
        ),

        -- Compute weighted centroids
        weighted_centroids AS (
            SELECT
                pg_id,
                (node1_weight * node1_long + node2_weight * node2_long + node3_weight * node3_long) AS centroid_long,
                (node1_weight * node1_lat + node2_weight * node2_lat + node3_weight * node3_lat) AS centroid_lat,
                (node1_weight * node1_wse + node2_weight * node2_wse + node3_weight * node3_wse) AS wse_weighted_average
            FROM triangle_vertices
        ),

        -- Combine original triangle data with computed centroids
        combined_table AS (
            SELECT
                t.*,
                wc.centroid_long,
                wc.centroid_lat,
                wc.wse_weighted_average
            FROM triangle t
            JOIN weighted_centroids wc ON t.pg_id = wc.pg_id
        )

        -- Replace the existing triangle table with the new data
        CREATE OR REPLACE TABLE triangle AS 
        SELECT * FROM combined_table;
        """
    )



