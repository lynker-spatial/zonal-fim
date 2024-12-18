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
    # data_conn.raw_sql(
    #     """
    #     -- Add idx mapping to nodes and save it
    #     CREATE OR REPLACE TABLE nodes_element_crosswalk AS
    #     SELECT 
    #         node_id,
    #         ROW_NUMBER() OVER () - 1 AS idx
    #     FROM tampa_nodes;

    #     -- Map elements (triangles) to node indices using idx
    #     CREATE OR REPLACE TABLE indexed_triangles AS
    #     SELECT
    #         el.pg_id,
    #         n1.idx AS node_id_1,
    #         n2.idx AS node_id_2,
    #         n3.idx AS node_id_3
    #     FROM elements el
    #     LEFT JOIN nodes_element_crosswalk n1 ON el.node_id_1 = n1.node_id
    #     LEFT JOIN nodes_element_crosswalk n2 ON el.node_id_2 = n2.node_id
    #     LEFT JOIN nodes_element_crosswalk n3 ON el.node_id_3 = n3.node_id
    #     WHERE n1.idx IS NOT NULL AND n2.idx IS NOT NULL AND n3.idx IS NOT NULL;

    #     -- Restore original node IDs for indexed triangles
    #     CREATE OR REPLACE TABLE original_triangles AS
    #     SELECT
    #         it.pg_id,
    #         n1.node_id AS node_id_1,
    #         n2.node_id AS node_id_2,
    #         n3.node_id AS node_id_3
    #     FROM indexed_triangles it
    #     LEFT JOIN nodes_element_crosswalk n1 ON it.node_id_1 = n1.idx
    #     LEFT JOIN nodes_element_crosswalk n2 ON it.node_id_2 = n2.idx
    #     LEFT JOIN nodes_element_crosswalk n3 ON it.node_id_3 = n3.idx;
    #     """
    # )

    # Compute WSE based on barycentric weights for all elements
    data_conn.raw_sql(
        """
        -- Map node indices to their coordinates and elevation (wse)
        CREATE OR REPLACE TEMP TABLE nodes_data AS 
        SELECT 
            node_id,
            long,
            lat,
            wse
        FROM nodes;

        -- Join triangle table with node data to get vertex coordinates
        CREATE OR REPLACE TEMP TABLE triangle_vertices AS 
        SELECT 
            t.pg_id,
            t.node_id_1, t.node_id_2, t.node_id_3,
            t.node1_weight, t.node2_weight, t.node3_weight,
            n1.long AS node1_long, n1.lat AS node1_lat, n1.wse AS node1_wse,
            n2.long AS node2_long, n2.lat AS node2_lat, n2.wse AS node2_wse,
            n3.long AS node3_long, n3.lat AS node3_lat, n3.wse AS node3_wse
        FROM original_triangles t
        JOIN nodes_data n1 ON t.node_id_1 = n1.node_id 
        JOIN nodes_data n2 ON t.node_id_2 = n2.node_id 
        JOIN nodes_data n3 ON t.node_id_3 = n3.node_id 
        ),

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
        FROM original_triangles t
        JOIN weighted_centroids wc ON t.pg_id = wc.pg_id;
        """
    )

    data_conn.con.close()
    return 



