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

    # 

