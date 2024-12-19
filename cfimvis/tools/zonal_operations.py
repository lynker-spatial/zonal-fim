# zonal_operations.py

import duckdb
import ibis
from ibis import _
import rasterio

def read_zonal_outputs(database_path: str, zonal_output_path: str) -> None:
    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql('LOAD spatial')
    data_conn.raw_sql(
        f"""
        CREATE OR REPLACE TABLE z_w AS 
        SELECT * FROM {zonal_output_path}
        """
    )
    data_conn.con.close()
    return

