# tools/mask_bounds.py
import duckdb
import ibis
from ibis import _
import geopandas as gpd
import os

ibis.options.interactive = True

def create_general_mask(database_path: str, triangles_path: str,
                        schisim_table_name: str, state_table_name: str,
                        levee_table_name: str, nwm_table_name: str,
                        water_table_name: str, dissolve: bool) -> None:
    
    """
    Creates a general coastal mask using spatial data from various input sources 
    and applies the mask to a triangles dataset. The result is stored in a DuckDB database.

    Parameters:
    -----------
    database_path : str
        Path to the DuckDB database file.
    triangles_path : str
        Path to the file containing the triangles dataset.
    schisim_table_name : str
        Name of the SCHISM table in the DuckDB database.
    state_table_name : str
        Name of the state boundaries table in the DuckDB database.
    levee_table_name : str
        Name of the levee data table in the DuckDB database.
    nwm_table_name : str
        Name of the National Water Model table in the DuckDB database.
    water_table_name : str
        Name of the water bodies table in the DuckDB database.
    dissolve : bool
        Whether to dissolve overlapping geometries at each stage.

    Returns:
    --------
    None

    Example:
    --------
    >>> create_general_mask(
    ...     database_path="example.duckdb",
    ...     triangles_path="triangles.parquet",
    ...     schisim_table_name="schisim",
    ...     state_table_name="state",
    ...     levee_table_name="levee",
    ...     nwm_table_name="nwm",
    ...     water_table_name="water",
    ...     dissolve=True
    ... )
    """
     
    # Create a single mask
    mask_conn = ibis.duckdb.connect(database_path)
    mask_conn.raw_sql('LOAD spatial')
    sch_b = mask_conn.table(schisim_table_name).execute()
    sch_b = sch_b.set_crs("EPSG:4326")
    state = mask_conn.table(state_table_name).execute()
    state = state.set_crs("EPSG:4326")
    levee =  mask_conn.table(levee_table_name).execute()
    levee = levee.set_crs("EPSG:4326")
    nwm =  mask_conn.table(nwm_table_name).execute()
    nwm = nwm.set_crs("EPSG:4326")
    water =  mask_conn.table(water_table_name).execute()
    water = water.set_crs("EPSG:4326")

    # Setup bounds
    sch_bb = gpd.GeoDataFrame(geometry=sch_b.buffer(3))
    sch_bb = sch_bb.dissolve(by=None)
    state = gpd.overlay(state, sch_bb, how="intersection")
    levee = gpd.overlay(levee, sch_bb, how="intersection")
    nwm = gpd.overlay(nwm, sch_bb, how="intersection")
    water = gpd.overlay(water, sch_bb, how="intersection")

    # Implement 5 step masking logic
    # 1. "mask_schism_boundary_atlantic.shp","exterior"
    # 2. "mask_state_boundaries_conus.shp","exterior"
    # 3. "mask_levee_protected_area_conus.shp","interior"
    # 4. "mask_nwm_lakes_conus.shp","interior"
    # 5. "mask_water_polygon_conus.shp","interior" 
    step_1 = gpd.overlay(sch_bb, sch_b, how='difference')
    step_1 = step_1[['geometry']]
    step_2 = gpd.GeoDataFrame(geometry=state.buffer(5))
    step_2 = step_2.dissolve(by=None)
    step_2 = gpd.overlay(step_2, state, how='difference')
    step_2 = gpd.overlay(step_2, sch_b, how='difference')
    step_2 = gpd.overlay(step_1, step_2, how='union').dissolve(by=None)
    step_2 = step_2[['geometry']]
    if dissolve:
        step_3 = gpd.overlay(step_2, levee, how='union').dissolve(by=None)
        step_3 = step_3[['geometry']]
        step_4 = gpd.overlay(step_3, nwm, how='union').dissolve(by=None)
        step_4 = step_4[['geometry']]
        step_5 = gpd.overlay(step_4, water, how='union').dissolve(by=None)
        step_5 = step_5[['geometry']]
    else:
        step_3 = gpd.overlay(step_2, levee, how='union')
        step_3 = step_3[['geometry']]
        step_4 = gpd.overlay(step_3, nwm, how='union')
        step_4 = step_4[['geometry']]
        step_5 = gpd.overlay(step_4, water, how='union')
        step_5 = step_5[['geometry']]
    # Ensure crs is set
    step_5 = step_5.set_crs("EPSG:4326") 
    
    # Write to database
    directory = 'temp'
    if not os.path.exists(directory):
        os.makedirs(directory)
    temp_file = os.path.join(directory, 'mask.parquet')
    step_5.to_parquet(temp_file, index=False)
    mask_conn.raw_sql(f"""
                        CREATE OR REPLACE TABLE step_5 AS 
                        SELECT * FROM '{temp_file}'
                        """)
    mask_conn.raw_sql(f"""
                        CREATE OR REPLACE TABLE triangles AS 
                        SELECT * FROM '{triangles_path}'
                        """)
    
    # Mask triangles that are fully overlap with costal mask
    mask_conn.raw_sql("""
                        CREATE OR REPLACE TABLE triangles_masked AS
                        SELECT 
                            t.* -- Keep pg_id from the table
                        FROM 
                            triangles AS t, step_5 AS s
                        WHERE 
                            ST_Intersects(t.geometry, s.geometry) -- Keep triangles that intersect
                            AND NOT ST_Within(t.geometry, s.geometry) -- Exclude those completely within step_5
                            OR NOT ST_Intersects(t.geometry, s.geometry); -- Keep triangles that do not intersect at all
                        """)
    # Cleanup
    os.remove(temp_file)
    os.rmdir(directory)
    mask_conn.close()
    return