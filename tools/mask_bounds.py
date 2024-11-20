# tools/mask_bounds.py
import duckdb
import ibis
from ibis import _
ibis.options.interactive = True

def create_general_mask(database_path, traingles_path, dissolve):
    # Create a single mask
    mask_conn = ibis.duckdb.connect(database_path)
    mask_conn.raw_sql('LOAD spatial')
    sch_b = mask_conn.table('exterior_mask_schism_boundary_atlantic_atlgulf').execute()
    sch_b = sch_b.set_crs("EPSG:4326")
    state = mask_conn.table('exterior_mask_state_boundaries_conus_atlgulf').execute()
    state = state.set_crs("EPSG:4326")
    levee =  mask_conn.table('interior_mask_levee_protected_area_conus_atlgulf').execute()
    levee = levee.set_crs("EPSG:4326")
    nwm =  mask_conn.table('interior_mask_nwm_lakes_conus_atlgulf').execute()
    nwm = nwm.set_crs("EPSG:4326")
    water =  mask_conn.table('interior_mask_water_polygon_conus_atlgulf').execute()
    water = water.set_crs("EPSG:4326")

    step_1 = gpd.GeoDataFrame(geometry=sch_b.buffer(3))
    step_1 = step_1.dissolve(by=None)
    step_1 = gpd.overlay(step_1, sch_b, how='difference')
    step_2 = gpd.GeoDataFrame(geometry=state.buffer(5))
    step_2 = step_2.dissolve(by=None)
    step_2 = gpd.overlay(step_2, state, how='difference')
    step_2 = gpd.overlay(step_2, sch_b, how='difference')
    step_2 = gpd.overlay(step_1, step_2, how='union').dissolve(by=None)
    if dissolve:
        step_3 = gpd.overlay(step_2, levee, how='union').dissolve(by=None)
        step_4 = gpd.overlay(step_3, nwm, how='union').dissolve(by=None)
        step_5 = gpd.overlay(step_4, water, how='union').dissolve(by=None)
    else:
        step_3 = gpd.overlay(step_2, levee, how='union')
        step_4 = gpd.overlay(step_3, nwm, how='union')
        step_5 = gpd.overlay(step_4, water, how='union')
    # Ensure crs is set
    step_5 = step_5.set_crs("EPSG:4326") 

    # write to database
    step_5.to_parquet('atlgulf/mask.parquet', index=False)
    mask_conn.raw_sql(f"""
                        CREATE OR REPLACE TABLE step_5 AS 
                        SELECT * FROM 'atlgulf/mask.parquet'
                        """)
    mask_conn.raw_sql(f"""
                        CREATE OR REPLACE TABLE triangles AS 
                        SELECT * FROM '{traingles_path}'
                        """)
    
    # Mask triangles that are fully overlap with costal mask
    mask_conn.raw_sql("""
                        CREATE OR REPLACE TABLE triangles_masked AS
                        SELECT 
                            t.* -- Keep all columns from the triangles table
                        FROM 
                            triangles AS t, step_5 AS s
                        WHERE 
                            ST_Intersects(t.geometry, s.geometry) -- Keep triangles that intersect
                            AND NOT ST_Within(t.geometry, s.geometry) -- Exclude those completely within step_5
                            OR NOT ST_Intersects(t.geometry, s.geometry); -- Keep triangles that do not intersect at all
                        """)
    return