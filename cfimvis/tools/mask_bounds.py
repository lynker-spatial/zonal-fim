# tools/mask_bounds.py
import duckdb
import ibis
from ibis import _
import geopandas as gpd
import os
import rasterio
from rasterio.mask import mask

ibis.options.interactive = True

def create_general_mask(database_path: str,
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
    try:
        mask_conn.raw_sql('LOAD spatial')
    except: 
        mask_conn.raw_sql('INSTALL spatial')
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
    print('Completed step_1 ...')
    step_2 = gpd.GeoDataFrame(geometry=state.buffer(5))
    step_2 = step_2.dissolve(by=None)
    step_2 = gpd.overlay(step_2, state, how='difference')
    step_2 = gpd.overlay(step_2, sch_b, how='difference')
    step_2 = gpd.overlay(step_1, step_2, how='union').dissolve(by=None)
    step_2 = step_2[['geometry']]
    print('Completed step_2 ...')
    if dissolve:
        step_3 = gpd.overlay(step_2, levee, how='union').dissolve(by=None)
        step_3 = step_3[['geometry']]
        print('Completed step_3 ...')
        step_4 = gpd.overlay(step_3, nwm, how='union').dissolve(by=None)
        step_4 = step_4[['geometry']]
        print('Completed step_4 ...')
        step_5 = gpd.overlay(step_4, water, how='union').dissolve(by=None)
        step_5 = step_5[['geometry']]
        print('Completed step_5 ...')
    else:
        step_3 = gpd.overlay(step_2, levee, how='union')
        step_3 = step_3[['geometry']]
        print('Completed step_3 ...')
        step_4 = gpd.overlay(step_3, nwm, how='union')
        step_4 = step_4[['geometry']]
        print('Completed step_4 ...')
        step_5 = gpd.overlay(step_4, water, how='union')
        step_5 = step_5[['geometry']]
        print('Completed step_5 ...')
    # Ensure crs is set
    step_5 = step_5.set_crs("EPSG:4326") 
    
    # Write to database
    directory = 'temp'
    if not os.path.exists(directory):
        os.makedirs(directory)
    temp_file = os.path.join(directory, 'mask.parquet')
    step_5.to_parquet(temp_file, index=False)
    mask_conn.raw_sql(
        f"""
        CREATE OR REPLACE TABLE step_5 AS 
        SELECT * FROM '{temp_file}'
        """
    )
    # Cleanup
    print('Completed masking.')
    os.remove(temp_file)
    os.rmdir(directory)
    mask_conn.con.close()
    return

def mask_triangles(database_path: str, triangles_path: str) -> None:
    """
    Filters triangles based on the 'masked_coverage_fraction' table and updates the triangles table.

    Args:
        database_path (str): Path to the DuckDB database file.
        triangles_path (str): Path to the file containing triangle data to be filtered.

    Output:
        - A new or updated 'triangles' table containing only triangles whose 'pg_id' exists 
          in the 'masked_coverage_fraction' table.
    """
    data_conn = ibis.duckdb.connect(database_path)
    data_conn.raw_sql(
        f"""
        CREATE OR REPLACE TABLE triangles AS 
        SELECT t.* 
        FROM '{triangles_path}' t
        WHERE t.pg_id IN (SELECT pg_id FROM masked_coverage_fraction)
        """
    )
    data_conn.con.close()
    return

def filter_valid_elements(data_database_path: str, table_name:str) -> None:
    """
    Filters elements based on valid node IDs and creates a new table with valid elements.

    Args:
        data_database_path (str): Path to the DuckDB database file.
        table_name (str): Name of the table containing the nodes to which elevation data will be added.

    Output:
        - A new or updated table 'null_filtered_masked_elements' containing elements from 
          'masked_elements' where 'node_id_1', 'node_id_2', and 'node_id_3' are found in 
          the 'masked_nodes' table.
    """
    data_conn = ibis.duckdb.connect(data_database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
    except: 
        data_conn.raw_sql('INSTALL spatial')
        data_conn.raw_sql('LOAD spatial')
  
    # Also filter for null values in elevation (points outside domain)
    data_conn.raw_sql(
        f""" 
        CREATE OR REPLACE TABLE null_filtered_masked_elements AS
        SELECT *
        FROM masked_elements AS me
        WHERE node_id_1 IN (SELECT node_id FROM '{table_name}') 
            AND node_id_2 IN (SELECT node_id FROM '{table_name}') 
            AND node_id_3 IN (SELECT node_id FROM '{table_name}');
        """
    )

    data_conn.con.close()
    return

def mask_raster(mask_database_path: str, raster_path: str) -> None:
    """
    Applies a spatial mask to a raster (DEM) based on geometries from a database table.

    Args:
        mask_database_path (str): Path to the DuckDB database containing the mask geometries.
        raster_path (str): Path to the raster (DEM) file to be masked.

    Output:
        - A new masked raster file ('DEM_masked_4326.tif') saved to the 'data' directory. 
          The raster is masked using geometries from the 'step_5' table in the DuckDB database 
          and reprojected to the raster’s coordinate reference system (CRS) if necessary.
    """
    # Load mask
    mask_conn = ibis.duckdb.connect(mask_database_path)
    try:
        mask_conn.raw_sql('LOAD spatial')
    except: 
        mask_conn.raw_sql('INSTALL spatial')
        mask_conn.raw_sql('LOAD spatial')
    mask_gdf = mask_conn.table('step_5').execute()
    mask_gdf = mask_gdf.set_crs("EPSG:4326")

    # Load DEM
    with rasterio.open(raster_path) as src:
        raster_meta = src.meta
    if mask_gdf.crs != raster_meta['crs']:
        mask_gdf = mask_gdf.to_crs(raster_meta['crs'])

    # Mask DEM
    with rasterio.open(raster_path) as src:
        out_image, out_transform = mask(src, mask_gdf.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({"driver": "GTiff",
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "transform": out_transform,
                        "crs": "EPSG:4326"})
    # Save masked DEM
    with rasterio.open("data/DEM_masked_4326.tif", 'w', **out_meta) as dest:
        dest.write(out_image)
    mask_conn.con.close()
    return

def filter_nodes(database_path: str) -> None:
    """
    Filters nodes based on valid node IDs and updates the 'nodes' table.

    Args:
        database_path (str): Path to the DuckDB database file.

    Output:
        - A new or updated 'nodes' table containing only the nodes whose 'node_id' 
          is found in the 'valid_nodes_elevations' table.
    """
    data_conn = ibis.duckdb.connect(database_path)
    try:
        data_conn.raw_sql('LOAD spatial')
    except: 
        data_conn.raw_sql('INSTALL spatial')
        data_conn.raw_sql('LOAD spatial')

    data_conn.raw_sql(
        """
        CREATE OR REPLACE TABLE nodes AS
        SELECT *
        FROM nodes
        WHERE node_id IN (
            SELECT node_id
            FROM valid_nodes_elevations
        );
        """
    )

    data_conn.con.close()
    return