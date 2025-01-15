# test_mask.py

import argparse
import time
from cfimvis.tools import mask_bounds as mb
from cfimvis.tools import geometry_manipulation as gm
from cfimvis.tools import read_schisim as rs
from cfimvis.tools import barrycentric as bc
from cfimvis.tools import zonal_operations as zo
from cfimvis.code import bary_estimation as be
from cfimvis.code import bary_interpolation as bi 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build general mask')
    parser.add_argument('-o','--s3_path',help='s3_path', required=True, type=str)
    parser.add_argument('-i','--file_path',help='gr3_file_path',required=True,type=str)
    parser.add_argument('-c','--database_path',help='Path to the DuckDB database file.',required=True,type=str)
    parser.add_argument('-a','--mask_database_path',help='Path to the mask DuckDB database file.',required=True,type=str)
    parser.add_argument('-u','--triangles_path',help='Path to the file containing the triangles dataset.',required=True,type=str)
    parser.add_argument('-w','--zonal_path',help='Path to the zonal file containing the coverage fractions.',required=True,type=str)
    parser.add_argument('-p','--schisim_table_name',help='Name of the SCHISM table in the DuckDB database.',required=False,type=str,default='exterior_mask_schism_boundary_atlantic_buffer_atlgulf')
    parser.add_argument('-n','--state_table_name',help='Name of the state boundaries table in the DuckDB database.',required=False,type=str,default='exterior_mask_state_boundaries_conus_atlgulf')
    parser.add_argument('-x','--levee_table_name',help='Name of the levee data table in the DuckDB database.',required=False,type=str,default='interior_mask_levee_protected_area_conus_atlgulf')
    parser.add_argument('-t','--nwm_table_name',help='Name of the National Water Model table in the DuckDB database.',required=False,type=str, default="interior_mask_nwm_lakes_conus_atlgulf")
    parser.add_argument('-f','--water_table_name',help='Name of the water bodies table in the DuckDB database.',required=False,type=str, default="interior_mask_water_polygon_conus_atlgulf")
    parser.add_argument('-s','--dissolve',help='Whether to dissolve overlapping geometries at each stage.',required=False,type=bool, default=True)


    args = vars(parser.parse_args())
    s3_path = args['s3_path']
    file_path = args['file_path']
    database_path = args['database_path']
    mask_database_path = args['mask_database_path']
    triangles_path = args['triangles_path']
    zonal_path = args['zonal_path']
    schisim_table_name = args['schisim_table_name']
    state_table_name = args['state_table_name']
    levee_table_name = args['levee_table_name']
    nwm_table_name = args['nwm_table_name']
    water_table_name = args['water_table_name']
    dissolve = args['dissolve']
    print('\n')

    generate_mask = False
    preprocess = True

    if generate_mask:
        print('Creating single mask ...')
        mb.create_general_mask(database_path=mask_database_path, triangles_path=triangles_path, 
                            schisim_table_name=schisim_table_name, state_table_name=state_table_name, 
                            levee_table_name=levee_table_name, nwm_table_name=nwm_table_name, 
                            water_table_name=water_table_name, dissolve=dissolve) 
        print('Masking complete. \n')
    #_____________________________

    if preprocess:
        # Ingest coverage fraction data
        print('Ingesting zonal output file ...')
        gm.write_to_database(database_path, 'coverage_fraction', df_path=zonal_path) 
        zo.filter_masked(database_path)
        print('Added zonal output file to duckdb.\n')
        # Ingest and filter triangles data
        print("Ingesting triangles data ...")
        mb.mask_triangles(database_path, triangles_path)
        print('Added triangle elements to duckdb.\n')
        print('Ingest element to node crosswalk ...')
        point_df, elements_df, _ = rs.read_gr3(file_path) 
        gm.write_to_database(database_path, 'nodes', df=point_df) 
        gm.write_to_database(database_path, 'elements', df=elements_df)
        rs.mask_elements(database_path)
        print('Added elements crosswalk to duckdb.\n')
        print('Extracting elevation for nodes ...')
        gm.add_point_geo(database_path, 'nodes', 'lat', 'long') # -- needs to run many times --- maybe needed
        gm.extract_elevation(s3_path=s3_path, database_path=database_path)
        # gm.mask_nodes(database_path, 'nodes_elevation', 'masked_nodes_elevation') # -- needs to run many times
        print('Elevation extraction complete.\n')

    print('Reading gr3 file.')
    start_section_1 = time.time()
    point_df, _, _ = rs.read_gr3(file_path) # -- needs to run many times
    gm.write_to_database(database_path, 'nodes', df=point_df) # -- needs to run many times
    gm.mask_nodes(database_path, 'nodes', 'masked_nodes') # -- needs to run many times
    gm.add_elevation(database_path, 'masked_nodes', 'nodes_elevation')
    # gm.add_point_geo(database_path, 'masked_nodes', 'lat', 'long') # -- needs to run many times --- maybe needed
    end_section_1 = time.time()
    time_section_1 = end_section_1 - start_section_1
    print(f"Time taken for section 1: {time_section_1:.2f} seconds")
    print('gr3 reading process complete. \n')
    #_____________________________

    # print('Finding none overlapping nodes.')
    # gm.get_none_overlapping(s3_path=s3_path, database_path=database_path, point_gdf_table='nodes')
    # print('Found all none overlapping nodes. \n')
    
    # print('Masking elements, nodes, and DEM...')
    # mb.filter_valid_elements(data_database_path=database_path, mask_database_path=mask_database_path) 
    # # These will be used from now on -> nu+ll_filtered_masked_elements,  valid_nodes_elevations
    # mb.mask_raster(mask_database_path=mask_database_path, raster_path=s3_path)
    # print('Masking complete. \n')
    if preprocess:
        print('Extracting elevation for nodes and calculating barycentric ...')
        bc.compute_3d_barycentric(database_path=database_path, raster_path=s3_path, node_table_name='masked_nodes', 
                                    element_table_name='masked_elements')
        # Output triangle_weights 
        print('Completed barycentric. \n')


    # print('Barycentric interpolation...')
    # start_section_2 = time.time()
    # mb.filter_nodes(database_path=database_path) # -- needs to run many times
    # be.estimate(database_path=database_path)                  # -- needs to run many times
    # end_section_2 = time.time()
    # time_section_2 = end_section_2 - start_section_2
    # print(f"Time taken for section 2: {time_section_2:.2f} seconds")
    # # output triangle_barycentric
    # print('Completed barycentric interpolation. \n')

    # # print('Zonal operations...')
    # # zo.read_zonal_outputs(database_path=database_path, zonal_output_path=zonal_path)
    # # bi.interpolate(database_path=database_path, s3_path=s3_path)
    # # bi.make_depth_raster(dem_path="data/DEM_masked_4326.tif")
    # total_time = time_section_1 + time_section_2
    # total_hours = int(total_time // 3600)
    # total_minutes = int((total_time % 3600) // 60)
    # total_seconds = total_time % 60

    # # Print total time in hours, minutes, and seconds
    # print(f"Total time taken: {total_hours} hours, {total_minutes} minutes, {total_seconds:.2f} seconds")
    print('Script end.')

    