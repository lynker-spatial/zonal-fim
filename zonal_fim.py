# test_mask.py

import argparse
import time
import os
import cfimvis.tools.mask_bounds as mb
import cfimvis.tools.geometry_manipulation as gm
import cfimvis.tools.read_schisim as rs
import cfimvis.tools.barrycentric as bc
import cfimvis.tools.format_data as fd
import cfimvis.tools.zonal_operations as zo
import cfimvis.code.bary_estimation as be
import cfimvis.code.bary_interpolation as bi

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', 'yes', 'y', '1'}:
        return True
    elif value.lower() in {'false', 'f', 'no', 'n', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")
   
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build general mask')
    parser.add_argument('-z','--generate_mask',help='activates the generate_mask as a part of preprocess in case there is a change in masks', required=False, type=str_to_bool, default=False)
    parser.add_argument('-d','--preprocess',help='activates the preprocess if it is the first time this is being executed in the domain', required=False, type=str_to_bool, default=False)
    parser.add_argument('-e','--generate_wse',help='activates the generate_wse to produce a map of WSE', required=False, type=str_to_bool, default=False)
    parser.add_argument('-g','--generate_depth',help='activates the generate_depth to produce a map of Depth', required=False, type=str_to_bool, default=True)
    parser.add_argument('-b','--zarr_format',help='activates the zarr_format that replaces geotif writing with zarr', required=False, type=str_to_bool, default=True)
    parser.add_argument('-r','--execute',help='activates the execute that runs the pipeline for generating interpolated depth', required=False, type=str_to_bool, default=True)
    parser.add_argument('-o','--dem_path',help='dem_path', required=False, type=str, default='')
    parser.add_argument('-m','--depth_path',help='depth raster path to save the file to if zarr format is chosen it automatically converts tif extension to zarr  e.g., /data/raster_v1.tif', required=False, type=str, default='')
    parser.add_argument('-q','--wse_path',help='wse raster path to save the file to if zarr format is chosen it automatically converts tif extension to zarr  e.g., /data/raster_v1.tif', required=False, type=str, default='')
    parser.add_argument('-i','--file_path',help='gr3_file_path',required=False,type=str, default='')
    parser.add_argument('-k','--shape_file_folder_path',help='shape_file_folder_path for schisim elements',required=False,type=str, default='')
    parser.add_argument('-l','--output_folder_path',help='output_folder_path for schisim elements',required=False,type=str, default='')
    parser.add_argument('-c','--database_path',help='Path to the DuckDB database file.',required=False,type=str, default='')
    parser.add_argument('-a','--mask_database_path',help='Path to the mask DuckDB database file.',required=False,type=str, default='')
    parser.add_argument('-u','--triangles_path',help='Path to the file containing the triangles dataset.',required=False,type=str, default='')
    parser.add_argument('-w','--zonal_path',help='Path to the zonal file containing the coverage fractions.',required=False,type=str, default='')
    parser.add_argument('-p','--schisim_table_name',help='Name of the SCHISM table in the DuckDB database.',required=False,type=str,default='exterior_mask_schism_boundary_atlantic_buffer_atlgulf')
    parser.add_argument('-n','--state_table_name',help='Name of the state boundaries table in the DuckDB database.',required=False,type=str,default='exterior_mask_state_boundaries_conus_atlgulf')
    parser.add_argument('-x','--levee_table_name',help='Name of the levee data table in the DuckDB database.',required=False,type=str,default='interior_mask_levee_protected_area_conus_atlgulf')
    parser.add_argument('-t','--nwm_table_name',help='Name of the National Water Model table in the DuckDB database.',required=False,type=str, default="interior_mask_nwm_lakes_conus_atlgulf")
    parser.add_argument('-f','--water_table_name',help='Name of the water bodies table in the DuckDB database.',required=False,type=str, default="interior_mask_water_polygon_conus_atlgulf")
    parser.add_argument('-s','--dissolve',help='Whether to dissolve overlapping geometries at each stage.',required=False, type=str_to_bool, default=True)

    args = vars(parser.parse_args())
    generate_mask = args['generate_mask']
    preprocess = args['preprocess']
    generate_wse = args['generate_wse']
    generate_depth = args['generate_depth']
    zarr_format = args['zarr_format']
    execute = args['execute']
    dem_path = args['dem_path']
    depth_path = args['depth_path']
    wse_path = args['wse_path']
    file_path = args['file_path']
    shape_file_folder_path = args['shape_file_folder_path']
    output_folder_path = args['output_folder_path']
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

    if generate_mask:
        print('Creating single mask ...')
        fd.convert_elements_file(shape_file_folder_path=shape_file_folder_path,output_folder_path=output_folder_path)
        mb.create_general_mask(database_path=mask_database_path, triangles_path=triangles_path,
                                schisim_table_name=schisim_table_name, state_table_name=state_table_name,
                                levee_table_name=levee_table_name, nwm_table_name=nwm_table_name,
                                water_table_name=water_table_name, dissolve=dissolve)
        print('Masking complete. \n')
    # #_____________________________

    if preprocess:
        # Reproject raster
        # file_name, file_extension = os.path.splitext(dem_path)
        # output_dem_path = f"{file_name}_4326{file_extension}"
        # print('Reprojecting to EPSG:4326 ...')
        # fd.reproject_dem(dem_path, output_dem_path)
        # print('Reprojection complete.\n')
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
        gm.add_point_geo(database_path, 'nodes', 'lat', 'long')
        gm.extract_elevation(dem_path=dem_path, database_path=database_path) # output_dem_path
        print('Elevation extraction complete.\n')
        print('Storing DEM metadata ...')
        fd.store_metadata(dem_path=dem_path, database_path=database_path, table_name="dem_metadata") # output_dem_path
        print('DEM metadata storing complete.\n')

    if execute or preprocess:
        print('Reading gr3 file...')
        start_section_1 = time.time()
        if 'point_df' not in locals() and 'point_df' not in globals():
            point_df, _, _ = rs.read_gr3(file_path)
            gm.write_to_database(database_path, 'nodes', df=point_df)
        gm.mask_nodes(database_path, 'nodes', 'masked_nodes')
        gm.add_elevation(database_path, 'masked_nodes', 'nodes_elevation')
        end_section_1 = time.time()
        time_section_1 = end_section_1 - start_section_1
        print(f"Time taken for section 1: {time_section_1:.2f} seconds")
        print('gr3 reading process complete. \n')
    # # _____________________________

    if preprocess:
        print('Calculating barycentric ...')
        mb.filter_valid_elements(data_database_path=database_path)
        bc.compute_3d_barycentric(database_path=database_path, node_table_name='masked_nodes',
                                    element_table_name='null_filtered_masked_elements')
        # Output triangle_weights
        print('Completed barycentric. \n')
        print('Completed preprocessing. \n')

    if execute:
        print('Barycentric interpolation...')
        start_section_2 = time.time()
        be.estimate(database_path)                
        end_section_2 = time.time()
        time_section_2 = end_section_2 - start_section_2
        print(f"Time taken for section 2: {time_section_2:.2f} seconds")
        # output triangle_barycentric
        print('Completed barycentric interpolation. \n')

        print('Zonal Interpolation...')
        start_section_3 = time.time()
        bi.interpolate(database_path=database_path)
        end_section_3 = time.time()
        time_section_3 = end_section_3 - start_section_3
        print(f"Time taken for section 3: {time_section_3:.2f} seconds")
       
        print('\nWriting rasters...')
        start_section_4 = time.time()
        bi.make_wse_depth_rasters(database_path=database_path,
                                    generate_depth=generate_depth, output_depth_path=depth_path, 
                                    output_wse_path=wse_path, generate_wse=generate_wse, zarr_format=zarr_format)
        end_section_4 = time.time()
        time_section_4 = end_section_4 - start_section_4
        print(f"Time taken for section 4: {time_section_4:.2f} seconds")

        total_time = time_section_1 + time_section_2 + time_section_3 + time_section_4
        total_hours = int(total_time // 3600)
        total_minutes = int((total_time % 3600) // 60)
        total_seconds = total_time % 60

        # Print total time in hours, minutes, and seconds
        print(f"Total time taken: {total_hours} hours, {total_minutes} minutes, {total_seconds:.2f} seconds")
        print('Script end.')


   

