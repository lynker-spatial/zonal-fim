# test_mask.py

import argparse
from tools import mask_bounds as mb
from tools import geometry_manipulation as gm
from tools import read_schisim as rs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build general mask')
    parser.add_argument('-o','--s3_path',help='s3_path', required=True, type=str)
    parser.add_argument('-i','--file_path',help='gr3_file_path',required=True,type=str)
    parser.add_argument('-c','--database_path',help='Path to the DuckDB database file.',required=True,type=str)
    parser.add_argument('-a','--mask_database_path',help='Path to the mask DuckDB database file.',required=True,type=str)
    parser.add_argument('-u','--triangles_path',help='Path to the file containing the triangles dataset.',required=True,type=str)
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
    schisim_table_name = args['schisim_table_name']
    state_table_name = args['state_table_name']
    levee_table_name = args['levee_table_name']
    nwm_table_name = args['nwm_table_name']
    water_table_name = args['water_table_name']
    dissolve = args['dissolve']

    print('Creating single mask ...')
    mb.create_general_mask(database_path=mask_database_path, triangles_path=triangles_path, 
                           schisim_table_name=schisim_table_name, state_table_name=state_table_name, 
                           levee_table_name=levee_table_name, nwm_table_name=nwm_table_name, 
                           water_table_name=water_table_name, dissolve=dissolve) 
    print('Msking complete. \n')
    print('Reading gr3 file.')
    point_df, _, _ = rs.read_gr3(file_path)
    gm.write_to_database(database_path, 'nodes', point_df)
    gm.add_point_geo(database_path, 'nodes', 'lat', 'long')
    print('gr3 reading process complete. \n')
    print('Finding none overlapping nodes.')
    gm.get_none_overlapping(s3_path=s3_path, database_path=database_path, point_gdf_table='nodes')
    print('Script end.')