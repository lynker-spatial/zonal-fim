import os
import logging
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("schism_converter")

def read_hgrid_gr3(filename):
    with open(filename, 'r') as f:
        title = f.readline().strip()
        num_elements, num_nodes = map(int, f.readline().strip().split())
        
        # Read nodes data
        nodes = []
        for _ in range(num_nodes):
            parts = f.readline().strip().split()
            node_id = int(parts[0])  
            x = float(parts[1])      
            y = float(parts[2])      
            depth = float(parts[3])  
            nodes.append((node_id, x, y, depth))
        
        # Read elements data
        elements = []
        for _ in range(num_elements):
            parts = f.readline().strip().split()
            elem_id = int(parts[0])  # Element ID
            nodes_in_elem = list(map(int, parts[1:]))
            elements.append((elem_id, nodes_in_elem))
    
    return {
        "title": title,
        "nodes": nodes,
        "elements": elements
    }

def hgrid_elements_to_gdf(grid_data, output_path=None, crs="EPSG:4326", driver="GPKG", include_mean_wse=True):
    """
    Converts elements from a parsed SCHISM hgrid dataset to a GeoDataFrame 
    containing pg_id, num_nodes, node_id_1, node_id_2, node_id_3, and optionally 
    mean_wse, along with the geometry.
    
    Parameters:
    -----------
    grid_data : dict
        The dictionary returned by read_hgrid_gr3 containing 'nodes' and 'elements'.
    output_path : str, optional
        File path to save the generated spatial file.
    crs : str, optional
        Coordinate reference system for the GeoDataFrame (default WGS 84 / "EPSG:4326").
    driver : str, optional
        Fiona/GDAL driver to use when exporting the file (default "GPKG").
    include_mean_wse : bool, optional
        If True, calculates and appends the mean water surface elevation (depth/z-value)
        of the component nodes as a column (default True).
        
    Returns:
    --------
    gdf : geopandas.GeoDataFrame
        The generated spatial GeoDataFrame of element geometries.
    """
    nodes = grid_data.get("nodes", [])
    elements = grid_data.get("elements", [])
    
    if not nodes:
        logger.error("No node data found in grid_data.")
        return None
        
    if not elements:
        logger.error("No element data found in grid_data.")
        return None

    max_node_id = max(node[0] for node in nodes)
    
    coords = np.zeros((max_node_id + 1, 2), dtype=np.float64)
    wse_lookup = np.zeros((max_node_id + 1,), dtype=np.float64)
    
    logger.info("Initializing node lookup index mapping for element conversion...")
    for node_id, x, y, depth in nodes:
        coords[node_id] = (x, y)
        wse_lookup[node_id] = depth
    
    successful_elements = []
    failed_elements = []
    
    logger.info(f"Processing {len(elements)} elements...")
    
    for idx, (elem_id, nodes_in_elem) in enumerate(tqdm(elements, desc="Generating Element Polygons")):
        try:
            if not nodes_in_elem:
                raise ValueError("Element contains empty node list.")
            
            # Element type is specified by the first integer
            num_nodes_in_elem = nodes_in_elem[0]
            
            # Extract actual node IDs (excluding the element type indicator)
            node_ids = nodes_in_elem[1 : 1 + num_nodes_in_elem]
            
            if len(node_ids) < 3:
                raise ValueError(
                    f"Incomplete element node sequence: expected at least 3 nodes, found {len(node_ids)}."
                )
                
            # Perform boundary check on node IDs against lookup arrays
            for nid in node_ids:
                if nid <= 0 or nid > max_node_id:
                    raise ValueError(f"Node ID {nid} is out of valid range (1 - {max_node_id}).")
            
            elem_coords = coords[node_ids]
            poly = Polygon(elem_coords)
            
            elem_dict = {
                "pg_id": elem_id,
                "num_nodes": num_nodes_in_elem,
                "node_id_1": node_ids[0],
                "node_id_2": node_ids[1],
                "node_id_3": node_ids[2],
                "geometry": poly
            }
            
            if include_mean_wse:
                elem_dict["mean_wse"] = float(np.mean(wse_lookup[node_ids]))
                
            if len(node_ids) > 3:
                for i, nid in enumerate(node_ids[3:], start=4):
                    elem_dict[f"node_id_{i}"] = nid
            
            successful_elements.append(elem_dict)
            
        except Exception as e:
            failed_elements.append({
                "index": idx,
                "element_id": elem_id,
                "error": str(e),
                "nodes_in_elem": nodes_in_elem
            })
            
    # Log processing statistics
    logger.info(f"Geometry generation complete: {len(successful_elements)} elements successfully created.")
    
    if failed_elements:
        logger.warning(f"{len(failed_elements)} elements failed to process.")
        for i, fail in enumerate(failed_elements[:15]):
            logger.error(
                f"Failure #{i+1} at index {fail['index']} (ID: {fail['element_id']}): "
                f"{fail['error']} | Raw nodes: {fail['nodes_in_elem']}"
            )
        if len(failed_elements) > 15:
            logger.error(f"... and {len(failed_elements) - 15} additional element construction failures.")
            
    if successful_elements:
        df = pd.DataFrame(successful_elements)
        gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
    else:
        logger.warning("No elements were successfully processed. Returning empty GeoDataFrame.")
        gdf = gpd.GeoDataFrame(columns=["pg_id", "num_nodes", "node_id_1", "node_id_2", "node_id_3", "geometry"], crs=crs)
        
    if output_path and successful_elements:
        try:
            output_dir = os.path.dirname(output_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
                
            logger.info(f"Writing elements to: {output_path} with driver '{driver}'...")
            gdf.to_file(output_path, driver=driver)
            logger.info("File export completed successfully.")
        except Exception as e:
            logger.error(f"Failed to write file to path '{output_path}': {e}")
            
    return gdf

def hgrid_nodes_to_gdf(grid_data, output_path=None, crs="EPSG:4326", driver="GPKG"):
    """
    Converts node data from a parsed SCHISM hgrid dataset to a Point GeoDataFrame
    containing node_id, long, lat, and wse (water surface elevation / depth).
    
    Parameters:
    -----------
    grid_data : dict
        The dictionary returned by read_hgrid_gr3 containing 'nodes'.
    output_path : str, optional
        File path to save the generated spatial file.
    crs : str, optional
        Coordinate reference system for the GeoDataFrame (default WGS 84 / "EPSG:4326").
    driver : str, optional
        Fiona/GDAL driver to use when exporting the file (default "GPKG").
        
    Returns:
    --------
    gdf : geopandas.GeoDataFrame
        The generated spatial GeoDataFrame of Point geometries.
    """
    nodes = grid_data.get("nodes", [])
    
    if not nodes:
        logger.error("No node data found in grid_data.")
        return None
        
    logger.info(f"Converting {len(nodes)} nodes to Point geometries using vectorized mappings...")
    
    try:
        # Cast lists to NumPy structured array 
        node_array = np.array(nodes, dtype=np.float64)
        
        node_ids = node_array[:, 0].astype(np.int64)
        long_coords = node_array[:, 1]
        lat_coords = node_array[:, 2]
        wse_values = node_array[:, 3]
        
        if np.isnan(long_coords).any() or np.isnan(lat_coords).any():
            logger.warning("Warning: NaN values detected in node coordinates.")
        
        geometry = gpd.points_from_xy(long_coords, lat_coords)
        
        data_properties = {
            "node_id": node_ids,
            "long": long_coords,
            "lat": lat_coords,
            "wse": wse_values
        }
            
        gdf = gpd.GeoDataFrame(data_properties, geometry=geometry, crs=crs)
        
        logger.info(f"Successfully generated Point GeoDataFrame.")
        
    except Exception as e:
        logger.error(f"Failed to generate node Point geometries: {e}")
        return None
        
    if output_path and (gdf is not None) and (not gdf.empty):
        try:
            output_dir = os.path.dirname(output_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
                
            logger.info(f"Writing nodes to: {output_path} with driver '{driver}'...")
            gdf.to_file(output_path, driver=driver)
            logger.info("File export completed successfully.")
        except Exception as e:
            logger.error(f"Failed to write file to path '{output_path}': {e}")
            
    return gdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SCHISM hgrid.gr3 files to elements and nodes spatial formats.")
    parser.add_argument("--input", "-i", type=str, required=True, help="Path to raw hgrid.gr3 file.")
    parser.add_argument("--elements-output", "-eo", type=str, default=None, help="Path to save the elements spatial file (e.g. .gpkg).")
    parser.add_argument("--elements-output-pq", "-eopq", type=str, default=None, help="Path to save the elements Parquet file.")
    parser.add_argument("--nodes-output", "-no", type=str, default=None, help="Path to save the nodes spatial file (e.g. .gpkg).")
    parser.add_argument("--nodes-output-pq", "-nopq", type=str, default=None, help="Path to save the nodes Parquet file.")
    parser.add_argument("--crs", type=str, default="EPSG:4326", help="Coordinate reference system (default: EPSG:4326).")
    parser.add_argument("--driver", type=str, default="GPKG", help="Fiona/GDAL driver for vector export (default: GPKG).")
    parser.add_argument("--exclude-mean-wse", action="store_true", help="Do not calculate or include mean water surface elevation (WSE) for elements.")
    
    args = parser.parse_args()

    filename = args.input
    elements_output = args.elements_output
    elements_output_pq = args.elements_output_pq
    nodes_output = args.nodes_output
    nodes_output_pq = args.nodes_output_pq

    # Read raw grid file
    if os.path.exists(filename):
        grid_data = read_hgrid_gr3(filename)
        
        if elements_output or elements_output_pq:
            elements_gdf = hgrid_elements_to_gdf(
                grid_data=grid_data, 
                output_path=elements_output, 
                crs=args.crs, 
                driver=args.driver,
                include_mean_wse=not args.exclude_mean_wse
            )
            if elements_output_pq and (elements_gdf is not None) and (not elements_gdf.empty):
                try:
                    output_dir = os.path.dirname(elements_output_pq)
                    if output_dir and not os.path.exists(output_dir):
                        os.makedirs(output_dir, exist_ok=True)
                    logger.info(f"Writing elements to parquet: {elements_output_pq}...")
                    elements_gdf.to_parquet(elements_output_pq, index=False)
                    logger.info("Elements parquet export completed successfully.")
                except Exception as e:
                    logger.error(f"Failed to write parquet file to path '{elements_output_pq}': {e}")
        
        if nodes_output or nodes_output_pq:
            nodes_gdf = hgrid_nodes_to_gdf(
                grid_data=grid_data,
                output_path=nodes_output,
                crs=args.crs,
                driver=args.driver
            )
            if nodes_output_pq and (nodes_gdf is not None) and (not nodes_gdf.empty):
                try:
                    output_dir = os.path.dirname(nodes_output_pq)
                    if output_dir and not os.path.exists(output_dir):
                        os.makedirs(output_dir, exist_ok=True)
                    logger.info(f"Writing nodes to parquet: {nodes_output_pq}...")
                    nodes_gdf.to_parquet(nodes_output_pq, index=False)
                    logger.info("Nodes parquet export completed successfully.")
                except Exception as e:
                    logger.error(f"Failed to write parquet file to path '{nodes_output_pq}': {e}")
    else:
        logger.warning(f"File '{filename}' was not found. Please verify the path to run.")