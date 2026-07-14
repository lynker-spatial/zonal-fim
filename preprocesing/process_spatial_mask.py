import os
import argparse
import netCDF4 as nc
import numpy as np
import xarray as xr
from tqdm import tqdm

def process_spatial_join_nc(
    input_nc_path, 
    output_mask_nc_path, 
    output_csv_path,
    filter_mask_val=0,
    use_vectorized_write=True
):
    """
    Processes spatial join NetCDF data to generate a land mask and filtered node CSV.
    Optimized for execution speed and low memory footprint.
    
    Parameters:
    -----------
    input_nc_path : str
        Path to the source spatial join NetCDF file.
    output_mask_nc_path : str
        Path to save the generated land_mask.nc file.
    output_csv_path : str
        Path to save the filtered atl_orig_id.csv node list.
    filter_mask_val : int or float, optional
        The mask value to filter out (default is 0).
    use_vectorized_write : bool, optional
        If True, utilizes numpy.savetxt (highly optimized C backend). 
        If False, uses a python loop wrapped with a tqdm progress bar.
    """
    if not os.path.exists(input_nc_path):
        raise FileNotFoundError(f"Input NetCDF file not found: {input_nc_path}")

    print(f"Reading input dataset from {input_nc_path}...")
    with xr.open_dataset(input_nc_path) as ds:
        n = ds['Node'].values.astype(np.int64)
        x = ds['Longitude'].values
        y = ds['Latitude'].values
        m = ds['Mask'].values.astype(np.int32)

    if os.path.exists(output_mask_nc_path):
        os.remove(output_mask_nc_path)

    print(f"Writing land mask array to {output_mask_nc_path}....")
    with nc.Dataset(output_mask_nc_path, 'w', format='NETCDF4') as ncout:
        ncout.createDimension('nodes', len(n))
        ncm = ncout.createVariable('Mask', 'i4', ('nodes',))
        ncm[:] = m

    print("Filtering node IDs....")
    I = np.where(m != filter_mask_val)
    filtered_nodes = n[I]

    print(f"Saving {len(filtered_nodes)} node IDs to {output_csv_path}...")
    if use_vectorized_write:
        np.savetxt(output_csv_path, filtered_nodes, fmt='%d', newline='\n')
    else:
        with open(output_csv_path, 'w') as o:
            for i in tqdm(range(len(filtered_nodes)), desc="Writing node IDs"):
                o.write(str(filtered_nodes[i]) + '\n')

    print("Processing completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial join NetCDF data to generate a land mask and filtered node CSV.")
    parser.add_argument("--input", "-i", type=str, required=True, help="Path to the source spatial join NetCDF file.")
    parser.add_argument("--output-mask", "-om", type=str, required=True, help="Path to save the generated land mask NetCDF file.")
    parser.add_argument("--output-csv", "-oc", type=str, required=True, help="Path to save the filtered node ID CSV file.")
    parser.add_argument("--filter-mask-val", "-fv", type=int, default=0, help="The mask value to filter out (default is 0).")
    parser.add_argument("--use-loop-write", action="store_true", help="Disable vectorized C-level write and utilize a python loop with a progress bar.")

    args = parser.parse_args()

    for path in [args.output_mask, args.output_csv]:
        dir_name = os.path.dirname(path)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=True)

    process_spatial_join_nc(
        input_nc_path=args.input,
        output_mask_nc_path=args.output_mask,
        output_csv_path=args.output_csv,
        filter_mask_val=args.filter_mask_val,
        use_vectorized_write=not args.use_loop_write
    )