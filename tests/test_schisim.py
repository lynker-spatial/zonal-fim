import xarray as xr
from datetime import datetime
import pytest
import os
os.chdir("..")
import cfimvis.tools.read_schisim as rs

def test_nc_file_rows():
    # Generate URL with current date
    today = datetime.now().strftime("%Y%m%d")
    url = f"https://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.{today}/analysis_assim_coastal_atlgulf/nwm.t00z.analysis_assim_coastal.total_water.tm00.atlgulf.nc"
    
    # Open the dataset
    try:
        point_df = rs.read_netcdf(url, is_url=True, verbose=False)
    except Exception as e:
        pytest.fail(f"Failed to open dataset: {e}")
    
    # Check the wse variable
    variable_name = 'wse'  
    if variable_name not in point_df:
        pytest.fail(f"Variable '{variable_name}' not found in the dataset")
    
    # Assert the first dimension has appropriate nodes
    assert point_df[variable_name].shape[0] == 10481055, f"Expected 10481055 nodes, got {point_df[variable_name].shape[0]}"