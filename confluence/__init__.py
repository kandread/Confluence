"""

Implementation of Confluence framework for integrating and diagnosing
SWOT-estimated discharge.

"""

import netCDF4 as netcdf


def read_file(ncfile):
    """Read reach time series data from NetCDF file."""
    with netcdf.Dataset(ncfile) as f:
        rt = f.groups['Reach_Timeseries']
        Q = rt.variables['Q'][:]
        W= rt.variables['W'][:]
        H = rt.variables['H'][:]
    return Q, H, W


def integrate(Q):
    """Apply integrator algorithms to estimated discharge."""
