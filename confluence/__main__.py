import sys
import netCDF4 as netcdf
import integrator
from confluence import removeFlagged


def main(ncfile):
    """Driver function for the Confluence framework."""
    f = netcdf.Dataset(ncfile)
    H = f.variables['H'][:].data
    W = f.variables['W'][:].data
    S = f.variables['S'][:].data
    H, W, S = removeFlagged(H, W, S)
    n = f.variables['n'][:].data
    Ab = f.variables['Abase'][:].data
    if 'dA' in f.variables:
        dA = f.variables['dA'][:].data
    else:
        dA = None
    mf = integrator.MeanFlow(n, Ab, H, W, S, dA, None)
    n_est, Ab_est = mf.integrate()
    integrator.write(ncfile, Ab_est, n_est)


if __name__ == '__main__':
    ncfile = sys.argv[1]
    main(ncfile)
