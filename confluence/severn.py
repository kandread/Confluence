#!/usr/bin/env python

import glob
import netCDF4 as netcdf
import numpy as np


def read_data(filenames):
    """Read data from NetCDF files."""
    S = []
    W = []
    H = []
    Q = []
    for ncfile in filenames:
        with netcdf.Dataset(ncfile) as f:
            S.append(f['Reach_Timeseries']['S_90m'][:, 0].data)
            W.append(np.mean(f['XS_Timseries']['W'][:].data, axis=1))
            H.append(np.mean(f['XS_Timseries']['H_90m'][:].data, axis=1))
            Q.append(f['Reach_Timseries']['Q'][:, 0].data)
    Q, H, W, S = np.array(Q).T, np.array(H).T, np.array(W).T, np.array(S).T
    return Q, H, W, S


def main():
    """Main driver routine."""
    ncfiles = glob.glob("../Severn/*.nc")
    Q, H, W, S = read_data(ncfiles)
