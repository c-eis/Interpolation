#!/usr/bin/env python
import sys
from netCDF4 import Dataset
import numpy as np
from scipy.signal import medfilt
from scipy.ndimage.filters import median_filter

def smooth(file_old, file_new):
    data = Dataset(file_old, mode='r')          # open old file
    x = data.variables['x']                # save coordinates
    x_coor = x[:]
    y = data.variables['y']
    y_coor = y[:]
    z = data.variables['z']
    vel = z[:, :]
    data.close()                         # close old file
    #vel_new = medfilt(vel, [7, 7])
    vel_new = median_filter(vel, size=7, mode='reflect')
    data = Dataset(file_new, 'w')  # open new file
    data.createDimension('x', len(x_coor))
    data.createDimension('y', len(y_coor))
    x = data.createVariable('x', 'f8', ('x',))
    x[:] = x_coor
    y = data.createVariable('y', 'f8', ('y',))
    y[:] = y_coor
    z = data.createVariable('z', 'f4', ('y', 'x'))
    z[:, :] = vel_new
    data.close()                         # close new file


def main():
    smooth(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
