#!/usr/bin/env python
from __future__ import print_function
from time import *
from netCDF4 import Dataset
import numpy as np
import sys

# sys.argv[1]: input , sys.argv[2]:output


def add_edge(file_old, file_new, file_edge):
    '''Add edges from a prior field to data'''
    edge_data = Dataset(file_edge, mode='r')
    e = edge_data.variables['z']
    edge = e[:, :]
    edge_data.close()
    data = Dataset(file_old, mode='r')          # open old file
    z = data.variables['z']
    vel = z[:, :]
    data.close()                         # close old file
    n = vel.shape[0]-1
    m = vel.shape[1]-1
    vel[0, :] = edge[0, :]
    vel[:, 0] = edge[:, 0]
    vel[n, :] = edge[n, :]
    vel[:, m] = edge[:, m]
    data_new = Dataset(file_new, 'a')  # open new file
    z = data.variables['z']
    z[:, :] = vel
    data_new.close()                         # close new file


def main():
    add_edge(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
