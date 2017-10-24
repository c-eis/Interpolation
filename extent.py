#!/usr/bin/env python
import gdal
from gdalconst import GA_ReadOnly
import sys


def extent(filename):
    data = gdal.Open(filename, GA_ReadOnly)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    data = None
    return [minx, miny, maxx, maxy]


def main():
    ext = extent(sys.argv[1])
    print ext[0], ext[1], ext[2], ext[3]


if __name__ == "__main__":
    main()
