#tool to create viewsheds from a DEM raster

# library imports
import numpy as np
import rasterio as rio
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterstats
from osgeo import gdal

# https://gis.stackexchange.com/questions/434417/an-example-for-viewshed-analysis-in-python


ox=[292294] # obs x in srs units
oy=[6571861] # obs y in srs units
oz=[2] # obs height metres above DEM
src_ds = gdal.Open('data/ni_dtm/ni_dtm.tif')
ds = gdal.ViewshedGenerate(
    src_ds.GetRasterBand(1),
    "GTiff",
    "outtest.tif",
    ["INTERLEAVE=BAND"],
    ox[0],
    oy[0],
    oz[0],
    0,  # targetHeight
    255,  # visibleVal
    0,  # invisibleVal
    0,  # outOfRangeVal
    -1.0,  # noDataVal,
    0.85714,  # dfCurvCoeff
    gdal.GVM_Edge,
    0,  # maxDistance
    heightMode=gdal.GVOT_MIN_TARGET_HEIGHT_FROM_GROUND,
    options=["UNUSED=YES"],
)