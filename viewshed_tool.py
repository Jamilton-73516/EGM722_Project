#tool to create viewsheds from a DEM raster and perform analysis

# library imports
import numpy as np
import rasterio as rio
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterstats
import cartopy.crs as ccrs
from osgeo import gdal
from osgeo import osr
from shapely.geometry import Point, LineString, Polygon

#ox=[292890] # obs x in srs units
#oy=[428388] # obs y in srs units
#oz=[2] # obs height metres above DEM
#src_ds = gdal.Open('data/ni_dtm/ni_dtm.tif')

#ds = gdal.ViewshedGenerate(
    #src_ds.GetRasterBand(1),
    #"GTiff",
    #"test1.tif",
    #["INTERLEAVE=BAND"],
    #ox[0],
    #oy[0],
    #oz[0],
    #0,  # targetHeight
    #255,  # visibleVal
    #0,  # invisibleVal
    #0,  # outOfRangeVal
    #-1.0,  # noDataVal,
    #0.85714,  # dfCurvCoeff
    #gdal.GVM_Edge,
    #0,  # maxDistance
    #heightMode=gdal.GVOT_MIN_TARGET_HEIGHT_FROM_GROUND,
    #options=["UNUSED=YES"],
#)


# import shapefiles, transform to WGS84
landcover = gpd.read_file('data/landcover/landcover.shp')
landcover = landcover.to_crs(epsg=4326) #wgs84
buildings = gpd.read_file('data/buildings/buildings_binevenagh.shp')
buildings = buildings.to_crs(epsg=4346) #wgs84

# import temporary viewshed while viewshed tool is broken
viewshed = gpd.read_file('data/viewshed/trackviewshed_polygon_val1.shp')
viewshed = viewshed.to_crs(epsg=4326) #wgs84

#print(buildings.head())

print(landcover.loc[0])
print(landcover.loc[landcover['Shape_Area'] > 5000])



#test = "data/ni_dtm/ni_dtm.tif"
#ds = gdal.Open(test)
#print("Projection", ds.GetProjection())

#dataset = rio.open("data/ni_dtm/ni_dtm.tif")
#print(dataset.crs)


#plt.imshow(array)
#plt.colorbar()
#print(dsorig.GetGeoTransform())
#print(dsorig.GetProjection())



#dsRepro = gdal.Warp("NI_repro.tif", dsorig, dstSRS = "EPSG:4326")


