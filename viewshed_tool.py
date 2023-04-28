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
from osgeo import ogr
from shapely.geometry import Point, LineString, Polygon
import sys

def viewshedcreate():
    ox=[274368] # obs x in srs units
    oy=[433564] # obs y in srs units
    oz=[2] # obs height metres above DEM
    src_ds_dem = gdal.Open('data/ni_dtm/ni_dtm.tif')

    ds = gdal.ViewshedGenerate(
        src_ds_dem.GetRasterBand(1),
        "GTiff", #geotiff
        "test1.tif", #make sure this is the same as is used for the gpd dataframe
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
        heightMode=gdal.GVOT_NORMAL,
        options=["UNUSED=YES"],
    )
viewshedcreate()

def pixelcount():
    vs_pix = gdal.Open('test1.tif')
    band = vs_pix.GetRasterBand(1)
    Cols = vs_pix.RasterXSize
    Rows = vs_pix.RasterYSize
    data = band.ReadAsArray(0, 0, Cols, Rows).astype(float)

    class1 = np.where(data==255)
    print(np.sum(class1))
pixelcount()

def rast2poly():
    src_ds_rast = gdal.Open('test1.tif')
    srcband = src_ds_rast.GetRasterBand(1)
    dst_layername = 'layernametest'
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource('data/outputs/conversiontest.shp')

    sp_ref = osr.SpatialReference()
    sp_ref.SetFromUserInput('EPSG:29902')

    dst_layer = dst_ds.CreateLayer(dst_layername, srs = sp_ref)

    fld = ogr.FieldDefn("HA", ogr.OFTInteger)
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("HA")

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)

    del src_ds_rast
    del dst_ds
rast2poly()






# import shapefiles, transform to WGS84
#landcover = gpd.read_file('data/landcover/landcover.shp')
#landcover = landcover.to_crs(epsg=4326) #wgs84
#buildings = gpd.read_file('data/buildings/buildings_binevenagh.shp')
#buildings = buildings.to_crs(epsg=4346) #wgs84

# import temporary viewshed while viewshed tool is broken
#viewshed = gpd.read_file('data/viewshed/trackviewshed_polygon_val1.shp')
#viewshed = viewshed.to_crs(epsg=4326) #wgs84

#print(buildings.head())
#print(landcover.loc[0])
#print(landcover.loc[landcover['Shape_Area'] > 5000])
#print(landcover['Shape_Area'].sum() / 1000000)



#clipping landcover by viewshed
#landcover_clip = gpd.clip(landcover, viewshed, keep_geom_type=False)

#export geodataframe to shapefile
#landcover_clip.to_file('lcm_clip_by_viewshed.shp', crs="epsg:4326")










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


