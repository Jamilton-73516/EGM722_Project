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
import fiona #to support different drivers for geodataframe to file

def viewshedcreate():
    ox=[274368] # obs x in srs units
    oy=[433564] # obs y in srs units
    oz=[2] # obs height metres above DEM
    src_ds_dem = gdal.Open('data/ni_dtm/ni_dtm.tif')

    ds = gdal.ViewshedGenerate(
        src_ds_dem.GetRasterBand(1),
        "GTiff", #geotiff
        "data/outputs/viewshed_northcoast.tif", #make sure this is the same as is used for the gpd dataframe
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
    vs_pix = gdal.Open('data/outputs/viewshed_northcoast.tif')
    band = vs_pix.GetRasterBand(1)
    Cols = vs_pix.RasterXSize
    Rows = vs_pix.RasterYSize
    data = band.ReadAsArray(0, 0, Cols, Rows).astype(float)

    class1 = np.where(data==255)
    print(np.sum(class1))
pixelcount()

def rast2poly():
    src_ds_rast = gdal.Open('data/outputs/viewshed_northcoast.tif')
    srcband = src_ds_rast.GetRasterBand(1)
    dst_layername = 'layernametest'
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource('data/outputs/viewshed_poly.shp')

    sp_ref = osr.SpatialReference()
    sp_ref.SetFromUserInput('EPSG:29902')

    dst_layer = dst_ds.CreateLayer(dst_layername, srs = sp_ref)

    fld = ogr.FieldDefn("visible", ogr.OFTInteger)
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("visible")

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)

    del src_ds_rast
    del dst_ds
rast2poly()

def vs_visible_select():
    vs_visible = gpd.read_file('data/outputs/viewshed_poly.shp')
    vs_visible[vs_visible['visible']==255].to_file('data/outputs/viewshed_visible.shp')
vs_visible_select()

def lcm_clip():
    vs_mask = gpd.read_file('data/outputs/viewshed_visible.shp')
    vs_mask_tm65 = vs_mask.to_crs(epsg=29902)
    landcover = gpd.read_file('data/landcover/landcover.shp')
    landcover_tm65 = landcover.to_crs(epsg=29902)
    landcover_clipped = gpd.clip(landcover_tm65, vs_mask_tm65)
    landcover_clipped.to_file('data/outputs/landcover_clipped.shp')
lcm_clip()



