#tool to create viewsheds from a DEM raster and perform analysis

# library imports
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import os

#set pycharm display to show all columns of data
desired_width=320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)
pd.set_option('display.max_columns', None)

def viewshedcreate():
    '''
    Takes a raster DEM, returns a raster viewshed. Creates based on an observer point's X and Y location, and Y height
    '''
    ox=[274368] # observer X location in srs units
    oy=[433564] # observer Y location in srs units
    oz=[2] # observer height Y above DEM in metres
    src_ds_dem = gdal.Open('data/inputs/ni_dtm.tif') # DEM file location

    ds = gdal.ViewshedGenerate(
        src_ds_dem.GetRasterBand(1),
        "GTiff", #geotiff
        "data/outputs/viewshed_northcoast.tif", # viewshed raster output location
        ["INTERLEAVE=BAND"],
        ox[0],
        oy[0],
        oz[0],
        0,  # target height
        1,  # visible value
        0,  # invisible value
        0,  # out of range vlaue
        -1.0,  # no data value,
        0.85714,  # dfCurvCoeff
        gdal.GVM_Edge,
        0,  # maximum distance from observer to perform calculation
        heightMode=gdal.GVOT_NORMAL,
        options=["UNUSED=YES"],
    )
viewshedcreate()

def rast2poly():
    '''
    Takes viewshed raster, converts to a polygon and exports to shapefile
    '''
    src_ds_rast = gdal.Open('data/outputs/viewshed_northcoast.tif') # viewshed raster location
    srcband = src_ds_rast.GetRasterBand(1)
    dst_layername = 'layernametest'
    drv = ogr.GetDriverByName("ESRI Shapefile") # driver type, shapefile
    dst_ds = drv.CreateDataSource('data/outputs/viewshed_poly.shp') # viewshed polygon output location

    sp_ref = osr.SpatialReference()
    sp_ref.SetFromUserInput('EPSG:29902') # TM65 Irish National Grid

    dst_layer = dst_ds.CreateLayer(dst_layername, srs = sp_ref)

    fld = ogr.FieldDefn("visible", ogr.OFTInteger) # layer name for visible/not visible values
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("visible") # layer name for visible/not visible values

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)

    del src_ds_rast
    del dst_ds
rast2poly()

def vs_visible_select():
    '''
    Takes viewshed polygon shapefile, returns only the visible values i.e. those with 1, exports to shapefile
    '''
    vs_visible = gpd.read_file('data/outputs/viewshed_poly.shp') # viewshed polygon shapefile location
    vs_visible[vs_visible['visible']==1].to_file('data/outputs/viewshed_visible.shp') # select polygons with value of 1, exports to output location
vs_visible_select()

def lcm_clip():
    '''
    Takes the viewshed shapefile and clips landcover map shapefile
    '''
    vs_mask = gpd.read_file('data/outputs/viewshed_visible.shp') # viewshed shapefile location
    vs_mask_tm65 = vs_mask.to_crs(epsg=29902) # ensures it is set to TM65 Irish National Grid
    landcover = gpd.read_file('data/inputs/landcover.shp') # landcover map shapefile location
    landcover_tm65 = landcover.to_crs(epsg=29902) # ensures it is also set to TM65 Irish National Grid
    landcover_clipped = gpd.clip(landcover_tm65, vs_mask_tm65) # clips landcover by viewshed
    landcover_clipped.to_file('data/outputs/landcover_clipped.shp') # exports to output location
lcm_clip()

def lcm_stats():
    '''
    Finds area statistics of the clipped landcover map
    '''
    lcm = gpd.read_file('data/outputs/landcover_clipped.shp') # clipped landcover map location
    lcm = lcm.to_crs(epsg=29902) # ensures it is set to TM65 Irish National Grid
    lcm["lcm_area"] = lcm['geometry'].area #creates area geometry values for each attribute
    print(lcm.groupby(['category'])['lcm_area'].sum()/ 10**6) # prints list of landcover category area values in km2
    print(lcm.groupby(['landcover'])['lcm_area'].sum()/ 10**6) # prints list of landcover area values in km2
    print('Total Area (km2): {:.6f}'.format(lcm['lcm_area'].sum()/ 10**6)) # prints sum total of all area visible
lcm_stats()

def create_lcm_map():
    '''
    Creates a map plot of the visible landcover
    '''
    def generate_handles(labels, colors, edge='k', alpha=1):
        lc = len(colors)  # get the length of the color list
        handles = []
        for i in range(len(labels)):
            handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], edgecolor=edge, alpha=alpha))
        return handles


    #create scale bar
    def scale_bar(ax, location=(0.4, 0.05)): # location on plot
          x0, x1, y0, y1 = ax.get_extent()
          sbx = x0 + (x1 - x0) * location[0]
          sby = y0 + (y1 - y0) * location[1]
        #set scale bar segment positions
          ax.plot([sbx, sbx - 10000], [sby, sby], color='k', linewidth=9, transform=ax.projection)
          ax.plot([sbx, sbx - 5000], [sby, sby], color='k', linewidth=6, transform=ax.projection)
          ax.plot([sbx - 5000, sbx - 10000], [sby, sby], color='w', linewidth=6, transform=ax.projection)
        #set scale bar label positions
          ax.text(sbx, sby - 800, '20 km', transform=ax.projection, fontsize=8)
          ax.text(sbx - 5000, sby - 800, '10 km', transform=ax.projection, fontsize=8)
          ax.text(sbx - 10200, sby - 800, '0 km', transform=ax.projection, fontsize=8)

    outline = gpd.read_file(os.path.abspath('data/inputs/NI_outline.shp')) # loads outline shapefile
    lcm = gpd.read_file('data/outputs/landcover_clipped.shp') # loads landcover shapefile

    myFig = plt.figure(figsize=(15, 15)) # creates figure, 15 inches by 15 inches
    myCRS = ccrs.UTM(29) #set ref system to UTM section 29, for NI
    ax = plt.axes(projection=myCRS) #create axes object, set to the projection above


    outline_feature = ShapelyFeature(outline['geometry'], myCRS, edgecolor='k', facecolor='w', alpha=0.7, linewidth=0.7) # add outline to plot, set to above projection, black edge, white face, opacity 0.7, width 0.7

    xmin, ymin, xmax, ymax = lcm.total_bounds # sets total bounds of figure based on landcover
    ax.add_feature(outline_feature) # adds outline

    ax.set_extent([xmin-1000, xmax+1000, ymin-1000, ymax+1000], crs=myCRS)  # sets bounds +/-1000 of outline

    lcm_colors = ['seagreen', 'lime', 'gray', 'gold', 'orange', 'blue', 'saddlebrown'] # set landcover category colours, alphabetically

    lcm_cats = list(lcm.category.unique()) # get a list of unique names for the county boundaries
    lcm_cats.sort()  # sort the counties alphabetically by name

    for ii, name in enumerate(lcm_cats): # add landcover to map, iterates over landcover category
        feat = ShapelyFeature(lcm.loc[lcm['category'] == name, 'geometry'],  # first argument is the geometry
                              myCRS,  # second argument is the CRS
                              edgecolor='k',  # set outline to black
                              facecolor=lcm_colors[ii],  # set the face color to the corresponding colour from the list
                              linewidth=0.5,  # set line width
                              alpha=0.75)  # set opacity, 0 (transparent) to 1 (opaque)
        ax.add_feature(feat)  # add to map

    lcm_handles = generate_handles(lcm.category.unique(), lcm_colors, alpha=0.75) # generate handles for landcover categories, opacity set to 0.75

    # set landcover category names to begin with uppercase
    nice_names = [name.title() for name in lcm_cats]
    handles = lcm_handles
    labels = nice_names

    #set legend title, font sizes, location, frame, and opacity
    leg = ax.legend(handles, labels, title='Legend', title_fontsize=12,
                    fontsize=10, loc='upper left', frameon=True, framealpha=1)

    scale_bar(ax) #add scalebar to plot

    myFig.savefig('data/outputs/landcover_map.png', bbox_inches='tight', dpi=300) # Save map figure to desired location and filetype. bbox_inches='tight' crops the image to the axis. 300 dots per inch
create_lcm_map()

