#tool to create viewsheds from a DEM raster and perform analysis

# library imports
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from shapely.geometry import Point, LineString, Polygon
import sys
import os
import fiona #to support different drivers for geodataframe to file

#set pycharm display to show all columns of data
desired_width=320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)
pd.set_option('display.max_columns', None)

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

def lcm_stats():
    lcm = gpd.read_file('data/outputs/landcover_clipped.shp')
    lcm = lcm.to_crs(epsg=29902)
    lcm["lcm_area"] = lcm['geometry'].area #creates area values for each attribute
    print(list(lcm.columns))
    print(lcm.groupby(['category'])['lcm_area'].sum()/ 10**6)
    #maybe figure out a way to export this as a csv or something, just a quick summary of stuff

lcm_stats()

def point_select():
    pnt = gpd.read_file('data/inputs/buildings_binevenagh.shp')
    vs = gpd.read_file('data/outputs/viewshed_visible.shp')
    pnt = pnt.to_crs(epsg=29902)
    vs = vs.to_crs(epsg=29902)
    pnt_visible = gpd.sjoin(pnt, vs, how='inner', predicate='within')
    pnt_visible.to_file('data/outputs/buildings_visible.shp')
point_select()

def point_stats():
    pnt = gpd.read_file('data/outputs/buildings_visible.shp')
    pnt = pnt.to_crs(epsg=29902)
    print(pnt.groupby(['USE'])['USE'].count())
point_stats()

#map plotting
# def makeanicemap():
#     outline = gpd.read_file('data/inputs/NI_outline.shp')
#     outline = outline.to_crs(epsg=4326)
#     lcm = gpd.read_file('data/outputs/landcover_clipped.shp')
#     lcm = lcm.to_crs(epsg=4326)
#     pnt = gpd.read_file('data/outputs/buildings_visible.shp')
#     pnt = pnt.to_crs(epsg=4326)
#
#     myFig = plt.figure(figsize=(10, 10))
#     myCRS = ccrs.UTM(29)
#     ax = plt.axes(projection=myCRS)
#
#     outline_feature = ShapelyFeature(outline['geometry'], myCRS, edgecolor='k', facecolor='w')
#     xmin, ymin, xmax, ymax = outline.total_bounds
#     ax.add_feature(outline_feature)
#     ax.set_extent([xmin-5000, xmax+5000, ymin-5000, ymax+5000], crs=myCRS)
#     #myFig.savefig('maptest2.png', bbox_inches='tight', dpi=300)
# makeanicemap()

def generate_handles(labels, colors, edge='k', alpha=1):
    lc = len(colors)  # get the length of the color list
    handles = []
    for i in range(len(labels)):
        handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], edgecolor=edge, alpha=alpha))
    return handles


# create a scale bar of length 20 km in the upper right corner of the map
# adapted this question: https://stackoverflow.com/q/32333870
# answered by SO user Siyh: https://stackoverflow.com/a/35705477
def scale_bar(ax, location=(0.5, 0.2)):
    x0, x1, y0, y1 = ax.get_extent()
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    ax.plot([sbx, sbx - 10000], [sby, sby], color='k', linewidth=9, transform=ax.projection)
    ax.plot([sbx, sbx - 5000], [sby, sby], color='k', linewidth=6, transform=ax.projection)
    ax.plot([sbx-50000, sbx - 10000], [sby, sby], color='w', linewidth=6, transform=ax.projection)

    ax.text(sbx, sby-4500, '10 km', transform=ax.projection, fontsize=8)
    ax.text(sbx-12500, sby-4500, '5 km', transform=ax.projection, fontsize=8)
    ax.text(sbx-24500, sby-4500, '0 km', transform=ax.projection, fontsize=8)


# load the outline of Northern Ireland for a backdrop
outline = gpd.read_file(os.path.abspath('data/inputs/NI_outline.shp'))

# load the datasets
lcm = gpd.read_file('data/outputs/landcover_clipped.shp')
pnt = gpd.read_file('data/outputs/buildings_visible.shp')

# create a figure of size 10x10 (representing the page size in inches)
myFig = plt.figure(figsize=(10, 10))

myCRS = ccrs.UTM(29)  # create a Universal Transverse Mercator reference system to transform our data.
# NI is in UTM Zone 29, so we pass 29 to ccrs.UTM()

ax = plt.axes(projection=myCRS)  # finally, create an axes object in the figure, using a UTM projection,
# where we can actually plot our data.

# first, we just add the outline of Northern Ireland using cartopy's ShapelyFeature
outline_feature = ShapelyFeature(outline['geometry'], myCRS, edgecolor='k', facecolor='w', alpha=0.7, linewidth=0.7)

xmin, ymin, xmax, ymax = lcm.total_bounds
ax.add_feature(outline_feature)  # add the features we've created to the map.

# using the boundary of the shapefile features, zoom the map to our area of interest
ax.set_extent([xmin-1000, xmax+1000, ymin-1000, ymax+1000], crs=myCRS)  # because total_bounds
# gives output as xmin, ymin, xmax, ymax,
# but set_extent takes xmin, xmax, ymin, ymax, we re-order the coordinates here.

# pick colors, add features to the map
lcm_colors = ['firebrick', 'seagreen', 'royalblue', 'coral', 'violet', 'cornsilk', 'pink', 'lime']

# get a list of unique names for the county boundaries
lcm_cats = list(lcm.category.unique())
lcm_cats.sort()  # sort the counties alphabetically by name

# next, add the municipal outlines to the map using the colors that we've picked.
# here, we're iterating over the unique values in the 'CountyName' field.
# we're also setting the edge color to be black, with a line width of 0.5 pt.
# Feel free to experiment with different colors and line widths.
for ii, name in enumerate(lcm_cats):
    feat = ShapelyFeature(lcm.loc[lcm['category'] == name, 'geometry'],  # first argument is the geometry
                          myCRS,  # second argument is the CRS
                          edgecolor='k',  # outline the feature in black
                          facecolor=lcm_colors[ii],  # set the face color to the corresponding color from the list
                          linewidth=1,  # set the outline width to be 1 pt
                          alpha=0.25)  # set the alpha (transparency) to be 0.25 (out of 1)
    ax.add_feature(feat)  # once we have created the feature, we have to add it to the map using ax.add_feature()

# here, we're setting the edge color to be the same as the face color. Feel free to change this around,
# and experiment with different line widths.
# water_feat = ShapelyFeature(water['geometry'],  # first argument is the geometry
#                             myCRS,  # second argument is the CRS
#                             edgecolor='mediumblue',  # set the edgecolor to be mediumblue
#                             facecolor='mediumblue',  # set the facecolor to be mediumblue
#                             linewidth=1)  # set the outline width to be 1 pt
# ax.add_feature(water_feat)  # add the collection of features to the map

# river_feat = ShapelyFeature(rivers['geometry'],  # first argument is the geometry
#                             myCRS,  # second argument is the CRS
#                             edgecolor='royalblue',  # set the edgecolor to be royalblue
#                             linewidth=0.2)  # set the linewidth to be 0.2 pt
# ax.add_feature(river_feat)  # add the collection of features to the map

# # ShapelyFeature creates a polygon, so for point data we can just use ax.plot()
# town_handle = ax.plot(towns.geometry.x, towns.geometry.y, 's', color='0.5', ms=6, transform=myCRS)

# generate a list of handles for the county datasets
lcm_handles = generate_handles(lcm.category.unique(), lcm_colors, alpha=0.25)

# # note: if you change the color you use to display lakes, you'll want to change it here, too
# water_handle = generate_handles(['Lakes'], ['mediumblue'])
#
# # note: if you change the color you use to display rivers, you'll want to change it here, too
# river_handle = [mlines.Line2D([], [], color='royalblue')]  # have to make this a list

# update county_names to take it out of uppercase text
nice_names = [name.title() for name in lcm_cats]

# ax.legend() takes a list of handles and a list of labels corresponding to the objects you want to add to the legend
handles = lcm_handles #+ water_handle + river_handle + town_handle
labels = nice_names #+ ['Lakes', 'Rivers', 'Towns']

leg = ax.legend(handles, labels, title='Legend', title_fontsize=12,
                fontsize=10, loc='upper left', frameon=True, framealpha=1)

# add the text labels for the towns
# for ind, row in towns.iterrows():  # towns.iterrows() returns the index and row
#     x, y = row.geometry.x, row.geometry.y  # get the x,y location for each town
#     ax.text(x, y, row['TOWN_NAME'].title(), fontsize=8, transform=myCRS)  # use plt.text to place a label at x,y

# add the scale bar to the axis
scale_bar(ax)

# save the figure as map.png, cropped to the axis (bbox_inches='tight'), and a dpi of 300
myFig.savefig('data/outputs/maptest5.png', bbox_inches='tight', dpi=300)