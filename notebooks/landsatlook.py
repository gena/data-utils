# 2014.11.08 11:36:01 W. Europe Standard Time
#Embedded file name: landsatlook.py
import numpy as np
import requests, json
import ogr, osr
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from matplotlib.patches import Polygon
import os
import fiona
from shapely.geometry import shape
import datetime

def search_tiles(extent, max_cloud, start_date = None, end_date = None):
    if not start_date:
        start_date = '1980-01-01'
    if not end_date:
        t = datetime.datetime.now()
        end_date = '{:04d}-{:02d}-{:02d}'.format(t.year, t.month, t.day)
    print start_date + '...' + end_date
    u = 'http://landsatlook.usgs.gov/arcgis/rest/services/LandsatLook/ImageServer/query?where=(acquisitionDate%20>%3D%20date%27{5}%27%20AND%20%20acquisitionDate%20<%3D%20date%27{6}%27)%20AND%20(sensor+%3D+%27TM%27+OR+sensor+%3D+%27ETM%27+OR+sensor+%3D+%27LANDSAT_ETM%27+OR+sensor+%3D+%27OLI%27)+AND+%28cloudCover+<%3D+{4}%29&objectIds=&time=&geometry=%7B%22xmin%22%3A{0}%2C%22ymin%22%3A{1}%2C%22xmax%22%3A{2}%2C%22ymax%22%3A{3}%2C%22spatialReference%22%3A%7B%22wkid%22%3A102100%7D%7D&geometryType=esriGeometryEnvelope&inSR=102100&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=sceneID%2Csensor%2CacquisitionDate%2CdateUpdated%2Cpath%2Crow%2CPR%2CcloudCover%2CsunElevation%2CsunAzimuth%2CreceivingStation%2CsceneStartTime%2Cmonth%2Cyear%2COBJECTID%2CdayOfYear%2CdayOrNight%2CbrowseURL&returnGeometry=true&outSR=102100&returnIdsOnly=false&returnCountOnly=false&pixelSize=&orderByFields=year%2CdayOfYear&groupByFieldsForStatistics=&outStatistics=&f=json'
    xmin = extent[0]
    ymin = extent[1]
    xmax = extent[2]
    ymax = extent[3]
    u = u.format(xmin, ymin, xmax, ymax, max_cloud, start_date, end_date)
    r = requests.post(u)
    return json.loads(r.text)['features']


def download_tile_image(tile, extent, resolution = 15, path_prefix = ''):
    xmin = extent[0]
    ymin = extent[1]
    xmax = extent[2]
    ymax = extent[3]
    ratio = (xmax - xmin) / (ymax - ymin)
    w = (xmax - xmin) / resolution
    h = int(w / ratio)
    rasterid = tile['attributes']['OBJECTID']
    scene_id = tile['attributes']['sceneID']
    year = tile['attributes']['year']
    month = tile['attributes']['month']
    day = tile['attributes']['dayOfYear']
    cloud_cover = tile['attributes']['cloudCover']
    row = tile['attributes']['row']
    path = tile['attributes']['path']
    u = 'http://landsatlook.usgs.gov/arcgis/rest/services/LandsatLook/ImageServer/exportImage?f=image&format=jpgpng&bbox={0}%2C{1}%2C{2}%2C{3}&imageSR=102100&bboxSR=102100&size={4}%2C{5}&renderingRule=%7B%22rasterFunction%22%3A%22Stretch%22%2C%22rasterFunctionArguments%22%3A%7B%22StretchType%22%3A0%7D%2C%22variableName%22%3A%22Raster%22%7D&mosaicRule=%7B%22mosaicMethod%22%3A%22esriMosaicLockRaster%22%2C%22ascending%22%3Atrue%2C%22lockRasterIds%22%3A%5B{6}%5D%2C%22mosaicOperation%22%3A%22MT_FIRST%22%7D'
    u = u.format(xmin, ymin, xmax, ymax, w, h, rasterid)
    r = requests.post(u)
    name = path_prefix + 'landsat_' + str(year) + '_' + '{:04d}'.format(day) + '_' + '{:03d}'.format(cloud_cover) + '_row{:04d}'.format(row) + '_path{:04d}'.format(path)
    with open(name + '.json', 'w') as text_file:
        text_file.write(str(tile['attributes']))
    with open(name + '.jpw', 'w') as text_file:
        text_file.write(str(resolution) + '\n')
        text_file.write(str(0.0) + '\n')
        text_file.write(str(0.0) + '\n')
        text_file.write(str(-resolution) + '\n')
        text_file.write(str(xmin + resolution / 2.0) + '\n')
        text_file.write(str(ymax + resolution / 2.0) + '\n')
    f = open(name + '.jpg', 'wb')
    f.write(r.content)
    f.close()
    return name + '.jpg'


def transform_coord(coordinate, transform):
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(coordinate[0], coordinate[1])
    point.Transform(transform)
    return [point.GetX(), point.GetY()]


def transform_coord_geo_to_web(coordinate):
    return transform_coord(coordinate, create_transform_geo_to_web())


def create_transform_geo_to_web():
    inputEPSG = 4326
    outputEPSG = 3857
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outputEPSG)
    return osr.CoordinateTransformation(inSpatialRef, outSpatialRef)


def transform_coord_web_to_geo(coordinate):
    return transform_coord(coordinate, create_transform_web_to_geo())


def create_transform_web_to_geo():
    inputEPSG = 3857
    outputEPSG = 4326
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outputEPSG)
    return osr.CoordinateTransformation(inSpatialRef, outSpatialRef)


def plot_tiles(tiles, extent, latmin = 0, latmax = 0, lonmin = 0, lonmax = 0, latstep = 0, lonstep = 0):
    transform = create_transform_web_to_geo()
    fig, ax = plt.subplots(figsize=(20, 12))
    if latmin == 0 and latmax == 0:
        pt_geo = transform_coord(tiles[0]['geometry']['rings'][0][0], transform)
        latmin = pt_geo[1]
        latmax = pt_geo[1]
        lonmin = pt_geo[0]
        lonmax = pt_geo[0]
        for t in tiles:
            for pt in t['geometry']['rings'][0]:
                pt_geo = transform_coord(pt, transform)
                latmin = min(pt_geo[1], latmin)
                latmax = max(pt_geo[1], latmax)
                lonmin = min(pt_geo[0], lonmin)
                lonmax = max(pt_geo[0], lonmax)

        latstep = 0.5 * (latmax - latmin)
        lonstep = 0.5 * (lonmax - lonmin)
        latmin -= latstep * 4.0
        latmax += latstep * 4.0
        lonmin -= lonstep * 4.0
        lonmax += lonstep * 4.0
    m = Basemap(projection='merc', llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, lat_ts=latmin, resolution='c')
    m.drawcoastlines(color='#cccccc')
    m.drawparallels(np.arange(latmin, latmax, latstep), labels=[10, 0, 0, 0])
    m.drawmeridians(np.arange(m.lonmin, m.lonmax + 30, lonstep), labels=[0, 0, 0, 10])
    rings = []
    for t in tiles:
        verts = t['geometry']['rings']
        verts_merc = []
        for p in verts[0]:
            p_geo = transform_coord(p, transform)
            verts_merc.append(m(p_geo[0], p_geo[1]))
            rings.append(verts_merc)

    c = collections.PolyCollection(rings)
    c.set_alpha(0.01)
    ax.add_collection(c)
    xmin = extent[0]
    ymin = extent[1]
    xmax = extent[2]
    ymax = extent[3]
    s = [[xmin, ymin],
     [xmin, ymax],
     [xmax, ymax],
     [xmax, ymin]]
    s_merc = []
    for p in s:
        p_geo = transform_coord(p, transform)
        s_merc.append(m(p_geo[0], p_geo[1]))

    c = collections.PolyCollection([s_merc])
    c.set_color('r')
    c.set_alpha(0.5)
    ax.add_collection(c)
    plt.show()


def plot_point(ax, point, color = 'black'):
    pt = point.coords[0]
    ax.scatter([pt[0]], [pt[1]], marker='o', s=80, c=color, facecolor='red')


def plot_polygon(ax, poly, color = 'blue'):
    a = np.asarray(poly.exterior)
    ax.add_patch(Polygon(a, facecolor=color, alpha=0.3))
    ax.plot(a[:, 0], a[:, 1], color='black')


def plot_multipolygon(ax, geom, color = 'blue'):
    """ Can safely call with either Polygon or Multipolygon geometry
    """
    if geom.type == 'Polygon':
        plot_polygon(ax, geom, color)
    elif geom.type == 'MultiPolygon':
        for poly in geom.geoms:
            plot_polygon(ax, poly, color)


def plot_feature(ax, f, color = 'blue'):
    g = f['geometry']
    t = g['type']
    s = shape(g)
    if t == 'Polygon':
        plot_multipolygon(ax, s, color)
    if t == 'Point':
        plot_point(ax, s, color)
