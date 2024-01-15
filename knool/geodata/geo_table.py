from osgeo import gdal, osr, ogr
import geopandas as gpd
from shapely.wkt import loads
from fiona.crs import from_epsg
from . import geo_geom


def geom_ogr_to_geom_shpl(geom_ogr):
    wkt = geom_ogr.ExportToWkt()
    return loads(wkt)


def latlon_to_point(x, lat_name, lon_name):
    point = geo_geom.create_point(x[lat_name], x[lon_name])
    point = geom_ogr_to_geom_shpl(point)
    return point
