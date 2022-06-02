from osgeo import gdal, osr, ogr
import geopandas as gpd
from shapely.wkt import loads
from fiona.crs import from_epsg

def geom_shpl_from_ogr(geom_ogr):
    wkt=geom_ogr.ExportToWkt()
    return loads(wkt)
