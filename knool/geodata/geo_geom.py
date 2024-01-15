from osgeo import ogr


def create_point(lat, lon):
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(lon, lat)
    return point


def create_points(latlon_list):  # or latlon_array
    point_list = []
    for latlon in latlon_list:
        point = create_point(latlon[1], latlon[0])
        point_list.append(point)
    return point_list


def create_lineString(latlon_list):  # or latlon_array
    line = ogr.Geometry(ogr.wkbLineString)
    for latlon in latlon_list:
        line.AddPoint(latlon[1], latlon[0])
    return line


def create_polygon(latlon_list):  # or_latlon_array
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for latlon in latlon_list:
        ring.AddPoint(latlon[1], latlon[0])
    ring.AddPoint(latlon_list[0][1], latlon_list[0][0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


# or latlon_array [n_polygon, n_point, 2(latlon)]
def create_polygons(latlon_list):
    poly_list = []
    for polygon in latlon_list:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        ring.AddPoint(polygon[0][1], polygon[0][0])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        poly_list.append(poly)
    return poly_list


# or latlon_array [n_polygon, n_point, 2(latlon)]
def create_multiPolygon(latlon_list):
    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)
    for polygon in latlon_list:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        multipoly.AddGeometry(poly)
    return multipoly


def find_intersect_polygons(geom1, geom2_list):
    geom2_intersect0 = [geom2.Intersection(geom1) for geom2 in geom2_list]
    geom2_intersect_bool = [True if geom is not None and geom.Area(
    ) > 0 else False for geom in geom2_intersect0]
    geom2_intersect = [
        geom for geom in geom2_intersect0 if geom is not None and geom.Area() > 0]
    return geom2_intersect, geom2_intersect_bool
