from osgeo import ogr


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


def create_polygons(latlon_list):  # or latlon_array [n_polygon, n_point, 2(latlon)]
    poly_list = []
    for polygon in latlon_list:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        ring.AddPoint(polygon[0][1], polygon[0][0])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        poly_list.append(poly)
    return poly_list


def create_multiPolygon(latlon_list):  # or latlon_array [n_polygon, n_point, 2(latlon)]
    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)
    for polygon in latlon_list:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        multipoly.AddGeometry(poly)
    return multipoly


def find_intersect_polygons(poly1, poly2_list):
    poly2_intersect0 = [poly2.Intersection(poly1) for poly2 in poly2_list]
    poly2_intersect_bool = [True if poly is not None and poly.Area() > 0 else False for poly in poly2_intersect0]
    poly2_intersect = [poly for poly in poly2_intersect0 if poly is not None and poly.Area() > 0]
    return poly2_intersect, poly2_intersect_bool
