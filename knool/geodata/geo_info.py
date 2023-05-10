import os
from osgeo import osr, ogr, gdal
import numpy as np
from . import geo_transform, geo_info


def get_stereographic_proj4(lat0, lon0, false_e, false_n, lat_ts):
    # lat_ts:origin of latitude, southern hemis.: -70, northern hemis.: 70
    proj4 = (
        "+proj=stere +lat_0="
        + str(lat0)
        + " +lat_ts="
        + str(lat_ts)
        + " +lon_0="
        + str(lon0)
        + " +x_0="
        + str(false_e)
        + " +y_0="
        + str(false_n)
    )
    return proj4


def get_orthographic_proj4(lat0, lon0, false_e, false_n):
    proj4 = (
        "+proj=ortho +lat_0=" + str(lat0) + " +lon_0=" + str(lon0) + " +x_0=" + str(false_e) + " +y_0=" + str(false_n)
    )
    return proj4


def get_shifted_wgs84_proj4(lon0):
    proj4 = "+proj=longlat +pm=" + str(lon0) + " +datum=WGS84 +no_defs"
    return proj4


def get_coord_transform_epsg(source_epsg, target_epsg):
    source_ref = osr.SpatialReference()
    target_ref = osr.SpatialReference()
    source_ref.ImportFromEPSG(source_epsg)
    target_ref.ImportFromEPSG(target_epsg)
    if source_epsg == 4326:
        source_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if target_epsg == 4326:
        target_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    return osr.CoordinateTransformation(source_ref, target_ref)


# 注意！！！！本結果を用いてtrans.transformpointを実施してinfが出力された場合xとyを入れ替える。gdal未修正箇所。今後アップデートされるかも
def get_coord_transform_wkt(source_wkt, target_wkt):
    source_ref = osr.SpatialReference()
    target_ref = osr.SpatialReference()
    source_ref.ImportFromWkt(source_wkt)
    target_ref.ImportFromWkt(target_wkt)
    return osr.CoordinateTransformation(source_ref, target_ref)


# 注意！！！！本結果を用いてtrans.transformpointを実施してinfが出力された場合xとyを入れ替える。gdal未修正箇所。今後アップデートされるかも
def get_coord_transform_proj4(source_proj4, target_proj4):
    source_ref = osr.SpatialReference()
    target_ref = osr.SpatialReference()
    source_ref.ImportFromProj4(source_proj4)
    target_ref.ImportFromProj4(target_proj4)
    return osr.CoordinateTransformation(source_ref, target_ref)


def get_latlons_from_raster(raster, interval=1):
    # create the new coordinate system
    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(raster.GetProjectionRef())
    new_cs = osr.SpatialReference()
    new_cs.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(old_cs, new_cs)

    width = raster.RasterXSize
    height = raster.RasterYSize
    gt = raster.GetGeoTransform()
    # minx = gt[0] + width * gt[1] + height * gt[2]
    # miny = gt[3] + width * gt[4] + height * gt[5]
    xorder_vector = np.arange(0, width, interval) + 0.5
    yorder_vector = np.arange(0, height, interval) + 0.5

    loc_x = np.array([gt[0] + xorder_vector * gt[1] + yorder * gt[2] for yorder in yorder_vector])
    loc_y = np.array([gt[3] + xorder_vector * gt[4] + yorder * gt[5] for yorder in yorder_vector])
    cols = loc_x.shape[0]
    rows = loc_x.shape[1]
    loc_p = np.array([loc_x.reshape(cols * rows), loc_y.reshape(cols * rows)]).T  # .tolist()
    latlon = transform.TransformPoints(loc_p)
    latlon = np.array(latlon)[:, 0:2].reshape(cols, rows, 2).transpose((2, 0, 1))
    return latlon


def get_latlonval_vectors(lat_array, lon_array, val_array, mask):
    cols, rows = lat_array.shape()
    lat_vector = lat_array.shape(cols * rows)
    lon_vector = lon_array.shape(cols * rows)
    val_vector = val_array.shape(cols * rows)
    return lat_vector, lon_vector, val_vector


def get_raster_extent(raster):
    ext = []
    gt = raster.GetGeoTransform()

    xmin = gt[0]
    ymin = gt[3] + (gt[5] * raster.RasterYSize)
    xmax = gt[0] + (gt[1] * raster.RasterXSize)
    ymax = gt[3]

    ext = [xmin, ymin, xmax, ymax]
    return ext


def get_vector_extent(vector_ds):
    xmin = 9.9e33
    ymin = 9.9e33
    xmax = -9.9e33
    ymax = -9.9e33

    layer = vector_ds.GetLayer(0)
    for feature in layer:
        geom = feature.GetGeometryRef()
        xmin0, xmax0, ymin0, ymax0 = geom.GetEnvelope()
        if xmin0 < xmin:
            xmin = xmin0
        if ymin0 < ymin:
            ymin = ymin0
        if xmax0 > xmax:
            xmax = xmax0
        if ymax0 > ymax:
            ymax = ymax0
    ext = [xmin, ymin, xmax, ymax]
    return ext


def get_extent_from_corners(corners):
    ext = []
    ext.append(min([corners[0][0], corners[1][0], corners[2][0], corners[3][0]]))
    ext.append(min([corners[0][1], corners[1][1], corners[2][1], corners[3][1]]))
    ext.append(max([corners[0][0], corners[1][0], corners[2][0], corners[3][0]]))
    ext.append(max([corners[0][1], corners[1][1], corners[2][1], corners[3][1]]))
    return ext


def get_raster_corners(raster):
    ext = []
    gt = raster.GetGeoTransform()
    xarr = [0, raster.RasterXSize]
    yarr = [0, raster.RasterYSize]

    for px in xarr:
        for py in yarr:
            x = gt[0] + (px * gt[1]) + (py * gt[2])
            y = gt[3] + (px * gt[4]) + (py * gt[5])
            ext.append([x, y])
        yarr.reverse()
    return ext


def get_property_from_raster_with_proj(raster):
    prop = []
    prop.append(raster.RasterXSize)
    prop.append(raster.RasterYSize)
    prop.append(raster.GetGeoTransform())
    prop.append(raster.GetProjection())
    prop.append(osr.SpatialReference(wkt=raster.GetProjection()).GetAttrValue("AUTHORITY", 1))
    prop.append(raster.GetDriver())
    return prop


def get_property_from_raster_with_gcps(raster):
    prop = []
    prop.append(raster.RasterXSize)
    prop.append(raster.RasterYSize)
    prop.append(raster.GetGCPSpatialRef())
    prop.append(raster.GetGCPs())
    prop.append(raster.GetDriver())

    if not prop[2]:
        # fmt: off
        prop[2] = """
                    GEOGCS["WGS 84",
                    DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                        AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
                    PRIMEM["Greenwich",0,
                        AUTHORITY["EPSG","8901"]],
                    UNIT["degree",0.01745329251994328,
                        AUTHORITY["EPSG","9122"]],
                        AUTHORITY["EPSG","4326"]]"""
        # fmt: on
    return prop


def get_local_time_of_day(t_array, lon_array, day, rot_offset=0):
    lon2 = lon_array - rot_offset
    ltod = t_array + np.where(lon2 > 180, lon2 - 360, np.where(lon2 < -180, lon2 + 360, lon2)) / 360 - (day - 1)
    return ltod


def calc_distances(lons1, lats1, lons2, lats2):  # unit:km
    R = 6373.0
    lon1 = np.radians(lons1)
    lat1 = np.radians(lats1)
    lon2 = np.radians(lons2)
    lat2 = np.radians(lats2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c
    return distance


def calc_prime_vertical_radius(gdlat0, h=0, a=6378137.0, b=6356752.314245):
    gdlat = np.float64(np.radians(gdlat0))
    return a**2 / np.sqrt((a * np.cos(gdlat)) ** 2 + (b * np.sin(gdlat)) ** 2)


def transform_lla_to_ecef(lat0, lon0, alt, a=6378137.0, b=6356752.314245):
    f = (a - b) / a

    n = calc_prime_vertical_radius(lat0, a=a, b=b)
    lon = np.float64(np.radians(lon0))
    lat = np.float64(np.radians(lat0))

    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    x = (n + alt) * cosLat * np.cos(lon)
    y = (n + alt) * cosLat * np.sin(lon)
    z = ((1 - f) ** 2 * n + alt) * sinLat
    return x, y, z


def transform_ecef_to_lla(x, y, z, a=6378137.0, b=6356752.314245):
    f = (a - b) / a

    e_sq = f * (2 - f)
    eps = e_sq / (1.0 - e_sq)
    p = np.sqrt(x * x + y * y)
    q = np.arctan2((z * a), (p * b))

    sin_q = np.sin(q)
    cos_q = np.cos(q)

    sin_q_3 = sin_q * sin_q * sin_q
    cos_q_3 = cos_q * cos_q * cos_q

    phi = np.arctan2((z + eps * b * sin_q_3), (p - e_sq * a * cos_q_3))
    lam = np.arctan2(y, x)

    v = a / np.sqrt(1.0 - e_sq * np.sin(phi) * np.sin(phi))
    h = (p / np.cos(phi)) - v

    lat = np.degrees(phi)
    lon = np.degrees(lam)

    return lat, lon, h


def transform_ecef_to_enu(x0, y0, z0, lat0, lon0, h0, x, y, z):
    # x0, y0, z0, lat0, lon0, h0は１次元配列or任意の値で、x,y,zとは大きさが異なっても良い。
    # 本来、x0,y0,z0とlat0,lon0,h0は情報が重複しているので、本関数内でどちらかを計算するのが望ましいが
    # 何回もfoward,backward計算する場合があるので、同じ計算をしないよう両方の情報を引数としている。
    #
    if type(x).__module__ != "numpy":
        proc = 2
    elif len(lat0.shape) == 0:
        proc = 2
    else:
        proc = 1

    lat = np.float64(np.radians(lat0))
    lon = np.float64(np.radians(lon0))

    xyz0 = np.array([x0, y0, z0])
    xyz = np.array([x, y, z])

    if len(xyz0.shape) == 1:
        dxyz = (xyz.T - xyz0).T
    else:
        dxyz = xyz - xyz0

    slon = np.sin(lon)
    clon = np.cos(lon)
    slat = np.sin(lat)
    clat = np.cos(lat)
    zero_arr = np.zeros(x0.shape)
    R = np.array(
        [[-slon, clon, zero_arr], [-slat * clon, -slat * slon, clat], [clat * clon, clat * slon, slat]]
    )  # .transpose([0, 1, 2])
    # print(R.shape)
    if proc == 1:
        sx, sy, sz = np.einsum("jik,ik->jk", R, dxyz)
    else:
        sx, sy, sz = R.dot(dxyz)
    return sx, sy, sz


def transform_enu_to_ecef(x0, y0, z0, lat0, lon0, h0, sx, sy, sz):
    # x0, y0, z0, lat0, lon0, h0は１次元配列or任意の値で、sx,sy,szとは大きさが異なっても良い。
    # 本来、x0,y0,z0とlat0,lon0,h0は情報が重複しているので、本関数内でどちらかを計算するのが望ましいが
    # 何回もfoward,backward計算する場合があるので、同じ計算をしないよう両方の情報を引数としている。

    if type(sx).__module__ != "numpy":
        proc = 2
    elif len(sx.shape) == 0:
        proc = 2
    else:
        proc = 1

    lat = np.float64(np.radians(lat0))
    lon = np.float64(np.radians(lon0))
    xyz0 = np.array([x0, y0, z0])
    xyz = np.array([sx, sy, sz])
    slon = np.sin(lon)
    clon = np.cos(lon)
    slat = np.sin(lat)
    clat = np.cos(lat)
    zero_arr = np.zeros(x0.shape)
    R = np.array(
        [[-slon, -slat * clon, clat * clon], [clon, -slat * slon, clat * slon], [zero_arr, clat, slat]]
    )  # .transpose([0, 1, 2])

    if proc == 1:
        sx, sy, sz = np.einsum("jik,ik->jk", R, xyz) + xyz0
    else:
        sx, sy, sz = R.dot(xyz) + xyz0

    return sx, sy, sz


def transform_lla_to_rotated_enu(lat0, lon0, eaz, lat, lon):  # eaz: Earth Azimuth
    # R = 6373000.0
    x0, y0, z0 = transform_lla_to_ecef(lat0, lon0, 0)  # , a=R, b=R)
    x, y, z = transform_lla_to_ecef(lat, lon, 0)  # , a=R, b=R)
    p, q, r = transform_ecef_to_enu(x0, y0, z0, lat0, lon0, 0, x, y, z)
    rad = np.radians(eaz)
    q2 = np.cos(rad) * q + np.sin(rad) * p
    p2 = np.cos(rad) * p - np.sin(rad) * q
    return p2, q2, r


def geodetic_to_geocentric_latitude(gdlat0, h=0, a=6378137.0, b=6356752.314245):
    gdlat = np.radians(gdlat0)
    f = (a - b) / a
    e2 = 2 * f - f**2
    n = calc_prime_vertical_radius(gdlat)
    gclat = np.arctan((n * (1 - f) ** 2 + h) / (n + h) * np.tan(gdlat))
    gclat = np.degrees(gclat)
    return gclat


def calc_slant_range1(gdlat, theta_inc, theta_look, h):
    # gdlat: geodetic latitude of observed point on earth surface
    # theta_inc: incidence angle
    # h: s/c height
    # theta_look: sensor look angle (obtained from offnadir angle and attitude data)
    ea = theta_inc - theta_look
    re = h * np.sin(theta_look) / (np.sin(np.pi - theta_inc) - np.sin(theta_look))
    # re=calc_prime_vertical_radius(gdlat,a=a,b=b)
    # srange = re * (np.sqrt(((h + re) / re) ** 2 - np.cos(np.pi / 2 - theta_inc) ** 2) - np.sin(np.pi / 2 - theta_inc))
    srange = re * np.sin(ea) / np.sin(theta_look)
    return srange


def calc_slant_range2(gdlat, theta_inc, theta_look, h):
    # gdlat: geodetic latitude of observed point on earth surface
    # theta_inc: incidence angle
    # h: s/c height
    # theta_look: sensor look angle (obtained from offnadir angle and attitude data)
    ea = theta_inc - theta_look
    re = h * np.sin(theta_look) / (np.sin(np.pi - theta_inc) - np.sin(theta_look))
    # re=calc_prime_vertical_radius(gdlat,a=a,b=b)
    # srange = re * (np.sqrt(((h + re) / re) ** 2 - np.cos(np.pi / 2 - theta_inc) ** 2) - np.sin(np.pi / 2 - theta_inc))
    srange = re * np.sin(ea) / np.sin(theta_look)
    return srange


def calc_line_buffer_point(lat0, lon0, h0, lat, lon, h, distance, ori="right"):
    R = 6373000.0
    x0, y0, z0 = transform_lla_to_ecef(lat0, lon0, h0, a=R, b=R)
    x, y, z = transform_lla_to_ecef(lat, lon, h, a=R, b=R)
    p, q, r = transform_ecef_to_enu(x0, y0, z0, lat0, lon0, h0, x, y, z)
    sign = 1
    if ori == "left":
        sign = -1

    theta = distance / R
    tan_th = np.tan(theta)
    distance2 = R * tan_th
    coef = 1 / (q * q / (p * p) + 1)
    q2 = np.sqrt(distance2**2 * coef)
    q2 = sign * q2
    p2 = -q * q2 / p
    if type(lat0).__module__ != "numpy":
        r2 = 0
    else:
        r2 = np.zeros(p2.shape)

    buff_ecef = transform_enu_to_ecef(x0, y0, z0, lat0, lon0, h0, p2, q2, r2)
    lat, lon, h = transform_ecef_to_lla(*buff_ecef, a=R, b=R)

    return lat, lon, h


def count_geometry_in_polygon(shape, polygon):
    poly = ogr.Open(polygon, 1)
    poly_lyr = poly.GetLayer(0)
    new_field = ogr.FieldDefn("Count", ogr.OFTReal)
    new_field.SetWidth(10)
    poly_lyr.CreateField(new_field)

    shp = ogr.Open(shape, 0)
    shp_lyr = shp.GetLayer(0)

    for feature in poly_lyr:
        ext = feature.GetGeometryRef()
        shp_lyr.SetSpatialFilter(ext)
        count = shp_lyr.GetFeatureCount()
        feature.SetField("Count", count)
        poly_lyr.SetFeature(feature)
        feature.Destroy()
        shp_lyr.SetSpatialFilter(None)

    poly.Destroy()
    shp.Destroy()


def calc_ocean_area_in_polygon(wkt):
    try:
        geometry = ogr.CreateGeometryFromWkt(wkt)
        lon, lat, _ = geometry.Centroid().GetPoint()
        minX, maxX, minY, maxY = geometry.GetEnvelope()
        raster = gdal.Open(os.path.dirname(__file__) + os.sep + "data/gshhg_2.3.7_shp_i_rasterize0.05.tif")
        raster_clip = geo_transform.clip_rectangle(raster, minX, minY, maxX, maxY)

        srs_pre = geo_info.get_shifted_wgs84_proj4(0)
        srs_af = geo_info.get_orthographic_proj4(lat, lon, 0, 0)
        trans = geo_info.get_coord_transform_proj4(srs_pre, srs_af)
        geometry.Transform(trans)
        minX, maxX, minY, maxY = geometry.GetEnvelope()

        raster_reproj = geo_transform.warp(
            raster_clip, "/vsimem/output.tif", srs_af, 5000, 5000, minX, minY, maxX, maxY, "near"
        )
        vector = geo_transform.get_vector_from_raster(raster_reproj)
        layer = vector.GetLayer()
        for feat in layer:
            test = feat.GetGeometryRef()
            geometry = geometry.Difference(test)

        ocean_area = geometry.Area() / 1000 / 1000
        # trans=geo_info.get_coord_transform_proj4(srs_af, srs_pre)
        # geometry.Transform(trans)
        # geo_io.make_vector_from_geom("../../test_data/test2.shp",geometry)
    except:
        ocean_area = None
    return ocean_area
