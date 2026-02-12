import os
from osgeo import osr, ogr, gdal
import numpy as np
from . import geo_transform, geo_info
from datetime import datetime, timedelta, timezone


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


def get_space_oblique_mercator_proj4(inc_angle, ps_rev, asc_lon, false_e, false_n):
    proj4 = (
        "+proj=som +inc_angle=" + str(inc_angle) + " +ps_rev=" + str(ps_rev) + " +asc_lon=" + str(asc_lon) + " +lon0=" + str(250) +
        " +x_0=" + str(false_e) + " +y_0=" + str(false_n)
    )
    return proj4


def get_oblique_mercator_proj4(alpha, lonc, lat_0, false_e, false_n):
    proj4 = (
        "+proj=omerc +gamma=" + str(alpha) + " +lonc=" + str(lonc) + " +lat_0=" + str(lat_0) +
        " +x_0=" + str(false_e) + " +y_0=" + str(false_n)
    )
    return proj4


def get_orthographic_proj4(lat0, lon0, false_e, false_n):
    proj4 = (
        "+proj=ortho +lat_0=" + str(lat0) + " +lon_0=" + str(lon0) +
        " +x_0=" + str(false_e) + " +y_0=" + str(false_n)
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
    # if target_epsg == 4326:
    #     target_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

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

    loc_x = np.array([gt[0] + xorder_vector * gt[1] + yorder * gt[2]
                     for yorder in yorder_vector])
    loc_y = np.array([gt[3] + xorder_vector * gt[4] + yorder * gt[5]
                     for yorder in yorder_vector])
    cols = loc_x.shape[0]
    rows = loc_x.shape[1]
    loc_p = np.array([loc_x.reshape(cols * rows),
                     loc_y.reshape(cols * rows)]).T  # .tolist()
    latlon = transform.TransformPoints(loc_p)
    latlon = np.array(latlon)[:, 0:2].reshape(
        cols, rows, 2).transpose((2, 0, 1))
    return latlon


def get_extent_from_latlon(lat, lon, source_epsg, target_epsg):
    ext = []
    coord_transform = get_coord_transform_epsg(source_epsg, target_epsg)
    coord_array = np.array([lat.reshape(-1), lon.reshape(-1)]).T
    trans_array = np.array(
        coord_transform.TransformPoints(coord_array))[:, 0:2]

    xmin = np.min(trans_array[:, 0])
    ymin = np.min(trans_array[:, 1])
    xmax = np.max(trans_array[:, 0])
    ymax = np.max(trans_array[:, 1])

    ext = [xmin, ymin, xmax, ymax]
    return ext


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
    ext.append(
        min([corners[0][0], corners[1][0], corners[2][0], corners[3][0]]))
    ext.append(
        min([corners[0][1], corners[1][1], corners[2][1], corners[3][1]]))
    ext.append(
        max([corners[0][0], corners[1][0], corners[2][0], corners[3][0]]))
    ext.append(
        max([corners[0][1], corners[1][1], corners[2][1], corners[3][1]]))
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
    prop.append(osr.SpatialReference(
        wkt=raster.GetProjection()).GetAttrValue("AUTHORITY", 1))
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
    ltod = t_array + np.where(lon2 > 180, lon2 - 360,
                              np.where(lon2 < -180, lon2 + 360, lon2)) / 360 - (day - 1)
    return ltod


def calc_distances(lons1, lats1, lons2, lats2, a=6378137.0, b=6356752.314245):  # unit:m
    # R = 6373.0*1000
    lon1 = np.radians(lons1)
    lat1 = np.radians(lats1)
    lon2 = np.radians(lons2)
    lat2 = np.radians(lats2)
    dlon = lon2 - lon1
    dlon = np.where(dlon > np.pi, dlon-2*np.pi,
                    np.where(dlon < -np.pi, dlon+2*np.pi, dlon))
    dlat = lat2 - lat1
    alat = (lat2 + lat1)/2
    # a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * \
    #     np.cos(lat2) * np.sin(dlon / 2) ** 2
    # c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    # distance = R * c

    e2 = (a**2 - b**2) / (a**2)
    w = np.sqrt(1 - e2 * (np.sin(alat)**2))
    m = a * (1 - e2) / (w**3)  # 子午線曲率半径
    n = a / w                         # 卯酉線曲半径
    distance = np.sqrt((m * dlat)**2 + (n * dlon * np.cos(alat))**2)  # 距離計測

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
        [[-slon, clon, zero_arr], [-slat * clon, -slat *
                                   slon, clat], [clat * clon, clat * slon, slat]]
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

    if (proc == 1) & (type(x0).__module__ != "numpy"):
        proc = 3
    elif (proc == 1) & (len(x0.shape) == 0):
        proc = 3

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
        [[-slon, -slat * clon, clat * clon],
            [clon, -slat * slon, clat * slon], [zero_arr, clat, slat]]
    )  # .transpose([0, 1, 2])

    if proc == 1:
        sx, sy, sz = np.einsum("jik,ik->jk", R, xyz) + xyz0
    elif proc == 2:
        sx, sy, sz = R.dot(xyz) + xyz0
    else:
        sx, sy, sz = np.einsum("ji,ik->jk", R, xyz) + \
            np.array([xyz0 for i in range(sx.shape[0])]).T

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

def get_satellite_local_frame_basis(r1, v1, mode="LVLH"):
    """
    r1, v1: shape (3, n)  (each column is one sample)
    mode: LVLH/RTN or VNC/VNR
    Returns
      E: shape (3, 3, n)  where E[:, :, i] is the basis matrix for sample i
         (rows are ex/ey/ez for LVLH or eV/eN/eR for VNR)
    """

    if mode == "LVLH":
        R = r1 / np.linalg.norm(r1, axis=0)
        H = np.cross(r1, v1, axis=0)
        N = H / np.linalg.norm(H, axis=0)
        T = np.cross(N, R, axis=0)
        T /= np.linalg.norm(T, axis=0)
        N = np.cross(R, T, axis=0)
        N /= np.linalg.norm(N, axis=0)
        E = np.stack([T, N, R], axis=0)
    # if mode == "LVLH":
    #     n_r1 = np.linalg.norm(r1, axis=0)
    #     ez = r1 / n_r1 

    #     # v_perp = v1 - (v1·ez) ez
    #     dot_v_ez = np.sum(v1 * ez, axis=0)
    #     v_perp = v1 - ez * dot_v_ez

    #     n_vp = np.linalg.norm(v_perp, axis=0)
    #     ex = v_perp / n_vp

    #     ey = np.cross(ez, ex, axis=0)
    #     n_ey = np.linalg.norm(ey, axis=0) 
    #     ey = ey / n_ey 

    #     E = np.stack([ex, ey, ez], axis=0)

    elif mode == "VNC":
        nv = np.linalg.norm(v1, axis=0)
        eV = v1 / nv

        h = np.cross(r1, v1, axis=0)
        nh = np.linalg.norm(h, axis=0)
        eN = h / nh

        eR = np.cross(eV, eN, axis=0)
        nR = np.linalg.norm(eR, axis=0) 
        eR = eR / nR 

        # re-orthonormalize N to keep right-handedness
        eN = np.cross(eR, eV, axis=0)
        eN /= np.linalg.norm(eN, axis=0)

        E = np.stack([eV, eN, eR], axis=0) 

    return E

def transform_earth_to_satellite_local_frame(r1, v1, r2, mode, origin_at_satellite=True):
    """
    r1, v1, r2: array-like shape (3,)
      r1 = (3,n) satellite position in ECEF/ECI
      v1 = (3,n) satellite velocity in ECEF/ECI
      r2 = (3,n,m) target position in ECEF/ECI
    mode: LVLH/RTN or VNC/VNR
    origin:
      True  -> use d = r2 - r1 (typical "projection into satellite local frame")
      False -> use d = r2      (coordinates w.r.t. Earth center but expressed in the sat frame axes)
    Returns:
      p : (3,)  = (X,Y,Z) components in satellite local frame
      E : (3,3) = basis matrix whose rows are (ex, ey, ez)
    """
    E=get_satellite_local_frame_basis(r1,v1,mode)
    d = (r2 - r1[:, :, None]) if origin_at_satellite else r2
    p = np.einsum('abn,bnm->anm', E, d)
    return p, E

def transform_satellite_local_frame_to_earth(r1, v1, r2, mode, origin_at_satellite=True):
    """
    r1, v1, r2: array-like shape (3,)
      r1 = (3,n) satellite position in ECEF/ECI
      v1 = (3,n) satellite velocity in ECEF/ECI
      r2 = (3,n,m) target position in ECEF/ECI
    mode: LVLH/RTN or VNC/VNR
    origin:
      True  -> use d = r2 - r1 (typical "projection into satellite local frame")
      False -> use d = r2      (coordinates w.r.t. Earth center but expressed in the sat frame axes)
    Returns:
      p : (3,)  = (X,Y,Z) components in satellite local frame
      E : (3,3) = basis matrix whose rows are (ex, ey, ez)
    """
    E=geo_info.get_satellite_local_frame_basis(r1,v1,mode)
    Einv = np.transpose(E, (1, 0, 2)) 
    p = np.einsum('ijn,jnm->inm', Einv, r2)
    p = (r1[:, :, None] + p) if origin_at_satellite else p
    return p, Einv


def transform_ecef_to_eci(r, t):
    """
    ECEF -> ECI (GMST-only) using numpy.

    Parameters
    ----------
    r :  (3,N,M)
        ECEF position vectors [same unit]
    times :datetime, shape (N)
        UTC datetimes for each vector
    """
    N = r.shape[1]

    theta = _gmst_rad_from_datetimes(t)
    c = np.cos(theta).reshape(N, 1)
    s = np.sin(theta).reshape(N, 1)

    x, y, z = r[0], r[1], r[2]
    x_eci = c * x - s * y
    y_eci = s * x + c * y
    z_eci = z

    out = np.stack([x_eci, y_eci, z_eci], axis=0)
    return out


def transform_eci_to_ecef(r, times):
    """
    ECI -> ECEF (GMST-only) using numpy.

    Parameters
    ----------
    r :  (3,N,M)
        ECI position vectors [same unit]
    times :datetime, shape (N)
        UTC datetimes for each vector
    """
    N = r.shape[1]

    theta = _gmst_rad_from_datetimes(times)
    c = np.cos(theta).reshape(N, 1)
    s = np.sin(theta).reshape(N, 1)

    x, y, z = r[0], r[1], r[2]
    # ECEF = Rz(-theta) * ECI
    x_ecef = c * x + s * y
    y_ecef = -s * x + c * y
    z_ecef = z

    out = np.stack([x_ecef, y_ecef, z_ecef], axis=0)
    return out

def transform_ecef_to_eci_posvel(r, v, times, OMEGA_EARTH = 7.2921150e-5):
    """
    r: (3,N,M)  ECEF position
    v: (3,N,M)  ECEF velocity (time-derivative in ECEF)
    times: datetime array (N,)
    # Earth rotation rate rad/s

    return:
      r_eci: (3,N,M)
      v_eci: (3,N,M)
    """
    N = r.shape[1]

    # Earth rotation rate
    omega = np.array([0.0, 0.0, OMEGA_EARTH], dtype=np.float64).reshape(3, 1, 1)

    theta = _gmst_rad_from_datetimes(times)  # (N,)
    c = np.cos(theta).reshape(N, 1)       # (N,1)
    s = np.sin(theta).reshape(N, 1)       # (N,1)

    x, y, z = r[0], r[1], r[2]  # each (N,M)
    r_eci = np.empty_like(r)
    r_eci[0] = c * x - s * y
    r_eci[1] = s * x + c * y
    r_eci[2] = z

    vx, vy, vz = v[0], v[1], v[2]
    v_rot = np.empty_like(v)
    v_rot[0] = c * vx - s * vy
    v_rot[1] = s * vx + c * vy
    v_rot[2] = vz

    # v_eci = R v_ecef + omega × r_eci
    omega_cross_r = np.cross(omega, r_eci, axisa=0, axisb=0, axisc=0)  # (3,N,M)
    v_eci = v_rot + omega_cross_r

    return r_eci, v_eci


def transform_eci_to_ecef_posvel(r, v, times, OMEGA_EARTH = 7.2921150e-5):
    """
    r: (3,N,M)  ECI position
    v: (3,N,M)  ECI velocity (time-derivative in ECEF)
    times: datetime array (N,)
    # Earth rotation rate rad/s
    return:
      r_eci: (3,N,M)
      v_eci: (3,N,M)
    """
    N = r.shape[1]

      # rad/s
    omega = np.array([0.0, 0.0, OMEGA_EARTH], dtype=np.float64).reshape(3, 1, 1)

    theta = _gmst_rad_from_datetimes(times)  # (N,)
    c = np.cos(theta).reshape(N, 1)       # (N,1)
    s = np.sin(theta).reshape(N, 1)       # (N,1)

    x, y, z = r[0], r[1], r[2]  # each (N,M)
    r_ecef = np.empty_like(r)
    r_ecef[0] = c * x + s * y
    r_ecef[1] = - s * x + c * y
    r_ecef[2] = z


    # v_eci = R v_ecef + omega × r_eci
    omega_cross_r = np.cross(omega, r, axisa=0, axisb=0, axisc=0)  # (3,N,M)
    v_dum = v - omega_cross_r
    vx, vy, vz = v_dum[0], v_dum[1], v_dum[2]
    v_ecef = np.empty_like(v)
    v_ecef[0] = c * vx + s * vy
    v_ecef[1] = - s * vx + c * vy
    v_ecef[2] = vz


    return r_ecef, v_ecef


def intersect_ray_ellipsoid_ecef(sp_ecef, d_ecef, a=6378137.0, b=6356752.314245, eps=1e-12):
    """
    Ray-ellipsoid intersection in ECEF.

    Inputs
    ------
    sp_ecef : array-like, shape (3,N)
        Satellite position in ECEF [m]
    d_ecef : array-like, shape (3,N)
        Line-of-sight unit vector in ECEF
        Ray is r(t) = r_sat + t * d_hat, t>=0
    a, b : float
        Ellipsoid semi-major and semi-minor axes [m]
        Ellipsoid: (x^2+y^2)/a^2 + z^2/b^2 = 1

    Returns
    -------
    r_gnd_ecef : ndarray, shape (3,N)
        Intersection point on ellipsoid in ECEF [m]
    t : float
        Range parameter [m] along the ray
    """

    ax2 = a * a
    bz2 = b * b

    xs, ys, zs = sp_ecef
    dx, dy, dz = d_ecef

    A = (dx*dx + dy*dy) / ax2 + (dz*dz) / bz2
    B = 2.0 * ((xs*dx + ys*dy) / ax2 + (zs*dz) / bz2)
    C = (xs*xs + ys*ys) / ax2 + (zs*zs) / bz2 - 1.0

    D = B*B - 4.0*A*C

    sqrtD = np.sqrt(D)

    # # Two solutions
    t1 = (-B - sqrtD) / (2.0*A)
    t2 = (-B + sqrtD) / (2.0*A)

    # # We need the nearest intersection in front of the satellite (t>=0)
    ts= np.array([t1,t2])

    t = np.min(ts,axis=0)

    r_gnd = sp_ecef + t * d_ecef
    return r_gnd,t


def _to_unix_seconds_utc(dt):
    """datetime -> unix seconds (UTC). naiveはUTC扱い。"""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt.timestamp()

def _gmst_rad_from_datetimes(times):
    """
    GMST angle [rad] from datetime array (UTC).
    Vallado等でよく使われる近似式（UT1差・極運動・歳差章動は無視）。
    """
    times = np.asarray(times, dtype=object)
    unix = np.array([_to_unix_seconds_utc(t) for t in times], dtype=np.float64)

    # Unix epoch -> Julian Date
    JD = 2440587.5 + unix / 86400.0
    T = (JD - 2451545.0) / 36525.0

    # GMST in degrees
    gmst_deg = (280.46061837
                + 360.98564736629 * (JD - 2451545.0)
                + 0.000387933 * T**2
                - (T**3) / 38710000.0)

    theta = np.deg2rad(np.mod(gmst_deg, 360.0))
    return theta

def get_Rx(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def get_Ry(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def get_Rz(psi):
    c, s = np.cos(psi), np.sin(psi)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def get_transform_matrix(roll, pitch, yaw, order="ZYX"):
    """
    roll, pitch, yaw: radians
    order="ZYX" means: R = Rz(yaw) @ Ry(pitch) @ Rx(roll)
    This returns R_{b<-l} if angles are defined as body rotation relative to LVLH.
    """
    if order == "ZYX":
        return get_Rz(yaw) @ get_Ry(pitch) @ get_Rx(roll)
    elif order == "XYZ":
        return get_Rx(roll) @ get_Ry(pitch) @ get_Rz(yaw)
    else:
        raise ValueError("unsupported order")

def transform_ecef2eci2sat2sens(gp,sp,sv,ot,roll,pich,yaw):
    gp_eci=transform_ecef_to_eci(gp, ot)
    sp, sv =transform_ecef_to_eci_posvel(sp[:,:,None], sv[:,:,None], ot)
    sp=sp[:,:,0]
    sv=sv[:,:,0]
    gp_sat,E=transform_earth_to_satellite_local_frame(sp,sv,gp_eci,mode="LVLH")
    Rbl = get_transform_matrix(roll, pich, yaw, order="ZYX")
    gp_sens=np.einsum("ji,ikl->jkl", Rbl, gp_sat)
    return gp_sens

def transform_sens2sat2eci2ecef(gp_sens,sp,sv,ot,roll,pich,yaw,stype="eci"):
    Rbl = get_transform_matrix(roll, pich, yaw, order="ZYX")
    Rbl = np.linalg.inv(Rbl)
    gp_sat=np.einsum("ji,ikl->jkl", Rbl, gp_sens)
    if stype=="ecef":
        sp, sv =transform_ecef_to_eci_posvel(sp[:,:,None], sv[:,:,None], ot)
        sp=sp[:,:,0]
        sv=sv[:,:,0]
    gp_eci,E=transform_satellite_local_frame_to_earth(sp,sv,gp_sat,mode="LVLH")
    gp_ecef=transform_eci_to_ecef(gp_eci, ot)
    return gp_ecef

def calc_orbit_normal_from_inc_raan(inc_rad, raan_rad):
    si, ci = np.sin(inc_rad), np.cos(inc_rad)
    sO, cO = np.sin(raan_rad), np.cos(raan_rad)
    h = np.array([si * sO, -si * cO, ci], dtype=float)
    h = h/np.linalg.norm(h)
    return h

def simulate_orbit_nav_from_circular_eci(dtimes, sp0_eci, v_mean, r_orbit, inc_deg=98.8, raan_deg=0.0):
    """
    Inputs:
      times: array-like [s]
      t0: float [s]
      r0_eci: (3,) initial position at t0 in ECI [m] (should be consistent with the chosen plane)
      v_mean: mean speed [m/s]
      r_orbit: orbit radius [m]
      inc_deg: inclination [deg] (98 deg)
      raan_deg: RAAN Ω [deg] (choose 0 if you don't care)

    Returns:
      r_eci: (3, N) positions in ECI [m]
    """
    inc = np.deg2rad(inc_deg)
    raan = np.deg2rad(raan_deg)
    n = v_mean / r_orbit
    theta = n * dtimes  # (N,)
    h_hat = calc_orbit_normal_from_inc_raan(inc, raan)

    # Build in-plane basis from r0 and h_hat
    p_hat = sp0_eci/np.linalg.norm(sp0_eci,axis=0)
    q_hat = np.cross(h_hat, p_hat, axisa=0, axisb=0, axisc=0)
    q_hat = q_hat/np.linalg.norm(q_hat,axis=0)

    c = np.cos(theta)
    s = np.sin(theta)
    r_eci = r_orbit * (p_hat[:, None] * c[None, :] + q_hat[:, None] * s[None, :])
    v_eci = (r_orbit * n) * (-p_hat[:, None] * s[None, :] + q_hat[:, None] * c[None, :])
    return r_eci,v_eci


def get_raan_from_sp_and_inc(r0_eci, inc_rad, eps=1e-12):
    """
    Solve RAAN Ω such that the orbit plane with inclination inc passes through r0:
        h_hat(inc,Ω) · r0 = 0
    Returns two solutions (Ω1, Ω2) in radians.
    """
    x0, y0, z0 = np.asarray(r0_eci, dtype=float).reshape(3,)
    Rxy = np.hypot(x0, y0)

    si, ci = np.sin(inc_rad), np.cos(inc_rad)

    lam = np.arctan2(y0, x0)
    s = -(ci * z0) / (si * Rxy)
    s = np.clip(s, -1.0, 1.0)

    alpha = np.arcsin(s)
    Omega1 = lam + alpha
    Omega2 = lam + (np.pi - alpha)

    # Normalize to [-pi, pi)
    Omega1 = (Omega1 + np.pi) % (2*np.pi) - np.pi
    Omega2 = (Omega2 + np.pi) % (2*np.pi) - np.pi
    return Omega1, Omega2


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
    re = h * np.sin(theta_look) / \
        (np.sin(np.pi - theta_inc) - np.sin(theta_look))
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
    re = h * np.sin(theta_look) / \
        (np.sin(np.pi - theta_inc) - np.sin(theta_look))
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


def get_pixel_index_from_points(ds, coords, coord_sr, coord_type="EPSG"):
    # 地理情報を取得
    geo_transform = ds.GetGeoTransform()
    proj = ds.GetProjection()
    source_srs = osr.SpatialReference()
    if coord_type == "EPSG":
        source_srs.ImportFromEPSG(coord_sr)  # WGS84
    elif coord_type == "WKT":
        source_srs.ImportFromWkt(coord_sr)  # WGS84
    elif coord_type == "PROJ4":
        source_srs.ImportFromProj4(coord_sr)  # WGS84

    target_srs = osr.SpatialReference()
    target_srs.ImportFromWkt(proj)
    transform = osr.CoordinateTransformation(source_srs, target_srs)
    coords2 = np.array(transform.TransformPoints(coords))[:, 0:2].T

    px = (coords2[0]+0.5 - geo_transform[0]) / geo_transform[1]
    py = (coords2[1]+0.5 - geo_transform[3]) / geo_transform[5]
    idx_arr = np.array([py, px]).astype(np.int32)
    idx_arr[0] = np.where((idx_arr[0] < 0) | (
        idx_arr[0] > ds.RasterYSize), np.NaN, idx_arr[0])
    idx_arr[1] = np.where((idx_arr[1] < 0) | (
        idx_arr[1] > ds.RasterXSize), np.NaN, idx_arr[1])
    return idx_arr.T


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
        raster = gdal.Open(os.path.dirname(__file__) +
                           os.sep + "data/gshhg_2.3.7_shp_i_rasterize0.05.tif")
        raster_clip = geo_transform.clip_rectangle(
            raster, minX, minY, maxX, maxY)

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
