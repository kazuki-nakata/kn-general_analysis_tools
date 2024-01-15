import os
from osgeo import ogr
from datetime import datetime as dt
from datetime import timedelta
import pandas as pd
from . import satellite
from ...geodata import geo_info, geo_table
from skyfield.api import utc
import geopandas as gpd


def matchup(sensor1, sensor2, stime, etime, interval, pinterval, diff_min, aoi_list, trans_for, trans_back):
    # stime,etime -> datetime object. It is noted that calculation period exeeds the period from stime to etime.
    # About output data
    # sensor1 #[aoi][scene] #[aoi][start or end][scene]
    # sensor2 #[aoi][sen1_scene][sen2_scene] #[aoi][sen1_scene][start or end][sen2_scene]

    geom_s1_llist00, _ = sensor1.calc_intersect_scene_areas(
        stime, etime, interval, pinterval, aoi_list, trans_for)
    geom_s1_llist0, attr_s1_dlist = sensor1.union_adjacent_scenes(
        geom_s1_llist00)
    geom_s1_llist = sensor1.transform(geom_s1_llist0, trans_back)

    print("The number of polygon Before Union (Only first AOI):",
          len(geom_s1_llist00[0]))
    print("The number of polygon After Union (Only first AOI)::",
          len(geom_s1_llist[0]))
    #    sensor1.export(r"../../test_data/output/final1.shp",option=0)
    #    sensor1.export(r"../../test_data/output/final2.shp",option=0)
    srs_pre = geo_info.get_shifted_wgs84_proj4(0)
    time_none = {"StartTime": [0], "EndTime": [0]}
    geom_s2_llist = []
    attr_s2_dlist = []
    j = 0
    # iterate for regions
    for attr_s1_dict, geom_s1_list in zip(attr_s1_dlist, geom_s1_llist):
        geom_s2_list = []
        attr_s2_list = []
        j += 1
        print("aoi number:", j)
        # iterate for results in a region during a period
        for i, geom_s1 in enumerate(geom_s1_list):

            stime = attr_s1_dict["StartTime"][i] - timedelta(minutes=diff_min)
            etime = attr_s1_dict["StartTime"][i] + timedelta(minutes=diff_min)
            print(stime, etime)

            lon, lat, _ = geom_s1.Centroid().GetPoint()
            srs_af = geo_info.get_orthographic_proj4(lat, lon, 0, 0)
            trans_for = geo_info.get_coord_transform_proj4(srs_pre, srs_af)
            trans_back = geo_info.get_coord_transform_proj4(srs_af, srs_pre)

            result00, _ = sensor2.calc_intersect_scene_areas(
                stime, etime, interval, pinterval, [
                    geom_s1], [trans_for], scene_all=True
            )
            if len(result00[0]) > 0:
                result0, attr = sensor2.union_adjacent_scenes(result00)
                result = sensor2.transform(result0, [trans_back])
                geom_s2_list.append(result[0])
                attr_s2_list.append(attr[0])
            else:
                geom_s2_list.append([0])
                attr_s2_list.append(time_none)

        geom_s2_llist.append(geom_s2_list)
        attr_s2_dlist.append(attr_s2_list)
    return geom_s1_llist, attr_s1_dlist, geom_s2_llist, attr_s2_dlist


def make_gdf_from_matchup(geom_s1, attr_s1, geom_s2, attr_s2, aoi_names, limit, outfile, outformat):
    aoi_n = len(aoi_names)

    s2_column = []
    s2_column.extend(["StartTime" + str(i) for i in range(1, limit + 1)])
    s2_column.extend(["EndTime" + str(i) for i in range(1, limit + 1)])
    s2_column.extend(["geometry" + str(i) for i in range(1, limit + 1)])

    gdf_result = [0] * aoi_n
    for aoi_i in range(aoi_n):
        ams_n = len(geom_s1[aoi_i])
        stime = attr_s1[aoi_i]["StartTime"]
        etime = attr_s1[aoi_i]["EndTime"]
        geom_shpl = list(map(geo_table.geom_ogr_to_geom_shpl, geom_s1[aoi_i]))
        aoi = [aoi_names[aoi_i]] * ams_n
        gdf_amsr2 = gpd.GeoDataFrame(
            data=list(zip(aoi, stime, etime, geom_shpl)),
            columns=["AOI", "StartTime", "EndTime", "geometry"],
            crs="EPSG:4326",
        )

        s2_list = []
        for ams_i in range(ams_n):
            res_test = [None] * limit * 3
            i = 0
            for val1, val2, val3 in zip(
                attr_s2[aoi_i][ams_i]["StartTime"], attr_s2[aoi_i][ams_i]["EndTime"], geom_s2[aoi_i][ams_i]
            ):
                if val1 != 0:
                    res_test[i] = val1
                    res_test[i + limit] = val2
                    # geo_table.geom_ogr_to_geom_shpl(val3)
                    res_test[i + limit * 2] = val3.ExportToIsoWkt()
                else:
                    res_test[i] = None
                    res_test[i + limit] = None
                    res_test[i + limit * 2] = None

                i += 1
            s2_list.append(res_test)
        gdf_modis = gpd.GeoDataFrame(data=s2_list, columns=s2_column)
        gdf_result[aoi_i] = pd.concat([gdf_amsr2, gdf_modis], axis=1)

    gdf_all = pd.concat([gdf_result[i] for i in range(aoi_n)])
    gdf_all.to_file(outfile, driver=outformat)

    date_column = ["StartTime", "EndTime"]
    date_column.extend(["StartTime" + str(i) for i in range(1, limit + 1)])
    date_column.extend(["EndTime" + str(i) for i in range(1, limit + 1)])

    for datetime in date_column:
        gdf_all[datetime] = pd.to_datetime(gdf_all[datetime])
    gdf_all = gdf_all.reset_index().set_index("StartTime")


# fmt: off
def matchup_AMSR2_MODIS(stime, etime, season, shpfile, interval=5, pinterval=2, diff_min=30,
                        split_period="7d", outdir=None, out_interval=None, aoi_names=None):
    # stime -> e.g., "2012-01-01"
    # etime -> e.g., "2012-02-01"
    # season -> e.g., [9,5]
    # pinterval >= 2 (unit:point)
    # interval (unit:minute)
    vector = ogr.Open(shpfile)
    layer = vector.GetLayer()
    geom_list = []
    for feat in layer:
        geom_list.append(feat.GetGeometryRef().Clone())

    aoi_n = len(geom_list)
    sensor_name = "amsr2"
    sensor1 = satellite.Sensor(sensor_name)
    sensor_name = "modis"
    sensor2 = satellite.Sensor(sensor_name)

    trans_for = []
    trans_back = []
    for geom in geom_list:
        lon, lat, _ = geom.Centroid().GetPoint()
        srs_pre = geo_info.get_shifted_wgs84_proj4(0)
        srs_af = geo_info.get_orthographic_proj4(lat, lon, 0, 0)
        trans_for.append(geo_info.get_coord_transform_proj4(srs_pre, srs_af))
        trans_back.append(geo_info.get_coord_transform_proj4(srs_af, srs_pre))

    smonth = season[0]
    emonth = season[1]
    date_ds0 = pd.date_range(start=stime, end=etime, freq=split_period, tz=utc)

    if smonth < emonth:
        date_ds = date_ds0[(date_ds0.month >= smonth) & (date_ds0.month <= emonth)]
    else:
        date_ds = date_ds0[(date_ds0.month <= smonth) | (date_ds0.month >= emonth)]

    date_df = pd.DataFrame(zip(date_ds, date_ds + timedelta(weeks=1)), columns=["startDate", "endDate"])

    if outdir is not None:
        i = 0
        str_stime = pd.to_datetime(date_df.head(1)["startDate"].values[0]).strftime("%Y%m%d")
        os.makedirs(outdir, exist_ok=True)

    for index, row in date_df.iterrows():  # iterate for periods

        print("Processing...")
        print(row)

        os1_geom, os1_attr, os2_geom, os2_attr = matchup(
            sensor1, sensor2, row["startDate"], row["endDate"], interval, pinterval,
            diff_min, geom_list, trans_for, trans_back)

        if outdir is None:
            if index == 0:
                out_s1_geom = os1_geom
                out_s1_attr = os1_attr
                out_s2_geom = os2_geom
                out_s2_attr = os2_attr
            else:
                [out_s1_geom[j].extend(os1_geom[j]) for j in range(aoi_n)]
                [out_s1_attr[j]["StartTime"].extend(os1_attr[j]["StartTime"]) for j in range(aoi_n)]
                [out_s1_attr[j]["EndTime"].extend(os1_attr[j]["EndTime"]) for j in range(aoi_n)]

                [out_s2_geom[j].extend(os2_geom[j]) for j in range(aoi_n)]
                [out_s2_attr[j].extend(os2_attr[j]) for j in range(aoi_n)]

        else:
            print(i, i % out_interval)
            if i % out_interval == 0:
                out_s1_geom_sub = os1_geom
                out_s1_attr_sub = os1_attr
                out_s2_geom_sub = os2_geom
                out_s2_attr_sub = os2_attr
            else:
                [out_s1_geom_sub[j].extend(os1_geom[j]) for j in range(aoi_n)]
                [out_s1_attr_sub[j]["StartTime"].extend(os1_attr[j]["StartTime"]) for j in range(aoi_n)]
                [out_s1_attr_sub[j]["EndTime"].extend(os1_attr[j]["EndTime"]) for j in range(aoi_n)]
                [out_s2_geom_sub[j].extend(os2_geom[j]) for j in range(aoi_n)]
                [out_s2_attr_sub[j].extend(os2_attr[j]) for j in range(aoi_n)]
            i += 1
            if i % out_interval == 0:
                filename = str_stime + "_" + dt.strftime(row["endDate"], "%Y%m%d") + ".gpkg"
                outfile = os.path.join(outdir, filename)
                print("Export results as a geopackage (filename:", outfile)
                print("scene number for first geometry:", len(out_s1_geom_sub[0]))
                make_gdf_from_matchup(
                    out_s1_geom_sub, out_s1_attr_sub, out_s2_geom_sub, out_s2_attr_sub, aoi_names, 2, outfile, "GPKG"
                )
                print("finish")
                str_stime = dt.strftime(row["endDate"], "%Y%m%d")
                out_s1_geom_sub = []
                out_s1_attr_sub = []
                out_s2_geom_sub = []
                out_s2_attr_sub = []

    if outdir is None:
        return out_s1_geom, out_s1_attr, out_s2_geom, out_s2_attr
# fmt: on


def get_scene_name_MODIS(date0):
    if date0 is None:
        return None
    minute0 = int(date0.minute / 5) * 5
    date = dt(date0.year, date0.month, date0.day, date0.hour, minute0)
    ver = "061"
    year = date.year
    jday = (date - dt(year, 1, 1)).days + 1
    hour = date.hour
    minute = date.minute
    name = str(year) + str(jday).zfill(3) + "." + \
        str(hour).zfill(2) + str(minute).zfill(2) + "." + ver
    #    file_list = ["MYD021KM.A" + naming_rule, "MYD03.A" + naming_rule, "MYD35_L2.A" + naming_rule]
    return name
