import os
import numpy as np
import pandas as pd
from skyfield.api import load
from datetime import timedelta
from skyfield.api import utc
from ...geodata import geo_io, geo_info, geo_geom
from ...helpers.misc import import_config


class TLEsPerSat:
    def __init__(self, infile):
        self.tles = load.tle_file(infile)
        date_list = []
        num_list = []
        for i, tle in enumerate(self.tles):
            date_list.append(tle.epoch.utc_datetime())
            num_list.append(i)
        self.df0 = pd.DataFrame(num_list, columns=["num"]).assign(Date=pd.to_datetime(date_list)).set_index("Date")
        self.df = self.df0

    def set_period(self, sdate0, edate0):
        sdate = sdate0.replace(tzinfo=utc)
        edate = edate0.replace(tzinfo=utc)
        self.df = self.df0[(self.df0.index > sdate) & (self.df0.index < edate)]

    def _get_position(self, num, date_sf):
        geocentric = self.tles[num].at(date_sf)
        subpoint = geocentric.subpoint()
        lat = subpoint.latitude.degrees
        lon = subpoint.longitude.degrees
        ele = subpoint.elevation.m
        return lat, lon, ele

    def calc_position_at(self, date):
        # input:
        # date: datetime object -> e.g., dt(2018, 2, 1, 12, 15, 30, 2000, tzinfo=utc)
        df = self.df
        row1 = df.loc[df.index == df.index.unique()[df.index.unique().get_loc(pd.to_datetime(date), method="nearest")]]
        date1_sf = load.timescale().from_datetime(date)
        num = row1.num.values[0]
        return self._get_position(num, date1_sf)

    def calc_positions_between(self, fdate, ldate, interval):
        # input:
        # ldate, fdate: text -> e.g., 2011-07-01
        # interval: int (min) -> 30 (min)
        ds_date = pd.date_range(start=fdate, end=ldate, freq=str(interval) + "min", tz=utc)
        pos_list = []
        date_list = []
        for date in ds_date:
            lat, lon, ele = self.calc_position_at(date)
            pos_list.append([lat, lon])
            date_list.append(date)
        self.output = pos_list
        self.output_attribute = date_list
        self.output_type = "csv_list_line"
        return pos_list

    def calc_positions_faster_between(self, fdt, ldt, interval):
        mdt = fdt + (ldt - fdt) / 2
        mdt = mdt.replace(tzinfo=utc)

        df = self.df
        row1 = df.loc[df.index == df.index.unique()[df.index.unique().get_loc(pd.to_datetime(mdt), method="nearest")]]
        num = row1.num.values[0]

        date_ds = pd.date_range(start=fdt, end=ldt, freq=str(interval) + "min", tz=utc)
        date_list = date_ds.to_list()
        date_sf = load.timescale().from_datetimes(date_list)

        lat, lon, _ = self._get_position(num, date_sf)
        self.output = np.stack([lat, lon]).T
        self.output_attribute = date_list
        self.output_type = "csv_list_line"
        return self.output

    def calc_buff_positions_between(self, fdt, ldt, interval, distance, ori="right"):
        # input:
        # ldate, fdate: text -> e.g., 2011-07-01
        # interval: int (min) -> 30 (min)
        ldt = ldt + timedelta(minutes=interval)
        ds_date = pd.date_range(start=fdt, end=ldt, freq=str(interval) + "min", tz=utc)
        pos_list = []
        date_list = []
        for i, date in enumerate(ds_date):
            lat, lon, ele = self.calc_position_at(date)

            date2 = ds_date[i + 1]
            lat2, lon2, ele2 = self.calc_position_at(date2)

            lat3, lon3, ele3 = geo_info.calc_line_buffer_point(lat, lon, 0, lat2, lon2, 0, distance, ori)
            pos_list.append([lat3, lon3])
            date_list.append(date)
        self.output_type = "csv_list_line"
        self.output_attribute = date_list
        self.output = pos_list
        return pos_list

    def calc_buff_positions_faster_between(self, fdt, ldt, interval, distance, ori="right"):
        mdt = fdt + (ldt - fdt) / 2
        mdt = mdt.replace(tzinfo=utc)

        df = self.df
        row1 = df.loc[df.index == df.index.unique()[df.index.unique().get_loc(pd.to_datetime(mdt), method="nearest")]]
        num = row1.num.values[0]

        date_ds = pd.date_range(start=fdt, end=ldt + timedelta(minutes=interval), freq=str(interval) + "min", tz=utc)
        date_list = date_ds.to_list()
        date_sf = load.timescale().from_datetimes(date_list)

        lat, lon, _ = self._get_position(num, date_sf)

        lat2, lon2, _ = geo_info.calc_line_buffer_point(lat[0:-1], lon[0:-1], 0, lat[1:], lon[1:], 0, distance, ori)
        self.output = np.stack([lat2, lon2]).T
        self.output_attribute = date_list[0:-1]
        self.output_type = "csv_list_line"
        return self.output

    def _calc_buff_positions__left_right_between(self, fdate, ldate, interval, distance, ori="middle"):

        if ori == "middle":
            rpos_array = self.calc_buff_positions_faster_between(fdate, ldate, interval, distance / 2, ori="right")
            lpos_array = self.calc_buff_positions_faster_between(fdate, ldate, interval, distance / 2, ori="left")
        elif ori == "right":
            rpos_array = self.calc_buff_positions_faster_between(fdate, ldate, interval, distance, ori="right")
            lpos_array = self.calc_positions_faster_between(fdate, ldate, interval)
        elif ori == "left":
            lpos_array = self.calc_buff_positions_faster_between(fdate, ldate, interval, distance, ori="left")
            rpos_array = self.calc_positions_faster_between(fdate, ldate, interval)
        else:
            print("The mode is not currently supported.")

        return lpos_array, rpos_array

    def calc_buff_area_between(self, fdate, ldate, interval, distance, ori="middle"):
        lpos_array, rpos_array = self._calc_buff_positions__left_right_between(fdate, ldate, interval, distance, ori)

        result = np.concatenate([rpos_array, lpos_array[::-1]])
        self.output = result
        self.output_attribute = None
        self.output_type = "csv_list_polygon"
        return self.output

    def calc_scene_areas(self, fdate, ldate, interval, distance, pinterval, ori="middle"):
        lpos_array, rpos_array = self._calc_buff_positions__left_right_between(fdate, ldate, interval, distance, ori)

        row_num = rpos_array.shape[0]
        col_num = 2
        bit = 8
        pol_num = int((row_num - 1) / (pinterval - 1))

        lpos_sub = np.lib.index_tricks.as_strided(
            lpos_array.copy(), (row_num, pinterval, 1, col_num), (bit * col_num * (pinterval - 1), bit * col_num, 8, 8)
        )[0:pol_num]
        rpos_sub = np.lib.index_tricks.as_strided(
            rpos_array.copy(), (row_num, pinterval, 1, col_num), (bit * col_num * (pinterval - 1), bit * col_num, 8, 8)
        )[0:pol_num]
        rpos_sub2 = rpos_sub[:, ::-1, :, :]
        outlines_sub = np.concatenate([lpos_sub, rpos_sub2], axis=1)[:, :, 0, :]

        self.output = outlines_sub
        date_list = self.output_attribute
        stime_list = [date_obj for date_obj in date_list[: -(pinterval - 1) : (pinterval - 1)]]
        ftime_list = [date_obj for date_obj in date_list[(pinterval - 1) :: (pinterval - 1)]]
        self.output_attribute = {"StartTime": stime_list, "EndTime": ftime_list}
        self.output_type = "csv_list_polygons"
        return self.output

    def calc_intersect_scene_areas(
        self, fdate, ldate, interval, distance, pinterval, geom1_list0, trans_list, ori="middle"
    ):
        geom1_list = []
        for geom1 in geom1_list0:
            geom1_list.append(geom1.Clone())

        self.calc_scene_areas(fdate, ldate, interval, distance, pinterval, ori)
        keys = self.output_attribute.keys()
        geom2_llist = []
        attr2_dlist = []
        for geom, trans in zip(geom1_list, trans_list):
            lon, lat, _ = geom.Centroid().GetPoint()
            bool_list = (
                np.max(geo_info.calc_distances(self.output[:, :, 1], self.output[:, :, 0], lon, lat), axis=1) < 5000
            )
            geom2_list = geo_geom.create_polygons(self.output[bool_list])
            [geom.Transform(trans) for geom in geom2_list]
            attr2_dict = {}
            for key in keys:
                attr2_dict[key] = [attr for (attr, bool) in zip(self.output_attribute[key], bool_list) if bool]
            geom2_llist.append(geom2_list)
            attr2_dlist.append(attr2_dict)

        [geom1.Transform(trans) for (geom1, trans) in zip(geom1_list, trans_list)]

        result_geom_list = []
        result_geom2_list = []
        result_attr_list = []

        for geom1, geom2_list, attr2_dict, trans in zip(geom1_list, geom2_llist, attr2_dlist, trans_list):
            poly2_intersect_list, bool_list = geo_geom.find_intersect_polygons(geom1, geom2_list)
            geom2_sub_list = [geom for (geom, bool) in zip(geom2_list, bool_list) if bool]
            attr2_sub_dict = {}
            for key in keys:
                attr2_sub_dict[key] = [attr for (attr, bool) in zip(attr2_dict[key], bool_list) if bool]

            result_geom_list.append(poly2_intersect_list)
            result_geom2_list.append(geom2_sub_list)
            result_attr_list.append(attr2_sub_dict)

        self.output = result_geom_list
        self.output_attribute = result_attr_list
        self.output_type = "list_of_geom_list_polygons"
        # output-> geom[aoi][scene], time[aoi][Start or End][scene]
        return result_geom_list, result_geom2_list

    def transform(self, result, trans_list):
        for i, trans in enumerate(trans_list):
            [geom.Transform(trans) for geom in result[i]]
        self.output = result
        return result

    def export(self, infile, option=None):

        ext = os.path.splitext(infile)[::-1][0]

        if ext == ".shp":
            sign = 1
        elif ext == ".kml":
            sign = -1
        else:
            print("The file type is not supported currently.")
            return

        if self.output_type == "csv_list_line":
            geom = geo_geom.create_lineString([latlon[::sign] for latlon in self.output])
            geo_io.make_vector_from_geom(geom, None, outfile=infile)

        elif self.output_type == "csv_list_polygon":
            geom = geo_geom.create_polygon([latlon[::sign] for latlon in self.output])
            geo_io.make_vector_from_geom(geom, None, outfile=infile)

        elif self.output_type == "csv_list_polygons":
            geom_list = geo_geom.create_polygons([[latlon[::sign] for latlon in polygon] for polygon in self.output])
            attr = {}
            attr["StartTime"] = [
                date_obj.strftime("%Y-%m-%d %H:%M:%S") for date_obj in self.output_attribute["StartTime"]
            ]
            attr["EndTime"] = [date_obj.strftime("%Y-%m-%d %H:%M:%S") for date_obj in self.output_attribute["EndTime"]]
            geo_io.make_vector_from_geomList(geom_list, attr, outfile=infile)

        elif self.output_type == "list_of_geom_list_polygons":
            index = option
            geom_list = self.output[index]
            attr = {}
            attr["StartTime"] = [
                date_obj.strftime("%Y-%m-%d %H:%M:%S") for date_obj in self.output_attribute[index]["StartTime"]
            ]
            attr["EndTime"] = [
                date_obj.strftime("%Y-%m-%d %H:%M:%S") for date_obj in self.output_attribute[index]["EndTime"]
            ]
            if sign == -1:
                [geom.SwapXY() for geom in geom_list]
            geo_io.make_vector_from_geomList(geom_list, attr, outfile=infile)


class SatelliteCatalog:

    ITEM_NAME = [
        "OBJECT_NAME",
        "OBJECT_ID",
        "NORAD_CAT_ID",
        "OBJECT_TYPE",
        "OPS_STATUS_CODE",
        "OWNER",
        "LAUNCH_DATE",
        "LAUNCH_SITE",
        "DECAY_DATE",
        "PERIOD",
        "INCLINATION",
        "APOGEE",
        "PERIGEE",
        "RCS",
        "DATA_STATUS_CODE",
        "ORBIT_CENTER",
        "ORBIT_TYPE",
    ]

    def __init__(self, infile=True):
        if infile:
            infile = os.path.dirname(__file__) + os.sep + "data/two_line_elements/satcat.csv"
        self.df = pd.read_csv(infile)  # , index_col=0

    def get_lines(self, item, val):
        return self.df[self.df[item] == val]

    def get_ncid_from_name(self, name):
        return self.get_lines("OBJECT_NAME", name)["NORAD_CAT_ID"].values[0]


class Sensor:
    def __init__(self, sensor_name):
        self.params = import_config(
            config_path=os.path.dirname(__file__) + os.sep + "conf/" + sensor_name + "_params.yaml"
        )
        self.param_keys = self.params.keys()
        satcat = SatelliteCatalog()
        self.ncid = satcat.get_ncid_from_name(self.params["satellite"])
        tle_path = (
            os.path.dirname(__file__)
            + os.sep
            + "data"
            + os.sep
            + "two_line_elements"
            + os.sep
            + "sat"
            + str(self.ncid)
            + ".txt"
        )
        self.sat = TLEsPerSat(tle_path)

    def union_adjacent_scenes(self, geom_llist):
        attr_dlist = self.sat.output_attribute
        geom_result_llist = []
        attr_result_dlist = []
        for geom_list, attr_dict in zip(geom_llist, attr_dlist):
            attr1_list = attr_dict["StartTime"]
            attr2_list = attr_dict["EndTime"]
            i = 0
            for i in range(len(geom_list) - 1):
                if attr2_list[i] == attr1_list[i + 1]:
                    geom_list[i + 1] = geom_list[i].Union(geom_list[i + 1])
                    attr1_list[i + 1] == attr1_list[i]
                    geom_list[i] = 0
            attr_dict["StartTime"] = [attr1 for (attr1, geom) in zip(attr1_list, geom_list) if geom != 0]
            attr_dict["EndTime"] = [attr2 for (attr2, geom) in zip(attr2_list, geom_list) if geom != 0]
            geom_list = [geom for geom in geom_list if geom != 0]
            geom_result_llist.append(geom_list)
            attr_result_dlist.append(attr_dict)
            self.sat.output = geom_result_llist
            self.sat.output_attribute = attr_result_dlist
        return geom_result_llist, attr_result_dlist

    def calc_intersect_scene_areas(self, stime, etime, interval, pinterval, aoi_list, trans_list, scene_all=False):
        ori = self.params["obs_direction"]
        obs_width = self.params["obs_width"]
        result1, result2 = self.sat.calc_intersect_scene_areas(
            stime, etime, interval, obs_width, pinterval, aoi_list, trans_list, ori
        )
        if scene_all:
            result = result2
            self.sat.output = result2
        #            result,attribute=self.union_adjacent_scenes(result2)
        else:
            result = result1
        #            result,attribute=self.union_adjacent_scenes(result1)
        return result, self.sat.output_attribute

    def transform(self, result, trans_list):
        result = self.sat.transform(result, trans_list)
        return result

    def export(self, infile, option):
        self.sat.export(infile, option)
