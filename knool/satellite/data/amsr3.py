from osgeo import gdal, osr
import numpy as np
import os
from datetime import datetime as dt
from datetime import timedelta
from ...helpers.misc import import_config
from ...geodata import geo_io
from ...fortlib import pmw_processor

gdal.PushErrorHandler('CPLQuietErrorHandler')


class AMSR3_L1B:
    """
    Read AMSR3_L1B data. Just after obj=AMSR3_L1B(hdfpath), you need to perform obj.set_lat_range(latmin,latmax,freq).
    freq: 6G, 7G, 10uG, 10G, 18G, 23G, 36G, 89AG, 89BG, 165G, 183r7G, 183r3G
    """

    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.params = import_config(config_path=os.path.dirname(
            __file__) + os.sep + "conf/amsr3_params.yaml")

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR3_L1B"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2]
            self.subdsID[name] = i

        self.metadata = self.ds.GetMetadata_Dict()
        self.overlapscans = int(
            self.metadata["NC_GLOBAL#NumberOfScansOverlap"])
        self.startTime = dt.strptime(
            self.metadata["NC_GLOBAL#ObservationStartDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.endTime = dt.strptime(
            self.metadata["NC_GLOBAL#ObservationEndDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.clip_array = slice(None)

    def keys(self):
        print(self.subdsID)

    def _get_value_freq(self, freq, str_list1, str_list2):
        output = []
        for str1, str2 in zip(str_list1, str_list2):
            var_name = str1+freq[:-1]+str2
            ds = self.get_subds(self.subdsID[var_name])
            scale_factor = float(ds.GetMetadata_Dict()[
                                 var_name+"#scale_factor"])
            add_offset = float(ds.GetMetadata_Dict()[var_name+"#add_offset"])
            output.append(ds.ReadAsArray()[::-1][
                          self.clip_array]*scale_factor+add_offset)
        return output

    def _read_latlon(self, freq):
        str_list1 = ["Latitude_P", "Longitude_P"]
        str_list2 = ["", ""]
        return self._get_value_freq(freq, str_list1, str_list2)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_latlon(self):
        return self.lat, self.lon

    def set_lat_range(self, lamin, lamax, freq, overlap=False):
        """freq: 6G, 7G, 10uG, 10G, 18G, 23G, 36G, 89AG, 89BG, 165G, 183r7G, 183r3G"""
        if (freq == "6G") | (freq == "7G"):
            freq = "0"+freq
        self.freq = freq
        self.clip_array = slice(None)
        lat, lon = self._read_latlon(freq)
        lat_1d_min = np.min(lat, axis=1)

        print("original size=", lat.shape)
        ovs = self.overlapscans
        self.clip_array = (lat_1d_min >= lamin) & (lat_1d_min <= lamax)

        if (not overlap) & (ovs!=0):
            self.clip_array[:ovs] = False
            self.clip_array[-ovs:] = False

        lmax = self.clip_array.shape[0]
        bool1 = False
        for i in range(lmax):
            if (not bool1) & self.clip_array[i]:
                start = i + 1
                bool1 = True
            if bool1 & (not self.clip_array[i]):
                finish = i
                break
        if i == lmax - 1:
            finish = i
        et = self.endTime
        st = self.startTime
        self.startTime_split = (et - st) * \
            (start - ovs - 1) / (lmax - ovs * 2 - 1) + st
        self.endTime_split = (et - st) * (finish - ovs -
                                          1) / (lmax - ovs * 2 - 1) + st

        self.lat = lat[self.clip_array]
        self.lon = lon[self.clip_array]

        self.length, self.width = self.lat.shape

    def get_1d_datetime(self):
        length = self.length
        width = self.width
        stime = self.get_subds(self.subdsID["ScanTimeUTC"]).ReadAsArray()[::-1][
            self.clip_array]
        length = self.length
        width = self.width
        date_list = [dt(stime[i, 0], stime[i, 1], stime[i, 2], stime[i, 3], stime[i, 4], stime[i, 5],stime[i, 6]*1000)
                     for i in range(length)]
        time_array = np.array(date_list)
        return time_array

    def get_obs_time(self):
        """
        Export obs. time (jday-1 + hour/24+minute/24/60+second/60/60/24). Note the part of [jday-1]
        """
        length = self.length
        width = self.width

        stime = self.get_subds(self.subdsID["ScanTimeUTC"]).ReadAsArray()[::-1][
            self.clip_array]
        length = self.length
        width = self.width
        date_list = [dt(stime[i, 0], stime[i, 1], stime[i, 2], stime[i, 3], stime[i, 4], stime[i, 5], stime[i, 6]*1000)
                     for i in range(length)]
        tm0 = np.array(
            [
                float(dt.strftime("%j")) - 1 + dt.hour / 24 +
                dt.minute / 24 / 60 + dt.second / 24 / 60 / 60 + dt.microsecond / 24 / 60 / 60 /1.0E6
                for dt in date_list
            ]
        )
        # time_array = np.tile(tm0, (width, 1)).T

        tm0 = np.append(tm0, np.array([(tm0[-1] - tm0[-2]) + tm0[-1]]), axis=0)
        time_array = np.array([(tm0[1:] - tm0[0:-1])*150/360 * i / self.width + tm0[0:-1] for i in range(self.width)]).transpose(
            1, 0
        )

        return time_array

    def get_brightness_temperature(self):
        freq = self.freq
        if (freq == "165G") | (freq == "183r7G") | (freq == "183r3G"):
            str_list1 = ["Tb_Ch"]
            str_list2 = ["V"]
            output = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        else:
            str_list1 = ["Tb_Ch", "Tb_Ch"]
            str_list2 = ["H", "V"]
            output = self._get_value_freq(self.freq, str_list1, str_list2)
        return output

    def get_navigation_data(self):
        nav = self.get_subds(self.subdsID["NavigationData"]).ReadAsArray()[::-1][
            self.clip_array]
        return nav

    def get_attitude_data(self):
        att = self.get_subds(self.subdsID["AttitudeData"]).ReadAsArray()[::-1][
            self.clip_array]
        return att

    def get_satellite_position(self):
        """格納されている位置データはスキャン開始時刻のもの。それをスキャン方向に内挿するメソッド。"""
        sp0 = self.get_subds(self.subdsID["NavigationData"]).ReadAsArray()[::-1][
            self.clip_array, 0:3]
        sp0 = np.append(sp0, np.array([(sp0[-1] - sp0[-2]) + sp0[-1]]), axis=0)
        sp = np.array([(sp0[1:] - sp0[0:-1])*150/360 * i / self.width + sp0[0:-1] for i in range(self.width)]).transpose(
            2, 1, 0
        )
        # sp = np.tile(sp0, (self.width,1,1)).transpose([2,1,0])

        return sp

    def get_satellite_velocity(self):
        """格納されている位置データはスキャン開始時刻のもの。それをスキャン方向に内挿するメソッド。"""
        sv0 = self.get_subds(self.subdsID["NavigationData"]).ReadAsArray()[::-1][
            self.clip_array, 3:6]
        sv0 = np.append(sv0, np.array([(sv0[-1] - sv0[-2]) + sv0[-1]]), axis=0)
        sv = np.array([(sv0[1:] - sv0[0:-1])*150/360 * i / self.width + sv0[0:-1] for i in range(self.width)]).transpose(
            2, 1, 0
        )
        # sv = np.tile(sv0, (self.width,1,1)).transpose([2,1,0])

        return sv

    def get_earth_azimuth(self):
        str_list1 = ["EarthAzimuth_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def get_earth_incidence(self):
        str_list1 = ["EarthIncidence_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def get_sun_azimuth(self):
        str_list1 = ["SunAzimuth_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def get_sun_elevation(self):
        str_list1 = ["SunElevation_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def get_land_ratio(self):
        str_list1 = ["LandAreaPercent_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def get_mean_height(self):
        str_list1 = ["AreaMeanHeight_P"]
        str_list2 = [""]
        val = self._get_value_freq(self.freq, str_list1, str_list2)[0]
        return val

    def set_output_prop(self, gcp_x=20, gcp_y=10, lat=None, lon=None):
        if lat is None:
            lat, lon = self.get_latlon()
        length, width = lat.shape
        gcps = []
        for ai in np.linspace(0, length - 1, gcp_y):
            for aj in np.linspace(0, width - 1, gcp_x):
                i = int(ai)
                j = int(aj)
                gcps.append(gdal.GCP(float(lon[i][j]), float(
                    lat[i][j]), 0.0, j + 0.5, i + 0.5))

        source_ref = osr.SpatialReference()
        source_ref.ImportFromEPSG(4326)

        self.output_prop = [width, length, source_ref, gcps]

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(
            self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")




class AMSR3_L1R:
    """
    Read AMSR3_L1R data. Just after obj=AMSR3_L1R(hdfpath), you need to perform obj.set_lat_range(latmin,latmax,freq).
    freq: 6G, 7G, 10uG, 10G, 18G, 23G, 36G, 89AG, 89BG, 165G, 183r7G, 183r3G
    """

    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.params = import_config(config_path=os.path.dirname(
            __file__) + os.sep + "conf/amsr3_params.yaml")

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR3_L1R"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2]
            self.subdsID[name] = i

        self.metadata = self.ds.GetMetadata_Dict()
        self.overlapscans = int(
            self.metadata["NC_GLOBAL#NumberOfScansOverlap"])
        self.startTime = dt.strptime(
            self.metadata["NC_GLOBAL#ObservationStartDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.endTime = dt.strptime(
            self.metadata["NC_GLOBAL#ObservationEndDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.clip_array = slice(None)

    def keys(self):
        print(self.subdsID)

    def _get_value_freq(self, freq, str_list1, str_list2):
        output = []
        for str1, str2 in zip(str_list1, str_list2):
            var_name = str1+freq[:-1]+str2
            ds = self.get_subds(self.subdsID[var_name])
            scale_factor = float(ds.GetMetadata_Dict()[
                                 var_name+"#scale_factor"])
            add_offset = float(ds.GetMetadata_Dict()[var_name+"#add_offset"])
            output.append(ds.ReadAsArray()[::-1][
                          self.clip_array]*scale_factor+add_offset)
        return output

    def _get_value_freq2(self, freq, freq2, str_list1, str_list2):
        output = []
        for str1, str2 in zip(str_list1, str_list2):
            var_name = str1+freq[:-1]+"Ch"+freq2[:-1]+str2
            ds = self.get_subds(self.subdsID[var_name])
            scale_factor = float(ds.GetMetadata_Dict()[
                                 var_name+"#scale_factor"])
            add_offset = float(ds.GetMetadata_Dict()[var_name+"#add_offset"])
            output.append(ds.ReadAsArray()[::-1][
                          self.clip_array]*scale_factor+add_offset)
        return output

    def _read_latlon(self):
        freq="89G"
        str_list1 = ["Latitude_P", "Longitude_P"]
        str_list2 = ["o", "o"]
        return self._get_value_freq(freq, str_list1, str_list2)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_latlon(self):
        return self.lat, self.lon

    def set_lat_range(self, lamin, lamax, overlap=False):
        self.clip_array = slice(None)
        lat, lon = self._read_latlon()
        lat_1d_min = np.min(lat, axis=1)

        print("original size=", lat.shape)
        ovs = self.overlapscans
        self.clip_array = (lat_1d_min >= lamin) & (lat_1d_min <= lamax)

        if not overlap:
            self.clip_array[:ovs] = False
            self.clip_array[-ovs:] = False

        lmax = self.clip_array.shape[0]
        bool1 = False
        for i in range(lmax):
            if (not bool1) & self.clip_array[i]:
                start = i + 1
                bool1 = True
            if bool1 & (not self.clip_array[i]):
                finish = i
                break
        if i == lmax - 1:
            finish = i

        et = self.endTime
        st = self.startTime
        self.startTime_split = (et - st) * \
            (start - ovs - 1) / (lmax - ovs * 2 - 1) + st
        self.endTime_split = (et - st) * (finish - ovs -
                                          1) / (lmax - ovs * 2 - 1) + st

        self.lat = lat[self.clip_array]
        self.lon = lon[self.clip_array]

        self.length, self.width = self.lat.shape

    def get_obs_time(self):
        """
        Export obs. time (jday-1 + hour/24+minute/24/60+second/60/60/24). Note the part of [jday-1]
        """
        length = self.length
        width = self.width

        stime = self.get_subds(self.subdsID["ScanTimeUTC"]).ReadAsArray()[::-1][
            self.clip_array]
        length = self.length
        width = self.width
        date_list = [dt(stime[i, 0], stime[i, 1], stime[i, 2], stime[i, 3], stime[i, 4], stime[i, 5])
                     for i in range(length)]
        time_array = np.array(
            [
                float(dt.strftime("%j")) - 1 + dt.hour / 24 +
                dt.minute / 24 / 60 + dt.second / 24 / 60 / 60
                for dt in date_list
            ]
        )
        time_array = np.tile(time_array, (width, 1)).T

        return time_array

    def get_brightness_temperature(self,freq1,freq2):
        """freq1: FOV (6G, 7G, 10uG, 10G, 18G, 23G, 36G, 89AG, 89BG, 165G, 183r7G, 183r3G"""
        
        if (freq1 == "6G") | (freq1 == "7G"):
            freq1 = "0"+freq1
        
        if (freq2 == "6G") | (freq2 == "7G"):
            freq2 = "0"+freq2

        if (freq2 == "165G") | (freq2 == "183r7G") | (freq2 == "183r3G"):
            str_list1 = ["Tb_FOV"]
            str_list2 = ["V_P89o"]
            output = self._get_value_freq2(freq1, freq2, str_list1, str_list2)[0]
        else:
            str_list1 = ["Tb_FOV", "Tb_FOV"]
            str_list2 = ["H_P89o", "V_P89o"]
            output = self._get_value_freq2(freq1, freq2, str_list1, str_list2)
        return output

    def get_navigation_data(self):
        nav = self.get_subds(self.subdsID["NavigationData"]).ReadAsArray()[::-1][
            self.clip_array]
        return nav

    def get_satellite_position(self):
        """格納されている位置データはスキャン開始時刻のもの。それをスキャン方向に内挿するメソッド。"""
        sp0 = self.get_subds(self.subdsID["NavigationData"]).ReadAsArray()[::-1][
            self.clip_array, 0:3]
        sp0 = np.append(sp0, np.array([(sp0[-1] - sp0[-2]) + sp0[-1]]), axis=0)
        sp = np.array([(sp0[1:] - sp0[0:-1]) * i / self.width + sp0[0:-1] for i in range(self.width)]).transpose(
            2, 1, 0
        )
        return sp

    def get_earth_azimuth(self,freq):
        str_list1 = ["EarthAzimuth_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def get_earth_incidence(self,freq):
        str_list1 = ["EarthIncidence_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def get_sun_azimuth(self,freq):
        str_list1 = ["SunAzimuth_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def get_sun_elevation(self,freq):
        str_list1 = ["SunElevation_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def get_land_ratio(self,freq):
        str_list1 = ["LandAreaPercent_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def get_mean_height(self,freq):
        str_list1 = ["AreaMeanHeight_P"]
        str_list2 = [""]
        val = self._get_value_freq(freq, str_list1, str_list2)[0]
        return val

    def set_output_prop(self, gcp_x=20, gcp_y=10, lat=None, lon=None):
        if lat is None:
            lat, lon = self.get_latlon()
        length, width = lat.shape
        gcps = []
        for ai in np.linspace(0, length - 1, gcp_y):
            for aj in np.linspace(0, width - 1, gcp_x):
                i = int(ai)
                j = int(aj)
                gcps.append(gdal.GCP(float(lon[i][j]), float(
                    lat[i][j]), 0.0, j + 0.5, i + 0.5))

        source_ref = osr.SpatialReference()
        source_ref.ImportFromEPSG(4326)

        self.output_prop = [width, length, source_ref, gcps]

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(
            self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")


class Open:
    def __new__(cls, hdfpath, *args):
        product_name = os.path.basename(hdfpath).split("_")[2][0:8]
        print(product_name)
        if product_name == "S1BTBBGA":
            load_class = super().__new__(type(AMSR3_L1B(None)))
            load_class.__init__(hdfpath)
        else:
            print("the product is not supported")
        return load_class
