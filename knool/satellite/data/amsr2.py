from osgeo import gdal, osr
import numpy as np
import os
from datetime import datetime as dt
from datetime import timedelta
from ...helpers.misc import import_config
from ...geodata import geo_info, geo_io, geo_geom
from ...fortlib import pmw_processor


class AMSR2_L1B:
    """
    Read AMSR2_L1B data. Just after obj=AMSR2_L1B(hdfpath), you need to perform obj.set_lat_range(latmin,latmax,freq).
    freq: 6G, 10G, 18G, 36G, 89G_A, 89G_B
    """

    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.params = import_config(config_path=os.path.dirname(
            __file__) + os.sep + "conf/amsr2_params.yaml")

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR2_L1B"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2][2:]
            self.subdsID[name] = i

        self.metadata = self.ds.GetMetadata_Dict()
        self.overlapscans = int(self.metadata["OverlapScans"])
        self.startTime = dt.strptime(
            self.metadata["ObservationStartDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.endTime = dt.strptime(
            self.metadata["ObservationEndDateTime"], "%Y-%m-%dT%H:%M:%S.%fZ")
        self.clip_array = slice(None)

    def keys(self):
        print(self.subdsID)

    def _calc_latlon(self, freq):

        if freq == "89G_A":
            lat = self.get_subds(
                self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
            lon = self.get_subds(
                self.subdsID["Longitude_of_Observation_Point_for_89A"]).ReadAsArray()
        elif freq == "89G_B":
            lat = self.get_subds(
                self.subdsID["Latitude_of_Observation_Point_for_89B"]).ReadAsArray()
            lon = self.get_subds(
                self.subdsID["Longitude_of_Observation_Point_for_89B"]).ReadAsArray()
        elif (freq == "6G") | (freq == "7G") | (freq == "10G") | (freq == "18G") | (freq == "23G") | (freq == "36G"):
            freq_str = freq + "-"
            latA = self.get_subds(
                self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
            lonA = self.get_subds(
                self.subdsID["Longitude_of_Observation_Point_for_89A"]).ReadAsArray()
            coreg_1 = float(
                [
                    var.split(freq_str)[1]
                    for var in self.ds.GetMetadata()["CoRegistrationParameterA1"].split(",")
                    if freq_str in var
                ][0]
            )
            coreg_2 = float(
                [
                    var.split(freq_str)[1]
                    for var in self.ds.GetMetadata()["CoRegistrationParameterA2"].split(",")
                    if freq_str in var
                ][0]
            )
            lat, lon = pmw_processor.calc_latlon_amsr2_l1b(
                latA.T, lonA.T, coreg_1, coreg_2)
            lat = lat.T
            lon = lon.T
        else:
            print("not have the freq data of " + freq)

        return lat, lon

    def get_latlon(self):
        return self.lat, self.lon

    def set_lat_range(self, lamin, lamax, freq, overlap=False):
        """freq: 6G, 7G, 10G, 18G, 23G, 36G, 89G_A, 89G_B"""
        self.freq = freq
        self.clip_array = slice(None)
        lat, lon = self._calc_latlon(freq)
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
        st = self.startTime_split
        et = self.endTime_split
        if st is None:
            print("Run set_lat_range function before get_obs_time.")
        else:
            length = self.length
            width = self.width
            grad = (et - st).seconds / 60 / 60 / 24 / (length - 1)
            date_list = [st + timedelta(days=grad * i) for i in range(length)]
            time_array = np.array(
                [
                    float(dt.strftime("%j")) - 1 + dt.hour / 24 +
                    dt.minute / 24 / 60 + dt.second / 24 / 60 / 60
                    for dt in date_list
                ]
            )
            time_array = np.tile(time_array, (width, 1)).T
        return time_array

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_brightness_temperature(self):
        """
        Just before this method, you need to peform set_lat_range().
        output is a list of h and v brightness temparature for band selected in set_lat_range()
        """
        subdsIDs = {
            "10G": [1, 2],
            "18G": [3, 4],
            "23G": [5, 6],
            "36G": [7, 8],
            "6G": [9, 10],
            "7G": [11, 12],
            "89G_A": [13, 14],
            "89G_B": [15, 16],
        }
        return [self.get_subds(subdsID).ReadAsArray()[self.clip_array] * 0.01 for subdsID in subdsIDs[self.freq]]

    def get_navigation_data(self):
        nav = self.get_subds(self.subdsID["Navigation_Data"]).ReadAsArray()[
            self.clip_array]
        return nav

    def get_satellite_position(self):
        sp0 = self.get_subds(self.subdsID["Navigation_Data"]).ReadAsArray()[
            self.clip_array, 0:3]
        sp0 = np.append(sp0, np.array([(sp0[-1] - sp0[-2]) + sp0[-1]]), axis=0)
        sp = np.array([(sp0[1:] - sp0[0:-1]) * i / self.width + sp0[0:-1] for i in range(self.width)]).transpose(
            2, 1, 0
        )
        return sp

    def get_earth_azimuth(self):
        eaz = self.get_subds(self.subdsID["Earth_Azimuth"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_earth_incidence(self):
        eaz = self.get_subds(self.subdsID["Earth_Incidence"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_sun_elevation(self):
        eaz = self.get_subds(self.subdsID["Sun_Elevation"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_sun_azimuth(self):
        eaz = self.get_subds(self.subdsID["Sun_Azimuth"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_land_ocean_flag(self, freq="36G"):
        num = {"6G": 0, "10G": 1, "18G": 2, "36G": 3}
        loflag = self.get_subds(self.subdsID["Land_Ocean_Flag_6_to_36"]).ReadAsArray()[
            num[freq]][self.clip_array]
        self.output = loflag
        return loflag

    def set_output_prop(self, gcp_x=20, gcp_y=10):
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


class AMSR2_L1R:
    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.params = import_config(config_path=os.path.dirname(
            __file__) + os.sep + "conf/amsr2_params.yaml")

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR2_L1R"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2][2:]
            self.subdsID[name] = i

        self.clip_array = slice(None)

    def keys(self):
        print(self.subdsID)

    def set_lat_range(self, lamin, lamax):
        lat = self.get_subds(
            self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
        lat_1d_min = np.min(lat, axis=1)
        lat_1d_max = np.max(lat, axis=1)
        # print(lat_1d_min, lat_1d_max)
        self.clip_array = (lat_1d_min >= lamin) & (lat_1d_max <= lamax)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_brightness_temperature(self, subdsID):
        return self.get_subds(subdsID).ReadAsArray()[self.clip_array] * 0.01

    def get_latlon(self, resolution="low"):
        # Positions of all resample data in L1R product are ajusted so as to they is consistent with latlons of 89A.
        lat = self.get_subds(
            self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
        lon = self.get_subds(
            self.subdsID["Longitude_of_Observation_Point_for_89A"]).ReadAsArray()
        if resolution == "low":
            lat = lat[:, 0::2]
            lon = lon[:, 0::2]
        lat = lat[self.clip_array]
        lon = lon[self.clip_array]
        return lat, lon

    def get_earth_azimuth(self):
        eaz = self.get_subds(self.subdsID["Earth_Azimuth"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_satellite_position(self):
        sp0 = self.get_subds(self.subdsID["Navigation_Data"]).ReadAsArray()[
            self.clip_array, 0:3]
        sp0 = np.append(sp0, np.array([(sp0[-1] - sp0[-2]) + sp0[-1]]), axis=0)
        sp = np.array([(sp0[1:] - sp0[0:-1]) * i / 243 + sp0[0:-1]
                      for i in range(243)]).transpose(2, 1, 0)
        return sp

    def get_earth_incidence(self):
        eaz = self.get_subds(self.subdsID["Earth_Incidence"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_sun_elevation(self):
        eaz = self.get_subds(self.subdsID["Sun_Elevation"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_sun_azimuth(self):
        eaz = self.get_subds(self.subdsID["Sun_Azimuth"]).ReadAsArray()[
            self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_land_ocean_flag(self, resolution="36GHz"):
        num = {"6GHz": 0, "10GHz": 1, "18GHz": 2, "36GHz": 3}
        loflag = self.get_subds(self.subdsID["Land_Ocean_Flag_6_to_36"]).ReadAsArray()[
            num[resolution]][self.clip_array]
        self.output = loflag
        return loflag

    def set_output_prop(self, gcp_x=20, gcp_y=10):
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


class AMSR2_L2SIC:
    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR2_L1R"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2][2:]
            self.subdsID[name] = i

        tmp_ds = self.get_subds(1)
        self.clip_array = np.full((tmp_ds.RasterYSize), True)  # slice(None)

    def keys(self):
        print(self.subdsID)

    def set_lat_range(self, lamin, lamax):
        lat = self.get_subds(
            self.subdsID["Latitude_of_Observation_Point"]).ReadAsArray()
        lat_1d_min = np.min(lat, axis=1)
        lat_1d_max = np.max(lat, axis=1)
        # print(lat_1d_min, lat_1d_max)
        self.clip_array = (lat_1d_min >= lamin) & (lat_1d_max <= lamax)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_sic(self, l1r_overlap=False):  # have a Bug for overlap. Must Change!!!!
        sic = self.get_subds(self.subdsID["Geophysical_Data"]).ReadAsArray()[
            self.clip_array][:, :, 0] * 0.1
        if l1r_overlap:
            self.clip_array[0:30]
            if self.clip_array[0]:
                sic = np.insert(sic, 0, np.full((30, 243), 101), axis=0)
            if self.clip_array[-1]:
                sic = np.vstack([sic, np.full((30, 243), 101)])
        return sic

    def get_latlon(self, l1r_overlap=False):
        lat = self.get_subds(self.subdsID["Latitude_of_Observation_Point"]).ReadAsArray()[
            self.clip_array]
        lon = self.get_subds(
            self.subdsID["Longitude_of_Observation_Point"]).ReadAsArray()[self.clip_array]
        if l1r_overlap:
            if self.clip_array[0]:
                lat = np.insert(lat, 0, np.full((30, 243), 999), axis=0)
                lon = np.insert(lon, 0, np.full((30, 243), 999), axis=0)
            if self.clip_array[-1]:
                lat = np.vstack([lat, np.full((30, 243), 999)])
                lon = np.vstack([lon, np.full((30, 243), 999)])
        return lat, lon

    def set_output_prop(self, gcp_x=20, gcp_y=10):
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


class AMSR2_L3:
    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR2_L3"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2][2:]
            self.subdsID[name] = i

    def keys(self):
        print(self.subdsID)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_val(self, gain):
        sic = self.get_subds(self.subdsID["Geophysical_Data"]).ReadAsArray()[
            :, :, 0] * gain
        return sic

    def get_latlon(self):
        lat = self.get_subds(
            self.subdsID["Latitude_of_Observation_Point"]).ReadAsArray()
        lon = self.get_subds(
            self.subdsID["Longitude_of_Observation_Point"]).ReadAsArray()
        return lat, lon

    def set_output_prop(self):
        ds = self.get_subds(0)
        dict1 = ds.GetMetadata_Dict()
        # print(dict1["Resolution"])
        # if dict1["Resolution"] == "10km":
        res = int(dict1["Resolution"][0:2])  # km
        print(res)
        if dict1["Projection"] == "PS-N":
            geotrans, geoproj = geo_io.make_north_nsidc_geoinfo(res)
        if dict1["Projection"] == "PS-S":
            geotrans, geoproj = geo_io.make_south_nsidc_geoinfo(res)

        length, width = self.output.shape
        self.output_prop = [width, length, geotrans, geoproj]

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_from_array_and_prop(
            self.output, filepath, self.output_prop, dtype=dtype, no_data=no_data, file_type=file_type
        )
        print("Exported")


class Open:
    def __new__(cls, hdfpath, *args):
        product_name = os.path.basename(hdfpath).split("_")[3][0:8]
        if product_name == "L1SGRTBR":
            load_class = super().__new__(type(AMSR2_L1R(None)))
            load_class.__init__(hdfpath)
        elif product_name == "L1SGBTBR":
            load_class = super().__new__(type(AMSR2_L1B(None)))
            load_class.__init__(hdfpath)
        elif product_name == "L2SGSICL":
            load_class = super().__new__(type(AMSR2_L2SIC(None)))
            load_class.__init__(hdfpath)
        elif product_name[0:5] == "L3SGT":
            load_class = super().__new__(type(AMSR2_L3(None)))
            load_class.__init__(hdfpath)
        else:
            print("the product is not supported")
        return load_class
