from ...helpers import settings
from ...helpers.misc import import_config
import pyrttov
import sys
import os
import re
import glob
import matplotlib.pyplot as plt
import numpy as np


# from rttov_wrapper_f2py import rttov_load_inst, rttov_call_direct, rttov_drop_all


class _BaseConf:
    def __init__(self):
        conf_list = glob.glob(os.path.join(
            os.path.dirname(__file__), "conf", "*.yaml"))
        # conf_list = glob.glob(os.path.join("./conf", "*.yaml"))
        num = [i+1 for i in range(len(conf_list))]
        self.conf_files = dict(zip(num, conf_list))

    def keys(self):
        print(self.conf_files)

    def read_conf(self, conf_num):
        conf = import_config(config_path=self.conf_files[conf_num])
        return conf


BaseConf = _BaseConf()


class Profiles_tool:
    def __init__(self, nprofiles, nlevels, nsurfaces, conf, gas_units=1, mmr_hydro=1, mmr_aer=1):
        self.nprofiles = nprofiles
        self.nlevels = nlevels
        self.nsurfaces = nsurfaces
        self.nlayers = nlevels-1
        self.conf = conf
        # "unknown" Default initialisation, ppmv over moist air will be used
        # "ppmv_dry" Gas units of ppmv over dry air
        # "kg_per_kg" Gas units of kg/kg over moist air
        # "ppmv_wet" Gas units of ppmv over moist air
        # gas_units=pyrttov.gasUnitType('kg_per_kg')
        self.profiles = pyrttov.Profiles(
            self.nprofiles, self.nlevels, self.nsurfaces)
        self.profiles.GasUnits = gas_units  # 1 => specific humidity in kg/kg for moisture

    def set_profiles_from_obj(self, obj):
        for key in vars(obj).keys():
            exec("self.profiles."+key+"=obj."+key)

    def confirm_profiles(self):
        conf = self.conf
        j1 = self._confirm(conf["Profiles"]["General"])
        j2 = self._confirm(conf["Profiles"]["Atmosphere"])
        j3 = self._confirm(conf["Profiles"]["Surface"])
        if j1 & j2 & j3:
            print("OK")

    def _confirm(self, params):
        n2 = 0
        for key in params.keys():
            if params[key]:
                exec("self.tmp=self.profiles."+key)
                if self.tmp is None:
                    print(key + " is not set up")
                    n2 = n2+1
        if n2 == 0:
            return True
        else:
            return False

    def set_PHalf_geomspace(self, p_surface, p_top):
        # geometric spacing in log-pressure
        half_pressure_levels = np.geomspace(p_surface, p_top, self.nlevels)
        self.profiles.PHalf = self.expand_nprofiles(
            p_surface-half_pressure_levels)

    def set_general_profiles(self, date=[2015, 1, 1, 0, 0, 0], geom=[80, 0.0, 0.0], angle=[55.0, 0.0, 0.0, 0.0]):
        self.profiles.DateTimes = np.tile(np.array(date), (self.nprofiles, 1))
        self.profiles.SurfGeom = np.tile(np.array(geom), (self.nprofiles, 1))
        self.profiles.Angles = np.tile(np.array([angle]), (self.nprofiles, 1))

    def expand_nprofiles(self, data):
        if len(data.shape) == 1:
            output = np.tile(data, (self.nprofiles, 1))
        if len(data.shape) == 2:
            output = np.tile(data, (self.nprofiles, 1, 1))
        return output

    def expand_nsurfaces(self, data):
        if len(data.shape) == 1:
            output = np.tile(data, (self.nsurfaces, 1))
        if len(data.shape) == 2:
            n1, n2 = data.shape
            output = np.tile(data.reshape((n1, 1, n2)), (1, self.nsurfaces, 1))
        return output


class RTTOV_tools:
    def __init__(self, nprofiles, nlevels, nsurfaces, conf, coefs, sens="mw"):
        # sens: mw, ir, vis
        #    self.rttov_dir=os.path.dirname(__file__)
        self.rttov_dir = settings.PYRTTOV_PATH + os.sep + "../"
        self.sens = sens
        self.nprofiles = nprofiles
        self.nlevels = nlevels
        self.nsurfaces = nsurfaces
        self.nlayers = nlevels-1
        self.conf = conf

        self.rttov = pyrttov.Rttov()
        self._set_options(conf["Options"]["Settings"])
        self.rttov.FileCoef = coefs["FileCoef"]
        if "FileHydrotable" in coefs:
            self.rttov.FileHydrotable = coefs["FileHydrotable"]

        try:
            self.rttov.loadInst()
        except pyrttov.RttovError as e:
            sys.stderr.write("Error loading instrument: {!s}\n".format(e))
            sys.exit(1)

        self.nchannels = self.rttov.Nchannels  # should be 14 for AMSR2

        self.atlas = pyrttov.Atlas()
        if self.sens == "mw":
            self.atlas.AtlasPath = '{}/{}'.format(self.rttov_dir, "emis_data")
        elif self.sens == "vis":
            self.atlas.AtlasPath = '{}/{}'.format(self.rttov_dir, "brdf_data")
        elif self.sens == "ir":
            self.atlas.AtlasPath = '{}/{}'.format(
                self.rttov_dir, "iremis_data")

    def view_profiles_size(self):
        conf = self.conf
        self._view_profiles_size(conf["Profiles"]["General"])
        self._view_profiles_size(conf["Profiles"]["Atmosphere"])
        self._view_profiles_size(conf["Profiles"]["Surface"])

    def _view_profiles_size(self, params):
        tmp_profiles = pyrttov.Profiles

        p = r'\[.*\]'
        print("Variable: Nsize")
        for key in params.keys():
            if params[key]:
                # exec do not work for local var. (self is needed)
                exec("self.tmp=tmp_profiles."+key+".__doc__")
                nstr = re.findall(p, self.tmp)[0]
                print(key, ":", nstr)

    def _set_options(self, params):
        for key in params.keys():
            exec("self.rttov.Options."+key+"="+str(params[key]))

    def load_profiles(self, profiles):
        try:
            self.rttov.Profiles = profiles
        except pyrttov.RttovError as e:
            sys.stderr.write("Error setting profiles: {!s}\n".format(e))
            sys.exit(1)

    def load_emissivity_atlas(self, IncSea=False, IncLand=True, IncSeaIce=True, MaxDistance=0):
        surfemis = np.zeros(
            (5, self.nprofiles, self.nsurfaces, self.nchannels), dtype=np.float64)
        surfemis[:, :, :, :] = -1.
        self._initialize_atlas()
        surfemis[0, :, :, :] = self._get_emis_ref_from_atlas(
            IncSea=IncSea, IncLand=IncLand, IncSeaIce=IncSeaIce)
        self.rttov.SurfEmisRefl = surfemis

    def _initialize_atlas(self, camel_version=3, year=2007, ang_corr=False, atlas_id=-1):
        month = self.rttov.Profiles.DateTimes[0][1]
        if self.sens == "mw":
            self.atlas.loadMwEmisAtlas(month, inst=None, year=0, atlas_id=-1)
        elif self.sens == "vis":
            self.atlas.loadBrdfAtlas(month, inst=None, atlas_id=-1)
        elif self.sens == "ir":
            self.atlas.loadIrEmisAtlas(
                month, inst=None, camel_version=3, year=2007, ang_corr=False, atlas_id=-1)

    def _get_emis_ref_from_atlas(self, IncSea=True, IncLand=True, IncSeaIce=True, MaxDistance=0):
        self.atlas.IncSea = IncSea
        self.atlas.IncSeaIce = IncSeaIce
        self.atlas.IncLand = IncLand
        self.atlas.MaxDistance = MaxDistance

        try:
            emisref = self.atlas.getEmisBrdf(self.rttov)
        except pyrttov.RttovError as e:
            sys.stderr.write("Error calling atlas: {!s}\n".format(e))

        return emisref

    def run(self, mode="forward"):
        if mode == "forward":
            self.rttov.runDirect()
        elif mode == "jacobian":
            self.rttov.runK()
        elif mode == "tangent_linear":
            self.rttov.runTL()
        elif mode == "ajoint":
            self.rttov.runAD()

    # --------------runDirect Output------------------------------
    def get_surf_emissivity(self):
        return self.rttov.SurfEmis

    def get_radiance(self):
        return self.rttov.Rads

    def get_cloudy_BT(self):
        return self.rttov.Bt

    def get_clearSky_BT(self):
        return self.rttov.BtClear

    def get_delta_BT_cloud(self):
        return self.rttov.Bt - self.rttov.BtClear

    def get_tau_at_all_levels(self):
        return self.rttov.TauLevels

    def get_tau_at_surface_level(self):
        return self.rttov.EmisTermsTauSfc

    def get_upward_radiance(self):
        return self.rttov.EmisTermsRadUp

    def get_downward_radiance(self):
        return self.rttov.EmisTermsRadDown

    def get_yacobian_TQ(self):
        t = self.rttov.TK
        q = self.rttov.QK
        return t, q

    def get_yacobian_surfem(self):
        return self.rttov.SkinK

    def get_yacobian_skin(self):
        return self.rttov.SkinK

    def get_yacobian_nsurface(self):
        return self.rttov.NearSurfaceK

    def get_yacobian_scatt_params(self):
        clw = self.rttov.MWClwK
        ciw = self.rttov.MWCiwK
        snow = self.rttov.MWSnowK
        rain = self.rttov.MWRainK
        return clw, ciw, snow, rain

    def get_yacobian_sfraction(self):
        return self.rttov.SurfaceFractionK
