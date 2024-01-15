import sys
import os
import numpy as np
from osgeo import osr
import pandas as pd
from osgeo import gdal, ogr
from . import geo_transform, geo_io, geo_geom
from scipy import interpolate
from scipy.spatial import cKDTree as KDTree
from ..fortlib import grid_data


class Search:
    # 本クラスは画素値は入力として与えない。最初に画像Aの各ピクセルの座標値を入力する。
    # その後、指定したSRS・解像度のグリッドセルに対応する画像Aのピクセルインデックスが割り当てられる。
    # search_indexに座標値を入力すると、画像Aのピクセルインデックスを高速で取得することができる。
    def __init__(self, source_ref, target_ref, iarray, jarray, res):
        # source_ref and target_ref: ref=osr.SpatialReference() -> ref.ImportFromEPSG(epsg)
        self.target_ref = target_ref
        self.res = res
        self.build(source_ref, target_ref, iarray, jarray, res)

    def _2Dcoord_to_1Dcoord_and_1Dindex(self, iarray, jarray):
        # for latlon, iarray:lat,jarray:lon
        coord = np.array([iarray.reshape(-1), jarray.reshape(-1)]).T
        length = iarray.shape[0]
        width = jarray.shape[1]
        index_i, index_j = np.meshgrid(
            np.arange(0, length), np.arange(0, width))
        index = np.array([index_i.T.reshape(-1), index_j.T.reshape(-1)])
        return coord, index

    def build(self, source_ref, target_ref, iarray, jarray, res):
        # coord_array : For latlon, array(n,2). (n,0)->lat (n,1)->lon
        # output trans_array: array(n,2). (n,0) -> x axis (horizontal)
        coord_array, index_array = self._2Dcoord_to_1Dcoord_and_1Dindex(
            iarray, jarray)

        coord_transform = osr.CoordinateTransformation(source_ref, target_ref)
        self.trans_array = np.array(
            coord_transform.TransformPoints(coord_array))[:, 0:2]
        self.trans_array = np.int32(self.trans_array / res)

        imin = np.min(self.trans_array[:, 0])
        imax = np.max(self.trans_array[:, 0])
        jmin = np.min(self.trans_array[:, 1])
        jmax = np.max(self.trans_array[:, 1])
        self.nx = imax - imin + 1
        self.ny = jmax - jmin + 1
        self.offset_i = imin
        self.offset_j = jmin
        # ------i -> horizon
        self.trans_array[:, 0] = self.trans_array[:, 0] - self.offset_i
        self.trans_array[:, 1] = self.trans_array[:, 1] - self.offset_j

        self.map_field_i = np.full([self.nx, self.ny], -1)
        self.map_field_j = np.full([self.nx, self.ny], -1)

        # should change array(n,2) to array(2,n)
        trans_array_t = self.trans_array.T

        # self.map_field_i[trans_array_t.tolist()] = index_array[0]
        # self.map_field_j[trans_array_t.tolist()] = index_array[1]
        self.map_field_i[tuple(map(tuple, trans_array_t))] = index_array[0]
        self.map_field_j[tuple(map(tuple, trans_array_t))] = index_array[1]

    def search_index(self, source_ref, iarray, jarray):
        # coord_array : For latlon, array(n,2). (n,0)->lat (n,1)->lon
        # output coord_array: array(n,2). (n,0) -> x axis (horizontal)

        coord_array, index_array = self._2Dcoord_to_1Dcoord_and_1Dindex(
            iarray, jarray)
        index_array = index_array.T

        coord_transform = osr.CoordinateTransformation(
            source_ref, self.target_ref)
        trans_array = np.array(
            coord_transform.TransformPoints(coord_array))[:, 0:2]
        trans_array = np.int32(trans_array / self.res)
        trans_array[:, 0] = trans_array[:, 0] - self.offset_i
        trans_array[:, 1] = trans_array[:, 1] - self.offset_j

        trans_array2 = trans_array[
            (trans_array[:, 0] > 0)
            & (trans_array[:, 0] < self.nx)
            & (trans_array[:, 1] > 0)
            & (trans_array[:, 1] < self.ny)
        ].T
        index_array2 = index_array[
            (trans_array[:, 0] > 0)
            & (trans_array[:, 0] < self.nx)
            & (trans_array[:, 1] > 0)
            & (trans_array[:, 1] < self.ny)
        ].T
        #        out_ij = np.array([self.map_field_i[trans_array2.tolist()], self.map_field_j[trans_array2.tolist()]])
        out_ij = np.array(
            [self.map_field_i[tuple(map(tuple, trans_array2))],
             self.map_field_j[tuple(map(tuple, trans_array2))]]
        )
        bool_array = out_ij[0, :] != -1
        out_ij = out_ij[:, bool_array]
        index_array2 = index_array2[:, bool_array]

        return out_ij, index_array2


class Grid:
    # ある投影座標系・グリッド設定にデータの座標値を変換しピクセル値とともにデータフレームとしてスタックする。
    # 投影グリッドの陸海マスクデータを自動で出力するメソッドを追加。
    # スタックされたデータはpandasのdfとして格納される。
    def __init__(self, source_ref, target_ref, lt_x, lt_y, rb_x, rb_y, res, colname=["data"]):
        # source_ref and target_ref: ref=osr.SpatialReference() -> ref.ImportFromEPSG(epsg)
        # rt_x,rt_y,lb_x,lb_y: corner position in right top and left bottom cells.
        self.target_ref = target_ref
        self.res = res
        self.coord_transform = osr.CoordinateTransformation(
            source_ref, target_ref)
        self.nx = int((rb_x - lt_x) / res)
        self.ny = int((lt_y - rb_y) / res)
        self.rb_x = rb_x
        self.rb_y = rb_y
        self.lt_x = lt_x
        self.lt_y = lt_y
        self.colname_list = ["id", "lat", "lon", "grid_x", "grid_y"]
        self.colname_list.extend(colname)
        self.df = pd.DataFrame(
            index=[], columns=self.colname_list, dtype="float64")
        self.num = 0
        self.id_list = []
        self.concave_hull = []

    def get_geotrans(self):
        geotrans = (self.lt_x, self.res, 0.0, self.lt_y, 0.0, -self.res)
        return geotrans

    def make_empty_raster_ds(self, nodata=0, num_band=1, out_dtype=gdal.GDT_Float32, outfile="/vsimem/output.tif"):
        geotrans = self.get_geotrans()
        pixel_x = self.nx
        pixel_y = self.ny
        geoproj = self.target_ref.ExportToWkt()

        prop = []
        prop.append(pixel_x)
        prop.append(pixel_y)
        prop.append(geotrans)
        prop.append(geoproj)

        tmp_ds = geo_io.make_empty_raster(
            prop, nodata=nodata, num_band=num_band, out_dtype=out_dtype, outfile=outfile)
        return tmp_ds

    def get_land_data(self, nodata=0, outfile="/vsimem/output.tif", all_touched=False, shpfile=None):
        if shpfile is None:
            shpfile = os.path.join(os.path.dirname(
                __file__), "data", "GSHHS_i_L1L5_dissolve.shp")
        print(shpfile)
        tmp_ds = self.make_empty_raster_ds(nodata=nodata, outfile=outfile)
        s_ds = ogr.Open(shpfile, 0)
        poly_ds = geo_transform.reproject_vector(
            s_ds, self.target_ref, outfile="/vsimem/output.shp")
        poly_layer = poly_ds.GetLayer()
        if all_touched:
            gdal.RasterizeLayer(tmp_ds, [1], poly_layer, options=[
                                "ATTRIBUTE=level", "ALL_TOUCHED=TRUE"])
        else:
            gdal.RasterizeLayer(tmp_ds, [1], poly_layer, options=[
                                "ATTRIBUTE=level"])

        del poly_ds
        del s_ds
        return tmp_ds

    def _2Dcoord_to_1Dcoord_and_1Dindex(self, iarray, jarray):
        # for latlon, iarray:lat,jarray:lon
        coord = np.array([iarray.reshape(-1), jarray.reshape(-1)]).T
        length = iarray.shape[0]
        width = jarray.shape[1]
        index_i, index_j = np.meshgrid(
            np.arange(0, length), np.arange(0, width))
        index = np.array([index_i.T.reshape(-1), index_j.T.reshape(-1)]).T

        return coord, index

    # reference: https://soback.jp/detail/26737
    def _list_to_numpy(self, LoL, default=np.nan):
        cols = len(max(LoL, key=len))
        rows = len(LoL)
        AoA = np.empty(
            (
                rows,
                cols,
            )
        )
        AoA.fill(default)
        for idx in range(rows):
            AoA[idx, 0: len(LoL[idx])] = LoL[idx]
        return AoA

    def stack(self, iarray, jarray, varray_list, id, concave_hull=False):
        # coord_array : For latlon, array(n,2). (n,0)->lat (n,1)->lon
        # output coord_array: array(n,2). (n,0) -> x axis (horizontal)
        input_dim = iarray.ndim
        if input_dim == 2:
            coord_array, index_array = self._2Dcoord_to_1Dcoord_and_1Dindex(
                iarray, jarray)
            trans_array = np.array(
                self.coord_transform.TransformPoints(coord_array))[:, 0:2]
        else:
            coord_array = np.array([iarray, jarray]).T
            index_array = np.array(
                [np.arange(0, iarray.shape[0]), np.arange(0, iarray.shape[0])]).T
            trans_array = np.array(
                self.coord_transform.TransformPoints(coord_array))[:, 0:2]

        if concave_hull:
            coord_array_ch = self._get_concave_hull(iarray, jarray)
            trans_array_ch = np.array(
                self.coord_transform.TransformPoints(coord_array_ch))[:, 0:2]
            self.concave_hull.append(self._get_concave_hull_mask(
                trans_array_ch.astype(np.float64)))

        mask = (trans_array[:, 0] > self.lt_x) & (trans_array[:, 0] < self.rb_x) & (
            trans_array[:, 1] > self.rb_y) & (trans_array[:, 1] < self.lt_y)

        trans_array2 = trans_array[mask].T
        index_array2 = index_array[mask].T
        coord_array2 = coord_array[mask].T

        # trans_array3 = np.empty(trans_array2.shape)
        trans_array2[0, :] = (trans_array2[0, :] - self.lt_x) / self.res
        trans_array2[1, :] = (trans_array2[1, :] - self.rb_y) / self.res

        # trans_array3 = np.rint(trans_array2)
        num_array = np.full(coord_array2.shape[1], id).reshape(
            1, coord_array2.shape[1])
        varray_list2 = [num_array, coord_array2, trans_array2]

        if input_dim == 2:
            for varray in varray_list:
                varray = varray[tuple(map(tuple, index_array2))]
                varray = varray.reshape(1, varray.shape[0])
                varray_list2.append(varray)
        else:
            for varray in varray_list:
                varray = varray[mask]
                varray = varray.reshape(1, varray.shape[0])
                varray_list2.append(varray)

        test = np.concatenate(varray_list2, axis=0).T

        # gi,gj:grid coord, si,sj: data coord
        df = pd.DataFrame(
            data=test, columns=self.colname_list, dtype="float64")
        print("stack num=", len(df))
        self.df = pd.concat([self.df, df], ignore_index=True)

        id_list = self.id_list
        id_list.append(id)
        self.id_list = list(set(id_list))
        self.num = len(self.id_list)

    def _get_concave_hull_mask(self, coord):
        ds_ras = self.make_empty_raster_ds(nodata=0, num_band=1)
        geom = geo_geom.create_polygon(coord[:, ::-1])
        ds_poly = geo_io.make_vector_from_geom(
            geom, epsg=None, srs=self.target_ref, outfile="/vsimem/output.shp")
        ds_dst = geo_transform.rasterize_by_raster_with_proj(
            ds_ras, ds_poly, outfile="/vsimem/output2.tif", attr=False)
        output = ds_dst.ReadAsArray()
        output = output == 1
        ds_dst = None
        ds_poly = None
        ds_ras = None
        return output

    def _get_concave_hull(self, iarray, jarray):
        length, width = iarray.shape
        iarray2 = np.concatenate(
            [iarray[0, :], iarray[:, width - 1], iarray[length - 1, ::-1], iarray[::-1, 0]])
        jarray2 = np.concatenate(
            [jarray[0, :], jarray[:, width - 1], jarray[length - 1, ::-1], jarray[::-1, 0]])
        coord = np.array([iarray2, jarray2]).T
        return coord

    def get_image(self, var_name, id=0, interp="nearest", mode=1, param=1):
        """
        mode 1 -> kdtree mask, param = distance (float, unit: cell num)
        mode 2 -> mask by concave hull numpy array. param is note used in this mode.
        mode 3 -> fast computation for nearest neighbor using fortran library. interp is not used in this mode.
        """
        nx = self.nx
        ny = self.ny
        X = np.linspace(1, nx, nx)
        Y = np.linspace(1, ny, ny)
        X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
        df = self.df[self.df["id"] == id]
        grid_x = df["grid_x"].values.astype(np.float32) + 0.5
        grid_y = df["grid_y"].values.astype(np.float32) + 0.5
        grid_z = df[var_name].values.astype(np.float32)

        if mode == 1:
            tree = KDTree(np.c_[grid_x, grid_y])
            dist, _ = tree.query(np.c_[X.ravel(), Y.ravel()], k=1)
            dist = dist.reshape(X.shape)
            mask = dist < param
            Z = interpolate.griddata(
                list(zip(grid_x, grid_y)), grid_z, (X, Y), method=interp)[::-1, :]

        elif mode == 2:
            mask = self.concave_hull[id]
            if mask.sum() == 0:
                return np.full(mask.shape, np.nan)
            Z = interpolate.griddata(
                list(zip(grid_x, grid_y)), grid_z, (X, Y), method=interp)[::-1, :]

        elif mode == 3:
            mask = self.concave_hull[id]
            mask_grid = np.zeros([self.nx, self.ny])
            wsize, res, threshold = param[:]
            Z = grid_data.nearest_neighbor(
                grid_x, grid_y, grid_z, mask_grid, wsize, res, threshold).T[::-1, :]

        Z[~mask] = np.nan

        return Z

    def remove(self, varname, value):
        df = self.df
        self.df = df[df[varname] != value]

    def flush(self):
        self.df = pd.DataFrame(
            index=[], columns=self.colname_list, dtype="float64")
        self.num = 0
        self.concave_hull = []
        self.id_list = []

    def export_numpy(self):
        groups = self.df.groupby(["grid_i", "grid_j"])
        result = []
        for colname in self.colname_list[2:]:
            groups = self.df.groupby(["grid_i", "grid_j"])
            glist = groups["grid_x"].apply(list)
            result.append(self._list_to_numpy(glist.values))
        glist = glist.reset_index()
        grid_i = np.array(glist["grid_i"].values)
        grid_j = np.array(glist["grid_j"].values)
        return grid_i, grid_j, np.array(result)

    def export_output(self, filepath, no_data=9.9e33, file_type="GTiff", dtype=gdal.GDT_Float32):
        geotrans = self.get_geotrans()
        pixel_x = self.nx
        pixel_y = self.ny
        geoproj = self.target_ref.ExportToWkt()
        output_prop = [pixel_x, pixel_y, geotrans, geoproj]
        geo_io.make_raster_from_array_and_prop(
            self.output, filepath, output_prop, dtype, no_data, file_type)
        print("Exported")
