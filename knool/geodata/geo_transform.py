from osgeo import gdal, osr, ogr
import os
import numpy as np
from . import ogr2ogr
from . import gdal_merge as merge
from . import geo_info, geo_io

OGRTypes = {int: ogr.OFTInteger, str: ogr.OFTString, float: ogr.OFTReal}


def reverse_clip(in_maskFile, out_maskFile):
    #    options = ["", "-f", "ESRI Shapefile", "-sql",
    #    ''"SELECT ST_Difference(a.geometry, b.geometry) AS geometry FROM a,'b.shp'.b"'',
    #    wkt_pr, outmaskFile, maskFile]
    #    ogr2ogr difference.shp a.shp -dialect SQLite -sql
    geoproj = osr.SpatialReference()
    geoproj.ImportFromEPSG(3411)
    wkt_pr = geoproj.ExportToWkt()
    options = ["", "-f", "ESRI Shapefile", "-t_srs", wkt_pr, out_maskFile, in_maskFile]
    ogr2ogr.main(options)


def reproject(ds, outfile="/vsimem/output.tif", epsg_str="EPSG:4326", NODATA_VALUE=0):
    output_ds = gdal.Warp(outfile, ds, dstSRS=epsg_str, resampleAlg="bilinear", dstNodata=NODATA_VALUE)
    return output_ds


def geocode(
    ds,
    outfile="/vsimem/output.tif",
    src_epsg_str="EPSG:4326",
    dst_epsg_str="EPSG:4326",
    resample_alg="near",
    trans_bug_correction=False,
    NODATA_VALUE=None,
):
    # ds -> raster with GCPs
    # どの座標系のGCPでもOK。また、GCPに定義されている座標系の変換（再投影）にも対応できるようコーディング。
    # うまく再投影できない場合は、xとyが逆の可能性あり（gdalのバグ）。その場合はtrans_bug_correctionをTrueにする。
    gcps = ds.GetGCPs()
    src_epsg = int(src_epsg_str[5:])
    dst_epsg = int(dst_epsg_str[5:])
    target_ref = osr.SpatialReference()
    target_ref.ImportFromEPSG(dst_epsg)
    trans = geo_info.get_coord_transform_epsg(src_epsg, dst_epsg)

    gcps2 = []
    for gcp in gcps:
        if not trans_bug_correction:
            x, y = trans.TransformPoint(gcp.GCPX, gcp.GCPY)[0:2]
        else:
            x, y = trans.TransformPoint(gcp.GCPY, gcp.GCPX)[0:2]
        gcp2 = gdal.GCP(x, y, 0, gcp.GCPPixel, gcp.GCPLine)
        gcps2.append(gcp2)

    tmp_ds = geo_io.copy_dataset(ds)
    tmp_ds.SetGCPs(gcps2, target_ref)

    output_ds = gdal.Warp(
        outfile,
        tmp_ds,
        srcSRS=None,  # src_epsg_str,
        dstSRS=dst_epsg_str,
        resampleAlg=resample_alg,
        dstNodata=NODATA_VALUE,
        tps=True,
        #    outputBounds=bound,
    )
    return output_ds


def adjust_shape(ds, ds2, outfile="/vsimem/output.tif"):
    cols = ds2.RasterXsize
    rows = ds2.RasterYSize
    output_ds = gdal.Warp(outfile, ds, width=cols, height=rows)
    return output_ds


def clip_rectangle(ds, minX, minY, maxX, maxY, outfile="/vsimem/output.tif"):
    output_ds = gdal.Warp(outfile, ds, outputBounds=(minX, minY, maxX, maxY))
    return output_ds


def warp(ds, outfile, epsg_str, xres, yres, minX, minY, maxX, maxY, resample):
    output_ds = gdal.Warp(
        outfile, ds, dstSRS=epsg_str, xRes=xres, yRes=yres, outputBounds=(minX, minY, maxX, maxY), resampleAlg=resample
    )
    return output_ds


def edit_tifftag(ds, option_str=None, tag_dict=None, outfile="/vsimem/output.tif"):
    # For th edetails, see gdal official HP https://gdal.org/drivers/raster/gtiff.html
    # TFW=YES :
    # RPB=YES : RPC 情報が利用可能な場合
    # TILED=YES : デフォルトでは、ストライプ化された TIFF ファイルが作成されます。このオプションは、タイル化された TIFF ファイルの作成を強制するために使用できます。
    # BLOCKXSIZE=n : タイル幅を設定します。デフォルトは 256 です。
    # BLOCKYSIZE=n : タイルまたはストリップの高さを設定します。タイルの高さのデフォルトは 256 で、ストリップの高さのデフォルトは、1 つのストリップが 8K 以下になるような値です。
    # COMPRESS=[JPEG/LZW/PACKBITS/DEFLATE/CCITTRLE/CCITTFAX3/CCITTFAX4/LZMA/ZSTD/LERC/LERC_DEFLATE/LERC_ZSTD/WEBP/JXL/NONE]
    # For example, option_srt = "-of GTiff -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=240 -co BLOCKYSIZE=128"
    # For example, tag_dict={"TIFFTAG_YRESOLUTION":"1", "TIFFTAG_RESOLUTIONUNIT":"None"}
    if tag_dict is not None:
        for name in tag_dict.keys():
            ds.SetMetadataItem(name, tag_dict[name])

    output_ds = gdal.Translate(outfile, ds, options=option_str)

    return output_ds


def collocate(outfile, infile1, infile2, ulx, uly, lrx, lry):
    # fmt: off
    options = ["", "-o", outfile, "-separate", "-n", "0", "-ul_lr",
               str(ulx), str(uly), str(lrx), str(lry), infile1, infile2]
    # fmt: on
    print(options)
    merge.main(options)


def mosaic(infile_list, outfile):
    options = ["", "-o", outfile]
    for infile in infile_list:
        options.append(infile)
    merge.main(options)


def raster_to_polygon(raster, dst_layername):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_layername + ".shp")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    gdal.Polygonize(srcband, srcband, dst_layer, -1, ["8CONNECTED=8"])


def get_vector_from_raster(raster):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("GPKG")
    dst_ds = drv.CreateDataSource(r"/vsimem/output.gpkg")
    dst_layer = dst_ds.CreateLayer("output", srs=srs)
    gdal.Polygonize(srcband, srcband, dst_layer, -1, ["8CONNECTED=8"])
    vector = ogr.Open("/vsimem/output.gpkg")
    return vector


def raster_to_polygon_with_field(raster, dst_layername, field_name):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_layername + ".shp")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    fd = ogr.FieldDefn(field_name, ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex(field_name)
    gdal.Polygonize(srcband, srcband, dst_layer, dst_field, ["8CONNECTED=8"])


def rasterize_by_raster_with_gcps(
    raster_ds,
    poly_ds,
    outfile="/vsimem/output.tif",
    src_epsg_str="EPSG:4326",
    dst_epsg_str="EPSG:4326",
    res=0.1,
    trans_bug_correction=False,
    out_dtype=gdal.GDT_Float32,
):
    # 注意：このメソッドはラスタライズの際にGCP画像上に投影するわけではない。（GDALでは簡単なコードでそれを実装できない。）
    # GCP画像の範囲を読み込み、指定した分解能と座標系でラスタライズする。
    prop = geo_info.get_property_from_raster_with_gcps(raster_ds)

    poly_layer = poly_ds.GetLayer()

    num_band = 1
    gcps = prop[3]

    src_epsg = int(src_epsg_str[5:])
    dst_epsg = int(dst_epsg_str[5:])
    target_ref = osr.SpatialReference()
    target_ref.ImportFromEPSG(dst_epsg)
    geoproj_wkt = target_ref.ExportToWkt()
    trans = geo_info.get_coord_transform_epsg(src_epsg, dst_epsg)

    gcps_pixel = []
    gcps_line = []
    for gcp in gcps:
        if not trans_bug_correction:
            x, y = trans.TransformPoint(gcp.GCPX, gcp.GCPY)[0:2]
        else:
            x, y = trans.TransformPoint(gcp.GCPY, gcp.GCPX)[0:2]
        gcps_pixel.append(x)
        gcps_line.append(y)
    #        print(x,y)

    array1 = np.array(gcps_pixel)
    array2 = np.array(gcps_line)
    pmin = np.min(array1)
    pmax = np.max(array1)
    lmin = np.min(array2)
    lmax = np.max(array2)

    pixel_x = int((pmax - pmin) / res)
    pixel_y = int((lmax - lmin) / res)

    geotrans = (pmin - res / 2, res, 0.0, lmax + res / 2, 0.0, -res)

    driver = gdal.GetDriverByName("GTiff")
    tmp_ds = driver.Create(outfile, pixel_x, pixel_y, num_band, out_dtype)
    tmp_ds.SetGeoTransform(geotrans)
    tmp_ds.SetProjection(geoproj_wkt)

    band = tmp_ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.FlushCache()

    gdal.RasterizeLayer(tmp_ds, [1], poly_layer, options=["ATTRIBUTE=id", "ALL_TOUCHED=TRUE"])  # , burn_values=[1]
    tmp_ds.FlushCache()


def rasterize_by_raster_with_proj(raster_ds, poly_ds, outfile="/vsimem/output.tif"):
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.CreateCopy(outfile, raster_ds)
    poly_layer = poly_ds.GetLayer()
    gdal.RasterizeLayer(dst_ds, [1], poly_layer, options=["ATTRIBUTE=id", "ALL_TOUCHED=TRUE"])  # , burn_values=[1]
    return dst_ds


def reproject_vector(s_ds, t_srs, outfile="/vsimem/output.shp"):
    # s_ds = ogr.Open(source_shp, 0)
    s_lyr = s_ds.GetLayer(0)
    s_srs = s_lyr.GetSpatialRef()

    coord_trans = osr.CoordinateTransformation(s_srs, t_srs)

    t_driver = ogr.GetDriverByName("ESRI Shapefile")
    t_ds = t_driver.CreateDataSource(outfile)
    t_lyr_name = os.path.splitext(os.path.split(outfile)[1])[0]
    t_lyr = t_ds.CreateLayer(t_lyr_name, geom_type=ogr.wkbMultiPolygon)

    s_lyr_defn = s_lyr.GetLayerDefn()
    for i in range(0, s_lyr_defn.GetFieldCount()):
        field_defn = s_lyr_defn.GetFieldDefn(i)
        t_lyr.CreateField(field_defn)

    t_lyr_defn = t_lyr.GetLayerDefn()

    for feature in s_lyr:
        geom = feature.GetGeometryRef()
        geom.Transform(coord_trans)
        t_feature = ogr.Feature(t_lyr_defn)
        t_feature.SetGeometry(geom)

        for i in range(0, t_lyr_defn.GetFieldCount()):
            t_feature.SetField(t_lyr_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))

        t_lyr.CreateFeature(t_feature)

        t_feature.Destroy()
        feature.Destroy()

    # t_srs.MorphToESRI()
    # f_name = os.path.splitext(outfile)[0]
    # file = open(f_name + ".prj", "w")
    # file.write(t_srs.ExportToWkt())
    # file.close()

    s_ds.Destroy()
    # t_ds.Destroy()
    return t_ds
