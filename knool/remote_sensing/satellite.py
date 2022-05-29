import sys
import os
import numpy as np
import pandas as pd
from skyfield.api import load
from datetime import datetime as dt
from datetime import timedelta
from skyfield.api import utc
from osgeo import gdal, osr, ogr
from ..geodata_processor import geo_io, geo_info,geo_geom

class tles_for_asat():
    def __init__(self,infile):
        self.tles = load.tle_file(infile)
        date_list=[]
        num_list=[]
        for i, tle in enumerate(self.tles):
            date_list.append(tle.epoch.utc_datetime())
            num_list.append(i)
        self.df0=pd.DataFrame(num_list,columns=['num']).assign(Date=pd.to_datetime(date_list)).set_index('Date')
        self.df=self.df0
        
    def set_period(self,fdate,ldate):
        ftime = dt.strptime(fdate, '%Y-%m-%d %H:%M:%S').replace(tzinfo=utc)
        ltime = dt.strptime(ldate, '%Y-%m-%d %H:%M:%S').replace(tzinfo=utc)
        self.df=self.df0[(self.df0.index>ftime) & (self.df0.index<ltime)]

    def _get_position(self,num,date_sf):
        geocentric=self.tles[num].at(date_sf)
        subpoint = geocentric.subpoint()
        lat = subpoint.latitude.degrees
        lon = subpoint.longitude.degrees
        ele = subpoint.elevation.m
        return lat,lon,ele
    
    def calc_position_at(self,date):
        #input:
        # date: datetime object -> e.g., dt(2018, 2, 1, 12, 15, 30, 2000, tzinfo=utc)
        df=self.df
        row1=df.loc[df.index==df.index.unique()[df.index.unique().get_loc(pd.to_datetime(date), method='nearest')]]
        date1_sf=load.timescale().from_datetime(date)
        num=row1.num.values[0]
        return self._get_position(num,date1_sf)
    
    def calc_positions_between(self,fdate,ldate,interval):
        #input:
        # ldate, fdate: text -> e.g., 2011-07-01
        # interval: int (min) -> 30 (min)
        ds_date=pd.date_range(start=fdate, end=ldate, freq=str(interval)+"min",tz=utc)
        pos_list=[]
        date_list=[]
        for date in ds_date:
            lat, lon, ele = self.calc_position_at(date)
            pos_list.append([lat,lon])
            date_list.append(date)
        self.output=pos_list
        self.output_attribute= date_list
        self.output_type="csv_list_line"
        return pos_list

    def calc_positions_faster_between(self,fdate,ldate,interval):
        fdt = dt.strptime(fdate, '%Y-%m-%d %H:%M:%S')
        ldt = dt.strptime(ldate, '%Y-%m-%d %H:%M:%S')
        mdt=fdt+(ldt-fdt)/2
        mdt=mdt.replace(tzinfo=utc)
        
        df=self.df
        row1=df.loc[df.index==df.index.unique()[df.index.unique().get_loc(pd.to_datetime(mdt), method='nearest')]]
        num=row1.num.values[0]
 
        date_ds=pd.date_range(start=fdate, end=ldate, freq=str(interval)+"min",tz=utc)
        date_list=date_ds.to_list()
        date_sf=load.timescale().from_datetimes(date_list)
        
        lat, lon, _ =self._get_position(num,date_sf)
        self.output=np.stack([lat,lon]).T
        self.output_attribute= date_list
        self.output_type="csv_list_line"
        return self.output

    def calc_buff_positions_between(self,fdate,ldate,interval,distance,ori="right"):
        #input:
        # ldate, fdate: text -> e.g., 2011-07-01
        # interval: int (min) -> 30 (min)
        fdt = dt.strptime(fdate, '%Y-%m-%d %H:%M:%S')
        ldt = dt.strptime(ldate, '%Y-%m-%d %H:%M:%S')
        ldt = ldt + timedelta(minutes=interval)
        ds_date=pd.date_range(start=fdt, end=ldt, freq=str(interval)+"min",tz=utc)
        pos_list=[]
        date_list=[]
        for i, date in enumerate(ds_date):
            lat, lon, ele = self.calc_position_at(date)
            
            date2=ds_date[i+1]
            lat2, lon2, ele2 = self.calc_position_at(date2)

            lat3,lon3,ele3=geo_info.calc_line_buffer_point(lat,lon,0,lat2,lon2,0,distance,ori)
            pos_list.append([lat3,lon3])
            date_list.append(date)
        self.output_type="csv_list_line"
        self.output_attribute = date_list
        self.output=pos_list
        return pos_list

    def calc_buff_positions_faster_between(self,fdate,ldate,interval,distance,ori="right"):
        fdt = dt.strptime(fdate, '%Y-%m-%d %H:%M:%S')
        ldt = dt.strptime(ldate, '%Y-%m-%d %H:%M:%S')
        mdt=fdt+(ldt-fdt)/2
        mdt=mdt.replace(tzinfo=utc)
        
        df=self.df
        row1=df.loc[df.index==df.index.unique()[df.index.unique().get_loc(pd.to_datetime(mdt), method='nearest')]]
        num=row1.num.values[0]
        
        date_ds=pd.date_range(start=fdt, end=ldt+timedelta(minutes=interval), freq=str(interval)+"min",tz=utc)
        date_list=date_ds.to_list()
        date_sf=load.timescale().from_datetimes(date_list)
        
        lat, lon, _ =self._get_position(num, date_sf)
        
        lat2,lon2,_=geo_info.calc_line_buffer_point(lat[0:-1],lon[0:-1],0,lat[1:],lon[1:],0,distance,ori)
        self.output=np.stack([lat2,lon2]).T
        self.output_attribute= date_list[0:-1]
        self.output_type="csv_list_line"
        return self.output

    def _calc_buff_positions__left_right_between(self,fdate,ldate,interval,distance,ori="middle"):

        if ori=="middle":
            rpos_array=self.calc_buff_positions_faster_between(fdate,ldate,interval,distance/2,ori="right")
            lpos_array=self.calc_buff_positions_faster_between(fdate,ldate,interval,distance/2,ori="left")
        elif ori=="right":
            rpos_array=self.calc_buff_positions_faster_between(fdate,ldate,interval,distance,ori="right")
            lpos_array=self.calc_positions_faster_between(fdate,ldate,interval)
        elif ori=="left":
            lpos_array=self.calc_buff_positions_faster_between(fdate,ldate,interval,distance,ori="left")
            rpos_array=self.calc_positions_faster_between(fdate,ldate,interval)
        else:
            print("The mode is not currently supported.")
        
        return lpos_array,rpos_array  

    def calc_buff_area_between(self,fdate,ldate,interval,distance,ori="middle"):
        lpos_array,rpos_array=self._calc_buff_positions__left_right_between(fdate,ldate,interval,distance,ori)
        
        result=np.concatenate([rpos_array,lpos_array[::-1]])
        self.output=result
        self.output_attribute=None
        self.output_type="csv_list_polygon"
        return self.output   


    def calc_scene_areas(self,fdate,ldate,interval,distance,pinterval,ori="middle"):
        lpos_array,rpos_array=self._calc_buff_positions__left_right_between(fdate,ldate,interval,distance,ori)
        
        row_num=rpos_array.shape[0]
        col_num=2
        bit=8
        pol_num=int((row_num-1)/(pinterval-1))

        lpos_sub = np.lib.index_tricks.as_strided(lpos_array.copy(),(row_num,pinterval,1,col_num),(bit*col_num*(pinterval-1),bit*col_num,8,8))[0:pol_num]
        rpos_sub = np.lib.index_tricks.as_strided(rpos_array.copy(),(row_num,pinterval,1,col_num),(bit*col_num*(pinterval-1),bit*col_num,8,8))[0:pol_num]
        rpos_sub2=rpos_sub[:,::-1,:,:]
        outlines_sub=np.concatenate([lpos_sub,rpos_sub2],axis=1)[:,:,0,:]

        self.output=outlines_sub
        date_list=self.output_attribute
        stime_list=[date_obj.strftime('%Y-%m-%d %H:%M:%S') for date_obj in date_list[:-(pinterval-1):(pinterval-1)]]
        ftime_list=[date_obj.strftime('%Y-%m-%d %H:%M:%S') for date_obj in date_list[(pinterval-1)::(pinterval-1)]]
        print(len(stime_list),len(ftime_list),outlines_sub.shape)
        self.output_attribute={"StartTime":stime_list,"EndTime":ftime_list}
        self.output_type="csv_list_polygons"
        return self.output

    def calc_intersect_scene_areas(self,fdate,ldate,interval,distance,pinterval,vecfile,ori="middle"):

        self.calc_scene_areas(fdate,ldate,interval,distance,pinterval,ori)
        geom2_list=geo_geom.create_polygons(self.output)
        attr2_dict=self.output_attribute
        keys=attr2_dict.keys()

        vector=ogr.Open(vecfile)
        layer1=vector.GetLayer()
        geom1_list=[]
        for feat in layer1:
            geom1_list.append(feat.GetGeometryRef().Clone())

        result_geom_list=[]
        result_geom2_list=[]
        result_attr_list=[]
        for geom1 in geom1_list:
            poly2_intersect_list,bool_list=geo_geom.find_intersect_polygons(geom1,geom2_list)
            geom2_sub_list=[i for (i,v) in zip(geom2_list,bool_list) if v]
            for key in keys:
                attr2_dict[key]= [i for (i,v) in zip(attr2_dict[key],bool_list) if v]
            result_geom_list.append(poly2_intersect_list)
            result_geom2_list.append(geom2_sub_list)
            result_attr_list.append(attr2_dict)

        self.output=result_geom_list
        self.output_attribute=result_attr_list
        self.output_type="list_of_geom_list_polygons"
        return result_geom_list,result_geom2_list


    def export(self,infile,option=None):
    
        ext=os.path.splitext(infile)[::-1][0]

        if ext == ".shp": 
            sign=1
        elif ext == ".kml":
            sign=-1
        else:
            print("The file type is not supported currently.")
            return

        if self.output_type=="csv_list_line":
            geom=geo_geom.create_lineString([latlon[::sign] for latlon in self.output])
            geo_io.export_vector_from_geom(infile,geom,None)

        elif self.output_type=="csv_list_polygon":
            geom=geo_geom.create_polygon([latlon[::sign] for latlon in self.output])
            geo_io.export_vector_from_geom(infile,geom,None)

        elif self.output_type=="csv_list_polygons":
            geom_list=geo_geom.create_polygons([[latlon[::sign] for latlon in polygon] for polygon in self.output])
            geo_io.export_vector_from_geomList(infile,geom_list,self.output_attribute)

        elif self.output_type=="list_of_geom_list_polygons":
            index=option
            geom_list=self.output[index]
            geo_io.export_vector_from_geomList(infile,geom_list,self.output_attribute[index])
