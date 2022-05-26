import sys
import os
import numpy as np
import pandas as pd
from skyfield.api import load
from datetime import datetime as dt
from datetime import timedelta
from skyfield.api import utc
from osgeo import gdal, osr, ogr
from ..geodata_processor import geo_io, geo_info

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
        for date in ds_date:
            lat, lon, ele = self.calc_position_at(date)
            pos_list.append([lat,lon])
        self.output=pos_list
        self.output_type="line"
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
        self.output_type="line"
        return self.output

    def calc_buff_positions_between(self,fdate,ldate,interval,distance,ori="right"):
        #input:
        # ldate, fdate: text -> e.g., 2011-07-01
        # interval: int (min) -> 30 (min)
        ds_date=pd.date_range(start=fdate, end=ldate, freq=str(interval)+"min",tz=utc)
        pos_list=[]
       
        for i, date in enumerate(ds_date[:-1]):
            lat, lon, ele = self.calc_position_at(date)
            
            date2=ds_date[i+1]
            lat2, lon2, ele2 = self.calc_position_at(date2)

            lat3,lon3,ele3=geo_info.calc_line_buffer_point(lat,lon,0,lat2,lon2,0,distance,ori)
            pos_list.append([lat3,lon3])
            
        self.output_type="line"
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
        
        date_ds=pd.date_range(start=fdate, end=ldate, freq=str(interval)+"min",tz=utc)
        date_list=date_ds.to_list()
        date_sf=load.timescale().from_datetimes(date_list)
        
        lat, lon, _ =self._get_position(num, date_sf)
        
        lat2,lon2,_=geo_info.calc_line_buffer_point(lat[0:-1],lon[0:-1],0,lat[1:],lon[1:],0,distance,ori)
        self.output=np.stack([lat2,lon2]).T
        self.output_type="line"
        return self.output

    def calc_buff_area_between(self,fdate,ldate,interval,distance,ori="middle"):

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
        
        result=np.concatenate([rpos_array,lpos_array[::-1]])
        self.output=result
        self.output_type="polygon"
        return self.output   

    def find_scenes(self,fdate,ldate,interval,distance,ori="middle"):
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
        
        
        result=np.vstack([rpos_array,lpos_array[::-1]])
        print(result.shape)
        self.output=result
        self.output_type="polygon"
        return self.output

    def export(self,infile):
        if self.output_type=="line":
            geo_io.create_line_string(infile,self.output)
        elif self.output_type=="polygon":
            geo_io.create_polygon(infile,self.output)