#!/usr/bin/env python
import cdsapi
import os

c = cdsapi.Client()

year=2020
while year<=2020:

    month=1
    while month<=1:
     
        print(str(year),str(month).zfill(2))
        target="era5_"+str(year)+str(month).zfill(2)+".nc"
        
        dmax=31
        if month==2:
            dmax=28
            if year==1980 or year==1984 or year==1988 or year==1992 or year==1996 or year==2000 or year==2004 or year==2008 or year==2012 or year==2016 or year==2020:
                dmax=29
        if month==4 or month==6 or month==9 or month==11:
            dmax=30

        d=list(range(1,dmax+1))
        #print(d)
        dlist=[str(n).zfill(2) for n in d]
        print(dlist)

        vals=['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature', 'mean_sea_level_pressure', 'sea_surface_temperature', 'total_cloud_cover', 'total_column_water_vapour',]
        #vals=['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature', 'mean_sea_level_pressure', 'sea_surface_temperature', 'total_cloud_cover', 'total_column_cloud_liquid_water', 'total_column_water_vapour',]
        print(vals)

        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': vals,
                'year': str(year),
                'month': str(month).zfill(2),
                'day': dlist,
                'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
            },
            target)




        
        month+=1
    year+=1

    
