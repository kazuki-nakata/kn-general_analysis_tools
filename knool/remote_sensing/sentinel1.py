import os
from osgeo import gdal
import xml.etree.ElementTree as ET
import numpy as np
import math
import datetime
from datetime import datetime as dt
import astropy.time
import matplotlib.pyplot as plt
import io
import os
import zipfile
import traceback

class read_iw_slc():
    
    __light_speed = 299792458.0

    def __init__(self, zip_file):

        self.subswath_path = []
        self.zip_file = zip_file
        
        self.ads_header = []
        self.quality_information = []
        self.general_annotation = []
        self.image_annotation = []
        self.doppler_centroid = []
        self.antenna_pattern = []
        self.swath_timing = []
        self.geolocation_grid = []
        self.coordinate_conversion = []
        self.swath_merging  = []
        self.lines_per_burst = []
        self.samples_per_burst = []
        self.burst_list = []
        self.burst_count = []
        self.azimuth_time_interval = []
        self.azimuth_time_list = []
        self.radar_frequency = []
        self.azimuth_steering_rate = []
        self.range_sampling_rate = []
        self.slant_range_time = []         
        
        
        with zipfile.ZipFile(zip_file, 'r') as zip_data:
            infos = zip_data.infolist()
            xml_list = []
            
            for i,info in enumerate(infos):
                _, ext = os.path.splitext(info.filename)
                find = info.filename.find("annotation/s1a")
                if ext == ".xml" and find != -1:
                    xml_list.append(info.filename)
                elif ext == ".tiff":
                    self.subswath_path.append(info.filename)
            xml_list.sort()
            self.subswath_path.sort()
            
            for xml in xml_list:
                xml_data = io.BytesIO(zip_data.read(xml))
                tree = ET.parse(xml_data)
                root = tree.getroot()
                self.build(root)
                
    def build(self, root):
       
        
        self.ads_header.append(root.find("adsHeader"))
        self.quality_information.append(root.find("qualityInformation"))
        self.general_annotation.append(root.find("generalAnnotation"))
        
        image_annotation = root.find("imageAnnotation")        
        self.image_annotation = image_annotation
        
        self.doppler_centroid.append(root.find("dopplerCentroid"))
        self.antenna_pattern.append(root.find("antennaPattern"))
        
        swath_timing = root.find("swathTiming")
        self.swath_timing.append(swath_timing)
        
        self.geolocation_grid.append(root.find("geolocationGrid"))
        self.coordinate_conversion.append(root.find("coordinateConversion"))
        self.swath_merging.append(root.find("swathMerging"))
        self.lines_per_burst.append(int(swath_timing.find("linesPerBurst").text))
        self.samples_per_burst.append(int(swath_timing.find("samplesPerBurst").text))

        burst_list = swath_timing.find("burstList")
        self.burst_list.append(burst_list)
        self.burst_count.append(int(burst_list.get("count")))

        for value in root.iter('azimuthTimeInterval'):
            azimuth_time_interval = float(value.text)
        self.azimuth_time_interval.append(azimuth_time_interval)
        
        azimuth_time_list = []
        for burst in burst_list:
            azimuth_time_list.append(dt.strptime(burst.find("azimuthTime").text, '%Y-%m-%dT%H:%M:%S.%f'))
        
        self.azimuth_time_list.append(azimuth_time_list)
        
        for value in root.iter('radarFrequency'):
            radar_frequency = float(value.text)

        self.radar_frequency.append(radar_frequency)
        
        for value in root.iter('azimuthSteeringRate'):
            azimuth_steering_rate = float(value.text)
        
        self.azimuth_steering_rate.append(azimuth_steering_rate)

        for value in root.iter('rangeSamplingRate'):
            range_sampling_rate = float(value.text)
        self.range_sampling_rate.append(range_sampling_rate)
 
        slant_range_time = float(image_annotation.find("imageInformation").find("slantRangeTime").text)#.iter("slantRangeTime")
        self.slant_range_time.append(slant_range_time)

        
    def read_as_array(self,sw_number):
        with zipfile.ZipFile(self.zip_file, 'r') as zip_data:
            filename = self.subswath_path[sw_number]
            path = "/vsizip" + os.sep + self.zip_file + os.sep + filename
            print(path)
            slc = gdal.Open(path)
            
#            slc = gdal.Open(io.BytesIO(zip_data.read(self.subswath_path[sw_number])))
            
        return slc.ReadAsArray()
        
    def select_burst(self,sw_number):
        try:
            burst_count = self.burst_count[sw_number]
            burst_list = self.burst_list[sw_number]
            fvs_count = int(burst_list[0].find("firstValidSample").get("count"))
            lvs_count = int(burst_list[0].find("lastValidSample").get("count"))
  
            first_valid_sample_array = np.zeros((fvs_count,burst_count), dtype=np.int32)
            last_valid_sample_array = np.zeros((lvs_count, burst_count), dtype=np.int32)

            for i,burst in enumerate(burst_list):
                first_valid_sample_array[:,i] = burst.find("firstValidSample").text.split(" ")
                last_valid_sample_array[:,i] = burst.find("lastValidSample").text.split(" ")

            print(last_valid_sample_array.shape)
            print(first_valid_sample_array.shape)

        except Exception as e:
            print("error:\n" + traceback.format_exc())
            pass
        
        return first_valid_sample_array,last_valid_sample_array
    
    def calc_ks(self,sw_number):
        #------Set Parameters for "Doppler Centroid Rate Introduced by the Scanning of the Antenna (ks)"---------------------        
        burst_count = self.burst_count[sw_number]
        orbit_list = self.general_annotation[sw_number].find("orbitList")
        orbit_count = int(self.general_annotation[sw_number].find("orbitList").get("count"))
        azimuth_time_interval = self.azimuth_time_interval[sw_number]
        lines_per_burst = self.lines_per_burst[sw_number]
        azimuth_time_list = self.azimuth_time_list[sw_number]
        azimuth_steering_rate = self.azimuth_steering_rate[sw_number]
        radar_frequency = self.radar_frequency[sw_number]
        
        t_array = np.zeros((orbit_count),dtype=np.float64)
        v_array = np.zeros((orbit_count),dtype=np.float64) 
        
        for i, orbit in enumerate(orbit_list):
            text_time = orbit.find("time").text
            date_time = dt.strptime(text_time, '%Y-%m-%dT%H:%M:%S.%f')
            if i == 0:
                offset = int(astropy.time.Time(date_time).jd)
            t_array[i] = (astropy.time.Time(date_time).jd - offset)

            vx = float(orbit.find("velocity").find("x").text)
            vy = float(orbit.find("velocity").find("y").text)
            vz = float(orbit.find("velocity").find("z").text)
            v_array[i] = math.sqrt(vx**2 + vy**2 + vz**2)

         #--------------Fitting-------------------------
        coefs = np.polyfit(t_array,v_array,2)

        azt_array = np.zeros(burst_count,dtype=np.float64)
        vs_array = np.zeros(burst_count,dtype=np.float64)
        
        for i, time in enumerate(azimuth_time_list):
            t_plus = time + datetime.timedelta(seconds = azimuth_time_interval * float(lines_per_burst)/2)
            t = astropy.time.Time(time).jd - offset
            azt_array[i] = t
            vs_array[i] = coefs[0]*t**2+ coefs[1]*t + coefs[2]

#         #-----------create figure-----------------------
#         v_pred = coefs[0]*t_array**2+ coefs[1]*t_array + coefs[2] 
#         fig = plt.figure()
#         ax = fig.add_subplot(1,1,1)
#         ax.scatter(t_array,v_array)
#         ax.scatter(t_array,v_pred)
#         ax.scatter(azt_array,vs_array)
#         plt.show()

        ks_array = np.zeros(burst_count, dtype=np.float64)
        
        for i, vs in enumerate(vs_array):
            ks_array[i] = 2*vs*radar_frequency*azimuth_steering_rate/self.__light_speed*2*math.pi/360
        
        return ks_array

    
    def calc_ka(self,sw_number):
        #------Set Parameters for "Doppler FM Rate (ka)"---------------------
        burst_count = self.burst_count[sw_number]
        general_annotation = self.general_annotation[sw_number]        
        afmrl_count = int(general_annotation.find("azimuthFmRateList").get("count"))
        range_sampling_rate = self.range_sampling_rate[sw_number]
        slant_range_time = self.slant_range_time[sw_number]
        azimuth_time_list = self.azimuth_time_list[sw_number]
        samples_per_burst = self.samples_per_burst[sw_number]
        
        coef_polynomial_array0 = np.zeros((3,afmrl_count),dtype=np.float64)
        slant_range_time_origin_array0 = np.zeros(afmrl_count,dtype=np.float64)
        azimuth_time_fmrate_list = []

        coef_polynomial_array = np.zeros((3,burst_count),dtype=np.float64)
        slant_range_time_origin_array = np.zeros(burst_count,dtype=np.float64)

        for i,value in enumerate(general_annotation.find("azimuthFmRateList")):
            azimuth_time_fmrate_list.append(dt.strptime(value.find("azimuthTime").text, '%Y-%m-%dT%H:%M:%S.%f'))
            coef_polynomial_array0[0:3,i] = value.find("azimuthFmRatePolynomial").text.split(" ")
            slant_range_time_origin_array0[i] = float(value.find("t0").text)

        for i, time in enumerate(azimuth_time_list):
            delta_t = [abs((afmr_time-time).total_seconds()) for afmr_time in azimuth_time_fmrate_list] 
            i_t=np.argmin(delta_t)
            coef_polynomial_array[:,i] = coef_polynomial_array0[:,i_t]
            slant_range_time_origin_array[i] = slant_range_time_origin_array0[i_t]

        ka_array = np.zeros((samples_per_burst,burst_count),dtype=np.float64)
        for i in range(0,burst_count):
            c0,c1,c2=coef_polynomial_array[:,i]
            tau0=slant_range_time - slant_range_time_origin_array[i]
            delta_tau = 1/range_sampling_rate
            ka_array[:,i] = [c0+c1*(tau0+float(j)*delta_tau)+c2*(tau0+float(j)*delta_tau)**2 for j in range(0,samples_per_burst)]
        
        return ka_array

    
    def calc_kt(self, sw_number):
        ka_array = self.calc_ka(sw_number)
        ks_array = self.calc_ks(sw_number)
        
        samples_per_burst = self.samples_per_burst[sw_number]
        burst_count = self.burst_count[sw_number]
        
        kt_array = np.zeros((samples_per_burst,burst_count),dtype=np.float64)
        
        for i in range(0,burst_count):
            kt_array[:,i] = [ka_array[j,i]*ks_array[i]/(ka_array[j,i]-ks_array[i]) for j in range(0,samples_per_burst)]
            
        return kt_array
   

    def calc_dc_frequency(self, sw_number):
        #------Set Parameters for "Doppler Centroid Frequency (fnc)"---------------------
        # for value in doppler_centroid.find("dcEstimateList"):
        #     print(value.tag)
        burst_count = self.burst_count[sw_number]
        doppler_centroid = self.doppler_centroid[sw_number]
        afmrl_count = int(self.general_annotation[sw_number].find("azimuthFmRateList").get("count"))        
        azimuth_time_list = self.azimuth_time_list[sw_number]
        samples_per_burst = self.samples_per_burst[sw_number]
        range_sampling_rate = self.range_sampling_rate[sw_number]
        slant_range_time = self.slant_range_time[sw_number]
        
        dcrl_count = int(doppler_centroid.find("dcEstimateList").get("count"))
        coef_polynomial_array0 = np.zeros((3,afmrl_count),dtype=np.float64)
        slant_range_time_origin_array0 = np.zeros(afmrl_count,dtype=np.float64)
        
        dc_rate_time_list = []

        coef_polynomial_array = np.zeros((3,burst_count),dtype=np.float64)
        slant_range_time_origin_array = np.zeros(burst_count,dtype=np.float64)

        for i,value in enumerate(doppler_centroid.find("dcEstimateList")):
            dc_rate_time_list.append(dt.strptime(value.find("azimuthTime").text, '%Y-%m-%dT%H:%M:%S.%f'))
            coef_polynomial_array0[0:3,i] = value.find("dataDcPolynomial").text.split(" ")
            slant_range_time_origin_array0[i] = float(value.find("t0").text)
#            print(value.find("t0").text,dc_rate_time_list[i])

        for i, time in enumerate(azimuth_time_list):
            delta_t = [abs((dcr_time-time).total_seconds()) for dcr_time in dc_rate_time_list] 
            i_t=np.argmin(delta_t)
            coef_polynomial_array[:,i] = coef_polynomial_array0[:,i_t]
            slant_range_time_origin_array[i] = slant_range_time_origin_array0[i_t]
#            print(time)

        fnc_array = np.zeros((samples_per_burst,burst_count),dtype=np.float64)
        for i in range(0,burst_count):
            c0,c1,c2 = coef_polynomial_array[:,i]
            tau0 =  slant_range_time - slant_range_time_origin_array[i]
#            print(tau0,slant_range_time,slant_range_time_origin_array[i])
            delta_tau = 1/range_sampling_rate
            
            fnc_array[:,i] = [c0+c1*(tau0+float(j)*delta_tau)+c2*(tau0+float(j)*delta_tau)**2 for j in range(0,samples_per_burst)]

        return fnc_array

    
    def calc_rzd_azimuth_time(self, sw_number):
        #------Set Parameters for "Reference zero-Doppler Azimuth Time"---------------------
        samples_per_burst = self.samples_per_burst[sw_number]
        burst_count = self.burst_count[sw_number]
        ka_array = self.calc_ka(sw_number)
        fnc_array = self.calc_dc_frequency(sw_number)
        
        eta_ref_array = np.zeros((samples_per_burst,burst_count),dtype=np.float64)
        n = int(samples_per_burst/2)
        
        for i in range(0,burst_count):
            eta_ref_array[:,i] = [ -fnc_array[j,i]/ka_array[j,i] + fnc_array[n,i]/ka_array[n,i] for j in range(0,samples_per_burst)]
            
        return eta_ref_array

    
    def calc_zd_azimuth_time(self, sw_number):    
        #-----Set parameters "zero-doppler azimuth time at each line"-------------------------------
        lines_per_burst = self.lines_per_burst[sw_number]
        azimuth_time_interval = self.azimuth_time_interval[sw_number]
        
        eta_array = np.zeros(lines_per_burst,dtype=np.float64)
        
        for i in range(0,lines_per_burst):
            eta_array[i] = (float(i) - float(lines_per_burst)/2) * azimuth_time_interval

        return eta_array
    
    def calc_all_deramp_parameters(self, sw_number): 
        kt_array=self.calc_kt(sw_number)
        eta_array=self.calc_zd_azimuth_time(sw_number)
        eta_ref_array=self.calc_rzd_azimuth_time(sw_number)
        fvs,lvs=self.select_burst(sw_number)
        return eta_array,eta_ref_array,kt_array,fvs,lvs
        
    def read_deburst_parameters(self, sw_number):     
        #--------------Set parameters for burst stitching (deburst)--------------------------
        burst_count = self.burst_count[sw_number]
        azimuth_time_list = self.azimuth_time_list[sw_number]
        azimuth_time_interval = self.azimuth_time_interval[sw_number]
        lines_per_burst = self.lines_per_burst[sw_number]
        
        delta_p = np.zeros(burst_count,dtype=np.int32)
        for i in range(1,burst_count):
            delta_t = azimuth_time_list[i]-azimuth_time_list[i-1]
            delta_p[i] = int(delta_t.total_seconds()/azimuth_time_interval - lines_per_burst)
            print(delta_p[i])

        db_size = sum(lines_per_burst + delta_p )# + lines_per_burst
        
        fvs,lvs=self.select_burst(sw_number)    
        
        return fvs,lvs,delta_p,db_size

    def deramp(self,slc_array,eta_array,eta_ref_array,kt_array,fvs,lvs):
        deramp_slc = deramp_process(slc_array,eta_array,eta_ref_array,kt_array,fvs,lvs)
        return deramp_slc

    def deburst(self,slc_array,fvs,lvs,delta_p,db_size):
        deburst_slc = deburst_process(slc_array,delta_p,fvs,lvs,db_size)
        return deburst_slc