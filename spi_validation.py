# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:24:37 2023

@author: tojo1
"""

## import libraries
import os
import xarray as xr
import spei as si
import pandas as pd
import scipy.stats as scs
import numpy as np

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')

test_data = pd.read_csv('spi_test/spi_test_results/test_data0_SPI_M_01_06.csv', header=1)

file = xr.load_dataset('spi_3_years.nc', engine =  'netcdf4')


#file = xr.load_dataset('Monthly_precipitation_data_1981_2020.nc', engine =  'netcdf4')

#data = file['precip'][:360,35:40,160:185]

#time = data['time'][:]
#longitude = data['lon'][:]
#longitude_panda = longitude.to_dataframe()
#latitude = data['lat'][:]
#latitude_panda = file['lat'][:].to_dataframe()
#precip = data['tp'][:]
#precip = precip * 1000
#temp = file['t2m'][:]

memory = 1

#spi_values = np.empty(((len(time)-(memory-1)), len(latitude), len(longitude)))
#x = 0

#total_lats = len(latitude)
#total_lons = len(longitude)

#for lat_index in range(total_lats):
#    for lon_index in range(total_lons):
#        
#        # get the values for the lat/lon grid cell
#        precip_values = data[:,lat_index, lon_index]
#        precip_panda = precip_values.to_dataframe()
#        nan_count = precip_panda['precip'].isnull().sum()
#        if nan_count < 10:
#            precip_panda = precip_panda.rolling(memory, min_periods = memory).sum().dropna()
#    
#            spi = si.spi(precip_panda['precip'], dist=scs.gamma)
#            spi_values[:, lat_index, lon_index] = spi
#        
#        x = x + 1
#
#spi_array = xr.DataArray(
#    data = np.array(spi_values),
#    dims =('time','lat', 'lon',),
#    coords={'time': time,'lat': latitude, 'lon': longitude},
#    attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})

#spi_array[-1,:,:].plot()

#spi_array[:,3,10].plot()

#spi_array.to_netcdf('spi_test.nc')

#reference_data = xr.open_dataset('spg01_m_wld_20000101_20001201_m.nc', engine = 'netcdf4')

#reference = xr.load_dataset('esf01_f_euu_2018_2020.nc', engine =  'netcdf4')
SPI_ref = reference_data['spg01'][:,140:145,160:185]
spi_array_check = spi_array[228:240,:,:]
#SPI_ref['time'] = spi_array_check['time'][:12]

spi_diff = SPI_ref - spi_array_check


