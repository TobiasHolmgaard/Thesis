# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:07:58 2023

@author: tojo1
"""
#%%
## Cell 1 - calculate SPI value sin python

## import libraries
import os
import xarray as xr
import spei as si
import pandas as pd
import scipy.stats as scs
import numpy as np
import matplotlib.pyplot as plt

### Time to test with longer timeseries to see if that is causing the error...


## Set working directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')

E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_precip_1950_2022.nc', engine = 'netcdf4')

file = E_obs#[:,99:119,182:202]

## Load files Downloaded from Copernicus CDS
# file1 = xr.load_dataset('ERA5_precip_1950_2023_5_x_5_degree_box.grib', 
#                         engine = 'cfgrib', filter_by_keys={'stepType': 'avgas'})
# file2 = xr.load_dataset('ERA5_precip_1950_2023_5_x_5_degree_box.grib',
#                         engine = 'cfgrib', filter_by_keys={'stepType': 'avgad'})

# ## Merge files, as keys change by year 2000
# file = xr.merge([file1, file2])

#%%

## Define variables
time = file['time'][:]
longitude = file['longitude'][174:193]
longitude_panda = longitude.to_dataframe()
latitude = file['latitude'][101:114]
latitude_panda = latitude.to_dataframe()
precip = file['rr'][:,101:114,174:193]
#precip = precip * 1000 * 30

monthly = 'y'

## Set the lenght of the period for SPI calculation
memory = 60
if monthly == 'y':
    m = int(memory / 30)
    time = time.resample(time='M').mean(dim='time')['time']
    m1 = memory-1
else:
    m = memory - 1

## Create dataarray for results
if monthly == 'y':
    spi_values = np.empty((len(time-m), len(latitude), len(longitude)))
else:
    spi_values = np.empty((len(time-m), len(latitude), len(longitude)))

## Define lenght for looping through grid
total_lats = len(latitude)
total_lons = len(longitude)

## set counter
x = 0

error = []

locator = []

## Loop through every pixel to calculate index-values
for lat_index in range(total_lats):
    for lon_index in range(total_lons):
        locator.append([x, lat_index,lon_index])
        ## get the values for the lat/lon grid cell
        precip_values = precip[:,lat_index, lon_index]
        ## convert to dataframe, so library function works
        precip_panda = precip_values.to_dataframe()
        ## Roll to create a cross-monthyl dataset
        if monthly == 'y':
            precip_panda = precip_panda.resample('M').sum()
            
        precip_panda_roll = precip_panda.rolling(memory, min_periods = memory).sum()
            
        ## drop nan values
        precip_panda_roll = precip_panda_roll.dropna()
        
        if len(precip_panda_roll['rr']) > 0:
        
        ## Calculate SPI on chosen statistical distribution
            try:
                spi = si.spi(precip_panda_roll['rr'], dist=scs.gamma)
            except:
                error.append(x)
                continue
                
            #monthly_spi = spi.resample('M').mean()
            monthly_spi = spi
            ## Write values to datarray
            spi_values[:m,lat_index, lon_index] = np.nan
            if monthly == 'y':
                spi_values[m1:, lat_index, lon_index] = monthly_spi
            else:
                spi_values[m:, lat_index, lon_index] = spi
        else:
            continue
        ## Write precipitation data to csv for comparison with other software
        filename = 'spi_test/test_data' + str(x) + '.csv'
        precip_panda.to_csv(filename, sep = ',',columns = ['rr'])
        ## update counter
        print(x)
        x = x + 1

## Convert dataarray to Xarray
spi_array = xr.Dataset(
    {'SPI' : xr.DataArray(
        data = np.array(spi_values),
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': time},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})
##### USE THE SPI CALCULATOR FOR GENERATING COMPARISON VALUES #####


#%%
## Cell 2 - comparing SPI calculator and library

if monthly == 'y':
    m = m1

## Load data for comparison with software
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
test_data = pd.read_excel('spi_test/spi_test_results/test_data57_SPI_M_02.xlsx', header=1)
tester1 = test_data['spi2'][m:]

tester = spi_array['SPI'][m:,3,2]

times = time[m:]

## Plot the two datasets together
plt.plot(times[100:200], tester[100:200], label = 'Python', linewidth=0.8)
plt.plot(times[100:200], tester1[100:200], label = 'SPI_generator', linewidth=0.8)
plt.legend()
plt.title('Gamma distribution 30 days')
plt.show()

error = tester - tester1
error_sum =  error.sum().values
error_mean = error.mean().values
error_std = error.std().values

## Plot the difference between the two datasets as a graph
plot_text = 'sum: ' + str(np.round(error_sum,4)) + 'mean: ' + str(np.round(error_mean,4)) + 'std:' + str(np.round(error_std,4))
error[100:200].plot()
plt.title('Gamma distribution 30 days')
#plt.text(2, 6, 'Text Here', fontsize=12, color='red')
plt.text(2000, 1, plot_text, fontsize=10, color='blue',
         bbox=dict(facecolor='white', edgecolor='blue', boxstyle='round'))
plt.show()

error_sum =  error.sum().values
error_mean = error.mean().values
error_std = error.std().values

#### longer timeseries gives better results and lower error values
#### I suspect the data distrubution is given in the software already

#spi_array.to_netcdf('spi_3_years.nc')

#reference_data = xr.open_dataset('esf01_f_euu_20220101_20221201_m.nc', engine = 'netcdf4')

#%%

## Cell 3 - comparing python values with drought values from MARS - doesn't really make sense
## MARS is not only SPI, but is a mixed index with other parameters also included.

Mars = xr.load_dataset('SPI_mars_reference_2018_2020.nc', engine = 'netcdf4')
Mars_dry = Mars['spg01'][:,140:145,185:190]
tester_mars = spi_array[(816):(852-m),:20,:20]
tester_mars_resample = tester_mars.coarsen(latitude=4).mean().coarsen(longitude=4).mean()
tester_mars_resample = tester_mars_resample.rename(latitude='lat')
tester_mars_resample['lat'] = tester_mars_resample['lat'] - 0.125
tester_mars_resample = tester_mars_resample.reindex({'lat': tester_mars_resample['lat'].values[::-1]})
tester_mars_resample = tester_mars_resample.rename(longitude='lon')
tester_mars_resample['lon'] = tester_mars_resample['lon'] + 0.125
one_day = pd.Timedelta(hours=6)
tester_mars_resample['time'] = tester_mars_resample['time'] + one_day

error3d = Mars_dry - tester_mars_resample
Mars_error_sum = error3d.sum().values
Mars_error_mean = error3d.mean().values

#%%

## cell 4 comparing python values from E_obs with SPI from KNMI'

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

knmi_spi = xr.load_dataset('spi_909_timesteps.nc', engine = 'netcdf4')

longis = (1 - ((knmi_spi['x'][299] - knmi_spi['x'][:]) / knmi_spi['x'][299])) * (7.57697351779 - 3.20584717258) + 3.20584717258

latis = ((1 - ((knmi_spi['y'][349] - knmi_spi['y'][:]) / (knmi_spi['y'][349] - knmi_spi['y'][0]))) * (53.7226030601 - 50.5847875746)) + 50.5847875746

knmi_spi['x'] = longis

knmi_spi['y'] = latis

knmi_spi_rolling = knmi_spi['sp1'].rolling(x=20, center = True).mean().rolling(y=25, center=True).mean()

t = len(knmi_spi_rolling['time'])

knmi_coarser = xr.Dataset(
    {'SPI' : xr.DataArray(
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': time},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

knmi_coarser['SPI'][:,:,:] = np.nan

knmi_coarser = knmi_coarser.sel(time=slice(knmi_spi['time'][0],knmi_spi['time'][-1]))

for lat_index in range(total_lats):
    for lon_index in range(total_lons):
        y_1 = knmi_spi_rolling['y'].sel(y=latitude[lat_index], method='nearest').values
        x_1 = knmi_spi_rolling['x'].sel(x=longitude[lon_index], method='nearest').values
        selected_data = knmi_spi_rolling.sel(x = x_1, y = y_1, method = 'nearest').values
        knmi_coarser['SPI'][:,lat_index,lon_index] = selected_data

#knmi_coarser = knmi_coarser.sel(time=slice(knmi_spi['time'][0],knmi_spi['time'][-1]))

knmi_coarser['SPI'].plot(label = 'KNMI')
plt.title('KNMI SPI values and calculated values')
spi_array['SPI'][15706-16-t:15706-16,:,:].plot(label = 'My Calc')
plt.legend()
plt.show()

knmi_error = spi_array['SPI'][15706-16-t:15706-16,:,:] - knmi_coarser['SPI']

knmi_error.plot()
plt.title('KNMI SPI vs calculated SPI values distrubution')
plt.show()

knmi_error[:,7,7].plot()
plt.title('KNMI SPI vs calculated SPI for single cell data time series')
plt.show()

plt.plot(time[15706-16-t:15706-16], spi_array['SPI'][15706-16-t:15706-16,7,7], label = 'Python', linewidth=0.8)
plt.plot(time[15706-16-t:15706-16], knmi_coarser['SPI'][:,7,7], label = 'SPI_generator', linewidth=0.8)
plt.legend()
plt.title('Gamma distribution 30 days')
plt.show()
