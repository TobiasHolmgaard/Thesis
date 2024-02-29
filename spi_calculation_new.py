# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:09:01 2023

@author: tojo1
"""
## import libraries
import os
import xarray as xr
import spei as si
import pandas as pd
import scipy.stats as scs
import numpy as np
import matplotlib.pyplot as plt

## change directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_1980_2022_precip.nc', engine = 'netcdf4')

#%%
## Extract values from dataset
time = E_obs['time'][:]
longitude = E_obs['longitude'][180:185]
latitude = E_obs['latitude'][106:109]
precip = E_obs['rr'][:,106:109,180:185]


#%%
## calculate SPI with daily time series

### some error in this part to fix!!!! ###

## define perdiod for calculating values
days = 30
t = days - 1

## create array for resulting values
spi_daily_values = np.empty((len(time), len(latitude), len(longitude)))

## Loop through cell to treat them as individual time series
for lat_idx in range(len(latitude)):
    for lon_idx in range(len(longitude)):
        # Extract the time series for the current grid cell
        time_series = precip[:, lat_idx, lon_idx].to_dataframe()
        ## Caluclate fraction of nan-values
        nans = (len(time_series['rr']) - len(time_series['rr'].dropna())) / len(time_series['rr'])
        ## Select the correct data
        filteret_time_series = time_series['rr']
        ## Roll data for chosen period
        rolled_time_series =  filteret_time_series.rolling(days, min_periods = days).sum()
        rolled_time_series = rolled_time_series.dropna()
        ## calculate SPI if there is not too much data missing
        spi_calc = si.spi(rolled_time_series, dist=scs.gamma)
        spi_daily_values[:t,lat_idx,lon_idx] = np.nan
        spi_daily_values[t:,lat_idx,lon_idx] = spi_calc

## Write data to xarray dataarray
spi_daily_array = xr.Dataset(
    {'SPI' : xr.DataArray(
        data = np.array(spi_daily_values),
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': time},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})


#%%

# calculate monthly SPI-values from daily time series

## define perdiod for calculating values
months = 1
t = months - 1

## resample time variable for monthly values
time_months = time.resample(time="1M").max(dim="time")

## create array for resulting values
spi_monthly_values = np.empty((len(time_months), len(latitude), len(longitude)))

## Loop through cell to treat them as individual time series
for lat_idx in range(len(latitude)):
    for lon_idx in range(len(longitude)):
        ## Extract the time series for the current grid cell
        time_series = precip[:, lat_idx, lon_idx].to_dataframe()
        ## Caluclate fraction of nan-values
        nans = (len(time_series['rr']) - len(time_series['rr'].dropna())) / len(time_series['rr'])
        ## Select the correct data
        filteret_time_series = time_series['rr']
        ## convert daily values to monthly values
        monthly_time_series = filteret_time_series.resample('M').sum()
        ## Roll data for chosen period
        rolling_time_series =  monthly_time_series.rolling(months, min_periods = months).sum()
        rolling_time_series = rolling_time_series.dropna()
        ## calculate SPI if there is not too much data missing
        spi_calc = si.spi(rolling_time_series, dist=scs.gamma)
        spi_monthly_values[:t,lat_idx,lon_idx] = np.nan
        spi_monthly_values[t:,lat_idx,lon_idx] = spi_calc


## Write data to xarray dataarray
spi_monthly_array = xr.Dataset(
    {'SPI' : xr.DataArray(
        data = np.array(spi_monthly_values),
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': time_months},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

#%%

## Compare data to downloaded SPI-values from KNMI - only works for data covering Netherlands
## chose data area in the sections above before running this section

## change working directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
## open dataset
knmi_spi = xr.load_dataset('spi_909_timesteps.nc', engine = 'netcdf4')

knmi_time = knmi_spi['time']

## calculate coordinates in lat/lon for KNMI dataset
longis = (1 - ((knmi_spi['x'][299] - knmi_spi['x'][:]) / knmi_spi['x'][299])) * (7.57697351779 - 3.20584717258) + 3.20584717258
latis = ((1 - ((knmi_spi['y'][349] - knmi_spi['y'][:]) / (knmi_spi['y'][349] - knmi_spi['y'][0]))) * (53.7226030601 - 50.5847875746)) + 50.5847875746

## Give coordinates to KNMI dataset
knmi_spi['x'] = longis
knmi_spi['y'] = latis

## Roll data to get a coarser resolution
knmi_spi_rolling = knmi_spi['sp1'].rolling(x=20, center = True).mean().rolling(y=25, center=True).mean()

## create dataset to assign coarser KNMI resolution SPI values
knmi_coarser = xr.Dataset(
    {'SPI' : xr.DataArray(
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': knmi_time},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

## Loop through coarser dataset to assign values from rolled high-resolution data
for lat_index in range(len(latitude)):
    for lon_index in range(len(longitude)):
        y_1 = knmi_spi_rolling['y'].sel(y=latitude[lat_index], method='nearest').values
        x_1 = knmi_spi_rolling['x'].sel(x=longitude[lon_index], method='nearest').values
        selected_data = knmi_spi_rolling.sel(x = x_1, y = y_1, method = 'nearest').values
        knmi_coarser['SPI'][:,lat_index,lon_index] = selected_data

t = len(knmi_spi_rolling['time'])

knmi_error = spi_daily_array['SPI'][15706-16-t:15706-16,:,:] - knmi_coarser['SPI']

knmi_error.plot()
plt.title('Distribution of error values - KNMI SPI and calculated 80-22')
plt.xlabel('SPI error')
plt.ylabel('count')
plt.show()

knmi_error[:,1,1].plot()
plt.title('KNMI SPI and calculated SPI difference single cell 80-22')
plt.ylabel('SPI difference - daily')
plt.show()

mini = 800

plt.plot(time[15706-16-t:15706-16-mini], spi_daily_array['SPI'][15706-16-t:15706-16-mini,1,1], label = 'Python E_obs', linewidth=0.8)
plt.plot(time[15706-16-t:15706-16-mini], knmi_coarser['SPI'][:909-mini,1,1], label = 'SPI KNMI', linewidth=0.8)
plt.title('Data from 1980 to 2022')
plt.xticks(rotation = 15)
plt.legend()
plt.show()


#%%

## Monthly values for camparison
 
## Compare data to downloaded SPI-values from KNMI - only works for data covering Netherlands
## chose data area in the sections above before running this section

## change working directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
## open dataset
knmi_spi = xr.load_dataset('spi_909_timesteps.nc', engine = 'netcdf4')

knmi_spi_months = knmi_spi.resample(time="1M").mean(dim="time")

knmi_time_monthly = knmi_spi_months['time']

## calculate coordinates in lat/lon for KNMI dataset
longis = (1 - ((knmi_spi['x'][299] - knmi_spi['x'][:]) / knmi_spi['x'][299])) * (7.57697351779 - 3.20584717258) + 3.20584717258
latis = ((1 - ((knmi_spi['y'][349] - knmi_spi['y'][:]) / (knmi_spi['y'][349] - knmi_spi['y'][0]))) * (53.7226030601 - 50.5847875746)) + 50.5847875746

## Give coordinates to KNMI dataset
knmi_spi_months['x'] = longis
knmi_spi_months['y'] = latis

## Roll data to get a coarser resolution
knmi_spi_rolling_monthly = knmi_spi_months['sp1'].rolling(x=20, center = True).mean().rolling(y=25, center=True).mean()

## create dataset to assign coarser KNMI resolution SPI values
knmi_coarser_monthly = xr.Dataset(
    {'SPI' : xr.DataArray(
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': knmi_time_monthly},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

#spi_monthly_array_shifted = spi_monthly_array.assign_coords(time = spi_monthly_array.time + pd.Timedelta(days=15))

## Loop through coarser dataset to assign values from rolled high-resolution data
for lat_index in range(len(latitude)):
    for lon_index in range(len(longitude)):
        y_1 = knmi_spi_rolling_monthly['y'].sel(y=latitude[lat_index], method='nearest').values
        x_1 = knmi_spi_rolling_monthly['x'].sel(x=longitude[lon_index], method='nearest').values
        selected_data = knmi_spi_rolling_monthly.sel(x = x_1, y = y_1, method = 'nearest').values
        knmi_coarser_monthly['SPI'][:,lat_index,lon_index] = selected_data
    
t = len(knmi_spi_rolling_monthly['time'])

knmi_error_monthly = spi_monthly_array['SPI'][876-t:876,:,:] - knmi_coarser_monthly['SPI']

knmi_error_monthly.plot()
plt.title('Distribution of error values - KNMI SPI and calculated monthly 50-22')
plt.xlabel('SPI error')
plt.ylabel('count')
plt.show()

knmi_error_monthly[:,1,1].plot()
plt.title('KNMI SPI and calculated SPI monthly difference single cell 50-22')
plt.ylabel('SPI difference - daily')
plt.show()

mini = 0

plt.plot(time_months[876-t:876-mini], spi_monthly_array['SPI'][876-t:876-mini,1,1], label = 'Python E_obs', linewidth=0.8)
plt.plot(time_months[876-t:876-mini], knmi_coarser_monthly['SPI'][:31-mini,1,1], label = 'SPI KNMI', linewidth=0.8)
plt.title('Data from 1950 to 2022')
plt.xticks(rotation = 25)
plt.legend()
plt.show()


