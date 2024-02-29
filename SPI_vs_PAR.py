# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:56:29 2023

@author: tojo1
"""

import cartopy.crs as ccrs
#from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
import numpy.polynomial.polynomial as poly

#%%
months_number = 1

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')
print('SPI open')

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc', engine = 'netcdf4')

zero_value_mask = E_obs['rr'].mean(dim = 'time').copy()

for lat_idx in range(len(spi_monthly_array['latitude'])):
    for lon_idx in range(len(spi_monthly_array['longitude'])):
        zeroser = (E_obs['rr'][:,lat_idx,lon_idx] == 0).sum(dim='time')
        if zeroser > 0.22 * len(spi_monthly_array['time']):
            zero_value_mask[lat_idx, lon_idx] = np.nan
        else: 
            zero_value_mask[lat_idx, lon_idx] = 1
            
spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')

## Combined droght classifcation with destintive classes
# combined_index = xr.load_dataset('EDO_CDI_2012_2023.nc')

## PAR anomaly index for vegetation drought indication relative to std
# par_index = xr.load_dataset('PAR_EDO_2002_2022.nc', chunks={'lat': 250, 'lon': 250})
# print('PAR index open')

## Modelled soil moisture index from ECMWF LISFLOOD
# soil_moisture_index = xr.load_dataset('EDO_SM_1995_2022.nc')

years = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012,
       2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022]

folder_path = r"C:\Users\tojo1\Documents\Speciale\Data\drought_index\PAR_EDO"

# # List all files in the folder
# files = os.listdir(folder_path)

# par_index = xr.DataArray(np.nan, coords={'time': spi_monthly_array['time'][:], 'latitude': spi_monthly_array['latitude'], 'longitude': spi_monthly_array['longitude']}, dims=('time', 'latitude', 'longitude'))

# par_all_years = xr.DataArray(
#     par_index,
#     coords={'time' : spi_monthly_array['time'], 'latitude': spi_monthly_array['latitude'] , 'longitude' : spi_monthly_array['longitude']},  # Use your own coordinate names and values
#     dims=('time', 'latitude', 'longitude'),
#     attrs={'description': 'PAR anomaly'},)

# x = 624

# # Loop through each file in the folder and gather data in one xarray dataset
# for file_name in files:
#     # Do something with each file, for example, print the file name
#     print(file_name)
#     par_index = xr.load_dataset('PAR_EDO/' + file_name)
#     par_index_rolling = par_index.rolling(lon=6, center = True).mean().rolling(lat=6, center=True).mean()
#     print('rolling done')
#     par_index_rolling['lat'] = np.round(par_index['lat'].rolling(lat=6,center = True).mean().values, decimals=3)
#     par_index_rolling['lon'] = np.round(par_index['lon'].rolling(lon=6,center = True).mean().values, decimals=3)
#     print('lat and lon rounded')
#     par_index_mean = par_index_rolling.resample(time="1M").mean(dim= 'time')
#     condition = par_index_mean.isin(spi_monthly_array['SPI'])
#     print('condiion made')
#     par_index_rolling_selected = par_index_mean.where(par_index_mean['time'] == spi_monthly_array['time'], drop = True)
#     par_index_rolling_selected = par_index_rolling_selected.where(par_index_rolling_selected['lat'].isin(spi_monthly_array['latitude']), drop = True)
#     par_index_rolling_selected = par_index_rolling_selected.where(par_index_rolling_selected['lon'].isin(spi_monthly_array['longitude']), drop = True)
#     print('data selected')
#     par_all_years[x:x+len(par_index_rolling_selected['time']),11:187,34:362] = par_index_rolling_selected['fapan']
#     x = x+12
    
# par_index_diff = spi_monthly_array['SPI'] - par_all_years

# os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
# par_all_years.to_netcdf('PAR_values_scaled_for_SPI.nc')

#%% 
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')

par_all_years = xr.open_dataset('PAR_values_scaled_for_SPI.nc')

par_all_years = par_all_years.where(zero_value_mask == 1, drop = True)

par_index_diff = spi_monthly_array['SPI'] - par_all_years['__xarray_dataarray_variable__']

## plot the difference between SPI and PAR for a specific time
## set time index:
ind = 700

par_index_diff = par_index_diff.dropna(dim='latitude', how = 'all')
par_index_diff = par_index_diff.dropna(dim='longitude', how = 'all')

## plot the data
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1.00, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
par_index_diff[ind].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
ax.set_title('SPI and PAR anomaly difference ' + str(par_index_diff['time'][ind].values)[:7], fontsize=25)
plt.show()

#%%

## plot the time mean of the differences spatially

par_index_diff_mean = par_index_diff.mean(dim = 'time')

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1.00, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
par_index_diff_mean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': 'SPI - PAR'})
ax.set_title('SPI and PAR anomaly difference mean', fontsize=25)
plt.show()

## plot the sum of the differences

par_index_diff_sum = par_index_diff.sum(dim = 'time')

par_index_diff_sum = par_index_diff_sum.where(~par_index_diff_mean.isnull())

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-100, 100, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
par_index_diff_sum.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': 'SPI - PAR'})
ax.set_title('SPI and PAR anomaly difference sum', fontsize=25)
plt.show()

#%%

months = [5,6,7,8,9]

months_number = 1

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

## SPI and PAR correlation maps

selected_SPI_data = spi_monthly_array['SPI'][624:]
selected_SPI_data = selected_SPI_data.sel(time=selected_SPI_data['time.month'].isin(months))
selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]
selected_PAR_data = selected_PAR_data.sel(time=selected_PAR_data['time.month'].isin(months))

correlation = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

regressions = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

x = 0

for lat in range(len(spi_monthly_array['latitude'])):
    for lon in range(len(spi_monthly_array['longitude'])):
        spi_data = selected_SPI_data[:,lat,lon].values
        # spi_data[spi_data == np.inf] = np.nan
        par_data = selected_PAR_data[:,lat,lon]
        nan_count = np.isnan(par_data).sum()
        if np.all(np.isnan(par_data)) ==  True:
            regressions[lat, lon] = np.nan
        elif np.isnan(spi_data).any() ==  True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif np.any(np.isinf(spi_data)) == True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif (nan_count/len(par_data)) > 0.25:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        else:
            valid_indices = np.arange(len(par_data))[~np.isnan(par_data)]
            par_data = np.interp(np.arange(len(par_data)), valid_indices, 
                                          par_data[valid_indices])
            polifit = poly.polyfit(spi_data[:-3], par_data[:-3],1)
            regressions[lat, lon] = polifit[0]
            corrco = np.corrcoef(spi_data, par_data)
            correlation[lat, lon] = corrco[0,1]

    x = x+1 
    print('SPI and PAr correlation coefficient ' + str(round((x*100/(len(spi_monthly_array['latitude'] ))),1)) + ' %')
## Find correlations on of PC and the SPI-values

# corr_spi_par = correlation_map(spi_monthly_array['SPI'][624:,:,:],par_all_years[624:,:,:].values)

## correlation data as Xarray dataarray

regressions_xr = xr.DataArray(
    regressions,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

correlation_xr = xr.DataArray(
    correlation,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

correlation_xr_clean = correlation_xr.where(correlation_xr != 0)

correlation_xr_clean = correlation_xr_clean.dropna(dim='latitude', how = 'all')
correlation_xr_clean = correlation_xr_clean.dropna(dim='longitude', how = 'all')
        
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
correlation_xr_clean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': ''})
ax.set_title('SPI-' + str(months_number) + ' and PAR correlation coefficient', fontsize=22)
plt.show()

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
correlation_xr_clean.to_netcdf('PAR_and_SPI_' + str(months_number) + '_correlation.nc')

#%%

months_number = 3

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

## SPI and PAR correlation maps
selected_SPI_data = spi_monthly_array['SPI'][624:]
selected_SPI_data = selected_SPI_data.sel(time=selected_SPI_data['time.month'].isin(months))
selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]
selected_PAR_data = selected_PAR_data.sel(time=selected_PAR_data['time.month'].isin(months))
# selected_SPI_data = spi_monthly_array['SPI'][624:]
# selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]

correlation = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

regressions = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

x = 0

for lat in range(len(spi_monthly_array['latitude'])):
    for lon in range(len(spi_monthly_array['longitude'])):
        spi_data = selected_SPI_data[:,lat,lon].values
        # spi_data[spi_data == np.inf] = np.nan
        par_data = selected_PAR_data[:,lat,lon]
        nan_count = np.isnan(par_data).sum()
        if np.all(np.isnan(par_data)) ==  True:
            regressions[lat, lon] = np.nan
        elif np.isnan(spi_data).any() ==  True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif np.any(np.isinf(spi_data)) == True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif (nan_count/len(par_data)) > 0.25:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        else:
            valid_indices = np.arange(len(par_data))[~np.isnan(par_data)]
            par_data = np.interp(np.arange(len(par_data)), valid_indices, 
                                          par_data[valid_indices])
            polifit = poly.polyfit(spi_data[:-3], par_data[:-3],1)
            regressions[lat, lon] = polifit[0]
            corrco = np.corrcoef(spi_data, par_data)
            correlation[lat, lon] = corrco[0,1]

    x = x+1 
    print('SPI and PAr correlation coefficient ' + str(round((x*100/(len(spi_monthly_array['latitude'] ))),1)) + ' %')
## Find correlations on of PC and the SPI-values

# corr_spi_par = correlation_map(spi_monthly_array['SPI'][624:,:,:],par_all_years[624:,:,:].values)

## correlation data as Xarray dataarray

regressions_xr = xr.DataArray(
    regressions,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'regresiion coefficient'},)

correlation_xr = xr.DataArray(
    correlation,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

correlation_xr_clean = correlation_xr.where(correlation_xr != 0)

correlation_xr_clean = correlation_xr_clean.dropna(dim='latitude', how = 'all')
correlation_xr_clean = correlation_xr_clean.dropna(dim='longitude', how = 'all')
        
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
correlation_xr_clean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': ''})
ax.set_title('SPI-' + str(months_number) + ' and PAR correlation coefficient', fontsize=22)
plt.show()

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
correlation_xr_clean.to_netcdf('PAR_and_SPI_' + str(months_number) + '_correlation.nc')

#%%

months_number = 6

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

## SPI and PAR correlation maps
selected_SPI_data = spi_monthly_array['SPI'][624:]
selected_SPI_data = selected_SPI_data.sel(time=selected_SPI_data['time.month'].isin(months))
selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]
selected_PAR_data = selected_PAR_data.sel(time=selected_PAR_data['time.month'].isin(months))
# selected_SPI_data = spi_monthly_array['SPI'][624:]
# selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]

correlation = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

regressions = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

x = 0

for lat in range(len(spi_monthly_array['latitude'])):
    for lon in range(len(spi_monthly_array['longitude'])):
        spi_data = selected_SPI_data[:,lat,lon].values
        # spi_data[spi_data == np.inf] = np.nan
        par_data = selected_PAR_data[:,lat,lon]
        nan_count = np.isnan(par_data).sum()
        if np.all(np.isnan(par_data)) ==  True:
            regressions[lat, lon] = np.nan
        elif np.isnan(spi_data).any() ==  True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif np.any(np.isinf(spi_data)) == True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif (nan_count/len(par_data)) > 0.25:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        else:
            valid_indices = np.arange(len(par_data))[~np.isnan(par_data)]
            par_data = np.interp(np.arange(len(par_data)), valid_indices, 
                                          par_data[valid_indices])
            polifit = poly.polyfit(spi_data[:-3], par_data[:-3],1)
            regressions[lat, lon] = polifit[0]
            corrco = np.corrcoef(spi_data, par_data)
            correlation[lat, lon] = corrco[0,1]

    x = x+1 
    print('SPI and PAr correlation coefficient ' + str(round((x*100/(len(spi_monthly_array['latitude'] ))),1)) + ' %')
## Find correlations on of PC and the SPI-values

# corr_spi_par = correlation_map(spi_monthly_array['SPI'][624:,:,:],par_all_years[624:,:,:].values)

## correlation data as Xarray dataarray

regressions_xr = xr.DataArray(
    regressions,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

correlation_xr = xr.DataArray(
    correlation,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

correlation_xr_clean = correlation_xr.where(correlation_xr != 0)

correlation_xr_clean = correlation_xr_clean.dropna(dim='latitude', how = 'all')
correlation_xr_clean = correlation_xr_clean.dropna(dim='longitude', how = 'all')
        
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
correlation_xr_clean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': ''})
ax.set_title('SPI-' + str(months_number) + ' and PAR correlation coefficient', fontsize=22)
plt.show()

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
correlation_xr_clean.to_netcdf('PAR_and_SPI_' + str(months_number) + '_correlation.nc')

#%%

months_number = 12

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

## SPI and PAR correlation maps

selected_SPI_data = spi_monthly_array['SPI'][624:]
selected_SPI_data = selected_SPI_data.sel(time=selected_SPI_data['time.month'].isin(months))
selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]
selected_PAR_data = selected_PAR_data.sel(time=selected_PAR_data['time.month'].isin(months))

# selected_SPI_data = spi_monthly_array['SPI'][624:]
# selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]

correlation = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

regressions = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

x = 0

for lat in range(len(spi_monthly_array['latitude'])):
    for lon in range(len(spi_monthly_array['longitude'])):
        spi_data = selected_SPI_data[:,lat,lon].values
        # spi_data[spi_data == np.inf] = np.nan
        par_data = selected_PAR_data[:,lat,lon]
        nan_count = np.isnan(par_data).sum()
        if np.all(np.isnan(par_data)) ==  True:
            regressions[lat, lon] = np.nan
        elif np.isnan(spi_data).any() ==  True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif np.any(np.isinf(spi_data)) == True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif (nan_count/len(par_data)) > 0.25:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        else:
            valid_indices = np.arange(len(par_data))[~np.isnan(par_data)]
            par_data = np.interp(np.arange(len(par_data)), valid_indices, 
                                          par_data[valid_indices])
            polifit = poly.polyfit(spi_data[:-3], par_data[:-3],1)
            regressions[lat, lon] = polifit[0]
            corrco = np.corrcoef(spi_data, par_data)
            correlation[lat, lon] = corrco[0,1]

    x = x+1 
    print('SPI and PAr correlation coefficient ' + str(round((x*100/(len(spi_monthly_array['latitude'] ))),1)) + ' %')
## Find correlations on of PC and the SPI-values

# corr_spi_par = correlation_map(spi_monthly_array['SPI'][624:,:,:],par_all_years[624:,:,:].values)

## correlation data as Xarray dataarray

regressions_xr = xr.DataArray(
    regressions,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

correlation_xr = xr.DataArray(
    correlation,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

correlation_xr_clean = correlation_xr.where(correlation_xr != 0)

correlation_xr_clean = correlation_xr_clean.dropna(dim='latitude', how = 'all')
correlation_xr_clean = correlation_xr_clean.dropna(dim='longitude', how = 'all')
        
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
correlation_xr_clean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': ''})
ax.set_title('SPI-' + str(months_number) + ' and PAR correlation coefficient', fontsize=22)
plt.show()

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
correlation_xr_clean.to_netcdf('PAR_and_SPI_' + str(months_number) + '_correlation.nc')

#%%

months_number = 24

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

## SPI and PAR correlation maps

selected_SPI_data = spi_monthly_array['SPI'][624:]
selected_SPI_data = selected_SPI_data.sel(time=selected_SPI_data['time.month'].isin(months))
selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]
selected_PAR_data = selected_PAR_data.sel(time=selected_PAR_data['time.month'].isin(months))

# selected_SPI_data = spi_monthly_array['SPI'][624:]
# selected_PAR_data = par_all_years['__xarray_dataarray_variable__'][624:]

correlation = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

regressions = np.full((len(spi_monthly_array['latitude']), len(spi_monthly_array['longitude'])),np.nan)

x = 0

for lat in range(len(spi_monthly_array['latitude'])):
    for lon in range(len(spi_monthly_array['longitude'])):
        spi_data = selected_SPI_data[:,lat,lon].values
        # spi_data[spi_data == np.inf] = np.nan
        par_data = selected_PAR_data[:,lat,lon]
        nan_count = np.isnan(par_data).sum()
        if np.all(np.isnan(par_data)) ==  True:
            regressions[lat, lon] = np.nan
        elif np.isnan(spi_data).any() ==  True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif np.any(np.isinf(spi_data)) == True:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        elif (nan_count/len(par_data)) > 0.25:
            regressions[lat, lon] = np.nan
            correlation[lat, lon] = np.nan
        else:
            valid_indices = np.arange(len(par_data))[~np.isnan(par_data)]
            par_data = np.interp(np.arange(len(par_data)), valid_indices, 
                                          par_data[valid_indices])
            polifit = poly.polyfit(spi_data[:-3], par_data[:-3],1)
            regressions[lat, lon] = polifit[0]
            corrco = np.corrcoef(spi_data, par_data)
            correlation[lat, lon] = corrco[0,1]

    x = x+1 
    print('SPI and PAr correlation coefficient ' + str(round((x*100/(len(spi_monthly_array['latitude'] ))),1)) + ' %')
## Find correlations on of PC and the SPI-values

# corr_spi_par = correlation_map(spi_monthly_array['SPI'][624:,:,:],par_all_years[624:,:,:].values)

## correlation data as Xarray dataarray

regressions_xr = xr.DataArray(
    regressions,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

correlation_xr = xr.DataArray(
    correlation,
    coords={'latitude' : spi_monthly_array['latitude'], 'longitude' : spi_monthly_array['longitude']},  
    dims=('latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

correlation_xr_clean = correlation_xr.where(correlation_xr != 0)

correlation_xr_clean = correlation_xr_clean.dropna(dim='latitude', how = 'all')
correlation_xr_clean = correlation_xr_clean.dropna(dim='longitude', how = 'all')
        
fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-1, 1, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
correlation_xr_clean.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True,
                         cbar_kwargs={'label': ''})
ax.set_title('SPI-' + str(months_number) + ' and PAR correlation coefficient', fontsize=22)
plt.show()

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\drought_index')
correlation_xr_clean.to_netcdf('PAR_and_SPI_' + str(months_number) + '_correlation.nc')



