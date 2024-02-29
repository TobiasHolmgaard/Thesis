# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 16:15:45 2024

@author: tojo1
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 14:39:19 2023

@author: tojo1
"""

### Import libraries
import cartopy.crs as ccrs
#from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
from eofs.standard import Eof
from eofs.tools.standard import correlation_map, covariance_map
import spei as si
import scipy.stats as scs
import time
import numpy.ma as ma
import numpy.polynomial.polynomial as poly
from matplotlib.patches import PathPatch, Path
from scipy.ndimage import uniform_filter
from matplotlib.colorbar import Colorbar
import pandas as pd

start_time = time.time()

#%%
## Calculate SPI-values - import data
## change directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc', engine = 'netcdf4')

print('precip data loaded')
#%%
## Calculate SPI-values on daily time series for monthly SPI 

### OBS! Can be edited for monthly date imput...

time_spi = E_obs['time'][:]
longitude = E_obs['longitude']#[180:222]
latitude = E_obs['latitude']#[100:142]
precip = E_obs['rr']#[:,100:142,180:222]

E_obs.close()
# calculate monthly SPI-values from daily time series

### plot the spatial distributions of standard deviations for precipitation
std_precip = precip.std(dim='time')

mean_precip = precip.mean(dim = 'time')

std_norm_precip = std_precip/mean_precip

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(0, 120, 13)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_precip[15:185,70:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('standard deviation (mm)')
ax.set_title('Std for precipitation in gridcells 1950-2022', fontsize=25)
plt.show()

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(0, 2, 9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_norm_precip[15:185,70:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_norm_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('Std / mean (mm)')
ax.set_title('Std normalized precip in gridcells 1950-2022', fontsize=25)
plt.show()

zero_value_mask = E_obs['rr'].mean(dim = 'time').copy()

for lat_idx in range(len(latitude)):
    for lon_idx in range(len(longitude)):
        zeroser = (E_obs['rr'][:,lat_idx,lon_idx] == 0).sum(dim='time')
        if zeroser > 0.22 * len(time_spi):
            zero_value_mask[lat_idx, lon_idx] = np.nan
        else: 
            zero_value_mask[lat_idx, lon_idx] = 1

mean_precip_mask = mean_precip.where(zero_value_mask == 1, drop = True)

std_precip_mask = std_precip.where(zero_value_mask == 1, drop = True)

std_norm_precip_mask = std_norm_precip.where(zero_value_mask == 1, drop = True)

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(0, 120, 13)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_precip_mask.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('standard deviation (mm)')
ax.set_title('Std for precipitation in gridcells 1950-2022 - masked', fontsize=20)
plt.show()

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(0, 2, 9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_norm_precip_mask.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_norm_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('Std / mean (mm)')
ax.set_title('Std normalized precip in gridcells 1950-2022 - masked', fontsize=20)
plt.show()

months_numbers = [1,3,6,12,24]
## define perdiod for calculating values
# months_number = 24   

months_number = 12
# for months_number in months_numbers:
    
t = months_number - 1

## resample time variable for monthly values
#time_months = time.resample(time="1M").max(dim="time")

## create array for resulting values
spi_monthly_values = np.empty((len(time_spi), len(latitude), len(longitude)))

x = 0

print('time for SPI calculation')

# ## Loop through cell to treat them as individual time series
# for lat_idx in range(len(latitude)):
#     for lon_idx in range(len(longitude)):
#         ## Extract the time series for the current grid cell
#         time_series = precip[:, lat_idx, lon_idx].to_dataframe()
#         ## Caluclate fraction of nan-values
#         nans = (len(time_series['rr']) - len(time_series['rr'].dropna())) / len(time_series['rr'])
#         ## Select the correct data
#         filteret_time_series = time_series['rr']
#         ## convert daily values to monthly values
#         #monthly_time_series = filteret_time_series.resample('M').sum()
#         ## Roll data for chosen period
#         rolling_time_series = filteret_time_series.rolling(months_number, min_periods = months_number).sum()
#         if rolling_time_series.sum() == 0.0:
#             #print('zeros')
#             spi_monthly_values[:,lat_idx,lon_idx] = np.nan
#         elif len(rolling_time_series.dropna()) < (0.5 * len(rolling_time_series)):
#             spi_monthly_values[:,lat_idx,lon_idx] = np.nan
#             #print('missing data')
#         else:
#             rolling_time_series = rolling_time_series.dropna()
#             ## calculate SPI if there is not too much data missing
#             spi_calc = si.spi(rolling_time_series, dist=scs.lognorm)
#             spi_monthly_values[:t,lat_idx,lon_idx] = np.nan
#             spi_monthly_values[t:,lat_idx,lon_idx] = spi_calc
#         x = x+1
#     print('SPI: ' + str(x*100/(len(latitude)*len(longitude))) + '%')
        
# ## Write data to xarray dataarray
# spi_monthly_array = xr.Dataset(
#     {'SPI' : xr.DataArray(
#         data = np.array(spi_monthly_values),
#         dims =('time','latitude', 'longitude',),
#         coords={'latitude': latitude, 'longitude': longitude, 'time': time_spi},
#         attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

# print('SPI values loaded')

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
# spi_monthly_array.to_netcdf('SPI_' + str(months_number) + '_month_1950_2022.nc')
#%%
spi_monthly_array = xr.open_dataset('SPI_' + str(months_number) + '_month_1950_2022.nc')

spi_monthly_array = spi_monthly_array.where(zero_value_mask == 1, drop = True)

#%%
## Plot the spatial distribution of infinite values and convert into nan

count_infinite = np.isinf(spi_monthly_array['SPI']).sum(dim='time')

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(0,80,9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
count_infinite.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.gray_r,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
ax.set_title('SPI ' + str(months_number) + ' month - infinite values count', fontsize=25)
plt.show()

### plot precipitation and SPI fom the location with most infinite values
fig = plt.figure(figsize=(12, 8))
coordinate = np.unravel_index(count_infinite.argmax(), count_infinite.shape)
precip[:,coordinate[0],coordinate[1]].plot()
plt.title('Precipitation time-series for most infinite values', fontsize=25)
plt.xlabel('')
plt.show()

coordinate = np.unravel_index(count_infinite.argmax(), count_infinite.shape)
spi_monthly_array['SPI'][:,coordinate[0],coordinate[1]].plot()
plt.title('SPI time-series for most infinite values', fontsize = 25)
plot_text = 'Infinite values: ' + str(count_infinite.max().values)
plt.text(-8000, 8, plot_text, fontsize=10, color='black')
plt.xlabel('')
plt.ylabel('Standardized Precicpitation Index (SPI)')
plt.show()

### Make infinite values into nan

spi_monthly_array_1 = spi_monthly_array.copy()

# Create a new dataset with NaN values where infinite values exist
spi_monthly_array_1['SPI'] = spi_monthly_array_1['SPI'].where(spi_monthly_array_1['SPI'] != np.inf, drop = False)

#spi_monthly_array_1['SPI'] = spi_monthly_array_1['SPI'].fillna(spi_monthly_array_1['SPI'])

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

periods = [str(' MJJAS 1950-2022'), str(' NDJFM 1950-2022'), str(' All year 1950-2022')]

months_period = [[5,6,7,8,9],[1,2,3,11,12],[1,2,3,4,5,6,7,8,9,10,11,12]]

period = 'EOF NDJFM and SPI MJJAS'

num = 0

# for period in periods:
    
# months = months_period[num]

months_SPI = [5,6,7,8,9]
months_EOF = [1,2,3,11,12]
# months_EOF = [-1,1,2,3,4]
    
# Define number of EOF's wanted
neof = 4
# Read data from existing daatafile
ERA5_data = xr.load_dataset('ERA5_GP500_50_2023.grib', engine = 'cfgrib')
# Extract geopotential from datafile
z_values = ERA5_data['z'][:-9]

#z_values = z_values[60:300,:,:]

if len(months_EOF) < 12:
    z_values = z_values.sel(time=z_values['time.month'].isin(months_EOF))

if len(months_EOF) < 12:
    z_values = z_values.where(z_values['time.month'].isin(months_EOF))

z_values_norm = z_values / z_values.mean(dim = 'latitude').mean(dim = 'longitude')

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(z_values.coords['latitude'].values)).clip(0., 1.)
wgts = np.sqrt((coslat))[..., np.newaxis]
solver = Eof(z_values_norm.values, weights=wgts)

# Retrieve the leading EOF's, expressed as the covariance between the leading 
# PC time series and the input SLP anomalies at each grid point.
pcs = solver.pcs(npcs = neof)
eofs = solver.eofs(neofs = neof, eofscaling=0)

varfrac = solver.varianceFraction(neof)
errors = solver.northTest(neigs=neof, vfscaled=True)
eigenvalues = solver.eigenvalues(neigs=neof)



eof_res_array = xr.Dataset(
    {'modes' : xr.DataArray(
        data = np.array(eofs),
        dims =('EOF','latitude', 'longitude',),
        coords={'latitude': z_values['latitude'], 'longitude': z_values['longitude'], 'EOF': np.arange(1, neof + 1)},
        attrs={'description': 'EOF at 500 hPa lavel', 'units': 'Relative scale'})},
    {'PC_modes' : xr.DataArray(
        data = np.array(pcs),
        dims =('time', 'pcs'),
        coords={'time': z_values['time'],'pcs': np.arange(1, neof + 1)},
        attrs={'description': 'PCS at 500 hPa lavel', 'units': 'Relative scale'})})

# eof_res_array = eof_res_array.reindex(time=ERA5_data['time'])
# z_values_norm = z_values_norm.reindex(time=ERA5_data['time'])

for i in eof_res_array['EOF']:
    ## Find the highest value for the EOF-mode
    max_EOF = eof_res_array['modes'][i-1,:,:].where(eof_res_array['modes'][i-1,:,:] == eof_res_array['modes'][i-1,:,:].max(), drop = True)
    ## Find the pressure values for the high PC-values for the mode
    z_values_high_PC = z_values_norm.where(eof_res_array['PC_modes'][:,i-1] > 0.65, np.nan)
    ## Take the mean of the preassure for high PC-values
    z_values_high_PC_mean = z_values_high_PC.mean(dim = 'time')
    ## Find the mean pressure for the location with highest EOF-anomaly
    relative_pressure_value = z_values_high_PC_mean.sel(longitude = max_EOF['longitude'].values, latitude = max_EOF['latitude'].values).values
    ## Calculate the difference between the pressure value and the mean value for this mode
    test_value = relative_pressure_value - z_values_high_PC_mean.mean().values
    if test_value < 0:
        ## inverse the centers of action, if test value is negative
        eof_res_array['modes'][i-1] = eof_res_array['modes'][i-1] * (-1)
        pcs[:,i-1] = pcs[:,i-1] * (-1)
        eof_res_array['PC_modes'][:,i-1] = eof_res_array['PC_modes'][:,i-1] * (-1)
    
eof_res_array['modes'][3] = eof_res_array['modes'][3] * (-1)
pcs[:,3] = pcs[:,3] * (-1)
eof_res_array['PC_modes'][:,3] = eof_res_array['PC_modes'][:,3] * (-1)

# if months == months_period[1]:
    
# elif months == months_period[0]:
#     eof_res_array['modes'][2] = eof_res_array['modes'][2] * (-1)
#     pcs[:,2] = pcs[:,2] * (-1)
#     eof_res_array['PC_modes'][:,2] = eof_res_array['PC_modes'][:,2] * (-1)
    
#     eof_res_array['modes'][3] = eof_res_array['modes'][3] * (-1)
#     pcs[:,3] = pcs[:,3] * (-1)
#     eof_res_array['PC_modes'][:,3] = eof_res_array['PC_modes'][:,3] * (-1)



eof_res_array_GP = (eof_res_array['modes'] * z_values.mean()) / 9.80665

ERA5_data.close()

print('EOF 500 hPa calculation done')

#%%
## plotting SPI values for year 1976

plotting_data = spi_monthly_array['SPI'][months_SPI[0]+312:months_SPI[0]+len(months_SPI)+312].mean(dim = 'time')

fig = plt.figure(figsize=(12, 8))
clevs = np.linspace(-2, 2, 9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
ax.set_title('Mean SPI values ' + str(months_number) + '-month for MJJAS 1976', fontsize=25)
plt.show()

#%%
## calculate PC for MSLP correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

neof = 4
# Read data from existing daatafile
ERA5_datas = xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# Extract geopotential from datafile
sp_values = ERA5_datas['msl'][:-9]  

#sp_values = sp_values[60:300,:,:]

if len(months_EOF) < 12:
    sp_values = sp_values.sel(time=sp_values['time.month'].isin(months_EOF))

sp_values_norm = sp_values / sp_values.mean(dim = 'latitude').mean(dim = 'longitude')

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(z_values.coords['latitude'].values)).clip(0., 1.)
wgts = np.sqrt((coslat))[..., np.newaxis]
solver = Eof(sp_values_norm.values, weights=wgts)

# Retrieve the leading EOF's, expressed as the covariance between the leading 
# PC time series and the input SLP anomalies at each grid point.
eof_res = solver.eofsAsCovariance(neofs=neof)

eofs_sf = solver.eofs(neofs = neof, eofscaling=0)
pcs_sf = solver.pcs(npcs = neof)

varfrac_sf = solver.varianceFraction(neof)
errors_sf = solver.northTest(neigs=neof, vfscaled=True)
eigenvalues_sf = solver.eigenvalues(neigs=neof)

eof_res_array_sf = xr.Dataset(
    {'modes' : xr.DataArray(
        data = np.array(eofs_sf),
        dims =('EOF','latitude', 'longitude',),
        coords={'latitude': sp_values['latitude'], 'longitude': sp_values['longitude'], 'EOF': np.arange(1, neof + 1)},
        attrs={'description': 'EOF at 500 hPa lavel', 'units': 'Relative scale'})},
    {'PC_modes' : xr.DataArray(
        data = np.array(pcs_sf),
        dims =('time', 'pcs'),
        coords={'time': sp_values['time'],'pcs': np.arange(1, neof + 1)},
        attrs={'description': 'PCS at MSLP', 'units': 'Relative scale'})})

for i in eof_res_array_sf['EOF']:
    ## Find the highest value for the EOF-mode
    max_EOF_sf = eof_res_array_sf['modes'][i-1,:,:].where(eof_res_array_sf['modes'][i-1,:,:] == 
                                                        eof_res_array_sf['modes'][i-1,:,:].max(), 
                                                        drop = True)
    ## Find the pressure values for the high PC-values for the mode
    sp_values_high_PC = sp_values_norm.where(eof_res_array_sf['PC_modes'][:,i-1] > 0.65, np.nan)
    ## Take the mean of the preassure for high PC-values
    sp_values_high_PC_mean = sp_values_high_PC.mean(dim = 'time')
    ## Find the mean pressure for the location with highest EOF-anomaly
    relative_pressure_value_sf = sp_values_high_PC_mean.sel(longitude = max_EOF_sf['longitude'].values, 
                                                   latitude = max_EOF_sf['latitude'].values).values
    ## Calculate the difference between the pressure value and the mean value for this mode
    test_value_sf = relative_pressure_value_sf - sp_values_high_PC_mean.mean().values
    if test_value_sf < 0:
        ## inverse the centers of action, if test value is negative
        eof_res_array_sf['modes'][i-1] = eof_res_array_sf['modes'][i-1] * (-1)
        pcs_sf[:,i-1] = pcs_sf[:,i-1] * (-1)
        eof_res_array_sf['PC_modes'][:,i-1] = eof_res_array_sf['PC_modes'][:,i-1] * (-1)

eof_res_array_sf['modes'][0] = eof_res_array_sf['modes'][0] * (-1)
pcs_sf[:,0] = pcs_sf[:,0] * (-1)
eof_res_array_sf['PC_modes'][:,0] = eof_res_array_sf['PC_modes'][:,0] * (-1)

# eof_res_array_sf['modes'][1] = eof_res_array_sf['modes'][1] * (-1)
# pcs_sf[:,1] = pcs_sf[:,1] * (-1)
# eof_res_array_sf['PC_modes'][:,1] = eof_res_array_sf['PC_modes'][:,1] * (-1)

eof_res_array_sf['modes'][2] = eof_res_array_sf['modes'][2] * (-1)
pcs_sf[:,2] = pcs_sf[:,2] * (-1)
eof_res_array_sf['PC_modes'][:,2] = eof_res_array_sf['PC_modes'][:,2] * (-1)

eof_res_array_sf['modes'][3] = eof_res_array_sf['modes'][3] * (-1)
pcs_sf[:,3] = pcs_sf[:,3] * (-1)
eof_res_array_sf['PC_modes'][:,3] = eof_res_array_sf['PC_modes'][:,3] * (-1)  
    
# elif months == months_period[1]:
#     eof_res_array_sf['modes'][3] = eof_res_array_sf['modes'][3] * (-1)
#     pcs_sf[:,3] = pcs_sf[:,3] * (-1)
#     eof_res_array_sf['PC_modes'][:,3] = eof_res_array_sf['PC_modes'][:,3] * (-1)  
    
eof_res_array_sf_hPa = (eof_res_array_sf['modes'] * sp_values.mean()) / 100

ERA5_datas.close()

print('EOF MSLP calculation done')
#%%
## plot the EOF results for MSLP and Z500
### MAKE UNITS IN HPA!!!!
## For Z500
for i in eof_res_array['EOF']:
    fig = plt.figure(figsize=(11.5, 6.8))
    clevs = np.linspace(-0.01, 0.01, 21)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    eof_res_array['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('EOF' + str(eof_res_array['EOF'][i-1].values) + ' GP500 anomaly' + period, fontsize=25)
    plt.text(-5000000, -3900000, plot_text, fontsize=25, color='black')
    plt.show()

## for MSLP
for i in eof_res_array_sf['EOF']:
    fig = plt.figure(figsize=(11.5, 6.8))
    clevs = np.linspace(-0.01, 0.01, 21)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    eof_res_array_sf['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('EOF' + str(eof_res_array_sf['EOF'][i-1].values) + ' MSLP anomaly' + period, fontsize=25)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    plt.text(-5000000, -3900000, plot_text, fontsize=25, color='black')
    plt.show()

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=-20, 
                                                                    central_latitude=60)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-0.01, 0.01, 21)
    ax.coastlines()
    
    eof_res_array['modes'][i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                              transform=ccrs.PlateCarree(), add_colorbar=False)
    plot_text = str(round(varfrac[i]*100, 1)) + '%'
    subtitle = 'EOF' + str(i+1) + ' GP500 ' + plot_text
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-0.01, 0.01, 21)
    ax.coastlines()
    
    eof_res_array_sf['modes'][i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)
    plot_text = str(round(varfrac_sf[i]*100, 1)) + '%'
    ax.set_title('EOF' + str(i+ 1) + ' MSLP ' + plot_text, fontsize=20)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.89])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=eof_res_array['modes'][0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('') 
plt.tight_layout()
plot_text = 'Explained variance - GP500: ' + str(round(varfrac.sum()*100,1)) + '% - MSLP: ' + str(round(varfrac_sf.sum()*100,1)) + '%'
#plt.text(-36,0.0095, plot_text, fontsize=25, color='black')
plt.suptitle('EOF modes of variability in North Atlantic' + period, y=1.12, ha='center', fontsize = 45)
plt.title(plot_text, x=-25, y=1.05, fontsize=28, color ='gray')

plt.show()


print('EOF plots')

# num_subplots = 4

# # Create a 2x4 grid of subplots
# fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8), 
#                         subplot_kw={'projection': ccrs.Orthographic(central_longitude=-20, 
#                                                                     central_latitude=60)})

# # Flatten the 2D array of subplots into a 1D array
# axs = axs.flatten()

# # Plotting correlations for GP500
# for i, ax in enumerate(axs[:num_subplots]):
#     clevs = np.linspace(-500, 500, 21)
#     ax.coastlines()
   
#     eof_res_array_GP[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                               transform=ccrs.PlateCarree(), add_colorbar=False)
#     plot_text = str(round(varfrac[i]*100, 1)) + '%'
#     subtitle = 'PC' + str(i+1) + ' GP500 ' + plot_text
#     ax.set_title(subtitle, fontsize=20)
    
#     #plt.show()
# # Plotting correlations for MSLP
# for i, ax in enumerate(axs[num_subplots:]):
#     clevs = np.linspace(-10, 10, 21)
#     ax.coastlines()

#     eof_res_array_sf_hPa[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                                 transform=ccrs.PlateCarree(), add_colorbar=False)
#     plot_text = str(round(varfrac_sf[i]*100, 1)) + '%'
#     ax.set_title('EOF' + str(i+ 1) + ' MSLP ' + plot_text, fontsize=20)
# cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
# colorbar1 = Colorbar(ax=cax1, mappable=eof_res_array['modes'][0].plot.contourf(levels = clevs, 
#                                                                 cmap=plt.cm.RdBu_r, 
#                                                                 add_colorbar=False), orientation='vertical')
# cax1.xaxis.label.set_text('')
# cax1.title.set_text('')

# # Adjust layout for better spacing
# plt.tight_layout()

# plt.title('EOF modes of variability in North Atlantic' + period, fontsize = 48, loc='right', pad=33)

# # Show the plot
# plt.show()

# print('New EOF plots')


#%%
## Calculate correlation between SPI and PC's

lats = spi_monthly_array_1['latitude']
longs = spi_monthly_array_1['longitude']

### Calculating the correaltion for the 500 hPa level

eof_res_array1 = eof_res_array.reindex(time=ERA5_data['time'][months_EOF[0]-1:-9])

selected_SPI_data = spi_monthly_array_1

if months_number > 1:
    pcs_summed = np.cumsum(np.nan_to_num(eof_res_array1['PC_modes'].values), axis = 0, dtype=float)
    pcs_summed_roll = pcs_summed.copy()
    pcs_summed_roll[months_number:] = pcs_summed[:-months_number]
    pcs_summed_roll[:months_number] = 0
    div_num = np.transpose(np.tile(np.cumsum(pcs_summed[:,0] != 0) - np.cumsum(pcs_summed_roll[:,0] != 0),(4,1)))
    pcs_roll_mean = (pcs_summed - pcs_summed_roll) / div_num
    nan_mask = np.isnan(eof_res_array1['PC_modes'].values)
    # pcs = pcs_roll_mean
    pcs_roll_mean[nan_mask] = np.nan
    pd_pcs_roll = pd.DataFrame(pcs_roll_mean)
    pd_pcs_roll = pd_pcs_roll.dropna()
    pcs = pd_pcs_roll.values

if len(months_SPI) < 12:
    selected_SPI_data = spi_monthly_array_1.sel(time=spi_monthly_array_1['time.month'].isin(months_SPI))
    

## Find correlations on of PC and the SPI-values
if len(pcs) <= len(selected_SPI_data['time']):
    corr = correlation_map(pcs[t:len(pcs)],selected_SPI_data['SPI'][t:len(pcs)].values)
else:
    corr = correlation_map(pcs[t:len(selected_SPI_data['time'])],selected_SPI_data['SPI'][t:len(selected_SPI_data['time'])].values)
    
## correlation data as Xarray dataarray
corr_xr = xr.DataArray(
    corr,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

corr_xr = corr_xr.dropna(dim='latitude', how = 'all')
corr_xr = corr_xr.dropna(dim='longitude', how = 'all')

## plotting the correlations
for i in corr_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()

### Calculate the correlations for the MSLP

eof_res_array_sf_1 = eof_res_array_sf.reindex(time=ERA5_data['time'][months_EOF[0]-1:-9])

if months_number > 1:
    pcs_summed_sf = np.cumsum(np.nan_to_num(eof_res_array_sf_1['PC_modes'].values), axis = 0, dtype=float)
    pcs_summed_roll_sf = pcs_summed_sf.copy()
    pcs_summed_roll_sf[months_number:] = pcs_summed_sf[:-months_number]
    pcs_summed_roll_sf[:months_number] = 0
    div_num_sf = np.transpose(np.tile(np.cumsum(pcs_summed_sf[:,0] != 0) - np.cumsum(pcs_summed_roll_sf[:,0] != 0),(4,1)))
    pcs_roll_mean_sf = (pcs_summed_sf - pcs_summed_roll_sf) / div_num_sf
    nan_mask_sf = np.isnan(eof_res_array_sf_1['PC_modes'].values)
    # pcs = pcs_roll_mean
    pcs_roll_mean_sf[nan_mask_sf] = np.nan
    pd_pcs_roll_sf = pd.DataFrame(pcs_roll_mean_sf)
    pd_pcs_roll_sf = pd_pcs_roll_sf.dropna()
    pcs_sf = pd_pcs_roll_sf.values

## Find correlations on of PC and the SPI-values
if len(pcs) <= len(selected_SPI_data['time']):
    corr_sf = correlation_map(pcs_sf[t:len(pcs)],selected_SPI_data['SPI'][t:len(pcs)].values)
else:
    corr_sf = correlation_map(pcs_sf[t:len(selected_SPI_data['time'])],selected_SPI_data['SPI'][t:len(selected_SPI_data['time'])].values)
    

## correlation data as Xarray dataarray
corr_sf_xr = xr.DataArray(
    corr_sf,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

corr_sf_xr = corr_sf_xr.dropna(dim='latitude', how = 'all')
corr_sf_xr = corr_sf_xr.dropna(dim='longitude', how = 'all')


print('500 hPa correlations calculation done')

## plotting the
for i in corr_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_sf_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()

print('surface correlations calculation done')

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8.9), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    corr_xr[i,:,:].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.06, 0.02, 0.88])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and ' + str(months_number) + '-months SPI-value correlation coefficient' + period, fontsize = 41, loc='right', pad=33)

# Show the plot
plt.show()

lats = corr_xr['latitude']
longs = corr_xr['longitude']

selected_SPI_data = selected_SPI_data.dropna(dim='latitude', how = 'all')
selected_SPI_data = selected_SPI_data.dropna(dim='longitude', how = 'all')

#%%

## Limit data to only high correlations

## high correlation level
H = 0.3

high_corr = ma.masked_inside(corr_xr, -H, H)

high_corr_xr = xr.DataArray(
    high_corr,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

for i in high_corr_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI ' + str(months_number) + '-month 50-22 Hcorr', fontsize=25)
    plt.show()

high_corr_sf = ma.masked_inside(corr_sf_xr, -H, H)

high_corr_sf_xr = xr.DataArray(
    high_corr_sf,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

for i in corr_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 Hcorr', fontsize=25)
    plt.show()

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8.9), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
cax1 = fig.add_axes([1, 0.06, 0.02, 0.88])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and ' + str(months_number) + '-months SPI-value high correlation coefficients' + period, fontsize = 38, loc='right', pad=33)

# Show the plot
plt.show()

print('High correlations plottet')
#%%

# Calculate the regression coefficient for the correlations

exp_coef = pcs[t:len(selected_SPI_data['time']),:]/pcs[t:len(selected_SPI_data['time']),:].std()

regressions = np.empty((neof, len(lats), len(longs)))

x = 0 
for i in range(neof):
    for lat in range(len(lats)):
        for lon in range(len(longs)):
            calc_data = selected_SPI_data['SPI'][t:,lat,lon]
            if np.isnan(calc_data).any() ==  True:
                regressions[i, lat, lon] = np.nan
            else:
                polifit = poly.polyfit(exp_coef[:,i],calc_data,1)
                regressions[i, lat, lon] = polifit[1]
        print('summer 500 hPa regression ' + str(round((x*100/(len(lats)*4)),1)) + ' %')
        x = x+1
            
regressions_xr = xr.DataArray(
    regressions,
    coords={'PC':np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

print('500 hPa regressions calculated')

regressions_sf = np.empty((neof, len(lats), len(longs)))

exp_coef_sf = pcs_sf[t:len(selected_SPI_data['time']),:]/pcs_sf[t:len(selected_SPI_data['time']),:].std()

x = 0 
for i in range(neof):
    for lat in range(len(lats)):
        for lon in range(len(longs)):
            calc_data = selected_SPI_data['SPI'][t:,lat,lon]
            if np.isnan(calc_data).any() ==  True:
                regressions_sf[i, lat, lon] = np.nan
            else:
                polifit = poly.polyfit(exp_coef_sf[:,i],calc_data,1)
                regressions_sf[i, lat, lon] = polifit[1]
        print('summer surface regression ' + str(round((x*100/(len(lats)*4)),1)) + '%')
        x = x+1
            
regressions_sf_xr = xr.DataArray(
    regressions_sf,
    coords={'PC':np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

print('MSLP regressions calculated')

#%%

## plot the regression coefficients calculated above

for i in regressions_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' Z500 and SPI ' + str(months_number) + '-month 50-22 RegSlope', fontsize=25)
    plt.show()

print('look at those plots!')

for i in regressions_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 RegSlope', fontsize=25)
    plt.show()

print('More beautiful plots are made!')

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8.9), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    regressions_xr[i,:,:].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    regressions_sf_xr[i,:,:].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
cax1 = fig.add_axes([1, 0.06, 0.02, 0.88])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and ' + str(months_number) + '-months SPI-value regression slope' + period, fontsize = 45, loc='right', pad=33)

# Show the plot
plt.show()
print('A combined plot was also made')

#%%
### plotting figures og correlations with high correlationsions hatched

for i in corr_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform= ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform = ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False, alpha=0)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()


for i in corr_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,
                                    alpha=0)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()
print('Hatched in stead of masked high correlations')


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8.9), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})
# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i
                             ].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,
                                    alpha=0)
   
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    ax.coastlines()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], add_colorbar=False,
                                    alpha=0)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
cax1 = fig.add_axes([1, 0.06, 0.02, 0.88])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and ' + str(months_number) + '-months SPI-value correlation coefficients' + period, fontsize = 41, loc='right', pad=33)

plt.show()

#%%
## Plotting data with contour line for EOF-modes
for i in corr_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                            transform= ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu,
                                    transform = ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
    
    clevs_eof = np.linspace(-0.010, 0.010, 21)
    contour_plot = eof_res_array['modes'][i].plot.contour(ax=ax, levels=clevs_eof,
                             transform=ccrs.PlateCarree(), add_colorbar=False,
                             linewidths = 0.7, linestyles = 
                             ['dashed' if val < 0 else 'dotted' if val > 0 
                              else 'solid' for val in clevs_eof], 
                             colors ='k', alpha = 1)
    ax.clabel(contour_plot, fontsize = 7)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()


print('Combined plot with hathcing')
for i in corr_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                            transform=ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
    
    clevs_eof = np.linspace(-0.010, 0.010, 21)
    contour_plot = eof_res_array_sf['modes'][i].plot.contour(ax=ax, levels=clevs_eof,
                             transform=ccrs.PlateCarree(), add_colorbar=False,
                             linewidths = 0.7, linestyles = 
                             ['dashed' if val < 0 else 'dotted' if val > 0 
                              else 'solid' for val in clevs_eof], 
                             colors ='k', alpha =1)
    ax.clabel(contour_plot, fontsize = 7)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 corr', fontsize=25)
    plt.show()


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(25, 10), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=-20, 
                                                                    central_latitude=60)})
# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
    
    clevs_eof = np.linspace(-50, 50, 21)
    contour_plot = eof_res_array_GP[i].plot.contour(ax=ax, levels=clevs_eof,
                             transform=ccrs.PlateCarree(), add_colorbar=False,
                             linewidths = 0.7, linestyles = 
                             ['dashed' if val < 0 else 'dotted' if val > 0 
                              else 'solid' for val in clevs_eof], 
                             colors ='k', alpha =1)
    ax.clabel(contour_plot, fontsize = 8)
    plot_text = str(round(varfrac[i]*100, 1)) + '%'
    subtitle = 'PC' + str(i+1) + ' GP500 - ' + plot_text
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    ax.coastlines()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)
    
    clevs_eof = np.linspace(-10, 10, 21)
    contour_plot = eof_res_array_sf_hPa[i].plot.contour(ax=ax, levels=clevs_eof,
                             transform=ccrs.PlateCarree(), add_colorbar=False,
                             linewidths = 0.7, linestyles = 
                             ['dashed' if val < 0 else 'dotted' if val > 0 
                              else 'solid' for val in clevs_eof], 
                             colors ='k', alpha =1)
    ax.clabel(contour_plot, fontsize = 8)
    
    plot_text = str(round(varfrac_sf[i]*100, 1)) + '%'
    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP - ' + plot_text, fontsize=20)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs for NDJFM and SPI-12 MJJAS correlation coefficients and EOF analysis', fontsize = 52, loc='right', pad=33)

plt.show()    


#%%

for i in regressions_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    regressions_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI ' + str(months_number) + '-month 50-22 RegSlope', fontsize=25)
    plt.show()


for i in regressions_sf_xr['PC']:
    fig = plt.figure(figsize=(12, 8))
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    regressions_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI ' + str(months_number) + '-month 50-22 RegSlope', fontsize=25)
    plt.show()


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 8.9), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    regressions_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)
   
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    ax.coastlines()
    regressions_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_sf_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False,alpha=0)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
cax1 = fig.add_axes([1, 0.06, 0.02, 0.88])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')
cax1.xaxis.label.set_text('')
cax1.title.set_text('')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and ' + str(months_number) + '-months SPI-value regression slopes' + period, fontsize = 45, loc='right', pad=33)

plt.show()

print('Regressions are now also available with hatching')

#%%

corr_xr_edit = corr_xr
corr_xr_edit['latitude'] = corr_xr_edit['latitude'] + 0.125
corr_xr_edit['longitude'] = corr_xr_edit['longitude'] + 0.125

regressions_xr_edit = regressions_xr
regressions_xr_edit['latitude'] = regressions_xr_edit['latitude'] + 0.125
regressions_xr_edit['longitude'] = regressions_xr_edit['longitude'] + 0.125

eof_res_array['corr'] = corr_xr_edit
eof_res_array['reg'] = regressions_xr_edit
eof_res_array.to_netcdf('500hPa_correlation_analysis_' + str(months_number) + '_month_'+ period + '.nc')

corr_sf_xr_edit = corr_sf_xr
corr_sf_xr_edit['latitude'] = corr_sf_xr_edit['latitude'] + 0.125
corr_sf_xr_edit['longitude'] = corr_sf_xr_edit['longitude'] + 0.125

regressions_sf_xr_edit = regressions_sf_xr
regressions_sf_xr_edit['latitude'] = regressions_sf_xr_edit['latitude'] + 0.125
regressions_sf_xr_edit['longitude'] = regressions_sf_xr_edit['longitude'] + 0.125

eof_res_array_sf['corr'] = corr_sf_xr_edit
eof_res_array_sf['reg'] = regressions_sf_xr_edit
eof_res_array_sf.to_netcdf('MSLP_correlation_analysis_' + str(months_number) + '_month_'+ period + '.nc')

num = num + 1    

print('This calculation is done!')

#%%
end_time = time.time()
elapsed_time = end_time - start_time
elapsed_time_hour = elapsed_time/3600
print('calculated in ' + str(round(elapsed_time_hour,2)) + ' hours - ' + str(round(elapsed_time)) + ' seconds')