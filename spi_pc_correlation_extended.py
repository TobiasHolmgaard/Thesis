# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 22:20:28 2023

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

start_time = time.time()
#%%
## Calculate PC for Z500 correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

months = [5,6,7,8,9]

period = str('MJJAS 1950-2022')

# Define number of EOF's wanted
neof = 4
# Read data from existing daatafile
ERA5_data = xr.load_dataset('ERA5_GP500_50_2023.grib', engine = 'cfgrib')
# Extract geopotential from datafile
z_values = ERA5_data['z']

#z_values = z_values[60:300,:,:]

if len(months) < 12:
    z_values = z_values.sel(time=z_values['time.month'].isin(months))

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

ERA5_data.close()

print('EOF 500 hPa calculation done')
#%%
## calculate PC for MSLP correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

neof = 4
# Read data from existing daatafile
ERA5_datas = xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# Extract geopotential from datafile
sp_values = ERA5_datas['msl']

#sp_values = sp_values[60:300,:,:]

if len(months) < 12:
    sp_values = sp_values.sel(time=sp_values['time.month'].isin(months))

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

ERA5_datas.close()

print('EOF MSLP calculation done')
#%%
## plot the EOF results for MSLP and Z500

## For Z500
for i in eof_res_array['EOF']:
    clevs = np.linspace(-0.009, 0.009, 19)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    eof_res_array['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('EOF' + str(eof_res_array['EOF'][i-1].values) + ' GP500 anomaly MJJAS 1955-1974', fontsize=14)
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

## for MSLP
for i in eof_res_array_sf['EOF']:
    clevs = np.linspace(-0.009, 0.009, 19)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    eof_res_array_sf['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('EOF' + str(eof_res_array_sf['EOF'][i-1].values) + ' MSLP anomaly MJJAS 1955-1974', fontsize=14)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
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
    clevs = np.linspace(-0.009, 0.009, 21)
    ax.coastlines()
   
    eof_res_array['modes'][i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                              transform=ccrs.PlateCarree(), add_colorbar=False)
    plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    subtitle = 'PC' + str(i+1) + ' GP500' + plot_text
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-0.009, 0.009, 21)
    ax.coastlines()

    eof_res_array_sf['modes'][i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    ax.set_title('EOF' + str(i+ 1) + ' MSLP' + plot_text, fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=eof_res_array['modes'][0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('EOF modes of variability in North Atlantic all year 1950-2022', fontsize = 46, loc='right', pad=33)

# Show the plot
plt.show()


print('EOF plots')
#%%
## Calculate SPI-values - import data - dataset is large! - Change dataset!!!
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

clevs = np.linspace(0, 120, 13)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_precip[15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('standard deviation (mm)')
ax.set_title('Std for precipitation in gridcells 1950-2022', fontsize=12)
plt.show()

clevs = np.linspace(0, 2, 9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
std_norm_precip[15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.Greys,
                         transform=ccrs.PlateCarree(), add_colorbar=False)
colorbar = plt.colorbar(mappable=std_norm_precip.plot.contourf(ax=ax, levels = clevs, 
                                                          cmap=plt.cm.Greys, add_colorbar=False), 
                        ax=ax, orientation='vertical', pad=0.02)
colorbar.set_label('Std / mean (mm)')
ax.set_title('Std normalized precip in gridcells 1950-2022', fontsize=12)
plt.show()

## define perdiod for calculating values
months_number = 3
t = months_number - 1

## resample time variable for monthly values
#time_months = time.resample(time="1M").max(dim="time")

## create array for resulting values
spi_monthly_values = np.empty((len(time_spi), len(latitude), len(longitude)))

# x = 0

# print('time for SPI calculation')

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
#         else:
#             rolling_time_series = rolling_time_series.dropna()
#             ## calculate SPI if there is not too much data missing
#             spi_calc = si.spi(rolling_time_series, dist=scs.lognorm)
#             spi_monthly_values[:t,lat_idx,lon_idx] = np.nan
#             spi_monthly_values[t:,lat_idx,lon_idx] = spi_calc
#         x = x+1
#     print('SPI: ' + str(x*100/(len(latitude)*len(longitude))) + '%')
        # 
# ## Write data to xarray dataarray
# spi_monthly_array = xr.Dataset(
#     {'SPI' : xr.DataArray(
#         data = np.array(spi_monthly_values),
#         dims =('time','latitude', 'longitude',),
#         coords={'latitude': latitude, 'longitude': longitude, 'time': time_spi},
#         attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})

# print('SPI values loaded')

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
# spi_monthly_array.to_netcdf('SPI_3_month_1950_2022.nc')
#%%
spi_monthly_array = xr.open_dataset('SPI_3_month_1950_2022.nc')

#%%
## plotting SPI values for the first year


plotting_data = spi_monthly_array['SPI'][months[0]+312:months[0]+len(months)+312].mean(dim = 'time')

clevs = np.linspace(-2, 2, 9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)

ax.coastlines()

plotting_data[15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
#plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
ax.set_title('Mean SPI values 3-month for All year 1976', fontsize=14)
#plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
plt.show()

#%%
## Plot the spatial distribution of infinite values and convert into nan

count_infinite = np.isinf(spi_monthly_array['SPI']).sum(dim='time')

clevs = np.linspace(0,80,9)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)
ax.coastlines()
count_infinite.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.gray_r,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
ax.set_title('SPI 3 month - infinite values count', fontsize=16)
plt.show()

### plot precipitation and SPI fom the location with most infinite values
coordinate = np.unravel_index(count_infinite.argmax(), count_infinite.shape)
precip[:,coordinate[0],coordinate[1]].plot()
plt.title('Precipitation time-series for most infinite values', fontsize=14)
plt.xlabel('')
plt.show()

coordinate = np.unravel_index(count_infinite.argmax(), count_infinite.shape)
spi_monthly_array['SPI'][:,coordinate[0],coordinate[1]].plot()
plt.title('SPI time-series for most infinite values', fontsize = 14)
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


#%%
## Calculate correlation between SPI and PC's

lats = spi_monthly_array_1['latitude']
longs = spi_monthly_array_1['longitude']

## select month for analysis - summer/winter/full year
selected_SPI_data = spi_monthly_array_1.sel(time=spi_monthly_array_1['time.month'].isin(months))

### Calculating the correaltion for the 500 hPa level

if months_number > 1:
    pcs_summed = np.cumsum(pcs, axis = 0, dtype=float)
    pcs_summed_roll = pcs_summed.copy()
    pcs_summed_roll[months_number:] = pcs_summed[:-months_number]
    pcs_roll_mean = (pcs_summed - pcs_summed_roll) / months_number
    pcs = pcs_roll_mean

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

## plotting the correlations
for i in corr_xr['PC']:
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 3-month 50-22 corr', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

### Calculate the correlations for the MSLP

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

print('500 hPa correlations calculation done')

## plotting the
for i in corr_sf_xr['PC']:
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_sf_xr['PC'][i].values+1) + ' MSLP and SPI 3-month 50-22 corr', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('surface correlations calculation done')

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    corr_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    corr_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and 3-months SPI-value correlation coefficient all year 1950-2022', fontsize = 41, loc='right', pad=33)

# Show the plot
plt.show()

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
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 3-month 50-22 Hcorr', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

high_corr_sf = ma.masked_inside(corr_sf_xr, -H, H)

high_corr_sf_xr = xr.DataArray(
    high_corr_sf,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

for i in corr_sf_xr['PC']:
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 3-month 50-22 Hcorr', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    high_corr_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    high_corr_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and 3-months SPI-value high correlation coefficients all year 1950-2022', fontsize = 38, loc='right', pad=33)

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
        print('500 hPa ' + str(round((x*100/(len(lats)*4)),1)) + ' %')
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
        print('surface ' + str(round((x*100/(len(lats)*4)),1)) + '%')
        x = x+1
            
regressions_sf_xr = xr.DataArray(
    regressions_sf,
    coords={'PC':np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

print('MSLP regressions calculated')

#%%

## plot the regression coefficients calculated above

#scale = round(max(abs(regressions_xr.min()), regressions_xr.max(),abs(regressions_sf_xr.min()),regressions_sf_xr.max()),1)

for i in regressions_xr['PC']:
    clevs = np.linspace(-1, 1, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' Z500 and SPI 3-month 50-22 RegSlope', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('look at those plots!')

for i in regressions_sf_xr['PC']:
    clevs = np.linspace(-1, 1, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 3-month 50-22 RegSlope', fontsize=14)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('More beautiful plots are made!')

num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
   
    regressions_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=False)
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()

    regressions_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                                transform=ccrs.PlateCarree(), add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and 3-months SPI-value regression slope all year 1950-2022', fontsize = 44, loc='right', pad=33)

# Show the plot
plt.show()
print('A combined plot was also made')

#%%
### plotting figures og correlations with high correlationsions hatched

for i in corr_xr['PC']:
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=True)
    mask = high_corr_xr[i].notnull()
    hatches = np.full(corr_xr[i].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0.9
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False, linestyles = (['-']), linewidth = 0.0001)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 3-month 50-22 corr', fontsize=13)
    plt.show()


for i in corr_sf_xr['PC']:
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
    alpha = 0.5
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 3-month 50-22 corr', fontsize=13)
    plt.show()

print('Hatched in stead of masked high correlations')


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    corr_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_xr[i,15:185,105:345].notnull()
    hatches = np.full(corr_xr[i,15:185,105:345].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0.9
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i,15:185,105:345].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False, linestyles = (['-']), linewidth = 0.0001)
   
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    ax.coastlines()
    corr_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_sf_xr[i,15:185,105:345].notnull()
    hatches = np.full(corr_xr[i,15:185,105:345].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0.5
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i,15:185,105:345].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and 3-months SPI-value correlation coefficients all year 1950-2022', fontsize = 40, loc='right', pad=33)

print('Combined plot with hathcing')

#%%

for i in regressions_xr['PC']:
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
    alpha = 0.5
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 3-month 50-22 RegSlope', fontsize=13)
    plt.show()


for i in regressions_sf_xr['PC']:
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
    alpha = 0.5
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 3-month 50-22 RegSlope', fontsize=13)
    plt.show()


num_subplots = 4

# Create a 2x4 grid of subplots
fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
                        subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
                                                                    central_latitude=50.375)})

# Flatten the 2D array of subplots into a 1D array
axs = axs.flatten()

# Plotting correlations for GP500
for i, ax in enumerate(axs[:num_subplots]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    regressions_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_xr[i,15:185,105:345].notnull()
    hatches = np.full(corr_xr[i,15:185,105:345].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0.9
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_xr[i,15:185,105:345].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    add_colorbar=False, linestyles = (['-']), linewidth = 0.0001)
   
    subtitle = 'PC' + str(i+1) + ' GP500'
    ax.set_title(subtitle, fontsize=20)
    
    #plt.show()
# Plotting correlations for MSLP
for i, ax in enumerate(axs[num_subplots:]):
    clevs = np.linspace(-1, 1.00, 21)
    ax.coastlines()
    ax.coastlines()
    regressions_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                            transform=ccrs.PlateCarree(), add_colorbar=False)
    mask = high_corr_sf_xr[i,15:185,105:345].notnull()
    hatches = np.full(corr_xr[i,15:185,105:345].shape, '', dtype='U1')
    hatching_pattern = '....'
    # Set hatching pattern in the masked region
    hatches[mask] = hatching_pattern
    alpha = 0.5
    
    hatches = np.where(mask, hatching_pattern, hatches)
    high_corr_sf_xr[i,15:185,105:345].plot.contourf(hatch = hatching_pattern, ax=ax, levels=clevs, 
                                    cmap=plt.cm.RdBu_r,
                                    transform=ccrs.PlateCarree(), 
                                    hatches = [hatching_pattern], 
                                    alpha = alpha, add_colorbar=False)

    ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
#plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
                                                                cmap=plt.cm.RdBu_r, 
                                                                add_colorbar=False), orientation='vertical')

# Adjust layout for better spacing
plt.tight_layout()

plt.title('PCs and 3-months SPI-value regression slopes all year 1950-2022', fontsize = 41, loc='right', pad=33)


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
eof_res_array.to_netcdf('500hPa_correlation_analysis_3_month.nc')

corr_sf_xr_edit = corr_sf_xr
corr_sf_xr_edit['latitude'] = corr_sf_xr_edit['latitude'] + 0.125
corr_sf_xr_edit['longitude'] = corr_sf_xr_edit['longitude'] + 0.125

regressions_sf_xr_edit = regressions_sf_xr
regressions_sf_xr_edit['latitude'] = regressions_sf_xr_edit['latitude'] + 0.125
regressions_sf_xr_edit['longitude'] = regressions_sf_xr_edit['longitude'] + 0.125

eof_res_array_sf['corr'] = corr_sf_xr_edit
eof_res_array_sf['reg'] = regressions_sf_xr_edit
eof_res_array_sf.to_netcdf('MSLP_correlation_analysis_3_month.nc')

#%%
end_time = time.time()
elapsed_time = end_time - start_time
elapsed_time_hour = elapsed_time/3600
print('calculated in ' + str(round(elapsed_time_hour,2)) + ' hours - ' + str(round(elapsed_time)) + ' seconds')

#%%

# num_subplots = 4

# # Create a 2x4 grid of subplots
# fig, axs = plt.subplots(2, num_subplots, figsize=(12, 5.8), subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)})

# # Flatten the 2D array of subplots into a 1D array
# axs = axs.flatten()

# # Plotting correlations for GP500
# for i, ax in enumerate(axs[:num_subplots]):
#     clevs = np.linspace(-1, 1.00, 21)
#     ax.coastlines()
    
#     corr_xr[i,:,70:350].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                              transform=ccrs.PlateCarree(), add_colorbar=False)
#     ax.set_title(f'PC{corr_xr["PC"][i].values + 1} GP500 and SPI 3-month 50-22 corr', fontsize=10)
#     #plt.show()
# # Plotting correlations for MSLP
# for i, ax in enumerate(axs[num_subplots:]):
#     clevs = np.linspace(-1, 1.00, 21)
#     ax.coastlines()

#     corr_sf_xr[i,:,70:350].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                                 transform=ccrs.PlateCarree(), add_colorbar=False)

#     ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP and SPI 3-month 50-22 corr', fontsize=10)
# #plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
# cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
# colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, cmap=plt.cm.RdBu_r, add_colorbar=False), orientation='vertical')

# # Adjust layout for better spacing
# plt.tight_layout()

# plt.title('PCs and 3-months SPI values all year 1950-2022', fontsize = 30, loc='right', pad=30)

# # Show the plot
# plt.show()


# #%%

# num_subplots = 4

# # Create a 2x4 grid of subplots
# fig, axs = plt.subplots(2, num_subplots, figsize=(20, 9.8), 
#                         subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
#                                                                     central_latitude=50.375)})

# # Flatten the 2D array of subplots into a 1D array
# axs = axs.flatten()

# # Plotting correlations for GP500
# for i, ax in enumerate(axs[:num_subplots]):
#     clevs = np.linspace(-1, 1.00, 21)
#     ax.coastlines()
   
#     regressions_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                              transform=ccrs.PlateCarree(), add_colorbar=False)
#     subtitle = 'PC' + str(i+1) + ' GP500'
#     ax.set_title(subtitle, fontsize=20)
    
#     #plt.show()
# # Plotting correlations for MSLP
# for i, ax in enumerate(axs[num_subplots:]):
#     clevs = np.linspace(-2, 2.00, 21)
#     ax.coastlines()

#     regressions_sf_xr[i,15:185,105:345].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
#                                 transform=ccrs.PlateCarree(), add_colorbar=False)

#     ax.set_title(f'PC{corr_sf_xr["PC"][i].values + 1} MSLP', fontsize=20)
# #plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
# cax1 = fig.add_axes([1, 0.05, 0.02, 0.9])  # Adjust the position and size of the colorbar
# colorbar1 = Colorbar(ax=cax1, mappable=corr_xr[0].plot.contourf(ax=axs[0], levels = clevs, 
#                                                                 cmap=plt.cm.RdBu_r, 
#                                                                 add_colorbar=False), orientation='vertical')

# # Adjust layout for better spacing
# plt.tight_layout()

# plt.title('PCs and 3-months SPI-value regression slope all year 1950-2022', fontsize = 44, loc='right', pad=33)

# # Show the plot
# plt.show()

