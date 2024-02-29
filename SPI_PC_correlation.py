# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 09:57:34 2023

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

start_time = time.time()
#%%
## Calculate PC for Z500 correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

# Define number of EOF's wanted
neof = 4
# Read data from existing daatafile
ERA5_data = xr.load_dataset('ERA5_500hPa_GP_NDJFM_1950_2021.grib', engine = 'cfgrib')
# Extract geopotential from datafile
z_values = ERA5_data['z']

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
#%%
## calculate PC for MSLP correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

neof = 4
# Read data from existing daatafile
ERA5_datas = xr.load_dataset('ERA5_MSLP_NDJFM_1950_2021.grib', engine = 'cfgrib')
# Extract geopotential from datafile
sp_values = ERA5_datas['sp']

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
    ax.set_title('EOF' + str(eof_res_array['EOF'][i-1].values) + ' GP500 anomaly NDJFM 50-21', fontsize=16)
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
    ax.set_title('EOF' + str(eof_res_array_sf['EOF'][i-1].values) + ' MSLP anomaly NDJFM 50-21', fontsize=16)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

#%%
## Calculate SPI-values - import data - dataset is large! - Change dataset!!!
## change directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc', engine = 'netcdf4')

#%%
## Calculate SPI-values on daily time series for monthly SPI 

### OBS! Can be edited for monthly date imput...

time = E_obs['time'][:]
longitude = E_obs['longitude']#[180:222]
latitude = E_obs['latitude']#[100:142]
precip = E_obs['rr']#[:,100:142,180:222]

E_obs.close()
# calculate monthly SPI-values from daily time series

## define perdiod for calculating values
months = 1
t = months - 1

## resample time variable for monthly values
#time_months = time.resample(time="1M").max(dim="time")

## create array for resulting values
spi_monthly_values = np.empty((len(time), len(latitude), len(longitude)))

x = 0

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
        #monthly_time_series = filteret_time_series.resample('M').sum()
        ## Roll data for chosen period
        rolling_time_series = filteret_time_series.rolling(months, min_periods = months).sum()
        if rolling_time_series.sum() == 0.0:
            print('zeros')
            spi_monthly_values[:,lat_idx,lon_idx] = np.nan
        else:
            rolling_time_series = rolling_time_series.dropna()
            ## calculate SPI if there is not too much data missing
            spi_calc = si.spi(rolling_time_series, dist=scs.lognorm)
            spi_monthly_values[:t,lat_idx,lon_idx] = np.nan
            spi_monthly_values[t:,lat_idx,lon_idx] = spi_calc
        x = x+1
        print('SPI: ' + str(x*100/(len(latitude)*len(longitude)*len(time))) + '%')
        
## Write data to xarray dataarray
spi_monthly_array = xr.Dataset(
    {'SPI' : xr.DataArray(
        data = np.array(spi_monthly_values),
        dims =('time','latitude', 'longitude',),
        coords={'latitude': latitude, 'longitude': longitude, 'time': time},
        attrs={'description': 'Standardized Precipitation Index.', 'units': 'Relative scale'})})



#%%
## Calculate correlation between SPI and PC's

lats = spi_monthly_array['latitude']
longs = spi_monthly_array['longitude']

## select month for analysis - summer/winter/full year
selected_SPI_data = spi_monthly_array.sel(time=spi_monthly_array['time.month'].isin([1, 2, 3, 11, 12]))

### Calculating the correaltion for the 500 hPa level

## Find correlations on of PC and the SPI-values
corr = correlation_map(pcs,selected_SPI_data['SPI'][:360].values)

## correlation data as Xarray dataarray
corr_xr = xr.DataArray(
    corr,
    coords={'PC': [0,1,2,3], 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

## plotting the
for i in corr_xr['PC']:
    clevs = np.linspace(-1, 1, 19)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI NDJFM 50-21 corr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

### Calculate the correlations for the MSLP

## Find correlations on of PC and the SPI-values
corr_sf = correlation_map(pcs_sf[:360,:],selected_SPI_data['SPI'][:360].values)

## correlation data as Xarray dataarray
corr_sf_xr = xr.DataArray(
    corr_sf,
    coords={'PC': [0,1,2,3], 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

## plotting the
for i in corr_sf_xr['PC']:
    clevs = np.linspace(-1, 1, 19)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_sf_xr['PC'][i].values+1) + ' MSLP and SPI NDJFM 50-21 corr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

#end_time = time.time()
#elapsed_time = end_time - start_time
#print('calculated in ' + str(elapsed_time) + ' seconds')

#%%

## Limit data to only high correlations

## high correlation level
H = 0.3

high_corr = ma.masked_inside(corr_xr, -H, H)

high_corr_xr = xr.DataArray(
    high_corr,
    coords={'PC': [0,1,2,3], 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

for i in high_corr_xr['PC']:
    clevs = np.linspace(-1, 1, 19)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    high_corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI NDJFM 50-21 Hcorr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

high_corr_sf = ma.masked_inside(corr_sf_xr, -H, H)

high_corr_sf_xr = xr.DataArray(
    high_corr_sf,
    coords={'PC': [0,1,2,3], 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

for i in corr_sf_xr['PC']:
    clevs = np.linspace(-1, 1, 19)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    high_corr_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI NDJFM 50-21 Hcorr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()
