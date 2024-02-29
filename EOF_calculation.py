# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:41:17 2023

@author: tojo1
"""

import cartopy.crs as ccrs
#from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
from eofs.standard import Eof

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

# Define number of EOF's wanted
neof = 4
# Read data from existing daatafile
ERA5_data = xr.load_dataset('ERA5_500hPa_GP_NDJFM_1950_2021.grib', engine = 'cfgrib')
# Extract geopotential from datafile
z_values = ERA5_data['z']

# Compute anomalies by removing the time-mean.
#z_values = z_values - z_values.mean(dim='time')

z_values_norm = z_values / z_values.mean(dim = 'latitude').mean(dim = 'longitude')

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(z_values.coords['latitude'].values)).clip(0., 1.)
wgts = np.sqrt((coslat))[..., np.newaxis]
solver = Eof(z_values_norm.values, weights=wgts)

# Retrieve the leading EOF's, expressed as the covariance between the leading 
# PC time series and the input SLP anomalies at each grid point.
eof_res = solver.eofsAsCovariance(neofs=neof)

eofs = solver.eofs(neofs = neof, eofscaling=0)
pcs = solver.pcs(npcs = neof)

varfrac = solver.varianceFraction(neof)
#varfrac_corrected = varfrac[1:]/(1-varfrac[0])
errors = solver.northTest(neigs=neof, vfscaled=True)
eigenvalues = solver.eigenvalues(neigs=neof)

# Convert values back to xarray
#eof_res_array = xr.DataArray(eofs,coords={'y': z_values['latitude'].values,'x': z_values['longitude'].values,'neof': neof}, 
#dims=['EOF',"y", "x"],attrs={'description': 'S', 'units': 'Relative scale'})

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
        attrs={'description': 'PCS at 500 hPa level', 'units': 'Relative scale'})})
#eof_res_array = eof_res_array * -1


#### EOF 0 is the mean value of the preassure fields!!!!
# Plot preassure anomaly maps for the governing EOF's
for i in eof_res_array['EOF']:
    clevs = np.linspace(-0.009, 0.009, 19)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    eof_res_array['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('EOF' + str(eof_res_array['EOF'][i-1].values) + ' GP500 anomaly NDJFM 50-21', fontsize=16)
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

ERA5_data.close()

#%%
## for mean sea level preassure
neof = 4
# Read data from existing daatafile
ERA5_datas = xr.load_dataset('ERA5_MSLP_NDJFM_1950_2021.grib', engine = 'cfgrib')
# Extract geopotential from datafile
sp_values = ERA5_datas['sp']

# Compute anomalies by removing the time-mean.
#sp_values = sp_values - sp_values.mean(dim='time')

#sp_values_norm = sp_values / sp_values.std(dim = 'time')
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

# Convert values back to xarray
eof_res_array_sf = xr.DataArray(eofs,coords={'y': sp_values['latitude'].values,'x': sp_values['longitude'].values,'neof': neof}, 
dims=['EOF',"y", "x"])

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
#eof_res_array = eof_res_array * -1

# Plot preassure anomaly maps for the governing EOF's
for i in eof_res_array_sf['EOF']:
    clevs = np.linspace(-0.009, 0.009, 19)
    proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    eof_res_array_sf['modes'][i-1].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    ax.set_title('EOF' + str(eof_res_array_sf['EOF'][i-1].values) + ' MSLP anomaly NDJFM 59-14', fontsize=16)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

ERA5_datas.close()
