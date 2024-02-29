# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:47:55 2023

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

start_time = time.time()
#%%
## Calculate PC for Z500 correlation analysis
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

months = [1,2,3,11,12]

period = str('NDJFM 1950-2022')

# Define number of EOF's wanted
neof = 4
# Read data from existing daatafile
ERA5_data = xr.load_dataset('ERA5_GP500_50_2023.grib', engine = 'cfgrib')
# Extract geopotential from datafile
z_values = ERA5_data['z']

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
    ax.set_title('EOF' + str(eof_res_array['EOF'][i-1].values) + ' GP500 anomaly 50-22', fontsize=16)
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
    ax.set_title('EOF' + str(eof_res_array_sf['EOF'][i-1].values) + ' MSLP anomaly 50-22', fontsize=16)
    plot_text = str(round(varfrac_sf[i-1]*100, 1)) + '%'
    plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('EOF plots')
#%%
## Calculate SPI-values - import data - dataset is large! - Change dataset!!!
## change directory
os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc', engine = 'netcdf4')



#%%
### Import previously calculated SPI-value to reduce runtime

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
spi_monthly_array = xr.load_dataset('E_obs_MONTHLY_spi_1M_1950_2022.nc')

print('SPI values loaded')

#%%
## plotting SPI values for the first year

plotting_data = spi_monthly_array['SPI'][:len(months)].mean(dim = 'time')

clevs = np.linspace(-1, 1.00, 21)
proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
ax = plt.axes(projection=proj)

ax.coastlines()

plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                         transform=ccrs.PlateCarree(), add_colorbar=True)
#plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
ax.set_title('Mean SPI values for JFMND 1950', fontsize=16)
#plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
plt.show()


#%%
## Calculate correlation between SPI and PC's

lats = spi_monthly_array['latitude']
longs = spi_monthly_array['longitude']

## select month for analysis - summer/winter/full year
selected_SPI_data = spi_monthly_array.sel(time=spi_monthly_array['time.month'].isin(months))

### Calculating the correaltion for the 500 hPa level

## Find correlations on of PC and the SPI-values
corr = correlation_map(pcs[:len(selected_SPI_data['time'])],selected_SPI_data['SPI'][:len(pcs)].values)

## correlation data as Xarray dataarray
corr_xr = xr.DataArray(
    corr,
    coords={'PC': np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'correlation coefficient'},)

## plotting the
for i in corr_xr['PC']:
    clevs = np.linspace(-1, 1.00, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    corr_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 50-22 corr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

### Calculate the correlations for the MSLP

## Find correlations on of PC and the SPI-values
corr_sf = correlation_map(pcs_sf[:len(selected_SPI_data['time']),:],selected_SPI_data['SPI'][:len(pcs)].values)

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
    ax.set_title('PC' + str(corr_sf_xr['PC'][i].values+1) + ' MSLP and SPI 50-22 corr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('surface correlations calculation done')

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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 50-22 Hcorr', fontsize=16)
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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 50-22 Hcorr', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('High correlations plottet')
#%%

# Calculate the regression coefficient for the correlations

exp_coef = pcs[:len(selected_SPI_data['time']),:]/pcs[:len(selected_SPI_data['time']),:].std()

regressions = np.empty((neof, len(lats), len(longs)))

x = 0 
for i in range(neof):
    for lat in range(len(lats)):
        for lon in range(len(longs)):
            calc_data = selected_SPI_data['SPI'][:,lat,lon]
            if np.isnan(calc_data).any() ==  True:
                regressions[i, lat, lon] = np.nan
            else:
                polifit = poly.polyfit(exp_coef[:,i],calc_data,1)
                regressions[i, lat, lon] = polifit[1]
            print('500 hPa ' + str(x/(93264*4)))
            x = x+1
            
regressions_xr = xr.DataArray(
    regressions,
    coords={'PC':np.arange(0,neof), 'latitude' : lats, 'longitude' : longs},  # Use your own coordinate names and values
    dims=('PC', 'latitude', 'longitude'),
    attrs={'description': 'regression coefficient'},)

print('500 hPa regressions calculated')

regressions_sf = np.empty((neof, len(lats), len(longs)))

exp_coef_sf = pcs_sf[:len(selected_SPI_data['time']),:]/pcs_sf[:len(selected_SPI_data['time']),:].std()

x = 0 
for i in range(neof):
    for lat in range(len(lats)):
        for lon in range(len(longs)):
            calc_data = selected_SPI_data['SPI'][:,lat,lon]
            if np.isnan(calc_data).any() ==  True:
                regressions_sf[i, lat, lon] = np.nan
            else:
                polifit = poly.polyfit(exp_coef_sf[:,i],calc_data,1)
                regressions_sf[i, lat, lon] = polifit[1]
            print('surface ' + str(x/(93264*4)))
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
    clevs = np.linspace(-2, 2, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' Z500 and SPI 50-22 RegSlope', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('look at those plots!')

for i in regressions_sf_xr['PC']:
    clevs = np.linspace(-2, 2, 21)
    proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
    ax = plt.axes(projection=proj)
    ax.coastlines()
    #ax.set_global()
    regressions_sf_xr[i].plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True)
    #plot_text = str(round(varfrac[i-1]*100, 1)) + '%'
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 50-22 RegSlope', fontsize=16)
    #plt.text(-5100000, -4800000, plot_text, fontsize=12, color='black')
    plt.show()

print('More beautiful plots are made!')

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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 50-22 corr', fontsize=16)
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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 50-22 corr', fontsize=16)
    plt.show()

print('Hatched in stead of masked high correlations')

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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' GP500 and SPI 50-22 RegSlope', fontsize=16)
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
    ax.set_title('PC' + str(corr_xr['PC'][i].values+1) + ' MSLP and SPI 50-22 RegSlope', fontsize=16)
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
eof_res_array.to_netcdf('500hPa_correlation_analysis.nc')

corr_sf_xr_edit = corr_sf_xr
corr_sf_xr_edit['latitude'] = corr_sf_xr_edit['latitude'] + 0.125
corr_sf_xr_edit['longitude'] = corr_sf_xr_edit['longitude'] + 0.125

regressions_sf_xr_edit = regressions_sf_xr
regressions_sf_xr_edit['latitude'] = regressions_sf_xr_edit['latitude'] + 0.125
regressions_sf_xr_edit['longitude'] = regressions_sf_xr_edit['longitude'] + 0.125

eof_res_array_sf['corr'] = corr_sf_xr_edit
eof_res_array_sf['reg'] = regressions_sf_xr_edit
eof_res_array_sf.to_netcdf('MSLP_correlation_analysis.nc')
