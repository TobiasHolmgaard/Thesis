# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:00:34 2023

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
from matplotlib.colors import BoundaryNorm
from matplotlib.colorbar import ColorbarBase

#%%

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

ERA5_data = xr.load_dataset('ERA5_GP500_50_2023.grib', engine = 'cfgrib')
ERA5_data_sf = xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
ERA5_data_sf['msl'] = ERA5_data_sf['msl']/100

EOF_z500_summer = xr.load_dataset('500hPa_correlation_analysis_1_month_extended_summer_season.nc')
EOF_z500_winter = xr.load_dataset('500hPa_correlation_analysis_1_month_extended_winter_season.nc')
EOF_z500_year = xr.load_dataset('500hPa_correlation_analysis_1_month_extended_winter_season.nc')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')


os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
SPI1 = xr.load_dataset('SPI_1_month_1950_2022.nc')
SPI3 = xr.load_dataset('SPI_3_month_1950_2022.nc')
SPI6 = xr.load_dataset('SPI_6_month_1950_2022.nc')
SPI12 = xr.load_dataset('SPI_12_month_1950_2022.nc')
SPI24 = xr.load_dataset('SPI_24_month_1950_2022.nc')

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc', engine = 'netcdf4')

zero_value_mask = E_obs['rr'].mean(dim = 'time').copy()

for lat_idx in range(len(SPI1['latitude'])):
    for lon_idx in range(len(SPI1['longitude'])):
        zeroser = (E_obs['rr'][:,lat_idx,lon_idx] == 0).sum(dim='time')
        if zeroser > 0.22 * len(SPI1['time']):
            zero_value_mask[lat_idx, lon_idx] = np.nan
        else: 
            zero_value_mask[lat_idx, lon_idx] = 1
            
SPI1 = SPI1.where(zero_value_mask == 1, drop = True)
SPI3 = SPI3.where(zero_value_mask == 1, drop = True)
SPI6 = SPI6.where(zero_value_mask == 1, drop = True)
SPI12 = SPI12.where(zero_value_mask == 1, drop = True)
SPI24 = SPI24.where(zero_value_mask == 1, drop = True)

#%%

SPIs = [SPI1, SPI3, SPI6, SPI12, SPI24]
SN = [1,3,6,12,24]

start_month = 5
end_month = 9

length = end_month-start_month+1

drought_years = [1953, 1954, 1959, 1960, 1969, 1976, 1985, 1989, 1990, 1991, 2002,
                 2003, 2004, 2007, 2008, 2010, 2015, 2016, 2017, 2018]

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(11, 9))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        plotting_data = plotting_data.dropna(dim='longitude', how = 'all').dropna(dim='latitude', how = 'all')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=20, central_latitude=50.375)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=False)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 MJJAs: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 MJJAS: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-18, 18, 19)
        ERA5_data_sf_stand = ERA5_data_sf['msl'] - ERA5_data_sf['msl'].mean(dim='longitude').mean(dim='latitude')
        # data_plot = ERA5_data_sf_stand[start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time')
        # data_plot['longitude'] = data_plot['longitude'] - 0.125
        # data_plot['latitude'] = data_plot['latitude'] - 0.125
        # data_plot = data_plot.where(plotting_data, drop=True)
        # contour_plot = data_plot.plot.contour(ax=ax, levels=clevs_eof,
                                                                                                                            
        contour_plot = ERA5_data_sf_stand[start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 8)
        plt.xlim(-4300000, 1850000)
        plt.ylim(-2000000, 2500000)
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-4130000, 1800000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=12, color='black')
        ax.set_title('Mean SPI-' + str(SN[x]) + ' values for MJJAS ' + str(year), fontsize=28)
        plt.show()
        x=x+1
    y = y+1


fig = plt.figure(figsize=(11, 9))
cmap = plt.cm.RdBu
fig, ax = plt.subplots()
norm = BoundaryNorm(clevs, cmap.N)
# Create the colorbar
colorbar = ColorbarBase(ax=ax, cmap=cmap, norm=norm, orientation='horizontal')
colorbar.set_label('SPI-values')
ax.set_aspect(0.1)
plt.show()

fig = plt.figure(figsize=(11, 9))
cmap = plt.cm.RdBu
fig, ax = plt.subplots()
norm = BoundaryNorm(clevs, cmap.N)
# Create the colorbar
colorbar = ColorbarBase(ax=ax, cmap=cmap, norm=norm, orientation='vertical')
colorbar.set_label('SPI-values')
ax.set_aspect(10)
plt.show()

#%%

start_month = - 1
end_month = 3

length = end_month - start_month

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(11,9))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=20, central_latitude=50.375)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=False)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-18, 18, 19)
        ERA5_data_sf_stand = ERA5_data_sf['msl'] - ERA5_data_sf['msl'].mean(dim = 'longitude').mean(dim = 'latitude')
        contour_plot = ERA5_data_sf_stand[start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 8)
        plt.xlim(-4300000, 1850000)
        plt.ylim(-2000000, 2500000)
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-4130000, 1800000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=12, color='black')
        ax.set_title('Mean SPI-' + str(SN[x]) + ' values for NDJFM ' + str(year-1) + '/' + str(year), fontsize=28)
        plt.show()
        x=x+1
    y = y+1

fig = plt.figure(figsize=(11, 9))
cmap = plt.cm.RdBu
fig, ax = plt.subplots()
norm = BoundaryNorm(clevs, cmap.N)
# Create the colorbar
colorbar = ColorbarBase(ax=ax, cmap=cmap, norm=norm, orientation='horizontal')
colorbar.set_label('SPI-values', fontsize = 25)
ax.set_aspect(0.1)
plt.show()

fig = plt.figure(figsize=(11, 9))
cmap = plt.cm.RdBu
fig, ax = plt.subplots()
norm = BoundaryNorm(clevs, cmap.N)
# Create the colorbar
colorbar = ColorbarBase(ax=ax, cmap=cmap, norm=norm, orientation='vertical')
colorbar.set_label('SPI-values')
ax.set_aspect(10)
plt.show()

#%%

start_month = 0
end_month = 12

length = end_month - start_month

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(12, 6.8))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=20, central_latitude=50.375)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=True)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4 NDJFM: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-5, 5, 15)
        ERA5_data_sf_stand = ERA5_data_sf['msl'] - ERA5_data_sf['msl'].mean(dim='time')
        contour_plot = ERA5_data_sf_stand[start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 8)
        plt.xlim(-4300000, 1850000)
        plt.ylim(-2000000, 2500000)
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-4130000, 1800000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=11, color='black')
        ax.set_title('Mean SPI values ' + str(SN[x]) + '-month for All year ' + str(year), fontsize=26)
        plt.show()
        x=x+1
    y = y+1



#%%

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

ERA5_data = xr.load_dataset('ERA5_GP500_50_2023.grib', engine = 'cfgrib')
ERA5_data_sf = xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# ERA5_data_sf['msl'] = ERA5_data_sf['msl']/100
# ERA5_data = ERA5_data

EOF_z500_summer = xr.load_dataset('500hPa_correlation_analysis_1_month_extended_summer_season.nc')
EOF_z500_winter = xr.load_dataset('500hPa_correlation_analysis_1_month_extended_winter_season.nc')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')
# PI24= xr.load_dataset('ERA5_MSLP_50_2023.grib', engine = 'cfgrib')


os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
SPI1 = xr.load_dataset('SPI_1_month_1950_2022.nc')
SPI3 = xr.load_dataset('SPI_3_month_1950_2022.nc')
SPI6 = xr.load_dataset('SPI_6_month_1950_2022.nc')
SPI12 = xr.load_dataset('SPI_12_month_1950_2022.nc')
SPI24 = xr.load_dataset('SPI_24_month_1950_2022.nc')

SPI1 = SPI1.where(zero_value_mask == 1, drop = True)
SPI3 = SPI3.where(zero_value_mask == 1, drop = True)
SPI6 = SPI6.where(zero_value_mask == 1, drop = True)
SPI12 = SPI12.where(zero_value_mask == 1, drop = True)
SPI24 = SPI24.where(zero_value_mask == 1, drop = True)


SPIs = [SPI1, SPI3, SPI6, SPI12, SPI24]
SN = [1,3,6,12,24]

start_month = 5
end_month = 9

length = end_month-start_month+1

drought_years = [1953, 1954, 1959, 1960, 1969, 1976, 1985, 1989, 1990, 1991,
                 2003, 2007, 2008, 2010, 2015, 2016, 2017, 2018]

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(12, 6.8))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=True)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_summer['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-1500, 1500, 16)
        ERA5_data_stand = ERA5_data - ERA5_data.mean(dim='time')
        contour_plot = ERA5_data_stand['z'][start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 7)
    
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-5000000, 1930000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=11, color='black')
        ax.set_title('Mean SPI values ' + str(SN[x]) + '-month for MJJAS ' + str(year), fontsize=25)
        plt.show()
        x=x+1
    y = y+1


start_month = - 1
end_month = 3

length = end_month - start_month

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(12, 6.8))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=True)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-1500, 1500, 16)
        ERA5_data_stand = ERA5_data - ERA5_data.mean(dim='time')
        contour_plot = ERA5_data_stand['z'][start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 7)
    
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-5000000, 1930000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=11, color='black')
        ax.set_title('Mean SPI values ' + str(SN[x]) + '-month for NDJFM ' + str(year), fontsize=25)
        plt.show()
        x=x+1
    y = y+1

start_month = 0
end_month = 12

length = end_month - start_month

fig = plt.figure(figsize=(12, 6.8))

y = 0
for year in drought_years:
    x = 0
    for SPI in SPIs:
        fig = plt.figure(figsize=(12, 6.8))
        year_ind = (year - 1950) * 12
        year_ind_pc = (year - 1950) * (end_month - start_month + 1)
        plotting_data = SPI['SPI'][start_month - 1 + year_ind : end_month + year_ind].mean(dim = 'time')
        
        clevs = np.linspace(-2.5, 2.5, 11)
        proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
        ax = plt.axes(projection=proj)
        ax.coastlines()
        plotting_data.plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
                                 transform=ccrs.PlateCarree(), add_colorbar=True)
        
        if SN[x] == 12:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - length:year_ind_pc + length,3].mean()),2))
        elif SN[x] == 24:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc - 2*length:year_ind_pc + length,3].mean()),2))
        else:
            plot_text1 = 'PC 1: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,0].mean()),2))
            plot_text2 = 'PC 2: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,1].mean()),2))
            plot_text3 = 'PC 3: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,2].mean()),2))
            plot_text4 = 'PC 4: ' + str(round(float(EOF_z500_winter['PC_modes'][year_ind_pc:year_ind_pc + length,3].mean()),2))
        clevs_eof = np.linspace(-1500, 1500, 16)
        ERA5_data_stand = ERA5_data - ERA5_data.mean(dim='time')
        contour_plot = ERA5_data_stand['z'][start_month + year_ind - SN[x]: end_month + year_ind].mean(dim = 'time').plot.contour(ax=ax, 
                                                                                                                            levels=clevs_eof,
                                 transform=ccrs.PlateCarree(), add_colorbar=False,
                                 linewidths = 0.7, linestyles = 
                                 ['dashed' if val < 0 else 'dotted' if val > 0 
                                  else 'solid' for val in clevs_eof], 
                                 colors ='k', alpha =1)
        ax.clabel(contour_plot, fontsize = 7)
    
        plot_text = plot_text1 + '\n' + plot_text2 + '\n' + plot_text3 + '\n' + plot_text4
        plt.text(-5000000, 1930000, plot_text,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'), fontsize=11, color='black')
        ax.set_title('Mean SPI values ' + str(SN[x]) + '-month for All year ' + str(year), fontsize=25)
        plt.show()
        x=x+1
    y = y+1

