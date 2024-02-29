# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 14:56:29 2024

@author: tojo1
"""
import cartopy.crs as ccrs
#from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
from matplotlib.colorbar import Colorbar
from libpysal.weights import lat2W
import esda
import libpysal
import numpy.ma as ma

# %%

# os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

# SPI_periods = [1,3,6,12,24]

# for SPI in SPI_periods:
    
#     corr_xr_3_month = xr.open_dataset('500hPa_correlation_analysis_' + str(SPI) + '_month_extended_summer_season.nc')
    
#     # cor_for_moran = corr_xr_3_month['corr'][:,:,250:]
    
#     ## Make a moving window used for moran's I calculation
#     moving_window_sizes = [3,5,7,9,11,13]
    
#     PCs = [1,2,3,4]
    
#     ## drop rows and columns without data
#     cor_for_moran = corr_xr_3_month['corr'].dropna(dim='latitude', how = 'all')
#     cor_for_moran = cor_for_moran.dropna(dim='longitude', how = 'all')
#     non_nan_mask = cor_for_moran[0].notnull()
    
#     for PC in PCs:
#         t = 0
#         a = 0
#         PC_morans = []
#         Morans_significant = []
#         for moving_window_size in moving_window_sizes:
#             Morans = np.full((len(cor_for_moran['latitude']), len(cor_for_moran['longitude'])), np.nan)
#             # Morans_EI = np.full((len(cor_for_moran['latitude']), len(cor_for_moran['longitude'])), np.nan)
#             # Morans_ZS = np.full((len(cor_for_moran['latitude']), len(cor_for_moran['longitude'])), np.nan)
#             Morans_pv = np.full((len(cor_for_moran['latitude']), len(cor_for_moran['longitude'])), np.nan)
            
#             low = int(moving_window_size/2)
#             hi = round(moving_window_size/2)
#             x = 0
#             ## Calculate Moran's I for every gridcell
#             for lat in range(len(cor_for_moran['latitude'])):
#                 for lon in range(len(cor_for_moran['longitude'])):
            
#                     WeightMatrix = lat2W(cor_for_moran[PC-1,lat-low:lat+hi,lon-low:lon+hi].shape[0],
#                                           cor_for_moran[PC-1,lat-low:lat+hi,lon-low:lon+hi].shape[1])
                    
#                     moran_test = esda.Moran(cor_for_moran[PC-1,lat-low:lat+hi,lon-low:lon+hi].fillna(0), WeightMatrix)
#                     # print(moran_test.I)
#                     Morans[lat, lon] = moran_test.I
#                     # print(moran_test.EI_sim)
#                     # Morans_EI[lat, lon] = moran_test.EI
#                     # Morans_ZS[lat, lon] = moran_test.z_sim
#                     if moran_test.I != np.nan:
#                         Morans_pv[lat, lon] = moran_test.p_sim
                    
#                 x=x+1
#                 # print('Morans I: ' + str( x * 100 / (len(cor_for_moran['latitude']))) + '%')
            
#             PC_morans.append(xr.DataArray(
#                 Morans.copy(),
#                 coords={'latitude' : cor_for_moran['latitude'], 'longitude' : cor_for_moran['longitude']},  
#                 dims=('latitude', 'longitude'),
#                 attrs={'description': 'Morans I autocorrelation'},))
            
#             Morans_significant.append(xr.DataArray(
#                 Morans_pv.copy(),
#                 coords={'latitude' : cor_for_moran['latitude'], 'longitude' : cor_for_moran['longitude']},  
#                 dims=('latitude', 'longitude'),
#                 attrs={'description': 'P-value for Morans I significance'},))
           
#             # Morans_significant[a] =  Morans_significant[a].where(Morans_significant[a] < 0.05,1,np.nan)
            
#             # Morans_significant_cop = Morans_significant[a].copy()
#             # P_values = Morans_significant_cop.where(Morans_significant_cop < 0.05, drop = True)
            
#             fig = plt.figure(figsize=(12, 8))
#             clevs = np.linspace(-1, 1, 21)
#             proj = ccrs.Orthographic(central_longitude=17.5, central_latitude=50.375)
#             ax = plt.axes(projection=proj)
#             ax.coastlines()
#             PC_morans[a].where(non_nan_mask, drop=True).plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
#                                       transform=ccrs.PlateCarree(), add_colorbar=True,
#                                       cbar_kwargs={'label': ''})
#             masking = Morans_significant[a].copy()
#             masking = masking.where(np.isnan(PC_morans[a]) == False)
#             masking = ma.masked_where(masking > 0.01, masking)
            
#             mask_xr = xr.DataArray(
#                 masking,
#                 coords={'latitude' : cor_for_moran['latitude'], 'longitude' : cor_for_moran['longitude']},  
#                 dims=('latitude', 'longitude'),
#                 attrs={'description': 'P-value for Morans I significance'},)
#             mask = mask_xr.notnull()
#             hatches = np.full(PC_morans[a].shape, '', dtype='U1')
#             hatching_pattern = '....'
#             # Set hatching pattern in the masked region
#             hatches[mask] = hatching_pattern
#             hatches = np.where(mask, hatching_pattern, hatches)
#             mask_xr.where(non_nan_mask, drop=True).plot.contourf(ax=ax, levels=clevs, 
#                                             cmap=plt.cm.RdBu,
#                                             transform = ccrs.PlateCarree(), 
#                                             hatches = [hatching_pattern], 
#                                             add_colorbar=False,alpha =0)
#             ax.set_title('SPI-' + str(SPI) + ' PC' + str(PC) + ' ' + str(moving_window_sizes[a]) + ' x ' + str(moving_window_sizes[a]) + " gridcells - Mean: " + str(round(float(PC_morans[a].mean()), 3)), fontsize=17)
#             plt.show()
            
#             if moving_window_size < 8:
#                 b = 0
#             else:
#                 b = 1
#             if moving_window_size == 9:
#                 t = 0
#             elif moving_window_size == 3:
#                 t = 0
           
#             t = t+1
#             a = a+1
            
#             print('time for another Morans I calculation')
        
#         num_subplots = 3

#         # Create a 2x4 grid of subplots
#         fig, axs = plt.subplots(2, num_subplots, figsize=(15, 9.3), 
#                                 subplot_kw={'projection': ccrs.Orthographic(central_longitude=17.5, 
#                                                                             central_latitude=50.375)})

#         # Flatten the 2D array of subplots into a 1D array
#         axs = axs.flatten()

#         # Plotting correlations for GP500
#         for i, ax in enumerate(axs):
#             clevs = np.linspace(-1, 1.00, 21)
#             ax.coastlines()
           
#             PC_morans[i].where(non_nan_mask, drop=True).plot.contourf(ax=ax, levels=clevs, cmap=plt.cm.RdBu,
#                                       transform=ccrs.PlateCarree(), add_colorbar=False)
            
#             masking = Morans_significant[i].copy()
#             masking = masking.where(np.isnan(PC_morans[i]) == False)
#             masking = ma.masked_where(masking > 0.01, masking)
            
#             mask_xr = xr.DataArray(
#                 masking,
#                 coords={'latitude' : cor_for_moran['latitude'], 'longitude' : cor_for_moran['longitude']},  
#                 dims=('latitude', 'longitude'),
#                 attrs={'description': 'P-value for Morans I significance'},)
#             mask = mask_xr.notnull()
#             hatches = np.full(PC_morans[i].shape, '', dtype='U1')
#             hatching_pattern = '....'
#             # Set hatching pattern in the masked region
#             hatches[mask] = hatching_pattern
#             hatches = np.where(mask, hatching_pattern, hatches)
#             mask_xr.where(non_nan_mask, drop=True).plot.contourf(ax=ax, levels=clevs, 
#                                             cmap=plt.cm.RdBu,
#                                             transform = ccrs.PlateCarree(), 
#                                             hatches = [hatching_pattern], 
#                                             add_colorbar=False,alpha =0)
            
#             subtitle = str(moving_window_sizes[i]) + ' x ' + str(moving_window_sizes[i]) + " gridcells - Mean: " + str(round(float(PC_morans[i].mean()), 3))
#             ax.set_title(subtitle, fontsize=20)
       
#         #plt.title('PCs and 3-months SPI values MJJAS 1950-2022', fontsize = 30, loc='right', pad=150)
#         cax1 = fig.add_axes([1, 0.05, 0.03, 0.89])  # Adjust the position and size of the colorbar
#         colorbar1 = Colorbar(ax=cax1, mappable=PC_morans[0].plot.contourf(levels = clevs, 
#                                                                         cmap=plt.cm.RdBu, 
#                                                                         add_colorbar=False), orientation='vertical')
#         cax1.xaxis.label.set_text('')
#         cax1.title.set_text('')

#         # Adjust layout for better spacing
#         plt.tight_layout()
    
#         plt.title('SPI-' + str(SPI) + ' PC' + str(PC) + " Moran's I autocorrelation varying moving window size", fontsize = 33, loc='right', pad=40)
    
#         plt.show()
        
        
#%%

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

seasons = ['MJJAS', 'NDJFM', 'All year']

season = seasons[2]

SPI_periods = [12]

PCs = [1,2,3,4]

Morans = np.zeros((len(SPI_periods), len(PCs)))
Morans_p = np.zeros((len(SPI_periods), len(PCs)))
Morans_z = np.zeros((len(SPI_periods), len(PCs)))

# for season in seasons:
x = 0
for SPI in SPI_periods:
    y = 0
    corr_xr_3_month = xr.open_dataset('500hPa_correlation_analysis_12_month_EOF NDJFM and SPI MJJAS_shift.nc')
    for PC in PCs:
        
        Moran_corr = corr_xr_3_month['corr'][PC-1]
        
        latitude = corr_xr_3_month['latitude']
        longitude = corr_xr_3_month['longitude']
        
        WeightMatrix = lat2W(len(latitude), len(longitude))
        
        moran_test = esda.Moran(Moran_corr.fillna(0), WeightMatrix)
        
        Morans[x,y] = moran_test.I
        Morans_p[x,y] = moran_test.p_sim
        Morans_z[x,y] = moran_test.z_sim
        print(moran_test.I)
        print(moran_test.p_sim)
        print(moran_test.z_sim)
        
        y = y + 1
    x = x + 1
