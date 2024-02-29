# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 10:10:48 2023

@author: tojo1
"""

import os
import xarray as xr

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data')
## Load dataset
E_obs = xr.load_dataset('E_obs_1980_2010/E_obs_precip_1950_2022.nc', engine = 'netcdf4')

E_obs_monthly = E_obs.resample(time='M').sum(dim='time')

E_obs_monthly.to_netcdf('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc')

monthly_precip = xr.load_dataset('E_obs_1980_2010/E_obs_MONTHLY_precip_1950_2022.nc')


os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')

spi_monthly_array.to_netcdf('E_obs_MONTHLY_spi_1M_1950_2022.nc')

monthly_spi = xr.load_dataset('E_obs_MONTHLY_spi_1M_1950_2022.nc')
