# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:16:42 2023

@author: tojo1
"""

import os
import cdsapi
import zipfile
#import xarray as xr
#import netCDF4

#import xarray as xr
import cdsapi

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\pressure_data')

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'grib',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'geopotential',
        'pressure_level': '500',
        'year': [
            '1950', '1951',
            '1952', '1953', '1954',
            '1955', '1956', '1957',
            '1958', '1959', '1960',
            '1961', '1962', '1963',
            '1964', '1965', '1966',
            '1967', '1968', '1969',
            '1970', '1971', '1972',
            '1973', '1974', '1975',
            '1976', '1977', '1978',
            '1979', '1980', '1981',
            '1982', '1983', '1984',
            '1985', '1986', '1987',
            '1988', '1989', '1990',
            '1991', '1992', '1993',
            '1994', '1995', '1996',
            '1997', '1998', '1999',
            '2000', '2001', '2002',
            '2003', '2004', '2005',
            '2006', '2007', '2008',
            '2009', '2010', '2011',
            '2012', '2013', '2014',
            '2015', '2016', '2017',
            '2018', '2019', '2020',
            '2021', '2022', '2023',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'area': [
            85, -80, 20,
            40,
        ],
    },
    'ERA5_GP500_50_2023.grib')

# c.retrieve(
#     'reanalysis-era5-single-levels-monthly-means',
#     {
#         'format': 'grib',
#         'product_type': 'monthly_averaged_reanalysis',
#         'variable': 'total_precipitation',
#         'year': [
#             '1950', '1951', '1952',
#             '1953', '1954', '1955',
#             '1956', '1957', '1958',
#             '1959', '1960', '1961',
#             '1962', '1963', '1964',
#             '1965', '1966', '1967',
#             '1968', '1969', '1970',
#             '1971', '1972', '1973',
#             '1974', '1975', '1976',
#             '1977', '1978', '1979',
#             '1980', '1981', '1982',
#             '1983', '1984', '1985',
#             '1986', '1987', '1988',
#             '1989', '1990', '1991',
#             '1992', '1993', '1994',
#             '1995', '1996', '1997',
#             '1998', '1999', '2000',
#             '2001', '2002', '2003',
#             '2004', '2005', '2006',
#             '2007', '2008', '2009',
#             '2010', '2011', '2012',
#             '2013', '2014', '2015',
#             '2016', '2017', '2018',
#             '2019', '2020', '2021',
#             '2022', '2023',
#         ],
#         'area': [
#             55, 5, 50,
#             10,
#         ],
#         'time': '00:00',
#         'month': [
#             '01', '02', '03',
#             '04', '05', '06',
#             '07', '08', '09',
#             '10', '11', '12',
#         ],
#     },
#     'ERA5_precip_1950_2023_5_x_5_degree_box.grib')

### chose the data at the CDS website and it gives you the code, that can be used beneath

#c.retrieve('sis-biodiversity-era5-global',
#    {
#        'format': 'zip',
#        'temporal_aggregation': 'annual',
#        'variable': 'volumetric_soil_water',
#        'derived_variable': 'coldest_quarter',
#        'version': '1.0',
#    },
#    'download.zip')
#
#
#with zipfile.ZipFile('download.zip', 'r') as zip_ref:
#    zip_ref.extractall(r'C:\Users\tojo1\Documents\Speciale\Data')
#    
#os.remove('download.zip')


