"""
 -*- coding: utf-8 -*-
Created on Tue Sep 26 14:52:49 2023

@author: tojo1
"""
import os
#import asyncio
#import logging
import requests
import json
import xarray as xr

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\spi_data')
directory = os.getcwd()

parameter = 'spi1_daily'

api_token = 'eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjMwZGI1ZDhhOGRhYjQyNjI4MjM3YTZiYTk5OWU3MDNjIiwiaCI6Im11cm11cjEyOCJ9'

url = "https://api.dataplatform.knmi.nl/open-data/v1/datasets/spi1_daily/versions/1.0/files/"
header = {"Authorization": f"Bearer {api_token}"}

#%%

## for downloading several files in a folder

page_token = None
page_size = 10

parameters = {"pageSize": page_size, "begin":'2022-06-04 00:00:00'}

y = 0
while True:
    x=0
    parameters = {"pageSize": page_size}
    url1 = url
    
    if page_token == True:
        url1 = url + '/PageToken=' + next_page_token
        parameters['pageToken'] = next_page_token
        parameters = {"pageSize": page_size, 'nextPageToken': next_page_token}
        
    
    response0 = requests.get(url, headers=header, params=parameters).json()
    filenames = response0['files']
    
    for file in filenames:
        if filenames:
            # Choose the first filename in the list for download
            download_url = url + str(filenames[x]['filename']) + '/url'
            
            # Download the file and save it locally
            response1 = requests.get(download_url, headers=header, stream=True)
            if response1.status_code == 200:
                with open('knmi_data' + str(y) + str(x) + '.json', 'wb') as file:
                    for chunk in response1.iter_content(chunk_size=8192):
                        file.write(chunk)
                    
                #print("Data downloaded and saved")
            else:
                print(f"Failed to download data (Status code: {response1.status_code})")
        else:
            print("No filenames found for download.")
        infile = 'knmi_data' + str(y) + str(x) + '.json'
        with open(infile, 'r') as json_file:
            data_from_knmi = json.load(json_file)
            
        file_download_url = data_from_knmi['temporaryDownloadUrl']
    
        if 'file_download_url' in locals():
            # Fetch the SPI data from the temporary download URL
            response2 = requests.get(file_download_url)
    
            if response2.status_code == 200:
                # Successful response, save the data to a NetCDF file
                with open('spi_data' + str(y) + str(x) + '.nc', 'wb') as file:
                    for chunk in response2.iter_content(chunk_size=8192):
                        file.write(chunk)
    
                print("SPI data downloaded and saved ")
    
            else:
                print(f"Failed to download SPI data (Status code: {response2.status_code})")
        else:
            print("File download URL not found.")
        os.remove(infile)
        x=x+1
    
     # Check for a next page token
    if "nextPageToken" in response0:
        print('token!!')
        next_page_token = response0["nextPageToken"]
        page_token = True
    else:
        break
    y = y+1


### time to merge files on time in CDO - LINUX ###

#%%

spi_dataset = xr.open_dataset('all_spi_data.nc')


#%%
### Precipitation data download setup

os.chdir(r'C:\Users\tojo1\Documents\Speciale\Data\precipitation_data\knmi_precip')

parameter = 'Rd1'

api_token = 'eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjMwZGI1ZDhhOGRhYjQyNjI4MjM3YTZiYTk5OWU3MDNjIiwiaCI6Im11cm11cjEyOCJ9'

url = "https://api.dataplatform.knmi.nl/open-data/v1/datasets/Rd1/versions/5/files"
header = {"Authorization": f"Bearer {api_token}"}

#%%

## Download precipitation data


page_token = None
page_size = 10

y = 0
while True:
    x=0
    parameters = {"pageSize": page_size}
    url1 = url
    
    if page_token == True:
        url1 = url + '/PageToken=' + next_page_token
        parameters['pageToken'] = next_page_token
        parameters = {"pageSize": page_size, 'nextPageToken': next_page_token}
        
    
    response0 = requests.get(url, headers=header, params=parameters).json()
    filenames = response0['files']
    
    for file in filenames:
        if filenames:
            # Choose the first filename in the list for download
            download_url = url + str(filenames[x]['filename']) + '/url'
            
            # Download the file and save it locally
            response1 = requests.get(download_url, headers=header, stream=True)
            if response1.status_code == 200:
                with open('knmi_data' + str(y) + str(x) + '.json', 'wb') as file:
                    for chunk in response1.iter_content(chunk_size=8192):
                        file.write(chunk)
                    
                #print("Data downloaded and saved")
            else:
                print(f"Failed to download data (Status code: {response1.status_code})")
        else:
            print("No filenames found for download.")
        infile = 'knmi_data' + str(y) + str(x) + '.json'
        with open(infile, 'r') as json_file:
            data_from_knmi = json.load(json_file)
            
        file_download_url = data_from_knmi['temporaryDownloadUrl']
    
        if 'file_download_url' in locals():
            # Fetch the SPI data from the temporary download URL
            response2 = requests.get(file_download_url)
    
            if response2.status_code == 200:
                # Successful response, save the data to a NetCDF file
                with open('knmi_precip_data' + str(y) + str(x) + '.nc', 'wb') as file:
                    for chunk in response2.iter_content(chunk_size=8192):
                        file.write(chunk)
    
                print("SPI data downloaded and saved ")
    
            else:
                print(f"Failed to download SPI data (Status code: {response2.status_code})")
        else:
            print("File download URL not found.")
        os.remove(infile)
        x=x+1
    
     # Check for a next page token
    if "nextPageToken" in response0:
        print('token!!')
        next_page_token = response0["nextPageToken"]
        page_token = True
    else:
        break
    y = y+1
