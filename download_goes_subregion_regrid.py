"""
GOES Satellite Data Download and Regridding Script

Written by Declan Finney
Designed for use on JASMIN computing facility, for which some packages would need to be installed in a virtual environment.
Published Jun 2023

DESCRIPTION:
This script downloads a subregion relevant for the DCMEX campaign for all of July-August 2022
when the campaign occurred. It can be used to download a variety of GOES products, but doesn't work with
all those available through goes2go. See some products below and any notes regarding them not working.
The script has not been developed explicitly for use by others so may not generalize well.
However, it is shared as an example code which may make it easier to code up similar tasks.
It also provides the base script for producing the netcdf that is then used with the parallax correction code.

IMPORTANT:
Lat/lons calculated in this script assume the image is of mean sea level height. There is NO parallax adjustment.
See parallax_latlons_goes_cloud.py for a way to make that adjustment.

USAGE:
python download_goes_subregion_regrid.py [product] [month] [channel]

Examples:
- python download_goes_subregion_regrid.py
- python download_goes_subregion_regrid.py ABI-L2-ACHAC 7 2

SYNTAX NOTES:
Comments use "< >" to mark paths etc that the user will need to set to their preference.
"""

## Written by Declan Finney
## Designed for use on JASMIN computing facility, for which some packages would need to be installed in a virtual environment.
## Published Jun 2023

## some interesting products. (there are more) https://github.com/blaylockbk/goes2go/blob/main/goes2go/product_table.txt
#ABI-L1b-RadC,   Advanced Baseline Imager Level 1b CONUS
#ABI-L2-ACHTF,   Advanced Baseline Imager Level 2 Cloud Top Temperature Full Disk (there isn't a CONUS version, unfortunately)
#ABI-L2-ACMC,    Advanced Baseline Imager Level 2 Clear Sky Mask CONUS
#ABI-L2-AODC,    Advanced Baseline Imager Level 2 Aerosol Optical Depth CONUS
#ABI-L2-DSRC,    Advanced Baseline Imager Level 2 Downward Shortwave Radiation CONUS ...(at surface!)
#ABI-L2-RSRC,    Advanced Baseline Imager Level 2 Reflected Shortwave Radiation Top-Of-Atmosphere CONUS
#GLM-L2-LCFA,    Geostationary Lightning Mapper Level 2 Lightning Detection
#ABI-L2-CODC,    Advanced Baseline Imager Level 2 Cloud Optical Depth CONUS
#ABI-L2-CPSC,    Advanced Baseline Imager Level 2 Cloud Particle Size CONUS
#ABI-L2-ACTPC,   Advanced Baseline Imager Level 2 Cloud Top Phase CONUS
#ABI-L2-ACHAC,   Advanced Baseline Imager Level 2 Cloud Top Height CONUS
#ABI-L2-RRQPEF,  Advanced Baseline Imager Level 2 Rainfall Rate (Quantitative Precipitation Estimate) Full Disk
#ABI-L2-DSIC,    Advanced Baseline Imager Level 2 Derived Stability Indices CONUS # not got this working yet
#ABI-L2-CTPC,    Advanced Baseline Imager Level 2 Cloud Top Pressure CONUS
#ABI-L2-MCMIPC,  Advanced Baseline Imager Level 2 Cloud and Moisture Imagery CONUS ## dlf27oct2022 cant currently use this in script

from goes2go import GOES 
import numpy as np
import sys
import xarray as xr
import os
import shutil
import xesmf as xe

# Check for required dependencies
try:
    from goes2go import GOES
except ImportError:
    print("Error: goes2go package not found. Please install with: pip install goes2go")
    sys.exit(1)

try:
    import xesmf as xe
except ImportError:
    print("Error: xesmf package not found. Please install with: pip install xesmf")
    print("Note: xesmf may require additional system libraries (ESMF)")
    sys.exit(1)

## function to return a string with 0 preceeding number, if number <10
## This is not necessary. e.g. Can replace with str(<x>).zfill(2) for a 2 digit number for example.
def str0(x):
    if x<10:
        xstr = "0"+str(x)
    else:
        xstr = str(x)
    return xstr

###################################
### USER TO CHOOSE OPTIONS BELOW
### it is possible to pass prods, mons, channels as terminal arguments
###################################
# SATELLITE
sat = 16
# REGION
regionName="Magda"; lat1=33.0;lat2=35.0;lon1=-108.0;lon2=-106.0
# which products
prods = ["ABI-L2-CTPC","ABI-L2-DSIC","ABI-L2-ACHAC","ABI-L2-ACHTF"]
if len(sys.argv)>1:
    prods = [str(sys.argv[1])]
    print(prods)
## time period
yr = 2022
mons = np.array([7,8])
if len(sys.argv)>2:
    mons = [int(sys.argv[2])]
    print(mons)
days_mon = np.array([np.arange(31)+1,np.arange(31)+1]) ## this needs to be a list of days, with a list for each month in the array above. for all days in july and august, e.g. np.array([np.arange(31)+1,np.arange(31)+1])
# if "ABI-L1b-RadC" listed then need to choose channels. wont' be used for other prods
## I RECOMMEND ONLY RUNNING ONE CHANNEL AT A TIME BECAUSE THE REGRIDDING CODE IN THIS SCRIPT WILL NOT BE EFFICIENT WITH MORE
channels = [2]
if len(sys.argv)>3:
    channels = [int(sys.argv[3])]
    print(channels)

# currently downloading all hours in day
##ARGS..
make_regrid_ctl = True # do you want a regridded version of data
only_keep_regrid_ctl = False # do you only want the regridded version
verbose=False # will output some information as script runs
download_data= True # knows not to download if already have files anyway.
only_download_data = False # if you don't want to make any files, maybe? [old option]
remove_all_GOES_data_upon_completion = True # this removes the raw downloaded data from GOES. Recommended.
#####################################
#####################################

 

##loop through prods, mons, days, hours
for prod in prods:
    if make_regrid_ctl:
        make_regrid=True # some variables will turn off regridding in the following loops because it is not sensible for them. but will want to turn back on for other variables.
    else:
        make_regrid=False
    if only_keep_regrid_ctl:
        only_keep_regrid=True # some variables will turn off only saving regrid because they dont need regridding and will want to save the raw data
    else:
        only_keep_regrid=False
    for imon,mon in enumerate(mons):
        tmppath="/gws/nopw/j04/dcmex/users/dfinney/data/tmpdata_GOES/"+prod+"_"+str(mon)+"_"+str(channels[0])+"/"
        ## make sure tmppath exists
        if not os.path.exists(tmppath):
            os.makedirs(tmppath)
        n_vars=1 ## set this as default each time start a new prod. Will be modified if a product takes more than 1 variable.
        firstfile_for_prod = True
        ## open a log file to track files that failed
        logf_name = "<PATH_TO_STORE_LOGFILE>/GOES"+str(sat)+"/GOES"+str(sat)+regionName+"_"+prod+"_"+str(mon)+"_"+str(channels[0])+"log_failedfiles.txt"
        logf = open(logf_name, "w")

        for day in days_mon[imon]:
            for hr in np.arange(24):                  
                startdatetime =str(yr)+"-"+str0(mon)+"-"+str0(day)+' '+str0(hr)+":00"
                enddatetime = str(yr)+"-"+str0(mon)+"-"+str0(day)+' '+str0(hr)+":59"
                print("~~~~~~~~~~~~~~~~~~~~~~")
                print(sat,prod,startdatetime,enddatetime)

                ## EXAMPLE ##
                # First, create a GOES object to specify the satellite, data product, and domain
                if prod=="ABI-L1b-RadC":
                    G = GOES(satellite=sat, product=prod, domain='C', bands=channels)
                else:
                    G = GOES(satellite=sat, product=prod, domain='C')
                # and make list of files in time range without downloading. Try it, because sometimes the files might not exist (e.g. didn't get any for ABI-L2-ACHTF 2022/07/20 20:00-20:59
                try:
                    ds = G.df(start=startdatetime, end=enddatetime)
                except TypeError:
                    print("couldn't find files for some reason")
                    logf.write(str("no files for "+prod+" "+startdatetime+".\n"))
                    continue # just move to next hour
                if verbose:
                    print(ds[0:5])
                    print(ds[-5:])

                # Download data for a specified time range. GOES to a folder.
                if download_data:
                    try:
                        G.timerange(start=startdatetime, end=enddatetime,save_dir=tmppath)
                    except ValueError:
                        print("Couldn't download files for some reason for "+startdatetime+".\n")
                        logf.write(str("couldnt download for "+prod+" "+startdatetime))
                        continue # just move to next hour
                    if only_download_data:
                        continue
                else:
                    continue

                ### MAKE AND SAVE SUBSELECTED REGION
                ### loop though each file
                # ds contains all filenames so can go through each one loading processing and deleting.
                for ifile,filen in enumerate(ds["file"]):
                    if verbose:
                        print("--------------------------------------")
                        print(ifile, ifile/ds["file"].shape[0], filen)
                    fullfname = tmppath+ds.iloc[ifile]["file"]
                    ## try to open netcdf. if can't then skip to next file.
                    try:
                        ds_disk = xr.open_dataset(fullfname)
                    except ValueError:
                        # seems to be a problem with some files in the time variable. e.g. /home/users/dfinney/data/noaa-goes16/ABI-L2-AODC/2022/190/11/OR_ABI-L2-AODC-M4_G16_s20221901100220_e20221901105120_c20221901108249.nc
                        logf.write("ValueError for "+fullfname+"\n")
                        ds_disk = xr.open_dataset(fullfname, decode_times=False)
                    except OSError:
                        print(fullfname,"didn't open")
                        logf.write("OSError for "+fullfname+".\n")
                        continue                           
 
                    
                    variable_names = list(ds_disk.data_vars)
                    var1 = ds_disk[variable_names[0]] #dlf14oct2022 I'm hoping the first variable is normally the one I want 
                    
                    if prod=="ABI-L2-ACMC": ## then want the binary mask and the multi-level mask
                        var2 = ds_disk[variable_names[1]]
                        var3 = ds_disk[variable_names[2]]
                        n_vars = 3 ## set this to make a more general approach for the remainder of code

                    #if a different channel or something then need to redo latlon stuff
                    coord_list = list(ds_disk.coords)
                    if "lon" and "lat" in coord_list: # some prods already have lats and lons (I don't know if they're parallax corrected)
                        var_select = var1.sel(lat=slice(lat2,lat1), lon=slice(lon1,lon2))
                        flag_select = ds_disk["DQF"].sel(lat=slice(lat2,lat1), lon=slice(lon1,lon2))
                        ## combine vars into a dataset
                        ds_select = xr.Dataset(data_vars=dict(var1=var_select,flag=flag_select))
                        if make_regrid:
                            ## don't think there's any need to regrid because already on a latlon grid
                            make_regrid=False
                            ## also want to turn this off for this product so that we save the raw version
                            only_keep_regrid=False
                            
                    elif "flash_lon" and "flash_lat" in coord_list: # glm data
                        ## lat and lons are parallax corrected as described by https://doi.org/10.1175/JTECH-D-19-0100.1
                        ## I believe it just assumes a fixed cloud height in time, but which varies and is optimised for each location. I may be wrong.
                        ## this data is just a list so needs to be treated differently
                        ## dlf7nov2022 first check there's flashes in this file (which seemed to be come an issue from 2022-08-15 21:30Z
                        if ds_disk.flash_lat.shape[0]==0:
                            logf.write("No lightning flashes in "+fullfname+"\n")
                            continue
                        index = np.where((ds_disk.flash_lat>=lat1) & (ds_disk.flash_lat <=lat2) & (ds_disk.flash_lon>=lon1) & (ds_disk.flash_lon <=lon2))
                        if len(index[0]) >=1:
                            ds_select = ds_disk["flash_energy"][index[0]]
                            if make_regrid:
                                ## not bothering to regrid lightning because that is not on any sort of grid to begin with.
                                make_regrid=False
                            
                        else:
                            # don't try and save anything if there's no flashes in that time period
                            continue
                    else:
                        newsource = filen.rsplit('/',1)[1].rsplit('_')[0:3]
                        if ifile==0 or newsource!=source:
                            ###### This explains how to convert data to lat, lon coords
                            # https://makersportal.com/blog/2018/11/25/goes-r-satellite-latitude-and-longitude-grid-projection-algorithm\
                            ## It assumes the image is of mean sea level height. If you want to adjust each location by a fixed amount
                            ## then you could add an array to r_eq and r_pol.
                            # GOES-R projection info and retrieving relevant constants
                            proj_info = ds_disk.variables['goes_imager_projection']
                            lon_origin = proj_info.attrs["longitude_of_projection_origin"]
                            H = proj_info.attrs["perspective_point_height"] +  proj_info.attrs["semi_major_axis"]
                            r_eq = proj_info.attrs["semi_major_axis"]
                            r_pol = proj_info.attrs["semi_minor_axis"]

                            # Data info
                            lat_rad_1d = ds_disk.variables['x'][:]
                            lon_rad_1d = ds_disk.variables['y'][:]

                            # create meshgrid filled with radian angles
                            lat_rad,lon_rad = np.meshgrid(lat_rad_1d,lon_rad_1d)

                            # lat/lon calc routine from satellite radian angle vectors

                            lambda_0 = (lon_origin*np.pi)/180.0

                            a_var = np.power(np.sin(lat_rad),2.0) + (np.power(np.cos(lat_rad),2.0)*(np.power(np.cos(lon_rad),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(lon_rad),2.0))))
                            b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
                            c_var = (H**2.0)-(r_eq**2.0)

                            r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)

                            s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
                            s_y = - r_s*np.sin(lat_rad)
                            s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)

                            lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
                            lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
                            var_latlon = var1.assign_coords({"lat": (("y", "x"), lat), "lon": (("y", "x"), lon)})

                            ## find index for subselect region
                            index=np.where((var_latlon.lat.data>lat1) & (var_latlon.lat.data<lat2)
                                   & (var_latlon.lon.data>lon1) & (var_latlon.lon.data<lon2))
                            xy_minmax = [np.min(index[0]),np.max(index[0]), np.min(index[1]),np.max(index[1])]
                            

                            # set source to newsource
                            source= newsource

                            if firstfile_for_prod:
                                if make_regrid:  ## calculate the regrid grid
                                    ## use average lat and lon spacing
                                    dlat = np.abs(np.nanmean(var_latlon.lat.values[1:] - var_latlon.lat.values[:-1]))
                                    ##dlat seems much more robust, for ctt, got a dlon of e-10 so going to just set dlon to be same as dlat.
                                    #dlon = np.abs(np.nanmean(var_latlon.lon.values[1:] - var_latlon.lon.values[:-1]))
                                    dlon = dlat
                                    lat_rgd_vals = np.arange(lat1, lat2+dlat, dlat)
                                    lon_rgd_vals = np.arange(lon1, lon2+dlon, dlon)
                                    ds_out = xr.Dataset({"lat": (["lat"], lat_rgd_vals), "lon": (["lon"], lon_rgd_vals),})
                                    regridder = xe.Regridder(var_latlon[xy_minmax[0]:xy_minmax[1],xy_minmax[2]:xy_minmax[3]].transpose(), ds_out.transpose(), "bilinear")

                                    # when going through multiple channels there can be different grids so just going to remake regridder each time, but otherwise, on need to make the first time
                                    # dlf11nov2022 THIS IS PROBABLY VERY INEFFICIENT. Better to do one channel at a time, or rewrite this code.
                                    if (prod!="ABI-L1b-RadC") or (len(channels)<=1):
                                        firstfile_for_prod = False # set this to false, until start new prod because dont want a different regrid grid for every file.

                        else:
                            ## apply latlons to data and select region
                            var_latlon = var1.assign_coords({"lat": (("y", "x"), lat), "lon": (("y", "x"), lon)})

                        if verbose:
                            print("::::::::::::::")
                            print("source", source)
                        ## subselect region now made for whatever file has come through
                        var_select = var_latlon[xy_minmax[0]:xy_minmax[1],xy_minmax[2]:xy_minmax[3]]
                        ## also take the quality flags
                        flag_select = ds_disk["DQF"][xy_minmax[0]:xy_minmax[1],xy_minmax[2]:xy_minmax[3]]
                        ## if there's more vars...
                        if n_vars ==3:
                            v2_select = var2[xy_minmax[0]:xy_minmax[1],xy_minmax[2]:xy_minmax[3]]
                            v3_select =	var3[xy_minmax[0]:xy_minmax[1],xy_minmax[2]:xy_minmax[\
3]]

                        
                        ## combine vars into a dataset
                        if n_vars==3:
                            ds_select = xr.Dataset(data_vars=dict(var1=var_select,var2=v2_select, var3=v3_select,flag=flag_select))
                            if make_regrid:
                                var_select_rgd = regridder(var_select)
                                v2_select_rgd = regridder(v2_select)
                                v3_select_rgd = regridder(v3_select)
                                flag_select_rgd = regridder(flag_select)
                                ds_select_rgd = xr.Dataset(data_vars=dict(var1=var_select_rgd,var2=v2_select_rgd, var3=v3_select_rgd,flag=flag_select_rgd))
                        else:
                            ds_select = xr.Dataset(data_vars=dict(var1=var_select,flag=flag_select))
                            if make_regrid:
                                var_select_rgd = regridder(var_select)
                                flag_select_rgd = regridder(flag_select)
                                ds_select_rgd = xr.Dataset(data_vars=dict(var1=var_select_rgd,flag=flag_select_rgd))


                    # save file
                    if not only_keep_regrid:
                        outdir = "/gws/nopw/j04/dcmex/data/GOES"+str(sat)+"/"+regionName+'/'+prod+"/"+str(yr)+"/"+str(mon).zfill(2)+"/"+str(day).zfill(2)+"/"+str(hr).zfill(2)+"/"
                        outfilen = outdir+ds.iloc[ifile]["file"].rsplit('/',1)[1][:-3]+"_select.nc"
                        if not os.path.exists(outdir):
                            os.makedirs(outdir)
                        ## then save the netcdf    
                        ds_select.to_netcdf(outfilen)
                    if make_regrid:
                        outdir = "/gws/nopw/j04/dcmex/data/GOES"+str(sat)+"rgd/"+regionName+'/'+prod+"/"+str(yr)+"/"+str(mon).zfill(2)+"/"+str(day).zfill(2)+"/"+str(hr).zfill(2)+"/"
                        outfilen = outdir+ds.iloc[ifile]["file"].rsplit('/',1)[1][:-3]+"_select_rgd.nc"
                        if not os.path.exists(outdir):
                            os.makedirs(outdir)
                        ## then save the netcdf    
                        ds_select_rgd.to_netcdf(outfilen)


                # and delete raw downloaded file
                if remove_all_GOES_data_upon_completion: # try and delete both
                    directory = tmppath+"noaa-goes"+str(sat)
                    shutil.rmtree(directory,ignore_errors=True) #ignoring errors because if not it fails because of some .nfs00* files left. but I don't care if these are left around.


    ##close failed files log
    logf.close()

