## Written by Declan Finney
## Designed for use on JASMIN computing facility, for which some packages would need to be installed in a virtual environment.
## Published Jun 2023

## DESCRIPTION
## This script was developed to calculate parallax-adjusted lat and lons of points in the GOES image for cloud properties
## such as cloud optical depth. Parallax adjustment is needed if lat lon precision on the scale of the cloud height is important.
## i.e. you want to know the clou dposition with a few km, but the cloud is over 10km high.
## The approach I have taken uses the cloud top height field from GOES (which is very coarse, but better than nothing).
## It then adds that to the minor and major axes radii of earth in the function to compute lat-lons.
## This means that an array of major and minor axe radii are passed through the function.
## The output array size is 2D, as it is if constants are used, but a different surface height from the earth centre
## has been assumed for each point, instead of the smoothly varying height of mean sea level.

## WARNING / REQUEST
## An observed clout-top height based adjustment for parallax shift is not something I've managed to find from someone else,
## which is why I'm sharing. However, there may be errors, despite my best efforts. Many functions, e.g. for regridding, are not
## set up for time varying lat and lons, and can be a bit awkward with curvilinear coordinates. I've had to work round this.
## Please let me know if you think there are errors.

## OUTPUT
## as well as resaving the read in files with updated lat and lons, a file with the dta regridded to a regular latlon grid
## is also saved.

## SYNTAX
## I use "< >" to mark paths etc that the user will need to set to their preference.

## EXAMPLE
## Example cloud optical depth and height files are included with this script on github. Those files
## were produced using the provided download_goes_subregion_regrid.py script.


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
from glob import glob
import xesmf as xe
import pandas as pd
from goes2go import GOES
import os
from datetime import datetime, timedelta

def goes_lonlat_parallax_corrected(goesdataset, heights_ASL):
    # heights_ASL will need to be same grid as goesdataset

    ## haven't got the goes_imager_projection variables in my raw saved data, so need to obtain them
    G = GOES(satellite=16, product="ABI-L2-ACHAC", domain='C')

    ## There are some sat files I have with incorrect t.values. I tried to edit but this still fails, I presume, because there's an error in the GOES archives files too!
    ## So in such cases I take a time 20mins earlier, which would hopefully have similar projections.
    try:
        satInfo = G.nearesttime(heights_ASL.t.values.astype(str))
    except ValueError:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("Couldn't get satInfo for this time, so using 20mins earlier")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        satInfo = G.nearesttime((heights_ASL.t.values - np.timedelta64(20,'m')).astype(str))
        
    proj_info = satInfo.variables['goes_imager_projection']
    lon_origin = proj_info.attrs["longitude_of_projection_origin"]
    ##dlf25apr2023 even though H includes semi_major_axis, I don't think I need to add cloud heights because 
    # i think perspective_point_height must be height above the surface, so we just need the radius to the surface,
    # to get full H
    H = proj_info.attrs["perspective_point_height"] +  proj_info.attrs["semi_major_axis"]


    ### dlf18apr2023 adjusting these for the cloud height interpolated to the goesdataset grid
    ## I believe its roughly implementing a similar approach to
    # eq15 in https://doi.org/10.3390/rs12030365
    r_eq = heights_ASL.fillna(0) + proj_info.attrs["semi_major_axis"]
    r_pol = heights_ASL.fillna(0) + proj_info.attrs["semi_minor_axis"]
    
    # Data info
    lat_rad_1d = goesdataset['x'][:]
    lon_rad_1d = goesdataset['y'][:]
    
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

    return lon, lat


##################
## USER TO EDIT
################
# I don't recommend applying this to the L1b Rad variables !
# since I'm only adjusting cloud for parallax, the surface wont be adjusted according to its height and therefore could screw up regridding.
# that's because land points may end up underneath cloud points, and then these would be averaged when regridding.
# in the listed variables above I think the clear-sky masking should prevent there being an issue
prods = ["ABI-L2-ACHAC","ABI-L2-CTPC","ABI-L2-ACHTF","ABI-L2-ACMC","ABI-L2-ACTPC"]
spacing=0.01 ## for regrid size. I figure ~1km spacing is going to be accurate as any parallax adjustment given the coarsness of GOES CTH
## i consider nearest neighbour to be most honest regirdder for this purpose but I found linear was ok.
## there are some strange edge artifacts outside the data subregion, so people will need to code a correction or use their intelligence
## when using the data.
rgdmethod="nearest_s2d"

goes_path = "<PATH_TO_GOES_FILES_TO_READ_IN>"
save_path = "<PATH_TO_SAVE_DATA_WITH_NEW_LATLONS>"
save_path_rgd = "<PATH_TO_SAVE_REGRIDDED_DATA_WITH_NEW_LATLONS>"

for prodnm in prods:
    ### list all files to reproduce
    ## Will correct latlons of every file I have saved for that prod.
    files_orig = glob(goes_path+prodnm+"/*/*/*/*/*")

    for fn in files_orig:
        print("resaving", fn)

        prod = xr.open_dataset(fn)
        ## strip the file name for the non-product specific bit
        prod_subdir_path = fn.split(prodnm)[1].split("OR_")[0]
        if prodnm == "ABI-L2-ACHAC":
            cth = prod
        else:
            ## USER may need to adapt this file name to their own format
            fn_ending = fn.split(prodnm)[2].split("_c")[0]+"_c*_select.nc"
            fn_cth = glob(goes_path+"ABI-L2-ACHAC"+prod_subdir_path+"OR_ABI-L2-ACHAC"+fn_ending)
            ## if cth file doesn't exist then need to skip this file
            ## e.g. 2022/07/16/16:11 (at least in my archive)
            if len(fn_cth) >0:
                cth = xr.open_dataset(fn_cth[0])
            else:
                print("============================")
                print("equivalent cth file doesn't exist so skipping this instance")
                print("=============================")
                continue

        ## Check datetime in file matches file name. There's at least one file that doesn't
        ## if it doesn't match then replace datetime with the midpoint of dates from the filename
        fn_startdatetime = fn.split("_s")[1].split("_e")[0]
        fn_enddatetime  = fn.split("_s")[1].split("_e")[1].split("_c")[0]
        year = int(fn_startdatetime[:4])
        doy = int(fn_startdatetime[4:7])
        hour = int(fn_startdatetime[7:9])
        minute = int(fn_startdatetime[9:11])
        second = int(fn_startdatetime[11:13])
        hrdiff = int(fn_enddatetime[7:9]) - int(fn_startdatetime[7:9])
        mindiff = int(fn_enddatetime[9:11]) - int(fn_startdatetime[9:11])
        secdiff = int(fn_enddatetime[11:13]) - int(fn_startdatetime[11:13])
        middatetime = datetime(year=year, month=1, day=1, hour=hour, minute=minute, second=second) + timedelta(days=doy - 1) + 0.5*timedelta(hours=hrdiff,minutes=mindiff, seconds=secdiff)
        # compare to the nearest minute; wont check seconds as they can be slightly different
        if prod.t.values.astype('datetime64[m]') != np.datetime64(middatetime).astype('datetime64[m]'):
            prod.t.values = np.datetime64(middatetime)
            cth.t.values = np.datetime64(middatetime)
            print("no")



        # get correct lat lons for cth
        cth_lon_original = cth["lon"] ## keep to compare
        cth_lat_original = cth["lat"]
        cth["lon"], cth["lat"] = goes_lonlat_parallax_corrected(cth.drop(["lat","lon"]).var1, cth.drop(["lat","lon"]).var1)

        prod_lon_original = prod["lon"] ## keep to compare
        prod_lat_original = prod["lat"]
        ## regrid cth to prod grid, just setting x and y to be lon and lat for the regridder because it only works if its lon and lat.
        # but I can't pass the cth lon and lat because they have the time dimension and the regridder wont accept that.
        ds_in = xr.Dataset({ "lat": (["y", "x"],cth_lat_original.values ), "lon": (["y", "x"], cth_lon_original.values )})
        ds_out = xr.Dataset({ "lat": (["y", "x"],prod_lat_original.values ), "lon": (["y", "x"], prod_lon_original.values )})
        regridder_cth_prod = xe.Regridder(ds_in,ds_out, method="bilinear")
        cth_prodrgd = regridder_cth_prod(cth.fillna(0))
        ## for some reason this leaves x and y but doesn't pass the prod x and y values so correcting.
        cth_prodrgd["x"], cth_prodrgd["y"] = prod["x"], prod["y"]

        ## calculate correct lat and lons for prod
        prod["lon"], prod["lat"] = goes_lonlat_parallax_corrected(prod.drop(["lat","lon"]).var1, cth_prodrgd.var1)

        # for the COD product also save the regridded cth with it so that COD can be related to the height.
        if prodnm == "ABI-L2-CODC":
            prod["cth_rgd"] = cth_prodrgd.var1

        # first make directory
        outdir = save_path+prodnm+prod_subdir_path
        outfilen = outdir+"OR_"+fn.split("OR_")[1].split("_select.nc")[0]+"_select_pc.nc"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        ## then save the netcdf    
        prod.to_netcdf(outfilen)

        print("regridding", fn)

        ### USER WILL NEED TO ADAPT THESE HARD-CODED lats and lons to ones that suit their domain.
        # forcing this with COD values because CTH seem to vary, which is awkward
        lat1 =32.8967170715332 ;lat2 = 35.17671707153275;lon1 =-109.12045288085938 ;lon2 =-105.0604528808573
        

        ## design the output regular latlon grid for regridding
        lat_rgd_vals = np.arange(lat1, lat2+spacing, spacing)
        lon_rgd_vals = np.arange(lon1, lon2+spacing, spacing)
        ds_out = xr.Dataset({ "lat": (["lat"], lat_rgd_vals), "lon": (["lon"], lon_rgd_vals)})

        regridder_prod = xe.Regridder(prod, ds_out, method=rgdmethod)

        cth_prodrgd_rgd = regridder_prod(cth_prodrgd)
        ## where cth is zero, the data is not reliable as it has not been adjusted for surface altitude parallax so dropping it.
        prod_rgd = regridder_prod(prod).where(cth_prodrgd_rgd>0)

        # first make directory
        outdir = save_path_rgd+prodnm+prod_subdir_path
        outfilen = outdir+"OR_"+fn.split("OR_")[1].split("_select.nc")[0]+"_select_pcrgd.nc"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        ## then save the netcdf    
        prod_rgd.to_netcdf(outfilen)

