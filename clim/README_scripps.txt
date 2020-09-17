To obtain Roemmich-Gilson Argo Climatology from Scripps
Go to:
http://sio-argo.ucsd.edu/RG_Climatology.html

Download the latest temperature and salinity climatology netcdf files
and place them in the clim/ folder (where this README exists).

If you want to add new months, they provide individual files. Place them in the monthly_adds/ folder.
Then see 'create_monthly_climatologies.m', adding new months is pretty simple.

Note that this will need to be updated, since the Scripps website is occassionally updated.
We downloaded the 2004 - 2017 version, but a 2004 - 2018 is available as of Sep 2020.

After the climatology is downloaded, you need to run 
'scripps_clim_pyc_density.m' 
to get mat file of mid-pycnocline depths, or else ../proc/argo_get_initial_scvs.m will fail
