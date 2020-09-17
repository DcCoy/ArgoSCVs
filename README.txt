This folder contains all the scripts you need to find SCVs in the Argo dataset.
argo_scvs_master.m is the master script, use it to perform each step of the processing.

Note that we set up these scripts using parallel processing. If you are unable to parallel process,
find all instances of the following and comment them.

delete(gcp('nocreate'))
parpool(...)
parfor i = ...

Then, replace the parfor loop with a standard for-loop
i.e.
parfor i = 1:length(argo_lon)
becomes 
for i = 1:length(argo_lon)

Before you begin, open argo_scvs_master.m and view the header.
Each file contains information at the top of the script, but a quick description of each file is given below.

####################################################################
####################################################################
##################################
argo_scvs_master.m
##################################
Master script, use it to process and detect SCVs in the Argo record.
Make sure you have the argo data and Scripps climatology downloaded before you begin.
See the header for more information.

##################################
argo_load_good_data.m
##################################
Loads Argo profiles from Pacific, Atlantic, and Indian basins. Loops through each day and loads the netcdf
file provided by the GDAC server. Only allows data that has been flagged as good or probably good at each
pressure level. Some data flags are applied which reject obviously bad T/S/P data.

##################################
argo_qc_data.m
##################################
Script which applies several quality control routines to the Argo record. Rejects profiles that have been
previously flagged by other users (grey-listed, file included), as well as throwing away profiles with poor data resolution
or shallow casts.

##################################
argo_proc_data.m
##################################
Script to process Argo profiles that have passed quality control. Interpolates all data to 10-dbar grid.
Also checks that density is always increasing, rejecting profiles which don't fit this criteria.
Also calls ksr.m, a Kernel smoothing regression technique. See that file for more information.

##################################
argo_get_anom.m
##################################
Script to match each profile with a climatological profile provided by the Roemmich-Gilson Argo Climatology.
Then, calculates anomalies along isopycnals for each cast.

##################################
argo_get_thresh.m
##################################
Calculates the interquartile range (IQR) for spiciness and buoyancy frequency anomalies along each isopycnal.
Also calls find_nearby_floats, which is a simple routine that looks for other nearby/neartime casts used
to calculate the IQR.

##################################
argo_get_initial_scvs.m
##################################
Script to detect instances of spiciness and buoyancy frequency anomaly exceeding interquartile range thresholds
along isopycnals. Keeps those that pass this test as potential SCVs.

##################################
argo_make_datafile.m
##################################
Script that organizes initial detections into an easier to use datafile. 

##################################
initial_QC.m
##################################
Script that takes initial detections and checks that T/S data are not completely outside of interquartile range
thresholds, which would indicated temperature or salinity sensor errors.

##################################
get_spicy_scvs.m / get_minty_scvs
##################################
Script that takes initial detections that have passed 'initial_QC.m' and performs several tasks:
1) Fit a gaussian model to spice anomalies in pressure space
2) Calculates dynamic height from the deepest pressure, also for the climatological cast.
3) Performs a vertical mode decomposition on the climatological cast.
4) Extracts 1st horizontal baroclinic mode of climatological cast.
5) Fits BC1 mode to dynamic height anomaly profile and removes it.
6) Checks that there is a peak in dynamic height anomaly within the limits of the SCV as estimated
   by the gaussian model.
Also calls fit_gaussian.m, which fits a 1-term gaussian model to spiciness anomaly profiles that have passed IQR 
thresholds. This script calls findpeaks2.m, which is similar to findpeaks.m but doesn't require a matlab toolbox.

##################################
clean_individual_scvs
##################################
Script that makes the spicy and minty scv datafile easier to read for other users.

##################################
spicy_scvs_TS / minty_scvs_TS
##################################
Script that identifies Argo floats which have detected SCVs over consecutive casts, essentially constructing
SCV time-series.

##################################
clean_timeseries_scvs
##################################
Script that makes the spicy and minty timeseries datafile easier to read for other users.

####################################################################
####################################################################

Reach out to Daniel McCoy if you have any questions
dmccoy801@gmail.com
