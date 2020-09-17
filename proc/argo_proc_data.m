%// argo_proc_data.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to process Argo profiles that have passed quality control
%// The main steps in processing are:
%// 1. Remove duplicate pressure entries
%// 2. Copy temp/salinity profiles, smooth both using ksr filter
%// 3. Use smoothed T/S to calculate N2, normal T/S to get spice,sigma0
%// 4. Interpolate all data to 10dbar grid
%// 5. Check that density is increasing w/ depth (if not, sort and check for large differences)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%// Addpath to gsw_toolbox scripts
addpath scripts/

%// Clear unwanted variables
clearvars -except argo flag steps ARGO
close all
load('../datadir.mat');

%// Defaults
depth  = 10; %// Used in ksr filter
dbar2  = [0:2:2000]';  %// 2-dbar grid
dbar10 = [0:10:2000]'; %// 10-dbar grid

%// Initiate flag
flag_idx = zeros(1,length(argo));

%// Initiate sigma0, N2, spice
for i = 1:length(argo)
    argo(i).sigma0 = [];
    argo(i).N2     = [];
    argo(i).spice  = [];
end

%// Initiate parpool
%// Will take a long time to run without parallel processing...
delete(gcp('nocreate'))
parpool(10)

%// Track progress
progress = [0:10000:length(argo)];
parfor i = 1:length(argo) 
try
    ind = [];
    raw = argo(i);

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp(i)
    end

    %// Remove duplicate pressure entries (rare)
    [a,b]        = unique(argo(i).pres);
    argo(i).pres = argo(i).pres(b);
    argo(i).temp = argo(i).temp(b);
    argo(i).salt = argo(i).salt(b);
    a = []; b = [];

    %// Remove NaNs from P,T,S where they appear
    dat          = [argo(i).pres + argo(i).temp + argo(i).salt];
    argo(i).temp = argo(i).temp(~isnan(dat));
    argo(i).salt = argo(i).salt(~isnan(dat));
    argo(i).pres = argo(i).pres(~isnan(dat));
    dat = [];

    %// Interpolate temporary salinity/temperature to 2dbar grid
    smooth_temp = argo(i).temp;
    smooth_salt = argo(i).salt;
    smooth_temp = interp1(argo(i).pres,smooth_temp,dbar2);
    smooth_salt = interp1(argo(i).pres,smooth_salt,dbar2);

    %// Remove NaNs, apply smoothing filter
    smooth_pres = dbar2(~isnan(smooth_temp));
    smooth_temp = smooth_temp(~isnan(smooth_temp));
    smooth_salt = smooth_salt(~isnan(smooth_salt));
    rt          = ksr(smooth_pres,smooth_temp,depth,length(smooth_temp));
    rs          = ksr(smooth_pres,smooth_salt,depth,length(smooth_salt));
    
    %// Interpolate all data to 10dbar grid (including smoothed T/S)
    argo(i).temp = interp1(argo(i).pres,argo(i).temp,dbar10);
    argo(i).salt = interp1(argo(i).pres,argo(i).salt,dbar10);
    argo(i).pres = dbar10;
    smooth_temp  = interp1(rt.x,rt.f,dbar10);
    smooth_salt  = interp1(rs.x,rs.f,dbar10);
    rs = []; rt = [];

    %// Add absolute salinity and conservative temperature (also do this for smoothed T/S)
    dat1            = [argo(i).pres + argo(i).temp + argo(i).salt];
    dat2            = [argo(i).pres + smooth_temp + smooth_salt];
    salt_abs        = gsw_SA_from_SP(argo(i).salt(~isnan(dat1)),argo(i).pres(~isnan(dat1)),argo(i).lon,argo(i).lat);
    smooth_salt_abs = gsw_SA_from_SP(smooth_salt(~isnan(dat2)),argo(i).pres(~isnan(dat2)),argo(i).lon,argo(i).lat);
    theta           = gsw_CT_from_t(salt_abs,argo(i).temp(~isnan(dat1)),argo(i).pres(~isnan(dat1)));
    smooth_theta    = gsw_CT_from_t(smooth_salt_abs,smooth_temp(~isnan(dat2)),argo(i).pres(~isnan(dat2)));

    %// Add density and spiciness
    argo(i).sigma0 = gsw_sigma0(salt_abs,theta);
    argo(i).spice  = gsw_spiciness0(salt_abs,theta);

    %// Get N2 from smoothed data
    [N2_mid,pres_mid] = gsw_Nsquared(smooth_salt_abs,smooth_theta,argo(i).pres(~isnan(dat2)),argo(i).lat);

    %// Interpolate new variables to 10-dbar
    ndat           = isnan(N2_mid);
    argo(i).N2     = interp1(pres_mid(ndat==0),N2_mid(ndat==0),argo(i).pres);
    argo(i).spice  = interp1(argo(i).pres(~isnan(dat1)),argo(i).spice,argo(i).pres);
    argo(i).sigma0 = interp1(argo(i).pres(~isnan(dat1)),argo(i).sigma0,argo(i).pres);
    smooth_temp = []; smooth_salt = []; smooth_theta = []; smooth_salt_abs = []; N2_mid = []; pres_mid = [];

    %// Check that density is ascending with depth
    %// For profiles that aren't, sort the data and see if there are very small differences
    %// If profiles show difference > 0.1 kg/m3, then reject them (probably spikes in data)
    %// This routine catches very small density differences, usually right at the surface (0 - 100m)
    if issorted(argo(i).sigma0(~isnan(argo(i).sigma0)))==0
	sigma0_orig = argo(i).sigma0(~isnan(argo(i).sigma0));
	sigma0_sort = sort(sigma0_orig);
	if max(abs(sigma0_orig-sigma0_sort))>0.1
	    flag_idx(i) = 1; continue
	else
	    ind = find(isnan(argo(i).sigma0)==0);
	    argo(i).sigma0(ind) = sort(argo(i).sigma0(ind));
	end
    end

    %// Only keep levels with all data
    dat = [argo(i).sigma0 + argo(i).temp + argo(i).salt + argo(i).spice + argo(i).N2];
    argo(i).temp(isnan(dat)==1)   = NaN;
    argo(i).salt(isnan(dat)==1)   = NaN;
    argo(i).pres(isnan(dat)==1)   = NaN;
    argo(i).spice(isnan(dat)==1)  = NaN;
    argo(i).N2(isnan(dat)==1)     = NaN;
    argo(i).sigma0(isnan(dat)==1) = NaN;
catch
flag_idx(i) = 1;
end
end

%// Kill parpool
delete(gcp('nocreate'))

%// Add new flags to flag.total
ID               = [argo.ID];
ind              = find(flag_idx == 1);
flag.bad_density = ID(ind);
flag.total       = [flag.total flag.bad_density];
argo(ind)        = [];  
disp(['Removing ',num2str(length(ind)),' profiles due to bad density'])

%// Save data
fname = [datadir,'PROC_global_argo.mat'];
save(fname,'argo','flag','-v7.3')

clearvars -except argo steps flag 
