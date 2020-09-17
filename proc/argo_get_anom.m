%// argo_get_anom.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to calculate anomalies for each float 
%// against the Scripps gridded Argo climatology.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// NOTES
%// See ../clim/README_scripps.txt first!!! You must have the
%// climatology downloaded before you begin!

%// Addpath to gsw scripts
addpath scripts/

%// Clear unwanted variables
clearvars -except argo flag steps
close all
load('../datadir.mat');

%// Load Scripps climatology (updated to 2019)
disp('Loading SCRIPPS climatology')
load('../clim/RG_ArgoClim_TS_Climatologies.mat');
scripps = RG_ArgoClim_TS_Climatologies; clear RG_ArgoClim_TS_Climatologies;

%// Get Scripps lon and lat (fix lon);
disp('Grab lon/lat data')
sclim.lon = scripps.lon; sclim.lon(sclim.lon>360) = sclim.lon(sclim.lon>360)-360;
sclim.lat = scripps.lat;

%// Get Argo lon and lat (fix lon);
alon = [argo.lon]; alon(alon<0) = alon(alon<0)+360;
alat = [argo.lat];

%// Make matrix data (easier to vectorize)
disp('Get Argo CTD data')
temp   = [argo.temp];
salt   = [argo.salt];
pres   = [argo.pres]; %// Same for all floats
spice  = [argo.spice];
N2     = [argo.N2];
sigma0 = [argo.sigma0];

%// Start filling in climatological data for each float
disp('Make climatology structure for proc')
atime   = nan(1,length(argo));
lon_ind = nan(1,length(argo));
lat_ind = nan(1,length(argo));

%// Keep track of progress
progress = [0:100000:length(argo)];
for i = 1:length(argo)

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str((i/length(argo))*100),' %'])
    end

    %// Get closest longitude/latitude from Scripps
    atime(i)   = argo(i).time(2);
    lon_ind(i) = min(find(abs(sclim.lon-alon(i)) == min(abs(sclim.lon-alon(i)))));
    lat_ind(i) = min(find(abs(sclim.lat-alat(i)) == min(abs(sclim.lat-alat(i)))));

    %// Also build clim matrix for argo_add_properties
    clim(i).lon  = sclim.lon(lon_ind(i));
    clim(i).lat  = sclim.lat(lat_ind(i));
    clim(i).pres = scripps.pres;

    %// Check if temp/salinity data exists, if not get average from nearby grid cells
    check = squeeze(scripps.temp(lon_ind(i),lat_ind(i),:,atime(i)));
    if isnan(nanmean(check)) == 1 | length(check(~isnan(check))) < 11 %//less than 100m

	%// Temperature of nearest 8 cells
	try
	temp1 = squeeze(scripps.temp(lon_ind(i)-1,lat_ind(i),:,atime(i)));
	catch temp1 = nan(size(scripps.pres)); end
	try
	temp2 = squeeze(scripps.temp(lon_ind(i)+1,lat_ind(i),:,atime(i)));
	catch temp2 = nan(size(scripps.pres)); end
	try
	temp3 = squeeze(scripps.temp(lon_ind(i),lat_ind(i)-1,:,atime(i)));
	catch temp3 = nan(size(scripps.pres)); end
	try
	temp4 = squeeze(scripps.temp(lon_ind(i),lat_ind(i)+1,:,atime(i)));
	catch temp4 = nan(size(scripps.pres)); end
	try
	temp5 = squeeze(scripps.temp(lon_ind(i)+1,lat_ind(i)+1,:,atime(i)));
	catch temp5 = nan(size(scripps.pres)); end
	try
	temp6 = squeeze(scripps.temp(lon_ind(i)-1,lat_ind(i)+1,:,atime(i)));
	catch temp6 = nan(size(scripps.pres)); end
	try
	temp7 = squeeze(scripps.temp(lon_ind(i)-1,lat_ind(i)-1,:,atime(i)));
	catch temp7 = nan(size(scripps.pres)); end
	try
	temp8 = squeeze(scripps.temp(lon_ind(i)+1,lat_ind(i)-1,:,atime(i)));
	catch temp8 = nan(size(scripps.pres)); end

	%// Salinity of nearest 8 cells
	try
	salt1 = squeeze(scripps.salt(lon_ind(i)-1,lat_ind(i),:,atime(i)));
	catch salt1 = nan(size(scripps.pres)); end
	try
	salt2 = squeeze(scripps.salt(lon_ind(i)+1,lat_ind(i),:,atime(i)));
	catch salt2 = nan(size(scripps.pres)); end
	try
	salt3 = squeeze(scripps.salt(lon_ind(i),lat_ind(i)-1,:,atime(i)));
	catch salt3 = nan(size(scripps.pres)); end
	try
	salt4 = squeeze(scripps.salt(lon_ind(i),lat_ind(i)+1,:,atime(i)));
	catch salt4 = nan(size(scripps.pres)); end
	try
	salt5 = squeeze(scripps.salt(lon_ind(i)+1,lat_ind(i)+1,:,atime(i)));
	catch salt5 = nan(size(scripps.pres)); end
	try
	salt6 = squeeze(scripps.salt(lon_ind(i)-1,lat_ind(i)+1,:,atime(i)));
	catch salt6 = nan(size(scripps.pres)); end
	try
	salt7 = squeeze(scripps.salt(lon_ind(i)-1,lat_ind(i)-1,:,atime(i)));
	catch salt7 = nan(size(scripps.pres)); end
	try
	salt8 = squeeze(scripps.salt(lon_ind(i)+1,lat_ind(i)-1,:,atime(i)));
	catch salt8 = nan(size(scripps.pres)); end

	%// Average profiles from all nearby grid cells (with data)
	TEMP = [temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8];
	SALT = [salt1 salt2 salt3 salt4 salt5 salt6 salt7 salt8];
	for row = 1:length(TEMP)
	    nan_ind = find(isnan(TEMP(row,:))==0);
	    if length(nan_ind)>0
		ttemp(row) = sum(TEMP(row,nan_ind))/length(nan_ind);
		tsalt(row) = sum(SALT(row,nan_ind))/length(nan_ind);
	    else
		ttemp(row) = NaN; tsalt(row) = NaN;
	    end
	end
	ttemp = ttemp'; tsalt = tsalt';
	if isnan(nanmean(ttemp))==1 | isnan(nanmean(tsalt))==1 | length(ttemp(~isnan(ttemp))) < 11
	    ttemp = nan(size(scripps.pres));
	    tsalt = nan(size(scripps.pres));
	end
	clim(i).temp = ttemp;
	clim(i).salt = tsalt;

    else
	%// Data exists in region, grab T/S profile for that time of year
	clim(i).temp = squeeze(scripps.temp(lon_ind(i),lat_ind(i),:,atime(i)));
	clim(i).salt = squeeze(scripps.salt(lon_ind(i),lat_ind(i),:,atime(i)));
    end

    %// Fix orientation (if necessary)
    [a,b] = size(clim(i).temp);
    if b > a
	clim(i).temp = clim(i).temp';
	clim(i).salt = clim(i).salt';
    end
end

%// Establish temp,salt data from scripps clim for each float
disp('Get climatology CTD data')
slon  = [clim.lon]; slon(slon>360) = slon(slon>360)-360; %//switch back to normal lon
slat  = [clim.lat];
stemp = [clim.temp];
ssalt = [clim.salt];
spres = [clim.pres];
clear clim scripps sclim

%// Add other data
disp('Grab salt_abs (patience!)')
ssalt_abs = gsw_SA_from_SP(ssalt,spres,slon,slat);
disp('Grab theta')
stheta = gsw_CT_from_t(ssalt_abs,stemp,spres);
disp('Grab sigma0')
ssigma0 = gsw_sigma0(ssalt_abs,stheta);
disp('Grab spice')
sspice = gsw_spiciness0(ssalt_abs,stheta);
disp('Grab N2')
[sN2_mid,pres_mid] = gsw_Nsquared(ssalt_abs,stheta,spres,slat);
sdat = isnan(sN2_mid);

%// Interpolate N2 back to regular pressure
disp('Interpolate N2 back to regular pressure grid')
progress = [0:100000:length(argo)];
for i = 1:length(argo)

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str((i/length(argo))*100),' %'])
    end

    %// Try/catch will reject bad profiles
    try
	sN2(:,i) = interp1(pres_mid(sdat(:,i)==0,i),sN2_mid(sdat(:,i)==0,i),spres(:,1));
    catch
	sN2(:,i) = nan(58,1);
    end
end
clear sN2_mid pres_mid

disp('Finding anomalies along density surfaces')
%// Initiate bad climatology flag
flag.climatology = [];
flag_idx         = [];

%// Keep track of progress
progress = [0:100000:length(argo)];
for i = 1:length(argo)

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str((i/length(argo))*100),' %'])
    end

    %// Check that climatology exists
    if isnan(nanmean(stemp(:,i)))==1 | isnan(nanmean(ssalt(:,i)))==1
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue;
    end

    %// Check for inverted climatological density (at surface)
    clim_sigma0 = ssigma0(:,i);
    if issorted(clim_sigma0(~isnan(clim_sigma0)))==0
	ind = find(isnan(ssigma0(:,i))==0);
	ssigma0(ind,i) = sort(ssigma0(ind,i));
    end

    %// Find where argo(i).sigma0 is available
    aind         = find(isnan(argo(i).sigma0)==0);
    
    %// Find where all rows have data
    dat          = [ssigma0(:,i) + stemp(:,i) + ssalt(:,i) + sspice(:,i) + sN2(:,i)];

    %// Interpolate climatology to float sigma0 grid
    clim_sigma0  = ssigma0(:,i);
    clim_temp    = stemp(:,i);
    clim_temp    = interp1(clim_sigma0(~isnan(dat)),clim_temp(~isnan(dat)),argo(i).sigma0(aind));
    filler       = nan(length(argo(i).sigma0),1);
    filler(aind) = clim_temp; clim_temp = filler;
    clim_salt    = ssalt(:,i);
    clim_salt    = interp1(clim_sigma0(~isnan(dat)),clim_salt(~isnan(dat)),argo(i).sigma0(aind));
    filler       = nan(length(argo(i).sigma0),1);
    filler(aind) = clim_salt; clim_salt = filler;
    clim_spice   = sspice(:,i);
    clim_spice   = interp1(clim_sigma0(~isnan(dat)),clim_spice(~isnan(dat)),argo(i).sigma0(aind));
    filler       = nan(length(argo(i).sigma0),1);
    filler(aind) = clim_spice; clim_spice = filler;
    clim_N2      = sN2(:,i);
    clim_N2      = interp1(clim_sigma0(~isnan(dat)),clim_N2(~isnan(dat)),argo(i).sigma0(aind));
    filler       = nan(length(argo(i).sigma0),1);
    filler(aind) = clim_N2; clim_N2 = filler;
    clim_pres    = spres(:,i);
    clim_pres    = interp1(clim_sigma0(~isnan(dat)),clim_pres(~isnan(dat)),argo(i).sigma0(aind));
    filler       = nan(length(argo(i).sigma0),1);
    filler(aind) = clim_pres; clim_pres = filler;

    %// Build structure
    argo_anom(i).float      = argo(i).float;
    argo_anom(i).cycle      = argo(i).cycle;
    argo_anom(i).ID         = argo(i).ID;
    argo_anom(i).lon        = argo(i).lon;
    argo_anom(i).lat        = argo(i).lat;
    argo_anom(i).time       = argo(i).time;
    argo_anom(i).temp_anom  = argo(i).temp - clim_temp;
    argo_anom(i).salt_anom  = argo(i).salt - clim_salt;
    argo_anom(i).spice_anom = argo(i).spice - clim_spice;
    argo_anom(i).N2_anom    = argo(i).N2 - clim_N2;
    argo_anom(i).pres_anom  = argo(i).pres - clim_pres;
    argo_anom(i).sigma0     = argo(i).sigma0;

    %// Check for empty data, add to flags
    if length(argo_anom(i).temp_anom(~isnan(argo_anom(i).temp_anom))) == 0
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue
    elseif length(argo_anom(i).salt_anom(~isnan(argo_anom(i).salt_anom))) == 0
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue
    elseif length(argo_anom(i).spice_anom(~isnan(argo_anom(i).spice_anom))) == 0
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue
    elseif length(argo_anom(i).N2_anom(~isnan(argo_anom(i).N2_anom))) == 0
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue
    elseif length(argo_anom(i).pres_anom(~isnan(argo_anom(i).pres_anom))) == 0
	flag.climatology = [flag.climatology argo(i).ID]; flag_idx = [flag_idx i]; continue
    end
end

%// Save flags of floats with no climatology (for QC figures)
disp(['Removing ',num2str(length(flag.climatology)),' profiles due to no climatology'])
flag.total          = [flag.total flag.climatology];
argo_anom(flag_idx) = [];
argo(flag_idx)      = [];

%// Save data
fname = [datadir,'ANOM_global_argo.mat'];
save(fname,'argo_anom','flag','-v7.3')

clearvars -except argo argo_anom flag steps



