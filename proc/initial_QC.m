%// Script to remove bad SCV detections using the below criteria:
%// 1. Minimum pressure < 100 and Maximum pressure > 700
%// 2. T/S/Spice profile within 1.5*IQR at some point
%// 3. IQR exceedence must be localized with depth (not just at sfc)

disp('Checking detections for bad temp/salt anomaly or pressure values')
disp('Adding climatological reference cast to structure')
load('../datadir.mat');

% Addpath to gsw scripts
addpath scripts/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Load initial scv dataset
disp('Loading initial SCV data')
load([datadir,'global_initial_scv_data.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Load Scripps climatology
disp('Loading Scripps Climatology')
load('../clim/RG_ArgoClim_TS_Climatologies.mat');
scripps = RG_ArgoClim_TS_Climatologies; clear RG_ArgoClim_TS_Climatologies;

%// Get scripps lon/lat, fix lon
sclim.lat                = scripps.lat;
sclim.lon                = scripps.lon;
sclim.lon(sclim.lon>180) = sclim.lon(sclim.lon>180)-360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Get IDs so we can reassemble datasets later, only keep unique entries
mID       = [minty_scv.ID];
sID       = [spicy_scv.ID];

%// Group all SCV data together (for better IQR results)
scv_data = [spicy_scv minty_scv];

%// Establish flags
bad_temp_prof  = [];
bad_salt_prof  = [];
bad_spice_prof = [];
bad_spice_anom = [];
bad_reference  = [];

%// Scroll through to add reference profile to structure and flag bad profiles
for i = 1:length(scv_data)

    %// Track progress
    progress = [0:5000:length(scv_data)];
    if ismember(i,progress)==1
	disp([num2str((i/length(scv_data))*100),' percent complete'])
    end

    %// Rebuild climatological profile  
    stime                = scv_data(i).time(2);
    x                    = abs(sclim.lon-scv_data(i).lon);
    y                    = abs(sclim.lat-scv_data(i).lat);
    lon_ind              = min(find(x==min(x)));
    lat_ind              = min(find(y==min(y)));
    scv_data(i).ref.lon  = sclim.lon(lon_ind);
    scv_data(i).ref.lat  = sclim.lat(lat_ind);
    scv_data(i).ref.pres = scripps.pres;
    scv_data(i).ref.temp = squeeze(scripps.temp(lon_ind,lat_ind,:,stime));
    scv_data(i).ref.salt = squeeze(scripps.salt(lon_ind,lat_ind,:,stime));

    %// Flag shallow reference casts
    if max(scv_data(i).ref.pres(~isnan(scv_data(i).ref.temp))) < 700
	bad_reference = [bad_reference i];
	continue;
    elseif isempty(scv_data(i).ref.pres(~isnan(scv_data(i).ref.temp))) == 1
	bad_reference = [bad_reference i];
	continue;
    end

%// Add zero at surface (for VMD), use first t/s data
    scv_data(i).ref.pres = [0;scripps.pres];
    t                    = scv_data(i).ref.temp(~isnan(scv_data(i).ref.temp)); t = t(1);
    s                    = scv_data(i).ref.salt(~isnan(scv_data(i).ref.salt)); s = s(1);
    scv_data(i).ref.temp = [t; scv_data(i).ref.temp];
    scv_data(i).ref.salt = [s; scv_data(i).ref.salt];   

%// Interpolate to 10dbar grid (0 - 2000 dbar)
    pgrid                = [0:10:2000]'; %'
    scv_data(i).ref.temp = interp1(scv_data(i).ref.pres(~isnan(scv_data(i).ref.temp)),scv_data(i).ref.temp(~isnan(scv_data(i).ref.temp)),pgrid);
    scv_data(i).ref.salt = interp1(scv_data(i).ref.pres(~isnan(scv_data(i).ref.salt)),scv_data(i).ref.salt(~isnan(scv_data(i).ref.salt)),pgrid);
    scv_data(i).ref.pres = pgrid;
	
%// Add other properties
    scv_data(i).ref.salt_abs = gsw_SA_from_SP(scv_data(i).ref.salt,scv_data(i).ref.pres,scv_data(i).ref.lon,scv_data(i).ref.lat);
    scv_data(i).ref.theta    = gsw_CT_from_t(scv_data(i).ref.salt_abs,scv_data(i).ref.temp,scv_data(i).ref.pres);
    scv_data(i).ref.sigma0   = gsw_sigma0(scv_data(i).ref.salt_abs,scv_data(i).ref.theta);
    scv_data(i).ref.spice    = gsw_spiciness0(scv_data(i).ref.salt_abs,scv_data(i).ref.theta);
    [n2_mid,pres_mid]        = gsw_Nsquared(scv_data(i).ref.salt_abs,scv_data(i).ref.theta,scv_data(i).ref.pres,scv_data(i).ref.lat);
    sdat                     = isnan(n2_mid);
    scv_data(i).ref.N2       = interp1(pres_mid(~isnan(sdat)),n2_mid(~isnan(sdat)),scv_data(i).ref.pres);
    clear sdate n2_mid pres_mid 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// PROFILE CHECKS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %// Find densities below mid-pycnocline
    indpyc = find(scv_data(i).sigma0 >= scv_data(i).pyc_dens);

    %// Check that temp profile isn't constantly offset from climatology
    indp = find(scv_data(i).temp_anom(indpyc) > 0); 
    indn = find(scv_data(i).temp_anom(indpyc) < 0); 
    if isempty(indp) == 1 | isempty(indn) == 1 | [length(indp)/length(scv_data(i).temp_anom(indpyc))] > 0.9 | [length(indn)/length(scv_data(i).temp_anom(indpyc))] > 0.9
	bad_temp_prof = [bad_temp_prof i];
	continue
    end

    %// Check that salt profile isn't constantly offset from climatology
    indp = find(scv_data(i).salt_anom(indpyc) > 0); 
    indn = find(scv_data(i).salt_anom(indpyc) < 0); 
    if isempty(indp) == 1 | isempty(indn) == 1 | [length(indp)/length(scv_data(i).salt_anom(indpyc))] > 0.9 | [length(indn)/length(scv_data(i).salt_anom(indpyc))] > 0.9
	bad_salt_prof = [bad_salt_prof i];
	continue
    end

    %// Check that spice profile isn't constantly offset from climatology
    indp = find(scv_data(i).spice_anom(indpyc) > 0); 
    indn = find(scv_data(i).spice_anom(indpyc) < 0); 
    if isempty(indp) == 1 | isempty(indn) == 1 | [length(indp)/length(scv_data(i).spice_anom(indpyc))] > 0.9 | [length(indn)/length(scv_data(i).spice_anom(indpyc))] > 0.9
	bad_spice_prof = [bad_spice_prof i];
	continue
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// ANOMALY CHECKS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %// Check that spice exceeding threshold is localized in depth
    sp_test_high = scv_data(i).spice_anom - [scv_data(i).spice_limits(:,2) + scv_data(i).spice_IQR];
    sp_test_low  = scv_data(i).spice_anom - [scv_data(i).spice_limits(:,1) - scv_data(i).spice_IQR];
    sp_test_high = sp_test_high(~isnan(sp_test_high));
    sp_test_low  = sp_test_low(~isnan(sp_test_low));
    ind = find(sp_test_high > 0 | sp_test_low < 0);
    if ismember(1,ind)==1 %| ismember(length(sp_test_high),ind)==1
	di = diff(ind);
	indd = find(di > 1);
	if isempty(indd)==1
	    bad_spice_anom = [bad_spice_anom i];
	    continue
	end
    end
end

%// Group flags together, remove them
bad_casts = sort(unique([bad_temp_prof bad_salt_prof bad_spice_prof bad_spice_anom bad_reference]));
disp(['Removing ',num2str(length(bad_casts)),' bad SCV detections'])
scv_data(bad_casts) = [];

%// Rebuild minty/spicy SCV datasets with outliers removed
aID        = [scv_data.ID];
mind       = find(ismember(aID,mID)==1);
sind       = find(ismember(aID,sID)==1);
spicy_scvs = scv_data(sind);
minty_scvs = scv_data(mind);

clearvars -except spicy_scvs minty_scvs
save([datadir,'initial_spicy_scv_qc.mat'],'spicy_scvs','-v7.3');
save([datadir,'initial_minty_scv_qc.mat'],'minty_scvs','-v7.3');
