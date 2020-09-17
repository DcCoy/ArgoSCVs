%// argo_get_initial_scvs.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to detect instances of high/low spiciness coupled
%// with low N2 (anticyclonic spicy/minty scvs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%// Clear unwanted variables
clearvars -except N2 argo_anom argo_thresh flag steps
close all
load('../datadir.mat');

%// Default values
ts_IQR      = [1.5]; %// q1 - 1.5*IQR & q3 + 1.5*IQR
n2_IQR      = [1.5]; %// same
N2_range    = [0.2]; %// within +- 0.1 kg/m^3 of high/low spice anom

%// Load premade Scripps mid-pycnocline density matrix, fix longitudes
load('../clim/scripps_pycno_dens.mat');
slon(slon<0)   = slon(slon<0)+360;
slon(slon>360) = slon(slon>360)-360;

%// Grab ID/lon/lat of Argo casts (fix lon)
IDa      = [argo_anom.ID];
alon     = [argo_anom.lon]; alon(alon<0) = alon(alon<0)+360;
alat     = [argo_anom.lat];
sigma0   = [argo_anom.sigma0];

%// Get pycnocline density
disp('Find pycnocline density')
progress = [0:100000:length(argo_anom)]; 
for i = 1:length(argo_anom)

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str(floor((i/length(argo_anom))*100)),' %'])
    end

    %// Grab pycnocline density from scripps pre-made data  
    lon_ind          = min(find(abs(slon-alon(i)) == min(abs(slon-alon(i)))));
    lat_ind          = min(find(abs(slat-alat(i)) == min(abs(slat-alat(i)))));
    argo_pyc_dens(i) = pyc_dens(lon_ind,lat_ind,argo_anom(i).time(2));

    %// Estimate from cast if climatological routine fails
    if isnan(argo_pyc_dens(i))==1 | isempty(argo_pyc_dens(i))==1
	castN2           = N2(:,i);
	castSIG          = sigma0(:,i);
	pyc_ind          = find(castN2 == max(castN2));
	argo_pyc_dens(i) = castSIG(pyc_ind(1));
    end
end
clear slon slat progress alat alon pyc_dens pyc_pres pyc_ind

%// Get matrices of data
disp('Get matrices of data');
spice_anom = [argo_anom.spice_anom];
buoy_anom  = [argo_anom.N2_anom];
clear argo_anom

%// Get matrices of limits
disp('Get anomaly thresholds')
spice_limits = [argo_thresh.spice_limits];
buoy_limits  = [argo_thresh.N2_limits];

%// Get matrices of IQRs
disp('Get IQRs')
spice_IQR = [argo_thresh.spice_IQR];
buoy_IQR  = [argo_thresh.N2_IQR];
clear argo_thresh

%// Get upper/lower limits
disp('Define upper/lower limits')
spice_lowr = spice_limits(:,1:2:end) - ts_IQR*spice_IQR;
spice_high = spice_limits(:,2:2:end) + ts_IQR*spice_IQR;
buoy_lowr  = buoy_limits(:,1:2:end)  - n2_IQR*buoy_IQR;
buoy_high  = buoy_limits(:,2:2:end)  + n2_IQR*buoy_IQR;
clear spice_limits buoy_limits

%// Define tests for spicy
disp('Define spicy/minty tests')
ac_n_test     = [buoy_lowr  - buoy_anom];  % low buoy, anticyclonic
c_n_test      = [buoy_anom  - buoy_high];  % high buoy, cyclonic
spicy_sp_test = [spice_anom - spice_high]; % high spice
minty_sp_test = [spice_lowr - spice_anom]; % low spice
clear buoy_anom buoy_high buoy_lowr

%// Get indices of spicy/minty Anticyclonic SCVs
disp('FINDING ANTICYLCONIC SPICY SCVs')
[spicy_index]   = get_initial_scvs(spicy_sp_test,ac_n_test,argo_pyc_dens,sigma0,N2_range);
spicy_scv_index = IDa(spicy_index);
disp([num2str(length(spicy_index)),' POTENTIAL SPICY SCVs FOUND'])

disp('FINDING ANTICYLONIC MINTY SCVs')
[minty_index]   = get_initial_scvs(minty_sp_test,ac_n_test,argo_pyc_dens,sigma0,N2_range);
minty_scv_index = IDa(minty_index);
disp([num2str(length(minty_index)),' POTENTIAL MINTY SCVs FOUND'])

%// Save data
mname = [datadir,'global_initial_scv_detections.mat'];
save(mname,'spicy_scv_index','minty_scv_index','argo_pyc_dens','flag')
