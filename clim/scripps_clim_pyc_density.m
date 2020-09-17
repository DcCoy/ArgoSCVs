%// Script to find the density of the mid-pycnocline
%// Using peak in N2 as target depth.

close all ; clear all
warning off

%// Addpaths to gsw toolbox
addpath ../proc/scripts/gsw_matlab_03_06_11/
addpath ../proc/scripts/gsw_matlab_03_06_11/library/
addpath ../proc/scripts/gsw_matlab_03_06_11/pdf/
addpath ../proc/scripts/gsw_matlab_03_06_11/thermodynamics_from_t/

%// Load Scripps climatology
disp('Loading SCRIPPS climatology')
load('RG_ArgoClim_TS_Climatologies.mat');
scripps = RG_ArgoClim_TS_Climatologies; clear RG_ArgoClim_TS_Climatologies;

%// Grab lon/lat
slon = [scripps.lon]; slon(slon>360) = slon(slon>360)-360;
slat = [scripps.lat];

%// Fill matrices with NaNs
pyc_dens = nan(length(scripps.lon),length(scripps.lat),12);
pyc_salt = nan(size(pyc_dens));
pyc_temp = nan(size(pyc_dens));
pyc_pres = nan(size(pyc_dens));

%// Scroll through each lat/lon and find density of peak N2
for i = 1:length(scripps.lon)
disp(['---',num2str((i/length(scripps.lon))*100),' percent complete---'])
    for j = 1:length(scripps.lat)
	for t = 1:12
%// Grab data from climatology, calculate other variables
	    lon = slon(i);
	    lat = slat(j);

%// Try adding temperature, if not grab average from nearby profiles
	    temp = squeeze(scripps.temp(i,j,:,t));
	    if isnan(nanmean(temp))==1 & i > 1 & j > 1 & i < length(scripps.lon) & j < length(scripps.lat)
		disp('Averaging nearby climatologies')
		temp1 = squeeze(scripps.temp(i-1,j,:,t));
		temp2 = squeeze(scripps.temp(i+1,j,:,t));
		temp3 = squeeze(scripps.temp(i,j-1,:,t));
		temp4 = squeeze(scripps.temp(i,j+1,:,t));
		temp5 = squeeze(scripps.temp(i+1,j+1,:,t));
		temp6 = squeeze(scripps.temp(i-1,j+1,:,t));
		temp7 = squeeze(scripps.temp(i-1,j-1,:,t));
		temp8 = squeeze(scripps.temp(i+1,j-1,:,t));
		TEMP = [temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8];
		for row = 1:length(TEMP)
		    nan_ind = find(isnan(TEMP(row,:))==0);
		    if length(nan_ind)>0
			ttemp(row) = sum(TEMP(row,nan_ind))/length(nan_ind);
		    else
			ttemp(row) = NaN;
		    end
		end
		temp = ttemp;
		if isnan(nanmean(temp))==1 | length(temp(~isnan(temp))) < 11 %// less than 100m
		    pyc_dens(i,j,t) = NaN; disp('Failed'); continue
		end
	    elseif isnan(nanmean(temp))==1 | length(temp(~isnan(temp))) < 11 %// less than 100m
		pyc_dens(i,j,t) = NaN; disp('Failed'); continue
	    end

%// Try adding salinity, if not grab average from nearby profiles
	    salt = squeeze(scripps.salt(i,j,:,t));
	    if isnan(nanmean(salt))==1 & i > 1 & j > 1 & i < length(scripps.lon) & j < length(scripps.lat)
		salt1 = squeeze(scripps.salt(i-1,j,:,t));
		salt2 = squeeze(scripps.salt(i+1,j,:,t));
		salt3 = squeeze(scripps.salt(i,j-1,:,t));
		salt4 = squeeze(scripps.salt(i,j+1,:,t));
		salt5 = squeeze(scripps.salt(i+1,j+1,:,t));
		salt6 = squeeze(scripps.salt(i-1,j+1,:,t));
		salt7 = squeeze(scripps.salt(i-1,j-1,:,t));
		salt8 = squeeze(scripps.salt(i+1,j-1,:,t));
		SALT = [salt1 salt2 salt3 salt4 salt5 salt6 salt7 salt8];
		for row = 1:length(SALT)
		    nan_ind = find(isnan(SALT(row,:))==0);
		    if length(nan_ind)>0
			tsalt(row) = sum(SALT(row,nan_ind))/length(nan_ind);
		    else
			tsalt(row) = NaN;
		    end
		end
		salt = tsalt;
		if isnan(nanmean(salt))==1 | length(salt(~isnan(salt))) < 11 %// less than 100m
		    pyc_dens(i,j,t) = NaN; disp('Failed'); continue
		end
	    elseif isnan(nanmean(salt))==1 | length(salt(~isnan(salt))) < 11 %// less than 100m
		pyc_dens(i,j,t) = NaN; disp('Failed'); continue
	    end

%// Add remaining properties
	    pres              = scripps.pres;
	    salt_abs          = gsw_SA_from_SP(salt,pres,lon,lat);
	    theta             = gsw_CT_from_t(salt_abs,temp,pres);
	    sigma0            = gsw_sigma0(salt_abs,theta);
	    [N2_mid,pres_mid] = gsw_Nsquared(salt_abs,theta,pres,lat);
	    dat               = isnan(N2_mid);
	    N2                = interp1(pres_mid(dat==0),N2_mid(dat==0),pres);

%// Find max N2, grab density
	    ind = find(N2 == max(N2));
	    pyc_dens(i,j,t) = sigma0(ind+1);
	    pyc_pres(i,j,t) = pres(ind+1);
	    pyc_temp(i,j,t) = temp(ind+1);
	    pyc_salt(i,j,t) = salt(ind+1);
	end
    end
end

save('scripps_pycno_dens','pyc_dens','pyc_pres','slon','slat');

