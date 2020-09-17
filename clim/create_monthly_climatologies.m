% Run this script to update the monthly climatology
% of RG Scripps data using the latest file updates
% See README_scripps.txt

file_description = 'SCRIPPS ARGO climatology using data up to June 2017';
add_recent       = true;	% To add rencent monthly climatologies

indir = 'monthly_add';
 
% -----------------------------------------------
% First, loads most recent compete climatology from SCRIPPS
infile1 = 'RG_ArgoClim_Temperature_2017.nc';
infile2 = 'RG_ArgoClim_Salinity_2017.nc';
tmp1 = netcdf_load('dir',indir,'file',infile1);
tmp2 = netcdf_load('dir',indir,'file',infile2);

% Adds anomalies to means to recover monthly climatologies
temp           = bsxfun(@plus,tmp1.argo_temperature_anomaly,tmp1.argo_temperature_mean);
temp(temp<-10) = nan;
salt           = bsxfun(@plus,tmp2.argo_salinity_anomaly,tmp2.argo_salinity_mean);
salt(salt<0)   = nan;
time           = tmp1.time;
tmp            = size(temp);
size3d         = tmp(1:3);

% -----------------------------------------------
% If required adds on the latest monthly climatologies, one by one
% Note according to SCRIPPS website the anomalies are caluclated still
% from the previous climaotlogy years - no update of new months fields?
if (add_recent)
    disp(['Adding latest monthly climatologies']);
    files_to_add = {'201701','201702','201703','201704','201705','201706', ... 
                    '201707','201708','201709','201710','201711','201712', ...
		    '201801','201802','201803','201804','201805','201806', ... 
                    '201807','201808','201809','201810','201811','201812', ...
		    '201901','201902','201903','201904','201905','201906', ...
		    '201907','201908','201909','201910','201911','201912'};
    temp_add = nan([size3d length(files_to_add)]);
    salt_add = nan([size3d length(files_to_add)]);
    time_add = nan(length(files_to_add),1);
    
    for indf=1:length(files_to_add)
       clear tmp3
       infile               = ['RG_ArgoClim_' files_to_add{indf} '.nc'];
       tmp3                 = netcdf_load('dir',indir,'file',infile);
       temp_add(:,:,:,indf) = tmp1.argo_temperature_mean + tmp3.argo_temperature_anomaly;
       salt_add(:,:,:,indf) = tmp2.argo_salinity_mean + tmp3.argo_salinity_anomaly;
       time_add(indf,1)     = tmp3.time;
    end

    % Updates temp and salt adding the latest months
    % WARNING: the month order need to be preserved, to calculate monthly climatologies
    temp_new                           = nan([size3d size(temp,4)+size(temp_add,4)]); 
    salt_new                           = nan([size3d size(salt,4)+size(salt_add,4)]); 
    temp_new(:,:,:,1:size(temp,4))     = temp;
    temp_new(:,:,:,size(temp,4)+1:end) = temp_add;
    salt_new(:,:,:,1:size(salt,4))     = salt;
    salt_new(:,:,:,size(salt,4)+1:end) = salt_add;
    % substitutes and clean
    temp           = temp_new;
    salt           = salt_new;
    temp(temp<-10) = nan;
    salt(salt<0)   = nan;
    clear temp_new salt_new temp_add salt_add
    time           = cat(1,time,time_add);
    % A check
    if any(diff(time)==0);error(['Problem with time vector']);end
end

% -----------------------------------------------
tmonth = mod(time-0.5,12)+1;
month  = [1:12];

temp_clim = nan([size3d 12]);
salt_clim = nan([size3d 12]);

for indm=1:12
    imonth 	 	  = find(tmonth==indm);
    temp_sampled 	  = temp(:,:,:,imonth);
    temp_clim(:,:,:,indm) = squeeze(nanmean(temp_sampled,4));
    salt_sampled 	  = salt(:,:,:,imonth);
    salt_clim(:,:,:,indm) = squeeze(nanmean(salt_sampled,4));
end

% Ads back the NaN mask
nanmask             = tmp1.bathymetry_mask;
nanmask(nanmask~=1) = nan;
 
temp_clim = bsxfun(@times,temp_clim,nanmask);
salt_clim = bsxfun(@times,salt_clim,nanmask);

% Fills in Matlab strucutre
RG_ArgoClim_TS_Climatologies.description = file_description;
RG_ArgoClim_TS_Climatologies.lon         = tmp1.longitude;
RG_ArgoClim_TS_Climatologies.lat         = tmp1.latitude;
RG_ArgoClim_TS_Climatologies.pres        = tmp1.pressure;
RG_ArgoClim_TS_Climatologies.month       = month;
RG_ArgoClim_TS_Climatologies.temp        = temp_clim;
RG_ArgoClim_TS_Climatologies.salt        = salt_clim;

save RG_ArgoClim_TS_Climatologies RG_ArgoClim_TS_Climatologies;

