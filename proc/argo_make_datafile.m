%// argo_make_datafile.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This scripts takes the initial detection IDs and creates a much
%// easier to use datafile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%// Load all data, reduce to only detected SCVs and save data
clearvars -except steps ; close all
load('../datadir.mat');
disp('load final flags, scv indicies, and argo pycnocline density')
load([datadir,'global_initial_scv_detections.mat'])

%// Load raw data
disp('load raw')
load([datadir,'RAW_global_argo.mat'],'argo');
IDr = [argo.ID]; indr = find(ismember(IDr,flag.total)==1); argo(indr) = []; IDr(indr) = [];
argo_raw = argo; clear argo
ind             = find(ismember(IDr,spicy_scv_index)==1);
spicy_scv_raw   = argo_raw(ind);
ind             = find(ismember(IDr,minty_scv_index)==1);
minty_scv_raw   = argo_raw(ind);
clear argo_raw

%// Load proc data
disp('load proc')
load([datadir,'PROC_global_argo.mat'],'argo');
IDd = [argo.ID]; indd = find(ismember(IDd,flag.total)==1); argo(indd) = []; IDd(indd) = [];
argo_proc = argo; clear argo
ind             = find(ismember(IDd,spicy_scv_index)==1);
spicy_scv_proc  = argo_proc(ind);
spicy_scv_pycd  = argo_pyc_dens(ind);
ind             = find(ismember(IDd,minty_scv_index)==1);
minty_scv_proc  = argo_proc(ind);
minty_scv_pycd  = argo_pyc_dens(ind);
clear argo_proc argo_pyc_dens

%// Load anom data
disp('load anom')
load([datadir,'ANOM_global_argo.mat'],'argo_anom');
IDa = [argo_anom.ID]; inda = find(ismember(IDa,flag.total)==1); argo_anom(inda) = []; IDa(inda) = [];
ind              = find(ismember(IDa,spicy_scv_index)==1);
spicy_scv_anom   = argo_anom(ind);
ind              = find(ismember(IDa,minty_scv_index)==1);
minty_scv_anom   = argo_anom(ind);
clear argo_anom

%// Load thresh data
disp('load thresh')
load([datadir,'THRESH_global_argo.mat'],'argo_thresh');
IDt = [argo_thresh.ID]; indt = find(ismember(IDt,flag.total)==1); argo_thresh(indt) = []; IDt(indt) = [];
ind              = find(ismember(IDt,spicy_scv_index)==1);
spicy_scv_thresh = argo_thresh(ind);
ind              = find(ismember(IDt,minty_scv_index)==1);
minty_scv_thresh = argo_thresh(ind);
clear argo_thresh

%// Make new structures with all data
disp('get spicy scv data')
for i = 1:length(spicy_scv_proc)
	spicy_scv(i).float        = spicy_scv_proc(i).float;
	spicy_scv(i).cycle        = spicy_scv_proc(i).cycle;
	spicy_scv(i).ID           = spicy_scv_proc(i).ID;
	spicy_scv(i).lon          = spicy_scv_proc(i).lon;
	spicy_scv(i).lat          = spicy_scv_proc(i).lat;
	spicy_scv(i).time         = spicy_scv_proc(i).time;
	spicy_scv(i).data_mode    = spicy_scv_proc(i).data_mode;
	spicy_scv(i).pyc_dens     = spicy_scv_pycd(i);
	spicy_scv(i).temp         = spicy_scv_proc(i).temp;
	spicy_scv(i).salt         = spicy_scv_proc(i).salt;
	spicy_scv(i).pres         = spicy_scv_proc(i).pres;
	spicy_scv(i).spice        = spicy_scv_proc(i).spice;
	spicy_scv(i).N2           = spicy_scv_proc(i).N2;
	spicy_scv(i).sigma0       = spicy_scv_proc(i).sigma0;
	spicy_scv(i).temp_anom    = spicy_scv_anom(i).temp_anom;
	spicy_scv(i).salt_anom    = spicy_scv_anom(i).salt_anom;
	spicy_scv(i).spice_anom   = spicy_scv_anom(i).spice_anom;
	spicy_scv(i).N2_anom      = spicy_scv_anom(i).N2_anom;
	spicy_scv(i).pres_anom    = spicy_scv_anom(i).pres_anom;
	spicy_scv(i).spice_limits = spicy_scv_thresh(i).spice_limits;
	spicy_scv(i).N2_limits    = spicy_scv_thresh(i).N2_limits;
	spicy_scv(i).spice_IQR    = spicy_scv_thresh(i).spice_IQR;
	spicy_scv(i).N2_IQR       = spicy_scv_thresh(i).N2_IQR;
	spicy_scv(i).RAW_DATA     = spicy_scv_raw(i);          
end

disp('get minty scv data')
for i = 1:length(minty_scv_proc)
	minty_scv(i).float        = minty_scv_proc(i).float;
	minty_scv(i).cycle        = minty_scv_proc(i).cycle;
	minty_scv(i).ID           = minty_scv_proc(i).ID;
	minty_scv(i).lon          = minty_scv_proc(i).lon;
	minty_scv(i).lat          = minty_scv_proc(i).lat;
	minty_scv(i).time         = minty_scv_proc(i).time;
	minty_scv(i).data_mode    = minty_scv_proc(i).data_mode;
	minty_scv(i).pyc_dens     = minty_scv_pycd(i);
	minty_scv(i).temp         = minty_scv_proc(i).temp;
	minty_scv(i).salt         = minty_scv_proc(i).salt;
	minty_scv(i).pres         = minty_scv_proc(i).pres;
	minty_scv(i).spice        = minty_scv_proc(i).spice;
	minty_scv(i).N2           = minty_scv_proc(i).N2;
	minty_scv(i).sigma0       = minty_scv_proc(i).sigma0;
	minty_scv(i).temp_anom    = minty_scv_anom(i).temp_anom;
	minty_scv(i).salt_anom    = minty_scv_anom(i).salt_anom;
	minty_scv(i).spice_anom   = minty_scv_anom(i).spice_anom;
	minty_scv(i).N2_anom      = minty_scv_anom(i).N2_anom;
	minty_scv(i).pres_anom    = minty_scv_anom(i).pres_anom;
	minty_scv(i).spice_limits = minty_scv_thresh(i).spice_limits;
	minty_scv(i).N2_limits    = minty_scv_thresh(i).N2_limits;
	minty_scv(i).spice_IQR    = minty_scv_thresh(i).spice_IQR;
	minty_scv(i).N2_IQR       = minty_scv_thresh(i).N2_IQR;
	minty_scv(i).RAW_DATA     = minty_scv_raw(i);   
end

%// Save results
disp('Save results')
fname = [datadir,'global_initial_scv_data'];
save(fname,'spicy_scv','minty_scv')
