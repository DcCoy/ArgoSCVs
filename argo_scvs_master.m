%// Main script to load, quality control, and process Argo data, and find initial SCVs
%// Before running this script, you must download all Argo data from a GDAC server
%// See REAMDE_argo_download.txt for instructions
%// You will also need to download the Scripps Argo climatology
%// See clim/REAMDE_scripps.txt for instructions
%// After initial data is downloaded into pacific, atlantic, indian directories, run this code and choose where to begin.
%// Daniel McCoy
%// September 15, 2020

clear all; close all

%// Set data directory here
datadir = ['/data/project1/demccoy/data/argo/update/'];
save('datadir.mat','datadir);

%// Decide where to start the code
steps = zeros(1,12);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('        ARGO SCVs        ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('STEP 1: LOAD ALL DATA')
disp('STEP 2: QUALITY CONTROL')
disp('STEP 3: PROCESSING')
disp('STEP 4: GET ANOMALIES')
disp('Step 5: GET THRESHOLDS')
disp('Step 6: FIND SCV CANDIDATES')
disp('Step 7: MAKE INITIAL FILE')
disp('Step 8: QC INITIAL CANDIDATES')
disp('Step 9: FIND INDIVIDUAL SCVS')
disp('Step10: CLEAN INDIVIDUAL SCV FILE');
disp('Step11: FIND SCV TIME-SERIES');
disp('Step12: CLEAN TIME-SERIES DATA');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
choice = input('Where do we begin? :');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
steps(choice:end) = 1;

% Change to processing directory
cd proc

%// Load data and perform basic QC (obviously bad T,S,P)
disp('STEP 1: LOADING DATA')    
if steps(1) == 1
	argo_load_good_data
end

%// Quality control Argo profiles
disp('STEP 2: QUALITY CONTROL')
if steps(1) == 1 & steps(2) == 1
    	argo_qc_data
elseif steps(1) == 0 & steps(2) == 1
    	load([datadir,'RAW_global_argo.mat']);
    	argo_qc_data
end

%// Process QC'd Argo profiles
disp('STEP 3: PROCESSING')
if steps(2) == 1 & steps(3) == 1
    	argo_proc_data
elseif steps(2) == 0 & steps(3) == 1
    	load([datadir,'QC_global_argo.mat']);
    	argo_proc_data
end

%// Get climatologies for all floats
disp('STEP 4: CALCULATE ANOMALIES')
if steps(3) == 1 & steps(4) == 1
    	argo_get_anom
elseif steps(3) == 0 & steps(4) == 1
    	load([datadir,'PROC_global_argo.mat']);
    	argo_get_anom
end

%// Get IQR thresholds for all floats
disp('STEP 5: CALCULATE THRESHOLDS')
if steps(4) == 1 & steps(5) == 1
    	argo_get_thresh
elseif steps(4) == 0 & steps(5) == 1
    	load([datadir,'ANOM_global_argo.mat']);
    	argo_get_thresh
end

%// Get IQR thresholds for all floats
disp('STEP 6: FIND SCV CANDIDATES')
if steps(5) == 1 & step(6) == 1
    	argo_get_initial_scvs
elseif steps(5) == 0 & steps(6) == 1
    	load([datadir,'THRESH_global_argo.mat'],'flag');
    	load([datadir,'PROC_global_argo.mat']);
    	%// Just grab N2 for routine
    	ID  = [argo.ID];
    	ind = find(ismember(ID,flag.total)==1);
    	argo(ind) = [];
    	N2 = [argo.N2];
    	clear argo ID
    	load([datadir,'ANOM_global_argo.mat']);
    	load([datadir,'THRESH_global_argo.mat']);
    	IDa = [argo_anom.ID];
    	ind = find(ismember(IDa,flag.total)==1); %//apply flags from thresh routine
    	argo_anom(ind) = [];
    	argo_get_initial_scvs
end

%// Make initial datafile
disp('STEP 7: MAKING INITIAL DATA FILE')
if steps(7) == 1
    	argo_make_datafile
end

%// QC initial detections
if steps(8) == 1
	initial_QC
end

%// Find individual SCVs
if steps(9) == 1
	get_spicy_scvs
	get_minty_scvs
end

%// Clean individual SCV data
if steps(10) == 1
	clean_individual_scvs
end

%// Find SCV time-series from the same float
if steps(11) == 1
	spicy_scvs_TS
	minty_scvs_TS
end

%// Clean time-series SCV data
if steps(12) == 1
	clean_timeseries_scvs
end

