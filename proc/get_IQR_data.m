function [data_IQR data_p25 data_p75 flag_idx] = get_IQRs(anom,sigma0,float);
%// Function that takes matrix of anomalies (chosen by input)
%// and returns the IQR threshold (typically 1st/3rd quartile -/+ 1.5*IQR, respectively)
%// along isopycnal surfaces. Needs 'float' structure from 'find_nearby_floats.m')
%// that specifies, for each float, the indicies of nearby floats
%// Used in argo_get_thresh.m
%// Danny McCoy
%// JUNE 21 2019

%// Start parallel workers
delete(gcp('nocreate'));
parpool(12);

%// Specify minimum number of samples for IQR
N = 60;

%// Initiate matrices
disp('Initiate matrices')
flag_idx = zeros(1,length(anom));
data_IQR = nan(size(anom));
data_p25 = nan(size(anom));
data_p75 = nan(size(anom));

%// Convert float to cell array
float = struct2cell(float);
float = squeeze(float);

%// Check progress
progress = [0:50000:length(anom)];
%// Start loop
disp('Starting for-loop')
parfor i = 1:length(anom)

    %// Display progress
    if ismember(i,progress)==1
	  i
    end

    %// Check for nearby floats
    if length(float{i}) > N

	%// Grab spice anomaly, N2 anomaly, and their densities from nearby casts
	DATA      = [];
	DATA      = anom(:,float{i});
	DENS      = [];
	DENS      = sigma0(:,float{i});

	%// Find where sample cast has data
	ind       = [];
	ind       = find(isnan(sigma0(:,i))==0);

	%// Initialize data_grid matrix
	[a,~]     = size(anom);
	[b]       = length(float{i});
	data_grid = nan(a,b);

	%// Interpolate each cast to the sample cast densities
	for ii = 1:length(float{i})
	    
	    %// Grab individual cast data
	    data  = [];
	    data  = DATA(:,ii); 
	    dens = [];
	    dens  = DENS(:,ii);

	    %// Only allow unique densities
	    c = []; ia = []; ic = [];
	    [c,ia,ic] = unique(dens);
	    data      = data(ia);
	    dens      = dens(ia);

	    %// Attempt to interpolate to sample cast density
	    %// Keep original NaN levels, fill in matrix
	    data_int  = [];
	    try
		data_int         = interp1(dens(~isnan(data)),data(~isnan(data)),sigma0(ind,i));
		filler           = nan(length(sigma0(:,i)),1);
		filler(ind)      = data_int;
		data_grid(:,ii)  = filler;
	    catch
	    end
	end

	%// Grab interquartile range of values along sample cast's isopycnals
	%// Only allow levels that have at least N points
	[a,b] = size(data_grid);
	for ii = 1:b
	    ind = find(isnan(data_grid(:,ii))==1);
	    if length(ind) < N
		data_grid(:,ii) == NaN;
	    end
	end

	data_IQR(:,i) = iqr(data_grid'); %'
	data_p25(:,i) = prctile(data_grid',25); %'
	data_p75(:,i) = prctile(data_grid',75); %'
    else
	flag_idx(i) = 1; continue;
    end
end

delete(gcp('nocreate'))
    



  
