%// argo_get_thresh.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to add IQR thresholds to temp, salt, spice, 
%// and N2 anomaly profiles using surrounding float profile anomalies.
%// This takes a long time to run!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%// Clear unwanted variables
clearvars -except argo_anom flag steps
close all
load('../datadir.mat')

%// Default Parameters
search_circle = 2;   %// 2 degrees

%// Make lon/lat matrix
alon = [argo_anom.lon]'; %'
alat = [argo_anom.lat]'; %'
for i = 1:length(argo_anom)
    amnth(i) = argo_anom(i).time(2);
    aday(i)  = argo_anom(i).time(3);
    adate(i) = datenum(argo_anom(i).time);
end
amnth = amnth'; aday = aday'; adate = adate'; %'

%// Find nearby floats using parfoor routine
disp('Finding nearby/neartime floats')
[float] = find_nearby_floats(alon,alat,amnth,aday,adate,search_circle)
save('float_index.mat','float','-v7.3')
disp('Done')
clear alon alat search_circle amnth aday adate

%// Make spice/N2 anom matrices, also grab density levels
disp('Making preIQR_data.mat')
spice_anom = [argo_anom.spice_anom];
N2_anom    = [argo_anom.N2_anom];
sigma0     = [argo_anom.sigma0];
save('preIQR_data.mat','spice_anom','N2_anom','sigma0');
clearvars -except spice_anom N2_anom sigma0 float
disp('Done')

%// Load premade float index
disp('Load nearby float index')
load float_index

%// To limit memory usage, only load one variable at a time
%// i.e. load spice anomalies, calculate IQR, then clear all
%// then load N2 anomalies etc.
load('preIQR_data.mat','spice_anom','sigma0')
disp('...spiciness')
[spice_IQR spice_p25 spice_p75 spflags] = get_IQR_data(spice_anom,sigma0,float);
save('spice_iqr_stats.mat','spice_IQR','spice_p25','spice_p75','spflags','-v7.3');
clear spice_IQR spice_p25 spice_p75 spflags spice_anom
load('preIQR_data.mat','N2_anom')
disp('...buoyancy frequency')
[N2_IQR    N2_p25    N2_p75    n2flags] = get_IQR_data(N2_anom,sigma0,float);
save('N2_iqr_stats.mat','N2_IQR','N2_p25','N2_p75','n2flags','-v7.3');
clear N2_IQR N2_p25 N2_p75 n2flags N2_anom float
disp('DONE')

%// Load cleared data
load('spice_iqr_stats');
load('N2_iqr_stats');
load([datadir,'ANOM_global_argo.mat']);

%// Get all flags
flag_idx = spflags+n2flags;

%// Build argo_thresh matrix
%// Keep track of progress
progress = [0:100000:length(argo_anom)];
for i = 1:length(argo_anom)

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str((i/length(argo_anom))*100),' %'])
    end

    if flag_idx(i) == 0
	%// Save results in separate struct (if parallel processing)
	argo_thresh(i).float        = argo_anom(i).float;
	argo_thresh(i).cycle        = argo_anom(i).cycle;
	argo_thresh(i).lon          = argo_anom(i).lon;
	argo_thresh(i).lat          = argo_anom(i).lat;
	argo_thresh(i).time         = argo_anom(i).time;
	argo_thresh(i).ID           = argo_anom(i).ID;
	argo_thresh(i).spice_IQR    = spice_IQR(:,i);
	argo_thresh(i).spice_limits = single([spice_p25(:,i) spice_p75(:,i)]);
	argo_thresh(i).N2_IQR       = N2_IQR(:,i);
	argo_thresh(i).N2_limits    = single([N2_p25(:,i) N2_p75(:,i)]);
	argo_thresh(i).sigma0       = argo_anom(i).sigma0;
    end
end

%// Define IQR flag and remove those profiles
ID                      = [argo_anom.ID];
flag.IQR                = ID(flag_idx>0);
flag.total              = [flag.total flag.IQR];
argo_thresh(flag_idx>0) = [];
argo_anom(flag_idx>0)   = [];

%// Save data
fname = [datadir,'THRESH_global_argo.mat'];
save(fname,'argo_thresh','flag','-v7.3')

clearvars -except argo argo_anom argo_thresh flag steps

