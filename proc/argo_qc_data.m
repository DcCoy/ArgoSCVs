%// argo_qc_data.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to flag profiles that are unfit for the SCV analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Flags applied:
%// 1. Cycle-0: These are the first 'test' cast for each profile and are rejected
%// 2. Empty: Remove any empty casts (obviously dont want these)
%// 3. Grey-list: Known 'bad' floats, obtained from GDAC server
%// 4. Depths: Need at least 40 data points
%// 5. Resolution: Ignore any data that is >65 dbar from previous pressure (1000 dbar or less)
%// 6. Resolution: Ignore any data that is >101 dbar from previous pressure (1000 dbar or greater)
%// 7. Max Depth: Profile needs to have data that passed the above go to at least 700 dbar

clearvars -except argo steps
load('../datadir.mat');

%// Load greylist (5/29/2019) from GDAC, grab unique float entries
grey                              = xlsread('argo_greylist.xlsx');
[c,ia,ic]                         = unique(grey(:,1:2),'rows');
grey                              = grey(ia,:);
grey(find(isnan(grey(:,1))==1),:) = []; %// Remove NaN

%// Grab float number, start/end date
float_num  = grey(:,1);
start_date = grey(:,2); 
end_date   = grey(:,3);

%// Fix start/end date (matlab time)
yr 			    = num2str(start_date); yr = yr(:,1:4);
mn 			    = num2str(start_date); mn = mn(:,5:6);
dy                          = num2str(start_date); dy = dy(:,7:8);
start_date 		    = datenum(str2num(yr),str2num(mn),str2num(dy));
ind 			    = find(isnan(end_date)==1);
end_date(ind)               = 99999999; %// Filler for now
yr 			    = num2str(end_date); yr = yr(:,1:4);
mn 		            = num2str(end_date); mn = mn(:,5:6);
dy 			    = num2str(end_date); dy = dy(:,7:8);
end_date 		    = datenum(str2num(yr),str2num(mn),str2num(dy));
end_date(end_date > 800000) = NaN;
clear yr mn dy c ia ic ind grey

%// Quality control each Argo profile
flag.cycle      = [];
flag.empty      = [];
flag.grey       = [];
flag.deep       = [];
flag.shallow    = [];
flag.depths     = [];
flag.resolution = [];
flag_idx        = [];

% Should be able to run in a par-for loop
progress = [0:100000:length(argo)];
for i = 1:length(argo) 
    ind = [];

    %// Display progress of for-loop
    if ismember(i,progress)==1
	disp([num2str((i/length(argo))*100),' %'])
    end

    %// Check for 'cycle - 0' profiles (not usable, first 'test' cast)
    if argo(i).cycle == 0
	flag.cycle = [flag.cycle argo(i).ID]; flag_idx = [flag_idx i]; continue
    end

    %// Check for empty profiles
    if isempty(argo(i).pres) == 1 | isempty(argo(i).temp)==1 | isempty(argo(i).salt)==1
	flag.empty = [flag.empty argo(i).ID]; flag_idx = [flag_idx i]; continue
    end

    %// Check for grey-list
    if ismember(str2num(argo(i).float),float_num)==1
	ind        = find(float_num == str2num(argo(i).float));
	grey_start = start_date(ind);
	grey_end   = end_date(ind);
	argo_date  = datenum(argo(i).time);
	if isnan(grey_end)==1
	    if argo_date > grey_start
		flag.grey = [flag.grey argo(i).ID]; flag_idx = [flag_idx i]; continue
	    end
	elseif isnan(grey_end)==0
	    if grey_start <= argo_date & argo_date <= grey_end
		flag.grey = [flag.grey argo(i).ID]; flag_idx = [flag_idx i]; continue
	    end
	end	    
    end

    %// Check minimum # of data points (40)
    if length(unique(argo(i).pres)) < 40
	flag.depths = [flag.depths argo(i).ID]; flag_idx = [flag_idx i]; continue
    end

    %// Check resolution
    %// Check shallow data resolution (<65 dbar)
    ind = find(100 < argo(i).pres & argo(i).pres < 700);
    if max(diff(argo(i).pres(ind)) > 65);
	flag.resolution = [flag.resolution argo(i).ID]; flag_idx = [flag_idx i]; continue
    end

    %// Ignore data below the point of resolution > 65dbar
    %// If we remove too much, flag.shallow will catch it
    ind1 = find(100 < argo(i).pres & argo(i).pres < 1000);
    ind2 = find(diff(argo(i).pres(ind1))>65);
    if isempty(ind2)==0
	argo(i).pres([ind1(ind2(1))+1]:end) = [];
	argo(i).temp([ind1(ind2(1))+1]:end) = [];
	argo(i).salt([ind1(ind2(1))+1]:end) = [];
    end

    %// Ignore data below 1000 dbar with resolution > 105 dbar
    ind1 = find(argo(i).pres > 1000);
    if isempty(ind1)==0
	ind2 = find(diff(argo(i).pres(ind1))>105);
	if isempty(ind2)==0
	    argo(i).pres([ind1(ind2(1))+1]:end) = [];
	    argo(i).temp([ind1(ind2(1))+1]:end) = [];
	    argo(i).salt([ind1(ind2(1))+1]:end) = [];
	end
    end
	 
    %// Check maximum pressure (>700)
    if max(argo(i).pres) < 700
	flag.shallow = [flag.shallow argo(i).ID]; flag_idx = [flag_idx i]; continue
    end 

    %// Check minimum pressure (<100)
    if min(argo(i).pres) > 100
	flag.deep = [flag.deep argo(i).ID]; flag_idx = [flag_idx i]; continue
    end

    %// Check minimum # of data points (40) again (in case we removed a lot with the above)
    if length(unique(argo(i).pres)) < 40
	flag.depths = [flag.depths argo(i).ID]; flag_idx = [flag_idx i]; continue
    end
end	

%// Get total flagged profiles
flag.total = [flag.cycle flag.empty flag.grey flag.shallow flag.deep flag.depths flag.resolution];
disp(['Removing ',num2str(length(flag.cycle)),' cycle-0 profiles'])
disp(['Removing ',num2str(length(flag.empty)),' profiles due to empty data'])
disp(['Removing ',num2str(length(flag.grey)),' profiles due to grey-listing'])
disp(['Removing ',num2str(length(flag.shallow)),' profiles due to max depth'])
disp(['Removing ',num2str(length(flag.deep)),' profiles due to min depth'])
disp(['Removing ',num2str(length(flag.depths)),' profiles due small # of samples'])
disp(['Removing ',num2str(length(flag.resolution)),' profiles due to bad vertical resolution'])
disp(['Removing ',num2str(length(flag.total)),' total profiles due to the above QC'])

%// Remove flagged profiles
argo(flag_idx) = [];

%// Save data
fname = [datadir,'QC_global_argo.mat'];
save(fname,'argo','flag','-v7.3')

clearvars -except argo flag steps
	
