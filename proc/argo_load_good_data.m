%// argo_load_good_data.m
%// Daniel McCoy - May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// This script is used to load Argo data for SCV analysis. It goes through each
%// basin (Pacific, Atlantic, Indian) for each day and loads all the available 
%// Argo float profiles. The script only allows data that has been flagged either
%// 1 (good) or 2 (probably good) at each pressure level. Delayed mode or adjusted
%// real-time profiles also load the adjusted, QC'd data (same flags are applied). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Other Data Flags:
%// Obviously bad pressure    ( 0 to 10000 dbar)
%// Obviously bad temperature (-4 to 40 degC)
%// Obviously bad salinity    (10 to 45 PSU)

clearvars -except steps
close all
load('../datadir.mat')

%// Establish bad min/max oceanic values for p,t,s
%//p_lim = [0 10000];
%//t_lim = [-4 40];
%//s_lim = [10 45];

%// Keep track of total floats, for structures
count = 1;
for y = 1995:2020 %// all years
    disp(['YEAR ',num2str(y)])
    disp(['TOTAL ARGO PROFILES SO FAR: ',num2str(count-1)])
    for m = 1:12 %// all months
	for d = 1:31 %// all days
	    for b = 1:3 %// all basins

		%// Select directory for data (Pac --> Atl --> Ind)
		if b == 1
		    diri = ['/data/project1/data/argo/pub/outgoing/argo/geo/pacific_ocean/'];
		elseif b == 2
		    diri = ['/data/project1/data/argo/pub/outgoing/argo/geo/atlantic_ocean/'];
		elseif b == 3
		    diri = ['/data/project1/data/argo/pub/outgoing/argo/geo/indian_ocean/'];
		end
		
		%// Grab filename string
		if d < 10 & m < 10
		    fname = [diri,num2str(y),'/','0',num2str(m),'/',num2str(y),'0',num2str(m),'0',num2str(d),'_prof.nc'];
		elseif d < 10 & m >= 10
		    fname = [diri,num2str(y),'/',num2str(m),'/',num2str(y),num2str(m),'0',num2str(d),'_prof.nc'];
		elseif d >= 10 & m < 10
		    fname = [diri,num2str(y),'/','0',num2str(m),'/',num2str(y),'0',num2str(m),num2str(d),'_prof.nc'];
		elseif d >= 10 & m >= 10
		    fname = [diri,num2str(y),'/',num2str(m),'/',num2str(y),num2str(m),num2str(d),'_prof.nc'];
		end

		%// Load float meta, profile, and QC data for that day
		try
		    float = ncread(fname,'PLATFORM_NUMBER');
		    cast = ncread(fname,'CYCLE_NUMBER');
		    lat = ncread(fname,'LATITUDE');
		    lon = ncread(fname,'LONGITUDE');
		    time = ncread(fname,'JULD');
		    mode = ncread(fname,'DATA_MODE');
		    profile_mode = ncread(fname,'DIRECTION');
		    temp = ncread(fname,'TEMP');
		    temp_qc = ncread(fname,'TEMP_QC');
		    temp_adj = ncread(fname,'TEMP_ADJUSTED');
		    temp_adj_qc = ncread(fname,'TEMP_ADJUSTED_QC');
		    salt = ncread(fname,'PSAL');
		    salt_qc = ncread(fname,'PSAL_QC');
		    salt_adj = ncread(fname,'PSAL_ADJUSTED');
		    salt_adj_qc = ncread(fname,'PSAL_ADJUSTED_QC');
		    pres = ncread(fname,'PRES');
		    pres_qc = ncread(fname,'PRES_QC');
		    pres_adj = ncread(fname,'PRES_ADJUSTED');
		    pres_adj_qc = ncread(fname,'PRES_ADJUSTED_QC');

		    %// For all the profiles from that day, check for bad points and remove them
		    for f = 1:length(mode)

			%// Ignore descending profiles (only want ascending)
			if strmatch('D',profile_mode(f))==1
			    continue
			end

			%// Save meta data
			argo(count).float = [deblank(char(float(:,f))')]; %'
			argo(count).cycle = [cast(f)]; 
			argo(count).lon = lon(f);
			argo(count).lat = lat(f);
			argo(count).time = datevec(time(f)+datenum([1950 01 01 0 0 0]));
			argo(count).data_mode = mode(f);
						
			%// Check for adjusted data or delayed mode, only keep flags 1 or 2 (good or probably good)
			if strmatch(mode(f),['A';'D'])>0
			    indt = find(str2num(temp_adj_qc(:,f))<3);
			    inds = find(str2num(salt_adj_qc(:,f))<3);
			    indp = find(str2num(pres_adj_qc(:,f))<3);
			    indts = intersect(indt,inds);
			    indtsp = intersect(indts,indp);
			    if isempty(indtsp) == 1
				continue
			    elseif length(indtsp) < 25
				continue
			    end
			    argo(count).temp = temp_adj(indtsp,f);
			    argo(count).salt = salt_adj(indtsp,f);
			    argo(count).pres = pres_adj(indtsp,f);

			    %// Remove obviously bad points
			    argo(count).temp((argo(count).temp > 40)) = NaN;
			    argo(count).temp((argo(count).temp < -4)) = NaN;
			    argo(count).salt((argo(count).salt > 40)) = NaN;
			    argo(count).salt((argo(count).salt < 10)) = NaN;
			    argo(count).pres((argo(count).pres < 0))  = NaN;
			    argo(count).pres((argo(count).pres > 1e4))= NaN;

			    %// Remove any other levels where data is NaN
			    badline = [argo(count).pres + argo(count).temp + argo(count).salt];
			    ind = find(isnan(badline)==1);
			    argo(count).pres(ind) = [];
			    argo(count).salt(ind) = [];
			    argo(count).temp(ind) = [];

			%// If no delayed or adjusted, check for good real-time data flags (will remove grey-list floats later)
			elseif strmatch('R',mode(f))==1
			    indt = find(str2num(temp_qc(:,f))<3);
			    inds = find(str2num(salt_qc(:,f))<3);
			    indp = find(str2num(pres_qc(:,f))<3);
			    indts = intersect(indt,inds);
			    indtsp = intersect(indts,indp);
			    if isempty(indtsp) ==1
				continue
			    elseif length(indtsp) < 25
				continue
			    end
			    argo(count).temp = temp(indtsp,f);
			    argo(count).salt = salt(indtsp,f);
			    argo(count).pres = pres(indtsp,f);

			    %// Remove obviously bad points
			    argo(count).temp((argo(count).temp > 40)) = NaN;
			    argo(count).temp((argo(count).temp < -4)) = NaN;
			    argo(count).salt((argo(count).salt > 40)) = NaN;
			    argo(count).salt((argo(count).salt < 10)) = NaN;
			    argo(count).pres((argo(count).pres < 0))  = NaN;
			    argo(count).pres((argo(count).pres > 1e4))= NaN;

			    %// Remove any other levels where data is NaN
			    badline = [argo(count).pres + argo(count).temp + argo(count).salt];
			    ind = find(isnan(badline)==1);
			    argo(count).pres(ind) = [];
			    argo(count).salt(ind) = [];
			    argo(count).temp(ind) = [];
			end

			%// Assuming all is good, increase count and move onto next float
			count = count + 1;
		    end
		catch
		    %// Move on, no profiles found
		end
	    end    
	end
    end	
end

%// Add my own ID number and make data 'single' to save memory
for i = 1:length(argo)
    argo(i).ID    = i;
    argo(i).cycle = single(argo(i).cycle);
    argo(i).lon   = single(argo(i).lon);
    argo(i).lat   = single(argo(i).lat);
    argo(i).temp  = single(argo(i).temp);
    argo(i).salt  = single(argo(i).salt);
    argo(i).pres  = single(argo(i).pres);
end
		
%// Save data
fname = [datadir,'RAW_global_argo.mat'];
save(fname,'argo','-v7.3')

clearvars -except argo steps
