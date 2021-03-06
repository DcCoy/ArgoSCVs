%// Script to group SCVs by float ID and save them for later

clear all; close all
load('../datadir.mat');

%// Load good Spicy SCVs, remove duplicates
load([datadir,'good_spicy_scvs.mat'])

%// Set consecutive cast limit
N = 3;

%// Grab float,cast,date,location and ID
for i = 1:length(scv_data)
    floatnum(i) = str2num(scv_data(i).float);
    castnum(i)  = scv_data(i).cycle;
    sdate(i)    = datenum(scv_data(i).time);
end
sID  = [scv_data.ID];
slon = [scv_data.lon];
slat = [scv_data.lat];

%// Sort by floats
[a,b]    = sort(floatnum);
floatnum = floatnum(b);
castnum  = castnum(b);
sID      = sID(b);
slon     = slon(b);
slat     = slat(b);
sdate    = sdate(b);

%// Flip dimensions
floatnum  = floatnum';
castnum   = castnum';
sID       = sID'; 
slon      = slon'; 
slat      = slat'; 
sdate     = sdate';
floatlist = unique(floatnum);

%// Scroll through floatlist and group consecutive casts by ID
scv_list = [];
cnt      = [1];
for i = 1:length(floatlist)
    ind   = find(floatnum == floatlist(i));
    casts = castnum(ind);
    if length(casts) == 1
	scv_list(cnt).IDs = sID(ind);
	cnt = cnt + 1;
	continue
    else
	[a,b] = sort(castnum(ind));
	[c,ia,ic] = unique(a);
	tID = sID(ind(b(ia)));
	casts = castnum(ind(b(ia)));
	a = diff(casts);
	b = find([a;inf]>3);
	c = diff([0;b]);
        d = cumsum(c);
	idx_start = 1;
	for j = 1:length(d)
	    idx = idx_start:1:d(j);
	    scv_list(cnt).IDs = tID(idx);
	    idx_start = idx(end)+1;
	    cnt = cnt + 1;
	end
    end
end

%// Rebuild final SCV structure
scv_count = [1];
for i = 1:length(scv_list)

    %// Grab IDs of consecutive casts from float
    IDs = [scv_list(i).IDs];

    %// Scroll through and rebuild
    for kk = 1:length(IDs)
	ind                                          = find(ismember(sID,IDs(kk))==1);
	spicy_scv(scv_count).float{kk}               = scv_data(ind).float;
	spicy_scv(scv_count).ID(kk)                  = scv_data(ind).ID;
	spicy_scv(scv_count).cycle(kk)               = scv_data(ind).cycle;
	spicy_scv(scv_count).lon(kk)                 = scv_data(ind).lon;
	spicy_scv(scv_count).lat(kk)                 = scv_data(ind).lat;
	spicy_scv(scv_count).time{kk}                = scv_data(ind).time;
	spicy_scv(scv_count).temp{kk}                = scv_data(ind).temp;
	spicy_scv(scv_count).salt{kk}                = scv_data(ind).salt;
	spicy_scv(scv_count).pres{kk}                = scv_data(ind).pres;
	spicy_scv(scv_count).spice{kk}               = scv_data(ind).spice;
	spicy_scv(scv_count).N2{kk}                  = scv_data(ind).N2;
	spicy_scv(scv_count).sigma0{kk}              = scv_data(ind).sigma0;
	spicy_scv(scv_count).temp_anom{kk}           = scv_data(ind).temp_anom;
	spicy_scv(scv_count).salt_anom{kk}           = scv_data(ind).salt_anom;
	spicy_scv(scv_count).pres_anom{kk}           = scv_data(ind).pres_anom;
	spicy_scv(scv_count).N2_anom{kk}             = scv_data(ind).N2_anom;
	spicy_scv(scv_count).spice_anom{kk}          = scv_data(ind).spice_anom;
	spicy_scv(scv_count).N2_limits{kk}           = scv_data(ind).N2_limits;
	spicy_scv(scv_count).spice_limits{kk}        = scv_data(ind).spice_limits;
	spicy_scv(scv_count).spice_IQR{kk}           = scv_data(ind).spice_IQR;
	spicy_scv(scv_count).N2_IQR{kk}              = scv_data(ind).N2_IQR;
	spicy_scv(scv_count).RAW_DATA(kk)            = scv_data(ind).RAW_DATA;
	spicy_scv(scv_count).ref(kk)                 = scv_data(ind).ref;
	spicy_scv(scv_count).stats(kk)               = scv_data(ind).stats;
	spicy_scv(scv_count).limits(kk)              = scv_data(ind).limits;
	spicy_scv(scv_count).gauss(kk)               = scv_data(ind).gauss;
	spicy_scv(scv_count).dyn_height{kk}          = scv_data(ind).dyn_height;
	spicy_scv(scv_count).dyn_height_anom{kk}     = scv_data(ind).dyn_height_anom;
	spicy_scv(scv_count).dyn_height_anom_BC1{kk} = scv_data(ind).dyn_height_anom_BC1;
    end

    %// Increase count
    scv_count = scv_count + 1;
end

%// Save final SCVs and all casts
fname = [datadir,'final_spicy_scvs.mat'];
save(fname,'spicy_scv','-v7.3') 
clearvars -except spicy_scv

	
return
    

    
  







