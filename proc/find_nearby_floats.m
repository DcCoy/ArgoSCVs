function [float] = find_nearby_floats(argo_lon,argo_lat,argo_month,argo_day,argo_date,search_circle)
%// Function that takes arrays of Argo float longitude/latitude
%// and returns the index of nearby floats within a search circle radius.
%// For use in argo_get_thresh.m
%// Danny McCoy
%// June 2019

%// Start parallel processing
delete(gcp('nocreate'))
parpool(12)

%// Keep track of progress
progress = [0:100000:length(argo_lon)];
parfor i = 1:length(argo_lon)

    %// Report progress
    if ismember(i,progress)==1
	disp(i)
    end

    %// Find distance
    lon  = argo_lon(i)*ones(size(argo_lon));
    lat  = argo_lat(i)*ones(size(argo_lat));
    dist = sqrt((argo_lon - lon).^2 + (argo_lat - lat).^2);

    %// Re-do if longitude is near zero
    if lon <= 2 | lon >= 358
	lon = lon + 360;
	alon = argo_lon;
	alon(alon<=2) = alon(alon<=2)+360;
	dist = sqrt((alon - lon).^2 + (argo_lat - lat).^2);
    end

    %// Get distance index
    didx = find(dist <= search_circle);

    %// Find floats within +- 1.5 months (year to year)
    month = argo_month(i);
    day   = argo_day(i);
    yr_span =[1995:2020]'; %'
    mn_span = month*ones(size(yr_span));
    dy_span = day*ones(size(yr_span));
    date_span = datenum([yr_span mn_span dy_span]); 
    date_span = [date_span-45 date_span+45];

    %// Of those nearby floats, look for those within the date_span
    tidx = [];
    for ii = 1:length(date_span)
        di   = find(date_span(ii,1) < argo_date & argo_date < date_span(ii,2));
        tidx = [tidx ; di]; 
    end

    %// Find intersection of nearby and neartime profiles
    if isempty(didx) == 0 & isempty(tidx) == 0
	float(i).idx = intersect(didx,tidx);
    else
	float(i).idx = NaN;
    end
end

delete(gcp('nocreate'))
