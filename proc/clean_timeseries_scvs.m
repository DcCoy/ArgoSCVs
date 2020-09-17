% Script to make final time-series SCVs file for sharing

close all ; clear all
load('../datadir.mat');

% SPICY SCVs
load([datadir,'final_spicy_scvs.mat']);
scv_data = [];

% MINTY SCVs
load([datadir,'final_minty_scvs.mat']);
scv_data = [];

% Make pressure array to interpolate climatological dynamic height to regular 0:10:2000 grid
pres_array = [0:10:2000];

% Fill in dynamic height data with NaNs where no data exists (0:10:2000 grid);
% Also get modes on same grid
for i = 1:length(spicy_scv)
	for j = 1:length(spicy_scv(i).lon)
		tmpdyn                         = nan(length(pres_array),1);
		tmppmodes		       = nan(length(pres_array),size(spicy_scv(i).ref(j).pmodes,2));
		tmpwmodes		       = nan(length(pres_array),size(spicy_scv(i).ref(j).pmodes,2));
		idx                            = ismember(pres_array,spicy_scv(i).ref(j).dyn_pres);
		tmpdyn(idx==1)                 = spicy_scv(i).ref(j).dyn_height;
		spicy_scv(i).ref(j).dyn_height = tmpdyn;
		idx                            = ismember(pres_array,spicy_scv(i).ref(j).mode_pres);
		tmppmodes(idx==1,:)            = spicy_scv(i).ref(j).pmodes;
		tmpwmodes(idx==1,:)            = spicy_scv(i).ref(j).wmodes;
		spicy_scv(i).ref(j).pmodes     = tmppmodes;
		spicy_scv(i).ref(j).wmodes     = tmpwmodes;
	end
end

% Make pressure array to interpolate climatological dynamic height to regular 0:10:2000 grid
pres_array = [0:10:2000];

% Fill in dynamic height data with NaNs where no data exists (0:10:2000 grid);
% Also get modes on same grid
for i = 1:length(minty_scv)
	for j = 1:length(minty_scv(i).lon)
		tmpdyn                         = nan(length(pres_array),1);
		tmppmodes		       = nan(length(pres_array),size(minty_scv(i).ref(j).pmodes,2));
		tmpwmodes		       = nan(length(pres_array),size(minty_scv(i).ref(j).pmodes,2));
		idx                            = ismember(pres_array,minty_scv(i).ref(j).dyn_pres);
		tmpdyn(idx==1)                 = minty_scv(i).ref(j).dyn_height;
		minty_scv(i).ref(j).dyn_height = tmpdyn;
		idx                            = ismember(pres_array,minty_scv(i).ref(j).mode_pres);
		tmppmodes(idx==1,:)            = minty_scv(i).ref(j).pmodes;
		tmpwmodes(idx==1,:)            = minty_scv(i).ref(j).wmodes;
		minty_scv(i).ref(j).pmodes     = tmppmodes;
		minty_scv(i).ref(j).wmodes     = tmpwmodes;
	end
end


% Make separate structures for different data types

% meta        = Float ID, location, data mode, etc
% raw         = Raw T/S/P data
% profile     = Interpolated (0:10:2000dbar) profile data
% anomalies   = Anomalies along isopycnals (compared with Scripps climatology)
% thresholds  = Q1(25%), Q3(75%), and interquartile range (Q3 - Q1) for spice and N^2
% climatology = Scripps climatology, interpolated to 0:10:2000 grid
% mode_decomp = Baroclinic horizontal (pmodes) and vertical (wmodes) velocity modes
%               with first horizontal baroclinic mode (BC1) and modal amplitude (alpha) for easy loading
% gaussian    = Gaussian model results (spice_anomaly vs pressure) and R^2 fit results
%               along with the limits of the SCV core and entire SCV in both pressure and density space
% scv_stats   = Various SCV stats, such as hydrographic properties at the SCV core, scale height
%               vertical extent, deformation radius (scale length), coriollis frequency etc. 

for i = 1:length(spicy_scv)
	for j = 1:length(spicy_scv(i).lon);

		% Metadata
		scv_data.spicy_scvs(i).meta(j).ID    = spicy_scv(i).ID(j);
		scv_data.spicy_scvs(i).meta(j).float = spicy_scv(i).float{j};
		scv_data.spicy_scvs(i).meta(j).cycle = spicy_scv(i).cycle(j);
		scv_data.spicy_scvs(i).meta(j).lon   = spicy_scv(i).lon(j);
		scv_data.spicy_scvs(i).meta(j).lat   = spicy_scv(i).lat(j);
		scv_data.spicy_scvs(i).meta(j).time  = spicy_scv(i).time{j};

		% Profile data
		scv_data.spicy_scvs(i).profile(j).temp   = single(spicy_scv(i).temp{j});
		scv_data.spicy_scvs(i).profile(j).salt   = single(spicy_scv(i).salt{j});
		scv_data.spicy_scvs(i).profile(j).spice  = single(spicy_scv(i).spice{j});
		scv_data.spicy_scvs(i).profile(j).N2     = single(spicy_scv(i).N2{j});
		scv_data.spicy_scvs(i).profile(j).pres   = single(spicy_scv(i).pres{j});
		scv_data.spicy_scvs(i).profile(j).sigma0 = single(spicy_scv(i).sigma0{j});

		% Anomaly data
		scv_data.spicy_scvs(i).anomalies(j).temp  = single(spicy_scv(i).temp_anom{j}); 
		scv_data.spicy_scvs(i).anomalies(j).salt  = single(spicy_scv(i).salt_anom{j});
		scv_data.spicy_scvs(i).anomalies(j).spice = single(spicy_scv(i).spice_anom{j});
		scv_data.spicy_scvs(i).anomalies(j).N2    = single(spicy_scv(i).N2_anom{j});
		scv_data.spicy_scvs(i).anomalies(j).pres  = single(spicy_scv(i).pres_anom{j});

		% IQR Thresholds
		scv_data.spicy_scvs(i).thresholds(j).spice.Q1  = single(spicy_scv(i).spice_limits{j}(:,1));
		scv_data.spicy_scvs(i).thresholds(j).spice.Q3  = single(spicy_scv(i).spice_limits{j}(:,2));
		scv_data.spicy_scvs(i).thresholds(j).spice.IQR = single(spicy_scv(i).spice_IQR{j});
		scv_data.spicy_scvs(i).thresholds(j).N2.Q1     = single(spicy_scv(i).N2_limits{j}(:,1));
		scv_data.spicy_scvs(i).thresholds(j).N2.Q3     = single(spicy_scv(i).N2_limits{j}(:,2));
		scv_data.spicy_scvs(i).thresholds(j).N2.IQR    = single(spicy_scv(i).N2_IQR{j});

		% Climatology from Scripps
		scv_data.spicy_scvs(i).climatology(j).lon        = single(spicy_scv(i).ref(j).lon);
		scv_data.spicy_scvs(i).climatology(j).lat        = single(spicy_scv(i).ref(j).lat);
		scv_data.spicy_scvs(i).climatology(j).temp       = single(spicy_scv(i).ref(j).temp);
		scv_data.spicy_scvs(i).climatology(j).salt       = single(spicy_scv(i).ref(j).salt);
		scv_data.spicy_scvs(i).climatology(j).spice      = single(spicy_scv(i).ref(j).spice);
		scv_data.spicy_scvs(i).climatology(j).N2         = single(spicy_scv(i).ref(j).N2);
		scv_data.spicy_scvs(i).climatology(j).sigma0     = single(spicy_scv(i).ref(j).sigma0);
		scv_data.spicy_scvs(i).climatology(j).dyn_height = single(spicy_scv(i).ref(j).dyn_height);
		scv_data.spicy_scvs(i).climatology(j).pres       = single(spicy_scv(i).ref(j).pres);

		% Vertical Mode decomposition of climatological dynamic height anomaly profile
		scv_data.spicy_scvs(i).mode_decomp(j).pmodes   = single(spicy_scv(i).ref(j).pmodes);
		scv_data.spicy_scvs(i).mode_decomp(j).wmodes   = single(spicy_scv(i).ref(j).wmodes);
		scv_data.spicy_scvs(i).mode_decomp(j).alphaBC1 = single(spicy_scv(i).ref(j).VMD.alpha);
		scv_data.spicy_scvs(i).mode_decomp(j).pmodeBC1 = single(spicy_scv(i).ref(j).pmodes(:,1));

		% Gaussian model on 0:10:2000 grid along with R^2 of fit
		scv_data.spicy_scvs(i).gaussian(j).spice_anom        = single(spicy_scv(i).gauss(j).X);
		scv_data.spicy_scvs(i).gaussian(j).pres              = single(spicy_scv(i).gauss(j).Y);
		scv_data.spicy_scvs(i).gaussian(j).R2                = single(spicy_scv(i).gauss(j).R2);
		scv_data.spicy_scvs(i).gaussian(j).core_shallow_pres = single(spicy_scv(i).limits(j).core_plims(1));
		scv_data.spicy_scvs(i).gaussian(j).core_deep_pres    = single(spicy_scv(i).limits(j).core_plims(2));
		scv_data.spicy_scvs(i).gaussian(j).shallow_pres      = single(spicy_scv(i).limits(j).shallow_pres);
		scv_data.spicy_scvs(i).gaussian(j).deep_pres         = single(spicy_scv(i).limits(j).deep_pres);
		scv_data.spicy_scvs(i).gaussian(j).core_shallow_dens = single(spicy_scv(i).limits(j).core_dlims(1));
		scv_data.spicy_scvs(i).gaussian(j).core_deep_dens    = single(spicy_scv(i).limits(j).core_dlims(2));
		scv_data.spicy_scvs(i).gaussian(j).shallow_dens      = single(spicy_scv(i).limits(j).shallow_dens);
		scv_data.spicy_scvs(i).gaussian(j).deep_dens         = single(spicy_scv(i).limits(j).deep_dens);

		% Other statistics
		% Get index of core pressure
		idx = find(spicy_scv(i).pres{j} == round(spicy_scv(i).limits(j).core_pres/10)*10);
		scv_data.spicy_scvs(i).scv_stats(j).core_pres            = single(spicy_scv(i).limits(j).core_pres);
		scv_data.spicy_scvs(i).scv_stats(j).core_dens            = single(spicy_scv(i).limits(j).core_dens);
		scv_data.spicy_scvs(i).scv_stats(j).core_temp            = single(spicy_scv(i).temp{j}(idx));
		scv_data.spicy_scvs(i).scv_stats(j).core_salt            = single(spicy_scv(i).salt{j}(idx));
		scv_data.spicy_scvs(i).scv_stats(j).core_spice           = single(spicy_scv(i).spice{j}(idx));
		scv_data.spicy_scvs(i).scv_stats(j).core_N2              = single(spicy_scv(i).N2{j}(idx));
		scv_data.spicy_scvs(i).scv_stats(j).core_dyn_height_anom = single(spicy_scv(i).stats(j).Mag);
		scv_data.spicy_scvs(i).scv_stats(j).vertical_extent      = single(spicy_scv(i).stats(j).Height);
		scv_data.spicy_scvs(i).scv_stats(j).coriollis_freq       = single(spicy_scv(i).stats(j).f);
		scv_data.spicy_scvs(i).scv_stats(j).scale_height         = single(sqrt((spicy_scv(i).stats(j).Height^2)/8));
		scv_data.spicy_scvs(i).scv_stats(j).deform_radius        = single(([sqrt(spicy_scv(i).ref(j).N2(idx))/spicy_scv(i).stats(j).f])*sqrt((spicy_scv(i).stats(j).Height^2)/8)); 
	end
end

for i = 1:length(minty_scv)
	for j = 1:length(minty_scv(i).lon);

		% Metadata
		scv_data.minty_scvs(i).meta(j).ID    = minty_scv(i).ID(j);
		scv_data.minty_scvs(i).meta(j).float = minty_scv(i).float{j};
		scv_data.minty_scvs(i).meta(j).cycle = minty_scv(i).cycle(j);
		scv_data.minty_scvs(i).meta(j).lon   = minty_scv(i).lon(j);
		scv_data.minty_scvs(i).meta(j).lat   = minty_scv(i).lat(j);
		scv_data.minty_scvs(i).meta(j).time  = minty_scv(i).time{j};

		% Profile data
		scv_data.minty_scvs(i).profile(j).temp   = single(minty_scv(i).temp{j});
		scv_data.minty_scvs(i).profile(j).salt   = single(minty_scv(i).salt{j});
		scv_data.minty_scvs(i).profile(j).spice  = single(minty_scv(i).spice{j});
		scv_data.minty_scvs(i).profile(j).N2     = single(minty_scv(i).N2{j});
		scv_data.minty_scvs(i).profile(j).pres   = single(minty_scv(i).pres{j});
		scv_data.minty_scvs(i).profile(j).sigma0 = single(minty_scv(i).sigma0{j});

		% Anomaly data
		scv_data.minty_scvs(i).anomalies(j).temp  = single(minty_scv(i).temp_anom{j}); 
		scv_data.minty_scvs(i).anomalies(j).salt  = single(minty_scv(i).salt_anom{j});
		scv_data.minty_scvs(i).anomalies(j).spice = single(minty_scv(i).spice_anom{j});
		scv_data.minty_scvs(i).anomalies(j).N2    = single(minty_scv(i).N2_anom{j});
		scv_data.minty_scvs(i).anomalies(j).pres  = single(minty_scv(i).pres_anom{j});

		% IQR Thresholds
		scv_data.minty_scvs(i).thresholds(j).spice.Q1  = single(minty_scv(i).spice_limits{j}(:,1));
		scv_data.minty_scvs(i).thresholds(j).spice.Q3  = single(minty_scv(i).spice_limits{j}(:,2));
		scv_data.minty_scvs(i).thresholds(j).spice.IQR = single(minty_scv(i).spice_IQR{j});
		scv_data.minty_scvs(i).thresholds(j).N2.Q1     = single(minty_scv(i).N2_limits{j}(:,1));
		scv_data.minty_scvs(i).thresholds(j).N2.Q3     = single(minty_scv(i).N2_limits{j}(:,2));
		scv_data.minty_scvs(i).thresholds(j).N2.IQR    = single(minty_scv(i).N2_IQR{j});

		% Climatology from Scripps
		scv_data.minty_scvs(i).climatology(j).lon        = single(minty_scv(i).ref(j).lon);
		scv_data.minty_scvs(i).climatology(j).lat        = single(minty_scv(i).ref(j).lat);
		scv_data.minty_scvs(i).climatology(j).temp       = single(minty_scv(i).ref(j).temp);
		scv_data.minty_scvs(i).climatology(j).salt       = single(minty_scv(i).ref(j).salt);
		scv_data.minty_scvs(i).climatology(j).spice      = single(minty_scv(i).ref(j).spice);
		scv_data.minty_scvs(i).climatology(j).N2         = single(minty_scv(i).ref(j).N2);
		scv_data.minty_scvs(i).climatology(j).sigma0     = single(minty_scv(i).ref(j).sigma0);
		scv_data.minty_scvs(i).climatology(j).dyn_height = single(minty_scv(i).ref(j).dyn_height);
		scv_data.minty_scvs(i).climatology(j).pres       = single(minty_scv(i).ref(j).pres);

		% Vertical Mode decomposition of climatological dynamic height anomaly profile
		scv_data.minty_scvs(i).mode_decomp(j).pmodes   = single(minty_scv(i).ref(j).pmodes);
		scv_data.minty_scvs(i).mode_decomp(j).wmodes   = single(minty_scv(i).ref(j).wmodes);
		scv_data.minty_scvs(i).mode_decomp(j).alphaBC1 = single(minty_scv(i).ref(j).VMD.alpha);
		scv_data.minty_scvs(i).mode_decomp(j).pmodeBC1 = single(minty_scv(i).ref(j).pmodes(:,1));

		% Gaussian model on 0:10:2000 grid along with R^2 of fit
		scv_data.minty_scvs(i).gaussian(j).spice_anom        = single(minty_scv(i).gauss(j).X);
		scv_data.minty_scvs(i).gaussian(j).pres              = single(minty_scv(i).gauss(j).Y);
		scv_data.minty_scvs(i).gaussian(j).R2                = single(minty_scv(i).gauss(j).R2);
		scv_data.minty_scvs(i).gaussian(j).core_shallow_pres = single(minty_scv(i).limits(j).core_plims(1));
		scv_data.minty_scvs(i).gaussian(j).core_deep_pres    = single(minty_scv(i).limits(j).core_plims(2));
		scv_data.minty_scvs(i).gaussian(j).shallow_pres      = single(minty_scv(i).limits(j).shallow_pres);
		scv_data.minty_scvs(i).gaussian(j).deep_pres         = single(minty_scv(i).limits(j).deep_pres);
		scv_data.minty_scvs(i).gaussian(j).core_shallow_dens = single(minty_scv(i).limits(j).core_dlims(1));
		scv_data.minty_scvs(i).gaussian(j).core_deep_dens    = single(minty_scv(i).limits(j).core_dlims(2));
		scv_data.minty_scvs(i).gaussian(j).shallow_dens      = single(minty_scv(i).limits(j).shallow_dens);
		scv_data.minty_scvs(i).gaussian(j).deep_dens         = single(minty_scv(i).limits(j).deep_dens);

		% Other statistics
		% Get index of core pressure
		idx = find(minty_scv(i).pres{j} == round(minty_scv(i).limits(j).core_pres/10)*10);
		scv_data.minty_scvs(i).scv_stats(j).core_pres            = single(minty_scv(i).limits(j).core_pres);
		scv_data.minty_scvs(i).scv_stats(j).core_dens            = single(minty_scv(i).limits(j).core_dens);
		scv_data.minty_scvs(i).scv_stats(j).core_temp            = single(minty_scv(i).temp{j}(idx));
		scv_data.minty_scvs(i).scv_stats(j).core_salt            = single(minty_scv(i).salt{j}(idx));
		scv_data.minty_scvs(i).scv_stats(j).core_spice           = single(minty_scv(i).spice{j}(idx));
		scv_data.minty_scvs(i).scv_stats(j).core_N2              = single(minty_scv(i).N2{j}(idx));
		scv_data.minty_scvs(i).scv_stats(j).core_dyn_height_anom = single(minty_scv(i).stats(j).Mag);
		scv_data.minty_scvs(i).scv_stats(j).vertical_extent      = single(minty_scv(i).stats(j).Height);
		scv_data.minty_scvs(i).scv_stats(j).coriollis_freq       = single(minty_scv(i).stats(j).f);
		scv_data.minty_scvs(i).scv_stats(j).scale_height         = single(sqrt((minty_scv(i).stats(j).Height^2)/8));
		scv_data.minty_scvs(i).scv_stats(j).deform_radius        = single(([sqrt(minty_scv(i).ref(j).N2(idx))/minty_scv(i).stats(j).f])*sqrt((minty_scv(i).stats(j).Height^2)/8)); 
	end
end

% Save final dataset
fname = [datadir,'final_timeseries_scvs.mat'];
save(fname,'scv_data','-v7.3');
