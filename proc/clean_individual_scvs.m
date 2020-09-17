% Script to make final 'good' SCVs file for sharing
% These are individual SCVs, the time-series SCVs are in a separate file

close all ; clear all
load('../datadir.mat');

% SPICY SCVs
load([datadir,'good_spicy_scvs.mat']);
good_spicy_scvs = scv_data; clear scv_data

% MINTY SCVs
load([datadir,'good_minty_scvs.mat']);
good_minty_scvs = scv_data; clear scv_data

% Make pressure array to interpolate climatological dynamic height to regular 0:10:2000 grid
pres_array = [0:10:2000];

% Fill in dynamic height data with NaNs where no data exists (0:10:2000 grid);
% Also get modes on same grid
for i = 1:length(good_spicy_scvs)
	tmpdyn                      = nan(length(pres_array),1);
	tmppmodes		    = nan(length(pres_array),size(good_spicy_scvs(i).ref.pmodes,2));
	tmpwmodes		    = nan(length(pres_array),size(good_spicy_scvs(i).ref.pmodes,2));
	idx                         = ismember(pres_array,good_spicy_scvs(i).ref.dyn_pres);
	tmpdyn(idx==1)              = good_spicy_scvs(i).ref.dyn_height;
	good_spicy_scvs(i).ref.dyn_height = tmpdyn;
	idx                         = ismember(pres_array,good_spicy_scvs(i).ref.mode_pres);
	tmppmodes(idx==1,:)         = good_spicy_scvs(i).ref.pmodes;
	tmpwmodes(idx==1,:)         = good_spicy_scvs(i).ref.wmodes;
	good_spicy_scvs(i).ref.pmodes     = tmppmodes;
	good_spicy_scvs(i).ref.wmodes     = tmpwmodes;
end

for i = 1:length(good_minty_scvs)
	tmpdyn                      = nan(length(pres_array),1);
	tmppmodes		    = nan(length(pres_array),size(good_minty_scvs(i).ref.pmodes,2));
	tmpwmodes		    = nan(length(pres_array),size(good_minty_scvs(i).ref.pmodes,2));
	idx                         = ismember(pres_array,good_minty_scvs(i).ref.dyn_pres);
	tmpdyn(idx==1)              = good_minty_scvs(i).ref.dyn_height;
	good_minty_scvs(i).ref.dyn_height = tmpdyn;
	idx                         = ismember(pres_array,good_minty_scvs(i).ref.mode_pres);
	tmppmodes(idx==1,:)         = good_minty_scvs(i).ref.pmodes;
	tmpwmodes(idx==1,:)         = good_minty_scvs(i).ref.wmodes;
	good_minty_scvs(i).ref.pmodes     = tmppmodes;
	good_minty_scvs(i).ref.wmodes     = tmpwmodes;
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

% Clear scv_data
scv_data  = [];

for i = 1:length(good_spicy_scvs)

	% Metadata
	scv_data.spicy_scvs(i).meta.ID        = good_spicy_scvs(i).ID;
	scv_data.spicy_scvs(i).meta.float     = good_spicy_scvs(i).float;
	scv_data.spicy_scvs(i).meta.cycle     = good_spicy_scvs(i).cycle;
	scv_data.spicy_scvs(i).meta.lon       = good_spicy_scvs(i).lon;
	scv_data.spicy_scvs(i).meta.lat       = good_spicy_scvs(i).lat;
	scv_data.spicy_scvs(i).meta.time      = good_spicy_scvs(i).time;
	scv_data.spicy_scvs(i).meta.data_mode = good_spicy_scvs(i).data_mode;

	% Raw data
	scv_data.spicy_scvs(i).raw.temp = good_spicy_scvs(i).RAW_DATA.temp;
	scv_data.spicy_scvs(i).raw.salt = good_spicy_scvs(i).RAW_DATA.salt;
	scv_data.spicy_scvs(i).raw.pres = good_spicy_scvs(i).RAW_DATA.pres;

	% Profile data
	scv_data.spicy_scvs(i).profile.temp       = good_spicy_scvs(i).temp;
	scv_data.spicy_scvs(i).profile.cons_temp  = good_spicy_scvs(i).theta;
	scv_data.spicy_scvs(i).profile.salt       = good_spicy_scvs(i).salt;
	scv_data.spicy_scvs(i).profile.salt_abs   = single(good_spicy_scvs(i).salt_abs);
	scv_data.spicy_scvs(i).profile.spice      = good_spicy_scvs(i).spice;
	scv_data.spicy_scvs(i).profile.N2         = good_spicy_scvs(i).N2;
	scv_data.spicy_scvs(i).profile.dyn_height = single(good_spicy_scvs(i).dyn_height);
	scv_data.spicy_scvs(i).profile.pres       = single(good_spicy_scvs(i).pres);
	scv_data.spicy_scvs(i).profile.sigma0     = good_spicy_scvs(i).sigma0;

	% Anomaly data
	scv_data.spicy_scvs(i).anomalies.temp            = good_spicy_scvs(i).temp_anom; 
	scv_data.spicy_scvs(i).anomalies.salt            = good_spicy_scvs(i).salt_anom;
	scv_data.spicy_scvs(i).anomalies.spice           = good_spicy_scvs(i).spice_anom;
	scv_data.spicy_scvs(i).anomalies.N2              = good_spicy_scvs(i).N2_anom;
	scv_data.spicy_scvs(i).anomalies.pres            = single(good_spicy_scvs(i).pres_anom);
	scv_data.spicy_scvs(i).anomalies.dyn_height_init = single(good_spicy_scvs(i).dyn_height_anom);
	scv_data.spicy_scvs(i).anomalies.dyn_height_adj  = single(good_spicy_scvs(i).dyn_height_anom_BC1);

	% IQR Thresholds
	scv_data.spicy_scvs(i).thresholds.spice.Q1  = good_spicy_scvs(i).spice_limits(:,1);
	scv_data.spicy_scvs(i).thresholds.spice.Q3  = good_spicy_scvs(i).spice_limits(:,2);
	scv_data.spicy_scvs(i).thresholds.spice.IQR = single(good_spicy_scvs(i).spice_IQR);
	scv_data.spicy_scvs(i).thresholds.N2.Q1     = good_spicy_scvs(i).N2_limits(:,1);
	scv_data.spicy_scvs(i).thresholds.N2.Q3     = good_spicy_scvs(i).N2_limits(:,2);
	scv_data.spicy_scvs(i).thresholds.N2.IQR    = single(good_spicy_scvs(i).N2_IQR);

	% Climatology from Scripps
	scv_data.spicy_scvs(i).climatology.lon        = single(good_spicy_scvs(i).ref.lon);
	scv_data.spicy_scvs(i).climatology.lat        = single(good_spicy_scvs(i).ref.lat);
	scv_data.spicy_scvs(i).climatology.temp       = single(good_spicy_scvs(i).ref.temp);
	scv_data.spicy_scvs(i).climatology.salt       = single(good_spicy_scvs(i).ref.salt);
	scv_data.spicy_scvs(i).climatology.spice      = single(good_spicy_scvs(i).ref.spice);
	scv_data.spicy_scvs(i).climatology.N2         = single(good_spicy_scvs(i).ref.N2);
	scv_data.spicy_scvs(i).climatology.sigma0     = single(good_spicy_scvs(i).ref.sigma0);
	scv_data.spicy_scvs(i).climatology.dyn_height = single(good_spicy_scvs(i).ref.dyn_height);
	scv_data.spicy_scvs(i).climatology.pres       = single(good_spicy_scvs(i).ref.pres);

	% Vertical Mode decomposition of climatological dynamic height anomaly profile
	scv_data.spicy_scvs(i).mode_decomp.pmodes     = single(good_spicy_scvs(i).ref.pmodes);
	scv_data.spicy_scvs(i).mode_decomp.wmodes     = single(good_spicy_scvs(i).ref.wmodes);
	scv_data.spicy_scvs(i).mode_decomp.alphaBC1   = single(good_spicy_scvs(i).ref.VMD.alpha);
	scv_data.spicy_scvs(i).mode_decomp.pmodeBC1   = single(good_spicy_scvs(i).ref.pmodes(:,1));

	% Gaussian model on 0:10:2000 grid along with R^2 of fit
	scv_data.spicy_scvs(i).gaussian.spice_anom        = single(good_spicy_scvs(i).gauss.X);
	scv_data.spicy_scvs(i).gaussian.pres              = single(good_spicy_scvs(i).gauss.Y);
	scv_data.spicy_scvs(i).gaussian.R2                = single(good_spicy_scvs(i).gauss.R2);
	scv_data.spicy_scvs(i).gaussian.core_shallow_pres = single(good_spicy_scvs(i).limits.core_plims(1));
	scv_data.spicy_scvs(i).gaussian.core_deep_pres    = single(good_spicy_scvs(i).limits.core_plims(2));
	scv_data.spicy_scvs(i).gaussian.shallow_pres      = single(good_spicy_scvs(i).limits.shallow_pres);
	scv_data.spicy_scvs(i).gaussian.deep_pres         = single(good_spicy_scvs(i).limits.deep_pres);
	scv_data.spicy_scvs(i).gaussian.core_shallow_dens = single(good_spicy_scvs(i).limits.core_dlims(1));
	scv_data.spicy_scvs(i).gaussian.core_deep_dens    = single(good_spicy_scvs(i).limits.core_dlims(2));
	scv_data.spicy_scvs(i).gaussian.shallow_dens      = single(good_spicy_scvs(i).limits.shallow_dens);
	scv_data.spicy_scvs(i).gaussian.deep_dens         = single(good_spicy_scvs(i).limits.deep_dens);

	% Other statistics
	% Get index of core pressure
	idx = find(good_spicy_scvs(i).pres == round(good_spicy_scvs(i).limits.core_pres/10)*10);
	scv_data.spicy_scvs(i).scv_stats.core_pres            = single(good_spicy_scvs(i).limits.core_pres);
	scv_data.spicy_scvs(i).scv_stats.core_dens            = single(good_spicy_scvs(i).limits.core_dens);
	scv_data.spicy_scvs(i).scv_stats.core_temp            = single(good_spicy_scvs(i).temp(idx));
	scv_data.spicy_scvs(i).scv_stats.core_salt            = single(good_spicy_scvs(i).salt(idx));
	scv_data.spicy_scvs(i).scv_stats.core_spice           = single(good_spicy_scvs(i).spice(idx));
	scv_data.spicy_scvs(i).scv_stats.core_N2              = single(good_spicy_scvs(i).N2(idx));
	scv_data.spicy_scvs(i).scv_stats.core_dyn_height_anom = single(good_spicy_scvs(i).stats.Mag);
	scv_data.spicy_scvs(i).scv_stats.vertical_extent      = single(good_spicy_scvs(i).stats.Height);
	scv_data.spicy_scvs(i).scv_stats.coriollis_freq       = single(good_spicy_scvs(i).stats.f);
	scv_data.spicy_scvs(i).scv_stats.scale_height         = single(sqrt((good_spicy_scvs(i).stats.Height^2)/8));
	scv_data.spicy_scvs(i).scv_stats.deform_radius        = single(([sqrt(good_spicy_scvs(i).ref.N2(idx))/good_spicy_scvs(i).stats.f])*sqrt((good_spicy_scvs(i).stats.Height^2)/8)); 
end

for i = 1:length(good_minty_scvs)

	% Metadata
	scv_data.minty_scvs(i).meta.ID        = good_minty_scvs(i).ID;
	scv_data.minty_scvs(i).meta.float     = good_minty_scvs(i).float;
	scv_data.minty_scvs(i).meta.cycle     = good_minty_scvs(i).cycle;
	scv_data.minty_scvs(i).meta.lon       = good_minty_scvs(i).lon;
	scv_data.minty_scvs(i).meta.lat       = good_minty_scvs(i).lat;
	scv_data.minty_scvs(i).meta.time      = good_minty_scvs(i).time;
	scv_data.minty_scvs(i).meta.data_mode = good_minty_scvs(i).data_mode;

	% Raw data
	scv_data.minty_scvs(i).raw.temp = good_minty_scvs(i).RAW_DATA.temp;
	scv_data.minty_scvs(i).raw.salt = good_minty_scvs(i).RAW_DATA.salt;
	scv_data.minty_scvs(i).raw.pres = good_minty_scvs(i).RAW_DATA.pres;

	% Profile data
	scv_data.minty_scvs(i).profile.temp       = good_minty_scvs(i).temp;
	scv_data.minty_scvs(i).profile.cons_temp  = good_minty_scvs(i).theta;
	scv_data.minty_scvs(i).profile.salt       = good_minty_scvs(i).salt;
	scv_data.minty_scvs(i).profile.salt_abs   = single(good_minty_scvs(i).salt_abs);
	scv_data.minty_scvs(i).profile.spice      = good_minty_scvs(i).spice;
	scv_data.minty_scvs(i).profile.N2         = good_minty_scvs(i).N2;
        scv_data.minty_scvs(i).profile.dyn_height = single(good_minty_scvs(i).dyn_height);
	scv_data.minty_scvs(i).profile.pres       = single(good_minty_scvs(i).pres);
	scv_data.minty_scvs(i).profile.sigma0     = good_minty_scvs(i).sigma0;

	% Anomaly data
	scv_data.minty_scvs(i).anomalies.temp  	      = good_minty_scvs(i).temp_anom; 
	scv_data.minty_scvs(i).anomalies.salt  	      = good_minty_scvs(i).salt_anom;
	scv_data.minty_scvs(i).anomalies.spice 	      = good_minty_scvs(i).spice_anom;
	scv_data.minty_scvs(i).anomalies.N2    	      = good_minty_scvs(i).N2_anom;
	scv_data.minty_scvs(i).anomalies.pres            = single(good_minty_scvs(i).pres_anom);
	scv_data.minty_scvs(i).anomalies.dyn_height_init = single(good_minty_scvs(i).dyn_height_anom);
	scv_data.minty_scvs(i).anomalies.dyn_height_adj  = single(good_minty_scvs(i).dyn_height_anom_BC1);

	% IQR Thresholds
	scv_data.minty_scvs(i).thresholds.spice.Q1  = good_minty_scvs(i).spice_limits(:,1);
	scv_data.minty_scvs(i).thresholds.spice.Q3  = good_minty_scvs(i).spice_limits(:,2);
	scv_data.minty_scvs(i).thresholds.spice.IQR = single(good_minty_scvs(i).spice_IQR);
	scv_data.minty_scvs(i).thresholds.N2.Q1     = good_minty_scvs(i).N2_limits(:,1);
	scv_data.minty_scvs(i).thresholds.N2.Q3     = good_minty_scvs(i).N2_limits(:,2);
	scv_data.minty_scvs(i).thresholds.N2.IQR    = single(good_minty_scvs(i).N2_IQR);

	% Climatology from Scripps
	scv_data.minty_scvs(i).climatology.lon        = single(good_minty_scvs(i).ref.lon);
	scv_data.minty_scvs(i).climatology.lat        = single(good_minty_scvs(i).ref.lat);
	scv_data.minty_scvs(i).climatology.temp       = single(good_minty_scvs(i).ref.temp);
	scv_data.minty_scvs(i).climatology.salt       = single(good_minty_scvs(i).ref.salt);
	scv_data.minty_scvs(i).climatology.spice      = single(good_minty_scvs(i).ref.spice);
	scv_data.minty_scvs(i).climatology.N2         = single(good_minty_scvs(i).ref.N2);
	scv_data.minty_scvs(i).climatology.sigma0     = single(good_minty_scvs(i).ref.sigma0);
	scv_data.minty_scvs(i).climatology.dyn_height = single(good_minty_scvs(i).ref.dyn_height);
	scv_data.minty_scvs(i).climatology.pres       = single(good_minty_scvs(i).ref.pres);

	% Vertical Mode decomposition of climatological dynamic height anomaly profile
	scv_data.minty_scvs(i).mode_decomp.pmodes     = single(good_minty_scvs(i).ref.pmodes);
	scv_data.minty_scvs(i).mode_decomp.wmodes     = single(good_minty_scvs(i).ref.wmodes);
	scv_data.minty_scvs(i).mode_decomp.alphaBC1   = single(good_minty_scvs(i).ref.VMD.alpha);
	scv_data.minty_scvs(i).mode_decomp.pmodeBC1   = single(good_minty_scvs(i).ref.pmodes(:,1));

	% Gaussian model on 0:10:2000 grid along with R^2 of fit
	scv_data.minty_scvs(i).gaussian.spice_anom        = single(good_minty_scvs(i).gauss.X);
	scv_data.minty_scvs(i).gaussian.pres              = single(good_minty_scvs(i).gauss.Y);
	scv_data.minty_scvs(i).gaussian.R2                = single(good_minty_scvs(i).gauss.R2);
	scv_data.minty_scvs(i).gaussian.core_shallow_pres = single(good_minty_scvs(i).limits.core_plims(1));
	scv_data.minty_scvs(i).gaussian.core_deep_pres    = single(good_minty_scvs(i).limits.core_plims(2));
	scv_data.minty_scvs(i).gaussian.shallow_pres      = single(good_minty_scvs(i).limits.shallow_pres);
	scv_data.minty_scvs(i).gaussian.deep_pres         = single(good_minty_scvs(i).limits.deep_pres);
	scv_data.minty_scvs(i).gaussian.core_shallow_dens = single(good_minty_scvs(i).limits.core_dlims(1));
	scv_data.minty_scvs(i).gaussian.core_deep_dens    = single(good_minty_scvs(i).limits.core_dlims(2));
	scv_data.minty_scvs(i).gaussian.shallow_dens      = single(good_minty_scvs(i).limits.shallow_dens);
	scv_data.minty_scvs(i).gaussian.deep_dens         = single(good_minty_scvs(i).limits.deep_dens);

	% Other statistics
	% Get index of core pressure
	idx = find(good_minty_scvs(i).pres == round(good_minty_scvs(i).limits.core_pres/10)*10);
	scv_data.minty_scvs(i).scv_stats.core_pres            = single(good_minty_scvs(i).limits.core_pres);
	scv_data.minty_scvs(i).scv_stats.core_dens            = single(good_minty_scvs(i).limits.core_dens);
	scv_data.minty_scvs(i).scv_stats.core_temp            = single(good_minty_scvs(i).temp(idx));
	scv_data.minty_scvs(i).scv_stats.core_salt            = single(good_minty_scvs(i).salt(idx));
	scv_data.minty_scvs(i).scv_stats.core_spice           = single(good_minty_scvs(i).spice(idx));
	scv_data.minty_scvs(i).scv_stats.core_N2              = single(good_minty_scvs(i).N2(idx));
	scv_data.minty_scvs(i).scv_stats.core_dyn_height_anom = single(good_minty_scvs(i).stats.Mag);
	scv_data.minty_scvs(i).scv_stats.vertical_extent      = single(good_minty_scvs(i).stats.Height);
	scv_data.minty_scvs(i).scv_stats.coriollis_freq       = single(good_minty_scvs(i).stats.f);
	scv_data.minty_scvs(i).scv_stats.scale_height         = single(sqrt((good_minty_scvs(i).stats.Height^2)/8));
	scv_data.minty_scvs(i).scv_stats.deform_radius        = single(([sqrt(good_minty_scvs(i).ref.N2(idx))/good_minty_scvs(i).stats.f])*sqrt((good_minty_scvs(i).stats.Height^2)/8)); 
end

% Save final dataset
fname = [datadir,'final_individual_scvs.mat'];
save(fname,'scv_data','-v7.3');
