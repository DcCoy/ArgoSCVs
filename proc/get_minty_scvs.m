%// Script to analyze initial SCV detections by:
%// 1. Performing some QC (removing bad profiles)
%// 2. Fitting a gaussian model to spice anomalies in pressure space
%// 3. Calculating dynamic height starting from the deepest pressure
%// 4. Performing a vertical mode decomposition on the reference cast
%// 5. Fitting the first baroclinic horizontal velocity mode to dyn height anomalies
%// 6. Characterizing the structure of detected SCVs (height, max anomaly, etc).
%// NOTES:
%// FLAGGING PROFILES USING THE BELOW CRITERIA
%// 1. T/S/Spice anomaly profile must not be 100% outside IQR thresholds
%// 2. Spice threshold can't only be exceeded at the surface
%// 3. Reference cast must go to 1000dbar
%// 4. Gaussian model must succeed
%//    - Able to fit a peak in spiciness anomaly that coincides with negative N2 passing threshold (+- 0.1 kg/m^3)
%//    - Minimum peak in spiciness == 0.1 kg/m3 (visual inspection of good profiles)
%//    - Predicted H between 100 - 1200 dbar
%//    - Fit must be at least r^2 of 0.5 and NRMSE of < 0.5
%//    - Must have data at the core limits of the Gaussian
%// 7. Must be a maxima in dynamic height anomaly between limits

close all ; clear all
load('../datadir.mat');

% Addpath to gsw scripts
addpath scripts/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MINTYY SCVs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Load quality controlled SCV detections
%// Checked for huge salt/temp anomalies, rejects bad profiles
disp('Loading MINTY SCVs')
load([datadir,'initial_minty_scv_qc.mat']);
scv_data = minty_scvs;
choice   = 2;
clearvars -except scv_data choice datadir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Initiate flags for bad detections
disp('Initiating matrices to fill')
bad_lat       = zeros(1,length(scv_data)); %// FLAG #1
bad_model     = zeros(1,length(scv_data)); %// FLAG #2
bad_dh_peak   = zeros(1,length(scv_data)); %// FLAG #3
bad_dh_lims   = zeros(1,length(scv_data)); %// FLAG #4
bad_dh_points = zeros(1,length(scv_data)); %// FLAG #5
bad_other     = zeros(1,length(scv_data)); %// FLAG #6
good_scv      = zeros(1,length(scv_data)); %// FLAG #00

%// Initiate some constant-sized fields to fill
%// Final structure will be cleaned at end of the file
for i = 1:length(scv_data)

    %// Dynamics
    scv_data(i).salt_abs            = nan(size(scv_data(i).temp));
    scv_data(i).theta               = nan(size(scv_data(i).temp));
    scv_data(i).dyn_height          = [];
    scv_data(i).dyn_height_anom     = [];
    scv_data(i).dyn_height_anom_BC1 = [];
    scv_data(i).dyn_height_pres_BC1 = [];
    scv_data(i).dyn_pres            = [];
    
    %// Stats
    scv_data(i).stats.f             = [];
    scv_data(i).stats.N             = [];
    scv_data(i).stats.Height        = [];
    scv_data(i).stats.Peak          = [];
    scv_data(i).stats.Rd            = [];
    scv_data(i).stats.Vol           = [];
    scv_data(i).stats.Mag           = [];

    %// SCV limits
    scv_data(i).limits.core_pres    = [];
    scv_data(i).limits.core_plims   = [];    
    scv_data(i).limits.shallow_pres = [];
    scv_data(i).limits.deep_pres    = [];
    scv_data(i).limits.core_dens    = [];
    scv_data(i).limits.core_dlims   = [];  
    scv_data(i).limits.shallow_dens = [];
    scv_data(i).limits.deep_dens    = [];

    %// Reference fields
    scv_data(i).ref.dyn_height      = [];
    scv_data(i).ref.dyn_pres        = [];
    scv_data(i).ref.dyn_N2          = [];
    scv_data(i).ref.pmodes          = [];
    scv_data(i).ref.wmodes          = [];
    scv_data(i).ref.mode_pres       = [];
    scv_data(i).ref.BC1_data        = [];
    scv_data(i).ref.BC1_pres        = [];

    %// Model fields
    scv_data(i).gauss.Y             = [];
    scv_data(i).gauss.X             = [];
    scv_data(i).gauss.R2            = [];
    scv_data(i).gauss.NRMSE         = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start parallel proc')
%// Final prep before for-loop
%// No display, parallel on for optimization routine
opts1    = optimset('display','off','UseParallel',true);

%// Check progress matrix
progress = [0:1000:length(scv_data)];

%// Initiate parpool
delete(gcp('nocreate'))
parpool(12)

%// Go through each float group and get model, DH' for SCVs
parfor i = 1:length(scv_data)
try
    %// Show progress
    if ismember(i,progress)==1
	disp(i)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 1         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Reject SCVs near equator
    if abs(scv_data(i).lat)<5
	bad_lat(i) = 1;
	continue
    end

    %// Rename reference for easier to read script
    ref = [];
    ref = scv_data(i).ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Attempt to fit a gaussian structure to the spice anomaly structure
    [gauss_model] = [];
    [gauss_stats] = [];
    [gauss_lims]  = [];
    [gauss_model gauss_stats gauss_lims] = fit_gaussian(scv_data(i),choice);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 2         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Reject failed profiles
    if gauss_model.x == 999 | gauss_model.y == 999
	bad_model(i) = 1;
	continue
    end    

    %// Grab stats for plotting and saving for later
    fity     = []; fity     = gauss_model.y;
    fitx     = []; fitx     = gauss_model.x;
    fitR2    = []; fitR2    = gauss_stats.R2;
    fitNRMSE = []; fitNRMSE = gauss_stats.NRMSE;
    mn       = []; mn       = gauss_lims.p_peak;
    liml     = []; liml     = gauss_lims.p_lims(1);
    limh     = []; limh     = gauss_lims.p_lims(2);
    scv_data(i).gauss.Y     = fity;
    scv_data(i).gauss.X     = fitx;
    scv_data(i).gauss.R2    = fitR2;
    scv_data(i).gauss.NRMSE = fitNRMSE;    

    %// Find pressure index of the limits
    %// Lower limit
    indl = [];
    indl = find(abs([liml - scv_data(i).pres]) == min(abs([liml - scv_data(i).pres])));
    if length(indl) > 1
       indl = indl(end);
    end

    %// Upper limit
    indh = [];
    indh = find(abs([limh - scv_data(i).pres]) == min(abs([limh - scv_data(i).pres])));
    if length(indh) > 1
	indh = indh(1);
    end

    %// Peak
    indm = [];
    indm = find(abs([mn - scv_data(i).pres]) == min(abs([mn - scv_data(i).pres])));
    if length(indm) > 1
	indm = indm(1);
    end

%// Save limit levels
    scv_data(i).limits.core_pres    = gauss_lims.p_peak;
    scv_data(i).limits.core_plims   = gauss_lims.p_clims;    
    scv_data(i).limits.shallow_pres = gauss_lims.p_lims(1);
    scv_data(i).limits.deep_pres    = gauss_lims.p_lims(2);
    scv_data(i).limits.core_dens    = gauss_lims.d_peak;
    scv_data(i).limits.core_dlims   = gauss_lims.d_clims;
    scv_data(i).limits.shallow_dens = gauss_lims.d_lims(1);
    scv_data(i).limits.deep_dens    = gauss_lims.d_lims(2);

%// Define easier to use limit levels
    pl  = gauss_lims.p_lims(1);  %// lower limit (pressure)
    plc = gauss_lims.p_clims(1); %// lower core limit (pressure)
    ph  = gauss_lims.p_lims(2);  %// upper limit (pressure)
    phc = gauss_lims.p_clims(2); %// upper core limit (pressure)
    pm  = gauss_lims.p_peak;     %// core pressure
    dl  = gauss_lims.d_lims(1);  %// lower limit (density)
    dh  = gauss_lims.d_lims(2);  %// upper limit (density)
    dm  = gauss_lims.d_peak;     %// core density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Calculate dynamic height for both profiles starting from deepest p
    %// Add in salt abs and theta for SCV profiles
    scv_data(i).salt_abs = gsw_SA_from_SP(scv_data(i).salt,scv_data(i).pres,scv_data(i).lon,scv_data(i).lat);
    scv_data(i).theta    = gsw_CT_from_t(scv_data(i).salt_abs,scv_data(i).temp,scv_data(i).pres);

    %// Remove NaNs for dyn_height calculation
    ps  = [];
    ps  = scv_data(i).pres(~isnan([scv_data(i).salt_abs + scv_data(i).theta]));
    cts = [];
    cts = scv_data(i).theta(~isnan([scv_data(i).salt_abs + scv_data(i).theta]));
    sas = [];
    sas = scv_data(i).salt_abs(~isnan([scv_data(i).salt_abs + scv_data(i).theta])); 
    pr  = [];
    pr  = ref.pres(~isnan([ref.salt_abs + ref.theta]));
    ctr = [];
    ctr = ref.theta(~isnan([ref.salt_abs + ref.theta]));
    sar = [];
    sar = ref.salt_abs(~isnan([ref.salt_abs + ref.theta]));
    n2r = [];
    n2r = ref.N2(~isnan([ref.salt_abs + ref.theta]));
	
    %// Use intersect to get max comparable density, limit matrices to identical levels
    psr  = [];
    psr  = intersect(ps,pr);
    inds = [];
    inds = find(ismember(ps,psr)==1);
    indr = [];
    indr = find(ismember(pr,psr)==1);
    cts  = cts(inds); sas = sas(inds); ps = ps(inds);
    ctr  = ctr(indr); sar = sar(indr); pr = pr(indr); n2r = n2r(indr);

    %// Calculate dynamic height for both casts starting at deepest point
    scv_data(i).dyn_height = gsw_geo_strf_dyn_height(sas,cts,ps,max(psr));
    scv_data(i).dyn_pres   = ps;
    ref.dyn_height         = gsw_geo_strf_dyn_height(sar,ctr,pr,max(psr));
    ref.dyn_pres           = pr;
    ref.dyn_N2             = n2r;

    %// Define dynamic height anomaly along isobars
    scv_data(i).dyn_height_anom = scv_data(i).dyn_height - ref.dyn_height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Decompose reference profile into vertical/horizontal modes
    %// Perform decomposition
    [ref.wmodes,ref.pmodes,~,~] = dynmodes(ref.N2(~isnan(ref.N2)),ref.pres(~isnan(ref.N2)),1);
	
    %// Grab pressure levels of decomposition, add zero level (surface)
    ref.mode_pres = ref.pres(~isnan(ref.N2));	
    ref.mode_pres = [0;ref.mode_pres];

    %// Interpolate 1st baroclinic mode to pressure of SCV cast
    ref.BC1_data = ref.pmodes(:,1);
    ref.BC1_pres = ref.mode_pres;
    ref.BC1_data = interp1(ref.mode_pres(~isnan(ref.BC1_data)),ref.BC1_data(~isnan(ref.BC1_data)),scv_data(i).dyn_pres);
    ref.BC1_pres = scv_data(i).dyn_pres;

    %// Create function that describes residuals between projected BC1 and dyn_height_anom
    %// Exclude data inbetween SCV limits for better fit to first mode
    dat = [];
    dat = [ref.BC1_data + scv_data(i).dyn_height_anom];
    x_o = [];
    x_o = ref.BC1_data(~isnan(dat));
    x_p = [];
    x_p = ref.BC1_pres(~isnan(dat));
    x_f = [];
    x_f = scv_data(i).dyn_height_anom(~isnan(dat));

    %// Remove values between upper/lower limits of SCV to avoid bad fit
    ind = [];
    ind = find(pl < x_p & x_p < ph);
    if ind(end) == length(x_p)
	ind = ind(1:end-1);
    end
    x_o(ind) = [];
    x_f(ind) = [];

    %// Remove mixed layer depths (Lynne Talley method, first density greater than 0.03 from sfc value
    ind      = [];
    mld_dens = scv_data(i).sigma0(~isnan(scv_data(i).sigma0));
    mld_pres = scv_data(i).pres(~isnan(scv_data(i).sigma0));
    ind      = find(mld_dens > mld_dens(1)+0.03);
    mld_pres = mld_pres(ind(1));
    ind      = find(x_p < mld_pres);
    x_o(ind) = [];
    x_f(ind) = [];

    %// f simply evaluates a given alpha (modal amplitude) and returns the 
    %// difference between the input DHanom profile and the projected 1st mode
    %// We want to restrict our solutions such that the bottom of the projected
    %// profile is equal to the bottom of the DHanom profile
    %// SO let alpha2 = DHanom(end) - alpha*BT1(end)
    f = [];
    f = @(alpha) (alpha*x_o - x_f + (x_f(end) - alpha*x_o(end))); 
    x0  = 0.05; %// First guess

    %// Solve for best modal amplitude
    alpha = [];
    alpha = lsqnonlin(f,x0,[-1],[1],opts1);

    %// Redfine x_o and x_f with full profile
    x_o = ref.BC1_data(~isnan(dat));
    x_p = ref.mode_pres(~isnan(dat));
    x_f = scv_data(i).dyn_height_anom(~isnan(dat));

    %// Fix dynamic height anomaly by removing projected 1st mode, add back in barotopic mode
    scv_data(i).dyn_height_anom_BC1 = [x_f] - [x_o*alpha + (x_f(end) - alpha*x_o(end))]; 
    scv_data(i).dyn_height_pres_BC1 = scv_data(i).dyn_pres(~isnan(dat));

    %// Save VMD results
    ref.VMD.x_f      = x_f;
    ref.VMD.x_o      = x_o;
    ref.VMD.x_p      = x_p;
    ref.VMD.alpha    = alpha;
    scv_data(i).ref  = ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 3         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Check for 'peak' in dynamic height around core of SCV
    plidx  = find(abs([pl-scv_data(i).dyn_height_pres_BC1]) == min(abs([pl-scv_data(i).dyn_height_pres_BC1])));
    phidx  = find(abs([ph-scv_data(i).dyn_height_pres_BC1]) == min(abs([ph-scv_data(i).dyn_height_pres_BC1])));
    scvidx = find(pl <= scv_data(i).dyn_height_pres_BC1 & scv_data(i).dyn_height_pres_BC1 <= ph);
    dh_low   = scv_data(i).dyn_height_anom_BC1(plidx);
    dh_high  = scv_data(i).dyn_height_anom_BC1(phidx);    
    dh_peak  = max(scv_data(i).dyn_height_anom_BC1(scvidx));
    if isempty([dh_low+dh_high+dh_peak])==1 | dh_low >= dh_peak | dh_high >= dh_peak
	bad_dh_peak(i) = 1;
	continue
    end
    
%// Find core of SCV in density space in order to calculate background N^2
    indpkr = [];
    indpkr = find(abs([pm-ref.dyn_pres]) == min(abs([pm-ref.dyn_pres])));
    if length(indpkr)>1
	indpkr = indpkr(end);
    end
	
%// Generate statistics for each SCV
    scv_data(i).stats.f      = gsw_f(scv_data(i).lat);
    scv_data(i).stats.N      = sqrt(ref.dyn_N2(indpkr));
    scv_data(i).stats.Height = ph - pl;
    scv_data(i).stats.Peak   = gauss_lims.d_peak;
    scv_data(i).stats.Rd     = [scv_data(i).stats.N*scv_data(i).stats.Height]/scv_data(i).stats.f;
    scv_data(i).stats.Vol    = (4/3)*pi*scv_data(i).stats.Height*[scv_data(i).stats.Rd^2];
    scv_data(i).stats.Mag    = dh_peak;

%// Assuming all went well, generate plots and move onto next SCV
    good_scv(i) = 1;
    close all
	
catch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 6         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bad_other(i) = 1;
end
end

%// Kill parpool (if any)
delete(gcp('nocreate'))

%// Keep bad scvs
scv_orig       = scv_data;
bad_lat        = scv_data(find(bad_lat==1));
bad_dh_peak    = scv_data(find(bad_dh_peak==1));
bad_dh_lims    = scv_data(find(bad_dh_lims==1));
bad_dh_points  = scv_data(find(bad_dh_points==1));
bad_model      = scv_data(find(bad_model==1));
bad_other      = scv_data(find(bad_other==1));

%// Save good scvs
scv_data = scv_data(find(good_scv==1));

%// Remove duplicate profiles
%// Grab float,cast,date,location and ID
for i = 1:length(scv_data)
    floatnum(i) = str2num(scv_data(i).float);
    castnum(i)  = scv_data(i).cycle;
    ID(i)       = scv_data(i).ID;
end

%// Duplicate IDs
[c,ia,ic]       = unique(ID);
ID              = ID(ia);
castnum         = castnum(ia);
floatnum        = floatnum(ia);
scv_data        = scv_data(ia);

%// Duplicate floatnum/castnum combo
[c,ia,ic]       = unique(floatnum);
floatlist       = floatnum(ia);
good_idx        = [];
for i = 1:length(floatlist)
    ind         = find(floatnum==floatlist(i));
    casts       = castnum(ind);
    [c,ia,ic]   = unique(casts);
    good_idx    = [good_idx ind(ia)];
end
scv_data = scv_data(good_idx);

%// Fix dynamic height info
for i = 1:length(scv_data)
    ind                             = ismember(scv_data(i).pres,scv_data(i).dyn_pres);
    dh                              = nan(size(scv_data(i).pres));
    dh(ind==1)                      = scv_data(i).dyn_height;
    scv_data(i).dyn_height          = dh;
    dh                              = nan(size(scv_data(i).pres));
    dh(ind==1)                      = scv_data(i).dyn_height_anom;
    scv_data(i).dyn_height_anom     = dh;
    ind                             = ismember(scv_data(i).pres,scv_data(i).dyn_height_pres_BC1);
    dh                              = nan(size(scv_data(i).pres));
    dh(ind==1)                      = scv_data(i).dyn_height_anom_BC1;
    scv_data(i).dyn_height_anom_BC1 = dh;
end
scv_data = rmfield(scv_data,{'dyn_pres','dyn_height_pres_BC1'});

disp(['# of SCVs: ',num2str(length(scv_data))])
pause(2)
disp('...')

%// Save 'good' SCVs (passed gaussian and DH' tests)
save([datadir,'good_minty_scvs'],'scv_data','-v7.3');

