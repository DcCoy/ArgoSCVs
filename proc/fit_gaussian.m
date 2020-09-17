function [gauss_model gauss_stats gauss_lims] = fit_gaussian(scv_data,choice)
%// This script will look for instances of spice/N2 anomaly exceeding 
%// interquartile range thresholds. It will then attempt to fit a 1-term
%// gaussian model to the spice anomalies in pressure space. The function
%// returns the model, the model limits in pressure/density space (including core)
%// and the goodness-of-fit stats

%// FLAGS:
%// FLAG-1: No identifiable peaks in spiciness anomaly
%// FLAG-2: No matching IQR signals in spice/N2 within 0.1 kg/m3
%// FLAG-3: Bad goodness-of-fit results 
%//         R2    > 0.5
%//         NRMSE < 0.5 
%// FLAG-4: SCV height below 150m or above 1200m
%// FLAG-5: No data at SCV limits
%// FLAG-6: No N2 signal within SCV core

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// For spicy, keep things how they are
%// For minty, flip the sign of spice anomaly and thresholds and run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Flags trigger a return
%// Initiate flag values
%// Flag values (code calling this function will reject profile)
gauss_model.x      = 999;
gauss_model.y      = 999;
gauss_model.A      = 999;
gauss_model.zo     = 999;
gauss_model.h      = 999;
gauss_model.ps     = 999;
gauss_stats.R2     = 999;
gauss_stats.NRMSE  = 999;
gauss_lims.p_lims  = 999;
gauss_lims.p_clims = 999;
gauss_lims.p_peak  = 999;
gauss_lims.d_lims  = 999;
gauss_lims.d_clims = 999;
gauss_lims.d_peak  = 999;

%// Goodness-of-fit thresholds
R2_thresh    = 0.5;
NRMSE_thresh = 0.5;

try    
    %// If minty, swap spice anom and thresholds (easy fix for code to run fine)
    if choice == 2
	scv_data.spice_anom   = -scv_data.spice_anom;
	scv_data.spice_limits = -scv_data.spice_limits;	
    end

    %// Find data below mid-pycnocline density
    ind = find(scv_data.sigma0 >= scv_data.pyc_dens);

    %// Limit data to below mid-pycnocline density
    pres      = scv_data.pres(ind);
    dens      = scv_data.sigma0(ind);
    sp_anom   = scv_data.spice_anom(ind); %smooth(scv_data.spice_anom(ind));
    n2_anom   = scv_data.N2_anom(ind);
    sp_thresh = scv_data.spice_limits(ind,:);
    sp_IQR    = scv_data.spice_IQR(ind,:);
    n2_thresh = scv_data.N2_limits(ind,:);
    n2_IQR    = scv_data.N2_IQR(ind,:);

    %// Find peaks of smoothed data
    [init_pks,init_locs] = findpeaks2(double(sp_anom));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 1         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag profiles with no initial peaks
    if isempty(init_locs)==1
	return
    end

    %// Check that peak in spiciness matches N2 anomaly signal
    %// Find where spice and N2 anomaly exceeds IQR threshold
    if choice == 1 %// spicy
	sp_exceed = sp_anom - [sp_thresh(:,2)+1.5*sp_IQR];
    elseif choice == 2 %// minty
	sp_exceed = sp_anom - [sp_thresh(:,1)+1.5*sp_IQR];
    end
    n2_exceed = [n2_thresh(:,1)-1.5*n2_IQR] - n2_anom;

    %// Find index of spice anomaly above IQR
    eind = find(sp_exceed>0);

    %// Test that each peak in spice anomaly is above IQR
    test           = ismember(dens(init_locs),dens(eind));
    good_init_locs = init_locs(test>0); %// where it is!

    %// Now look through each peak and check that N2 anomaly test is also triggered
    %// Also check that peak spice anomaly is greater than 0.10 kg/m3
    bad_locs = [];
    for ii = 1:length(good_init_locs)
	peakdens     = dens(good_init_locs(ii));
	peakspice    = sp_anom(good_init_locs(ii));
	peakiqr      = sp_IQR(good_init_locs(ii));
	if choice == 1
	    peaklim = sp_thresh(good_init_locs(ii),2);
	else
	    peaklim = sp_thresh(good_init_locs(ii),1);
	end
	ind          = find(peakdens-0.2 <= dens & dens <= peakdens+0.2); %changed from 0.1
	n2_check     = n2_exceed(ind);
	n2_check_ind = find(n2_check>0);
	%// Reject weak peaks
	if isempty(n2_check_ind)==1 | peakspice <= 0.10 | peakspice < [peaklim + 1.5*peakiqr]
	    bad_locs = [bad_locs ii];
	end
    end
    good_init_locs(bad_locs) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 2         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag profiles with no good peaks above IQR
    if sum(test) == 0 | isempty(good_init_locs)==1
	return
    end

    %// Select the larger peaks (if more than one exists)
    if length(good_init_locs)>1
	sigs       = dens(good_init_locs);
	anoms      = sp_anom(good_init_locs);
	indr       = find(anoms == max(anoms));
	sig_spike  = sigs(indr);
	anom_spike = anoms(indr);
    else
	sig_spike  = dens(good_init_locs);
	anom_spike = sp_anom(good_init_locs);
    end

    %// Find pressure of spike
    ind        = find(dens == sig_spike);
    if length(ind) > 1
	pres_spike = pres(ind(end));
    else
	pres_spike = pres(ind);
    end

    %// Set Gaussian magnitude == max above IQR
    A    = double(anom_spike);

    %// Find best H according to least squares regression
    %// Allow height to vary by 50:10:850 (roughly 150 - 1200 dbar)
    %// Allow A to vary by +- 20%
    %// Allow core pressure to range by +- 20% of height
    pdrng = [-0.2:0.05:0.2];
    arng  = [0.8:0.05:1.2];
    hrng  = [50:10:850];
    lse   = nan(length(pdrng),length(arng),length(hrng)); 
    R2    = nan(size(lse));
    NRMSE = nan(size(lse));

    hcnt = 0;
    for h = [50:10:850]
	hcnt = hcnt + 1; 

	acnt = 0;
	for ad = [0.8:0.05:1.2]
	    acnt = acnt + 1;
	  
	    pdcnt = 0;
	    for pd = [-0.2:0.05:0.2]
		pdcnt = pdcnt + 1;

		%// Get anomaly profile to compare model to
		if choice == 1
		    sa   = double(scv_data.spice_anom);
		elseif choice == 2
		    sa   = double(-scv_data.spice_anom);
		end

		%// Center Gaussian model  around pres_spike + pd*Height
		%// Where Height = 4*sqrt(h^2/2) --> Half height = 2*sqrt(h^2/2);
		zo = [];
		zo   = double(scv_data.pres - [pres_spike + pd*(4)*sqrt(h^2/2)]);

		%// Reduce to only data, ignore NaNs
		dat = sa + zo;
		sa = sa(~isnan(dat));
		zo = zo(~isnan(dat));

		%// Get modified gaussian amplitude (+- 20%)
		AA = A*ad;
		
		%// Generate model using update amplitude (AA), center (zo), and height (h)
		gauss = AA*exp((-(zo.^2))/(h.^2));

		%// Check R2
		pl = [pres_spike + pd*(4)*sqrt(h^2/2)] - 2*sqrt((h^2)/2); pl  = round(pl/10)*10;
		ph = [pres_spike + pd*(4)*sqrt(h^2/2)] + 2*sqrt((h^2)/2); ph  = round(ph/10)*10;

		%// Calculate R2 (must be > 0.5);
		zp = [zo + pres_spike + pd*(4)*sqrt(h^2/2)];
		dataX  = scv_data.spice_anom(pl <= scv_data.pres & scv_data.pres <= ph);
		dataY  = scv_data.pres(pl <= scv_data.pres & scv_data.pres <= ph);
		modelX = gauss(pl <= zp & zp <= ph);
		modelY = zp(pl <= zp & zp <= ph);
		if length(dataX) < length(modelX) | length(modelX) < length(dataX)
		    [c,~,~] = intersect(dataY,modelY);
		    ind     = find(min(c) <= dataY & dataY <= max(c));
		    dataX   = dataX(ind);   dataY = dataY(ind);
		    ind     = find(min(c) <= modelY & modelY <= max(c));
		    modelX  = modelX(ind); modelY = modelY(ind); 
		end
		R2(pdcnt,acnt,hcnt) = corr2(dataX,modelX).^2;

		%// Calculate NRMSE (must be < 0.5);
		RMSE                   = sqrt(sum((dataX - modelX).^2)/length(dataX));
		NRMSE(pdcnt,acnt,hcnt) = RMSE/(max(dataX) - min(dataX)); 

		if R2(pdcnt,acnt,hcnt) < R2_thresh
		    lse(pdcnt,acnt,hcnt) = NaN;
		elseif NRMSE(pdcnt,acnt,hcnt) > NRMSE_thresh
		    lse(pdcnt,acnt,hcnt) = NaN;
		else
		    lse(pdcnt,acnt,hcnt) = sum([dataX - modelX].^2);
		end
	    end
	end
    end

    %// Find best zo,A,H combo according to lse
    [minlse,idxlse] = min(lse(:));
    [a,b,c] = ind2sub(size(lse),idxlse);

    %// Update parameters
    A          = A*arng(b);
    H          = hrng(c);
    pres_spike = pres_spike + pdrng(a)*(4)*sqrt(H^2/2);
    zo   = double(scv_data.pres - [pres_spike]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 3         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag bad 'goodness-of-fit' models
    if isnan(lse(a,b,c)==1)
	return
    end

    %// Generate final model
    gauss = A*exp((-(zo.^2))/(H.^2));

    %// Get pressure limits
    pliml  = pres_spike - 2*sqrt((H^2)/2); pliml  = round(pliml/10)*10;
    plimh  = pres_spike + 2*sqrt((H^2)/2); plimh  = round(plimh/10)*10;
    plimcl = pres_spike - sqrt((H^2)/2);   plimcl = round(plimcl/10)*10;
    plimch = pres_spike + sqrt((H^2)/2);   plimch = round(plimch/10)*10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 4         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag too small/too big H (min == 150m max == 1000m)
    if [plimh - pliml] < 150 | [plimh - pliml] > 1200 | pliml < 100
	return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 5         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag profiles without data within limits
    pp = scv_data.pres(~isnan(scv_data.spice_anom));
    if min(pp) > pliml | max(pp) < plimh
	return
    end

    %// Get density limits
    %// Return NaNs if not available (most likely deep)
    indl   = find(scv_data.pres == pliml); 
    if isempty(indl)
	dliml = NaN;
    else
	dliml  = scv_data.sigma0(indl); 
    end
    indh = find(scv_data.pres == plimh);
    if isempty(indh)
	dlimh = NaN;
    else
	dlimh  = scv_data.sigma0(indh);
    end
    indcl = find(scv_data.pres == plimcl);
    if isempty(indcl) == 1
	dlimcl = NaN;
    else
	dlimcl = scv_data.sigma0(indcl);
    end
    indch = find(scv_data.pres == plimch);
    if isempty(indch)
	dlimch = NaN;
    else
	dlimch = scv_data.sigma0(indch);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%        FLAG 6         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %// Flag profiles without N2 signal within core
    indn2 = find(n2_exceed(dlimcl <= dens & dens <= dlimch) > 0);
    n2sum = sum(n2_anom(dlimcl <= dens & dens <= dlimch));
    if isempty(indn2) == 1 | n2sum > 0
	return
    end

    %// Restore pressure of model
    zp = zo + pres_spike;

    %// Grab NRMSE from routine
    NRMSE = NRMSE(a,b,c);

    %// Grab R2 fit from routine
    R2 = R2(a,b,c);

    %// Generate model and results
    gauss_model.x      = gauss;
    gauss_model.y      = zp;
    gauss_model.A      = A;
    gauss_model.zo     = zo;
    gauss_model.h      = H;
    gauss_model.ps     = pres_spike;

    gauss_stats.R2     = R2;
    gauss_stats.NRMSE  = NRMSE;

    gauss_lims.p_lims  = [pliml  plimh];
    gauss_lims.p_clims = [plimcl plimch];
    gauss_lims.p_peak  = pres_spike;
    gauss_lims.d_lims  = [dliml  dlimh];
    gauss_lims.d_clims = [dlimcl dlimch];
    gauss_lims.d_peak  = sig_spike;

catch
    %// Apply flags over again, just to be safe
    gauss_model.x       = 999;
    gauss_model.y       = 999;
    gauss_stats.R2      = 999;
    gauss_stats.NRMSE   = 999;
    gauss_lims.p_lims   = 999;
    gauss_lims.p_clims  = 999;
    gauss_lims.p_peak   = 999;
    gauss_lims.d_lims   = 999;
    gauss_lims.d_clims  = 999;
    gauss_lims.d_peak   = 999;
end

%// Finally, fix orientation of model if minty
if choice == 2
    gauss_model.x = -gauss_model.x;
end
