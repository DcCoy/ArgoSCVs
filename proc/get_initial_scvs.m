function [scv_index] = get_initial_scvs(sp_test,n_test,pyc_dens,sigma0,N2_range);
%// Script to find SCVs by looking for spice anomalies that exceed their
%// respective thresholds at densities where n2 anomalies also exceed
%// their thresholds (within N2_range)

%// Define index for SCVs, and progress-checker
idx      = zeros(1,length(sp_test));
progress = [0:100000:length(sp_test)];

%// Loop through and check for t/s/n2 thresholds being passed
for i = 1:length(sp_test)

    %// Keep track of progress
    if ismember(i,progress)==1
	disp([num2str(floor((i/length(sp_test))*100)),'% complete'])
    end

    %// Get data below pycnocline
    ind = find(sigma0(:,i) >= pyc_dens(i));
    if isempty(ind) == 1 | length(ind) < 2
	continue
    end

    %// Find where sp_test is positive
    indsp  = find(sp_test(ind,i) > 0);

    %// Find where n_test is positive
    indn  = find(n_test(ind,i) > 0);

    %// Continue if empty or only one point exceeds thresh
    if isempty(indsp)==1 | isempty(indn)==1
	continue
    elseif length(indsp) < 2
	continue
    end

    %// Get densities where thresholds are passed
    sp_dens = sigma0(ind(indsp),i);
    n2_dens = sigma0(ind(indn),i);

    %// Check that n2_dens is within 0.1kg/m3 of ts_dens
    for ii = 1:length(n2_dens)
	scv_test = abs([sp_dens - n2_dens(ii)]);
	ind_test = find(scv_test <= N2_range);
	if isempty(ind_test) == 0
	    idx(i) = 1;
	end
    end
end

%// Record scv_index
scv_index = find(idx == 1);

    
