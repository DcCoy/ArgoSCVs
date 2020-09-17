% Find peaks function without the need of a matlab toolbox
function [pks, locs] = findpeaks2(data)
    pks  = zeros(numel(data), 1);
    locs = zeros(numel(data), 1);
    count = 0;
    for i = 2:(numel(data) - 1)
        if (data(i - 1) < data(i)) && (data(i) > data(i + 1))
            count = count + 1;
            pks(count)  = data(i);
            locs(count) = i;
        end
    end
    pks  = pks(1:count);
    locs = locs(1:count);
end
