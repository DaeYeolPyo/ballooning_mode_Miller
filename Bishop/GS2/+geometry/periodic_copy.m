function a = periodic_copy(obj, a0, ext, nperiod)
    if nperiod <= 1
        return;
    end

    ntheta = 2*obj.nth;
    itot_min = -nth + (-nperiod + 1) * ntheta;
    itot_max =  nth + ( nperiod - 1) * ntheta;

    itot_grid = itot_min:itot_max;
    a = nan(1, numel(itot_grid));

    offset = 1 - itot_min;

    % central period, k = 0
    for i = -nth:nth
        itot = i;
        a(itot + offset) = a0(i + nth + 1);
    end

    % negative periods
    for k = -nperiod+1:-1
        for i = -nth:nth
            itot = i + k * ntheta;
            a(itot + offset) = a0(i + nth + 1) + k * ext;
        end
    end

    % positive periods
    for k = 1:nperiod-1
        for i = -nth:nth
            itot = i + k * ntheta;
            a(itot + offset) = a0(i + nth + 1) + k * ext;
        end
    end
end