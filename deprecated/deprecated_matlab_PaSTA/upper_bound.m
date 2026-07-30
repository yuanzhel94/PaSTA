function idx = upper_bound(x, val)
    lo = 1;
    hi = numel(x) + 1;

    while lo < hi
        mid = floor((lo + hi) / 2);
        if x(mid) <= val
            lo = mid + 1;
        else
            hi = mid;
        end
    end

    idx = lo - 1;
end