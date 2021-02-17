function F = rescale(F, minv, maxv)
%F = RESCALE(F, min, max) Rescales values of matrix F so to fit the
%interval [minv, maxv]. If F contains some infinite or NaN values, they are
%ignored. Also, minv and maxv must be finite, not NaN and such that minv <
%maxv.
    
    % Check for validity of arguments
    if isnan(minv) || isinf(minv)
        error("minv must be finite and not NaN.");
    elseif isnan(maxv) || isinf(maxv)
        error("maxv must be finite and not NaN.");
    elseif minv >= maxv
        error("minv must be stricly less than maxv.");
    end
    
    % Compute the target interval
    intv = maxv - minv;
    % If F has no inf or nans
    if all(~isinf(F), 'all') && all(~isnan(F), 'all')
        % Normalize F
        minF = min(F, [], 'all');
        maxF = max(F, [], 'all');
        F = (F - minF) / (maxF - minF);
        F = F * intv + minv;
    % Otherwise, normalize all the valid elements
    else
        valids = find(~isinf(F) && ~isnan(F));
        minF = min(F(valids), [], 'all');
        maxF = max(F(valids), [], 'all');
        F(valids) = (F(valids) - minF) / (maxF - minF);
        F(valids) = F(valids) * intv + minv;
    end
end

