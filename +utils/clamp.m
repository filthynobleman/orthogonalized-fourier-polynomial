function F = clamp(F, minv, maxv)
%F = CLAMP(F, minv, maxv) Clamps the values of the given matrix F between
%the interval [minv, maxv]. Notice F is allowed to contain infinite values
%and NaN values. The latters will be ignored in clamping operation. Also,
%notice minv and maxv are allowed to be infinite, but they cannot be NaN
%and it must hold minv < maxv

    % Check validity or arguments
    if isnan(minv)
        error("minv must be definite.");
    elseif isnan(maxv)
        error("maxv must be definite.");
    elseif minv >= maxv
        error("minv must be strictly less than maxv.");
    end
    
    % Clamp the matrix
    F = min(max(F, minv), maxv);
end

