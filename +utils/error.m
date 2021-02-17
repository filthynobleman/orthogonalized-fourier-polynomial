function err = error(target, f, varargin)
%err = ERROR(target, f, tol) Computes the mean squared relative error
%between target and f. The denominator is the maximum between target and
%tol.
%
%err = ERROR(target, f) Same as full arguments, but tolerance is computed
%by the function as 10^(-5) times the maximum between the mean of the
%abolute values of target and 1.
%
%err = ERROR(__, thresh, perc) After the error computation, the function
%excludes from f all the points where the local error is greater than
%thresh*err, if the number of such points is less than or equals to
%perc*size(f, 1). If the exclusion takes place, the error is recomputed. No
%further iterations are performed.
%WARNING: drastical reduction in performance if more than one function are
%given.

    % Retrieve the arguments
    if nargin == 5
        thresh = varargin{2};
        perc = varargin{3};
    elseif nargin == 4
        thresh = varargin{1};
        perc = varargin{2};
    end
    if nargin == 2 || nargin == 4
        tol = 1e-5 * max(mean(abs(target), 'all'), 1);
    else
        tol = varargin{1};
    end
    % Compute the error
    num = (target - f);
    den = max(tol, abs(target));
    n = numel(target);
    err = norm(num ./ den) / sqrt(n);
    % If the thresh and perc are not given, compute the error and return
    if nargin < 4
        return;
    end
    % Initialize the error and the number of considered elements
    errtot = 0;
    elmtnum = 0;
    % For each function...
    for i = 1:size(f, 2)
        % Compute the points where the error is bounded by thresh*err
        num = target(:, i) - f(:, i);
        den = max(tol, abs(target(:, i)));
        incl = find(abs(num ./ den) <= thresh * err);
        % If the number of points exceeding such a threshold is too large,
        % does not exclude any point. Simply sum the local error to the
        % total error and consider each point in the error computation
        if size(f, 1) - incl > perc * size(f, 1)
            errtot = errtot + sum((num ./ den).^2);
            elmtnum = elmtnum + size(f, 1);
        % Otherwise, exclude the points with a too large error and
        % recompute the local error. Add the local error to the total one
        % and consider in the computation only the non-exluded points
        else
            floc = f(incl, i);
            num = target(incl, i) - floc;
            den = max(tol, abs(target(incl, i)));
            errtot = errtot + sum((num ./ den).^2);
            elmtnum = elmtnum + length(incl);
        end
        % The real error is the total error divided by the square root of
        % the number of considered elements.
        err = sqrt(errtot / elmtnum);
    end
    
end

