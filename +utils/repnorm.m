function res = repnorm(varargin)
%res = REPNORM(X) Returns the norm of each row of the matrix X.
%
%res = REPNORM(X, p) Returns the norm p of each row of the matrix X.
%
%res = REPNORM(--, dim) Select the dimension of the matrix on which
%computing the norm.

    % Check we have at least one argument
    if nargin == 0
        error("repnorm requires at least one argument.");
    end
    % Check we have at most 3 arguments
    if nargin > 3
        error("repnorm accepts at most three arguments.");
    end
    X = varargin{1};
    % Get the size of X
    [n, m] = size(X);
    % Initialize p and dim
    p = 2;
    dim = "rows";
    % If the number of arguments is two, check if it is a string or a
    % scalar
    if nargin >= 2
        if isscalar(varargin{2})
            p = varargin{2};
            if nargin == 3
                if ischar(varargin{2})
                    dim = convertCharsToStrings(varargin{2});
                elseif isstring(varargin{2})
                    dim = varargin{2};
                else
                    error("The third argument must be a string or a character array.");
                end
            end
        elseif ischar(varargin{2})
            dim = convertCharsToStrings(varargin{2});
            if nargin == 3
                error("If the second argument is a character array, repnorm does not accept a third argument.");
            end
        elseif isstring(varargin{2})
            dim = varargin{2};
            if nargin == 3
                error("If the second argument is a string, repnorm does not accept a third argument.");
            end
        else
            error("The second argument must be a scalar, a string or a character array.");
        end
    end
    
    % Determine the axis
    if dim == "rows"
        axis = 2;
    elseif dim == "cols"
        axis = 1;
    else
        error("The dimension must be 'rows' or 'cols'.");
    end
    
    % Compute the norm
    res = sum(X.^p, axis).^(1/p);
end