function M = init(varargin)
%M = INIT(filename) Loads a triangle mesh from the given file.
%
%M = INIT(--, opts) Loads a triangle mesh from the given file and executes
%the given options.
%
%M = INIT('Cache', cachename) Loads a triangle mesh from the given cache
%file. This function assumes that the cache file contains only the triangle
%mesh or, at least, that the triangle mesh name is the first one returned
%by the function fieldnames() on the loaded data.
%
%
% Acceptable options:
%   - 'Laplace': Computes the Laplace-Beltrami operator, in form of a mass
%   matrix and a stiffness matrix.
%   - 'LaplaceBasis', k: Computes the Lalace-Beltrami operator, as in the
%   'Laplace' option, and its first k eigenfunctions. Notice the
%   initialization actually computes the first k + 1 eigenfunctions, since
%   the first one is always constant.
%   WARNING: In this last setting, k is expected to be an integer value. To
%   ensure this, the function will cast it to int64, so any floating point
%   precision is lost and the number is rounded.

    if nargin == 1
        % If only one argument is given, load the mesh and return
        if endsWith(varargin{1}, '.off')
            M = mesh.io.load_off(varargin{1});
        elseif endsWith(varargin{1}, '.ply')
            M = mesh.io.load_ply(varargin{1});
        elseif endsWith(varargin{1}, '.obj')
            M = mesh.io.load_obj(varargin{1});
        else
            error("Format not supported yet.");
        end
        
        % Be sure the vertex numbering in the triangles starts from 0
        if min(min(M.TRIV)) == 0
            M.TRIV = M.TRIV + 1;
        end
        
    else
        % The first argument MUST be a string or a character array
        if ~isstring(varargin{1}) && ~ischar(varargin{1})
            error(strcat("The first argument for mesh.init() ", ...
                         "must be a string or a character array."));
        end
        
        % Parse the arguments
        if lower(varargin{1}) == "cache"
            args = parse_args(varargin{:});
        else
            args = parse_args(varargin{2:end});
        end
        
        % If the first argument is the string "Cache"
        if args.usecache
            % Check validity of arguments
            if ~isstring(args.cache) && ~ischar(args.cache)
                error(strcat("If mesh.init() gets 'Cache' as argument", ...
                             ", the second argument must be a string", ...
                             " or a character array."));
            end
            % Load the data and return
            data = load(args.cache);
            fields = fieldnames(data);
            M = data.(fields{1});
            return;
        end
        % Otherwise, the first argument is the filename
        M = mesh.io.load_off(varargin{1});
        % Now, check for the other options
        % Parse the arguments
        args = parse_args(varargin{2:end});
        if args.laplace
            % If Laplace-Beltrami operator is requested, compute it with
            % the desired order
            [M.S, M.A, M.Al] = mesh.proc.FEM_order_p(M, args.laplaceorder);
        end
        if args.laplacebasis > 0
            % If requested, compute the laplace basis
            try
                k = int64(args.laplacebasis);
            catch
                error(strcat("The length of the Laplace-Beltrami ", ...
                             "basis must be castable to int64."));
            end
            M = mesh.proc.compute_laplace_basis(M, double(k));
        end
    end
    
    % Alias the vertices
    M.X = M.VERT(:, 1);
    M.Y = M.VERT(:, 2);
    M.Z = M.VERT(:, 3);
    M.TRIV = M.TRIV;
    M.n = M.n;
    M.m = M.m;
end

function args = parse_args(varargin)
    i = 1;
    args.usecache = false;
    args.cache = "";
    args.laplaceorder = 0;
    args.laplacebasis = 0;
    args.laplace = false;
    while i <= nargin
        if lower(varargin{i}) == "cache"
            args.usecache = true;
            args.cache = varargin{i + 1};
            break;
        elseif lower(varargin{i}) == "laplacebasis"
            args.laplace = true;
            args.laplacebasis = varargin{i + 1};
            args.laplaceorder = 1;
            i = i + 1;
        elseif lower(varargin{i}) == "laplace"
            args.laplace = true;
            args.laplaceorder = 1;
        elseif lower(varargin{i}) == "laplaceorder"
            args.laplace = true;
            args.laplaceorder = varargin{i + 1};
            i = i + 1;
        else
            error("'%s' is not a valid option.", varargin{i});
        end
        
        i = i + 1;
    end
end

