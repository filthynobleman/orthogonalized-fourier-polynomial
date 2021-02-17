function [I1, I2, I3, I4] = calc_FEM_test_funcs(Order, ICache)

    if nargin < 2
        ICache = "Cache/FEM/";
    end
    [I1, I2, I3, I4, found] = check_cache(Order, ICache);
    if found
        return;
    end

    syms u v;
    U(1) = u^0;
    V(1) = v^0;
    for i = 1:Order
        U(i + 1) = u * U(i);
        V(i + 1) = v * V(i);
    end
    
    UV = V' * U;
    UV = fliplr(UV);
    Vars = [];
    for i = 0:Order
        Vars = [Vars; diag(UV, Order - i)];
    end
    Vars = Vars';
    
    M = [];
    M = [subs(subs(Vars, u, 0), v, 0);
         subs(subs(Vars, u, 1), v, 0);
         subs(subs(Vars, u, 0), v, 1)];
    Space = linspace(0, 1, Order + 1);
    for i = 2:Order
        M = [M; subs(subs(Vars, u, Space(i)), v, 0)];
    end
    for i = 2:Order
        M = [M; subs(subs(Vars, u, 1 - Space(i)), v, Space(i))];
    end
    for i = 2:Order
        M = [M; subs(subs(Vars, u, 0), v, 1 - Space(i))];
    end
    for i = 2:Order
        for j = 2:Order - i + 1
            M = [M; subs(subs(Vars, u, Space(i)), v, Space(j))];
        end
    end
    
    % Compute coefficients for each test function
    C = zeros(length(Vars));
    for i = 1:length(Vars)
        rhs = zeros(length(Vars), 1);
        rhs(i) = 1;
        C(:, i) = M \ rhs;
    end
    H = C' * Vars';
    
    
    I1 = zeros(length(Vars));
    I2 = zeros(length(Vars));
    I3 = zeros(length(Vars));
    I4 = zeros(length(Vars));
    for i = 1:length(Vars)
        for j = i:length(Vars)
            f = H(i) * H(j);
            I4(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), u) * diff(H(j), u);
            I1(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), v) * diff(H(j), v);
            I2(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), v) * diff(H(j), u) + diff(H(i), u) * diff(H(j), v);
            I3(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
        end 
    end
    
    
    I1 = I1 + I1';
    I1(1:length(Vars)+1:end) = diag(I1) ./ 2;
    
    I2 = I2 + I2';
    I2(1:length(Vars)+1:end) = diag(I2) ./ 2;
    
    I3 = I3 + I3';
    I3(1:length(Vars)+1:end) = diag(I3) ./ 2;
    
    I4 = I4 + I4';
    I4(1:length(Vars)+1:end) = diag(I4) ./ 2;
    
    cachename = sprintf("Cache/FEM/Order%d.mat", Order);
    save(cachename, 'I1', 'I2', 'I3', 'I4');
    
end


function [I1, I2, I3, I4, found] = check_cache(p, cachedir)
    if nargin == 1
        cachedir = "Cache/FEM/";
    end
    cachename = sprintf("%s/Order%d.mat", cachedir, p);
    if exist(cachename, 'file')
        Cache = load(cachename);
        I1 = Cache.I1;
        I2 = Cache.I2;
        I3 = Cache.I3;
        I4 = Cache.I4;
        found = true;
    else
        I1 = [];
        I2 = [];
        I3 = [];
        I4 = [];
        found = false;
    end
end

