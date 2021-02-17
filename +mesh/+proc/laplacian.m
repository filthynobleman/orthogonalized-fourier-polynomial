function varargout = laplacian(M, order)
%[S, A] = LAPLACIAN(M)
%Computes the stifness and mass matrix for the Laplace-Beltrami operator,
%using Finite Element discretization.
%
%[S, A, Al] = LAPLACIAN(M)
%Computes the stifness and mass matrix for the Laplace-Beltrami operator,
%using Finite Element discretization. Furthermore, the function returns the
%lumped mass matrix, which is the diagonal of the mass.
%
%[__] = LAPLACIAN(M, order)
%The Finite Element dicretization is computed using hat functions of the
%given order.

    % Infer missing arguments
    if nargin < 2
        order = 1;
    end
    
    % Parameters are ok?
    if order ~= 1 && order ~= 2 && order ~= 3
        error("The order must be integer in [1, 3]. %d given.", order);
    end
    
    if order == 1
        [W, Sc, Sl] = calc_LB_FEM(M);
    elseif order == 2
        [W, Sc, Sl] = mesh.proc.calc_LB_FEM_quad(M);
    else
        [W, Sc, Sl] = calc_LB_FEM_cubic(M);
    end
    
    varargout{1} = W;
    varargout{2} = Sc;
    if nargout == 2
        return;
    end
    varargout{3} = Sl;
end

function [W, Sc, Sl] = calc_LB_FEM(M)
    % Stiffness (p.s.d.)
    angles = zeros(M.m,3);
    for i=1:3
        a = mod(i-1,3)+1;
        b = mod(i,3)+1;
        c = mod(i+1,3)+1;
        ab = M.VERT(M.TRIV(:,b),:) - M.VERT(M.TRIV(:,a),:);
        ac = M.VERT(M.TRIV(:,c),:) - M.VERT(M.TRIV(:,a),:);
        %normalize edges
        ab = ab ./ (sqrt(sum(ab.^2,2))*[1 1 1]);
        ac = ac ./ (sqrt(sum(ac.^2,2))*[1 1 1]);
        % normalize the vectors
        % compute cotan of angles
        angles(:,a) = cot(acos(sum(ab.*ac,2)));
        %cotan can also be computed by x/sqrt(1-x^2)
    end

    indicesI = [M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,3);M.TRIV(:,2);M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,1);M.TRIV(:,3)];
    values   = [angles(:,3);angles(:,1);angles(:,2);angles(:,1);angles(:,3);angles(:,2)]*0.5;
    W = sparse(indicesI, indicesJ, -values, M.n, M.n);
    W = W-sparse(1:M.n,1:M.n,sum(W));

    % Mass
    areas = mesh.proc.tri_areas(M);
    indicesI = [M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,3);M.TRIV(:,2);M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,1);M.TRIV(:,3)];
    values   = [areas(:); areas(:); areas(:); areas(:); areas(:); areas(:)]./12;
    Sc = sparse(indicesI, indicesJ, values, M.n, M.n);
    Sc = Sc+sparse(1:M.n, 1:M.n, sum(Sc));
    
    % Lumped mass
    Sl = spdiags(sum(Sc,2), 0, M.n, M.n);
    
    boundary = calc_boundary_edges(M.TRIV);
    boundary = unique(boundary(:));
%     E = [M.TRIV(:, [1 2]);
%          M.TRIV(:, [2 3]);
%          M.TRIV(:, [3 1])];
%     M.Adj = sparse(E(:, 1), E(:, 2), 1, M.n, M.n);
%     M.Adj = double(logical((M.Adj + M.Adj') > 0));
%     bbound = zeros(M.n, 1);
%     bbound(boundary) = 1;
%     bbound = logical(M.Adj * bbound);
%     bbound = find(bbound);
%     [W, Sc, Sl] = dirichlet_bc(W, Sc, Sl, [boundary; bbound], M.n);
    [W, Sc, Sl] = dirichlet_bc(W, Sc, Sl, boundary, M.n);
    
    W = (W + W')/2;
    Sc = (Sc + Sc')/2;
end

function [W, A, AL] = dirichlet_bc(W_full, A_full, AL_full, boundary, n)
    boundary = unique(boundary(:));
    l = size(boundary,1);
    inside = setdiff(1:n, boundary);
    W = sparse(n, n);
    W(boundary, boundary) = W_full(boundary, boundary);% eye(l);
    W(inside, inside) = W_full(inside, inside);
    A = sparse(n, n);
    A(inside,inside) = A_full(inside, inside);
    AL = sparse(n, n);
    AL(inside,inside) = AL_full(inside, inside);
end

function bd = calc_boundary_edges(triangles)

if isstruct(triangles)
    triangles = triangles.TRIV;
end

[c,d,~,~] = get_boundary(triangles);
bd = zeros(length(d),2);

for i=1:length(c)
    t = triangles(c(i),:);
    v = true(1,3);
    v(d(i)) = false;
    v = t(v);
    bd(i,1) = v(1);
    bd(i,2) = v(2);
end

function [ c,d,I,v ] = get_boundary( tri )
%GET_BOUNDARY determines the boundary edges of a triangular mesh 
%   [c,d] = get_boundary(tri) takes as input a list tri of consistently oriented triangles
%   returns the indices c of the triangles the boundary edges belong to and
%   the (local) indices d (in {1,2,3}) of the vertices opposing the boundary
%   edge
%   One gets the global indices of those vertices via F(sub2ind(size(F),c,d))
%   Via 
%   d1 = mod(d+1,3); d1(d1==0) = 3;
%   d2 = mod(d+2,3); d2(d2==0) = 3;
%   one gets the local indices of the boundary vertices.


if size(tri,1)<size(tri,2)
    tri=tri';
end
m = size(tri,1);

%% Check for opposing halfedges

% Matrix of directed edges
I = [tri(:,1) tri(:,2);
     tri(:,2) tri(:,3);
     tri(:,3) tri(:,1)];

b = not(ismember(I(:,[2 1]),I,'rows'));
b = find(b);



% Triangle indices
c = mod(b,m);
c(c==0) = m;

% vertex opposing boundary edge
d = floor((b-1)/m);
d(d==0)=3;


% Directed boundary edges
I=I(b,:);

% Boundary vertices
[~,~,v] = find(I);
v = unique(v);

end

end

function [Stiff, Mass, LumpedMass] = calc_LB_FEM_cubic(M)

Ia = 1/13440 .* ...
    [76 11 11 18 0 27 27 0 18 36;
    11 76 11 0 18 18 0 27 27 36;
    11 11 76 27 27 0 18 18 0 36;
    18 0 27 540 -189 -135 -54 -135 270 162;
    0 18 27 -189 540 270 -135 -54 -135 162;
    27 18 0 -135 270 540 -189 -135 -54 162;
    27 0 18 -54 -135 -189 540 270 -135 162;
    0 27 18 -135 -54 -135 270 540 -189 162;
    18 27 0 270 -135 -54 -135 -189 540 162;
    36 36 36 162 162 162 162 162 162 1944];

Ib = 1/80 .* ...
    [34 -7 0 -54 27 -3 -3 3 3 0;
    -7 34 0 27 -54 3 3 -3 -3 0;
    0 0 0 0 0 0 0 0 0 0;
    -54 27 0 135 -108 0 0 0 0 0;
    27 -54 0 -108 135 0 0 0 0 0;
    -3 3 0 0 0 135 -27 27 27 -162;
    -3 3 0 0 0 -27 135 -135 27 0;
    3 -3 0 0 0 27 -135 135 -27 0;
    3 -3 0 0 0 27 27 -27 135 -162;
    0 0 0 0 0 -162 0 0 -162 324];

Ic = 1/80 .* ...
    [34 0 -7 3 3 -3 -3 27 -54 0;
    0  0 0 0 0 0 0 0 0 0;
    -7 0 34 -3 -3 3  3 -54 27 0;
    3 0 -3 135 -27 27 27 0 0 -162;
    3 0 -3 -27 135 -135 27 0 0 0;
    -3 0 3 27 -135 135 -27 0 0 0;
    -3 0 3 27 27 -27 135 0 0 -162;
    27 0 -54 0 0 0 0 135 -108 0;
    -54 0 27 0 0 0 0 -108 135 0;
    0  0 0 -162 0 0 -162 0 0 324];

Id = 1/80 .* ...
    [68 -7 -7 -51 30 -6 -6 30 -51 0;
    -7 0 7 24 -57 57 -24 0 0 0;
    -7 7 0 0 0 -24 57 -57 24 0;
    -51 24 0 135 -108 27 27 -27 135 -162;
    30 -57 0 -108 135 -135 27 -27 -27 162;
    -6 57 -24 27 -135 135 54 27 27 -162;
    -6 -24 57 27 27 54 135 -135 27 -162;
    30 0 -57 -27 -27 27 -135 135 -108 162;
    -51 0 24 135 -27 27 27 -108 135 -162;
    0 0 0 -162 162 -162 -162 162 -162 324];

q = size(Ib, 1);

Ad = calc_indices_adj_cubic(M);
Nedges = nnz(triu(Ad));
Ntot = M.n + 2*Nedges + M.m;

TRIE = zeros(M.m, 3);
TRIE(:, 1) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n + ...
    (M.TRIV(:, 1) > M.TRIV(:, 2));
TRIE(:, 2) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n + ...
    (~(M.TRIV(:, 1) > M.TRIV(:, 2)));

TRIE(:, 3) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n + ...
    (M.TRIV(:, 2) > M.TRIV(:, 3));
TRIE(:, 4) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n + ...
    (~(M.TRIV(:, 2) > M.TRIV(:, 3)));

TRIE(:, 5) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n + ...
    (M.TRIV(:, 3) > M.TRIV(:, 1));
TRIE(:, 6) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n + ...
    (~(M.TRIV(:, 3) > M.TRIV(:, 1)));

TRIE(:, 7) = (1:M.m)' + M.n + 2*Nedges;
TRItot = [M.TRIV, TRIE];

P1 = M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 1), :);
P2 = M.VERT(M.TRIV(:, 3), :) - M.VERT(M.TRIV(:, 1), :);

%     Row-wise dot product, repmat to match I#
p11 = dot(P1, P1, 2);
p11b = repmat(p11', q, 1);
p11b = p11b(:);

p22 = dot(P2, P2, 2);
p22b = repmat(p22', q, 1);
p22b = p22b(:);

p12 = dot(P1, P2, 2);
p12b = repmat(p12', q, 1);
p12b = p12b(:);

pre2 = vecnorm(cross(P1,P2, 2), 2, 2);
pre2b = repmat(pre2', q, 1);
pre2b = pre2b(:);

I1b = repmat(Ib, M.m, 1);
I2b = repmat(Ic, M.m, 1);
I3b = repmat(Id, M.m, 1);
I4b = repmat(Ia, M.m, 1);

Alocal2 = (p11b .* I2b + p22b .* I1b - p12b .* I3b) ./ pre2b;
Blocal2 = I4b .* pre2b;

va = Alocal2';
va = va(:);

vb = Blocal2';
vb = vb(:);

idx_rows = repmat(1:M.m, q*q, 1);
idx_rows = idx_rows(:);

idx_cols = repmat(1:q, q, 1);
idx_cols = idx_cols(:);
idx_cols = repmat(idx_cols, M.m, 1);
rows = TRItot(sub2ind(size(TRItot), idx_rows, idx_cols));

idx_cols = repmat(1:q, 1, M.m * q)';
cols = TRItot(sub2ind(size(TRItot), idx_rows, idx_cols));

A = sparse(rows, cols, va, Ntot, Ntot);
B = sparse(rows, cols, vb, Ntot, Ntot);

Stiff = A;
Mass = B;
% LumpedMass = spdiag(sum(B, 2));
LumpedMass = spdiags(sum(B, 2), 0, Ntot, Ntot);
end

function Ad = calc_indices_adj_cubic(M)
    indicesI = [M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,3); M.TRIV(:,2); M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,1); M.TRIV(:,3)];
    Ad = sparse(indicesI, indicesJ, ones(M.m,6), M.n, M.n);
    [indicesI, indicesJ, ~] = find(triu(Ad));
    Nedges = nnz(triu(Ad));
    Ad = sparse(indicesI, indicesJ, 1:2:2*Nedges, M.n, M.n);
    Ad = Ad + Ad';
end

