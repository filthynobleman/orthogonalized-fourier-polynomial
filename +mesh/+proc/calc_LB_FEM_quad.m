function [Stiff, Mass, LumpedMass] = calc_LB_FEM_quad(M)

% temp = load('C:\Users\fmagg\Google Drive\Master\Thesis\Code\quadratic-fem\I1.MAT');
% I1 = temp.Expression1;
I1 = [ 1/2,  1/6, 0, -2/3,    0,    0;
       1/6,  1/2, 0, -2/3,    0,    0;
         0,    0, 0,    0,    0,    0;
      -2/3, -2/3, 0,  4/3,    0,    0;
         0,    0, 0,    0,  4/3, -4/3;
         0,    0, 0,    0, -4/3,  4/3];
% temp = load('C:\Users\fmagg\Google Drive\Master\Thesis\Code\quadratic-fem\I2.MAT');
% I2 = temp.Expression1;
I2 = [ 1/2, 0,  1/6,    0,    0, -2/3;
         0, 0,    0,    0,    0,    0;
       1/6, 0,  1/2,    0,    0, -2/3;
         0, 0,    0,  4/3, -4/3,    0;
         0, 0,    0, -4/3,  4/3,    0;
      -2/3, 0, -2/3,    0,    0,  4/3];
% temp = load('C:\Users\fmagg\Google Drive\Master\Thesis\Code\quadratic-fem\I3.MAT');
% I3 = temp.Expression1;
I3 = [   1,  1/6,  1/6, -2/3,    0, -2/3;
       1/6,    0, -1/6, -2/3,  2/3,    0;
       1/6, -1/6,    0,    0,  2/3, -2/3;
      -2/3, -2/3,    0,  4/3, -4/3,  4/3;
         0,  2/3,  2/3, -4/3,  4/3, -4/3;
      -2/3,    0, -2/3,  4/3, -4/3,  4/3];
% temp = load('C:\Users\fmagg\Google Drive\Master\Thesis\Code\quadratic-fem\I4.MAT');
% I4 = temp.Expression1;
I4 = [  1/60, -1/360, -1/360,     0, -1/90,     0;
      -1/360,   1/60, -1/360,     0,     0, -1/90;
      -1/360, -1/360,   1/60, -1/90,     0,     0;
           0,      0,  -1/90,  4/45,  2/45,  2/45;
       -1/90,      0,      0,  2/45,  4/45,  2/45;
           0,  -1/90,      0,  2/45,  2/45,  4/45];

% [I1, I2, I3, I4] = mesh.proc.calc_FEM_test_funcs(2);


Ad = calc_indices_adj(M);
Nedges = nnz(triu(Ad));
Ntot = M.n + Nedges;

TRIE = zeros(M.m, 3);
% for i=1:M.m
%    TRIE(i, 1) = Ad(M.TRIV(i, 1), M.TRIV(i, 2)) + M.n;
%    TRIE(i, 2) = Ad(M.TRIV(i, 2), M.TRIV(i, 3)) + M.n;
%    TRIE(i, 3) = Ad(M.TRIV(i, 3), M.TRIV(i, 1)) + M.n;
% end
TRIE(:, 1) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n;
TRIE(:, 2) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n;
TRIE(:, 3) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n;
TRItot = [M.TRIV, TRIE];

% A = spalloc(Ntot,Ntot,0);
% B = spalloc(Ntot,Ntot,0); 
P1 = M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 1), :);
P2 = M.VERT(M.TRIV(:, 3), :) - M.VERT(M.TRIV(:, 1), :);
P11 = dot(P1, P1, 2);
P22 = dot(P2, P2, 2);
P12 = dot(P1, P2, 2);
pre = cross(P1, P2, 2);
pre = sqrt(sum(pre.^2, 2));
Alocal = reshape(P11, 1, 1, M.m) .* repmat(I2, 1, 1, M.m) + ...
         reshape(P22, 1, 1, M.m) .* repmat(I1, 1, 1, M.m) - ...
         reshape(P12, 1, 1, M.m) .* repmat(I3, 1, 1, M.m);
Alocal = Alocal ./ reshape(pre, 1, 1, M.m);
Blocal = repmat(I4, 1, 1, M.m) .* reshape(pre, 1, 1, M.m);
AI = reshape(TRItot(:, repmat(1:6, 1, 6))', 36*M.m, 1);
AJ = reshape(TRItot(:, reshape(repmat(1:6, 6, 1), 1, 36))', 36*M.m, 1);
BI = reshape(TRItot(:, repmat(1:6, 1, 6))', 36*M.m, 1);
BJ = reshape(TRItot(:, reshape(repmat(1:6, 6, 1), 1, 36))', 36*M.m, 1);
A = sparse(AI, AJ, reshape(Alocal, 1, 36*M.m), Ntot, Ntot);
B = sparse(BI, BJ, reshape(Blocal, 1, 36*M.m), Ntot, Ntot);

% for i=1:M.m
% %     P1 = M.VERT(M.TRIV(i, 2), :) - M.VERT(M.TRIV(i, 1), :);
% %     P2 = M.VERT(M.TRIV(i, 3), :) - M.VERT(M.TRIV(i, 1), :);
% %     pre = norm(cross(P1,P2));
% %     Alocal = (P11(i)*I2 + P22(i)*I1 - P12(i)*I3)./pre(i); 
% %     Blocal = I4*pre(i);
% %     for l=1:6
% %         for m=1:6
% %             A(TRItot(i,l), TRItot(i,m)) = A(TRItot(i,l), TRItot(i,m)) + Alocal(l,m, i);
% %             B(TRItot(i,l), TRItot(i,m)) = B(TRItot(i,l), TRItot(i,m)) + Blocal(l,m, i);
% %         end
% %     end
%     A(TRItot(i, 1:6), TRItot(i, 1:6)) = A(TRItot(i, 1:6), TRItot(i, 1:6)) + Alocal(1:6, 1:6, i);
%     B(TRItot(i, 1:6), TRItot(i, 1:6)) = B(TRItot(i, 1:6), TRItot(i, 1:6)) + Blocal(1:6, 1:6, i);
% end

Stiff = A;
Mass = B;
LumpedMass = sparse(1:length(B), 1:length(B), sum(B, 2));

end

function Ad = calc_indices_adj(M)

indicesI = [M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,3); M.TRIV(:,2); M.TRIV(:,1)];
indicesJ = [M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,1); M.TRIV(:,3)];
Ad = sparse(indicesI, indicesJ, ones(M.m,6), M.n, M.n);
[indicesI, indicesJ, ~] = find(triu(Ad));
Nedges = nnz(triu(Ad));
Ad = sparse(indicesI, indicesJ, 1:Nedges, M.n, M.n);
Ad = Ad + Ad';

end