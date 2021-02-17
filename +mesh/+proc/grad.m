function [G] = grad(M)
  % GRAD Compute the numerical gradient operator for triangle or tet meshes
  %
  % G = grad(V,F)
  %
  % Inputs:
  %   V  #vertices by dim list of mesh vertex positions
  %   F  #faces by simplex-size list of mesh face indices
  % Outputs:
  %   G  #faces*dim by #V Gradient operator
  %
  % Example:
  %   L = cotmatrix(V,F);
  %   G  = grad(V,F);
  %   dblA = doublearea(V,F);
  %   GMG = -G'*repdiag(diag(sparse(dblA)/2),size(V,2))*G;
  %
  %   % Columns of W are scalar fields
  %   G = grad(V,F);
  %   % Compute gradient magnitude for each column in W
  %   GM = squeeze(sqrt(sum(reshape(G*W,size(F,1),size(V,2),size(W,2)).^2,2)));
  %

  V = M.VERT;
  F = M.TRIV;
  dim = size(V,2);
  ss = size(F,2);
  switch ss
  case 2
    % Edge lengths
    len = abs(V(F(:,2))-V(F(:,1)));
    % Gradient is just staggered grid finite difference
    G = sparse(repmat(1:size(F,1),2,1)',F,[1 -1]./len, size(F,1),size(V,1));
  case 3
    % append with 0s for convenience
    if size(V,2) == 2
      V = [V zeros(size(V,1),1)];
    end
  
    % Gradient of a scalar function defined on piecewise linear elements (mesh)
    % is constant on each triangle i,j,k:
    % grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
    % grad(Xijk) = Xj * (Vi - Vk)^R90 / 2A + Xk * (Vj - Vi)^R90 / 2A + 
    %             -Xi * (Vi - Vk)^R90 / 2A - Xi * (Vj - Vi)^R90 / 2A
    % where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
    % i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
    % 90 degrees
    %
    % renaming indices of vertices of triangles for convenience
    i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
    % #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v32 = V(i3,:) - V(i2,:);  v13 = V(i1,:) - V(i3,:); v21 = V(i2,:) - V(i1,:);
  
    % area of parallelogram is twice area of triangle
    % area of parallelogram is || v1 x v2 || 
    n  = cross(v32,v13,2); 
    
    % This does correct l2 norm of rows, so that it contains #F list of twice
    % triangle areas
    dblA = utils.repnorm(n, 'rows');
    
    % now normalize normals to get unit normals
    u = n ./ vecnorm(n, 2, 2);
  
    % rotate each vector 90 degrees around normal
    %eperp21 = bsxfun(@times,normalizerow(cross(u,v21)),normrow(v21)./dblA);
    %eperp13 = bsxfun(@times,normalizerow(cross(u,v13)),normrow(v13)./dblA);
    eperp21 = bsxfun(@times,cross(u,v21),1./dblA);
    eperp13 = bsxfun(@times,cross(u,v13),1./dblA);
  
    %g =                                              ...
    %  (                                              ...
    %    repmat(X(F(:,2)) - X(F(:,1)),1,3).*eperp13 + ...
    %    repmat(X(F(:,3)) - X(F(:,1)),1,3).*eperp21   ...
    %  );
  
    G = sparse( ...
      [0*size(F,1)+repmat(1:size(F,1),1,4) ...
      1*size(F,1)+repmat(1:size(F,1),1,4) ...
      2*size(F,1)+repmat(1:size(F,1),1,4)]', ...
      repmat([F(:,2);F(:,1);F(:,3);F(:,1)],3,1), ...
      [eperp13(:,1);-eperp13(:,1);eperp21(:,1);-eperp21(:,1); ...
       eperp13(:,2);-eperp13(:,2);eperp21(:,2);-eperp21(:,2); ...
       eperp13(:,3);-eperp13(:,3);eperp21(:,3);-eperp21(:,3)], ...
      3*size(F,1), size(V,1));
  
    %% Alternatively
    %%
    %% f(x) is piecewise-linear function:
    %%
    %% f(x) = ? ?i(x) fi, f(x ? T) = ?i(x) fi + ?j(x) fj + ?k(x) fk
    %% ?f(x) = ...                 = ??i(x) fi + ??j(x) fj + ??k(x) fk 
    %%                             = ??i fi + ??j fj + ??k) fk 
    %%
    %% ??i = 1/hjk ((Vj-Vk)/||Vj-Vk||)^perp = 
    %%     = ||Vj-Vk|| /(2 Aijk) * ((Vj-Vk)/||Vj-Vk||)^perp 
    %%     = 1/(2 Aijk) * (Vj-Vk)^perp 
    %% 
    %m = size(F,1);
    %eperp32 = bsxfun(@times,cross(u,v32),1./dblA);
    %G = sparse( ...
    %  [0*m + repmat(1:m,1,3) ...
    %   1*m + repmat(1:m,1,3) ...
    %   2*m + repmat(1:m,1,3)]', ...
    %  repmat([F(:,1);F(:,2);F(:,3)],3,1), ...
    %  [eperp32(:,1);eperp13(:,1);eperp21(:,1); ...
    %   eperp32(:,2);eperp13(:,2);eperp21(:,2); ...
    %   eperp32(:,3);eperp13(:,3);eperp21(:,3)], ...
    %  3*m,size(V,1));
  
    if dim == 2
      G = G(1:(size(F,1)*dim),:);
    end
  
  
    % Should be the same as:
    % g = ... 
    %   bsxfun(@times,X(F(:,1)),cross(u,v32)) + ...
    %   bsxfun(@times,X(F(:,2)),cross(u,v13)) + ...
    %   bsxfun(@times,X(F(:,3)),cross(u,v21));
    % g = bsxfun(@rdivide,g,dblA);

  case 4
    % really dealing with tets
    T = F;
    % number of dimensions
    assert(dim == 3);
    % number of vertices
    n = size(V,1);
    % number of elements
    m = size(T,1);
    % simplex size
    assert(size(T,2) == 4);
  
    % f(x) is piecewise-linear function:
    %
    % f(x) = ? ?i(x) fi, f(x ? T) = ?i(x) fi + ?j(x) fj + ?k(x) fk + ?l(x) fl
    % ?f(x) = ...                 = ??i(x) fi + ??j(x) fj + ??k(x) fk + ??l(x) fl
    %                             = ??i fi + ??j fj + ??k fk + ??l fl
    %
    % ??i = 1/hjk = Ajkl / 3V * (Facejkl)^perp
    %     = Ajkl / 3V * (Vj-Vk)x(Vl-Vk)
    %     = Ajkl / 3V * Njkl / ||Njkl||
    % 
  
    % get all faces
    F = [ ...
      T(:,1) T(:,2) T(:,3); ...
      T(:,1) T(:,3) T(:,4); ...
      T(:,1) T(:,4) T(:,2); ...
      T(:,2) T(:,4) T(:,3)];
    % compute areas of each face
    A = utils.doublearea(V,F)/2;
    N = normr(utils.normals(V,F));
  
    % compute volume of each tet
    vol = volume(V,T);
  
    G = sparse( ...
      [0*m + repmat(1:m,1,4) ...
       1*m + repmat(1:m,1,4) ...
       2*m + repmat(1:m,1,4)], ...
      repmat([T(:,4);T(:,2);T(:,3);T(:,1)],3,1), ...
      repmat(A./(3*repmat(vol,4,1)),3,1).*N(:), ...
      3*m,n);
  
  end
end