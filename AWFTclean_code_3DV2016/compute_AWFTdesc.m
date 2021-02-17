function desc = compute_AWFTdesc(shape,functions,param)
%
% desc = compute_AWFTdesc(shape,functions_names,parameters)
%    computes the anisotropic windowed Fourier tranform (AWFT) of a shape
% input:
%       shape = struct with shape.X, shape.Y, shape.Z, shape.TRIV: 
%       functions = cell where functions{i} = string with the name
%       params,
%
% output:
%    desc: AWFT descriptor on the input shape
%

freqs0 = param.freqs0;
taus = param.taus; 
angles = param.angles;
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
alpha =  param.alpha;
k = param.n_eigen;

% compute isotropic LBO
options.alpha = 0;
options.angle = 0;
[shape.Phi, shape.Lambda, shape.W, shape.A] = calc_ALB([shape.X,shape.Y,shape.Z], shape.TRIV, options);
totA= 1./sum(full(diag(shape.A)));
L1.A = totA.*shape.A;
L1.evecs = sqrt(totA).*shape.Phi;
L1.evals = totA.*shape.Lambda;
L1.n = length(shape.X);

% compute curvature
curv_options.curvature_smoothing = options.curv_smooth;
curv_options.verb = 0;
[ ~,  ~, k2, k1,Cmean,Cgauss,Normal] = compute_curvature([shape.X,shape.Y,shape.Z], shape.TRIV, curv_options);

% empty initialization
desc = [];

% prepare function
fs = zeros(length(shape.X),length(functions));
for s = 1 : length(functions)
    if strcmp(functions{s},'HKS')
        [HKS] = calc_HKS(L1);
        fs(:,s)= HKS(:,80);
    elseif strcmp(functions{s},'WKS')
       [WKS] = calc_WKS(L1, 100, 6.0);
        fs(:,s)= WKS(:,80);
    elseif strcmp(functions{s},'Constant')
        fs(:,s) = ones(L1.n,1);
    elseif strcmp(functions{s},'Fiedler')
        fs(:,s) = L1.evecs(:,2);
    elseif strcmp(functions{s},'Cmean')
        fs(:,s) = Cmean;
    elseif strcmp(functions{s},'Cgauss')
        fs(:,s) = Cgauss;
    elseif strcmp(functions{s},'normal1')
        fs(:,s) = Normal(1,:)';
    elseif strcmp(functions{s},'normal2')
        fs(:,s) = Normal(2,:)';
    elseif strcmp(functions{s},'normal3')
        fs(:,s) = Normal(3,:)';
    elseif strcmp(functions{s},'x')
        fs(:,s) = shape.X;
    elseif strcmp(functions{s},'y')
        fs(:,s) = shape.Y;
    elseif strcmp(functions{s},'z')
        fs(:,s) = shape.Z;
    elseif strcmp(functions{s},'ShapeIndex')
        fs(:,s)  = -(2/pi).*atan((k1+k2)./(k1-k2));
    elseif strcmp(functions{s}(1:6),'geovec')
        [gv] = compute_geovec(L1);
        fs(:,s) = gv(:,str2num(functions{s}(7:end)));
    else 
        disp('[!] not available descriptor')
    end
end

% loop over the ghats
for ell = 1:length(taus) 
    for an=1:length(angles)
        for al = 1:length(alpha)
            options.alpha = alpha(al);
            if angles(an) == -1
                % isotropic case
                options.alpha = 0;
                options.angle = 0;
                [evecs, evals, ~, ~] = calc_ALB([shape.X,shape.Y,shape.Z], shape.TRIV, options);
            else
                % anisotropic case
                options.angle = pi*angles(an)/180;
                [evecs, evals, ~, ~] = calc_ALB([shape.X,shape.Y,shape.Z], shape.TRIV, options);
            end
            
            % current ghat
            tau = sum(full(diag(shape.A)))*taus(ell);
            freq0 = freqs0(ell);
            ghat = exp(-(tau*abs(evals(freq0)-evals(1:k))));


            % conversion to single precision
            evecs = single(evecs);
            ghat = single(ghat);

            % area term
            A = single(full(shape.A));
            rho = sqrt(sum(A(:)));

            % translation
            traslhat = bsxfun(@times,ghat,evecs');
            nn = 1./sqrt(diag(traslhat'*traslhat));
            traslhat = bsxfun(@times, nn', traslhat);
            trasl =  evecs * traslhat;

            %%
%             idx = 900;
%             figure
%             trisurf(shape.TRIV,shape.X,shape.Y,shape.Z,double(trasl(:,idx)),'SpecularStrength',0.15);
%             axis equal; axis off, view([0 0]); shading interp; colormap(jet); colorbar;
%             title(['angle ',num2str(angles(an)),' \ ',num2str(taus(ell)),' \ ',num2str(alpha(al))]);
            %%
            for s = 1:length(functions)

                f = single(fs(:,s));

                % WFT atom: modulation * translation
                q = size(f,2);
                M = cell(q,1);
                for i = 1:q
                    M{i} = rho * bsxfun(@times,evecs,f(:,i))' * A * trasl;
                end

                % copy the data
                tmp_ = cat(1,M{:});
                Sf = tmp_';
                
                % TWP
                desc_ = sum( repmat(((evals.^2)./(evals'*evals))', size(Sf,1),1).*((Sf.^2)),2);

                % stack and normalize the desc
                desc = [desc,desc_./norm(desc_)];
            end
        end
    end
end

desc = double(desc);