function desc = compute_WFTdesc(shape,functions,param)
%
% desc = compute_WFTdesc(shape,functions_names,parameters)
%    computes the windowed Fourier tranform (WFT) of a shape
% input:
%       shape = struct with shape.X, shape.Y, shape.Z, shape.TRIV: 
%       functions = cell where functions{i} = string with the name
%       params,
%
% output:
%    desc: WFT descriptor on the input shape
%

freqs0 = param.freqs0;
taus = param.taus; 
options.curv_smooth = param.curv_smooth;
options.n_eigen = param.n_eigen;
k = param.n_eigen;

% compute isotropic LBO
options.alpha = 0;
options.angle = 0;
[shape.Phi, shape.Lambda, shape.W, shape.A] = calc_ALB([shape.X,shape.Y,shape.Z], shape.TRIV, options);

% normalize data
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
    freq0 = freqs0(ell);
    tau = sum(diag(full(shape.A)))*taus(ell);
    ghat = exp(-(tau*abs(shape.Lambda(freq0)-shape.Lambda(1:k))));
    ghat = ghat./norm(ghat,2);
         
    % area term
    A = single(full(shape.A));
    rho = sqrt(sum(A(:)));
    
    % conversion to single precision
    Phi = single(shape.Phi);
    ghat = single(ghat);

    for s = 1:length(functions)
        
        f = single(fs(:,s));
        
        % translation
        traslhat = bsxfun(@times,ghat,Phi');
        nn = 1./sqrt(diag(traslhat'*traslhat));
        traslhat = bsxfun(@times, nn', traslhat);
        trasl =  Phi * traslhat;

        % WFT atom: modulation * translation
        q = size(f,2);
        M = cell(q,1);
        for i = 1:q
            M{i} = rho * bsxfun(@times,Phi,f(:,i))' * A * trasl;
        end

        % copy the data
        tmp_ = cat(1,M{:});
        Sf = tmp_';
        desc_ = sum( repmat(((shape.Lambda.^2)./(shape.Lambda'*shape.Lambda))', size(Sf,1),1).*((Sf.^2)),2);

        % stack and normalize the desc
        desc = [desc,desc_./norm(desc_)];
    end
end

desc = double(desc);