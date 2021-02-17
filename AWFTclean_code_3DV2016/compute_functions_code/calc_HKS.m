function [HKS,t] = calc_HKS(M, t)

if nargin==1
    % based on J.Sun's code
    tmin = abs(4*log(10) / M.evals(end));
    tmax = abs(4*log(10) / M.evals(2));
    nstep = 99; %100;
    stepsize = (log(tmax) - log(tmin)) / nstep;
    logts = log(tmin):stepsize:log(tmax);
    t = exp(logts);
    %t = 2.^(1: 1/16 : 18-1/16);
end

HKS = zeros(M.n, length(t));

% skipping the first freq. as the authors do in their code
for i=1:length(t)
    HKS(:,i) = sum(...
        (M.evecs(:,2:end)).^2 .* repmat(exp(-t(i)*M.evals(2:end))',M.n,1), ...
        2);
end

end
