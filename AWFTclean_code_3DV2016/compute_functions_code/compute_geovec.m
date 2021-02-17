function [g0] = compute_geovec(shape)

        rho2 =  sum(sum(full(diag(shape.A))));
        Lambda = rho2.*abs(shape.evals);
        lambdamin = 0;
        lambdamax = Lambda(end);
        ncoefs = 150; 
        delta = (lambdamax - lambdamin)/(ncoefs-1);
        LAMBDA = lambdamin:delta:lambdamax;
        sigma = delta;
        B0 = basisfunc(Lambda, LAMBDA, sigma);
        g0 =(rho2.*shape.evecs.^2)*B0;
        g = g0(:,1:16);
end