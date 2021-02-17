function I = basisfunc(lambda, LAMBDA, sigma)

I = exp(-bsxfun(@minus, lambda(:), LAMBDA(:)').^2/2/sigma^2);


