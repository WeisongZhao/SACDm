function y=generate_rsf(gama,n)
sigma = gama./sqrt(8*log(2));
kernelRadius = min(ceil(sigma * sqrt(-2 * log(0.0002)))+1,floor(n/2));
ii=-kernelRadius:kernelRadius;
rsf_x=1/2*(erf((ii+0.5)./(sqrt(2).*sigma)) - erf((ii-0.5)./(sqrt(2).*sigma)));
kernel = rsf_x'* rsf_x;
y=kernel./sum(kernel(:));
end