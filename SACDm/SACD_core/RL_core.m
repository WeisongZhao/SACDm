function data_decon=RL_core(data,kernel,iteration,rule,gpu)

if nargin < 3 || isempty(iteration)
    iteration = 10;
end
if nargin < 4 || isempty(rule)
    rule = 1;
end
if nargin < 5 || isempty(gpu)
    gpu = 0;
end
%initialization
kernel=kernel./sum(sum(kernel));
[dx,dy]=size(data);
B=floor(min(dx,dy)/6);
data=padarray(data, [B,B] ,'replicate');
kernelm=kernel;
kernel=zeros(size(data,1),size(data,2));
kernel(1:size(kernelm,1),1:size(kernelm,2)) = kernelm;
[~,idx] = max(kernel(:));
[x,y,~] = ind2sub(size(kernel), idx);
kernel = circshift(kernel, 1-[x,y,1]);
if gpu==0
    yk=data;
    xk=zeros(size(data));
    vk=zeros(size(data));
else
    kernel=gpuArray(kernel);
    data=gpuArray(data);
    yk=data;
    xk=gpuArray.zeros(size(data));
    vk=xk;
end
otf = fftn(kernel,size(data));
% otf=psf2otf(kernel,size(data));
%%
rliter = @(estimate, data, otf)fftn(data./ max(ifftn(otf .* fftn(estimate)), 1e-6));
for iter = 1:iteration
    xk_update = xk;
    xk= yk .* real(ifftn(conj(otf) .* rliter(yk,data,otf)))...
        ./real(ifftn(fftn(ones(size(data))) .* otf));
    xk=max(xk,1e-6);
    vk_update = vk;
    vk=max(xk - yk,1e-6);
    if iter == 1
        alpha = 0;
    else
        alpha = sum(sum(vk_update.*vk))./...
            (sum(sum(vk_update.*vk_update))+eps);
        alpha = max(min(alpha,1),1e-6);
    end
    yk = xk + alpha * (xk-xk_update);
    yk=max(yk,1e-6);
    yk(isnan(yk)) = 1e-6;
end
yk(yk < 0) = 0;
data_decon = yk(B+1:end-B, B+1:end-B);