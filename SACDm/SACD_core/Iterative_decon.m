function data_de=Iterative_decon(data,kernel,iteration,rule,gpu)

if nargin < 3 || isempty(iteration)
    iteration=10;
end
if nargin < 4 || isempty(rule)
    rule=1;
end
if nargin < 5 || isempty(gpu)
    gpu = 0;
end
if ndims(data)==3
    for i=1:size(data,3)
        data_de(:,:,i)=real(RL_core(data(:,:,i),kernel,iteration,rule,gpu));
    end
else
    data_de=real(RL_core(data,kernel,iteration,rule,gpu));
end