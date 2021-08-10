function SACDresult = SACDm(imgstack, varargin)

%pixel     pixel size {example:65*10^-9}
%wavelength    wavelength {example:532*10^-9}
%iter1      Maximum deconvolution iterations {example:20}
%mag      magnification times  {example:8}
%scale    PSF^2 {example:2}
%order     Auto-correlation cumulant order  {example:2}

params.pixel = 65;
params.wavelength = 525;
params.NA = 1.3;
params.mag = 2;
params.iter1 = 10;
params.iter2 = 7;

params.subfactor = 0.8;
params.ACorder = 2;
params.scale = [];
params.psf = [];
params.resolution = [];
params.ifregistration = false;
params.ifbackground = false;
params.backgroundfactor = 2;
params.ifsparsedecon = false;
params.fidelity = 100;
params.tcontinuity = 0.1;
params.sparsity = 1;

warning('off');
addpath('./Utils');
addpath('./Sparse');
addpath('./SACD_core');
addpath('./Register');

if nargin > 2
    params =  read_params(params, varargin);
end
if isempty(params.scale)
    params.scale = params.ACorder;
end
%% psf calculation
tic
disp(['SACD reconstruction start...'])
disp(['Data size is ' num2str(size(imgstack,1)) ...
    ' X ' num2str(size(imgstack,2)) ' X ' num2str(size(imgstack,3))]);
params.pixel = params.pixel * 10^-9;
params.wavelength = params.wavelength * 10^-9;
disp(['PSF calculation...'])
if isempty(params.psf)&&isempty(params.resolution)
    psf = kernel(params.pixel, ...
        params.wavelength, params.NA, 0, min(size(imgstack,1),size(imgstack,2)));
    psfv2 = kernel(params.pixel/params.mag, ...
        params.wavelength, params.NA, 0, min(size(imgstack,1),size(imgstack,2)));
elseif ~isempty(params.psf)
    psf = params.psf;
    psfv2 = fourierInterpolation(params.psf,[params.mag,params.mag,1],'lateral');
elseif ~isempty(params.resolution)
    params.resolution = params.resolution * 10^-9;
    psf = generate_rsf(params.resolution/params.pixel,min(size(imgstack,1),size(imgstack,2)));
    psfv2 = generate_rsf(params.mag * params.resolution/params.pixel,...
        min(size(imgstack,1),size(imgstack,2)));
end
%% background subtration
t = size(imgstack,3);
imgstack = single(imgstack);
for i = 1 : t
    stack(:,:,i) = imgstack(:,:,i) - min(min(imgstack(:,:,i)));
end
if params.ifbackground
    disp(['Background estimation...'])
    BACK=background_estimation(stack./params.backgroundfactor);
    BACK(BACK < 0) = 0;
    stack = stack - BACK;
    stack(stack < 0) = 0;
    clear BACK
end
if params.ifregistration
    disp(['Drift estimation...'])
    for i = 2:t
        stack(:,:,i)=register(stack(:,:,i),stack(:,:,1));
    end
end
%% pre deconvolution
disp(['Pre RL deconvolution...'])
for i = 1:t
    % datadecon(:,:,i)=datasub(:,:,i);
    datadecon(:,:,i)=deconvlucy(stack(:,:,i),psf,params.iter1);
end
%% Fourier interpolation
disp(['Fourier interpolation...'])
datadeconl=fourierInterpolation(datadecon,[params.mag,params.mag,1],'lateral');
datadeconl(datadeconl < 0) = 0;
%% AC cumulant
disp(['Cumulant calculation...'])
datadeconl=abs(datadeconl - params.subfactor * mean(datadeconl,3));
cum = abs(cumulant(datadeconl,params.ACorder));
if params.ifsparsedecon
    disp(['Post sparsity-continuity reconstruction...'])
    cum = SparseHessian_core(cum,params.fidelity,params.tcontinuity,params.sparsity,100,0);
end
disp(['Post RL deconvolution...'])
SACDresult = deconvlucy(cum,psfv2.^params.scale, params.iter2);
ttime = toc;
disp('SACD reconstruction completed.');
disp(['Total time cost: ', num2str(ttime/60) ' mins'])
