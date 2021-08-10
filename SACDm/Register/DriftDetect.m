function [ShiftInteger, Shift] =  DriftDetect( Im1 , Im2, varargin )
%   [ShiftInteger, Shift] =  DriftDetect( Im1 , Im2 )
%   Detects the drift between two squared images of the same size

assert(...
    (size(Im1,1) == size(Im2,1)) && ...
    (size(Im1,2) == size(Im2,2)), ...
    'size(Im1) ~= size(Im2)');

assert(...
    (size(Im1,3) == 1) && ...
    (size(Im2,3) == 1) ,...
    '( size(Im1,3) || size(Im2,3) ) ~= 1');

params.maxIterations = 200;
params.roiRadius = 5;

if (nargin > 2)
    params =  read_params(params, varargin);
end

[h, w] = size(Im1);

% Crosscorrelation in the frequency domain
Nr = size(Im1,1)*2-1;
Nc = size(Im1,2)*2-1;
corr = fftshift(ifft2(fft2(Im1,Nr,Nc).*conj(fft2(Im2,Nr,Nc))));

[~,SND] = max(corr(:));
[IJ,JI] = ind2sub(size(corr),SND);
ShiftInteger(1) = h-JI;
ShiftInteger(2) = w-IJ;

%  gaussianfitting, centered according to integerShift
corr = real(corr);
corr_Size = size(corr);

shiftedCorr =  ShiftImage(corr,ShiftInteger);
sub_corr = shiftedCorr( ...
    floor( corr_Size(2)/2 ) - params.roiRadius : floor(corr_Size(2)/2) + params.roiRadius    , ...
    floor( corr_Size(1)/2 ) - params.roiRadius : floor(corr_Size(1)/2) + params.roiRadius   );

%   Gaussian Fitting
GaussianParameters =  G2DFit(sub_corr, params.maxIterations, 1 );

Shift(1) = params.roiRadius  - GaussianParameters.ux + ShiftInteger(1) + 1;
Shift(2) = params.roiRadius  - GaussianParameters.uy + ShiftInteger(2) + 1;

end