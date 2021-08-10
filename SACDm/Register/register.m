function [toreg]=register(toreg,GT)

params.driftCorrection = false;
params.roiRadius = 5;
params.maxIterations = 200;
[ox,oy]=size(toreg);
[~ , Shift]=  DriftDetect (GT, toreg, ...
    'maxIterations', params.maxIterations, 'roiRadius', params.roiRadius);  %Find shift
[xF,yF] = meshgrid(-floor(ox/2):-floor(ox/2)+ox-1,-floor(oy/2):-floor(oy/2)+oy-1); %Define shift in frequency domain
fftB = fftshift(fft2(toreg));
fftB=fftB.*exp(-1i*2*pi.*((-xF*Shift(1))/ox+(-yF*Shift(2))/oy))+eps; %Perform the shift
B=abs(real(ifft2(ifftshift(fftB))));
toreg=B./max(max(B));
fprintf(strcat('', 'SubPixel Shift:', '\n', ...
    num2str(Shift(1)),' pixels in x;' , '\n', ...
    num2str(Shift(2)),' pixels in y.', '\n' ));