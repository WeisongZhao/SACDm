function [ shiftedIm ] =  ShiftImage( Im, Shift )
%   [ shiftedIm ] = ShiftImage( Im, Shift )
%   The function shifts image of a given quantity in pixels in the frequency domain

assert(size(Im,3)==1,'size(Im,3)~=1');

assert(size(Shift,1) == 1 && size(Shift,2) == 2 ,'size(Shift,1) ~= 1 || size(Shift,2) ~= 2');

[h, w] = size(Im);

%Define shift in frequency domain
[xF,yF] = meshgrid(-floor(h/2):-floor(h/2)+h-1,-floor(w/2):-floor(w/2)+w-1);

fftIm = fftshift(fft2(Im));

%Perform the shift and antitransform
shiftedIm = real(ifft2(ifftshift(fftIm.*exp(-1i*2*pi.*((xF*Shift(1))/h+(yF*Shift(2))/w)))));

end