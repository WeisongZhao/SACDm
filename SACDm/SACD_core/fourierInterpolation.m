function [ img ] = fourierInterpolation( img, itp_fac, mirrorMode )
% USAGE: [ img ] = fourierInterpolation( img, itp_fac, padding )
% Interpolation of a 2D or 3D input image using zero padding in the Fourier
% domain. The input data can be mirrored along the lateral/axial or both
% dimensions to make the borders periodic, which reduces artifacts.
%  img     - 2D or 3D input image
%  itp_fac - Interpolation factor along each dimension [ipX,ipY] / [ipX,ipY,ipZ]
%            If a single number is given, the same factor is used for all
%            dimensions. The output is of size itp_fac.*size(img)
%  mirror - Specifies wether to use periodic mirroring of the input data or
%            not (example see fSOFI publication). The padding prevents
%            artifacts from non-periodic borders and is essential if a low
%            number of pixels is available along a specific dimension.
%            Possible values: 'none','lateral','axial','both'
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2017

itp_fac = itp_fac(:).';

if ~( numel(itp_fac) == 1 || numel(itp_fac) == ndims(img))
    error('%i interpolation factors specified. Give either one for all dimension or one per dimension!', numel(itp_fac))
end

% If all interpolation factors are 1, skip the interpolation
if sum(itp_fac == ones(1,numel(itp_fac))) == numel(itp_fac)
    return
end

if numel(itp_fac) == 1
    itp_fac = repmat(itp_fac,[1,ndims(img)]);
end
noip = (itp_fac==1); % for interpolation factors of 1, we perform neither padding nor interpolation

if nargin < 3 || isempty(mirrorMode)
    mirrorMode = 'none';
end

input_sz = size(img); % Input size

% Starting index to cut out from upsampled periodic image
sz = input_sz;
sz = sz - ~mod(sz,2);
idx = ceil(sz/2)+1 + (itp_fac-1).*floor(sz/2);

switch ndims(img)
    case 2
        switch mirrorMode
            case 'none'
                %                 img = interpft(img, itp_fac*size(img,1), 1);% Interpolate x-dir
                %                 img = interpft(img, itp_fac*size(img,2), 2);% Interpolate y-dir
                
                newsz = round(itp_fac.*size(img));
                img = fInterp_2D(img, newsz);
                return
                
            case 'lateral'
                padsize = [size(img,1)/2,size(img,2)/2];
                padsize(noip) = 0;
                img = padarray(img,ceil(padsize),'symmetric','pre');
                img = padarray(img,floor(padsize),'symmetric','post');
                
                % Fourier interpolation
                %                 img = interpft(img, itp_fac*size(img,1)-(itp_fac-1), 1);% Interpolate x-dir
                %                 img = interpft(img, itp_fac*size(img,2)-(itp_fac-1), 2);% Interpolate y-dir
                
                newsz = round(itp_fac.*size(img)-(itp_fac-1));
                img = fInterp_2D(img, newsz);
                % Cut out relevant part of padded image
                img = get_valid_part(img);
                
                return
                
            case 'axial'
                error('Padding ''axial'' only possible for 3D data.')
            case 'both'
                error('Padding ''both'' only possible for 3D data.')
            otherwise
                error('Unknown padding option ''%s''.', mirrorMode);
        end
        
    case 3
        switch mirrorMode
            case 'none'
                %                 img = interpft3D(img, itp_fac*size(img,1), 1);% Interpolate x-dir
                %                 img = interpft3D(img, itp_fac*size(img,2), 2);% Interpolate y-dir
                %                 img = interpft3D(img, itp_fac*size(img,3), 3);% Interpolate z-di
                
                newsz = round(itp_fac.*size(img));
                img = fInterp_3D(img, newsz);
                return
                
            case 'lateral'
                padsize = [size(img,1)/2,size(img,2)/2, 0];
                padsize(noip) = 0;
                img = padarray(img,ceil(padsize),'symmetric','pre');
                img = padarray(img,floor(padsize),'symmetric','post');
                
                % Fourier interpolation
                %img = interpft3D(img, itp_fac*size(img,1)-(itp_fac-1), 1);% Interpolate x-dir
                %img = interpft3D(img, itp_fac*size(img,2)-(itp_fac-1), 2);% Interpolate y-dir
                %img = interpft3D(img, itp_fac*size(img,3), 3);% Interpolate z-dir
                
                newsz = round([itp_fac(1)*size(img,1)-(itp_fac(1)-1), itp_fac(2)*size(img,2)-(itp_fac(2)-1), itp_fac(3)*size(img,3)]);
                img = fInterp_3D(img, newsz);
                % Cut out relevant part of padded image
                img = get_valid_part(img);
                return
                
            case 'axial'
                padsize = [0,0, size(img,3)/2];
                padsize(noip) = 0;
                img = padarray(img,ceil(padsize),'symmetric','pre');
                img = padarray(img,floor(padsize),'symmetric','post');
                
                % Fourier interpolation
                %                 img = interpft3D(img, itp_fac*size(img,1), 1);% Interpolate x-dir
                %                 img = interpft3D(img, itp_fac*size(img,2), 2);% Interpolate y-dir
                %                 img = interpft3D(img, itp_fac*size(img,3)-(itp_fac-1), 3);% Interpolate z-dir
                
                newsz = round([itp_fac(1)*size(img,1), itp_fac(2)*size(img,2), itp_fac(3)*size(img,3)-(itp_fac(3)-1)]);
                img = fInterp_3D(img, newsz);
                % Cut out relevant part of padded image
                img = get_valid_part(img);
                return
                
            case 'both'
                padsize = size(img)/2;
                padsize(noip) = 0;
                img = padarray(img,ceil(padsize),'symmetric','pre');
                img = padarray(img,floor(padsize),'symmetric','post');
                
                % Fourier interpolation
                %                 img = interpft3D(img, itp_fac*size(img,1)-(itp_fac-1), 1);% Interpolate x-dir
                %                 img = interpft3D(img, itp_fac*size(img,2)-(itp_fac-1), 2);% Interpolate y-dir
                %                 img = interpft3D(img, itp_fac*size(img,3)-(itp_fac-1), 3);% Interpolate z-dir
                
                newsz = round(itp_fac.*size(img)-(itp_fac-1));
                img = fInterp_3D(img, newsz);
                % Cut out relevant part of padded image
                img = get_valid_part(img);
                return
                
            otherwise
                error('Unknown padding option ''%s''.', mirrorMode);
        end
end

%% Nested functions
    function img = get_valid_part(img)
        % What part of the image to cut out depends on the dimensions
        % padding was applied to before interpolation
        %
        % doip(iDim): interpolation (with padding) was performed,
        % noip(iDim): no interpolation (with padding) was performed
        doip = ~noip;
        
        switch ndims(img)
            case 2
                if noip(1) && noip(2)
                    return
                else
                    if noip(1) && doip(2)
                        img = img(:, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1);
                    else
                        if doip(1) && noip(2)
                            img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, :);
                        else
                            if doip(1) && doip(2)
                                img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1);
                            end
                        end
                    end
                end
                return
                
            case 3
                switch mirrorMode
                    case 'lateral'
                        if noip(1) && noip(2)
                            return
                        else
                            if noip(1) && doip(2)
                                img = img(:, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,:);
                            else
                                if doip(1) && noip(2)
                                    img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, :,:);
                                else
                                    if doip(1) && doip(2)
                                        img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,:);
                                    end
                                end
                            end
                        end
                    case 'axial'
                        if doip(3)
                            img = img(:,:, idx(3):idx(3)+itp_fac(3)*input_sz(3)-1);
                        else
                            return
                        end
                        return
                        
                    case 'both'
                        if noip(3) %% No z-interpolation (with padding)
                            if noip(1) && noip(2)
                                return
                            else
                                if noip(1) && doip(2)
                                    img = img(:, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,:);
                                else
                                    if doip(1) && noip(2)
                                        img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, :,:);
                                    else
                                        if doip(1) && doip(2)
                                            img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,:);
                                        end
                                    end
                                end
                            end
                        else %% With z-interpolation
                            if noip(1) && noip(2)
                                img = img(:,:,idx(3):idx(3)+itp_fac*input_sz(3)-1);
                            else
                                if noip(1) && doip(2)
                                    img = img(:, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,idx(3):idx(3)+itp_fac(3)*input_sz(3)-1);
                                else
                                    if doip(1) && noip(2)
                                        img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, :,idx(3):idx(3)+itp_fac(3)*input_sz(3)-1);
                                    else
                                        if doip(1) && doip(2)
                                            img = img(idx(1):idx(1)+itp_fac(1)*input_sz(1)-1, idx(2):idx(2)+itp_fac(2)*input_sz(2)-1,idx(3):idx(3)+itp_fac(3)*input_sz(3)-1);
                                        end
                                    end
                                end
                            end
                        end
                end
        end
    end

end




function img_ip = fInterp_3D(img, newsz)
% Fourier interpolation of 3D image 'img' to new size 'newsz' = [nx,ny,nz]
% This is similar to performing MATLABs interpft along all dimensions
% individually, but faster.

% Fourier interpolation
sz = size(img);

%  If necessary, increase ny by an integer multiple to make ny > m.
if sum(newsz == 0) >= 1
    img_ip = [];
    return
end

isgreater = newsz >= sz;
incr = zeros(3,1);
for iDim = 1:3
    if isgreater(iDim)
        incr(iDim) = 1;
    else
        incr = floor(sz(iDim)/newsz(iDim)) + 1;
        newsz(iDim) = incr(iDim)*newsz(iDim);
    end
end

img_ip = zeros(newsz);
nyqst = ceil((sz+1)/2);
img = newsz(1)/sz(1)*newsz(2)/sz(2)*newsz(3)/sz(3)* fftn(img);%  multiplicative factor conserves the counts at the original positions

% zero padding, need to copy all 8 edges of the image cube
% note: xl:'x low', xh:'x high'
img_ip(1:nyqst(1),1:nyqst(2),1:nyqst(3)) = img(1:nyqst(1),1:nyqst(2),1:nyqst(3)); % xl, yl, zl
img_ip(end-(sz(1)-nyqst(1))+1:end, 1:nyqst(2), 1:nyqst(3)) = img(nyqst(1)+1:sz(1),1:nyqst(2),1:nyqst(3)); % xh, yl, zl
img_ip(1:nyqst(1),end-(sz(2)-nyqst(2))+1:end,1:nyqst(3)) = img(1:nyqst(1), nyqst(2)+1:sz(2) ,1:nyqst(3)); % xl, yh, zl
img_ip(1:nyqst(1),1:nyqst(2),end-(sz(3)-nyqst(3))+1:end) = img(1:nyqst(1),1:nyqst(2), nyqst(3)+1:sz(3)); % xl, yl, zh
img_ip(end-(sz(1)-nyqst(1))+1:end, end-(sz(2)-nyqst(2))+1:end, 1:nyqst(3)) = img(nyqst(1)+1:sz(1), nyqst(2)+1:sz(2),1:nyqst(3)); % xh, yh, zl
img_ip(end-(sz(1)-nyqst(1))+1:end, 1:nyqst(2), end-(sz(3)-nyqst(3))+1:end) = img(nyqst(1)+1:sz(1),1:nyqst(2), nyqst(3)+1:sz(3)); % xh, yl, zh
img_ip(1:nyqst(1), end-(sz(2)-nyqst(2))+1:end, end-(sz(3)-nyqst(3))+1:end) = img(1:nyqst(1),nyqst(2)+1:sz(2),nyqst(3)+1:sz(3)); % xl, yh, zh
img_ip(end-(sz(1)-nyqst(1))+1:end, end-(sz(2)-nyqst(2))+1:end, end-(sz(3)-nyqst(3))+1:end) = img(nyqst(1)+1:sz(1),nyqst(2)+1:sz(2),nyqst(3)+1:sz(3)); % xh, yh, zh

rm = rem(sz,2);
if rm(1) == 0 && newsz(1) ~= sz(1)
    img_ip(nyqst(1), :, :) = img_ip(nyqst(1), :, :)/2;
    img_ip(nyqst(1) + newsz(1)-sz(1), :, :) = img_ip(nyqst(1), :, :);
end
if rm(2) == 0 && newsz(2) ~= sz(2)
    img_ip(:, nyqst(2), :) = img_ip(:, nyqst(2), :)/2;
    img_ip(:, nyqst(2) + newsz(2)-sz(2), :) = img_ip(:, nyqst(2), :);
end
if rm(3) == 0 && newsz(3) ~= sz(3)
    img_ip(:, :, nyqst(3)) = img_ip(:, :, nyqst(3))/2;
    img_ip(:, :, nyqst(3) + newsz(3)-sz(3)) = img_ip(:, :, nyqst(3));
end

img_ip = real(ifftn(img_ip));
% Skip points if neccessary
img_ip = img_ip(1:incr(1):newsz(1), 1:incr(2):newsz(2), 1:incr(3):newsz(3));
end



function img_ip = fInterp_2D(img, newsz)
% Fourier interpolation of 2D image 'img' to new size 'newsz' = [nx,ny]
% This is similar to performing MATLABs interpft along all dimensions
% individually, but faster.

% Fourier interpolation
sz = size(img);

%  If necessary, increase ny by an integer multiple to make ny > m.
if sum(newsz == 0) >= 1
    img_ip = [];
    return
end

isgreater = newsz >= sz;
incr = zeros(2,1);
for iDim = 1:2
    if isgreater(iDim)
        incr(iDim) = 1;
    else
        incr = floor(sz(iDim)/newsz(iDim)) + 1;
        newsz(iDim) = incr(iDim)*newsz(iDim);
    end
end

img_ip = zeros(newsz);
nyqst = ceil((sz+1)/2);
img = newsz(1)/sz(1)*newsz(2)/sz(2)* fft2(img);%  multiplicative factor conserves the counts at the original positions

% zero padding, need to copy all 4 edges of the image plane
% note: xl:'x low', xh:'x high'
img_ip(1:nyqst(1),1:nyqst(2)) = img(1:nyqst(1),1:nyqst(2)); % xl, yl
img_ip(end-(sz(1)-nyqst(1))+1:end, 1:nyqst(2)) = img(nyqst(1)+1:sz(1),1:nyqst(2)); % xh, yl
img_ip(1:nyqst(1),end-(sz(2)-nyqst(2))+1:end) = img(1:nyqst(1), nyqst(2)+1:sz(2)); % xl, yh
img_ip(end-(sz(1)-nyqst(1))+1:end, end-(sz(2)-nyqst(2))+1:end) = img(nyqst(1)+1:sz(1), nyqst(2)+1:sz(2)); % xh, yh, zl

rm = rem(sz,2);
if rm(1) == 0  && newsz(1) ~= sz(1)
    img_ip(nyqst(1), :) = img_ip(nyqst(1), :)/2;
    img_ip(nyqst(1) + newsz(1)-sz(1), :) = img_ip(nyqst(1), :);
end
if rm(2) == 0  && newsz(2) ~= sz(2)
    img_ip(:, nyqst(2)) = img_ip(:, nyqst(2))/2;
    img_ip(:, nyqst(2) + newsz(2)-sz(2)) = img_ip(:, nyqst(2));
end

img_ip = real(ifft2(img_ip));
% Skip points if neccessary
img_ip = img_ip(1:incr(1):newsz(1), 1:incr(2):newsz(2));
end


function y = interpft3D(x,ny,dim)
% !! Patched Mathworks interpft function which works on 3D data !!
%
%
%INTERPFT 1-D interpolation using FFT method.
%   Y = INTERPFT(X,N) returns a vector Y of length N obtained
%   by interpolation in the Fourier transform of X.
%
%   If X is a matrix, interpolation is done on each column.
%   If X is an array, interpolation is performed along the first
%   non-singleton dimension.
%
%   INTERPFT(X,N,DIM) performs the interpolation along the
%   dimension DIM.
%
%   Assume x(t) is a periodic function of t with period p, sampled
%   at equally spaced points, X(i) = x(T(i)) where T(i) = (i-1)*p/M,
%   i = 1:M, M = length(X).  Then y(t) is another periodic function
%   with the same period and Y(j) = y(T(j)) where T(j) = (j-1)*p/N,
%   j = 1:N, N = length(Y).  If N is an integer multiple of M,
%   then Y(1:N/M:N) = X.
%
%   Example:
%      % Set up a triangle-like signal signal to be interpolated
%      y  = [0:.5:2 1.5:-.5:-2 -1.5:.5:0]; % equally spaced
%      factor = 5; % Interpolate by a factor of 5
%      m  = length(y)*factor;
%      x  = 1:factor:m;
%      xi = 1:m;
%      yi = interpft(y,m);
%      plot(x,y,'o',xi,yi,'*')
%      legend('Original data','Interpolated data')
%
%   Class support for data input x:
%      float: double, single
%
%   See also INTERP1.

%   Robert Piche, Tampere University of Technology, 10/93.
%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.15.4.5 $  $Date: 2010/08/23 23:11:51 $

error(nargchk(2,3,nargin,'struct'));

if nargin==2,
    [x,nshifts] = shiftdim(x);
    if isscalar(x), nshifts = 1; end % Return a row for a scalar
elseif nargin==3,
    perm = [dim:max(length(size(x)),dim) 1:dim-1];
    x = permute(x,perm);
end

siz = size(x);
[m,n,l] = size(x);
if ~isscalar(ny)
    error(message('MATLAB:interpft:NonScalarN'));
end

%  If necessary, increase ny by an integer multiple to make ny > m.
if ny > m
    incr = 1;
else
    if ny==0, y=[]; return, end
    incr = floor(m/ny) + 1;
    ny = incr*ny;
end
a = fft(x,[],1);
nyqst = ceil((m+1)/2);
b = [a(1:nyqst,:,:) ; zeros(ny-m,n,l) ; a(nyqst+1:m,:,:)];
if rem(m,2) == 0
    b(nyqst,:,:) = b(nyqst,:,:)/2;
    b(nyqst+ny-m,:,:) = b(nyqst,:,:);
end
y = ifft(b,[],1);
if isreal(x), y = real(y); end
y = y * ny / m;
y = y(1:incr:ny,:,:);  % Skip over extra points when oldny <= m.

if nargin==2,
    y = reshape(y,[ones(1,nshifts) size(y,1) siz(2:end)]);
elseif nargin==3,
    y = ipermute(y,perm);
end

end

