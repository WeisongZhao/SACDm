function GaussianParameters  =  G2DFit( data, iterationsNumber, pixelSize)
% G2DFit_gaussian2DFittingSupervised Fits 2D data to a gaussian function.
% Unsupervised version, estimating the rotation angle (theta).
%   data = data to fit
%   iterationsNumber = maximum number of iterations
%   pixelSize = conversion factor from pixel # to units
% Unupervised version, estimating the rotation angle.
% it starts from a initial hypothesis of center and sigma that can be
% passed as argument, or automatically estimated by the algorithm.

% arguments processing
if (iterationsNumber < 1)
    iterationsNumber = 1;
end
if (pixelSize <= 0)
    pixelSize = 1;
end

% constants
coincidenceLevel = 0.01;
minimumCondNumber = 1e-20;

%initialization
data = double(data);
[height, width] = size(data);
backgroudLevel = max( min(min(data)), 0.0 );
coeff = zeros(6,1,'single'); oldCoeff = coeff;

% iterative procedure
for currentIteration = 1 : iterationsNumber
    
    % A = [  N,   sx,   sy,   sx2,   sy2, sxy;
    %       sx,  sx2,  sxy,   sx3,  sxy2, sx2y;
    %       sy,  sxy,  sy2,  sx2y,   sy3, sxy2;
    %      sx2,  sx3, sx2y,   sx4, sx2y2, sx3y;
    %      sy2, sxy2,  sy3, sx2y2,   sy4, sxy3
    %      sxy, sx2y, sxy2, sx3y, sxy3, sx2y2 ];
    % b = [ lg;  xlg;  ylg;  x2lg;  y2lg];
    A = zeros(6,6,'single');
    b = zeros(6,1,'single');
    
    for yC = 1:height
        for xC = 1:width
            value = max( data(yC,xC) - backgroudLevel, 0.0 );
            x = (xC - 1) * pixelSize;
            y = (yC - 1) * pixelSize;
            if (currentIteration>1)
                oldValue = exp(sum(coeff'.*[1 x y x^2 y^2 x*y]));
            else
                oldValue = value;
            end
            %if (value > backgroudLevel)
            if (value > 0)
                dA = [  1,   x,     y,     x^2,     y^2,     x*y;
                    x,   x^2,   x*y,   x^3,     x*y^2,   x^2*y;
                    y,   x*y,   y^2,   x^2*y,   y^3,     x*y^2;
                    x^2, x^3,   x^2*y, x^4,     x^2*y^2, x^3*y;
                    y^2, x*y^2, y^3,   x^2*y^2, y^4,     x*y^3
                    x*y, x^2*y, x*y^2, x^3*y,   x*y^3,   x^2*y^2];
                db = [  log(value); x*log(value); y*log(value);
                    x^2*log(value); y^2*log(value); x*y*log(value)];
                A = A + (dA * oldValue^2);
                b = b + (db * oldValue^2);
            end
        end
    end
    
    condNum = rcond(A);
    if (isnan(condNum) || condNum < minimumCondNumber)
        GaussianParameters.sx = -1;
        GaussianParameters.sy = -1;
        GaussianParameters.ux = -1;
        GaussianParameters.uy = -1;
        GaussianParameters.fwhmX = -1;
        GaussianParameters.fwhmY = -1;
        GaussianParameters.peak = -1;
        GaussianParameters.kValue = -1;
        GaussianParameters.rSquare = -1;
        GaussianParameters.background = -1;
        GaussianParameters.iterations = -1;
        GaussianParameters.theta = -1;
        return;
    end
    
    coeff = A\b;
    
    % Evaluation of the model
    tailDistanceSqr = -4.5 / min( coeff(4), coeff(5) );
    hypX = -0.5 * coeff(2) / coeff(4);
    hypY = -0.5 * coeff(3) / coeff(5);
    modelData = zeros(height, width,'single');
    meanError = 0.0; divisor = 0.0;
    for yC = 1:height
        for xC = 1:width
            x = (xC - 1) * pixelSize;
            y = (yC - 1) * pixelSize;
            modelData(yC,xC) = exp(sum(coeff'.*[1 x y x^2 y^2 x*y]));
            if ( (x-hypX)^2 + (y-hypY)^2 >= tailDistanceSqr )
                meanError = meanError + ( data(yC, xC) - modelData(yC,xC) );
                divisor = divisor + 1;
            end
        end
    end
    assert( meanError == 0 || divisor>0 );
    if (divisor > 0)
        meanError = meanError ./ divisor;
    end
    backgroudLevel = meanError;
    if (backgroudLevel <= 0)
        backgroudLevel = 0;
    end
    
    % Stopping criterion: each coeff doesn't change for more than 0.01.
    if ( max( abs(coeff - oldCoeff) ) <= coincidenceLevel )
        break;
    end
    
    oldCoeff = coeff;
    
end

if ( coeff(2)<=0 || coeff(3)<=0 || coeff(4)>=0 || coeff(5)>=0 )
    GaussianParameters.sx = -1;
    GaussianParameters.sy = -1;
    GaussianParameters.ux = -1;
    GaussianParameters.uy = -1;
    GaussianParameters.fwhmX = -1;
    GaussianParameters.fwhmY = -1;
    GaussianParameters.peak = -1;
    GaussianParameters.kValue = -1;
    GaussianParameters.rSquare = -1;
    GaussianParameters.background = -1;
    GaussianParameters.iterations = -1;
    GaussianParameters.theta = -1;
    return;
end

theta = 0.5 * atan( coeff(6) / (coeff(5) - coeff(4) ) );
% thetaDeg = theta * 180.0 / pi;
sin2theta = sin(theta)^2;
cos2theta = cos(theta)^2;
s2x = ( 0.5 * cos(2*theta) / ( coeff(5) * sin2theta - coeff(4) * cos2theta ) );
s2y = ( 0.5 * cos(2*theta) / ( coeff(4) * sin2theta - coeff(5) * cos2theta ) );
ux = coeff(2) * ( s2y * sin2theta + s2x * cos2theta ) + 0.5 * coeff(3) * (s2y - s2x) * sin(2*theta);
uy = coeff(3) * ( s2y * cos2theta + s2x * sin2theta ) + 0.5 * coeff(2) * (s2y - s2x) * sin(2*theta);
sx = s2x ^ 0.5;
sy = s2y ^ 0.5;
fwhmX = sx * ( 2.0 * sqrt ( 2.0 * log (2.0) ) );
fwhmY = sy * ( 2.0 * sqrt ( 2.0 * log (2.0) ) );
peak = exp(sum(coeff'.*[1 ux uy ux^2 uy^2 ux*uy]));
kValue = (2*pi*sx*sy)* peak;

% r-Square
data = data - backgroudLevel;
data = data .* (data>0);
rss = sum(sum((data - modelData).^2));
meanData = mean(mean(data));
tss = sum(sum((data - meanData).^2));
rSquare = 1 - (rss / tss);

% Returning the parameters
GaussianParameters.sx = sx;
GaussianParameters.sy = sy;
GaussianParameters.ux = ux;
GaussianParameters.uy = uy;
GaussianParameters.fwhmX = fwhmX;
GaussianParameters.fwhmY = fwhmY;
GaussianParameters.peak = peak;
GaussianParameters.kValue = kValue;
GaussianParameters.rSquare = rSquare;
GaussianParameters.background = backgroudLevel;
GaussianParameters.iterations = currentIteration;
GaussianParameters.theta = theta;
end