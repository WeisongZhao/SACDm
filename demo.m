addpath(genpath('./SACDm'));
imgstack = imreadstack('XXXXXXXXX.tif');
SRimg = SACDm(imgstack);
imshow(SRimg,[])
