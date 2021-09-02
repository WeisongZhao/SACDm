clear;clc;close all;
%% Involve SACDm
addpath(genpath('./SACDm'));
%% Read data
imgstack = imreadstack('561 scmos-30ms-C1_2020-09-13_2-ROI.tif');
%% SACD recon
SRimg = SACDm(imgstack,'pixel',65,'NA',1.3,'wavelength',561);
%% Visualization
LRimg = imfilter(mean(double(imgstack),3),generate_rsf(2));
LRimg = LRimg./max(LRimg(:));
figure(1);imshow(LRimg,'colormap',hot)
background = 0.02;order = 2;
SRimg2vis = SRimg.^0.5;
SRimg2vis(SRimg2vis < order * background * max(SRimg2vis(:))) = 0;
figure(2);imshow(SRimg2vis,[],'colormap',hot)
