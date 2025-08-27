% project mouse whole cortex
% figure 4. cortico-basal ganglia projection
% subdivide SNr

%% 1. assign each cube in SNr an unique value
clear, clc, close all
addpath('/home/admin/Documents/YanLab/Task_Local/Matlab_shared');
%% Read configure files
config.maskpath = './structure_381_SNr.nrrd';
config.resolution = 50;
config.metric = 'LengthOfAxon';
config.normalize = 'false';
%% Generate mask
mask = nrrdread(config.maskpath);
[row, col, slice] = size(mask);
tmp_left = ones(row, col, slice/2, 'uint8');
tmp_right = ones(row, col, slice/2, 'uint8')*0;
MASK_HEMI = cat(3, tmp_left, tmp_right);
mask = mask.*MASK_HEMI;
% Assign ID to each cube
mask = uint16(mask);
I = find(1 == mask);
label = 1:1:length(I);
mask(I) = label';
func_matrix2tiff(mask, 16, './SNr_anno.tif');