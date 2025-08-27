


clear, clc, close all
addpath('/home/admin/Documents/YanLab/Task_Local/Matlab_shared');


%% Read configure files
config.maskpath = './structure_262_RT_50um.nrrd';
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

func_matrix2tiff(mask, 16, '/tmp/RT_anno.tif');



%% Read path of neurons
fid = fopen('./Hemisphere_Path_PT_CT_RT2.csv');
fdata = textscan(fid, '%s %s', 'Delimiter', ' ');
fclose(fid);

Hemisphere = fdata{1,1};
Path_SWC = fdata{1,2};


%% Calculate projection pattern
NumNeuron = length(Path_SWC);
NumCube = length(I);

matrix_projection = zeros(NumNeuron, NumCube);

parpool(8)
parfor iNeuron = 1:1:NumNeuron
    tmpPathSWC = Path_SWC{iNeuron,1};
    swc_allen = dlmread(tmpPathSWC(1, 2:end-1));
    
    tmpHemisphere = Hemisphere{iNeuron,1};
    if strcmp(tmpHemisphere(1, 2:end-1), 'Right')
        swc_allen(:,5) = 456*25 - swc_allen(:,5);
    end
    
    [tmp_data, type] = func_calProjectionPattern(swc_allen, config.resolution, mask, ...
                                                 config.metric, config.normalize);
    
    matrix_projection(iNeuron, :) = tmp_data';
    
    if(0 == rem(iNeuron, 100))
        disp(strcat(num2str(iNeuron),'/',num2str(NumNeuron)));
    end
end


%% remove cubes with few projections
sum_cube = sum(matrix_projection, 1);
map_cube = 1:1:size(matrix_projection, 2);

indexCube = find(sum_cube > 100); % TODO PFC projection 

map_cube_PFC = map_cube(indexCube);
map_cube_nonPFC = map_cube(find(sum_cube < 100));

matrix_projection_PFC = matrix_projection(:,indexCube);


%% Write projection
dlmwrite('./PT_CT_RT_projection_strength_100.txt', matrix_projection_PFC, 'precision', 16, 'delimiter',' ');



%% Step 04 write TIFF files
T = dlmread('./label_cluster_real_K6_2.csv');

map_cube_PFC = label';

NumCluster = 7;

out_mask = zeros(size(mask), 'uint8');
color_code = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
for i = 1:1:NumCluster
    out_mask(I(i==T)) = color_code(i);
end

%out_mask(I(map_cube_nonPFC)) = 6;

func_matrix2tiff(out_mask, 8, './label_cluster_real_K6_2.tif');


%% save each channel as a seperate image
for i = 1:1:5
    out_mask = zeros(size(mask), 'uint16');
    out_mask(I(map_cube_PFC(i==T))) = 65535;
    func_matrix2tiff(out_mask(:,:,:), 16, strcat('./channel/', num2str(i), '.tif'));
end
%% Crop mask 

template = nrrdread('/home/admin/Documents/YanLab/resource_structure_allen_50um/average_template_50.nrrd');
mask_CP_ACB = nrrdread('./structure_262_RT_50um.nrrd');
mask_CP_ACB(:,:,115:end) = 0;

template_CP_ACB = template;
template_CP_ACB(mask_CP_ACB == 1) = 0;

func_matrix2tiff(template_CP_ACB(:, :, :), 16, './template_RT.tiff');



