% project mouse whole cortex
% cortical columns

close all 
clear

p = genpath('/home/admin/Documents/YanLab/Task_Local/Matlab_shared');
addpath(p);

%% volumes of 255 cortical columns
CC_anno_K255 = ...
    nrrdread('/run/media/admin/LeGAO/CodeMatlab/Isocortex_ROI/result/label_PFC_volume_filled_K255_R20_10um.nrrd');

CC255_volume = zeros(255,1);
for iRegion = 1:1:255
    num_pixel = length(find(CC_anno_K255 == iRegion));
    CC255_volume(iRegion, 1) = num_pixel;
end

mean = mean(CC255_volume);
sem = std(CC255_volume) / sqrt(255);


statistics_volume = readtable('/run/media/admin/LeGAO/CodeMatlab/Isocortex_ROI/result/statistics_volumes.csv');

volume = statistics_volume.NumberOfVoxels(2:end) * 25*25*25 / 1000/1000/1000;

histogram(volume);
xlabel('Volume of cortical columns (mm^3)');
ylabel('Count');

mean(volume)
std(volume)

%% volumes of 43 cortical subregions in ABA
ABA_anno = nrrdread('/home/admin/Documents/YanLab/database_neuron/ABA_annotation/annotation_10_2017.nrrd');

length(find(ABA_anno == 417))

ABA_43 = ...
    readtable('/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_cortical_column/ABA_43.csv');

volume_43 = zeros(43,1);
for iRegion = 1:1:43
    tmpIndex = find(ABA_43.RegionID == iRegion);
    tmpStructureID = ABA_43.StuctureID(tmpIndex);
    
    tmpLength = length(tmpStructureID);
    num_pixel = 0;
    for iSubregion = 1:1:tmpLength
        num_pixel = num_pixel + length(find(ABA_anno == tmpStructureID(iSubregion)));
    end
    volume_43(iRegion, 1) = num_pixel;
end

volume_43_real = volume_43/2*10*10*10 / 1000 / 1000/1000;
mean(volume_43_real)
std(volume_43_real) / sqrt(43)

