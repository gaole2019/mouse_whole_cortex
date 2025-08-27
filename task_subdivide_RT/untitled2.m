% project mouse whole cortex
% subdivide RT
% generate mesh of subdomains using FIJI


p = genpath('/home/admin/Documents/YanLab/Task_Local/Matlab_shared');
addpath(p);


fileID = fopen("/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/K7/itksnap_label.txt");
C = textscan(fileID, '%d %d %d %d %d %d %d %s', "Delimiter","\t", "CommentStyle","#");
fclose(fileID);

structureID = C{1};
structureID = structureID(2:end,1);

NumStructure = length(structureID);

%anno_structures = nrrdread("/run/media/admin/LeGAO/Project_Human_Hypothalamus/brain_sample_3/HE_images/HE_raw_clean_ds0.2_align/v3_anno_structures.nrrd");

%%label_structures = readtable("v3_anno_structures.label", "FileType","text", "CommentStyle","#");
%%label_structures.Properties.VariableNames = {'IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'};

anno_structures = ...
    nrrdread("/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/K7/label_cluster_real_K7.nrrd");


for iStructure = 1:1:NumStructure
    tmp_structure = zeros(size(anno_structures), "uint8");
    tmp_structure(anno_structures == structureID(iStructure)) = 255;

func_matrix2tiff(...
    tmp_structure, 8, ...
    strcat('/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/K7/structures_mask/', ...
           num2str(structureID(iStructure)), '.tiff'));
end


