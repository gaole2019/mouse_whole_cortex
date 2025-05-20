# project mouse whole cortex
# figure 4. cortico-basal ganglia projections
# joint projections of IT neurons

#
# target_acronym <-
#   c("STR", "GPe", "STN", "GPi", "SNr")
#
# proj_IT_target_acronym <-
#   c(paste0("l_", target_acronym),
#     paste0("r_", target_acronym))
#
#
# proj_IT_BG <- shared_proj_Combined[
#   index_axon_IT,
#   snp::getStructureFromAcronym(proj_IT_target_acronym)]
#
# colnames(proj_IT_BG) <- proj_IT_target_acronym
#
#
# threshold_region <- c(10000, 1000, 100, 100, 1000,
#                       10000, 1000, 100, 100, 1000)

encode_neuron <- matrix(nrow = nrow(proj_IT_BG_new), ncol = 1)
for (iNeuron in 1:nrow(proj_IT_BG_new)) {
  encode <- 0
  for (iTarget in 1:ncol(proj_IT_BG_new)) {
    encode <- encode +
      ifelse(proj_IT_BG_new[iNeuron, iTarget] > threshold_region[iTarget], 1, 0) * 2^(iTarget - 1)
  }
  encode_neuron[iNeuron, 1] <- encode
}

aaaa <- sort(table(encode_neuron), decreasing = T)

barplot(aaaa[2:15])


vector_names <- matrix(nrow = length(aaaa), ncol = 1)
for (iPatttern in 1:length(aaaa)) {
  binary_pattern <- as.integer(intToBits(as.numeric(names(aaaa)[iPatttern])))[1:10]
  tmp_name <- ""
  for (iTarget in 1:ncol(proj_IT_BG_new)) {
    if (binary_pattern[iTarget] == 1) {
      tmp_name <- paste0(tmp_name, "_", colnames(proj_IT_BG_new)[iTarget])
    }
  }

  vector_names[iPatttern, 1] <- tmp_name
}




bbbb <- proj_IT_BG_new
bbbb$l_STR <- bbbb$l_STR > 1000
bbbb$l_GPe <- bbbb$l_GPe > 1000
bbbb$l_STN <- bbbb$l_STN > 100
bbbb$l_GPi <- bbbb$l_GPi > 100
bbbb$l_SNr <- bbbb$l_SNr > 1000

bbbb$r_STR <- bbbb$r_STR > 1000
bbbb$r_GPe <- bbbb$r_GPe > 1000
bbbb$r_STN <- bbbb$r_STN > 100
bbbb$r_GPi <- bbbb$r_GPi > 100
bbbb$r_SNr <- bbbb$r_SNr > 1000


bbbb[bbbb == TRUE] <- 1
bbbb[bbbb == FALSE] <- 0
bbbb <- as.data.frame(bbbb)

library(ComplexHeatmap)
m = make_comb_mat(bbbb)
m = m[comb_size(m) > 20 & comb_degree(m) > 0]
UpSet(m, comb_order = order(comb_size(m), decreasing = T), pt_size = unit(4, "mm"))





# IT neurons projecting to iSTR and iGPe ----

table(info_tmp$Region[index_axon_IT[which(encode_neuron == 3)]])


library(rgl)
mesh_lr_GPe <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1022_GPe.obj")
mesh_lr_STR <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/485_STRd.obj")

mesh_l_GPe <- getMeshLeft(mesh_lr_GPe)
mesh_l_STR <- getMeshLeft(mesh_lr_STR)


which(proj_IT_BG$l_STR[which(encode_neuron == 3)] > 1000 &
        proj_IT_BG$l_GPe[which(encode_neuron == 3)] > 4000)
# [1]  97 131 163 171 185 281 351 352 356 383 456 465 533 539

index_axon_IT_BG_3 <-
  index_axon_IT[which(encode_neuron == 3)][42]


nl_axon_IT_BG_3 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_IT_BG_3],
    neuronnames = shared_NeuronName[index_axon_IT_BG_3])

#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])
library(rgl)
open3d()
snp::plot_morphology(
  neuronList = nl_axon_IT_BG_3,
  hemisphere = info_tmp$Hemisphere[index_axon_IT_BG_3],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")

mesh_root <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")
shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "blue", alpha = 0.1)
shade3d(mesh_l_STR, color = "red", alpha = 0.1)



rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)



# IT neurons projecting to iSTR and rSTR (33) ----

table(info_tmp$Region[index_axon_IT[which(encode_neuron == 33)]])


which(proj_IT_BG$l_STR[which(encode_neuron == 33)] > 40000 &
        proj_IT_BG$r_STR[which(encode_neuron == 33)] > 40000)


index_axon_IT_BG_33 <-
  index_axon_IT[which(encode_neuron == 33)][2480]

info_tmp[index_axon_IT_BG_33,]

nl_axon_IT_BG_33 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_IT_BG_33],
    neuronnames = shared_NeuronName[index_axon_IT_BG_33])

#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])
library(rgl)
open3d()
snp::plot_morphology(
  neuronList = nl_axon_IT_BG_33,
  hemisphere = info_tmp$Hemisphere[index_axon_IT_BG_33],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")

mesh_root <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")
shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "blue", alpha = 0.1)
shade3d(mesh_lr_STR, color = "red", alpha = 0.1)



rgl.viewpoint(fov = 0, userMatrix = userMatrix_3d_2)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_horizontal)



# IT neurons. iSTR, cSTR, cGPe  (35) ----

table(info_tmp$Region[index_axon_IT[which(encode_neuron == 35)]])


which(proj_IT_BG$l_GPe[which(encode_neuron == 35)] > 4000)


index_axon_IT_BG_35 <-
  index_axon_IT[which(encode_neuron == 35)][117]

info_tmp[index_axon_IT_BG_35,]

nl_axon_IT_BG_35 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_IT_BG_35],
    neuronnames = shared_NeuronName[index_axon_IT_BG_35])

#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])
library(rgl)
open3d()
snp::plot_morphology(
  neuronList = nl_axon_IT_BG_35,
  hemisphere = info_tmp$Hemisphere[index_axon_IT_BG_35],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")

mesh_root <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")
shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_lr_GPe, color = "blue", alpha = 0.1)
shade3d(mesh_lr_STR, color = "red", alpha = 0.1)



rgl.viewpoint(fov = 0, userMatrix = userMatrix_3d_2)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_horizontal)

