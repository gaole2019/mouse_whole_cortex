# project mouse whole cortex
# figure 4 cortico-basal ganglia pathway
# PT neurons


index_axon_PT <-
  df_All_index_AxonType$Index[df_All_index_AxonType$AxonType %in% c("PT")]

# 1. set threshold for target regions ----
threshold_region <- c(1000, 1000, 100, 100, 1000,
                      1000, 1000, 100, 100, 1000)
names(threshold_region) <-
  c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
    "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")


encode_neuron <- matrix(nrow = nrow(proj_PT_BG_new), ncol = 1)
for (iNeuron in 1:nrow(proj_PT_BG_new)) {
  encode <- 0
  for (iTarget in 1:ncol(proj_PT_BG_new)) {
    encode <- encode +
      ifelse(proj_PT_BG_new[iNeuron, iTarget] > threshold_region[iTarget], 1, 0) * 2^(ncol(proj_PT_BG_new) - iTarget)
  }
  encode_neuron[iNeuron, 1] <- encode
}

aaaa <- sort(table(encode_neuron), decreasing = T)

barplot(aaaa)

#
# vector_names <- matrix(nrow = length(aaaa), ncol = 1)
# for (iPatttern in 1:length(aaaa)) {
#   binary_pattern <- as.integer(intToBits(as.numeric(names(aaaa)[iPatttern])))[1:10]
#   tmp_name <- ""
#   for (iTarget in 1:ncol(proj_PT_BG_new)) {
#     if (binary_pattern[iTarget] == 1) {
#       tmp_name <- paste0(tmp_name, "_", colnames(proj_PT_BG_new)[iTarget])
#     }
#   }
#
#   vector_names[iPatttern, 1] <- tmp_name
# }



# 2. plot combination matrix ----
bbbb <- proj_PT_BG_new
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
m = m[comb_size(m) > 19 & comb_degree(m) > 0]
UpSet(m = m,
      comb_order = order(comb_size(m), decreasing = T),
      pt_size = unit(4, "mm"))



motif_number <- sort(comb_size(m), decreasing = T)

motif_number_code <- strtoi(x = names(motif_number), 2)


df_BG_PT <-
  data.frame(
    Index = index_axon_PT,
    Layer = factor(info_tmp$Layer[index_axon_PT]),
    Region = factor(info_tmp$Region[index_axon_PT]),
    AxonType = factor(info_tmp$AxonType[index_axon_PT]),
    AxonEncode = factor(encode_neuron, levels = motif_number_code),
    XYZ_UV_CC[index_axon_PT,])
head(df_BG_PT)



library(dplyr)
library(tidyverse)
library(tidyr)

ttt <-
  df_BG_PT %>%
  group_by(Region, AxonEncode) %>%
  summarize(Count = n())


library(ggplot2)

ttt$AxonEncode <-
  factor(x = ttt$AxonEncode, levels = motif_number_code)


# 4. laminar distribution of neurons with each projection pattern ----
ggplot(data = ttt,
       mapping = aes(x = AxonEncode, y = Region)) +
  geom_point(mapping = aes(size = Count))


  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("L1" = "red",
                               "L2/3" = "red",
                               "L4" = "green",
                               "L5" = "blue",
                               "L6a" = "orange",
                               "L6b" = "orange"))

# 3. example neurons ----
## STN projecting neurons (pattern #9) ----

which(proj_PT_BG_new$l_STN > 1000 &
      encode_neuron == 640 &
        info_tmp$ProjectID[index_axon_PT] == "ISO")

library(nat)
nl_neurite_PT_BG_640 <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT[c(1687,1781)]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT[c(1687,1781)]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT[c(1687,1781)]]

open3d()
snp::plot_morphology(neuronList = nl_neurite_PT_BG_640[1],
                     hemisphere = hemisphere[1], color = "blue")

plot3d(x = nl_neurite_PT_BG_640)


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "orange", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_GPi, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_horizontal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)


## GPe projecting neurons (pattern #10) ----
which(proj_PT_BG_new$l_GPe > 2000 &
        encode_neuron == 864 &
        info_tmp$ProjectID[index_axon_PT] == "ISO")

library(nat)
nl_neurite_PT_BG_864 <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT[c(1547,1573)]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT[c(1547,1573)]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT[c(1547,1573)]]

open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_BG_864[1],
  hemisphere = hemisphere[1], color = "blue")


## STN and GPe co-projecting neurons (pattern #5) ----
which(proj_PT_BG_new$l_GPe > 2000 &
        proj_PT_BG_new$l_STN > 1000 &
        encode_neuron == 992 &
        info_tmp$ProjectID[index_axon_PT] == "ISO")

library(nat)
nl_neurite_PT_BG_992 <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT[c(1598,1639)]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT[c(1598,1639)]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT[c(1598,1639)]]

open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_BG_992[1],
  hemisphere = hemisphere[1], color = "blue")


# 4. new type of direct pathways ----
## GPi projecting neurons (pattern #3) ----
which(proj_PT_BG_new$l_GPi > 1000 &
        encode_neuron == 576 &
        info_tmp$ProjectID[index_axon_PT] == "ISO")

library(nat)
nl_neurite_PT_BG_576 <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT[c(1644,1675)]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT[c(1644,1675)]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT[c(1644,1675)]]

open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_BG_576,
  hemisphere = hemisphere,
  color = "blue")


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPi, color = "red", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)




## GPi and SNr projecting neurons (pattern #4) ----
which(proj_PT_BG_new$l_GPi > 1000 &
        proj_PT_BG_new$l_SNr > 2000 &
        encode_neuron == 736 &
        info_tmp$ProjectID[index_axon_PT] == "ISO")

vector_index_736 <- c(1908,1909)

library(nat)
nl_neurite_PT_BG_736 <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT[vector_index_736]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT[vector_index_736]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT[vector_index_736]]


open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_BG_736[2],
  hemisphere = hemisphere[2],
  color = "blue")

shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPi, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)



# backup -----

proj_PT_BG <- shared_proj_Combined[
  index_axon_PT,
  snp::getStructureFromAcronym(c("l_STR", "l_GPe", "l_GPi", "l_STN", "l_SNr"))]

colnames(proj_PT_BG) <- c("l_STR", "l_GPe", "l_GPi", "l_STN", "l_SNr")


# 1. select neurons projecting to major targets. note some neurons do not project to major targets ----
index_PT_minor <- which(rowSums(proj_PT_BG) < 0)

if(length(index_PT_minor) == 0) {
  index_axon_PT_BG <- index_axon_PT
} else {
  index_axon_PT_BG <- index_axon_PT[-index_PT_minor]
}


proj_PT_BG <- shared_proj_Combined[
  index_axon_PT_BG,
  snp::getStructureFromAcronym(c("l_STR", "l_GPe", "l_GPi", "l_STN", "l_SNr"))]

colnames(proj_PT_BG) <- c("l_STR", "l_GPe", "l_GPi", "l_STN", "l_SNr")




# 3. assign neuron projections as different patterns ----

projection_pattern <- matrix(nrow = nrow(proj_PT_BG), ncol = 1)
for (iNeuron in 1:nrow(proj_PT_BG)) {
  bit1 <- ifelse(proj_PT_BG$l_SNr[iNeuron] > 1000, 1, 0)
  bit2 <- ifelse(proj_PT_BG$l_GPi[iNeuron] > 500, 1, 0)
  bit3 <- ifelse(proj_PT_BG$l_STN[iNeuron] > 100, 1, 0)
  bit4 <- ifelse(proj_PT_BG$l_GPe[iNeuron] > 1000, 1, 0)
  bit5 <- ifelse(proj_PT_BG$l_STR[iNeuron] > 1000, 1, 0)

  projection_pattern[iNeuron, 1] <- bit5*16 + bit4*8 + bit3*4 + bit2*2 + bit1
}


p_STR <- length(which(proj_PT_BG$l_STR > 1000)) / nrow(proj_PT_BG)
p_GPe <- length(which(proj_PT_BG$l_GPe > 1000)) / nrow(proj_PT_BG)
p_STN <- length(which(proj_PT_BG$l_STN > 100)) / nrow(proj_PT_BG)
p_SNr <- length(which(proj_PT_BG$l_SNr > 1000)) / nrow(proj_PT_BG)

p_vector <- c(p_STR, p_GPe, p_STN, p_SNr)

p_motif_pattern <- matrix(nrow = 16, ncol = 1)

counter <- 0
for(iPattern in 0:15) {
  tmpFlag <- rev(as.integer(intToBits(iPattern))[1:4])

  p_motif <- ifelse(tmpFlag[1] == 1, p_vector[1], 1-p_vector[1]) *
    ifelse(tmpFlag[2] == 1, p_vector[2], 1-p_vector[2]) *
    ifelse(tmpFlag[3] == 1, p_vector[3], 1-p_vector[3]) *
    ifelse(tmpFlag[4] == 1, p_vector[4], 1-p_vector[4])

  counter <- counter + 1
  p_motif_pattern[counter, 1] <- p_motif
}

p_STR*(1-p_GPe)*(1-p_STN)*(1-p_BG) * 3366





barplot(sort(table(projection_pattern), decreasing = T))


sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 8)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 10)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 9)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 11)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 15)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 14)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 13)]]))
sort(table(info_tmp$Region[index_axon_PT_BG[which(projection_pattern == 12)]]))




mesh_root <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")

library(rgl)
mesh_lr_GPe <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1022_GPe.obj")
mesh_lr_STN <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/470_STN.obj")
mesh_lr_GPi <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1031_GPi.obj")
mesh_lr_SNr <- rgl::readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/381_SNr.obj")

mesh_l_GPe <- getMeshLeft(mesh_lr_GPe)
mesh_l_STN <- getMeshLeft(mesh_lr_STN)
mesh_l_GPi <- getMeshLeft(mesh_lr_GPi)
mesh_l_SNr <- getMeshLeft(mesh_lr_SNr)





# PT neurons project to four targets (22) ----

which(proj_PT_BG$l_GPe[which(projection_pattern == 22)] < 100 &
        proj_PT_BG$l_STN[which(projection_pattern == 22)] > 1000 &
        proj_PT_BG$l_GPi[which(projection_pattern == 22)] > 500 &
        proj_PT_BG$l_SNr[which(projection_pattern == 22)] < 100)
#[1] 101 130 149 275 300 349


index_axon_PT_BG_22 <-
  index_axon_PT_BG[which(projection_pattern == 22)][229]

info_tmp[index_axon_PT_BG_22,]


nl_axon_PT_BG_22 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_PT_BG_22],
    neuronnames = shared_NeuronName[index_axon_PT_BG_22])
#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_PT_BG_22,
  hemisphere = info_tmp$Hemisphere[index_axon_PT_BG_22],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "red", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_GPi, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)




# PT neurons project to four targets (10) ----

which(proj_PT_BG$l_GPe[which(projection_pattern == 10)] < 100 &
        proj_PT_BG$l_STN[which(projection_pattern == 10)] > 1000 &
        proj_PT_BG$l_SNr[which(projection_pattern == 10)] < 100)
# [1]  97 131 163 171 185 281 351 352 356 383 456 465 533 539

index_axon_PT_BG_10 <-
  index_axon_PT_BG[which(projection_pattern == 10)][539]


nl_axon_PT_BG_10 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_PT_BG_10],
    neuronnames = shared_NeuronName[index_axon_PT_BG_10])
#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_PT_BG_10,
  hemisphere = info_tmp$Hemisphere[index_axon_PT_BG_10],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "red", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)




# PT neurons project to four targets (13) ----

which(proj_PT_BG$l_STN[which(projection_pattern == 13)] < 100 &
        proj_PT_BG$l_GPe[which(projection_pattern == 13)] > 5000 &
        proj_PT_BG$l_SNr[which(projection_pattern == 13)] > 2000)


index_axon_PT_BG_13 <-
  index_axon_PT_BG[which(projection_pattern == 13)][303]


nl_axon_PT_BG_13 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_PT_BG_13],
    neuronnames = shared_NeuronName[index_axon_PT_BG_13])
#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_PT_BG_13,
  hemisphere = info_tmp$Hemisphere[index_axon_PT_BG_13],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "red", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)





# PT neurons project (9) ----


which(proj_PT_BG$l_GPe[which(projection_pattern == 9)] < 100 &
        proj_PT_BG$l_STN[which(projection_pattern == 9)] < 100 &
        proj_PT_BG$l_SNr[which(projection_pattern == 9)] > 5000)


index_axon_PT_BG_9 <-
  index_axon_PT_BG[which(projection_pattern == 9)][248]

info_tmp[index_axon_PT_BG_9,]


nl_axon_PT_BG_9 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_PT_BG_9],
    neuronnames = shared_NeuronName[index_axon_PT_BG_9])
#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_PT_BG_9,
  hemisphere = info_tmp$Hemisphere[index_axon_PT_BG_9],
  flip = TRUE, soma_size = 100, color = "blue")
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")


#shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "red", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)




# PT neurons project (12) ----


which(proj_PT_BG$l_GPe[which(projection_pattern == 12)] > 2000 &
        proj_PT_BG$l_STN[which(projection_pattern == 12)] < 100 &
        proj_PT_BG$l_SNr[which(projection_pattern == 12)] < 100)

#[1]  64 155 202 209 215 222 230 237 259 268 277 285 286


index_axon_PT_BG_12 <-
  index_axon_PT_BG[which(projection_pattern == 12)][237]

info_tmp[index_axon_PT_BG_12,]

nl_axon_PT_BG_12 <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_axon_PT_BG_12],
    neuronnames = shared_NeuronName[index_axon_PT_BG_12])
#
# nl_dendrite_PT_BG_15 <-
#   nat::read.neurons(
#     paths = shared_SWCPath_dendrite_allen[index_axon_PT_BG_15],
#     neuronnames = shared_NeuronName[index_axon_PT_BG_15])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_PT_BG_12,
  hemisphere = info_tmp$Hemisphere[index_axon_PT_BG_12],
  flip = TRUE,
  color = "blue",
  soma_size = 100)
#
# snp::plot_morphology(
#   neuronList = nl_dendrite_PT_SSbfd_VPM,
#   hemisphere = info_tmp$Hemisphere[index_full_PT_SSpbfd_VPM],
#   flip = TRUE, soma_size = 100, color = "red")


#shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_GPe, color = "red", alpha = 0.1)
shade3d(mesh_l_STN, color = "red", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.1)


rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)
rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)


# PT neurons projecting to SNr ----

index_axon_PT_SNr <- index_axon_PT_BG[which(projection_pattern %in% c(9, 11, 15, 13))]



df_region <- data.frame(Region = info_tmp$Region[index_axon_PT],
                        FlagSNr = projection_pattern)












proj_PT_normalized <-
  apply(proj_PT_BG, 2, function(x) x / max(x)) # normalized to the maximal value


proj_PT_region <-
  aggregate(x = proj_PT, by = list(Region = info_tmp$Region[index_axon_PT]), FUN = median)

rownames(proj_PT_region) <- proj_PT_region$Region

proj_PT_region_normalized <-
  apply(as.matrix(proj_PT_region[,-1]), 2, function(x) x/max(x)) # normalized to the maximal value

rownames(proj_PT_region_normalized) <- proj_PT_region$Region
colnames(proj_PT_region_normalized) <- proj_PT_target_acronym


library(ComplexHeatmap)
ht <- Heatmap(matrix = proj_PT_normalized, cluster_rows = T, cluster_columns = F)
draw(ht)


library(ComplexHeatmap)
Heatmap(
  #matrix = proj_L5PT_region_density,
  matrix = proj_PT_normalized,
  #matrix = proj_L5PT_region_normalized,
  name = "Normalized\naxon length",
  right_annotation =
    rowAnnotation(Ratio_contra = anno_barplot(x = region_ratio_contra, axis_param = c(side = "top")),
                  Mean_target = anno_barplot(df_summarize_numOfTarget$mean_target, axis_param = c(side = "top"))),

  column_split = substring(proj_PT_target_acronym, first = 1, last = 2),

  #
  # column_split =
  #   res_region$Refine_coarse[match(info_tmp$Region[index_axon_L5PT],
  #                                  res_region$Region_fine)],

  cluster_column_slices = F,
  show_column_names = T,
  column_names_side = "top",


  row_names_side = "left",

  cluster_columns = F,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",

  clustering_distance_rows = "spearman",
  clustering_method_rows = "ward.D2",
  cluster_rows = T)





