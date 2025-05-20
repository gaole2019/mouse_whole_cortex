# project mouse whole cortex
# figure 5. cortico-thalamic projections of PT neurons
# VISp

# 1. select neurons with extensive ptojections to selected thalamic nuclei ----
index_PT_VISp <-
  which(rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_LGd", "l_LP", "l_RT"))]) > 1000 &
        1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "PT"] &
        info_tmp$Region == "VISp")


Target_acronym <- paste0("l_", c("LGd", "LP", "RT", "PG"))

proj_PT_VISp <-
  shared_proj_Combined[index_PT_VISp,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_PT_VISp) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_PT_VISp]

colnames(proj_PT_VISp) <- Target_acronym



# 2. calculate their projection preference ----
ratio_LGd <-
  proj_PT_VISp[,"l_LGd"] / (proj_PT_VISp[,"l_LGd"] + proj_PT_VISp[,"l_LP"])

order_ratio_LGd <- order(ratio_LGd, decreasing = T)


VISp_label_pattern_proj <-
  matrix(nrow = length(index_PT_VISp), ncol = 1, data = "LGd and LP")

for (iNeuron in 1:length(index_PT_VISp)) {
  if (ratio_LGd[iNeuron] > 0.736)#quantile(ratio_LGd, 0.75))
    VISp_label_pattern_proj[iNeuron, 1] <- "LGd prefer"
  if (ratio_LGd[iNeuron] < 0.219)#quantile(ratio_LGd, 0.25))
    VISp_label_pattern_proj[iNeuron, 1] <- "LP prefer"
}


# 3. plot heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_PT_VISp),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = VISp_label_pattern_proj,

      Ratio = anno_points(x = ratio_LGd),


      col = list(AxonProj = c("LGd prefer" = "red",
                              "LGd and LP" = "green",
                              "LP prefer" = "blue"))
    ),
  name = "Axon length (um)",
  heatmap_legend_param	=
    list(title_gp = gpar(fontsize = 8),
         labels_gp = gpar(fontsize = 8),
         title_position = "leftcenter-rot"),

  cluster_rows = F,
  row_names_side = "left",
  row_names_rot = 0,
  row_names_gp = gpar(fontsize = 8),

  column_order = order_ratio_LGd,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)


# 4. example PT neurons preferentially projecting to LGd ----
mesh_LGd <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/170_LGd.obj")

mesh_LP <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/218_LP.obj")


mesh_SCm <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/294_SCm.obj")

mesh_PG <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/931_PG.obj")


mesh_l_LGd <- getMeshLeft(mesh_LGd)
mesh_l_LP <- getMeshLeft(mesh_LP)


index_demo_PT_LGd <-
  index_PT_VISp[which(VISp_label_pattern_proj == "LGd prefer" &
                      proj_PT_VISp$l_LGd > 10000)]

demo_PT_LGd <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/233666/swc_allen_axon/169.swc")

open3d()
plot3d(demo_PT_LGd, color = "blue", WithNodes = F, soma = 100)
shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_LP, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_SCm, color = "orange", alpha = 0.1)


# LGd and LP ----

index_demo_PT_LGd_LP <-
  index_PT_VISp[which(VISp_label_pattern_proj == "LGd and LP" &
                        proj_PT_VISp$l_LGd > 5000 &
                        proj_PT_VISp$l_LP > 5000)]

demo_PT_LGd_LP <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/233666/swc_allen_axon/126.swc")

open3d()
plot3d(demo_PT_LGd_LP, color = "blue", WithNodes = F, soma = 100)
shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_LP, color = "red", alpha = 0.1)

shade3d(mesh_SCm, color = "orange", alpha = 0.1)

shade3d(mesh_root, color = "white", alpha = 0.1)

# LP dominant ----
index_demo_PT_LP <-
  index_PT_VISp[which(VISp_label_pattern_proj == "LP prefer" &
                        proj_PT_VISp$l_LP > 5000 &
                        proj_PT_VISp$l_PG > 0)]

demo_PT_LP <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/221226/swc_allen_axon/015.swc")
demo_PT_LP$d$Z <- 11400 - demo_PT_LP$d$Z

open3d()
plot3d(demo_PT_LP, color = "blue", WithNodes = F, soma = 100)
shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_LP, color = "red", alpha = 0.1)

shade3d(mesh_PG, color = "orange", alpha = 0.1)

shade3d(mesh_root, color = "white", alpha = 0.1)
