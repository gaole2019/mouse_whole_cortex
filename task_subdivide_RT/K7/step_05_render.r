# project mouse whole cortex
# subdivide RT
# render segmented structures


library(rgl)
library(stringr)


anno_structure_label <- read.table(
  file = "/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/K7/itksnap_label.txt", 
  comment.char = "#", 
  col.names = c('IDX', 'R', 'G', 'B', 'A', 'VIS', 'MSH', 'LABEL'))
anno_structure_label <- anno_structure_label[-1,]

anno_structure_label$ColorHex <-
  rgb(red = anno_structure_label$R, 
      green = anno_structure_label$G, 
      blue = anno_structure_label$B, 
      maxColorValue = 255)


NumStructure <- nrow(anno_structure_label)

mesh_list <- list()
for(iStructure in 1:NumStructure) {
  tmp_mesh <- 
    readOBJ(con = 
              paste0("/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/K7/structures_mesh/", 
                     anno_structure_label$IDX[iStructure], 
                     ".tiff_smoothed.obj"))
  #tmp_mesh$vb[1,] <- tmp_mesh$vb[1,] * 0.7
  #tmp_mesh$vb[2,] <- tmp_mesh$vb[2,] * 0.7
  #tmp_mesh$vb[3,] <- tmp_mesh$vb[3,] * 1.0
  mesh_list[[iStructure]] <- tmp_mesh
}


mesh_l_CP <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh_left_hemisphere/485_STRd.obj")

mesh_RT <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/262_RT.obj")

mesh_root <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")

#
open3d()
for(iStructure in 1:NumStructure) {
  shade3d(mesh_list[[iStructure]], 
          alpha = 1, 
          color = anno_structure_label$ColorHex[iStructure], 
          specular = "black")
}

rgl::shade3d(mesh_root, color = "white", alpha = 0.1)
rgl::shade3d(mesh_RT, color = "blue", alpha = 0.1)


# plot one by one 
for(iStructure in 1:NumStructure) {
  open3d()
  par3d("windowRect" = c(970, 109, 1500, 625))
  
  shade3d(mesh_list[[iStructure]], 
          alpha = 1, 
          color = anno_structure_label$ColorHex[iStructure], 
          specular = "black")
  
  rgl::shade3d(mesh_l_CP, color = "white", alpha = 0.1)
  rgl::shade3d(mesh_l_ACB, color = "white", alpha = 0.1)
  
  #rgl.viewpoint(fov = 0, userMatrix = userMatrix_3d_2)
  rgl.viewpoint(fov = 0, userMatrix = userMatrix_LM, zoom = 0.62)
  
  rgl.snapshot(
    paste0("/run/media/admin/LeGAO/Project_Cortex/task_subdivide_striatum/rendering/results/", 
           stringr::str_pad(string = iStructure, width = 2, side = "left", pad = "0"), 
           "_sagittal.png"))
  
  rgl.viewpoint(fov = 0, userMatrix = userMatrix_AP, zoom = 0.62)
  rgl.snapshot(
    paste0("/run/media/admin/LeGAO/Project_Cortex/task_subdivide_striatum/rendering/results/", 
           stringr::str_pad(string = iStructure, width = 2, side = "left", pad = "0"), 
           "_coronal.png"))
  
  rgl.close()
}

view3d(userMatrix = userMatrix_horizontal, fov = 0)
view3d(userMatrix = userMatrix_coronal, fov = 0)
view3d(userMatrix = userMatrix_sagittal, fov = 0)
view3d(userMatrix = userMatrix_3d_2, fov = 0)


rgl::shade3d(mesh_root, color = "white", alpha = 0.1)



