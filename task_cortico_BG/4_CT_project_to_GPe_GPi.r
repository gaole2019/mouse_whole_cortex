# project mouse whole cortex
# figure 4. cortico-basal ganalia projections
# CT neurons project to GPe/i


# 1. example CT neuron projecting to GPe/i ----
index_example_CT <-
  which(rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_GPe", "l_GPi"))]) > 2000 &
          1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("CT")] &
          info_tmp$Region == "AId")


table(info_tmp$Region[index_example_CT])


nl_neuron_CT_1 <-
  nat::read.neuron(f = SWCPath_neurite[index_example_CT[1]])
nl_neuron_CT_1$d$Z <-
  ifelse(test = nl_neuron_CT_1$d$Z[1]>5700,
         yes = 11400 - nl_neuron_CT_1$d$Z,
         no = nl_neuron_CT_1$d$Z)

nl_neuron_CT_2 <-
  nat::read.neuron(f = SWCPath_neurite[index_example_CT[2]])
nl_neuron_CT_2$d$Z <-
  ifelse(test = nl_neuron_CT_2$d$Z[1]>5700,
         yes = 11400 - nl_neuron_CT_2$d$Z,
         no = nl_neuron_CT_2$d$Z)

nl_neuron_CT_3 <-
  nat::read.neuron(f = SWCPath_neurite[index_example_CT[3]])
nl_neuron_CT_3$d$Z <-
  ifelse(test = nl_neuron_CT_3$d$Z[1]>5700,
         yes = 11400 - nl_neuron_CT_3$d$Z,
         no = nl_neuron_CT_3$d$Z)



mesh_GPe <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/1022_GPe.obj")
mesh_GPi <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/1031_GPi.obj")

mesh_AId <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/104_AId.obj")

mesh_l_GPe <- getMeshLeft(mesh = mesh_GPe)
mesh_l_GPi <- getMeshLeft(mesh = mesh_GPi)

mesh_l_AId <- getMeshLeft(mesh = mesh_AId)



open3d()
plot3d(x = nl_neuron_CT_1, color = "red", WithNodes = F, soma = 200)
plot3d(x = nl_neuron_CT_2, color = "green", WithNodes = F, soma = 200)
plot3d(x = nl_neuron_CT_3, color = "blue", WithNodes = F, soma = 200)

shade3d(mesh_l_GPe, color = "orange", alpha = 0.1)
shade3d(mesh_l_GPi, color = "orange", alpha = 0.1)

shade3d(mesh_l_AId, color = "white", alpha = 0.1)


shade3d(mesh_root, color = "white", alpha = 0.1)
