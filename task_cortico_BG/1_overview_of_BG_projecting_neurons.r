# project mouse whole cortex
# figure 4. cortico-basal ganglia projection
# overview of neurons projecting to BG

df_All_index_AxonType <-
  read.csv(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_fntdist/df_All_index_AxonType.csv", row.names = 1)

# Note that CP and ACB were combined as STR in the following analysis.
target_acronym <-
  c("CP", "ACB", "GPe", "STN", "GPi", "SNr")

target_acronym <-
  c(paste0("l_", target_acronym), paste0("r_", target_acronym))


# IT
proj_IT_BG <-
  shared_proj_Combined[df_All_index_AxonType$Index[df_All_index_AxonType$AxonType %in% c("L23IT", "L4IT", "L5IT", "L6IT")],
  #index_axon_IT,
  snp::getStructureFromAcronym(target_acronym)]

colnames(proj_IT_BG) <- target_acronym


proj_IT_BG_new <-
  data.frame(l_STR = proj_IT_BG$l_CP + proj_IT_BG$l_ACB,
             l_GPe = proj_IT_BG$l_GPe,
             l_STN = proj_IT_BG$l_STN,
             l_GPi = proj_IT_BG$l_GPi,
             l_SNr = proj_IT_BG$l_SNr,

             r_STR = proj_IT_BG$r_CP + proj_IT_BG$r_ACB,
             r_GPe = proj_IT_BG$r_GPe,
             r_STN = proj_IT_BG$r_STN,
             r_GPi = proj_IT_BG$r_GPi,
             r_SNr = proj_IT_BG$r_SNr)

colnames(proj_IT_BG_new) <-
  c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
    "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")



# PT
proj_PT_BG <- shared_proj_Combined[
  df_All_index_AxonType$Index[df_All_index_AxonType$AxonType %in% c("PT")],
  snp::getStructureFromAcronym(target_acronym)]

colnames(proj_PT_BG) <- target_acronym


proj_PT_BG_new <-
  data.frame(l_STR = proj_PT_BG$l_CP + proj_PT_BG$l_ACB,
             l_GPe = proj_PT_BG$l_GPe,
             l_STN = proj_PT_BG$l_STN,
             l_GPi = proj_PT_BG$l_GPi,
             l_SNr = proj_PT_BG$l_SNr,

             r_STR = proj_PT_BG$r_CP + proj_PT_BG$r_ACB,
             r_GPe = proj_PT_BG$r_GPe,
             r_STN = proj_PT_BG$r_STN,
             r_GPi = proj_PT_BG$r_GPi,
             r_SNr = proj_PT_BG$r_SNr)

colnames(proj_PT_BG_new) <-
  c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
    "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")



# CT
proj_CT_BG <- shared_proj_Combined[
  index_axon_CT,
  snp::getStructureFromAcronym(target_acronym)]

colnames(proj_CT_BG) <- target_acronym



proj_CT_BG_new <-
  data.frame(l_STR = proj_CT_BG$l_CP + proj_CT_BG$l_ACB,
             l_GPe = proj_CT_BG$l_GPe,
             l_STN = proj_CT_BG$l_STN,
             l_GPi = proj_CT_BG$l_GPi,
             l_SNr = proj_CT_BG$l_SNr,

             r_STR = proj_CT_BG$r_CP + proj_CT_BG$r_ACB,
             r_GPe = proj_CT_BG$r_GPe,
             r_STN = proj_CT_BG$r_STN,
             r_GPi = proj_CT_BG$r_GPi,
             r_SNr = proj_CT_BG$r_SNr)

colnames(proj_CT_BG_new) <-
  c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
    "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")


###
rowSum_IT_BG <- rowSums(proj_IT_BG_new)
rowSum_PT_BG <- rowSums(proj_PT_BG_new)
rowSum_CT_BG <- rowSums(proj_CT_BG_new)

df_neuronClass_BG <-
  data.frame(NeuronClass = factor(c(rep("IT", length(index_axon_IT)),
                             rep("PT", length(index_axon_PT)),
                             rep("CT", length(index_axon_CT))),levels = c("IT", "PT", "CT")),
             Flag_BG = factor(c(ifelse(rowSum_IT_BG > 5000, 1, 0),
                         ifelse(rowSum_PT_BG > 5000, 1, 0),
                         ifelse(rowSum_CT_BG > 5000, 1, 0))))

# figure. ----
library(ggplot2)
ggplot(data = df_neuronClass_BG, aes(x = NeuronClass,fill = Flag_BG)) +
  geom_bar(position = "dodge")



# threshold of axon lengths to determine whether or not project to a specific target
# "l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr"
threshold_region <- c(1000, 1000, 100, 100, 1000,
                      1000, 1000, 100, 100, 1000)
names(threshold_region) <-
  c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
    "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")


# IT neurons ----
proportion_IT <- matrix(nrow = ncol(proj_IT_BG_new), ncol = 1)
rownames(proportion_IT) <- colnames(proj_IT_BG_new)

for (iTarget in 1:ncol(proj_IT_BG_new)) {
  proportion_IT[iTarget, 1] <-
    length(which(proj_IT_BG_new[,iTarget] > threshold_region[iTarget])) / nrow(proj_IT_BG_new)
}



# PT neurons ----
proportion_PT <- matrix(nrow = ncol(proj_PT_BG_new), ncol = 1)
rownames(proportion_PT) <- colnames(proj_PT_BG_new)

for (iTarget in 1:ncol(proj_PT_BG_new)) {
  proportion_PT[iTarget, 1] <-
    length(which(proj_PT_BG_new[,iTarget] > threshold_region[iTarget])) / nrow(proj_PT_BG_new)
}


# CT neurons ----
proportion_CT <- matrix(nrow = ncol(proj_CT_BG_new), ncol = 1)
rownames(proportion_CT) <- colnames(proj_CT_BG_new)

for (iTarget in 1:ncol(proj_CT_BG_new)) {
  proportion_CT[iTarget, 1] <-
    length(which(proj_CT_BG_new[,iTarget] > threshold_region[iTarget])) / nrow(proj_CT_BG_new)
}



df_neuronClass_single <-
  data.frame(NeuronClass = factor(c(rep("IT", 10),
                                    rep("PT", 10),
                                    rep("CT", 10)),
                                  levels = c("IT", "PT", "CT")),
             Target = factor(x = rep(c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
                                       "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr"), 3),
                             levels = c("l_STR", "l_GPe", "l_STN", "l_GPi", "l_SNr",
                                        "r_STR", "r_GPe", "r_STN", "r_GPi", "r_SNr")),
             Proportion = c(proportion_IT,
                            proportion_PT,
                            proportion_CT))

# figure ----
ggplot(data = df_neuronClass_single,
       mapping = aes(x = Target, y = Proportion, fill = NeuronClass)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_classic()






# 3. figure. flatmap showing soma distribution of neurons projecting to BG ----
soma_UV <-
  read.csv("/run/media/admin/LeGAO/Project_Cortex/shared_files_20240402/soma_XYZUV.csv")

# IT
U <- soma_UV$U[index_axon_IT]
V <- soma_UV$V[index_axon_IT]

Color <- rep("white", length(index_axon_IT))
Color[rowSum_IT_BG > 5000] <- "blue"
Radius <- 6
filename = "/tmp/IT_BG.svg"

# PT
U <- soma_UV$U[index_axon_PT]
V <- soma_UV$V[index_axon_PT]

Color <- rep("white", length(index_axon_PT))
Color[rowSum_PT_BG > 5000] <- "blue"
Radius <- 6
filename = "/tmp/PT_BG.svg"

# CT
U <- soma_UV$U[index_axon_CT]
V <- soma_UV$V[index_axon_CT]

Color <- rep("white", length(index_axon_CT))
Color[rowSum_CT_BG > 5000] <- "blue"
Radius <- 6
filename = "/tmp/CT_BG.svg"



conn <- file("/home/admin/shared/flatmap_Isocortex_template.svg")
svg_PFC <- readLines(conn)
close(conn)
fileConn <- file(filename, "w")
NumElement <- length(U)
svg_element <- c()
if (length(Color) == 1) {
  Color <- rep(Color, NumElement)
}
if (length(Radius) == 1) {
  Radius <- rep(Radius, NumElement)
}
for (iElement in 1:NumElement) {
  svg_element <-
    c(svg_element,
      sprintf("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke: none; fill: %s\"/>",
              U[iElement], V[iElement], Radius[iElement], Color[iElement]))
}

out <- c(svg_PFC[1:50], svg_element, svg_PFC[51:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)


