# project mouse whole cortex
# figure 5. cortico-thalamic projections
# RT projecting PT and CT neurons


# 1. flatmap. soma distribution of RT-projecting PT neurons. ----
shared_soma_XYZUV <-
  read.csv("/run/media/admin/LeGAO/Project_Cortex/shared_files_20240402/soma_XYZUV.csv")

U <- soma_UV$U[index_RT_PT]
V <- soma_UV$V[index_RT_PT]
Color <- "red"
Radius <- 6
filename = "/tmp/iRT_PT.svg"


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
    c(
      svg_element,
      sprintf(
        "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke: none; fill: %s\"/>",
        U[iElement],
        V[iElement],
        Radius[iElement],
        Color[iElement]
      )
    )
}

out <- c(svg_PFC[1:50], svg_element, svg_PFC[51:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)



# 2. flatmap. soma distribution of RT-projecting CT neurons. ----
U <- shared_soma_XYZUV$U[index_RT_CT]
V <- shared_soma_XYZUV$V[index_RT_CT]
Color <- "red"
Radius <- 6
filename = "/tmp/iRT_CT.svg"


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
    c(
      svg_element,
      sprintf(
        "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke: none; fill: %s\"/>",
        U[iElement],
        V[iElement],
        Radius[iElement],
        Color[iElement]
      )
    )
}

out <- c(svg_PFC[1:50], svg_element, svg_PFC[51:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)
