body_surface_area <- read.csv("surface_area.csv")
rownames(body_surface_area) <- body_surface_area$X
body_surface_area[body_surface_area$ID==389,"head_width"] = 946 # measured again
# sample 356 should be F. fusca but is F. sanguinea so it is mislabeled
body_surface_area <- body_surface_area[body_surface_area$ID!=356,]
body_surface_area <- body_surface_area[,colnames(body_surface_area)!="X"]
body_surface_area$view[body_surface_area$view=="top"] = "Thorax dorsal view"
body_surface_area$view[body_surface_area$view=="side"] = "Thorax lateral view"
body_surface_area$view[body_surface_area$view=="head"] = "Head dorsal view"
usethis::use_data(body_surface_area, overwrite = TRUE)


