#' @export
get_body_area_normalizer <- function(measurement_data, head_width_data){
  measurement_data$species <- head_width_data$species[match(measurement_data$ID, head_width_data$chromatogram_ID)]
  measurement_data$head_width <- head_width_data$head_width[match(measurement_data$ID, head_width_data$chromatogram_ID)]
  species_coeffs <- list()
  for(species in c("F. sanguinea", "F. fusca")){
    species_coeffs[[species]] <- numeric(3)
    for(view in c("Thorax dorsal view", "Head dorsal view", "Thorax lateral view")){
      mask <- measurement_data$view==view & measurement_data$species==species
      species_coeffs[[species]] <- species_coeffs[[species]] +
        lm(area~poly(head_width,2, raw=TRUE), data=measurement_data[mask,])$coefficients
    }
  }
  standarized_area <- sum(species_coeffs[["F. sanguinea"]]*c(1,1150,1150^2))
  # fusca_standard_head_width <- optim(par=1.15, function(x) abs(standarized_area-sum(species_coeffs[["F. fusca"]]*c(1,x,x^2))),
  #                                    method="Brent", lower=1, upper=1.3)$par
  # print(sprintf("fusca head width equivalent of 1.15 mm: %f", fusca_standard_head_width))
  normalizer <- function(species, head_width){
    areas <- mapply(function(sp,hw) sum(species_coeffs[[sp]]*c(1,hw,hw^2)), species, head_width)
    return(standarized_area/areas)
  }
  return(normalizer)
}
