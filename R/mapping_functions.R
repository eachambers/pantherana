import_range_maps <- function(path) {
  yava <- maptools::readShapePoly(paste0(path, "/yavapaiensis/data_0.shp"))
  yava <- fortify(yava)
    
  magn <- maptools::readShapePoly(paste0(path, "/magnaocularis_58656/species_58656.shp"))
  magn <- fortify(magn)
  
  forr <- maptools::readShapePoly(paste0(path, "/forreri/data_0.shp"))
  forr <- fortify(forr)
  
  spnov <- maptools::readShapePoly(paste0(path, "/species b (sp.4).shp"))
  spnov <- fortify(spnov)
  
  omil <- maptools::readShapePoly(paste0(path, "/new_omiltemana_wpapa/species_58687.shp"))
  omil <- fortify(omil)
  
  spec <- maptools::readShapePoly(paste0(path, "/spectabilis_58722/species_58722.shp"))
  spec <- fortify(spec)
  
  lenca <- maptools::readShapePoly(paste0(path, "/lenca_handmade/species_58653.shp"))
  lenca <- fortify(lenca)

  macro1 <- maptools::readShapePoly(paste0(path, "/redlist_species_data_63e03aaf-8763-4637-a167-0cd1cb2745e6/data_0.shp"))
  macro1 <- fortify(macro1)
  
  brown <- maptools::readShapePoly(paste0(path, "/redlist_species_data_bbde5bab-e44a-45d1-afd3-6dc21316d07b/data_0.shp"))
  brown <- fortify(brown)
  
  berl <- maptools::readShapePoly(paste0(path, "/redlist_species_data_f8832c88-565b-431b-a5d2-1be9e0f33f54/data_0.shp"))
  berl <- fortify(berl)
  
  neo <- maptools::readShapePoly(paste0(path, "/neovolcanica/data_0.shp"))
  neo <- fortify(neo)
  
  return(list(yava = yava, magn = magn, forr = forr, spnov = spnov, omil = omil, spec = spec, lenca = lenca))
}
