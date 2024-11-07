import_range_maps <- function(path) {
  # foothills species
  yava <- readOGR(paste0(path, "/yavapaiensis/data_0.shp"))
  magn <- readOGR(paste0(path, "/magnaocularis/data_0.shp"))
  omil <- readOGR(paste0(path, "/new_omiltemana_wpapa/species_58687.shp"))
  
  tls <- read_tsv(here("data", "type_localities_rec.txt"))
  aten_short <- tls %>% 
    filter(Species == "Aten_short_other" | Species == "Atenquique_short")
  aten_long <- tls %>% 
    filter(Species == "Atenquique_long")
  
  # CENTAM species
  lenca <- readOGR(paste0(path, "/lenca_handmade/species_58653.shp"))
  spnov <- readOGR(paste0(path, "/spnov/species b (sp.4).shp"))
  
  # forreri
  forr <- readOGR(paste0(path, "/forreri/data_0.shp"))
  
  # ATL_MXPL
  berl <- readOGR(paste0(path, "/berneo/species_58561.shp"))
  spec <- readOGR(paste0(path, "/spectabilis/species_58722.shp"))
  macro1 <- readOGR(paste0(path, "/new_macroglossa/data_0.shp"))
  macro2 <- readOGR(paste0(path, "/taylori_58732/species_58732.shp"))
  
  return(list(yava = yava, magn = magn, omil = omil, lenca = lenca, spnov = spnov, 
              forr = forr, berl = berl, spec = spec, macro1 = macro1, macro2 = macro2,
              aten_long = aten_long, aten_short = aten_short))
}

mapping_colors <- function(path_to_tls) {
  tls <- read_tsv(path_to_tls)
  
  # foothills
  yava <- tls %>% filter(Species == "yavapaiensis") %>% pull(Color)
  magn <- tls %>% filter(Species == "magnaocularis") %>% pull(Color)
  omil <- tls %>% filter(Species == "omiltemana") %>% pull(Color)
  aten_long <- tls %>% filter(Species == "Atenquique_long") %>% pull(Color)
  aten_short <- tls %>% filter(Species == "Atenquique_short") %>% pull(Color)
  # CENTAM
  lenca <- tls %>% filter(Species == "lenca") %>% pull(Color)
  spnov <- "#82ccc8"
  # forreri
  forr <- tls %>% filter(Species == "forreri") %>% pull(Color)
  # ATL_MXPL
  berl <- tls %>% filter(Species == "berlandieri") %>% pull(Color)
  spec <- tls %>% filter(Species == "spectabilis") %>% pull(Color)
  macro <- tls %>% filter(Species == "macroglossa") %>% pull(Color)
  
  return(list(yava = yava, magn = magn, omil = omil, lenca = lenca, spnov = spnov, 
              forr = forr, berl = berl, spec = spec, macro = macro,
              aten_long = aten_long, aten_short = aten_short))
}
