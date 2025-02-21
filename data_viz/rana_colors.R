#' Retrieve colors for plotting Rana analyses
#'
#' @param dataset_name specify dataset name ("ATL_MXPL", "PACMX", "CENTAM", or "forreri"); if NULL 
#' @param K_value K value to retrieve colors for
#' @param analysis type of analysis; options are "admixture" (default) or "hhsd"
#'
#' @return
#' @export
retrieve_kcols <- function(K_value = NULL, dataset_name = NULL, analysis = "admixture") {
  if (analysis == "admixture") {
    if (dataset_name == "ATL_MXPL") {
      if (K_value == 3) kcols = c(X1 = "#80a4bc", X3 = "#984625", X2 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#e0895a", X3 = "#80a4bc", X2 = "#984625", X4 = "#e0b693")
      if (K_value == 5) kcols = c(X5 = "#80a4bc", X1 = "#6f82b7", X4 = "#e0b693", X3 = "#984625", X2 = "#e0895a")
      if (K_value == 6) kcols = c(X4 = "#80a4bc", X5 = "#6f82b7", X1 = "#e0b693", X3 = "#984625", 
                                  X2 = "#e0895a", X6 = "#e97490")
      if (K_value == 7) kcols = c(X1 = "#80a4bc", X7 = "#6f82b7", X4 = "#e0b693", X6 = "#e97490", 
                                  X2 = "#e0895a", X5 = "#984625", X3 = "#c8b656")
    }
    
    if (dataset_name == "CENTAM") {
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "#924c62", X3 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#e0895a", X2 = "#73806d", X3 = "#82ccc8", X4 = "#924c62")
      if (K_value == 5) kcols = c(X1 = "#c8b656", X2 = "#e0895a", X3 = "#e0b693", X4 = "#924c62", X5 = "#73806d")
      if (K_value == 6) kcols = c(X1 = "#e0b693", X2 = "#c8b656", X3 = "#924c62", X4 = "#73806d", X5 = "#82ccc8", X6 = "#e0895a")
      if (K_value == 7) kcols = c(X1 = "#c8b656", X2 = "#e0895a", X3 = "#924c62", X4 = "#73806d", X5 = "#82ccc8", X6 = "#ba94a6", X7 = "#e0b693")
    }
    
    if (dataset_name == "PACMX") {
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "#ba94a6", X3 = "#e0895a")
      if (K_value == 4) kcols = c(X1 = "#ba94a6", X2 = "#73806d", X3 = "#6f82b7", X4 = "#e0895a")
      if (K_value == 5) kcols = c(X1 = "#73806d", X2 = "gray74", X3 = "#e0895a", X4 = "#ba94a6", X5 = "#6f82b7")
      if (K_value == 6) kcols = c(X1 = "#73806d", X2 = "#ba94a6", X3 = "#a8a2ca", X4 = "#e0895a", X5 = "#6f82b7", X6 = "#054051")
      if (K_value == 7) kcols = c(X1 = "#73806d", X2 = "#a8a2ca", X3 = "gray74", X4 = "#ba94a6", X5 = "#6f82b7", X6 = "#054051", X7 = "#e0895a")
      if (K_value == 8) kcols = c(X1 = "gray74", X2 = "#054051", X3 = "#e0895a", X4 = "gray30", X5 = "#73806d", X6 = "#428b9b", X7 = "#6f82b7", X8 = "#ba94a6")
      if (K_value == 9) kcols = c(X1 = "gray74", X2 = "#054051", X3 = "#ba94a6", X4 = "#6f82b7", X5 = "#e0895a", X6 = "#428b9b", X7 = "gray30", X8 = "#73806d", X9 = "#a8a2ca")
      if (K_value == 10) kcols = c(X1 = "#a8a2ca", X2 = "#ba94a6", X3 = "#73806d", X4 = "gray74", X5 = "#e0895a", X6 = "#e0b693", X7 = "#054051", X8 = "#428b9b", X9 = "gray30", X10 = "#6f82b7")
    }
    
    if (dataset_name == "forreri") {
      if (K_value == 2) kcols = c(X1 = "#73806d", X2 = "gray74")
      if (K_value == 3) kcols = c(X1 = "#73806d", X2 = "gray30", X3 = "gray74")
      if (K_value == 4) kcols = c(X1 = "gray30", X2 = "#73806d", X3 = "#eeb18f", X4 = "gray74")
      if (K_value == 5) kcols = c(X1 = "gray74", X2 = "#73806d", X3 = "#428b9b", X4 = "#d2a3a6", X5 = "gray30")
      if (K_value == 6) kcols = c(X1 = "gray30", X2 = "#eeb18f", X3 = "#d2a3a6", X4 = "gray74", X5 = "#428b9b", X6 = "#73806d")
    }
  }
  
  if (analysis == "hhsd") {
    mxpl = data.frame(pop = c("neov", "berl", "spec", "chic"),
                      color = c("#e97490", "#984625", "#80a4bc", "#6f82b7"),
                      assembly = "mxpl")
    forr = data.frame(pop = c("forr", "flor", "hill", "arce", "miad"),
                      color = c("gray74", "gray30", "#3b808f", "#d2a3a6", "#73806d"),
                      assembly = "forreri")
    foot = data.frame(pop = c("yava", "magn", "aten_short", "aten_long", "omil_GE", "omil_OA"),
                      color = c("#054051", "#f7cd5e", "#c0a06f", "#428b9b", "#a8a2ca", "#ba94a6"),
                      assembly = "foothills")

    kcols = bind_rows(mxpl, forr, foot)
  }

  return(kcols)
}

# "#80a4bc" = spectabilis
# "#6f82b7" = chichicuahutla
# "#e0895a" = macroglossa
# "#e0b693" = macroglossa2
# "#a8a2ca" = omiltemana oaxaca
# "#ba94a6" = omiltemana guerrero
# "#984625" = berlandieri
# "#e97490" = berlandieri neovolcanica
# "#054051" = yavapaiensis / foothills
# "#f7cd5e" = magnaocularis
# "#428b9b" = atenquique long (sp 8)
# "#c0a06f" = atenquique short (sp 7)
# "#82ccc8" = sp nov CR/PM
# "#73806d" = forreri miadis / southern forreri
# "#eeb18f" = miadis
# "gray30" = forreri floresi
# "gray74" = forreri forreri
# "#647777" = forreri4
# "#3b808f" = forreri hillisi
# "#924c62" = lenca
# "#c8b656" = taylori
# "#d2a3a6" = Arcelia form