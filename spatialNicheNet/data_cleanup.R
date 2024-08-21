sc_clean = function(){
  library(forcats)
  sc = readRDS('data/liver_mouseStSt_guilliams2022_withZonation.rds')
  sc$annot_fine_zonated = fct_collapse(sc$annot_fine_zonated, Fibroblast = c("Fibroblast 1","Fibroblast 2"))
  return(sc)
}

spatial_clean = function(){
  library(forcats)
  library(dplyr)
  ## RESOLVE A1-1
  a1 <- readRDS('data/A1-1.rds')
  levels(a1$annotationSave)[match("Kupffer cells",levels(a1$annotationSave))] <- "KCs"
  levels(a1$annotationSave)[match("central_vein_EC45",levels(a1$annotationSave))] <- "Central Vein Endothelial cells"
  levels(a1$annotationSave)[match("portal_vein_EC45",levels(a1$annotationSave))] <- "Portain Vein Endothelial cells"
  levels(a1$annotationSave)[match("capsular_fibroblasts45",levels(a1$annotationSave))] <- "Capsule fibroblasts"
  levels(a1$annotationSave)[match("HepatocytesCentral",levels(a1$annotationSave))] <- "Hepatocytes_central"
  levels(a1$annotationSave)[match("HepatocytesPortal",levels(a1$annotationSave))] <- "Hepatocytes_portal"
  levels(a1$annotationSave)[match("LSEC Portal",levels(a1$annotationSave))] <- "LSECs_portal"
  levels(a1$annotationSave)[match("LSEC Central",levels(a1$annotationSave))] <- "LSECs_central"
  levels(a1$annotationSave)[match("VSMC",levels(a1$annotationSave))] <- "VSMCs"
  levels(a1$annotationSave)[match("stellateAll",levels(a1$annotationSave))] <- "Stellate cells"
  levels(a1$annotationSave)[match("stellate PeriCentral",levels(a1$annotationSave))] <- "Stellate cells_central"
  levels(a1$annotationSave)[match("stellate PeriPortal",levels(a1$annotationSave))] <- "Stellate cells_portal"
  levels(a1$annotationSave)[match("FibroblastsCentral",levels(a1$annotationSave))] <- "Fibroblast 1"
  levels(a1$annotationSave)[match("fibroblastPortal",levels(a1$annotationSave))] <- "Fibroblast 2"
  levels(a1$annotationSave)[match("Portal LAM macrophages",levels(a1$annotationSave))] <- "MoMac1"
  levels(a1$annotationSave)[match("Capsule and Central Vein Mac",levels(a1$annotationSave))] <- "MoMac2"
  levels(a1$annotationSave)[match("LECs",levels(a1$annotationSave))] <- "Lymphatic Endothelial cells"
  a1@meta.data[['annotationSave']] = fct_collapse(a1@meta.data[['annotationSave']], Fibroblast = c("Fibroblast 1","Fibroblast 2"))

  ## RESOLVE A1-2
  a2 <- readRDS('data/A1-2.rds')
  levels(a2$annotationSave)[match("Kupffer cells",levels(a2$annotationSave))] <- "KCs"
  levels(a2$annotationSave)[match("central_vein_EC45",levels(a2$annotationSave))] <- "Central Vein Endothelial cells"
  levels(a2$annotationSave)[match("portal_vein_EC45",levels(a2$annotationSave))] <- "Portain Vein Endothelial cells"
  levels(a2$annotationSave)[match("capsular_fibroblasts45",levels(a2$annotationSave))] <- "Capsule fibroblasts"
  levels(a2$annotationSave)[match("HepatocytesCentral",levels(a2$annotationSave))] <- "Hepatocytes_central"
  levels(a2$annotationSave)[match("HepatocytesPortal",levels(a2$annotationSave))] <- "Hepatocytes_portal"
  levels(a2$annotationSave)[match("LSEC Portal",levels(a2$annotationSave))] <- "LSECs_portal"
  levels(a2$annotationSave)[match("LSEC Central",levels(a2$annotationSave))] <- "LSECs_central"
  levels(a2$annotationSave)[match("VSMC",levels(a2$annotationSave))] <- "VSMCs"
  levels(a2$annotationSave)[match("stellateAll",levels(a2$annotationSave))] <- "Stellate cells"
  levels(a2$annotationSave)[match("stellate PeriCentral",levels(a2$annotationSave))] <- "Stellate cells_central"
  levels(a2$annotationSave)[match("stellate PeriPortal",levels(a2$annotationSave))] <- "Stellate cells_portal"
  levels(a2$annotationSave)[match("FibroblastsCentral",levels(a2$annotationSave))] <- "Fibroblast 1"
  levels(a2$annotationSave)[match("fibroblastPortal",levels(a2$annotationSave))] <- "Fibroblast 2"
  levels(a2$annotationSave)[match("Portal LAM macrophages",levels(a2$annotationSave))] <- "MoMac1"
  levels(a2$annotationSave)[match("Capsule and Central Vein Mac",levels(a2$annotationSave))] <- "MoMac2"
  a2@meta.data[['annotationSave']] = fct_collapse(a2@meta.data[['annotationSave']], Fibroblast = c("Fibroblast 1","Fibroblast 2",'FibroblastAll'))
  levels(a2$annotationSave)[match("LECs",levels(a2$annotationSave))] <- "Lymphatic Endothelial cells"

  ## RESOLVE C2
  c2 <- readRDS('data/C2.rds')
  levels(c2$annotationSave)[match("Kupffer cells",levels(c2$annotationSave))] <- "KCs"
  levels(c2$annotationSave)[match("central_vein_EC45",levels(c2$annotationSave))] <- "Central Vein Endothelial cells"
  levels(c2$annotationSave)[match("portal_vein_EC45",levels(c2$annotationSave))] <- "Portain Vein Endothelial cells"
  levels(c2$annotationSave)[match("capsular_fibroblasts45",levels(c2$annotationSave))] <- "Capsule fibroblasts"
  levels(c2$annotationSave)[match("HepatocytesCentral",levels(c2$annotationSave))] <- "Hepatocytes_central"
  levels(c2$annotationSave)[match("HepatocytesPortal",levels(c2$annotationSave))] <- "Hepatocytes_portal"
  levels(c2$annotationSave)[match("LSEC Portal",levels(c2$annotationSave))] <- "LSECs_portal"
  levels(c2$annotationSave)[match("LSEC Central",levels(c2$annotationSave))] <- "LSECs_central"
  levels(c2$annotationSave)[match("VSMC",levels(c2$annotationSave))] <- "VSMCs"
  levels(c2$annotationSave)[match("stellateAll",levels(c2$annotationSave))] <- "Stellate cells"
  levels(c2$annotationSave)[match("stellate PeriCentral",levels(c2$annotationSave))] <- "Stellate cells_central"
  levels(c2$annotationSave)[match("stellate PeriPortal",levels(c2$annotationSave))] <- "Stellate cells_portal"
  levels(c2$annotationSave)[match("Portal LAM macrophages",levels(c2$annotationSave))] <- "MoMac1"
  levels(c2$annotationSave)[match("Capsule and Central Vein Mac",levels(c2$annotationSave))] <- "MoMac2"
  c2@meta.data[['annotationSave']] = fct_collapse(c2@meta.data[['annotationSave']], Fibroblast = c("FibroblastsCentral","fibroblastPortal",'FibroblastAll'))
  levels(c2$annotationSave)[match("LECs",levels(c2$annotationSave))] <- "Lymphatic Endothelial cells"

  ## RESOLVE D2-1
  d2 <- readRDS('data/D2-1.rds')
  levels(d2$annotationSave)[match("Kupffer cells",levels(d2$annotationSave))] <- "KCs"
  levels(d2$annotationSave)[match("central_vein_EC45",levels(d2$annotationSave))] <- "Central Vein Endothelial cells"
  levels(d2$annotationSave)[match("portal_vein_EC45",levels(d2$annotationSave))] <- "Portain Vein Endothelial cells"
  levels(d2$annotationSave)[match("capsular_fibroblasts45",levels(d2$annotationSave))] <- "Capsule fibroblasts"
  levels(d2$annotationSave)[match("HepatocytesCentral",levels(d2$annotationSave))] <- "Hepatocytes_central"
  levels(d2$annotationSave)[match("HepatocytesPortal",levels(d2$annotationSave))] <- "Hepatocytes_portal"
  levels(d2$annotationSave)[match("LSEC Portal",levels(d2$annotationSave))] <- "LSECs_portal"
  levels(d2$annotationSave)[match("LSEC Central",levels(d2$annotationSave))] <- "LSECs_central"
  levels(d2$annotationSave)[match("VSMC",levels(d2$annotationSave))] <- "VSMCs"
  levels(d2$annotationSave)[match("stellateAll",levels(d2$annotationSave))] <- "Stellate cells"
  levels(d2$annotationSave)[match("stellate PeriCentral",levels(d2$annotationSave))] <- "Stellate cells_central"
  levels(d2$annotationSave)[match("stellate PeriPortal",levels(d2$annotationSave))] <- "Stellate cells_portal"
  levels(d2$annotationSave)[match("FibroblastsCentral",levels(d2$annotationSave))] <- "Fibroblast 1"
  levels(d2$annotationSave)[match("fibroblastPortal",levels(d2$annotationSave))] <- "Fibroblast 2"
  levels(d2$annotationSave)[match("Portal LAM macrophages",levels(d2$annotationSave))] <- "MoMac1"
  levels(d2$annotationSave)[match("Capsule and Central Vein Mac",levels(d2$annotationSave))] <- "MoMac2"
  d2@meta.data[['annotationSave']] = fct_collapse(d2@meta.data[['annotationSave']], Fibroblast = c("Fibroblast 1","Fibroblast 2",'FibroblastAll'))
  levels(d2$annotationSave)[match("LECs",levels(d2$annotationSave))] <- "Lymphatic Endothelial cells"

  return(c(a1, a2, c2, d2))
}