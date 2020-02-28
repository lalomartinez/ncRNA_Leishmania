library(UpSetR)
library(png)
library(ggplot2)
library(ggpubr)

data <- read.table("~/Desktop/Leishmania_info_06-01-20/ncRNAs_Leishmania_spp/clustering_nc/Leishmania_nc.clstr.strain_binary_matrix",
                   header=T,sep="\t",row.names = 1)


set=c("LaetL147","LamaM2269","LaraLEM1108","LbraM2903","LbraM2904","LdonAG83","LdonBHU1220","LdonBPK282A1",
      "LdonPasteur","LenrLEM3045","LgerLEM452","LinfJPCM5","LinfTR01","LmajFriedlin","LmajLV39c5","LmajSD75",
      "LmexU1103","LpanL13","LpanPSC1","LperLEM1537","LperPAB4377","LspLD974","LspLEM2494","LtarParTarII","LtroL590","LturLEM423")


lv_sets= list(c("LdonAG83","LdonBHU1220","LdonBPK282A1","LdonPasteur","LinfJPCM5","LinfTR01"))
viannia_sets= list(c("LbraM2903","LbraM2904","LpanL13","LpanPSC1","LperLEM1537","LperPAB4377"))
flist= list(set)

setEPS()
postscript("~/Desktop/Leishmania_info_06-01-20/Imagenes_ncRNAs_paper/Shared_ncRNAS_cluster.eps", width = 12.0, height = 9.0)
  upset(data, sets = set,
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
      mainbar.y.label = "Shared ncRNAs clusters",
      queries = list(list(query = intersects, params = lv_sets, active = T, color ="#5B7CA7"),
                     list(query = intersects, params = viannia_sets, active=T, color="#DB9942"),
                     list(query = intersects, params = flist, active=T, color="#985A71")
                     )
      )

dev.off()

