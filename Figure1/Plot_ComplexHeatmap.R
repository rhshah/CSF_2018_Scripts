library(ComplexHeatmap)
library("beyonce")
#sample_order = scan("../MergedPlotWithcBio/PNET_PatientOrder_PID.txt",what = "character")
clinicalData = read.delim(
  "../ClinicalAnnotations_onlyPassed_2018April20.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
#nclinicalData = clinicalData[order(clinicalData$Histologic_Grade,clinicalData$Location_of_IMPACT_testing,clinicalData$Somatostatin_analog_Therapy,clinicalData$Targeted_therapy,clinicalData$Chemotherapy),]
nclinicalData = clinicalData
sampleOrder = as.vector(nclinicalData$Tumor_Sample_Barcode)
#nclinicalData<-clinicalData[match(sample_order, clinicalData$Sample_Id),]
#nd_clinicalData <- subset(nclinicalData, select = c(Location_of_IMPACT_testing,Tumor_Differentiation,Histologic_Grade))
rownames(nclinicalData) <- 1:nrow(nclinicalData)

df = data.frame(
  Histology = nclinicalData$Histology
)
col = list(
  Histology =  c(
    "Anaplastic oligodendroglioma, IDH-mutant and 1p/19q co-deleted"='#DC7633',
    "Oligodendroglioma, IDH-mutant and 1p/19q co-deleted"='#F7DC6F',
    "Anaplastic Astrocytoma, IDH-mutant"='#7B241C',
    "Anaplastic Astrocytoma, IDH-WT"='#F5B7B1',
    "Diffuse astrocytoma, IDH-mutant"='#D2B4DE',
    "Diffuse astrocytoma, IDH-WT"='#633974',
    "Glioblastoma, IDH-mutant"='#AED6F1',
    "Glioblastoma, IDH-WT"='#1F618D'
   

  )
)
ha = HeatmapAnnotation(
  df = df,
  col = col,
  Mutation_Burden = anno_barplot(
    nclinicalData$Mutations.Mb,
    axis = TRUE,
    ylim = c(0,max(nclinicalData$Mutations.Mb)),
    gp = gpar(fill = "#ACBAC7")
  ),
  annotation_name_side = "right",
  annotation_name_offset = unit(5, "mm"),
  annotation_name_gp = gpar(fontface = "bold"),
  gap = unit(c(2,2,2), "mm"),
  annotation_height = c(3, 3, 20),
  #height = unit(10,"cm"),
  
  show_annotation_name = TRUE
)
zero_row_mat = matrix(nrow = 0,
                      ncol = length(nclinicalData$Tumor_Sample_Barcode))
#colnames(zero_row_mat) = letters[1:90]
ht = Heatmap(zero_row_mat,
             top_annotation = ha,
             column_title = "only annotations")
pdf(file = "annoHeatmap_byCustomSort_v2_April2018.pdf",
    width = 20,
    height = 10)
draw(ht)
dev.off()
#mat = read.delim("../portal_processed_data/TopMutatedGenes/InputForComplexHeatmap_Transposed.txt",header = TRUE,stringsAsFactors=FALSE, sep = "\t")
#mat = read.delim("../MergedPlotWithcBio/Updated_Data_October262016/TestcBio.txt",header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat = read.delim(
  "../ComplexHeatMap_onlyPassed_2018April20.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)
mat = mat[match(sampleOrder, mat$Case_ID),]
mat[is.na(mat)] = ""
rownames(mat) = sampleOrder
#rownames(mat) = mat[, 1]
mat = mat[, -1]
#mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x,
              y,
              w - unit(0.25, "mm"),
              h - unit(0.25, "mm"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "#008000", col = NA))
  },
  InFrame = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "#708090", col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "black", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x,
              y,
              w - unit(0.25, "mm"),
              h - unit(0.25, "mm"),
              gp = gpar(fill = "lightcyan4", col = NA))
  },
  Multi = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "#808000", col = NA))
  },
  PRESENT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "#556B2F", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x,
              y,
              w - unit(0.25, "mm"),
              h - unit(0.25, "mm"),
              gp = gpar(fill = "red", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x,
              y,
              w - unit(0.25, "mm"),
              h - unit(0.25, "mm"),
              gp = gpar(fill = "blue", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.25, "mm"), h * 0.33, gp = gpar(fill = "purple", col = NA))
  }
  
)
rowOrder = c(
  "TERT",
  "TP53",
  "IDH1",
  "CDKN2A.B",
  "EGFR",
  "ATRX",
  "PTEN",
  "CIC",
  "NF1",
  "PIK3CA",
  "CDK4",
  "PIK3R1",
  "RB1",
  "PDGFRA",
  "SOX2",
  "NOTCH1",
  "KIT",
  "MDM2",
  "SETD2",
  "PTPN11",
  "FUBP1",
  "KMT2D",
  "FAT1",
  "KDR"
)
col = c(
  "MUT" = "#008000",
  "TRUNC" = "black",
  "InFrame" = "#708090",
  "Multi" = "#808000",
  "PRESENT"="#556B2F",
  "AMP" = "red",
  "HOMDEL" = "blue",
  "DEL" = "lightcyan4",
  "Fusion" = "purple"
 
)
ht_list = oncoPrint(
  mat,
  get_type = function(x)
    strsplit(x, ";")[[1]],
  alter_fun = alter_fun,
  col = col,
  row_order = rowOrder,
  column_order = sampleOrder,
  remove_empty_columns = F,
  top_annotation = ha,
  top_annotation_height = unit(20, "cm"),
  heatmap_legend_param = list(
    title = "Alternations",
    at = c("MUT", "TRUNC", "InFrame", "Multi","PRESENT","AMP", "HOMDEL", "DEL", "Fusion"),
    labels = c(
      "Missense Mutation",
      "Truncating Mutations",
      "In-Frame Mutations",
      "Multiple Mutations",
      "Mutations below threshold",
      "Amplification",
      "Deep deletion",
      "Intragenic Deletion",
      "EGFR vIII"
    ),
    ncol = 1,
    title_position = "topcenter"
  ),
  
  row_barplot_width = unit(8, "cm")
  #barplot_ignore = c("LOH"),
)
decorate_annotation("Mutation_Burden", {
  grid.text("Mutations/Mb", unit(10, "mm"), just = "bottom", rot = 90)
  grid.lines(c(0, 1), unit(c(0.2, 0.2), "native"), gp = gpar(lty = 2, col = "blue"))
})
pdf(file = "OverallOncoPrint_byCustomSort_onlypass_v2_April2018.pdf",
    width = 22,
    height = 20)
draw(ht_list,row_sub_title_side = "left")
dev.off()
