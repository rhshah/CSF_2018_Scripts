library(xlsx)
library("plyr")
library(dplyr)
library(ggplot2)
library(BBmisc)
library("reshape2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("scales")
library(beyonce)
library(stringr)
theme_mine <- function(base_size = 12,
                       base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text.x = element_text(
        size = 10,
        hjust = 0,
        vjust = 0.5,
        angle = 330,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 14,
        hjust = 1,
        face = "bold"
      ),
      axis.ticks.x =  element_blank(),
      axis.ticks.y =  element_blank(),
      axis.title.x =  element_blank(),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal"
    )
}
makePlot1 <- function(dataDF, outfile, dflen) {
  pal <- beyonce_palette(90, type = "continuous")
  #print (dataDF)
  print(outfile)
  datm = unique(melt(dataDF, id.vars = c("Gene_aa", "Chr", "Pos", "Ref", "Alt")))
  #print(datm)
  datm$Gene_aa = factor(datm$Gene_aa, levels = dataDF$Gene_aa)
  datm$variable = factor(datm$variable, levels = unique(datm$variable))
  datm$value = as.numeric(datm$value)
  base_size <- 9
  p <- ggplot(datm, aes(x = factor(Gene_aa), y = variable))
  p + geom_tile(aes(fill = value)) + geom_point(aes(size = value), colour = "black", show_guide =
                                                  TRUE) +
    scale_fill_gradient2(
      limits = c(0, 1),
      low = muted("red"),
      mid = "white",
      high = "steelblue",
      guide = guide_colorbar(
        title = "Variant Allele Frequency (VAF): ",
        title.theme = element_text(
          angle = 0,
          size = 16,
          face = "bold"
        ),
        title.position = "top",
        label.theme = element_text(
          angle = 0,
          size = 12,
          face = "bold"
        )
      )
    ) +
    xlab("") + ylab("") +
    coord_equal() +
    scale_size_continuous(breaks = seq(0, 1, 0.05), limits = c(0, 1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limit = (rev(levels(datm$variable))), expand = c(0, 0)) +
    theme_mine()
  
  ggsave(outfile, width = 25, height = 15)
  
  
}


makePlot <- function(dataDF, testDF, outfile, dflen) {
  pal <- beyonce_palette(90, type = "continuous")
  #print (dataDF)
  print(outfile)
  datm = unique(melt(dataDF, id.vars = c("Gene_aa", "Chr", "Pos", "Ref", "Alt")))
  dataT = unique(melt(textDF, id.vars = c("Gene_aa", "Chr", "Pos", "Ref", "Alt")))
  #print(datm)
  if (dflen <= 10) {
    fsize = 10
    
  }
  else if (dflen > 11 && dflen <= 20) {
    fsize = 8
    
  }
  else if (dflen > 20 && dflen <= 30) {
    fsize = 4
    
  }
  else{
    fsize = 3
    
  }
  #datm$Gene_aa = factor(datm$Gene_aa, levels = dataDF$Gene_aa)
  #datm$variable = factor(datm$variable, levels = unique(datm$variable))
  datm$value = as.numeric(datm$value)
  base_size <- 9
  p <- ggplot(datm, aes(x = factor(Gene_aa), y = variable))
  p + geom_tile(aes(fill = value)) + geom_tile(aes(fill = value), colour = "black", show_guide =
                                                 FALSE) +
    geom_text(data = dataT, aes(label = value), size = fsize) +
    scale_fill_gradient2(
      limits = c(0, 1),
      low = muted("red"),
      mid = "white",
      high = "steelblue",
      guide = guide_colorbar(
        title = "Variant Allele Frequency (VAF): ",
        title.theme = element_text(
          angle = 0,
          size = 16,
          face = "bold"
        ),
        title.position = "top",
        label.theme = element_text(
          angle = 0,
          size = 12,
          face = "bold"
        )
      )
    ) +
    xlab("") + ylab("") +
    coord_equal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limit = (rev(levels(datm$variable))), expand = c(0, 0)) +
    theme_mine()
  
  ggsave(outfile, width = 25, height = 15)
  
  
}
checkAF <- function(aVals) {
  aVals = as.numeric(aVals)
  aVF = paste("C", round(aVals[4],3), sep = ":")
  
  if (aVals[1] >= 20) {
    if (aVals[3] >= 2 && aVals[4] >= 0.01) {
      if (aVals[3] >= 8 && aVals[4] >= 0.02) {
        aVF = aVF
      }
      else{
        aVF = paste("P", round(aVals[4],3), sep = ":")
      }
    }
    else{
      aVF =  paste("N", round(aVals[4],3), sep = ":")
    }
    
  }
  else{
    aVF =  paste("N", round(aVals[4],3), sep = ":")
  }
  return(aVF)
}


outDir = getwd()
titleDF = read.delim(
  #"/Users/shahr2/Documents/MSKCC/CMO/CSF_Analysis/5500_CM/Proj_05500_CM_Zp_title.txt",
  #"/Users/rshah22/Documents/Projects/MSKCC/CSF//MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF3/Proj_05500_CSF3_Zp_title_v1.txt",
  #"/Users/rshah22/Documents/Projects/MSKCC/CSF//MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF2/Proj_05500_CSF2_Zp_title_v1.txt",
  "/Users/rshah22/Documents/Projects/MSKCC/CSF//MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF1/Proj_05500_CSF1_Zp_title_v1.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)
titleDF = titleDF[with (titleDF,!grepl("Normal", Class)), ]
dataDF = read.delim(
  #"/Users/shahr2/Documents/MSKCC/CMO/CSF_Analysis/5500_CM/annotated_exonic_variants.txt",
  #"/Users/rshah22/Documents/Projects/MSKCC/CSF/MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF3/annotated_variants_proj05500_csf3.txt",
  #"/Users/rshah22/Documents/Projects/MSKCC/CSF/MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF2/annotated_variants_proj05500_csf2.txt",
  "/Users/rshah22/Documents/Projects/MSKCC/CSF/MergedAnalysis_3March2017/June26_2017/NewRuns/Proj_05500_CSF1/annotated_variants_proj05500_csf1.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)
dataDF <- dataDF[dataDF$Call_Confidence == "HIGH", ]
uPID <- unique(dataDF$PatientID)
#gDF <- group_by(dataDF,PID)
allOutDF = list()
allTextDF = list()
for (pId in uPID) {
  pData = dataDF[dataDF$PatientID  == pId, ]
  numrows = nrow(pData)
  titleData = titleDF[titleDF$Patient_ID == pId, ]
  sampleIDs = sort(unique(titleData$Sample_ID))

  for (sampleID in sampleIDs){
    cat(sampleID,"\n")
    outDF = data.frame(
      Sample=character(numrows),
      Chr = character(numrows),
      Pos = integer(numrows),
      Ref = character(numrows),
      Alt = character(numrows),
      Gene = character(numrows),
      cDNAchange=character(numrows),
      AAchange = character(numrows),
      TotalDepth = double(numrows),
      VariantAlleleDepth = double(numrows),
      VariantAlleleFrequency = double(numrows),
      StatusOfMutation = character(numrows),
      stringsAsFactors = FALSE
    )
    for (i in 1:nrow(pData)) {
      outDF$Sample[i] = sampleID
      outDF$Gene[i] = pData[i, "Gene"]
      outDF$AAchange[i] = pData[i, "AAchange"]
      outDF$cDNAchange[i] = pData[i,"cDNAchange"]
      outDF$Chr[i] = as.character(pData[i, "Chrom"])
      outDF$Pos[i] = pData[i, "Start"]
      outDF$Ref[i] = pData[i, "Ref"]
      outDF$Alt[i] = pData[i, "Alt"]
      aVals = unlist((str_match_all(
      pData[i, sampleID][[1]], "[+-]?([0-9]*[.])?[0-9]+"
    )))[c(1, 2, 3, 4)]
      #cat(aVals,"\n")
      aVF = checkAF(aVals)
      outDF$TotalDepth[i] = as.character(aVals[1])
      outDF$VariantAlleleDepth[i] = as.character(aVals[3])
      outDF$VariantAlleleFrequency[i] = as.character(aVals[4])
      if(str_detect(aVF, "C:")){
        status="Called"
      }
      if(str_detect(aVF, "P:")){
        status="Present"
      }
      if(str_detect(aVF, "N:")){
        status="Absent"
      }
      outDF$StatusOfMutation[i] = status
    }
    colnames(outDF) <- c("Sample","Chr", "Pos", "Ref", "Alt","Gene","cDNAchange","AAchange","Total Depth(DP)","Variant Allele Depth(AD)","Variant Allele Frequency(VAF)","Mutation Status")
    allOutDF[[sampleID]] = unique(outDF)
  }
}
for (i in 1:length(allOutDF))
{
  samplename = names(allOutDF[i]);
  outname = paste(samplename ,"v1" ,".txt",sep = "_")
  write.table(unique(allOutDF[[i]]), file = outname
              ,quote = FALSE,sep = "\t",row.names=FALSE)
  
}