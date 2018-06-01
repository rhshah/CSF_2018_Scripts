library(ggplot2)
library(ggthemes)
library(ggtech)
library(beyonce)
library(reshape2)
library(dplyr)
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
p38<- read.delim("./Mutation_Table/s_IM_GBM_38.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p38$s_IM_GBM_38_PL = gsub("[^0-9.]", "", p38$s_IM_GBM_38_PL)
p38$s_IM_GBM_38_CSFb = gsub("[^0-9.]", "", p38$s_IM_GBM_38_CSFb)

m38DF = melt(p38,measure.vars = c("s_IM_GBM_38_CSFb","s_IM_GBM_38_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m38DF$value <- as.numeric(m38DF$value)
p38_t1 <- t.test(
  x = m38DF$value[m38DF$variable == "s_IM_GBM_38_PL"],
  y = m38DF$value[m38DF$variable == "s_IM_GBM_38_CSFb"]
)
p29<- read.delim("./Mutation_Table/s_IM_GBM_29.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p29$s_IM_GBM_29_CSF = gsub("[^0-9.]", "", p29$s_IM_GBM_29_CSF)
p29$s_IM_GBM_29_PL = gsub("[^0-9.]", "", p29$s_IM_GBM_29_PL)
m29DF = melt(p29,measure.vars = c("s_IM_GBM_29_CSF","s_IM_GBM_29_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m29DF$value <- as.numeric(m29DF$value)
p29_t1 <- t.test(
  x = m29DF$value[m29DF$variable == "s_IM_GBM_29_PL"],
  y = m29DF$value[m29DF$variable == "s_IM_GBM_29_CSF"]
)
p40<- read.delim("./Mutation_Table/s_IM_GBM_40.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p40$s_IM_GBM_40_PL = gsub("[^0-9.]", "", p40$s_IM_GBM_40_PL)
p40$s_IM_GBM_40_CSFa = gsub("[^0-9.]", "", p40$s_IM_GBM_40_CSFa)
m40DF = melt(p40,measure.vars = c("s_IM_GBM_40_CSFa","s_IM_GBM_40_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m40DF$value <- as.numeric(m40DF$value)
p40_t1 <- t.test(
  x = m40DF$value[m40DF$variable == "s_IM_GBM_40_PL"],
  y = m40DF$value[m40DF$variable == "s_IM_GBM_40_CSFa"]
)
p46<- read.delim("./Mutation_Table/s_IM_GBM_46.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p46$s_IM_GBM_46_PL = gsub("[^0-9.]", "", p46$s_IM_GBM_46_PL)
p46$s_IM_GBM_46_CSFa = gsub("[^0-9.]", "", p46$s_IM_GBM_46_CSFa)
m46DF = melt(p46,measure.vars = c("s_IM_GBM_46_CSFa","s_IM_GBM_46_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m46DF$value <- as.numeric(m46DF$value)
p46_t1 <- t.test(
  x = m46DF$value[m46DF$variable == "s_IM_GBM_46_PL"],
  y = m46DF$value[m46DF$variable == "s_IM_GBM_46_CSFa"]
)

p55<- read.delim("./Mutation_Table/s_IM_GBM_55.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p55$s_IM_GBM_55_PL = gsub("[^0-9.]", "", p55$s_IM_GBM_55_PL)
p55$s_IM_GBM_55_CSF = gsub("[^0-9.]", "", p55$s_IM_GBM_55_CSF)
m55DF = melt(p55,measure.vars = c("s_IM_GBM_55_CSF","s_IM_GBM_55_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m55DF$value <- as.numeric(m55DF$value)
p55_t1 <- t.test(
  x = m55DF$value[m55DF$variable == "s_IM_GBM_55_PL"],
  y = m55DF$value[m55DF$variable == "s_IM_GBM_55_CSF"]
)

p73<- read.delim("./Mutation_Table/s_IM_GBM_73.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p73$s_IM_GBM_73_PL = gsub("[^0-9.]", "", p73$s_IM_GBM_73_PL)
p73$s_IM_GBM_73_CSFa = gsub("[^0-9.]", "", p73$s_IM_GBM_73_CSFa)
m73DF = melt(p73,measure.vars = c("s_IM_GBM_73_CSFa","s_IM_GBM_73_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m73DF$value <- as.numeric(m73DF$value)
p73_t1 <- t.test(
  x = m73DF$value[m73DF$variable == "s_IM_GBM_73_PL"],
  y = m73DF$value[m73DF$variable == "s_IM_GBM_73_CSFa"]
)

p63<- read.delim("./Mutation_Table/s_IM_GBM_63.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
p63$s_IM_GBM_63_PL = gsub("[^0-9.]", "", p63$s_IM_GBM_63_PL)
p63$s_IM_GBM_63_CSFa = gsub("[^0-9.]", "", p63$s_IM_GBM_63_CSFa)
m63DF = melt(p63,measure.vars = c("s_IM_GBM_63_CSFa","s_IM_GBM_63_PL"),id.vars = c("Gene_aa","Chr","Pos", "Ref","Alt"))
m63DF$value <- as.numeric(m63DF$value)
p63_t1 <- t.test(
  x = m63DF$value[m63DF$variable == "s_IM_GBM_63_PL"],
  y = m63DF$value[m63DF$variable == "s_IM_GBM_63_CSFa"]
)

# Merge only once
#mDF = rbind(m29DF,m38DF,m40DF,m46DF,m55DF,m73DF,m63DF)
#write.table(mDF,row.names = FALSE,file="csf_plasma.txt",quote = FALSE,sep = "\t")
mDF = read.delim("./csf_plasma.txt",header = T,sep = '\t')
sample_order <- c("s_IM_GBM_29_CSF","s_IM_GBM_29_PL",
                  "s_IM_GBM_38_CSFb","s_IM_GBM_38_PL",
                  "s_IM_GBM_40_CSFa","s_IM_GBM_40_PL",
                  "s_IM_GBM_46_CSFa","s_IM_GBM_46_PL",
                  "s_IM_GBM_55_CSF","s_IM_GBM_55_PL",
                  "s_IM_GBM_63_CSFa","s_IM_GBM_63_PL",
                  "s_IM_GBM_73_CSFa","s_IM_GBM_73_PL"
                  )
mDF$variable <- factor(mDF$variable, levels = sample_order)
p = ggplot(mDF,aes(x=variable,y=value,fill=SampleType)) + 
  geom_boxplot(position = "identity") +
  scale_y_continuous(breaks = seq(0,1,0.05),limits = c(0,1)) +
  scale_fill_manual(values = c("#606060","#A9A9A9")) +
  ylab("Variant Allele Frequency")
p = p + 
  geom_text(data = mDF,aes(x = "s_IM_GBM_29_CSF", y = 0.95), label = paste("P = ", format.pval(p29_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_38_CSFb", y = 0.95), label = paste("P = ", format.pval(p38_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_40_CSFa", y = 0.95), label = paste("P = ", format.pval(p40_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_46_CSFa", y = 0.95), label = paste("P = ", format.pval(p46_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_55_CSF", y = 0.95), label = paste("P = ", format.pval(p55_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_63_CSFa", y = 0.95), label = paste("P = ", format.pval(p63_t1$p.value, 2)),size=5) +
  geom_text(data = mDF,aes(x = "s_IM_GBM_73_CSFa", y = 0.95), label = paste("P = ", format.pval(p73_t1$p.value, 2)),size=5)
pdf("csf_plasma_pval_v2.pdf",width = 12,height = 8)
p + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()


p = ggplot(mDF,aes(x=variable,y=value,fill=SampleType)) + 
  geom_boxplot(position = "identity") +
  scale_y_continuous(breaks = seq(0,1,0.05),limits = c(0,1)) +
  scale_fill_manual(values = c("#606060","#A9A9A9")) +
  ylab("Variant Allele Frequency") 
p = p + facet_grid(~Pt_ID,scales = "free",switch = 'x') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) 

pdf("csf_plasma_v2.pdf",width = 10,height = 8)
p
dev.off() 

#not done
#p63<- read.delim("../s_IM_GBM_63.txt",header = T,sep = "\t",stringsAsFactors = FALSE)
#p63$s_IM_GBM_63_PZ = gsub("[^0-9.]", "", p63$s_IM_GBM_63_PZ)
#p63$s_IM_GBM_63 = gsub("[^0-9.]", "", p63$s_IM_GBM_63)