library(ggplot2)
library(reshape2)
library(ggsignif)
library(ggpubr)
theme_mine <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.text.x = element_text(size=16),
      strip.text.y = element_text(size=16),
      strip.background = element_rect(colour="black", fill="white"),
      #axis.text.x = element_text(size=10,face="bold",angle=45,hjust=1,vjust=1),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=10,hjust=1,face="bold"),
      axis.ticks.x =  element_line(colour = "black"), 
      axis.ticks.y =  element_line(colour = "black"), 
      axis.title.x =  element_blank(),
      axis.title.y = element_text(size=16,angle=90,face="bold",vjust=2),
      legend.position = "top", 
      panel.grid.major.y = element_line(colour = "lightgray"), 
      panel.grid.minor = element_blank(),
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line = element_line(colour = "black")
    )
}
clinicalData = read.delim(
  "../../June26_2017/ClinicalAnnotations_onlyPassed_2018April20.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
melt_clindata = melt(clinicalData,measure.vars = 'Mutations.Mb')
mskdata = read.delim(
  "../../June26_2017/cBio_MSKIMPACT_10000/study_view_clinical_data_withMB.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
melt_mskdata = melt(mskdata,measure.vars = 'Mutation_Burden')
cdata = merge(melt_mskdata,melt_clindata,by=c('variable','value'),all=TRUE)
cdata$value = as.numeric(cdata$value)
#cdata[cdata$variable=="Mutation_Burden"]<-"MSK-IMPACT"
#cdata[cdata$variable=="Mutations.Mb"]<-"CSF-ctDNA"

p = ggplot(cdata,aes(x=variable,y=value,fill=variable)) + 
  geom_boxplot(position = "identity") + 
  #scale_y_continuous(breaks = c(10,20,seq(50,450,50)),limits = c(0,450)) +
  scale_y_log10(breaks=c(seq(0,10,1),50,100,150,200,400))+
  #scale_y_log10(breaks = c(seq(0,20,5),100,200,300,400,500))+
  scale_fill_manual(values = c("#A9A9A9","#606060"),labels=c("MSK-IMPACT","CSF ctDNA")) +
  #scale_x_discrete(labels=c("MSK-IMPACT","CSF ctDNA"))+
  ylab("Mutations/Mb") + 
  geom_signif(comparisons = list(c("Mutation_Burden", "Mutations.Mb")), 
              map_signif_level=TRUE)
p + theme_mine()
#ggboxplot(cdata, x = "variable", y = "value",fill = "variable", 
#          palette = c("#A9A9A9","#606060"), outlierColor="red",
#          theme = theme_mine()) +
#  stat_compare_means(method = "t.test")
ggsave("mutations_mb_till10_v1_2018April20.pdf",width = 6,height = 10)
