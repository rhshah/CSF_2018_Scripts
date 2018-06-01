library(ggplot2)
library(reshape2)
library(ggsignif)
theme_mine <- function(base_size = 12, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.text.x = element_text(size=16),
      strip.text.y = element_text(size=16),
      strip.background = element_rect(colour="black", fill="white"),
      axis.text.x = element_text(size=10,face="bold",angle=45,hjust=1,vjust=1),
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

rd<-read.delim(file="../Plot_MSKIMPACT_vs_CSF_2018April25.txt",sep = "\t",header = T,stringsAsFactors = F)
rd2 = rd
rd2$Gene <- factor(rd$Gene,levels = rd$Gene[order(rd$MSK.IMPACT,decreasing = T)])
mrd = melt(rd2,measure.vars = c("MSK.IMPACT","CSF.ctDNA"))
ggplot(mrd,aes(x=Gene,y=value,fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1))+
  scale_fill_manual(values = c("#A9A9A9","#606060"))+
  xlab("Genes") + ylab("Frequecy") +
  theme_mine()
ggsave("Plot_Frequecy_MSK_CSF_v1_2018April25.pdf",height = 5,width = 10)