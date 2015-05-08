########
# ggplot personalized Theme
# Quentin Schorpp
# 16.04.2015
########

library(ggplot2)
library(grid)
library(extrafont)

mytheme = 
  theme_bw() + 
  theme(plot.title = element_text(size=11,face="bold", family="Times New Roman"),
        axis.title.x = element_text(size=8, face="bold", family="Times New Roman"),
        axis.title.y = element_text(size=8, face="bold", family="Times New Roman"),
        axis.text.x = element_text(size=7),        
        axis.text.y = element_text(size=7),        
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),  
        strip.background = element_rect(color = "grey", fill="black", size=0.1),
        strip.text.x = element_text(size=8,  colour="white"),
        strip.text.y = element_text(size=8,  colour="white"),
        panel.border = element_rect(colour="black", size=0.2, fill=NA),        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin = unit(0, "lines"),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8, family="Times New Roman"),        
        legend.background=element_blank(),
        legend.key=element_blank())
        
        