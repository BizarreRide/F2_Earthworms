########
# ggplot personalized Theme
# Quentin Schorpp
# 16.04.2015
########

library(ggplot2)
library(grid)

mytheme = 
  theme_bw() + 
  theme(strip.background = element_rect(color = "grey", fill="black", size=0.1),
        strip.text.x = element_text(size=8,  colour="white", face="italic"),
        strip.text.y = element_text(size=8,  colour="white", face="italic"),
        axis.text.x = element_text(size=7),
        axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=8, family="Times New Roman"),
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),
        plot.title = element_text(size=11,face="bold", family="Times New Roman"),
        panel.margin = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=0.2, fill=NA),
        legend.key=element_blank(),
        legend.background=element_blank(),
        legend.text=element_text(size=8,face="italic", family="Times New Roman"),
        legend.title=element_text(size=8))