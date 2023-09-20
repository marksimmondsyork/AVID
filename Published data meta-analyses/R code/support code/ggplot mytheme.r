# my defult theme for plotting with ggplot

library(ggplot2)

theme_set(theme_bw(base_size=16))
theme_update(panel.grid.major=element_line(colour="#CCCCCC"),panel.border=element_rect(colour="black",fill=NA))
