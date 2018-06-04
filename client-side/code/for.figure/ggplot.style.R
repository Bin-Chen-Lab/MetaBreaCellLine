### A customerized ggplot style ###
require(ggplot2)

ggplot.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=55, face="bold"),
                                                 axis.text  = element_text( size=55, face="bold"),
                                                 plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                                 axis.line.x = element_line(colour = "black",size = 3),
                                                 axis.line.y = element_line(colour = "black",size = 3)
                                                 ) 