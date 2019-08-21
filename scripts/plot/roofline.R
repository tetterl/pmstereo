#!/usr/bin/Rscript
library(ggplot2)
library(MASS)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
    stop("usage: ./roofline.R outputpdf", call.=FALSE)
}
#library(gdata)
set.seed(16)
# ridge point
ridge_x <- 2
ridge_y <- 2
#origin of 45 degree angled line
orig_x <- ridge_x-ridge_y
# plot limits
x_lim <- 16
y_lim <- 16

# plot without data input
ggplot() +
    geom_point() +
    #xlim(0, x_lim)+
    #ylim(0,y_lim)+
    labs(title="Performance P=W/T [ops/cycle]")+
    geom_segment(aes(x = orig_x, xend = ridge_x, y = 0, yend = ridge_y)) +
    geom_segment(aes(x = ridge_x, xend = x_lim, y = ridge_y, yend = ridge_y)) +
    ylab("")+
    coord_fixed()+ # such that x/y = 1
    # Note: frations is a method in the MASS package that converts to a fraction
    scale_x_continuous(trans = 'log2', limits = c(1/8,x_lim), labels=fractions)+
    scale_y_continuous(trans = 'log2', limits = c(1/16,y_lim), labels=fractions)+
    xlab("Operational Intensity I=W/Q [ops/bytes]") 

ggsave(args[1], width=6,height=4)
