#!/usr/bin/Rscript
library(ggplot2)
library(scales)
library(directlabels)
library(wesanderson)
#library(gdata)
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
    stop("usage: ./peakperf.R csv-file output-pdf", call.=FALSE)
} 
data <- read.csv(file=args[1], header=TRUE, sep=",")

ggplot(data, aes(x=n,y=perf, group=type, colour=type))+
  geom_line(aes(linetype=bcast)) +
  #scale_colour_manual(values=cbPalette, guide='none')+
  ylim(0,16)+
  ylab(element_blank())+
  xlim(25,500)+
  labs(title="Performance with -03 -mavx2 -mfma [flops/cycles]")+
  theme(plot.title = element_text(size=10))+
  # peak perf and cache markings
  geom_vline(xintercept=340, linetype="dashed", color = "red")+
  # takes care of end of line labels
  geom_dl(aes(label = type), method = list(dl.combine("first.points", "last.points"), cex = 0.8))
ggsave(args[2], width=6,height=4)

