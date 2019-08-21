#!/usr/bin/Rscript
library(ggplot2)
library(scales)
library(directlabels)
library(wesanderson)
data_run_size_dep <- read.csv(file="run-size-dependency.txt", header=TRUE, sep=",")
data_run_size_dep_mcost <- read.csv(file="run-size-dependency-mcost.txt", header=TRUE, sep=",")
data_run_size_dep_mcost$byte_per_pixel
# data per pixel in bytes for simd_v3 (only data used in mcost)
# uint8_t* i; rgba, 4 byte
# float* g: 4 byte
# 2 views
data_run_size_dep_mcost[data_run_size_dep_mcost$kernel_name == "simd_v3", "byte_per_pixel"]  = (4+4)*2
# data per pixel in bytes for simd_v17 (only data used in mcost)
# float* i_f; rgba, 16 byte
# float* g: 4 byte
# 2 views
data_run_size_dep_mcost[data_run_size_dep_mcost$kernel_name == "simd_v17", "byte_per_pixel"]  = (16+4)*2
# data per pixel in bytes for simd_v21 (only data used in mcost)
# float* ir_f; 4 byte
# float* ig_f; 4 byte
# float* ib_f; 4 byte
# float* g: 4 byte
# 2 views
data_run_size_dep_mcost[data_run_size_dep_mcost$kernel_name == "simd_v21", "byte_per_pixel"]  = (12+4)*2
L1=32*1024
L2=256*1024
L3=30*1024*1024


data_run_bar_chart <- read.csv(file="run-bar-chart-validation.txt", header=TRUE, sep=",")
data_run_bar_chart$order
data_run_bar_chart[data_run_bar_chart$kernel_name == "baseline", "order"]  = 0
data_run_bar_chart[data_run_bar_chart$kernel_name == "compressed_planes", "order"]  = 1
data_run_bar_chart[data_run_bar_chart$kernel_name == "removed_inside_and_other", "order"]  = 2
data_run_bar_chart[data_run_bar_chart$kernel_name == "split_up_cost", "order"]  = 3
data_run_bar_chart[data_run_bar_chart$kernel_name == "weight_precomputation", "order"]  = 4
data_run_bar_chart[data_run_bar_chart$kernel_name == "fast_exp", "order"]  = 5
data_run_bar_chart[data_run_bar_chart$kernel_name == "weight_precomputation_fast_exp", "order"]  = 6
data_run_bar_chart[data_run_bar_chart$kernel_name == "remove_if", "order"]  = 7
data_run_bar_chart[data_run_bar_chart$kernel_name == "inline_l1_norm", "order"]  = 8
data_run_bar_chart[data_run_bar_chart$kernel_name == "rgba", "order"]  = 9
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v0", "order"]  = 10
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v1", "order"]  = 11
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v2", "order"]  = 12
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v3", "order"]  = 13
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v4", "order"]  = 14
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v5", "order"]  = 15
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v6", "order"]  = 16
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v7", "order"]  = 17
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v8", "order"]  = 18
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v9", "order"]  = 19
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v10", "order"]  = 20
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v11", "order"]  = 21
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v12", "order"]  = 22
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v13", "order"]  = 23
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v14", "order"]  = 24
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v15", "order"]  = 25
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v16", "order"]  = 26
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v17", "order"]  = 27
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v18", "order"]  = 28
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v21", "order"]  = 29
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v22", "order"]  = 30
data_run_bar_chart[data_run_bar_chart$kernel_name == "simd_v23", "order"]  = 31

# shorten long names
data_run_bar_chart$kernel_name <- as.character(data_run_bar_chart$kernel_name)

data_run_bar_chart[data_run_bar_chart == "removed_inside_and_other"]  = "removed_inside"
data_run_bar_chart[data_run_bar_chart == "weight_precomputation"]  = "precomp"
data_run_bar_chart[data_run_bar_chart == "weight_precomputation_fast_exp"]  = "precomp_fast_exp"

# Select a subset of kernels
final_kernels <- c("baseline", "weight_precomp", "rgba", "simd_v7", "simd_v10", "simd_v17")
data_run_bar_chart_final_subset <- data_run_bar_chart[data_run_bar_chart$kernel_name %in% final_kernels,]

model_string <- "XeonE5_2680v3@2.5GHz, L1: 32K, L2: 256K, L3: 30M"

#Plot run-size dependency for normal scale 
runSizeDependencyPlot <- function(data, log, outputFile){
  p<-ggplot(data, aes(x=rows*cols/1000000,y=time/60, group=kernel_name, colour=kernel_name))+
    geom_line(aes(linetype=binary_name)) +
    ylim(0,125)+
    ylab(element_blank())+
    xlab("MegaPixels")+
    xlim(0,3.5)+
    labs(title="Time [min]")+
    theme(plot.title = element_text(size=10))+
    geom_dl(aes(label = kernel_name), method = list(dl.combine("first.points", "last.points"), cex = 0.8))
    if(log == 1){
      p = p + scale_y_continuous(trans = "log10")
      p = p + scale_x_continuous(trans = "log10")
    }
  ggsave(outputFile, width=12,height=8)
}

#Plot run-size dependency for mcost
fac=1024*1024
runSizeDependencyMcostPlot <- function(data, log, outputFile){
  p<-ggplot(data, aes(x=rows*cols*byte_per_pixel/fac,y=avg_mcost_time, group=kernel_name, colour=kernel_name))+
    geom_point()+
    geom_line() +
    xlim(-10,115)+
    geom_vline(xintercept = L3/fac, linetype="dashed",
               color = "black", size=.5)+
    annotate("text",x=L3/fac, y=13000, label="L3 cache", size=6, angle=90, vjust=-0.5, hjust=0, color="black")+
    ylab(element_blank())+
    xlab("problem size / MB used by cost function (~cache usage)")+
    labs(title="average cycles for cost function")+
    theme(legend.position="none",
          plot.title = element_text(size=18),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14))+
    geom_dl(aes(label = paste("", kernel_name, "", sep=" ")), method = list(dl.combine("first.points", "last.points"), cex = 1.3))
  if(log == 1){
    #p = p + scale_y_continuous(trans = "log10")
    p = p + scale_x_continuous(trans = "log10")
  }
  ggsave(outputFile, width=9,height=6)
}

# Bar chart ordered by descending time
runBarChart <- function(data, width, outputFile){
  ggplot(data, aes(x=reorder(kernel_name, order), y=time/60))+
    geom_bar(stat="identity", width=0.5)+
    #ylim(0,125)+
    ylab(element_blank())+
    xlab("Kernel")+
    #xlim(0,3.5)+
    labs(title="Time [min]")+
    theme(plot.title = element_text(size=10), axis.text.x= element_text(angle=55, hjust=1, vjust=1))
    #geom_dl(aes(label = kernel_name), method = list(dl.combine("first.points", "last.points"), cex = 0.8))
  ggsave(outputFile, width=width,height=4)
}
#bar chart with speedup shown
runBarChartSpeedup <- function(data, width, outputFile){
  baseline_time<- data[data$kernel_name == 'baseline',]$time
  ggplot(data, aes(x=reorder(kernel_name, order), y=baseline_time/time))+
    geom_bar(stat="identity", width=0.5)+
    #ylim(0,125)+
    ylab(element_blank())+
    xlab("Kernel")+
    #xlim(0,3.5)+
    labs(title="Speedup")+
    theme(plot.title = element_text(size=10), axis.text.x= element_text(size=8,angle=40, hjust=1, vjust=1), axis.title.x = element_text(margin = margin(t = -15)))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
    #geom_dl(aes(label = kernel_name), method = list(dl.combine("first.points", "last.points"), cex = 0.8))
  ggsave(outputFile, width=width,height=2.5)
}
#runSizeDependencyPlot(data_run_size_dep, 0, "run-size-dependency.pdf")
#runSizeDependencyPlot(data_run_size_dep, 0, "run-size-dependency.png")
#runSizeDependencyPlot(data_run_size_dep, 1, "run-size-dependency-loglog.pdf")
#runSizeDependencyPlot(data_run_size_dep, 1, "run-size-dependency-loglog.png")

runSizeDependencyMcostPlot(data_run_size_dep_mcost, 0, "run-size-dependency-mcost.pdf")
runSizeDependencyMcostPlot(data_run_size_dep_mcost, 0, "run-size-dependency-mcost.png")
runSizeDependencyMcostPlot(data_run_size_dep_mcost, 1, "run-size-dependency-mcost-loglog.pdf")
runSizeDependencyMcostPlot(data_run_size_dep_mcost, 1, "run-size-dependency-mcost-loglog.png")

runBarChart(data_run_bar_chart, 8, "run-bar-chart.pdf")
runBarChart(data_run_bar_chart, 8, "run-bar-chart.png")

runBarChart(data_run_bar_chart_final_subset, 4, "run-bar-chart-subset.pdf")
runBarChart(data_run_bar_chart_final_subset, 4, "run-bar-chart-subset.png")

#Speedup plots
runBarChartSpeedup(data_run_bar_chart, 10, "run-bar-chart-speedup.pdf")
runBarChartSpeedup(data_run_bar_chart, 10, "run-bar-chart-speedup.png")

runBarChartSpeedup(data_run_bar_chart_final_subset, 4, "run-bar-chart-speedup-subset.pdf")
runBarChartSpeedup(data_run_bar_chart_final_subset, 4, "run-bar-chart-speedup-subset.png")
