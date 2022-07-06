#!/usr/bin/env Rscript

.libPaths(R.home("library"))

# load libraries
library(optparse)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)
require(cowplot)
require(grid)
library(RColorBrewer)
library(dplyr)
library(svglite)
library(gridExtra)      # provides side-by-side plotting
library(zoo)
#library(ggpubr)
#library(gplots)


# function to read the input file
CHARACTER_command_args <- commandArgs(trailingOnly=TRUE)
#CHARACTER_command_args <- "/Users/jjuma/Work/Sam_Oyola/AMR_PT/Camp-005-output/bedtools/coverage/Camp.coverage.gz"
gzipped = gzfile(CHARACTER_command_args[1], 'rt')
data <- read.csv(gzipped, header=F, sep = "\t")
data <- data %>% dplyr::rename(chrom="V1", locus="V2", coverage="V3")
data$locus <- as.numeric(data$locus)
data$chrom <- as.character(data$chrom)
data$chrom <- as.factor(data$chrom)


for (i in unique(data$chrom)){
  sample <- basename(paste0(substr(CHARACTER_command_args[1], 1, nchar(CHARACTER_command_args[1])-12)))
  d <- data[data$chrom == i,]
  windowed <- d %>% 
    group_by(chrom) %>% 
    do(
      data.frame(
        start = rollapply(.$locus, width=200, by=200, FUN=min, align="left"),
        end = rollapply(.$locus, width=200, by=200, FUN=max, align="left"),
        coverage = rollapply(.$coverage, width=200, by=200, FUN=median, align="left")
      )
    )
  windowed$coverage <- windowed$coverage + 1
  
  a <- ggplot(data = windowed, aes(x=end,y=coverage)) +
    geom_ribbon(aes(ymin=0, ymax=coverage), fill="#D55E00", data=) +
    theme_bw() +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(trans=log10_trans(),
                       breaks=10^c(0:10),
                       labels=trans_format('log10', math_format(10^.x)),
                       expand=c(0, 0)) +
    expand_limits(y=1) +
    ylab(bquote('log'[10]~'(Coverage+1)')) +
    xlab('Position (bp)') +
    ggtitle('coverage')
  
  # if (max(d$locus) > 1000000){
  #   a <- ggplot(data = d, aes(x=locus, y=depth, colour=chrom)) +
  #     geom_bar(stat="identity", width=0.6) +
  #     theme_bw() +
  #     scale_x_continuous(breaks=seq(0, max(d$locus), 400000)) +
  #     scale_y_continuous(trans=log10_trans(),
  #                        breaks=10^c(0:10),
  #                        labels=trans_format('log10', math_format(10^.x)),
  #                        expand=c(0, 0)) +
  #     expand_limits(y=1) +
  #     labs(title = paste(sample,'coverage'),
  #          x="Position in bp",
  #          y=bquote('log'[10]~'(Coverage+1)'),
  #          fill="Depth of coverage (log10)") +
  #     guides(fill=FALSE, color=FALSE)
  # } else {
  #   a <- ggplot(data = d, aes(x=locus, y=depth, colour=chrom)) +
  #     geom_bar(stat="identity", width=0.6) +
  #     theme_bw() +
  #     scale_x_continuous(breaks=seq(0, max(d$locus), 400)) +
  #     scale_y_continuous(trans=log10_trans(),
  #                        breaks=10^c(0:10),
  #                        labels=trans_format('log10', math_format(10^.x)),
  #                        expand=c(0, 0)) +
  #     expand_limits(y=1) +
  #     labs(title = paste(sample,'coverage'),
  #          x="Position in bp",
  #          y=bquote('log'[10]~'(Coverage+1)'),
  #          fill="Depth of coverage (log10)") +
  #     guides(fill=FALSE, color=FALSE)
  #   
  # }
  
  
  # save the plot as TIFF format (or any other)
   ggsave(file=paste0(substr(CHARACTER_command_args[1],1,nchar(CHARACTER_command_args[1])-12), ".chr.", i, ".pdf"), 
          a, width = 16.8, height = 8.4,
           units = "in", device = "pdf")
}



