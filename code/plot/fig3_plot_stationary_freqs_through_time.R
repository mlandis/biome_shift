library(ggplot2)
library(gginnards)
source("biome_shift_util.R")

# filesystem
fp = "/Users/mlandis/projects/gh_biome_shift/"
plot_fp = paste(fp, "code/plot/", sep="")
plot_fn = paste(plot_fp, "fig3_stationary_freqs_through_time.csv", sep="")
col_fn = paste(plot_fp, "biome_region_colors.txt", sep="")

# get colors for biome+region states
col_dat = read.csv(col_fn, sep=",", stringsAsFactors=F)
st_colors = as.vector(col_dat$color)
names(st_colors) = col_dat$str

# read in stationary freq data (generated through RevBayes; script not provided)
freqs = read.csv(plot_fn, sep=",", header=T)
freqs$state = col_dat$str[ freqs$s ]
freqs$state = factor( freqs$state, ordered=T, levels=col_dat$str )

# set dimensions
tmax = 70
freqs$t0[ freqs$t0==80 ] = tmax
freqs$dt = freqs$t0 - freqs$t1
freqs$mt = freqs$t0 - freqs$dt/2

# plot state frequencies for each time interval
p = ggplot( freqs, aes(x=mt, y=pi_s, fill=state, width=dt))
p = p + geom_bar(stat="identity",position="fill")
p = p + xlab("Age (Ma)")
p = p + ylab("Lineage states\n(proportions)")
p = p + scale_x_continuous(trans="reverse", limits=c(tmax,0) )
p = p + coord_cartesian(ylim=c(-0.15,1))
p = p + scale_fill_manual(values=st_colors, guide="legend")
p = p + theme_classic()
p = p + theme(legend.position="top",
                legend.text = element_text(size=8),
                legend.key.size = unit(0.1, "cm"))
p = p + guides( fill=guide_legend(title="Biome+Region",
                                    title.position="top",
                                    title.hjust=0.5, nrow=3, ncol=6, byrow =T,
                                    override.aes = list(size=3)) )

# add epoch info to timescale
p = fig3_add_epoch_times(p, max_age=72, dy=0.1) 

# save file
pdf(  file=paste(plot_fp, "fig/fig3_stationary_freqs_through_time.pdf", sep=""), width=8, height=4)
print(p)
dev.off()
