#
library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(bayestestR)
library(data.table)
library(dplyr)
library(gginnards)
library(ggpubr)
library(ggridges)

source("biome_shift_util.R")

# filesystem
fp       = "/Users/mlandis/projects/gh_biome_shift/"
plot_fp  = paste0(fp, "code/plot/")
out_fp   = paste0(fp, "output/")
atlas_fp = paste0(fp, "data/atlas/")
atlas_fn = paste0(atlas_fp, "epoch_names.txt")
col_fn   = paste0(plot_fp, "biome_region_colors.txt")
plot_fn  = paste0(plot_fp, "fig/fig8_event_epochs.pdf")
base_fn  = c("run_1.paleo", "run_1.modern", "run_1.null")

# get colors and names for biome+region states
dat_col = read.csv(col_fn, stringsAsFactors=F)
n_states = nrow(dat_col)
st_lbl = dat_col$str
st_colors = as.vector(dat_col$color)
names(st_colors) = st_lbl

# get graphs
graphs = make_graphs( atlas_fp, atlas_fn )

# get datasets
f_burn = 0
thinby = 1

# epochs as time bins
epoch_times = c(65, 56, 48, 34, 23, 16, 5.3, 0.0)

# labels for events/triplets
event_lbl = c("Warm+EAs to Cold+EAs", "Cold+EAs to Warm+EAs",
              "Warm+NAm to Cold+NAm", "Cold+NAm to Warm+NAm",
              "Warm+EAs to Warm+NAm", "Warm+NAm to Warm+EAs",
              "Cold+EAs to Cold+NAm", "Cold+NAm to Cold+EAs")

trip_lbl = c("Biome First", "Biome Flight", "Biome Reversal", "Region First", "Region Flight", "Region Reversal")

# plotting settings
linewidth = c(upper=0.5, mean=1, lower=0.5)
linetype = c(upper=6, mean=1, lower=6)
linealpha = c(upper=0.5, mean=1, lower=0.5)

# plot coord offset for event/triplet types
dx_event = (0:7-3.5)*0.1
dx_trip = (0:5-2.5)*0.1

# plot legends
my_guide_legend_1 = guide_legend(title="Biome shift or dispersal event",
                                 title.position="top", title.hjust=0.5,
                                 nrow=4, byrow = T)
my_guide1 = guides( color=my_guide_legend_1)
my_guide_legend_2 = guide_legend(title="Event series type",
                                 title.position="top", title.hjust=0.5,
                                 nrow=2, ncol=3, byrow =T)
my_guide2 = guides( color=my_guide_legend_2)



# process data
d_raw = list()
d_param = list()
d_trip = list()
d_state_ep = list()
d_trip_ep = list()
p_event = list()
p_trip = list()
for (i in 1:length(base_fn)) {
    # get output
    history_fn        = paste0(out_fp, base_fn[i],".history.tsv")
    model_fn          = paste0(out_fp, base_fn[i], ".model.log")
    # read triplet data
    d_raw[[i]]        = make_dat_triplet_raw(history_fn, n_states, f_burn, thinby )
    # read parameters
    d_param[[i]]      = read.table(model_fn, sep="\t", header=T)
    if (grepl("null", base_fn[i], fixed=T)) { d_param[[i]] = cbind( d_param[[i]], data.frame("w_adj.1."=1, "w_adj.2."=0, "w_adj.3."=0)) }
    # generate triplet table with simulated subroot states
    d_trip[[i]]       = make_triplet_df(d_raw[[i]], d_param[[i]], graphs[[1]], subroot=T)
    # build LSTT table for states + epochs
    d_state_ep[[i]]   = make_state_lstt_epoch( d_trip[[i]], epoch_times, event_lbl )
    # build LSTT table for triplets + epochs
    d_trip_ep[[i]]    = make_triplet_lstt_epoch( d_trip[[i]], epoch_times, trip_lbl )
    # adjust x for event/trip types
    d_state_ep[[i]]$x = d_state_ep[[i]]$Var1 + dx_event[d_state_ep[[i]]$Var2] - 1
    d_trip_ep[[i]]$x  = d_trip_ep[[i]]$Var1 + dx_trip[d_trip_ep[[i]]$Var2] - 1
    # make plots for events/triplets
    p_event[[i]] = make_event_plot(d_state_ep[[i]][ d_state_ep[[i]]$Var1 > 1, ])
    p_trip[[i]]  = make_trip_plot(d_trip_ep[[i]][ d_trip_ep[[i]]$Var1 > 1, ])
}

# make legends for events/triplets
leg_event = get_legend( make_event_plot(d_state_ep[[1]][ d_state_ep[[1]]$Var1 > 1, ], my_guide1) )
leg_trip  = get_legend( make_trip_plot(d_trip_ep[[1]][ d_trip_ep[[1]]$Var1 > 1, ], my_guide2) )

# combine plots
pspg = plot_grid(leg_event, leg_trip,
                 p_event[[1]], p_trip[[1]],
                 p_event[[2]], p_trip[[2]],
                 p_event[[3]], p_trip[[3]],
                 ncol=2, nrow=4, rel_heights = c(.75,2,2,2),
                 hjust=0, vjust=rep(-0.5, 8),
                 labels=c("","","(A) Paleobiome","(D) Paleobiome","(B) Modern Biome","(E) Modern Biome","(C) Null Biome","(F) Null Biome"), scale=0.9)

# save plot
CairoPDF( file=plot_fn, width=16, height=12)
print(pspg)
dev.off()

