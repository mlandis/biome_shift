library(RevGadgets)
library(ggtree)
library(ggplot2)
library(gginnards)
library(deeptime)
library(cowplot)
source("biome_shift_util.R")

# filesystem
fp = "/Users/mlandis/projects/gh_biome_shift/"
col_fn = paste(fp, "code/plot/biome_region_colors.txt",sep="")
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/", sep="")
base_fn = c("run_1.paleo", "run_1.modern", "run_1.null")

# manage colors and state labels
dat_col = read.csv(col_fn)
st_lbl = c("Trop+SEAs","Trop+EAs","Trop+Eur","Trop+NAm","Trop+CAm","Trop+SAm",
           "Warm+SEAs","Warm+EAs","Warm+Eur","Warm+NAm","Warm+CAm","Warm+SAm",
           "Cold+SEAs","Cold+EAs","Cold+Eur","Cold+NAm","Cold+CAm","Cold+SAm")
           #"--")
names(st_lbl) = c( as.vector(dat_col$state) ) #, "--")
st_colors = c( as.vector(dat_col$color), "gray" )
names(st_colors) = st_lbl

# generate pie-state plots
summary_statistic = "PieState"
ppl = list()
for (i in 1:length(base_fn)) {

    fn = base_fn[i]
    phy_fn = paste(fn, ".ase.tre",  sep="")
    stree_fn = paste(out_fp, phy_fn, sep="")
  
  
    zz=plot_ancestral_states(tree_file=stree_fn,
                          summary_statistic=summary_statistic,
                          state_labels=st_lbl,
                          state_colors=st_colors,
                          node_label_size=0,
                          node_size_range=c(0.5,2.0),
                          node_label_nudge_x=0.5,
                          tip_node_size=0.75,
                          tip_label_size=2*0.9,
                          tip_label_offset=0.5,
                          xlim_visible=c(0,70),
                          show_posterior_legend=T,
                          node_pie_diameter=5*1.2,
                          tip_pie_diameter=4*1.2,
                          pie_nudge_x=0.3,
                          pie_nudge_y=0.4,
                          alpha=1)
    
    ppl[[i]] = zz
}

# get number of nodes to set ylim
n_node = sum(!is.na(zz$data$label))

# generate plots for the three analysis types
plot_fn = paste0( plot_fp, "fig5_anc_", c("paleo","modern","null"), ".pdf")
for (i in 1:length(plot_fn)) {

    # get simple RevGadgets plot
    p2 = ppl[[i]]
    
    # tweak coordinate system
    x_height = max(p2$data$x)
    x_new_root = 75
    x_offset = x_new_root - x_height
    x_extra = 0
    x_lbl = 8
    dy = 4
    
    # add ggplot features
    p2 = p2 + theme_tree2()
    p2 = p2 + coord_cartesian(xlim = c(-(x_offset+x_extra), x_height+x_lbl), ylim=c(0,n_node+1), expand=TRUE)
    p2 = p2 + scale_x_continuous(breaks = seq(0,x_new_root,10)-x_offset+x_new_root%%10, labels = rev(seq(0,x_new_root,10)))
    p2 = p2 + labs(x="Age (Ma)")
    p2 = p2 + theme(legend.position=c(0.0, 0.5), axis.line = element_line(colour = "black"))
    p2 = p2 + scale_colour_manual( name="Biome+Region",
                                   drop=F,
                                   values=st_colors,
                                   breaks=names(st_colors),
                                   limits=names(st_colors))
    
    # add legend
    p2 = p2 + guides( colour=guide_legend(title="Biome+Region", override.aes = list(size = 4)) )
    my_guide = cowplot::get_legend(p2)
    p2 = p2 + guides( colour=FALSE )

    # add timescale
    p2 = fig5_add_epoch_times(p2, x_new_root,  dy_bars=-7, dy_text=-3) #add_epoch_times(p2, max_age=x_new_root, x_offset=x_offset, dy=4)
    
    # plot legend beside figure
    pg = plot_grid( NULL, my_guide, p2, align="hv", nrow=1, rel_widths = c(0.75,0.75,6))
    
    # print figure to file
    pdf(plot_fn[i], height=8, width=20*(1/3), useDingbats=F)
    print(pg)
    dev.off()
}
