library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(gginnards)
library(data.table)

source("biome_shift_util.R")

# filesystem
fp = "~/projects/gh_biome_shift/"
plot_fp = paste(fp, "code/plot/fig/", sep="")
col_fn = paste(fp, "code/plot/biome_region_colors.txt",sep="")
atlas_fp = paste(fp, "data/atlas/", sep="")
atlas_fn = c("tropical", "warm", "cold")
times_fn = paste(atlas_fp, "epoch_names.txt", sep="")

# MCMC processing
f_burn = 0 # how many iterations to discard
thinby = 1 # sampling frequency (e.g. thinby=10 is 1 sample every 10 iterations)

# get colors and names for biome+region states
dat_col = read.csv(col_fn, stringsAsFactors=F) # color data
n_states = nrow(dat_col)
st_lbl = c("Trop+SEAs","Trop+EAs","Trop+Eur","Trop+NAm","Trop+CAm","Trop+SAm",
           "Warm+SEAs","Warm+EAs","Warm+Eur","Warm+NAm","Warm+CAm","Warm+SAm",
           "Cold+SEAs","Cold+EAs","Cold+Eur","Cold+NAm","Cold+CAm","Cold+SAm")
st_colors = as.vector(dat_col$color)
names(st_colors) = st_lbl
match_cols = c(mismatch="#BBBBBB",match="#444444")

# get epoch times
epoch_times = read.csv(times_fn, sep=",", header=T)
epoch_times[,c("start_age","end_age")] = round(epoch_times[,c("start_age","end_age")]) 
epoch_times$index = epoch_times$index + 1 # base-1 indexing

# get all biome graphs
n_atlas = length(atlas_fn)
atlas_list = list()
for (i in 1:n_atlas) {
    atlas_list[[i]] = list()
    afn = paste0(atlas_fp, atlas_fn[i])
    atlas_files = list.files(afn, full.names = T )
    for (j in 1:length(atlas_files)) {
        gij = read.csv( atlas_files[[j]], sep=",", header=F)
        atlas_list[[i]][[j]] = gij
    }
}

# create plots for three analysis types (paleo, modern, null) and
# for two presentation times (LSTT states and LSTT biome-matches)
base_fn = c("run_1.paleo", "run_1.modern", "run_1.null")

d_state = list()
d_match = list()
pp_state = list()
pp_match = list()

for (i in 1:length(base_fn)) {
    
    # filename
    fn = paste0(fp, "output/", base_fn[i], ".history.tsv")
    
    # gather state LSTT and biome-match LSTT data
    d_state[[i]] = make_lstt_dat(fn, n_states, f_burn, thinby)
    d_match[[i]] = make_match_bins( d_state[[i]], epoch_times, atlas_list)
    
    # biome+region LSTT plot
    p1s = ggplot(data = d_state[[i]], aes(fill=State, x=age2, y=prop))
    p1s = p1s + geom_bar(stat="identity", position="fill")
    p1s = p1s + ylab("Lineage states\n(proportions)")
    p1s = p1s + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
    p1s = p1s + coord_cartesian(ylim=c(-0.15,1))
    p1s = p1s + scale_fill_manual(values=st_colors, guide=FALSE) #"legend")
    p1s = p1s + theme_classic()
    p1s = p1s + theme(legend.position="top",
                    legend.justification="center",
                    legend.text = element_text(size=8),
                    legend.key.size = unit(0.1, "cm"))
    p1s = fig6_add_epoch_times(p1s, max_age = 72)
    if (i == 1) {
      p1s = p1s + guides( fill=guide_legend(title="Biome+Region",
                                    title.position="top",
                                    title.hjust=0.5, nrow=3, ncol=6, byrow =T,
                                    override.aes = list(size=3)) )
    }
    if (i == 3) {
      p1s = p1s + xlab("Age (Ma)")
    }
    else {
      p1s = p1s + xlab("")
    }
    pp_state[[i]] = p1s
    
    # biome match LSTT plot
    p1m = ggplot(data = d_match[[i]], aes(fill=match, x=age2, y=prop))
    p1m = p1m + geom_bar(stat="identity", position="fill")
    p1m = p1m + xlab("")
    p1m = p1m + ylab("Lineage biome matches\n(proportions)")
    p1m = p1m + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
    p1m = p1m + coord_cartesian(ylim=c(-0.15,1))
    p1m = p1m + scale_fill_manual(values=match_cols, guide=FALSE, labels=c("lineage's region lacks biome", "lineage's region contains biome"))
    p1m = p1m + theme_classic()
    p1m = p1m + theme(legend.position="top",
                legend.justification ="center",
                legend.text = element_text(size=10),
                legend.key.size = unit(0.1, "cm"))
    p1m = fig6_add_epoch_times(p1m, max_age = 72)
    if (i == 1) {
      p1m = p1m + guides( fill=guide_legend(title="Lineage state match?",
                                    title.position="top",
                                    title.hjust=0, nrow=3, byrow=T,
                                    override.aes = list(size=3)) )
    }
    if (i == 3) {
      p1m = p1m + xlab("Age (Ma)")
    }
    else {
      p1m = p1m + xlab("")
    }
    
    pp_match[[i]] = p1m
}

# combine subfigures into one plot
ppg = plot_grid(pp_state[[1]], pp_match[[1]],
                pp_state[[2]], pp_match[[2]],
                pp_state[[3]], pp_match[[3]],
                ncol=2, nrow=3, rel_heights = c(2.55,2,2),
                hjust=0, vjust=c(7.5,7.5,-1,-1,-1,-1)+0.5,
                labels=c("(A) Paleobiome","(D) Paleobiome","(B) Modern Biome","(E) Modern Biome","(C) Null Biome","(F) Null Biome"),
                label_size=12)

CairoPDF( file=paste(plot_fp, "fig6_lstt_states.pdf", sep=""), width=16, height=12)
print(ppg)
dev.off()

