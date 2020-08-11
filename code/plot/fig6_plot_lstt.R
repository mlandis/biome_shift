library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(gginnards)
library(data.table)

source("biome_shift_util.R")

# filesystem
fp = "/Users/mlandis/projects/gh_biome_shift/"
plot_fp = paste(fp, "code/plot/fig/", sep="")
col_fn = paste(fp, "code/plot/biome_region_colors.txt",sep="")
atlas_fp = paste(fp, "data/atlas/", sep="")
atlas_fn = c("tropical", "warm", "cold")
times_fn = paste(atlas_fp, "epoch_names.txt", sep="")

# MCMC processing
f_burn = 0.5
thinby = 1

# get colors and names for biome+region states
dat_col = read.csv(col_fn, stringsAsFactors=F) # color data
n_states = nrow(dat_col)
st_lbl = c("Trop+SEAs","Trop+EAs","Trop+Eur","Trop+NAm","Trop+CAm","Trop+SAm",
           "Warm+SEAs","Warm+EAs","Warm+Eur","Warm+NAm","Warm+CAm","Warm+SAm",
           "Cold+SEAs","Cold+EAs","Cold+Eur","Cold+NAm","Cold+CAm","Cold+SAm")
st_colors = as.vector(dat_col$color)
names(st_colors) = st_lbl

# get epoch times
epoch_times = read.csv(times_fn, sep=",", header=T)
epoch_times[,c("start_age","end_age")] = round(epoch_times[,c("start_age","end_age")]) 
epoch_times$index = epoch_times$index + 1 # base-1 indexing

# get all biome graphs
n_atlas = length(atlas_fn)
atlas_list = list()
for (i in 1:n_atlas) {
    atlas_list[[i]] = list()
    afn = paste(atlas_fp, atlas_fn[i], sep="")
    atlas_files = list.files(afn, full.names = T )
    for (j in 1:length(atlas_files)) {
        gij = read.csv( atlas_files[[j]], sep=",", header=F)
        atlas_list[[i]][[j]] = gij
    }
}

# process files
base_fn1 = "run_1.paleo"
base_fn2 = "run_1.modern"
base_fn3 = "run_1.null"

d1 = make_lstt_dat( paste0(fp, "output/", base_fn1, ".history.tsv"), n_states, f_burn, thinby)
d2 = make_lstt_dat( paste0(fp, "output/", base_fn2, ".history.tsv"), n_states, f_burn, thinby)
d3 = make_lstt_dat( paste0(fp, "output/", base_fn3, ".history.tsv"), n_states, f_burn, thinby)

m1 = make_match_bins( d1, epoch_times, atlas_list)
m2 = make_match_bins( d2, epoch_times, atlas_list)
m3 = make_match_bins( d3, epoch_times, atlas_list)

age_oligo = 34
mf1_lt34 = 1 - sum(m1[ m1$match=="match"&m1$age<=age_oligo, ]$count) / sum(m1[m1$age<=age_oligo, ]$count)
mf2_lt34 = 1 - sum(m2[ m2$match=="match"&m2$age<=age_oligo, ]$count) / sum(m2[m2$age<=age_oligo, ]$count)
mf3_lt34 = 1 - sum(m3[ m3$match=="match"&m3$age<=age_oligo, ]$count) / sum(m3[m3$age<=age_oligo, ]$count)

age_oligo = 34
mf1_gt34 = 1 - sum(m1[ m1$match=="match"&m1$age>age_oligo, ]$count) / sum(m1[m1$age>age_oligo, ]$count)
mf2_gt34 = 1 - sum(m2[ m2$match=="match"&m2$age>age_oligo, ]$count) / sum(m2[m2$age>age_oligo, ]$count)
mf3_gt34 = 1 - sum(m3[ m3$match=="match"&m3$age>age_oligo, ]$count) / sum(m3[m3$age>age_oligo, ]$count)

age_oligo = 0
mf1_all = 1 - sum(m1[ m1$match=="match"&m1$age>age_oligo, ]$count) / sum(m1[m1$age>age_oligo, ]$count)
mf2_all = 1 - sum(m2[ m2$match=="match"&m2$age>age_oligo, ]$count) / sum(m2[m2$age>age_oligo, ]$count)
mf3_all = 1 - sum(m3[ m3$match=="match"&m3$age>age_oligo, ]$count) / sum(m3[m3$age>age_oligo, ]$count)

print("Treewide mismatch")
c(P=mf1_all, M=mf2_all,  N=mf3_all) #,  mf4_lt34)

print("Before Oligocene mismatch")
c(P=mf1_gt34, M=mf2_gt34,  N=mf3_gt34) #

print("After Oligocene mismatch")
c(P=mf1_lt34, M=mf2_lt34,  N=mf3_lt34) #

# age offset
d1$age2 = d1$age+0.5
d2$age2 = d2$age+0.5
d3$age2 = d3$age+0.5

m1$age2 = m1$age+0.5
m2$age2 = m2$age+0.5
m3$age2 = m3$age+0.5

# Lineage-states through time
my_guide = guides( fill=guide_legend(title="Biome+Region",
                                    title.position="top",
                                    title.hjust=0.5, nrow=3, ncol=6, byrow =T,
                                    override.aes = list(size=3)) )

#############
#############
# first plot is just the counts

# LSTT, dataset 1
p1s = ggplot(data = d1, aes(fill=State, x=age2, y=prop))
p1s = p1s + geom_bar(stat="identity", position="fill") # position_nudge(x = 0.5))
p1s = p1s + xlab("") #Age (Ma)")
p1s = p1s + ylab("Lineage states\n(proportions)")
p1s = p1s + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
#p1s = p1s + ylim(-0.5, 1)
p1s = p1s + coord_cartesian(ylim=c(-0.15,1))
p1s = p1s + scale_fill_manual(values=st_colors, guide=FALSE) #"legend")
p1s = p1s + theme_classic()
p1s = p1s + theme(legend.position="top",
                legend.justification="center",
                legend.text = element_text(size=8),
                legend.key.size = unit(0.1, "cm"))
p1s = p1s + my_guide
p1s = fig6_add_epoch_times(p1s, max_age = 72)
#p1s


# LSTT, dataset 2
p2s = ggplot(data = d2, aes(fill=State, x=age2, y=prop))
p2s = p2s + geom_bar(stat="identity", position="fill")
p2s = p2s + xlab("") #Age (Ma)")
p2s = p2s + ylab("Lineage states\n(proportions)")
p2s = p2s + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
p2s = p2s + coord_cartesian(ylim=c(-0.15,1))
#p3 = p3 + scale_y_continuous( breaks=seq(-10,10,2), limits=c(-5,5) )
p2s = p2s + scale_fill_manual(values=st_colors, guide="legend")
p2s = p2s + theme_classic()
p2s = p2s + theme(legend.position="top",
                legend.text = element_text(size=8),
                legend.key.size = unit(0.1, "cm"))
p2s = p2s + guides( fill=FALSE)
p2s = fig6_add_epoch_times(p2s, max_age = 72)


p3s = ggplot(data = d3, aes(fill=State, x=age2, y=prop))
p3s = p3s + geom_bar(stat="identity", position="fill")
p3s = p3s + xlab("Age (Ma)")
p3s = p3s + ylab("Lineage states\n(proportions)")
p3s = p3s + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
p3s = p3s + coord_cartesian(ylim=c(-0.15,1))
#p3 = p3 + scale_y_continuous( breaks=seq(-10,10,2), limits=c(-5,5) )
p3s = p3s + scale_fill_manual(values=st_colors, guide="legend")
p3s = p3s + theme_classic()
p3s = p3s + theme(legend.position="top",
                legend.text = element_text(size=8),
                legend.key.size = unit(0.1, "cm"))
p3s = p3s + guides( fill=FALSE)
p3s = fig6_add_epoch_times(p3s, max_age = 72)


# Lineage biome-matching through time

match_cols = c(mismatch="#BBBBBB",match="#444444")
my_guide_match = guides( fill=guide_legend(title="Lineage state match?",
                                    title.position="top",
                                    title.hjust=0, nrow=3, byrow=T,
                                    override.aes = list(size=3)) )


# LMTT, dataset 1
p1m = ggplot(data = m1, aes(fill=match, x=age2, y=prop))
p1m = p1m + geom_bar(stat="identity", position="fill")
p1m = p1m + xlab("")
p1m = p1m + ylab("Lineage biome matches\n(proportions)")
p1m = p1m + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
p1m = p1m + coord_cartesian(ylim=c(-0.15,1))
p1m = p1m + scale_fill_manual(values=match_cols, guide="legend", labels=c("lineage's region lacks biome", "lineage's region contains biome"))
p1m = p1m + theme_classic()
p1m = p1m + theme(legend.position="top",
                legend.justification ="center",
                legend.text = element_text(size=10),
                legend.key.size = unit(0.1, "cm"))
p1m = p1m + my_guide_match
p1m = fig6_add_epoch_times(p1m, max_age = 72)

# LMTT, dataset 2

p2m = ggplot(data = m2, aes(fill=match, x=age2, y=prop))
p2m = p2m + geom_bar(stat="identity", position="fill")
p2m = p2m + xlab("")
p2m = p2m + ylab("Lineage biome matches\n(proportions)")
p2m = p2m + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
p2m = p2m + coord_cartesian(ylim=c(-0.15,1))
p2m = p2m + scale_fill_manual(values=match_cols, guide="legend")
p2m = p2m + theme_classic()
p2m = p2m + guides( fill=FALSE)
p2m = fig6_add_epoch_times(p2m, max_age = 72)

p3m = ggplot(data = m3, aes(fill=match, x=age2, y=prop))
p3m = p3m + geom_bar(stat="identity", position="fill")
p3m = p3m + xlab("Age (Ma)")
p3m = p3m + ylab("Lineage biome matches\n(proportions)")
p3m = p3m + scale_x_continuous(trans="reverse", limits=c(70,0), breaks=seq(0,70,10) )
p3m = p3m + coord_cartesian(ylim=c(-0.15,1))
#p3 = p3 + scale_y_continuous( breaks=seq(-10,10,2), limits=c(-5,5) )
p3m = p3m + scale_fill_manual(values=match_cols, guide="legend")
p3m = p3m + theme_classic()
p3m = p3m + guides( fill=FALSE)
p3m = fig6_add_epoch_times(p3m, max_age = 72)

ppg = plot_grid(p1s,p1m,p2s,p2m,p3s,p3m,
                ncol=2, nrow=3, rel_heights = c(2.55,2,2),
                hjust=0, vjust=c(7.5,7.5,-1,-1,-1,-1)+0.5,
                labels=c("(A) Paleobiome","(D) Paleobiome","(B) Modern Biome","(E) Modern Biome","(C) Null Biome","(F) Null Biome"),
                label_size=12)

CairoPDF( file=paste(plot_fp, "fig6_lstt_states.pdf", sep=""), width=16, height=12)
print(ppg)
dev.off()
