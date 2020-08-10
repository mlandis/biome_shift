library(ggplot2)
library(cowplot)
library(reshape2)
library(patchwork)
source("biome_shift_util.R")


# filesystem
fp = "/Users/mlandis/projects/gh_biome_shift/"
atlas_fp = paste(fp, "data/atlas/", sep="")
atlas_fn = paste(atlas_fp, "epoch_names.txt", sep="")
out_fn = paste(fp, "code/plot/fig/fig2_paleobiome_maps.pdf", sep="")

# area names
area_names = c("SE As", "E As", "Eur", "N Am", "C Am", "S Am")

# read in names of epochs
atlas = read.table(atlas_fn,header=T,sep=",")

# biome graph prefixes
graph = c("null", "land", "tropical", "warm", "cold")
graph2 = c("Null", "Land", "Tropical", "Warm Temperate", "Cold Temperate")
graph_expr = c( expression(atop("Uninformative", paste("A" [U]))),
                expression(atop("Geographical", paste("A" [G], "(m)", sep=""))),
                expression(atop("Tropical", paste("A" [T], "(m)", sep=""))),
                expression(atop("Warm Temperate", paste("A" [W], "(m)", sep=""))),
                expression(atop("Cold Temperate", paste("A" [C], "(m)", sep=""))))
n_graphs = length(graph)

# colors associated with Null, Land, Tropical, Warm Temp., Cold Temp. graphs
graph_color = c("black","burlywood4","firebrick3","forestgreen","dodgerblue")

# construct data tables for each biome (i) and each timeslice (j)
dat = list()
for (i in 1:n_graphs) {
    dat[[ i ]] = list()
     
    files = list.files(path=paste(atlas_fp, graph[i], "/", sep=""), full.names = T)
    n_files = length(files)
    for (j in 1:n_files) {
        tmp = read.table(files[j], sep=",")
        rownames(tmp) = area_names
        colnames(tmp) = area_names
        dat[[i]][[j]] = melt(as.matrix(tmp))
        colnames(dat[[i]][[j]])=c("From","To","Area")
    }
}

# plot everything
idx = 1
plotlist=list()
for (j in 1:n_files) {
    for (i in 1:n_graphs) {
        
        # plot column title if first row
        col_title = " "
        if (j == 1) {
            col_title = graph_expr[i]
        }
        
        # plot row title if first column
        row_title = " "
        if (i == 1) {
            start_age = sprintf(atlas$start_age[j], fmt = '%#.1f') 
            end_age = sprintf(atlas$end_age[j], fmt = '%#.1f') 
            row_title = paste(atlas$name[j], "\n", start_age, " to ", end_age, " Ma\n", "m = ", j, sep="")
        }
        
        # generate matrix plot for biome (i) and time (j)
        p = ggplot(dat[[i]][[j]], aes(From,To))
        p = p + geom_tile(aes(fill=Area), colour="white", size=0.5)
        p = p + scale_y_discrete(limits = rev(levels(dat[[i]][[j]]$From)))
        p = p + scale_fill_gradientn(limits = c(-0.25,2),
                                     colours=c("white",graph_color[i]),
                                     breaks=c(0,1,2), labels=c("Marg.","Subdom.","Dom."))
        p = p + xlab(col_title) + ylab(row_title)
        p = p + scale_x_discrete(position = "top") 
        p = p + theme(legend.position="none", 
                      axis.title.y = element_text( size=18 ),
                      axis.title.x = element_text( size=18 ),
                      axis.text.y = element_text(size=14),
                      axis.text.x = element_text(size=14, angle=90, hjust=0),
                      axis.text.x.top = element_text(vjust = 0.5),
                      panel.grid.major = element_blank(),
                      panel.grid.minor =  element_blank(),
                      panel.background = element_blank(),
                      panel.border =  element_blank(),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
         
        plotlist[[ idx ]] = p
        idx = idx + 1
    }
}

# aggregate all matrices into single plot object
pg = wrap_plots(plotlist, ncol=n_graphs, nrow=n_files, byrow=T)

# save to file
save_plot(out_fn, pg, nrow=n_files, ncol=n_graphs, base_asp=1.05)
