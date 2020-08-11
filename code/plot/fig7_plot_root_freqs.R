library(HDInterval)
library(ggplot2)
library(reshape2)
library(ggridges)
library(dplyr)
library(bayestestR)

source("biome_shift_util.R")

# filesystem
fp = '/Users/mlandis/projects/gh_biome_shift/'
out_fp = paste(fp, 'output/', sep='')
plot_fp = paste(fp, "code/plot/fig/", sep="")
plot_fn = paste(plot_fp, "fig7_root_freqs.pdf", sep="")
col_fn = paste(fp, "code/plot/biome_region_colors.txt",sep="")
fn = paste(out_fp, c('run_1.paleo.model.log', 'run_1.modern.model.log', 'run_1.null.model.log'), sep='')

# get colors and names for biome+region states
dat_col = read.csv(col_fn, stringsAsFactors=F) # color data
st_lbl = c("Trop+SEAs","Trop+EAs","Trop+Eur","Trop+NAm","Trop+CAm","Trop+SAm",
           "Warm+SEAs","Warm+EAs","Warm+Eur","Warm+NAm","Warm+CAm","Warm+SAm",
           "Cold+SEAs","Cold+EAs","Cold+Eur","Cold+NAm","Cold+CAm","Cold+SAm")
st_colors = as.vector(dat_col$color)
names(st_colors) = st_lbl
st_shape = c( rep(22, 6), rep(21, 6), rep(24, 6) )
names(st_shape) = st_lbl
n_states = length(st_lbl)

# process files
nn = c("Paleo", "Modern", "Null")
df0 = data.frame(rf=NULL, biome=NULL, region=NULL, prob=NULL)
x = list()
for (i in 1:length(fn)) {
    
    xtmp = read.csv(fn[i], sep='\t', stringsAsFactors=F)
    x[[ nn[i] ]] = xtmp
    
    for (j in 1:18) {
        strtok = strsplit( st_lbl[j], split="\\+" )[[1]]
        rfj = paste("rf_simplex.",j,".",sep="")
        xtmp[[rfj]] = sort( xtmp[[rfj]] )
        hpd95 = hdi(xtmp[[rfj]], ci=0.95)
        hpd90 = hdi(xtmp[[rfj]], ci=0.80)
        df1 = data.frame(Model=nn[i], Biome=strtok[1], Region=strtok[2], State=st_lbl[j],
                         Mean=mean(xtmp[[rfj]]),
                         lower95=hpd95$CI_low, upper95=hpd95$CI_high,
                         lower90=hpd90$CI_low, upper90=hpd90$CI_high)
        #if (df1$upper90 > df1$upper95) {
        #    print(i)
        #    hist(xtmp[[rfj]])
        #}
        df0 = rbind(df0, df1)
    }
}

m = df0
m$State = factor(m$State, ordered=T, levels=rev(st_lbl))   
m$Model = factor(m$Model, ordered=T, levels=c("Null","Modern","Paleo"))   
m$y = c( rev(sort(rep(3:1,18))) + ((rep(18:1,3)-9.5)/18)*0.7 )
                 
p = ggplot(m)
p = p + geom_vline(xintercept = 1/18, linetype=2, color="gray")

p = p + geom_segment(data=m, mapping=aes(x=lower90, xend=upper90, y=y, yend=y, color=State), size=1.25, alpha=0.5)
p = p + geom_segment(data=m, mapping=aes(x=lower95, xend=upper95, y=y, yend=y, color=State), size=0.65, alpha=0.5)
#p = p + geom_line(data=m,  mapping=aes(x=x, y=value, color=Var2), linetype=2)
p = p + geom_point(data=m, mapping=aes(x=Mean, y=y, color=State),size=2)
p = p + geom_point(data=m, mapping=aes(x=Mean, y=y),size=0.5, color="white")

p = p + ylab("Biome structure") 
p = p + xlab(expression(paste("Posterior root stationary probability, ", pi,"(", italic(m)[root] ,')',sep="")))
p = p + scale_color_manual( name="Biome+Region", values=st_colors, breaks=names(st_colors) )
p = p + scale_shape_manual( name="Biome+Region", values=st_shape, breaks=names(st_colors) )
p = p + guides(shape = guide_legend(override.aes = list(size = 0.5)))
p = p + xlim(0.0,0.175)
p = p + scale_y_continuous( breaks=c(1,2,3), labels=c("Null","Modern","Paleo") )
p = p + theme_classic()
p = p + theme(axis.text.y = element_text(angle=90, hjust=0.5, size=10),
              legend.position = "top",
              legend.key.size = unit(0, 'lines'))

my_guide_legend = guide_legend(title="Biome+Region",
                                 title.position="top", title.hjust=0.5,
                                 nrow=3, ncol=6, byrow =T)

p = p + guides( color=my_guide_legend) #, fill=my_guide_legend_1) #my_guide_legend_1)


pdf(plot_fn, height=8, width=6)
print(p)
dev.off()
