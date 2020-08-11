library(HDInterval)
library(reshape2)
library(ggplot2)
library(data.table)
library(kdensity)
library(cowplot)
library(Cairo)
library(gginnards)
library(dplyr)

source("biome_shift_util.R")

# filesystem
fp = '/Users/mlandis/projects/gh_biome_shift/'
out_fp = paste(fp, 'output/sim/', sep='')
plot_fn = paste(fp, 'code/plot/fig/fig4_sim_results.pdf', sep='')

# get true parameter values
true_fn = paste(fp, 'output/sim/true_params.txt', sep='')
true_params = read.table( true_fn, sep=',', header=T)

# experiment settings
sim_name = c('sim.vstrong.paleo',
             'sim.strong.paleo',
             'sim.med.paleo',
             'sim.weak.paleo',
             'sim.null')
model_name = c("Paleo","Paleo","Paleo","Paleo","Null")
strength_name = c("V.Strong","Strong","Moderate","Weak","Null")
param_name = c('w_biome', 'w_land', 'w_null', 'w_not_biome', 'y_subdominant')

n_sim = length(sim_name)
n_rep = 1
f_burn = 0.10

# collect results
m_list = list()
m_idx = 1
# for each simulation setting
for (i in 1:n_sim) {
    # for each replicate per sim. setting
    for (j in 1:n_rep) {
        # get posterior
        fn_ij = paste( out_fp, sim_name[i], '.', j, '.paleo.model.log', sep='')
        cat("Processing", fn_ij, "\n")
        x = read.csv(fn_ij, sep='\t')
        
        # extract post-burnin samples for parameters of interest
        n_burn = max(1, f_burn * nrow(x))
        x = x[n_burn:nrow(x), param_name]
        
        # get posterior settings
        mdl = model_name[i]
        str = strength_name[i]
        
        # summarize posterior in various ways
        print("a")
        avg = apply(x, 2, median)
        print("b")
        hpd = apply(x, 2, hdi)
        print("c")
        cover = is_covered(hpd, true_params, mdl, str)
        print("d")
        bf = apply(x, 2, compute_bf)
        print("e")
        y = t(rbind(avg, hpd))
        print("f")
        
        # reformat posterior stats
        m = melt(y)
        names(m) = c("Parameter","Statistic","Value")
        m$Model = mdl
        m$Strength = str
        m$Replicate = j
        m$Covered = cover
        m$BayesFactor = bf
        
        # store stats
        m_list[[m_idx]] = m
        m_idx = m_idx + 1
    }
}

# format label names for aggregated results
m = rbindlist(m_list)
m$Strength[ m$Strength=="Medium" ] = "Moderate" # change name
m$Condition = m$Strength
m$Strength =  factor( m$Strength, ordered=T, levels=c("V.Strong","Strong","Moderate","Weak","Null"))
m$Model = factor( m$Model, ordered=T, levels=c("Paleo","Null"))
m$Condition=factor( m$Condition,
                    ordered=T,
                    levels=c("Null","Weak","Moderate","Strong","V.Strong"))

# format label names for true parameter values
m_true = melt(true_params)
colnames(m_true) = c("Model","Strength","Parameter","Value")
m_true$Condition = m_true$Strength
m_true = m_true[ m_true$Condition %in% levels(m$Condition), ]
m_true$Strength =  factor( m_true$Strength, ordered=T, levels=c("V.Strong","Strong","Moderate","Weak","Null"))
m_true$Model = factor( m_true$Model, ordered=T, levels=c("Paleo","Null"))
m_true$Condition=factor( m_true$Condition,
                    ordered=T,
                    levels=c("Null","Weak","Moderate","Strong","V.Strong"))
m_true$Statistic = "True"

# classify Bayes Factors
m$BayesFactorCategory = make_bf_cat( m$BayesFactor )
m$BayesFactorCategory = factor(m$BayesFactorCategory, ordered=T, levels=rev(levels(m$BayesFactorCategory)))
m$Statistic = as.vector(m$Statistic)
m$Statistic[ m$Statistic == "avg" ] = "Median"
m$Statistic[ m$Statistic == "upper" ] = "Upper (HPD95)"
m$Statistic[ m$Statistic == "lower" ] = "Lower (HPD95)"
m$Statistic = sapply( m$Statistic, tools::toTitleCase)
m$Statistic = factor(m$Statistic, ordered=T, levels=c("Lower (HPD95)","Median","Grand median","Upper (HPD95)", "True"))

# auxiliary data frames, for plotting points
m_med = m[ m$Statistic%in%c("Median","Lower (HPD95)","Upper (HPD95)")&m$Parameter=='w_biome', ]
m_gmed =  m[ m$Statistic=="Median"&m$Parameter=="w_biome", ] %>% group_by(Strength, Condition) %>% summarise(Value=mean(Value))
m_gmed$Statistic = "Mean-median"
m_bf = m[ m$Parameter=='w_biome'&m$Statistic=='Median',]
m_bf2 = m_bf %>% group_by(Strength, BayesFactorCategory) %>% summarise( Count=n() )

# compute info about HPD coverage
m_hpd = m[ m$Statistic%in%c('Lower (HPD95)','Upper (HPD95)')&m$Parameter=='w_biome',]
m_hpd$Statistic = factor(m_hpd$Statistic, ordered=T, levels=rev(levels(m_hpd$Statistic)))
f_covered = rep(0,5)
names(f_covered) = levels(m_hpd$Condition)
for (i in 1:5) {
    s = names(f_covered)[i]
    m_tmp = m_hpd[ m_hpd$Statistic=="Lower (HPD95)" & m_hpd$Condition==s, ]
    f_covered[i] = sum(m_tmp$Covered=="Covered") / nrow(m_tmp)
}
f_cov_1 = sum( (2/sort(m_hpd[ m_hpd$Condition=="Null"&m_hpd$Statistic=="Upper (HPD95)", ]$BayesFactor)) > 0.025 )/100
f_covered[1] = f_cov_1
f_covered = format(f_covered,nsmall=2)
m_covered = data.frame( Condition=names(f_covered), Frequency=f_covered )

# Plot A: posterior estimates, medians, coverage
dodge <- position_dodge(width=0.7)  
p = ggplot(m_med, aes(x=Condition, y=Value, color=Strength))
p = p + geom_point(mapping=aes(shape=Statistic), alpha=0.5, position=dodge)
p = p + geom_point(data=m_true[ m_true$Parameter=='w_biome',], mapping=aes(x=Condition, y=Value, color=Strength, shape=Statistic), size=4)
p = p + geom_point(data=m_gmed, mapping=aes(x=Condition, y=Value, color=Strength, shape=Statistic), size=4)
p = p + scale_color_manual( values=c("forestgreen","firebrick3","dodgerblue","goldenrod1","gray"))
p = p + scale_shape_manual( values=c("Upper (HPD95)"=2,"Median"=1,"Lower (HPD95)"=6,"Mean-median"=16, "True"=15),
                            breaks=c("Median","Upper (HPD95)","Lower (HPD95)","Mean-median","True"))
p = p + geom_text( data=m_covered, mapping=aes(x=Condition, y=1.05, label=Frequency), color="black" )
p = p + annotate( geom="text", x=3, y=1.125, label="Coverage frequency")
p = p + theme( axis.text.x = element_text(angle=90),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"))
p = p + scale_y_continuous( breaks=seq(0,1,0.2), limits=c(0,1.15))               
lbl2 = expression(paste("Posterior estimates of ", w[B], sep=""))
p = p + ylab(lbl2)
p1 = p


# Plot B: Bayes Factors
p = ggplot(m_bf, aes(x=Condition, fill=BayesFactorCategory), stat='count')
p = p + scale_fill_manual( name="Bayes Factor\nCategory",
                           values=rev( gray.colors(n=6, start=0, end=0.9)),
                           labels=c("None (<1)",
                                    "Insubstantial (1-3)",
                                    "Substantial (3-10)",
                                    "Strong (10-30)",
                                    "V. Strong (30-100)",
                                    "Decisive (>100"))
p = p + geom_bar()
lbl2 = expression(paste("Model support for ", w[B], " > 0", sep=""))
p = p + ylab(lbl2)
p = p + theme( axis.text.x = element_text(angle=90),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"))
p = p + scale_y_continuous( breaks=seq(0,100,25), limits=c(0,115))  
p2 = p

# Combine plots A and B
pg =  plot_grid( p1, p2, labels=c("(A)","(B)"), align="v", rel_widths = c(2,1.5), nrow=1)

# Save to file
CairoPDF( plot_fn, height=4.5, width=10)
print(pg)
dev.off()
