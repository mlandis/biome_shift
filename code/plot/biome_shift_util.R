# Fig 2

.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}


# Fig 3

fig3_add_epoch_times <- function( p, max_age, x_offset=0, dy=0.1 ) {

    dy2 = 0.2
    max_x = 0
    max_y = -0.005
    epoch_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")

    x_pos = max_x-c(max_age, 65, 56, 48, 34, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2

    box_bg = geom_rect( xmin=0, xmax=70, ymin=max_y-dy, ymax=max_y, fill="white", alpha=1.0, size=0)
    p = append_layers(p, box_bg, position = "top")

    for (k in 2:(length(x_pos))) {
        box_col = "gray92"
        if (k %% 2 == 0) box_col = "white"
        box = geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=max_y-dy2, ymax=y_pos[k], fill=box_col, size=0 )
        p = append_layers(p, box, position = "top")
    }
    for (k in 1:length(epoch_names)) {
        p = p + annotate("text", x=-x_pos_mid[k], y=max_y-dy, label=epoch_names[k], hjust=0.5, size=2.5)
    }
    return(p)

}


# Fig 4
make_bf = function(x) {
    k = kdensity(x=x, kernel='beta', bw=0.01)
    posterior = k(0)
    prior = 2 # equal for all points in simplex under flat dirichlet
    bf = prior / posterior
    return(bf)
}

is_covered = function(hpd, tp, mdl, str) {
    x = tp[ tp$model==mdl & tp$strength==str, ]
    print(x)
    cn = colnames(hpd)
    print(cn)
    hpd_low = hpd[c("lower"),cn]
    hpd_up = hpd[c("upper"),cn]
    covered = (x[cn] > hpd_low & x[cn] < hpd_up)
    
    return( sapply( covered, function(y) { if (y) { "Covered" } else { "Not covered" } } ))
}

lnL_beta = function(par, ...) {
    z = list(...)[[1]]
    negLnL = -sum( dbeta( z, par[1], par[2], log=T) )
    return(negLnL)
}

lnL_dirichlet = function(par, ...) {
    z = list(...)[[1]]
    negLnL = -sum(log(ddirichlet( z, alpha=par)))
    return(negLnL)
}

fit_beta_mle = function(x) {
    o = optim( par=c(1,1), fn=lnL_beta, ...=x, method="L-BFGS-B", lower=c(0,0), upper=c(Inf,Inf) )
    opar = c(shape1=o$par[1], shape2=o$par[2])
    return(opar)
}

fit_dirichlet_mle = function(x) {
    o = optim( par=c(1,1,1), fn=lnL_dirichlet, ...=x, method="L-BFGS-B", lower=c(0,0,0), upper=c(Inf,Inf,Inf) )
    opar = c(shape1=o$par[1], shape2=o$par[2])
    return(opar)
}

compute_bf = function(x, bw=0.05) {
    k = kdensity(x=x, kernel='beta')
    mlike_prior = 2 # dbeta(x=0,shape1=1, shape2=2)
    mlike_posterior = k(0) 
    bf = mlike_prior / mlike_posterior
    return(bf)
}

make_bf_cat = function(x) {
    cat = c(-Inf, 1, 3, 10, 30, 100)
    cat_lbl = c("None", "Insubstantial", "Substantial", "Strong", "V. Strong", "Decisive")
    idx = sapply(x, function(y) { idx = which(y>cat); return(length(idx)) })
    ret = factor(cat_lbl[idx], ordered=T, levels=rev(cat_lbl))
    return(ret)
}

# Fig 5

fig5_add_epoch_times <- function( p, max_age, dy_bars, dy_text ) {
    
    max_x = max(p$data$x)
    max_y = max(p$data$y)
    epoch_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
 
    x_pos = max_x-c(75, 65, 56, 48, 33.9, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2 

    for (k in 2:(length(x_pos))) {
        box_col = "gray92"
        if (k %% 2 == 0) box_col = "white"
        box = geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=dy_bars, ymax=y_pos[k], fill=box_col )
        p = append_layers(p, box, position = "bottom")
    }
    for (k in 1:length(epoch_names)) {
        p = p + annotate( geom="text", label=epoch_names[k], x=x_pos_mid[k], y=dy_text, hjust=0.5, size=2.5)
    }
    return(p)

}

# FIG 6


fig6_add_epoch_times <- function( p, max_age, x_offset=0, dy=0.1 ) {
    
    dy2 = 0.2
    max_x = 0
    max_y = 0
    epoch_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")

    x_pos = max_x-c(max_age, 65, 56, 48, 34, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2 

    for (k in 2:(length(x_pos))) {
        box_col = "gray92"
        if (k %% 2 == 0) box_col = "white"
        box = geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=0-dy2, ymax=y_pos[k], fill=box_col )
        p = append_layers(p, box, position = "bottom")
    }
    for (k in 1:length(epoch_names)) {
        #print(x_pos_mid[k])
        #print(epoch_names[k])
        p = p + annotate("text", x=-x_pos_mid[k], y=0-dy, label=epoch_names[k], hjust=0.5, size=2.5)
    }
    return(p)

}

make_bioregion = function(s) {
    if      (s==01) c(1,1)   
    else if (s==02) c(1,2)
    else if (s==03) c(1,3)
    else if (s==04) c(1,4)
    else if (s==05) c(1,5)
    else if (s==06) c(1,6)
    else if (s==07) c(2,1)
    else if (s==08) c(2,2)
    else if (s==09) c(2,3)
    else if (s==10) c(2,4)
    else if (s==11) c(2,5)
    else if (s==12) c(2,6)
    else if (s==13) c(3,1)
    else if (s==14) c(3,2)
    else if (s==15) c(3,3)
    else if (s==16) c(3,4)
    else if (s==17) c(3,5)
    else if (s==18) c(3,6)

}

get_match_idx = function( s_idx, time_idx, epoch_times, atlas_list ) {

    epoch_idx = get_epoch_idx( time_idx, epoch_times )
    br_idx = make_bioregion(s_idx)
    b_idx = br_idx[1]
    r_idx = br_idx[2]
    ret = rep(FALSE, length(epoch_idx))
    for (i in 1:length(epoch_idx)) {
        t = epoch_idx[i]
        if (atlas_list[[b_idx]][[t]][r_idx,r_idx] > 0) {
            ret[i] = TRUE
        }
    }
    return(ret)
}

get_epoch_idx = function( time_idx, epoch_times ) {
    ret = rep(0L, length(time_idx))
    for (i in 1:length(time_idx)) {
        ret[i] = rev(which( time_idx[i] < epoch_times$start_age ))[1]
    }
    return(ret)
}

make_graphs = function(atlas_fp, atlas_fn) {
    
    graph = c("null", "land", "tropical", "warm", "cold")
    n_graphs = length(graph)
    
    atlas = read.table(atlas_fn,header=T,sep=",")
    n_epochs = nrow(atlas)
    
    area_names = c("SE As","E As", "Eur", "N Am", "C Am", "S Am")
    dat = list()
    
    for (i in 1:n_epochs) {
        dat[[ i ]] = list()
    }
    for (i in 1:n_graphs) {
        files = list.files(path=paste(atlas_fp, graph[i], "/", sep=""), full.names = T)
        n_files = length(files)
        for (j in 1:n_files) {
            tmp = read.table(files[j], sep=",")
            rownames(tmp) = area_names
            colnames(tmp) = area_names
            dat[[j]][[i]] = as.matrix(tmp) + 1
        }
    }
    
    return(dat)
}


make_state = function(b1) {
    
    # get anc states
    b = rbind( as.data.frame(b1, stringsAsFactors=F))
    
    # sort joint events by transition time
    b = b[order(b$transition_time, decreasing=T),]
    
    # get the first valid start state for each character
    if (b$transition_type[1]=="no_change") {
        
        x1 = b$branch_start_time
        x2 = b$branch_end_time
        
    } else if (b$transition_type[1] == "anagenetic") {
        
        x1 = c( b$branch_start_time[1], b$transition_time )
        x2 = c( b$transition_time, b$branch_end_time[1] )
        
        # add event for last event segment
        last_row = nrow(b)
        b = rbind(as.data.frame(b, stringsAsFactors=F), 
                  as.data.frame(b[last_row,], stringsAsFactors=F))
        b[last_row+1, c("start_state","transition_time")] = b[last_row, c("end_state","branch_end_time")]
    }
    
    b = cbind(b, x1=x1, x2=x2, s1=b$start_state)
    
    return(b)
    
}


make_match_bins = function( d, epoch_times, atlas_list ) {
    
    ret = data.frame(NULL, stringsAsFactors=FALSE)

    ages = unique(d$age)
    for (i in 1:length(ages)) {
        di = d[ d$age==ages[i], ]
        match_idx = c()
        for (j in 1:nrow(di)) {
            dij = di[j,]
            match_idx[ dij$StateValue ] = get_match_idx( dij$StateValue, dij$age, epoch_times, atlas_list )
        }
        di$match = match_idx
        
        d_match = di[ di$match, ]
        d_mismatch = di[ !di$match, ]
        
        v_match = c( age=ages[i], count=sum(d_match$count), Support=1, prop=sum(d_match$prop), match=TRUE)
        v_mismatch = c( age=ages[i], count=sum(d_mismatch$count), Support=1, prop=sum(d_mismatch$prop), match=FALSE)
        
        ret = rbind(ret, v_match, v_mismatch)
        
    }
    
    colnames(ret) = c("age","count","Support","prop","match")
    ret$match = c("mismatch", "match")[ret$match+1]
    ret$match = factor( ret$match, ordered=T, levels=c("mismatch","match") )
    ret$age2 = as.numeric(ret$age + 0.5)
    
    return(ret)
   
}


make_lstt_dat = function(fn, n_states, f_burn=0, thinby=1) { 
    
    # read in data file
    stoch = read.csv(fn, sep="\t", stringsAsFactors=F)
    stoch = stoch[ stoch$transition_type != "cladogenetic", ]
    stoch$transition_time[ stoch$transition_type=="no_change" ] = stoch$branch_start_time[ stoch$transition_type=="no_change" ]
    
    # filter MCMC samples
    iterations = unique(stoch$iteration)
    n_burn = max(1, f_burn*length(iterations))
    iterations = iterations[n_burn:length(iterations)]
    iterations = iterations[ seq(1, length(iterations), by=thinby) ]
    
    # get branch indices
    branches = 1:max(unique(stoch$parent_index), na.rm=T)
    
    # STAGE 1, process posterior stochastic mappings
    stoch_bg_biome = data.frame(stringsAsFactors = F)
    branch_bg_biome = list()
    bidx = 1
    for (i in 1:length(iterations)) {
        
        # get biome and biogeography stochastic mappings per iteration
        it = iterations[i]
        cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
        sample = stoch[ stoch$iteration==it, ]
        
        # loop over branches
        for (j in 1:length(branches)) {
            
            # get biome and biogeography stochastic mappings per branch
            nd_idx = branches[j]
            branch = sample[ sample$node_index==nd_idx, ]
            
            # interleave biome and biogeography stochastic mappings
            branch_bg_biome[[bidx]] = as.data.frame( make_state(branch), stringsAsFactors=F )
            bidx = bidx + 1
        }
    }
    stoch_bg_biome = rbindlist( branch_bg_biome )

    # STAGE 2, put data in time bins
    max_time = 90
    bin_width = 1
    n_bins = max_time / bin_width
    
    # state bins will contain the number of samples in each state for each time
    state_bins = array(0, dim=c(n_states, n_bins))
    
    # create time bins, 0 to 1, 1 to 2, 2 to 3, etc.
    ages = seq(0.0, max_time, by=bin_width)

    # put counts into matrix
    for (i in 1:length(iterations)) {
        
        # get biome and biogeography stochastic mappings per iteration
        it = iterations[i]
        cat("Stage 2, processing iteration ",it," / ", max(iterations), "\n", sep="")
        sample = stoch_bg_biome[ stoch_bg_biome$iteration==it, ]
        
        # loop over branches
        for (j in 1:length(branches)) {
            branch_sample = sample[ sample$node_index==j, ]
            
            for (k in 1:nrow(branch_sample)) {
                s_idx = branch_sample$start_state[k] + 1
                start_age = branch_sample$x1[k]
                int_start_age = floor(start_age)
                end_age = branch_sample$x2[k]
                int_end_age = floor(end_age)
                age_bins = int_start_age:int_end_age * bin_width
                time_idx = age_bins + 1
                n_time = length(time_idx)
                
                bin_weights = rep(1, n_time)
                if (n_time==0) {
                    error("zero events??")
                } else if (n_time==1) {
                    #print( c(n_time, start_age, end_age) ) 
                    bin_weights[1] = (start_age-end_age) / bin_width
                } else if (n_time>1) {
                    bin_weights[1] = (start_age-int_start_age)
                    bin_weights[n_time] = bin_width - (end_age-int_end_age)
                }

                state_bins[ s_idx, time_idx ] = state_bins[ s_idx, time_idx ] + bin_weights
            }
        }
    }
    
    # get state bins
    dat_plot_2 = matrix(nrow=0, ncol=5)
    for (i in 1:dim(state_bins)[1]) {
        for (j in 1:dim(state_bins)[2]) {
            c_ij = state_bins[i,j]
            dat_plot_2 = rbind(dat_plot_2, c( ages[j], c_ij, i, st_lbl[i], 1))
        }
    }
    
    # convert into something ggplotable
    d2 = data.frame(dat_plot_2,  stringsAsFactors=FALSE)
    colnames(d2) = c("age","count","StateValue", "State","Support")
    d2$age = as.numeric(d2$age)
    d2$age2 = as.numeric(d2$age + 0.5)
    d2$count = as.numeric(d2$count)
    d2$prop = as.numeric(d2$count / length(iterations) )
    d2$StateValue = as.numeric(d2$StateValue)
    d2$State = factor(d2$State, ordered=T, levels=st_lbl)
    
    # return
    return(d2)
}