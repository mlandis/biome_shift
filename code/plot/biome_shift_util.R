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



# fig 8

make_graphs = function(atlas_fp, atlas_fn) {
    
    graph = c("null", "land", "tropical", "warm", "cold")
    n_graphs = length(graph)
    
    atlas = read.table(atlas_fn, header=T, sep=",")
    n_epochs = nrow(atlas)
    
    area_names = c("SE As","E As", "Eur", "N Am", "C Am", "S Am")
    dat = list()
    
    for (i in 1:n_epochs) {
        dat[[ i ]] = list()
    }
    for (i in 1:n_graphs) {
        files = list.files(path=paste0(atlas_fp, graph[i], "/"), full.names = T)
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


make_pi = function(params) {
    
    col_names = paste("rf_simplex.",1:18,".",sep="")
    pi = as.vector(params[,col_names])
    return(pi)
}



# connectivity
make_Q_epoch = function(params, graphs) {
    
    Q_epoch = list()
    for (i in 1:length(graphs)) {
        Q_epoch[[i]] = make_Q(params, graphs[[i]])
    }
    return(Q_epoch)
}

# probabilities of next event type
make_P = function(Q) {
    P = Q
    for (i in 1:nrow(Q)) {
        P[i,] = Q[i,] / -Q[i,i]
        P[i,i] = 0
    }
    return(P)
}

# connectivity
make_Q = function(params, graph) {
    
    na=6
    nb=3
    ns = na * nb
    nt = length(graphs)
    
    # clock
    clock = params$rate_biome
    
    # raw dispersal rates
    dr = matrix(1, nrow=na, ncol=na)
    dr = dr * params$dr_base * clock
    
    # raw biome shift rates
    br = matrix(0, nrow=nb, ncol=nb)
    br[1,2] = params$biome_rates.1. 
    br[2,3] = params$biome_rates.2.
    br[2,1] = params$biome_rates.3.
    br[3,2] = params$biome_rates.4.
    br = br * params$br_base * clock

    # graph weight params
    w_null = params$w_null
    w_land = params$w_land
    w_biome = params$w_biome
    w_not_biome = 1 - w_biome
    
    # connectivity params
    y = c( params$w_adj.1., params$w_adj.2., params$w_adj.3.)
    names(y) = c("y_marg","y_weak","y_strong")
    
    #cat( w_null, w_land, w_biome, y, "\n")
   
    Q = matrix(NA, nrow=ns, ncol=ns)
    #rownames(Q) <- colnames(Q) <- abnames
    
    for (bi in 1:nb) {
        for (ai in 1:na) {
            for (bj in 1:nb) {
                for (aj in 1:na) {
                    si = (bi-1) * na + ai
                    sj = (bj-1) * na + aj
        
                    # cat( "(",bi,",",ai,") -> (",bj,",",aj,"); ", si, "->", sj, "\n")
                    
                    biome_diff = F
                    region_diff = F
                    
                    if (si==sj) {
                        # do nothing
                    } else if (ai != aj && bi != bj) {
                        # do nothing
                    } else if (ai != aj) {
                        region_diff = T
                    } else if (bi != bj) {
                        biome_diff = T
                    }
                    
                    g_land = graph[[2]]
                    g_biome = graph[[2+bj]]
                    
                    Q[si,sj] = 0
                    if ( biome_diff ) {
                        # biome shift
                        y_biome_idx = g_biome[aj,aj]
                        Q[si,sj] = br[bi,bj] * (w_not_biome + w_biome * y[ y_biome_idx ] )
                    } else if ( region_diff ) {
                        # dispersal
                        y_land_idx = g_land[ai,aj]
                        y_biome_idx = g_biome[ai,aj]
                        Q[si,sj] = dr[ai,aj] * (w_null + w_land*y[y_land_idx] + w_biome*y[y_biome_idx])
                    }
                }
            }
        }
    }
    
    for (i in 1:ns) {
        Q[i,i] = -sum(Q[i,])
    }
    
    return(Q)
}

make_prob_ij = function(Q) {
    ns = ncol(Q)
    pij = matrix(0, ns, ns)
    
    # Prob( j | i ) -- forwards in time
    for (i in 1:ns) {
        for (j in 1:ns) {
            pij[i,j] = (Q[i,j] / -Q[i,i])
        }
    }
    diag(pij) = 0
    
    return(pij)
}

make_prob_ji = function(Q, pi) {
    
    ns = ncol(Q)
    
    pij = make_prob_ij(Q)
    pji = matrix(0, nrow=ns, ncol=ns)
    #pi = expm(Q * 100)
    
    # Prob( i | j ) -- backwards in time
    for (i in 1:ns) {
        for (j in 1:ns) {
            print(pi[i])
            print(pi[j])
            print(pij[i,j])
            pji[i,j] = pij[i,j] * as.numeric(pi[i]) / as.numeric(pi[j])
        }
    }
    return(pji)
}


darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

lighten <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col*factor
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
    col
}


pause = function() {
    if (interactive()) {
        invisible(readline(prompt = "Press <Enter> to continue..."))
    } else {
        cat("Press <Enter> to continue...")
        invisible(readLines(file("stdin"), 1))
    }
}

get_first_lbl = function() {
    lbl = c("Biome First","Biome Flight","Biome Reversal","Region First","Region Flight","Region Reversal")
    return(lbl)
}


get_state_lbl = function() {
    lbl = c()
    for (i1 in c("Trop","Warm","Cold")) {
        for (j1 in c("SEAs","EAs","Eur","NAm","CAm","SAm")) {
            for (i2 in c("Trop","Warm","Cold")) {
                for (j2 in c("SEAs","EAs","Eur","NAm","CAm","SAm")) {
                    if (i1 != i2 && j1 != j2) {
                        next
                    } else if (i1 == i2 && j1 == j2) {
                        next
                    } else {
                        lbl = c( lbl, paste0( i1, "+", j1, " to ", i2, "+", j2))
                    }
                }
            }
        }
    }
    return(lbl)
}

make_biome_name = function(s) {
    if (s==1) {
        "Trop"
    } else if (s==2) {
        "Warm"
    } else if (s==3) {
        "Cold"
    } else {
        error("Unknown biome state")
    }
}

make_region_name = function(s) {
    if (s==1) {
        "SEAs"
    } else if (s==2) {
        "EAs"
    } else if (s==3) {
        "Eur"
    } else if (s==4) {
        "NAm"
    } else if (s==5) {
        "CAm"
    } else if (s==6) {
        "SAm"
    } else {
        error("Unknown biome state")
    }
}

find_epoch = function(x, epoch_times=c(65, 56, 48, 34, 23, 16, 5.3, 0.0)) {
    z = sapply( x, function(y) { sum(y <= epoch_times)+1 })
    z[ z>length(epoch_times) ] = length(epoch_times)
    return(z)
}


make_triplets = function(x, params, graph, subroot=F) {
    
    # get root node and state
    idx_root = max(x$node_index)
    x_root = x[ x$parent_index==idx_root, ]$start_state[1]
    
    # get root state Q matrix
    Q_root = make_Q( params, graph )
    
    # turn into probability matrix for next event being type ij
    P_root = make_P( Q_root )
    
    # get root state probabilities
    rf_lbl = paste("rf_simplex.",1:ncol(Q_root),".",sep="")
    rf_prob = params[,c(rf_lbl)]
    
    # compute the probability of the previous state being X_subroot=i given you end in state X_root=j
    
    # Define:
    #   p1 : P( X_subroot=i | X_root=j, Q(m_root) ) -- want, but don't have
    #   p2 : P( X_root=j | X_subroot=i, Q(m_root) ) -- prob event of j | i
    #   p3 : pi( X_subroot=i | Q(m_root) ) -- stationary prob i
    #   p4 : pi( X_root=j | Q(m_root) -- stationary prob j
    #   p1 = p2 * p3 / p4 -- Bayes
    
    # fill in our probabilities for starting in i given we end in j
    rev_prob = rep(0, 18)
    for (i in 1:length(rf_prob)) {
        rev_prob[i] = P_root[i, x_root] * rf_prob[i] / rf_prob[x_root]
    }
    rev_prob = unlist(rev_prob)
    rev_prob = rev_prob / sum(rev_prob)
    
    # sample subroot state
    x_subroot = -1
    if (subroot) {
        x_subroot = sample(1:18, size=1, prob=rev_prob, replace=TRUE)
    }
    state_root = c(idx_root, -1, x_subroot, x_root)
    count_root = 0
    
    # start recursion
    events = make_triplet_recursion(x, idx_root, state_root, count_root)
    
    # return count mtx
    return(events)
}


make_triplet_df = function(d,params,graph,subroot=F) {
    # loop over iterations
    iterations = sort(unique(d$iteration))
    n_it = length(iterations)
    
    evt = NULL
    df_list = list()
    df_list_idx = 1
    for ( i in 1:n_it ) {
        # get iteration
        it = iterations[i]
        # get data for iteration
        dtmp = d[d$iteration==it, ]
        # get parameters for iteration
        partmp = params[ params$Iteration==it, ]
        # make event triplets for iteration
        evttmp = data.frame(it, make_triplets( dtmp, partmp, graph, subroot ))
        # add event triplets to data frame
        #evt = rbind(evt, evttmp)
        df_list[[df_list_idx]] = evttmp
        df_list_idx = df_list_idx + 1
    }
    evt = rbindlist(df_list)
    
    # return data frame
    df = data.frame(evt); colnames(df) = c("iteration","node_index","S1","S2","S3","age"); rownames(df)=NULL
    return(df)
}


make_triplet_recursion = function(x, idx, anc_state, count) {
    events = matrix(NA, nrow=0, ncol=5)
    
    # get node
    x_nd = x[ x$node_index==idx, ]
    
    # collect events if something happened on branch
    nr = nrow(x_nd)
    
    if (nr > 0 && x_nd$transition_type != "no_change" && !is.na(x_nd$parent_index[1]) ) {

        # loop over events
        for (i in 1:nr) {
            
            # get prev and curr state
            state_0 = anc_state[3]
            state_1 = x_nd$start_state[i]
            state_2 = x_nd$end_state[i]
            age = x_nd$transition_time[i]
            
            # update triplet vector
            anc_state = c(idx, state_0, state_1, state_2, age)
            
            # append to event dataframe
            events = rbind(events, anc_state)
        }
    }
    
    # recurse to children
    ch1_events = NULL
    ch1_idx = x_nd$child1_index[1]
    #print(ch1_idx)
    if (!is.na(ch1_idx)) {
        ch1_events = make_triplet_recursion(x, ch1_idx, anc_state, count)
    }
    ch2_events = NULL
    ch2_idx = x_nd$child2_index[1]
    #print(ch2_idx)
    if (!is.na(ch2_idx)) {
        ch2_events = make_triplet_recursion(x, ch2_idx, anc_state, count)
    }
    
    # collect events
    events = rbind(events, ch1_events, ch2_events)
    
    return(events)
}


make_triplet_table = function(df) {
    
    iterations = sort(unique(df$iteration))
    n_it = length(iterations)
    df = df[ df$S1>0, ]
    lbl = get_first_lbl()
    
    df_tzz = data.frame(NULL, stringsAsFactors = F)
    for (i in 1:n_it) {
        it = iterations[i]
        df_it = df[df$iteration==it,]
        zz = apply(df_it[c("S1","S2","S3")], 1, function(z) { make_tx3_type( z[1], z[2], z[3]) } )
        fzz = factor(zz, levels=lbl)
        tzz = table(fzz)
        df_tzz = rbind( df_tzz, tzz, df[df$iteration==it, c("age")] )
        colnames(df_tzz)=names(tzz)
    }

    return(df_tzz)
}

make_tx2_type = function(s1,s2) {
    x1 = make_bioregion(s1)
    x2 = make_bioregion(s2)
    val = paste0( make_biome_name(x1[1]), "+", make_region_name(x1[2]), " to ", make_biome_name(x2[1]), "+", make_region_name(x2[2]))
    return(val)
}

make_tx3_type = function(s1, s2, s3) {
    x1 = make_bioregion(s1)
    x2 = make_bioregion(s2)
    x3 = make_bioregion(s3)
    
    biome_12_match  = ( x1[1] == x2[1] )
    biome_23_match  = ( x2[1] == x3[1] )
    biome_13_match  = ( x1[1] == x3[1] )
    
    region_12_match = ( x1[2] == x2[2] )
    region_23_match = ( x2[2] == x3[2] )
    region_13_match = ( x1[2] == x3[2] )
    
    s = ""
    # RegionFirst    A,X -> A,Y -> B,Y    
    if (biome_12_match && !region_12_match && !biome_23_match && region_23_match) {
        s = "Region First"
    }
    # BiomeFirst     A,X -> B,X -> B,Y
    else if (!biome_12_match && region_12_match && biome_23_match && !region_23_match) {
        s = "Biome First"
    }
    # BiomeRev       A,X -> B,X -> A,X
    else if (!biome_12_match && region_12_match && biome_13_match && region_23_match) {
        s = "Biome Reversal"
    }
    # BiomeRun       A,X -> B,X -> C,X
    else if (!biome_12_match && !biome_23_match && !biome_13_match && region_12_match && region_23_match) {
        s = "Biome Flight"
    }
    # RegionRev      A,X -> A,Y -> A,X
    else if (biome_12_match && biome_23_match && !region_12_match && region_13_match) {
        s = "Region Reversal"
    }
    # RegionRun             A,X -> A,Y -> A,Z
    else if (biome_12_match && biome_23_match && !region_12_match && !region_23_match && !region_13_match) {
        s = "Region Flight"
    }
    else {
        s = "Unknown!"
    }

    return(s)
}

make_state = function(b1) {
    
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


make_trip_plot = function(dat, the_guide=NULL) {
    
    line_colors2 = c( "darkorchid", "gold", "dodgerblue", "firebrick3", "olivedrab", "violetred1")
    lbl_types_long = c("Biome First\n AX,BX,BY",
                     "Biome Flight\n  AX,BX,CX",
                     "Biome Reversal\n      AX,BX,AX",
                     "Region First\n  AX,AY,BX",
                     "Region Flight\n   AX,AY,BZ",
                     "Region Reversal\n      AX,AY,AX")
    lbl_types_short = c("Biome First",
                     "Biome Flight",
                     "Biome Reversal",
                     "Region First",
                     "Region Flight",
                     "Region Reversal")
    
    names(line_colors2)=lbl_types_short
    
    ppb = ggplot(dat, aes(color=Var2))
    
    
    dat_rect = data.frame( center=c(1,3,5,7) )
    dat_rect$left = dat_rect$center - 0.5
    dat_rect$right = dat_rect$center + 0.5
    ppb = ppb + geom_rect(data=dat_rect, mapping=aes(xmin=left, xmax=right, ymin=0, ymax=1), fill="gray92", color=NA)
    
    ppb = ppb + geom_segment(data=dat, mapping=aes(x=x, xend=x, y=lower80, yend=upper80), size=1.25, alpha=0.5)
    ppb = ppb + geom_segment(data=dat, mapping=aes(x=x, xend=x, y=lower95, yend=upper95), size=0.65, alpha=0.5)
    ppb = ppb + geom_line(data=dat,  mapping=aes(x=x, y=value, color=Var2), linetype=2)
    ppb = ppb + geom_point(data=dat, mapping=aes(x=x, y=value, color=Var2),size=2)
    ppb = ppb + geom_point(data=dat, mapping=aes(x=x, y=value),size=0.5, color="white")
    
    ppb = ppb + scale_color_manual(values=line_colors2, labels=lbl_types_long)
    ppb = ppb + xlab("Time interval") + ylab("Proportion of events")
    ppb = ppb + scale_x_continuous(breaks=1:7, labels=c("Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent"))
    ppb = ppb + ylim(0,1)
    ppb = ppb + theme_classic()
    ppb = ppb + theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="top",
                      legend.justification="center",
                      legend.title = element_text(size=12),
                      legend.text = element_text(size=8),
                      legend.key.height = unit(0.1, "cm"),
                      legend.key.width = unit(0.5,"cm"),
                      legend.spacing.x = unit(0.3, 'cm'),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (!is.null(the_guide)) {
        ppb = ppb + the_guide
    } else {
        ppb = ppb + guides(color=FALSE)
    }
    return(ppb)
}

make_event_plot = function(dat, the_guide=NULL) {
    
    lbl_types = c("Warm+EAs to Cold+EAs", "Cold+EAs to Warm+EAs",
              "Warm+NAm to Cold+NAm", "Cold+NAm to Warm+NAm",
              "Warm+EAs to Warm+NAm", "Warm+NAm to Warm+EAs",
              "Cold+EAs to Cold+NAm", "Cold+NAm to Cold+EAs")
    
    line_colors = c(lighten("dodgerblue",1.4), darken("dodgerblue",1.2),
                    lighten("gold", 1.0), darken("gold",1.25),
                    lighten("firebrick", 1.4), darken("firebrick",1.1),
                    lighten("mediumorchid1", 1.4), darken("mediumorchid",1.1))
    names(line_colors)=lbl_types
    
    dat_rect = data.frame( center=c(1,3,5,7) )
    dat_rect$left = dat_rect$center - 0.5
    dat_rect$right = dat_rect$center + 0.5
    
    
    ppb = ggplot(dat, aes(color=Var2))
    ppb = ppb + geom_rect(data=dat_rect, mapping=aes(xmin=left, xmax=right, ymin=0, ymax=1), fill="gray92", color=NA)
    ppb = ppb + geom_segment(data=dat, mapping=aes(x=x, xend=x, y=lower80, yend=upper80), size=1.25, alpha=0.5)
    ppb = ppb + geom_segment(data=dat, mapping=aes(x=x, xend=x, y=lower95, yend=upper95), size=0.65, alpha=0.5)
    ppb = ppb + geom_line(data=dat,  mapping=aes(x=x, y=value, color=Var2), linetype=2)
    ppb = ppb + geom_point(data=dat, mapping=aes(x=x, y=value, color=Var2),size=2)
    ppb = ppb + geom_point(data=dat, mapping=aes(x=x, y=value),size=0.5, color="white")
    
    ppb = ppb + scale_color_manual(values=line_colors)
    ppb = ppb + xlab("Time interval") + ylab("Proportion of events")
    ppb = ppb + scale_x_continuous(breaks=1:7, labels=c("Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent"))
    ppb = ppb + ylim(0,1)
    ppb = ppb + theme_classic()
    ppb = ppb + theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="top",
                      legend.justification="center",
                      legend.title = element_text(size=12),
                      legend.text = element_text(size=8),
                      legend.key.height = unit(0.1, "cm"),
                      legend.key.width = unit(0.5,"cm"),
                      legend.spacing.x = unit(0.3, 'cm'),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
    if (!is.null(the_guide)) {
        ppb = ppb + the_guide
    } else {
        ppb = ppb + guides(color=FALSE)
    }
    #ppb = ppb + theme(axis.text.x = element_text(angle = 90, vjust=0.5))
    return(ppb)
}

make_dat_triplet_raw = function(fn, n_states, f_burn=0.0, thinby=1) { 
    
    # files
    #out_fp = paste(fp, "output_trim/", sep="")
    #fn = paste(out_fp, base_fn, ".history.tsv", sep="")
    
    # read in data file
    stoch = read.csv(fn, sep="\t", stringsAsFactors=F)
    stoch = stoch[ stoch$transition_type != "cladogenetic", ]
    stoch$transition_time[ stoch$transition_type=="no_change" ] = stoch$branch_start_time[ stoch$transition_type=="no_change" ]
    
    # filter MCMC samples
    iterations = unique(stoch$iteration)
    n_burn = max(1, f_burn*length(iterations))
    iterations = iterations[n_burn:length(iterations)]
    iterations = iterations[ seq(1, length(iterations), by=thinby)  ]
        
    # get branch indexing
    branches = 1:max(unique(stoch$parent_index), na.rm=T)
    
    # STAGE 1, gather data
    # loop over iterations
    df_list = list()
    df_list_idx = 1
   
    for (i in 1:length(iterations)) {
        
        # get biome and biogeography stochastic mappings per iteration
        it = iterations[i]
        cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
        sample = stoch[ stoch$iteration==it, ]
        
        # loop over branches
        br_list = list()
        br_list_idx = 1
        
        dtmp = data.frame(stringsAsFactors = F)
        for (j in 1:length(branches)) {
            
            # get biome and biogeography stochastic mappings per branch
            nd_idx = branches[j]
            branch = sample[ sample$node_index==nd_idx, ]
            
            # interleave biome and biogeography stochastic mappings
            branch_bg_biome = branch
            branch_bg_biome$start_state = branch_bg_biome$start_state + 1
            branch_bg_biome$end_state = branch_bg_biome$end_state + 1
            
            dtmp = as.data.frame(branch_bg_biome, stringsAsFactors=F)
            br_list[[br_list_idx]] = dtmp
            br_list_idx = br_list_idx + 1
            
        }
        df_br = rbindlist(br_list)
        df_list[[df_list_idx]] = df_br
        df_list_idx = df_list_idx + 1
        
    }
    stoch_bg_biome = rbindlist(df_list)
    
    return(stoch_bg_biome)
}


make_state_lstt_epoch = function(df,  epoch_times, match_lbl) {
    
    iterations = sort(unique(df$iteration))
    n_it = length(iterations)
    df = df[ df$S1>0, ]
    lbl = get_state_lbl()
    n_events = length(lbl)
    n_epochs = length(epoch_times) # - 1
    df$age = ceiling(df$age)
    df$epoch = find_epoch(df$age)
    
    # matrix of (times x event_types x iterations)
    # we then use this to compute HPDs for events
    x = array(0, dim=c(n_epochs, n_events, n_it))
    dimnames(x) = list(NULL,lbl,NULL)
    
    for (i in 1:n_it) {
        it = iterations[i]
        df_it = df[df$iteration==it,]
        
        zz = apply(df_it[,c("S2","S3")], 1, function(z) { make_tx2_type( z[1], z[2]) } )
        fzz = factor(zz,levels=lbl)
        
        for (j in 1:length(zz)) {
            k = df_it$epoch[j]
            x[k, fzz[j], i] = x[k, fzz[j], i] + 1 #/ epoch_duration[k]
            #cat(it, df_it$age[j], fzz[j], i, x[df_it$age[j], fzz[j], i], "\n")
        }
    }
    
    # normalize by number of events in iteration per time interval, and by time interval length
    x_norm = x
    for (i in 1:n_epochs) {
        for (j in 1:n_it) {
            n_evt = sum(x_norm[i,,j]);# print(n_evt)
            if (n_evt > 0) {
                x_norm[i,,j] = x_norm[i,,j] / n_evt
            }
        }
    }

    x_norm = x_norm[,match(match_lbl,lbl),]
    
    # then summarize as posteriors
    ret = list(lower80=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]), 
               upper80=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]),
               lower95=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]), 
               upper95=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]),
               mean=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]))
    
    for (i in 1:dim(x_norm)[1]) {
        for (j in 1:dim(x_norm)[2]) {
            hpd80 = hdi(x_norm[i,j,], ci=0.80)
            hpd95 = hdi(x_norm[i,j,], ci=0.95)
            if (hpd95$CI_high == 0.0) hpd95$CI_high = 0.01
            if (hpd80$CI_high == 0.0) hpd80$CI_high = 0.01
            ret$lower80[i,j] = hpd80$CI_low
            ret$upper80[i,j] = hpd80$CI_high
            ret$lower95[i,j] = hpd95$CI_low
            ret$upper95[i,j] = hpd95$CI_high
            ret$mean[i,j] = mean(x_norm[i,j,])
        }
    }
    
    ret$x = x_norm
    
    # format for ggplot
    colnames(ret$lower80)=match_lbl
    colnames(ret$upper80)=match_lbl
    colnames(ret$lower95)=match_lbl
    colnames(ret$upper95)=match_lbl
    colnames(ret$mean)=match_lbl
    
    bb1 = ret
    bbm1 = melt(bb1$mean); bbm1$stat="mean"
    bbu951 = melt(bb1$upper95); bbu951$stat="upper95"
    bbl951 = melt(bb1$lower95); bbl951$stat="lower95"
    bbu801 = melt(bb1$upper80); bbu801$stat="upper80"
    bbl801 = melt(bb1$lower80); bbl801$stat="lower80"
    
    bba1 = bbm1
    bba1$upper95 = bbu951$value
    bba1$lower95 = bbl951$value
    bba1$upper80 = bbu801$value
    bba1$lower80 = bbl801$value
    
    ret = bba1
    
    ret$Var2 = factor(ret$Var2, ordered=T, levels=match_lbl)
    
    return(ret)
}

make_triplet_lstt_epoch = function(df, epoch_times, match_lbl) {
    

    iterations = sort(unique(df$iteration))
    n_it = length(iterations)
    df = df[ df$S1>0, ]
    lbl = get_first_lbl()
    n_events = length(lbl)
    n_epochs = length(epoch_times)
    df$age = ceiling(df$age)
    df$epoch = find_epoch(df$age)
    
    # matrix of (times x event_types x iterations)
    # we then use this to compute HPDs for events
    x = array(0, dim=c(n_epochs, 6, n_it))
    dimnames(x) = list(NULL,lbl,NULL)
    
    for (i in 1:n_it) {
        it = iterations[i]
        df_it = df[df$iteration==it,]
        
        zz = apply(df_it[,c("S1","S2","S3")], 1, function(z) { make_tx3_type( z[1], z[2], z[3]) } )
        fzz = factor(zz,levels=lbl)
        for (j in 1:length(zz)) {
            k = df_it$epoch[j]
            x[k, fzz[j], i] = x[k, fzz[j], i] + 1 #/ epoch_duration[k]
        }
    }
    
    # normalize by number of events in iteration per time interval, and by time interval length
    x_norm = x
    for (i in 1:n_epochs) {
        for (j in 1:n_it) {
            n_evt = sum(x_norm[i,,j])
            if (n_evt > 0) {
                x_norm[i,,j] = x_norm[i,,j] / n_evt
            }
        }
    }

    # then summarize as posteriors
    ret = list(lower80=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]), 
               upper80=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]),
               lower95=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]), 
               upper95=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]),
               mean=matrix(0, nrow=dim(x_norm)[1], ncol=dim(x_norm)[2]))
    
    for (i in 1:dim(x_norm)[1]) {
        for (j in 1:dim(x_norm)[2]) {
            hpd80 = hdi(x_norm[i,j,], ci=0.80)
            hpd95 = hdi(x_norm[i,j,], ci=0.95)
            if (hpd95$CI_high == 0.0) hpd95$CI_high = 0.01
            if (hpd80$CI_high == 0.0) hpd80$CI_high = 0.01
            ret$lower80[i,j] = hpd80$CI_low
            ret$upper80[i,j] = hpd80$CI_high
            ret$lower95[i,j] = hpd95$CI_low
            ret$upper95[i,j] = hpd95$CI_high
            ret$mean[i,j] = mean(x_norm[i,j,])
        }
    }
    
    ret$x = x_norm
    
    # format for ggplot
    colnames(ret$lower80)=lbl
    colnames(ret$upper80)=lbl
    colnames(ret$lower95)=lbl
    colnames(ret$upper95)=lbl
    colnames(ret$mean)=lbl
    
    bb1 = ret
    bbm1 = melt(bb1$mean); bbm1$stat="mean"
    bbu951 = melt(bb1$upper95); bbu951$stat="upper95"
    bbl951 = melt(bb1$lower95); bbl951$stat="lower95"
    bbu801 = melt(bb1$upper80); bbu801$stat="upper80"
    bbl801 = melt(bb1$lower80); bbl801$stat="lower80"
    
    bba1 = bbm1
    bba1$upper95 = bbu951$value
    bba1$lower95 = bbl951$value
    bba1$upper80 = bbu801$value
    bba1$lower80 = bbl801$value
    
    ret = bba1
    ret$Var2 = factor( as.vector(ret$Var2), ordered=T, levels=match_lbl)
    return(ret)
}