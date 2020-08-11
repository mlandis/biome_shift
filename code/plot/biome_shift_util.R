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
    # print(x)
    # print(covered)
    # if (covered) {
    #     return("Covered")
    # } else {
    #     return("Not covered")
    # }
    # return(covered)
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