.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}


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