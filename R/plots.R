
## Plot for each iteration of the R² of the linear model between Alpha.XX and levoglucosan

levo.plot <- function(wd, itn, bab, baf, baw, bcp, fcp, wcp, abs, afs, aws, r2b, r2f, r2w) {

    png( paste0(wd, "plots/preplots/", "iter_", itn, ".png"), units = "cm", width = 18, height = 16, res = 150 )

    plot( bab, bcp [as.character(bab)], type = "n",
         xlim = c(.8,5.2), ylim = c(0, 1),
         axes = FALSE,
         xlab = bquote(Parameter~value~"for"~alpha[BC]~and~alpha[FF]),
         ylab = bquote(R^2~of~linear~regression~with~Levoglucosan), mgp = c(2.5,1,0) )

    lines(1:5, r2b, col = "deepskyblue3", lty = 2, lwd = 0.75, pch = 19)
    points(match(bab, abs), bcp[as.character(bab)], col = "red", pch = 20)

    lines(1:5, r2f, col = "forestgreen", lty = 2, lwd = .75, pch = 19)
    points(match(baf, afs), fcp[as.character(baf)], col = "red", pch = 20)

    ## plotting WB anyway with FF.set x axis: the axis label is fixed later
    lines(1:5, r2w, col = "darkorange", lty = 2, lwd = .75, pch = 19)
    points(match(baw, aws), wcp[match(baw, aws)], col = "red", pch = 20)

    axis(side = 2, lwd = 0, lwd.ticks = 1)
    axis(side = 1, at = 1:5, labels = sprintf("%.3f", afs), lwd = 0, lwd.ticks = 1)
    axis(side = 3, at = 1:5, labels = sprintf("%.3f", aws), lwd = 0, lwd.ticks = 1)
    mtext(bquote(Parameter~value~"for"~alpha[WB]), side = 3, line = 2.2)

    abline(v = 1:5, lty = 2, lwd = .5)

    legend("bottomleft", legend = c(bquote(Best~alpha), bquote(alpha[BC]), bquote(alpha[FF]), bquote(alpha[WB])), lty = c(NA, rep(2,3)), pch = c(20, NA, NA, NA), col = c("red", "deepskyblue3", "forestgreen", "darkorange"), bty = "n")

    graphics::box(bty = "o")
    
    dev.off()

}

fitplot <- function(wd, sample, babc, baff, bawb) {
    
    png( paste0(wd, "plots/fit_plots/", sample$name[1], ".png"), units = "cm", width = 18, height = 16, res = 150 )
    
    par(mar = c(4.5, 4, .5, .5))

    plot(sample$wavelength, sample$abs, pch = 19, cex = .3, ylim = range(0, sample$abs + sample$u.abs, sample$abs + sample$u.abs, na.rm = TRUE),
         xlab = "Wavelength (nm)", ylab = bquote(Absorption~coefficient~b[abs]~(Mm^-1)),
         mgp = c(2.5,1,0), type = "n")       

    with(sample, segments(wavelength - u.wavelength, abs, wavelength + u.wavelength, abs))
    with(sample, segments(wavelength, abs - u.abs, wavelength, abs + u.abs))
    points(sample$wavelength, sample$abs, pch = 19, cex = .5)
    with(sample, lines(wavelength, a * wavelength ^ (-babc) + b * wavelength ^ (-alpha.brc[1]), col = "red", lwd = .75, lty = 1) ) ## BC + BrC fit
    with(sample, lines(wavelength, a * wavelength ^ (-babc), col = "black", lty = 1, lwd = .75)) ## BC fit
    with(sample, lines(wavelength, b * wavelength ^ (-alpha.brc[1]), col = "darkorange", lty = 1, lwd = .75))  ## BrC fit
    with(sample, lines(wavelength, a.p * wavelength ^ (-baff) + b.p * wavelength ^(-bawb), col = "blue", lty = 2, lwd = .75)) ## FF + WB fit
    with(sample, lines(wavelength, a.p * wavelength ^ (-baff), col = "magenta4", lty = 2, lwd = .75)) ## FF fit
    with(sample, lines(wavelength, b.p * wavelength ^(-bawb), col = "forestgreen",  lty = 2, lwd = .75))  ## WB fit
    legend("topright", legend = c("Sample data", "BC+BrC fit", "BC share", "BrC share", "FF+WB fit", "FF share", "WB share"),
           lty = c(NA, rep(1, 3), rep(2, 3)), pch = c(19, rep(NA, 6)),
           col = c("black", "red", "black", "darkorange", "blue", "magenta4", "forestgreen"),
           pt.cex = .5, bty = "n")
    abline(v = seq(0,1000, 100), lwd = .25, lty = 2, col = gray(.75))

    dev.off()
}

opt.app.timeseries.plot <- function(wd, data, lambda.short, lambda.long) {

    png( paste0(wd, "plots/opt_app_plots/short_lambda.png"), units = "cm", width = 20, height = 16, res = 150 )
    
    dum <- dplyr::bind_rows(data) %>%
        dplyr::filter(wavelength == lambda.short)

    par(mar = c(4.5, 4, .5, .5))
    with(dum, plot(1:length(unique(dum$name)), bc.ff, col = "black", type = "o", lwd = .75, pch = 19, cex = .75,
                   ylim = range(bc.ff, brc, bc.wb, na.rm = TRUE),
                   xaxt = "n", xlab = "",
                   ylab = bquote(Absorption~coefficient~"@"~.(lambda.short)~nm~(Mm^-1)),
                   mgp = c(2.5,1,0))
         )  ## BC FF

    with(dum, lines(1:length(unique(dum$name)), bc.wb, col = "forestgreen", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## WB

    with(dum, lines(1:length(unique(dum$name)), brc, col = "darkorange", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## brc

    axis(1, at = 1:length(unique(dum$name)), labels = unique(dum$name), las = 2)
    grid()

    legend("topright", legend = c(bquote(BC[FF]~"@"~.(lambda.short)~nm), c(bquote(BC[WB]~"@"~.(lambda.short)~nm)), "BrC"),
           lty = 1, pch = 19, col = c("black", "forestgreen", "darkorange"),
           pt.cex = .5, bty = "n")

    mtext("Sample name", 1, 3.2)

    dev.off()

    ## long wavelength
    png( paste0(wd, "plots/opt_app_plots/long_lambda.png"), units = "cm", width = 20, height = 16, res = 150 )
    
    dum <- dplyr::bind_rows(data) %>%
        dplyr::filter(wavelength == lambda.long)

    par(mar = c(4.5, 4, .5, .5))
    with(dum, plot(1:length(unique(dum$name)), bc.ff, col = "black", type = "o", lwd = .75, pch = 19, cex = .75,
                   ylim = range(bc.ff, brc, bc.wb, na.rm = TRUE),
                   xaxt = "n", xlab = "",
                   ylab = bquote(Absorption~coefficient~"@"~.(lambda.long)~nm~(Mm^-1)),
                   mgp = c(2.5,1,0))
         )  ## BC FF


    grid(nx = NULL, ny = NULL, col = "lightgray", lty = 2, lwd = .75)
    
    with(dum, lines(1:length(unique(dum$name)), bc.wb, col = "forestgreen", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## WB

    with(dum, lines(1:length(unique(dum$name)), brc, col = "darkorange", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## brc

    axis(1, at = 1:length(unique(dum$name)), labels = unique(dum$name), las = 2)

    legend("topright", legend = c(bquote(BC[FF]~"@"~.(lambda.long)~nm), c(bquote(BC[WB]~"@"~.(lambda.long)~nm)), "BrC"),
           lty = 1, pch = 19, col = c("black", "forestgreen", "darkorange"),
           pt.cex = .5, bty = "n")

    mtext("Sample name", 1, 3.2)

    dev.off()

}


mass.app.timeseries.plot <- function(wd, data) {

    png( paste0(wd, "plots/mass_app_plots/mass_app_time_series.png"), units = "cm", width = 20, height = 16, res = 150 )
    
    dum <- dplyr::bind_rows(data) %>%
        dplyr::select(matches("ec.ff|oc.ff|oc.nc|ec.wb|oc.wb|name")) %>% ## keep columns needed for plot
        dplyr::distinct(.) ## remove duplicates

    par(mar = c(4.5, 4, .5, .5))
    with(dum, plot(1:length(unique(dum$name)), ec.ff, col = "black", type = "o", lwd = .75, pch = 19, cex = .75,
                   ylim = range(ec.ff, oc.ff, ec.wb, oc.wb, oc.nc, na.rm = TRUE),
                   xaxt = "n", xlab = "",
                   ylab = bquote(Concentration~(mu*g^-3)),
                   mgp = c(2.5,1,0))
         )  ## EC_FF


    grid(nx = NULL, ny = NULL, col = "lightgray", lty = 2, lwd = .75)
    
    with(dum, lines(1:length(unique(dum$name)), ec.wb, col = "forestgreen", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## EC_WB

    with(dum, lines(1:length(unique(dum$name)), oc.ff, col = "darkorange", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## brc

    with(dum, lines(1:length(unique(dum$name)), oc.wb, col = "red", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## brc

    with(dum, lines(1:length(unique(dum$name)), oc.nc, col = "deepskyblue3", lty = 1, lwd = .75, cex = .75, pch = 19, type = "o")) ## brc

    axis(1, at = 1:length(unique(dum$name)), labels = unique(dum$name), las = 2)

    legend("topright", legend = c(bquote(EC[FF]), bquote(EC[WB]), bquote(OC[FF]), bquote(OC[WB]), bquote(OC[NC])),
           lty = 1, pch = 19, col = c("black", "forestgreen", "darkorange", "red", "deepskyblue3"),
           pt.cex = .5, bty = "n")

    mtext("Sample name", 1, 3.2)

    dev.off()

}


swipe.plot <- function(wd, sabs, sabss, sabcs){

    png( paste0(wd, "plots/fit_plots/swipe_brc.png"), units = "cm", width = 18, height = 16, res = 150 )

    par(mar = c(4.5, 4, .5, .5))

    plot(spline(sabcs, sabs), ylim = range(sabs + sabss, sabs - sabss, na.rm = TRUE),
         xlab = bquote(alpha[BC]), ylab = bquote(alpha[BrC]~average~across~all~samples),
         type = "n", col = "red", mgp = c(2.5, 1, 0))

    ## Adding grid
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = 2, lwd = .75)

    spline(rev(sabcs), rev(sabs - sabss)) %>%
        data.frame(.) %>%
        arrange(desc(x)) %>%        
        rbind.data.frame(spline(sabcs, c(sabs + sabss)), .) %>%
        polygon(., col = rgb(matrix(col2rgb("red")/255,ncol=3), alpha = 0.5, maxColorValue = 1), border = rgb(matrix(col2rgb("red")/255, ncol = 3), alpha = 0.5, maxColorValue = 1))

    lines(spline(sabcs, sabs), ylim = range(sabs + sabss, sabs - sabss, na.rm = TRUE),
         type = "l", col = "red")
  
    legend("topleft", legend = c(bquote(alpha[BrC]~average), bquote(alpha[BrC]~1*sigma)),
           lty = 1, lwd = c(1, 6), col = c("red", rgb(matrix(col2rgb("red")/255,ncol=3), alpha = 0.5, maxColorValue = 1)), bty = "n")
    
    dev.off()
    
}

brc.variability.plot <- function(wd, data){

    dum <- dplyr::bind_rows(data) %>%
        dplyr::select(matches("name|alpha.brc")) %>% ## keep columns needed for plot
        dplyr::distinct(.) ## remove duplicates

    png( paste0(wd, "plots/fit_plots/brc_variability.png"), units = "cm", width = 18, height = 16, res = 150 )

    par(mar = c(4.5, 4, .5, .5))
    
    plot(1:nrow(dum), dum$alpha.brc, ylim = range(dum$alpha.brc + dum$u.alpha.brc, dum$alpha.brc - dum$u.alpha.brc, na.rm = TRUE), xlab = "", ylab = bquote(alpha[BrC]~"\U00B1"~1*sigma), type = "p", col = "red", mgp = c(2.5, 1, 0), pch = 19, xaxt = "n")

    segments(1:length(dum), dum$alpha.brc - dum$u.alpha.brc, 1:length(dum), dum$alpha.brc + dum$u.alpha.brc, col = "red", lwd = .75)
    
    axis(1, at = 1:length(unique(dum$name)), labels = unique(dum$name), las = 2)
    mtext("Sample name", 1, 3.2)

    ## Adding grid
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = 2, lwd = .75)
    
    dev.off()
}