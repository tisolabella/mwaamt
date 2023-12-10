library(jsonlite)
library(dplyr)
library(minpack.lm)
library(nlme)

source("datahandle.R")
source("set_resolution.R")

## ~~~~~~~~~~~~~~~~ ##
## Wrapper function ##
## ~~~~~~~~~~~~~~~~ ##

config.file <- "files/m_paper_config.json"

run <- function(configuration.file.path) {
    
    ## Get the starting date and time
    start.time <- Sys.time()
    
    ## Suppress warnings
    ## options(warn=-1)
      
    ## ~~~~~~~~~~~~~~ ##
    ## Read the input ##
    ## ~~~~~~~~~~~~~~ ##    
    
    message("\n---> Opening configuration file...\n")
    tryCatch({
        cfg <- jsonlite::fromJSON(config.file)
        names(cfg) <- gsub(" ", ".", names(cfg), fixed = TRUE) ## replace space with dot from the names in the config file
    }, error = function(e) {
        print(e)
    })
    
    ## Read the data with the method from datahandle.py
    data <- data.read.mw(config.file)
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Preliminary fitting and correlations ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    ## Check for booked presets:
    
    tryCatch({
        levo.booked <- 'levoglucosan' %in% cfg$presets
    }, error = function(e) {
        message(paste("Error: keyword not present in the configuration file", e$message))
    })
    
    ## Fill in with other presets
    ## TODO implement user-defined fits and correlations (for the GUI)
    ## The values in the cfg file are the default to be used if there is
    ## no preliminary analysis
    
    alpha.BC <- cfg$alpha.BC
    alpha.ff <- cfg$alpha.FF
    alpha.wb <- cfg$alpha.WB
    best.alpha.BC <- alpha.BC
    best.alpha.ff <- alpha.ff
    best.alpha.wb <- alpha.wb
    
    ## Carry out the Levoglucosan correlation adjustment
    
    if (levo.booked == TRUE) {
        message(paste('---> Performing correlation maximisation with "levoglucosan" preset...\n'))

        ## iteration start
        for (iter.number in seq_len(cfg$iterations)) {
            message(paste("\n-----> Iteration number:", iter.number, "\n"))

            ## Set the search resolution
            res <- set.resolution( c(best.alpha.BC, best.alpha.ff, best.alpha.wb), iteration = iter.number - 1)
            names(res) <- c("alpha.BC","alpha.ff", "alpha.wb")
            alpha.BC.set <- res$alpha.BC
            alpha.ff.set <- res$alpha.ff
            alpha.wb.set <- res$alpha.wb
            rm(res)

            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Do the alpha.BC iteration ##
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##

            ## correlations
            r2.alpha.bc <- NULL
            BC.correlation.pairs <- NULL

            for (alpha.BC in alpha.BC.set) {
                
                alpha.BC <- round(alpha.BC, 3) ## rounding alpha.bc
                ## List to save the BrC(lambda.short) and levo
                BrC.set <- NULL
                levo.set <- NULL
                
                ## Loop over all the samples
                for (sample in data) {

                    BrC <- tryCatch({
                        ## fitting equation 1 of the paper
                        fitres <- nlsLM(abs ~  A * wavelength^(-alpha.BC) + B * wavelength^(-alpha.2),
                                        data = sample, 
                                        start = list("A" = 1e3, "B" = 1e10, "alpha.2" = 3),
                                        lower = c(1e2, 1e2, 1),
                                        upper = c(1e15, 1e15, 10),
                                        weights = 1/(sample$u.abs ^ 2),
                                        control = list(maxiter = 1000)
                                        )
                        
                        B <- coef(fitres)["B"] %>%
                            round(.,0)
                        alpha.BrC <- coef(fitres)["alpha.2"] %>%
                            round(.,3)
                        lambda.short <- min(sample$wavelength)

                        ## equation 3: estimate BrC
                        BrC <- as.numeric(B * lambda.short ^ (-alpha.BrC))
                        
                    }, error = function(e) {
                        message(paste("!! Fit error for alpha.BC =" , round(alpha.BC, 4), "in sample", sample$name[1], ":", e$message))
                        return(NA)
                    })

                    BrC.set <- c(BrC.set, BrC)
                    levo.set <- c(levo.set, sample$add.meas[1])
                    
                } ## end of the sample loop
                
                ## Calculate regression and append R^2
                r2 <- tryCatch({
                    regression.res <- lm(BrC.set ~ levo.set)
                    summary(regression.res)$r.squared
                }, error = function(e) {
                    message(paste("!! Regression error for for alpha.BC =", round(alpha.BC,4), "in sample", sample$name[1], ":", e$message))
                    return(NA)
                })
                
                r2.alpha.bc <- c(r2.alpha.bc, r2)
                BC.correlation.pairs <- c(BC.correlation.pairs, r2)
                
            } ## end of the BC scan loop

            ## add names to BC correlation pairs
            names(BC.correlation.pairs) <- alpha.BC.set

            ## For stability, use the best only if its significantly
            ## better than the previous best:

            if (max(r2.alpha.bc, na.rm = TRUE) - BC.correlation.pairs[as.character(best.alpha.BC)] < cfg$threshold) {
                ## Do not change the value
            } else {
                ## 
                best.alpha.BC <- alpha.BC.set[which.max(r2.alpha.bc)] %>%
                    round(., 3)
            }

            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Do the alpha.ff iteration ##
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##

            ## List for the correlations
            r2.alpha.ff <- NULL
            ff.correlation.pairs <- NULL
            
            for (alpha.ff in alpha.ff.set) {
                
                alpha.ff <- round(alpha.ff, 3)
                ## List to save the BC.wb(lambda.long) and levo
                BC.wb.set <- NULL
                levo.set <- NULL
                
                ## Iterate over all the samples
                for (sample in data) {

                    BC.wb <- tryCatch({
                        ## fitting equation 1 of the paper using best.alpha.BC
                        typefit <- nlsLM(abs ~  A * wavelength^(-best.alpha.BC) + B * wavelength^(-alpha.2),
                                         data = sample, 
                                         start = list("A" = 1e3, "B" = 1e10, "alpha.2" = 3),
                                         lower = c(1e2, 1e2, 1),
                                         upper = c(1e15, 1e15, 10),
                                         weights = 1/sample$u.abs^2,
                                         control = list(maxiter = 1000)
                                         )

                        ## fitting equation 2 of the paper using best.alpha.BC
                        sourcefit <- nlsLM(abs ~  A.p * wavelength^(-alpha.ff) + B.p * wavelength^(-alpha.wb),
                                           data = sample, 
                                           start = list("A.p" = 1e3, "B.p" = 1e10),
                                           lower = c(1e2, 1e2),
                                           upper = c(1e15, 1e15),
                                           weights = 1/sample$u.abs^2,
                                           control = list(maxiter = 1000)
                                           )
                        
                        A <- coef(typefit)["A"] %>%
                            round(., 0)
                        A.p <- coef(sourcefit)["A.p"] %>%
                            round(., 0)
                        
                        ## Apportion BC WB at the longest wavelength using eq. 3 of the paper
                        lambda.long <- max(sample$wavelength)
                        BC.wb <- (A - A.p) * lambda.long ^ (-best.alpha.BC)
                        
                    }, error = function(e) {
                        message(paste("!! Fit error for alpha.BC =" , round(alpha.BC, 4), "in sample", sample$name[1], ":", e$message))
                        return(NA)
                    })

                    BC.wb.set <- c(BC.wb.set, BC.wb)
                    levo.set <- c(levo.set, sample$add.meas[1])
                    
                } ## end of the sample loop

                ## Calculate regression and append R^2
                r2 <- tryCatch({
                    regression.res <- lm(BC.wb.set ~ levo.set)
                    summary(regression.res)$r.squared
                }, error = function(e) {
                    message(paste("!! Regression error for for alpha.FF =" ,alpha.ff, "in sample", sample$name[1], ":", e$message))
                    return(NA)
                })
                r2.alpha.ff <- c(r2.alpha.ff, r2)
                ff.correlation.pairs <- c(ff.correlation.pairs, r2)
                
            } ## end of the FF scan loop
            
            ## add names to BC correlation pairs
            names(ff.correlation.pairs) <- alpha.ff.set

            if (max(r2.alpha.ff, na.rm = TRUE) - ff.correlation.pairs[as.character(best.alpha.ff)] < cfg$threshold) {
                ## Do not change the value
            } else {
                best.alpha.ff <- alpha.ff.set[which.max(r2.alpha.ff)] %>%
                    round(., 3)
            }
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Do the alpha.wb iteration ##
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##

            ## List for the correlations
            r2.alpha.wb <- NULL
            wb.correlation.pairs <- NULL
            
            for (alpha.wb in alpha.wb.set) {

                alpha.wb <- round(alpha.wb, 3)
                
                ## List to save the BC.wb(lambda.long) and levo
                BC.wb.set <- numeric()
                levo.set <- numeric()
                
                ## Iterate over all the samples
                for (sample in data) {

                    BC.wb <- tryCatch({
                        ## fitting equation 1 of the paper using best.alpha.BC
                        typefit <- nlsLM(abs ~  A * wavelength^(-best.alpha.BC) + B * wavelength^(-alpha.2),
                                         data = sample, 
                                         start = list("A" = 1e3, "B" = 1e10, "alpha.2" = 3),
                                         lower = c(1e2, 1e2, 1),
                                         upper = c(1e15, 1e15, 10),
                                         weights = 1/sample$u.abs^2,
                                         control = list(maxiter = 1000)
                                         )

                        ## fitting equation 2 of the paper using best.alpha.BC
                        sourcefit <- nlsLM(abs ~  A.p * wavelength^(-best.alpha.ff) + B.p * wavelength^(-alpha.wb),
                                           data = sample, 
                                           start = list("A.p" = 1e3, "B.p" = 1e10),
                                           lower = c(1e2, 1e2),
                                           upper = c(1e15, 1e15),
                                           weights = 1/sample$u.abs^2,
                                           control = list(maxiter = 1000)
                                           )
                        
                        A <- coef(typefit)["A"] %>%
                            round(., 0)
                        A.p <- coef(sourcefit)["A.p"] %>%
                            round(., 0)

                        ## Apportion BC WB at the longest wavelength using eq. 3 of the paper
                        lambda.long <- max(sample$wavelength)
                        BC.wb <- (A - A.p) * lambda.long ^ (-best.alpha.BC)

                    }, error = function(e) {
                        message(paste("!! Fit error for alpha.ff =" ,alpha.ff, "in sample", sample$name[1], ":", e$message))
                        return(NA)
                    })

                    BC.wb.set <- c(BC.wb.set, BC.wb)
                    levo.set <- c(levo.set, sample$add.meas[1])
                    
                } ## end of the sample loop
                
                ## Calculate regression and append R^2
                r2 <- tryCatch({
                    regression.res <- lm(BC.wb.set ~ levo.set)
                    summary(regression.res)$r.squared
                }, error = function(e) {
                    message(paste("!! Regression error for for alpha.WB =" ,alpha.wb, "in sample", sample$name[1], ":", e$message))
                    return(NA)
                })
                r2.alpha.wb <- c(r2.alpha.wb, r2)
                wb.correlation.pairs <- c(wb.correlation.pairs, r2)
                
            } ## end of the WB scan loop

            ## add names to BC correlation pairs
            names(wb.correlation.pairs) <- alpha.wb.set

            if (max(r2.alpha.wb, na.rm = TRUE) - wb.correlation.pairs[as.character(best.alpha.wb)] < cfg$threshold) {
                ## Do not change the value
            } else {
                best.alpha.wb <- alpha.wb.set[which.max(r2.alpha.wb)] %>%
                    round(., 3)
            }

            ## From now on, use alpha.XX for the best values
            ## or the default values if no optimization was done
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## Save plots for the correlation vs parameter values ##
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

            ##     tryCatch({
            ##         if (cfg$plots) {
            ##             cat("\t\tsaving plots\n")
            ##             png(paste0(cfg$`working directory`, "/plots/preplots/iter", iter.number + 1, ".png"), width = 800, height = 600)
            ##             par(mfrow=c(1,1))
            ##             plot(best.alpha.BC, BC.correlation.pairs[as.character(best.alpha.BC)], type = "n", xlim = c(min(alpha.BC.set), max(alpha.BC.set)), ylim = c(0, 1), xlab = "Parameter value for alpha.BC and alpha.ff", ylab = "Levoglucosan analysis R^2")
            ##             lines(alpha.BC.set, r2.alpha.bc, col = "black", lty = 1, lwd = 2)
            ##             points(best.alpha.BC, BC.correlation.pairs[as.character(best.alpha.BC)], col = "red", pch = 20)
            ##             lines(alpha.ff.set, r2.alpha.ff, col = "red", lty = 2, lwd = 2)
            ##             points(best.alpha.ff, ff.correlation.pairs[as.character(best.alpha.ff)], col = "red", pch = 20)
            ##             axis(side = 1, labels = FALSE, at = c(best.alpha.BC, best.alpha.ff), tck = -0.02)
            ##             mtext(side = 1, text = c("alpha.BC", "alpha.ff"), line = 1, at = c(best.alpha.BC, best.alpha.ff))
            ##             axis(side = 2, las = 2, at = seq(0, 1, by = 0.1))
            ##             dev.off()
            ##         }
            ##     }, error = function(e) {
            ##         print(paste(missing.keyword, e$message))
            ##     })
            ## }
            
        }
        
        alpha.BC <- best.alpha.BC
        alpha.ff <- best.alpha.ff
        alpha.wb <- best.alpha.wb
        message(paste("\nThe best parameters for the the correlation with levoglucosan are",
                      "alpha.BC =", round(alpha.BC, 3),
                      ", alpha.ff =", round(alpha.ff, 3),
                      ", alpha.wb =", round(alpha.wb, 3), "\n"))
    }

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ALPHA BrC variation with ALPHA BC (swipe) ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    
    if (cfg[['alpha.bc.swipe']] == TRUE) {
        message("\n ----> Performing alpha bc swipe...\n")
        sw.alpha.BC.set <- seq(0.8, 1.2, length.out = 20)

        sw.alpha.BrC.set <- NULL
        sw.alpha.BrC.sd.set <- NULL

        ## alpha.bc loop
        for (abc in sw.alpha.BC.set) {

            BrC <- NULL

            ## sample loop
            for (sample in data) {

                fitres <- tryCatch({
                    ## fitting equation 1 of the paper
                    nlsLM(abs ~  A * wavelength^(-abc) + B * wavelength^(-alpha.brc),
                          data = sample, 
                          start = list("A" = 1e3, "B" = 1e10, "alpha.brc" = 3),
                          lower = c(1e2, 1e2, 1),
                          upper = c(1e15, 1e15, 10),
                          weights = 1/sample$u.abs^2,
                          control = list(maxiter = 1e3)
                          ) %>%
                        coefficients() %>%
                        extract2(3)
                }, error = function(e) {
                    message(paste("!! Fit error during alpha.BC swipe for alpha.BC =" , round(alpha.BC,4), "in sample", sample$name[1], ":", e$message))
                    return(NA)
                })

                BrC <- c(BrC, fitres)
            } ## end of sample loop

            sw.alpha.BrC.set <- c( sw.alpha.BrC.set, mean(BrC, na.rm = TRUE) )
            sw.alpha.BrC.sd.set <- c( sw.alpha.BrC.sd.set, sd(BrC, na.rm = TRUE) )
            
        } ## end of alpha.BC loop

        rm(fitres, BrC, sample)
    }


    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## FITTING PROCEDURE WITH OPTIMIZED PARAMETERS ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    message('\n ---> Fitting the experimental data...\n')

    failed.fit.count <- 0
    failed.fit <- list()

    for (sample in data) {
        
        ## AAE fit
        tryCatch({

            aae.fitres <- nlsLM(abs ~  A * wavelength^(-aae),
                                data = sample, 
                                start = list("A" = 1e5, "aae" = 1),
                                lower = c(1e2, .5),
                                upper = c(1e15, 3),
                                weights = 1/sample$u.abs^2,
                                control = list(maxiter = 1e3)
                                )
            
            ## Save first aae
            data[[sample$name[1]]]$aae <- coef(aae.fitres)["aae"]
            data[[sample$name[1]]]$u.aae <- sqrt(vcov(aae.fitres)[2, 2])
            
            ## Repeat fit with fixed aae to decrease uncertainty
            aae.fitres.sec <- nlsLM(abs ~  A * wavelength^(- coef(aae.fitres)["aae"]),
                                    data = sample,
                                    start = list("A" = 1e5),
                                    lower = c(1e2),
                                    upper = c(1e15),
                                    weights = 1/sample$u.abs^2,
                                    control = list(maxiter = 1e3)
                                    )

            ## Save the A scale parameter
            data[[sample$name[1]]]$scale <- coef(aae.fitres.sec)["A"]
            data[[sample$name[1]]]$u.scale <- sqrt(vcov(aae.fitres.sec)[1, 1])
            
            ## Calculates the reduced chisquared for this fit
            ## in python non è esplicitato AAE per la stima del modellato sotto
            data[[sample$name[1]]]$chisq.aae.fit <- sum(residuals(aae.fitres.sec)^2 / sample$u.abs^2) / (length(sample$wavelength) - 2)
            
        }, error = function(e) {
            message(paste('!! Fit error for AAE in sample', sample$name[1], ":", e$message, '\n'))
            failed.fit.count <<- failed.fit.count + 1
            failed.fit <- c(failed.fit, list(c('aae', paste('sample', sample$name))))
        })
        
        ## Components fit
        tryCatch({
            ## fitting equation 1 of the paper
            ## fix alpha.bc to the best value
            fitres <- nlsLM(abs ~  A * wavelength^(-alpha.BC) + B * wavelength^(-alpha.2),
                            data = sample, 
                            start = list("A" = 1e3, "B" = 1e10, "alpha.2" = 3),
                            lower = c(1e2, 1e2, 1),
                            upper = c(1e15, 1e15, 10),
                            weights = 1/sample$u.abs^2,
                            control = list(maxiter = 1e3)
                            )
            
            ## Save alpha.BrC
            data[[sample$name[1]]]$alpha.brc <- coef(fitres)[3] %>%
                round(., 3)
            data[[sample$name[1]]]$u.alpha.brc <- sqrt(vcov(fitres)[3, 3]) %>%
                round(., 3)
            
            ## repeat components fit with fixed alpha.BrC to reduce the uncertainty 
            fitres.rep <- nlsLM(abs ~  A * wavelength^(-alpha.BC) + B * wavelength^(-coef(fitres)[3]),
                                data = sample, 
                                start = list("A" = 1e3, "B" = 1e10),
                                lower = c(1e2, 1e2),
                                upper = c(1e15, 1e15),
                                weights = 1/sample$u.abs^2,
                                control = list(maxiter = 1e3)
                                )
            
            data[[sample$name[1]]]$a <- coef(fitres.rep)[1]
            data[[sample$name[1]]]$u.a <- sqrt(vcov(fitres.rep)[1, 1])
            data[[sample$name[1]]]$b <- coef(fitres.rep)[2]
            data[[sample$name[1]]]$u.b <- sqrt(vcov(fitres.rep)[2, 2])
            
            ## Calculates the reduced chisquare for this fit
            data[[sample$name[1]]]$chisq.type.fit <- sum(residuals(fitres.rep)^2 / sample$u.abs^2) / (length(sample$wavelength) - 3)
            
        }, error = function(e) {
            message(paste('!! Fit error for BrC in sample', sample$name[1], ":", e$message, '\n'))
            failed.fit.count <<- failed.fit.count + 1
            failed.fit <- c(failed.fit, list(c('type', paste('sample', sample$name))))
        })
        

        ## Source fit 
        tryCatch({
            ## fitting equation 2 with alpha_FF and alpha_WB fixed to their best values
            source.fit <- nlsLM(abs ~  A.p * wavelength^(-alpha.ff) + B.p * wavelength^(-alpha.wb),
                                data = sample, 
                                start = list("A.p" = 1e3, "B.p" = 1e10),
                                lower = c(1e2, 1e2),
                                upper = c(1e15, 1e15),
                                weights = 1/sample$u.abs^2,
                                control = list(maxiter = 1e3)
                                )

            data[[sample$name[1]]]$a.p <- coef(source.fit)[1] ## save A.p
            data[[sample$name[1]]]$u.a.p <- sqrt(vcov(source.fit)[1, 1])
            data[[sample$name[1]]]$b.p <- coef(source.fit)[2] ## save B.p
            data[[sample$name[1]]]$u.b.p <- sqrt(vcov(source.fit)[2, 2])
            
            ## Calculates the reduced chisquare for the source fit
            data[[sample$name[1]]]$chisq.source.fit <- sum( residuals(source.fit)^2 / sample$u.abs^2) / (length(sample$wavelength) - 3)

        }, error = function(e) {
            cat('FIT.ERROR', e, '\n')
            failed.fit.count <<- failed.fit.count + 1
            failed.fit <- c(failed.fit, list(c('source', paste('sample', sample$name[1]))))
        })
    } ## end of sample loop for apportionment

    ## ~~~~~~~~~~~~~~~~~~~~~ ##
    ## OPTICAL APPORTIONMENT ##
    ## ~~~~~~~~~~~~~~~~~~~~~ ##

    message("\n ---> Performing optical apportionment... \n")

    for (sample in data) {

        A <- sample$a[1] %>%
            round(., 0)
        B <- sample$b[1] %>%
            round(., 0)
        A.p <- sample$a.p[1] %>%
            round(., 0)
        B.p <- sample$b.p[1] %>%
            round(., 0)
        alpha.BrC <- sample$alpha.brc[1] %>%
            round(., 3)
        w <- sample$wavelength
        abs <- sample$abs
        
        ## BC wood burning
        data[[sample$name[1]]]$bc.wb <- (A - A.p) / (w ^ alpha.BC)
        data[[sample$name[1]]]$bc.wb.frac <- (A - A.p) / (w ^ alpha.BC) / abs
        ## BC fossil fuel
        data[[sample$name[1]]]$bc.ff <- A.p / (w ^ alpha.ff)
        data[[sample$name[1]]]$bc.ff.frac <- A.p / (w ^ alpha.ff) / abs
        ## BrC
        data[[sample$name[1]]]$brc <- B / (w ^ alpha.BrC)
        data[[sample$name[1]]]$brc.frac  <- B / (w ^ alpha.BrC) / abs
        
    }

    ## ~~~~~~~~~~~~~~~~~~ ##
    ## MASS APPORTIONMENT ##
    ## ~~~~~~~~~~~~~~~~~~ ##

    if (cfg$mass.appo == TRUE) {
        
        ## Check if k1 and k2 are provided
        ## NOTA: 0 è un numero, quindi è provided. Vuoi forse dire NA?
        ## non capisco questo passaggio: 
        if (cfg$k1 == 0 | cfg$k2 == 0) {
            do.fit <- TRUE
        } else {
            k1 <- cfg$k1
            k2 <- cfg$k2
            do.fit <- FALSE
        }
        
        message("\n ---> Performing mass apportionment...\n")
        
        ## EC apportionment
        ## NOTA: hai settato il EC/OC.FF.FRAC ma non lo calcoli
        
        for (sample in data) {
            
            ## get longest wavelength
            i.l <- which.max(sample$wavelength)
            
            ## Apportion EC
            data[[sample$name[1]]]$ec.ff <- sample$ec * ( sample$bc.ff[i.l] / (sample$abs[i.l] - sample$brc[i.l]) )
            data[[sample$name[1]]]$ec.wb <- sample$ec * ( sample$bc.wb[i.l] / (sample$abs[i.l] - sample$brc[i.l]) )
        }
        
        ## OC apportionment
        
        ## empty list to hold fit coeff
        k1.list <- k2.list <- NULL
        k1.x <- k1.y <- NULL
        k2.x <- k2.y <- NULL
        
        for (sample in data) {
            
            ## get longest wavelength
            i.l <- which.max(sample$wavelength)
            
            if (sample$aae[1] < cfg$AAE.high & sample$aae[1] > cfg$AAE.low) {
                ## aae of the sample is within the cfg bounds
                k1.list <- c(k1.list, sample$name[[1]])
                k1.x <- c(k1.x, sample$bc.ff[i.l])
                k1.y <- c(k1.y, sample$oc[1])
            } else {
                ## aae of the sample is outside the cfg bounds
                k2.list <- c(k2.list, sample$name[[1]])
            }
        }
        
        ## perform the linear regression fit, if needed
        if (do.fit == TRUE) {
            fit.1 <- lm(k1.y ~ k1.x)
            k1 <- coef(fit.1)[2]
            u.k1 <- summary(fit.1)$coefficients[2,"Std. Error"]
            int1 <- coef(fit.1)[1]
        }
        
        ## compute OC.FF for all samples, based on the k1 coeff
        for (sample in data) {
            i.l <- which.max(sample$wavelength)
            data[[sample$name[1]]]$oc.ff <- k1 * sample$bc.ff[i.l]
        }
        
        ## populate k2.x and k2.y for regression
        for (j in k2.list) {
            sample <- data[[j]]
            i.s <- which.min(sample$wavelength)
            k2.x <- c(k2.x, sample$brc[i.s])
            k2.y <- c(k2.y, sample$oc[1] - data[[sample$name[1]]]$oc.ff[1])
        }
        
        ## perform the linear regression fit, if needed
        if (do.fit == TRUE) {
            fit.2 <- lm(k2.y ~ k2.x)
            k2 <- coef(fit.2)[2]
            u.k2 <- summary(fit.2)$coefficients[2, "Std. Error"]
            int2 <- coef(fit.2)[1]
        }
        
        ## compute OC.wb and OC.nc for all samples, based on k2 coeff
        for (sample in data) {
            i.s <- which.min(sample$wavelength)
            data[[sample$name[1]]]$oc.wb <- k2 * sample$brc[i.s]
            data[[sample$name[1]]]$oc.nc <- sample$oc - data[[sample$name[1]]]$oc.ff - data[[sample$name[1]]]$oc.wb
        }
        
    } ## end of the mass apportionment loop
    
    ## Export fit parameters to file
    out.path <- paste0(cfg$working.directory, 'fit_results.csv')
    message(paste("\n ---> Writing fit results to", out.path, "\n"))

    write.table(data.frame(
        Name = sapply(data, function(x) x$name[1]),
        Scale = sapply(data, function(x) x$scale[1]),
        err.Scale = sapply(data, function(x) x$u.scale[1]),
        AAE = sapply(data, function(x) x$aae[1]),
        err.AAE = sapply(data, function(x) x$u.aae[1]),
        Red.chisq = sapply(data, function(x) x$chisq.aae.fit[1]),
        A = sapply(data, function(x) x$a[1]),
        err.A = sapply(data, function(x) x$u.a[1]),
        B = sapply(data, function(x) x$b[1]),
        err.B = sapply(data, function(x) x$u.b[1]),
        alpha.BrC = sapply(data, function(x) x$alpha.brc[1]),
        err.alpha.BrC = sapply(data, function(x) x$u.alpha.brc[1]),
        Red.chisq = sapply(data, function(x) x$chisq.type.fit[1]),
        A.p = sapply(data, function(x) x$a.p[1]),
        err.A.p = sapply(data, function(x) x$u.a.p[1]),
        B.p = sapply(data, function(x) x$b.p[1]),
        err.B.p = sapply(data, function(x) x$u.b.p[1]),
        Red.chisq = sapply(data, function(x) x$chisq.source.fit[1])
    ), file = out.path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ", ")

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Optical apportionment results export ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    out.path <- paste0(cfg$working.directory, 'opt_app_results.csv')
    message( paste("\n ---> Writing optical apportionment results to", out.path, "\n") )

    ## Usa il primo campione per scrivere l'intestazione
    cat('Name,', paste( paste(rep( c('BC.ff.frac', 'BC.wb.frac', 'BrC.frac', 'BC.ff', 'BC.wb', 'BrC'), each = length(data[[1]]$wavelength)),
                              rep(data[[1]]$wavelength, 6),
                              sep = "."), collapse = ","), file = out.path, append = FALSE)
    ## cat('\n,', paste(rep( c('BC.ff.frac', 'BC.wb.frac', 'BrC.frac', 'BC.ff', 'BC.wb', 'BrC'), each = length(data[[1]]$wavelength)), collapse = ","),
    ##         file = out.path, append = TRUE)
    cat('\n',file = out.path, append = TRUE)
    
    write.table(t(sapply(
        data, function(x) c(x$name[1], x$bc.ff.frac, x$bc.wb.frac, x$brc.frac,
                            x$bc.ff, x$bc.wb, x$brc))),
        out.path, row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE, sep =",")

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## mass apportionment results export ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    out.path <- paste0(cfg$working.directory, 'mass_app_results.csv')
    message(paste("\n ---> Writing mass apportionment results to", out.path, "\n"))
    
    cat(c('Name', "EC_FF", "EC_WB","OC_FF", "OC_WB","OC_NC\n"), file = out.path, append = FALSE)

    write.table(t(sapply(data, function(x) c(x$name[1], x$ec.ff[1], x$ec.wb[1], x$oc.ff[1], x$oc.wb[1], x$oc.nc[1])))
              , file = out.path, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE, sep =",")

    ## ~~~~~~~~~~~~~~ ##
    ## export logfile ##
    ## ~~~~~~~~~~~~~~ ##

    ## Get the endtime of the analysis
    end.time <- Sys.time()

    log <- NULL

    ## Prepare the lines to write
    log[1] <- paste0("Start time: ", format(start.time, "%Y-%m-%d %H:%M:%S"))
    log[2] <- paste0("End time: ", format(end.time, "%Y-%m-%d %H:%M:%S"))
    log[3] <- paste0("Input file: ", cfg$input.file)
    log[4] <- paste0("Output folder: ", cfg$working.directory)
    log[5] <- paste0("\n---------- PREPROCESSING ----------")
    log[6] <- paste0("Booked presets: ", cfg$presets)
    log[7] <- paste0("Best parameters: \n\talpha.BC = ", sprintf("%.3f", best.alpha.BC), "\talpha.FF = ", sprintf("%.3f", best.alpha.ff), "\talpha.WB = ", sprintf("%.3f", best.alpha.wb) )
    log[8] <- ifelse(cfg$plots & length(cfg$presets) > 0, paste0("Parameter optimization plots in: ", cfg$working.directory, "plots/preplots/"), "Parameter optimization plots not saved\n")
    log[9] <- paste0("\n---------- FIT ----------")
    log[10] <- paste0("Number of failed fits: ", failed.fit.count)
    log[11] <- paste0("Sample name of failed fits: ", failed.fit)
    log[12] <- paste0("Sample name of failed fits: ", failed.fit)

    log[13] <- sapply(data, function(d)  d$alpha.brc[1]) %>%
        mean(., na.rm = TRUE) %>%
        sprintf("%.4f", .) %>%
        paste0("Average alpha.BrC: ", .)

    log[14] <-sapply(data, function(d)  d$alpha.brc[1]) %>%
        sd(., na.rm = TRUE) %>%
        sprintf("%.4f", .) %>%
        paste0("Uncertainty (stdev) on alpha.BrC: ", .)

    log[15] <- ifelse(cfg$plots, paste0("Fit plots in: ", cfg$working.directory, "plots/fitplots/"), "Fit plots not saved")

    log[16] <- "\n---------- OPTICAL APPORTIONMENT ---------- "
    log[17] <- ifelse(cfg$plots, paste0("Optical apportionment plots in: ", cfg$working.directory, "plots/appoplots/"), "Optical apportionment plots not saved")

    ## Check for mass appo
    if (cfg$mass.appo == TRUE) {
        log[18] <- ifelse(do.fit, paste0("k1: ", sprintf("%.3f", k1), ", R-square: ", sprintf("%.3f", summary(fit.1)$r.squared) ), paste0("k1: ", sprintf("%.3f", k1)))
        log[19] <- ifelse(do.fit, paste0("k2: ", sprintf("%.3f", k2), ", R-square: ", sprintf("%.3f", summary(fit.2)$r.squared) ), paste0("k2: ", sprintf("%.3f", k2)))
        log[20] <- ifelse(cfg$plots == TRUE & do.fit == TRUE, paste0("Linear regression plots in: ", cfg$working.directory, "plots/mappoplots/"), "Linear regression plots not saved")
    }

    out.path <- paste0(cfg$working.directory, 'log.txt')
    message(paste("---> Writing log file in", out.path))

    cat("---------- GENERAL ----------\n", file = out.path, append = FALSE)
    write.table(log, file = paste0(cfg$working.directory, 'log.txt'), append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    


}
