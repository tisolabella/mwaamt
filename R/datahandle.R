## IMPORTS

library(jsonlite)
library(dplyr)
library(stringr)

## NUMERICAL CONSTANTS

default.abs.unc <- 0.05

## PRINTOUT STRINGS

input.file.key.not.present <- "CONFIGURATION FILE ERROR: the key 'input file' was not found in the configuration file"
input.file.no.abs <- "INPUT FILE ERROR: the input file does not contain either a column named ABS nor a column named Babs"
mass.appo.no.ec.oc <- "WARNING: the mass.appo flag has been set to True but the EC and OC values were not provided in the input  file. The mass apportionment will not be performed."
mass.appo.no.flag <- "WARNING: EC and OC were provided but the mass.appo flag has not been set to True. Mass apportionment will not be performed."

## ## CLASS DEFINITIONS
## property <- setClass(
##     "property",
##     slots = list(
##         data.type = "character",
##         wavelength = "numeric",
##         u.wavelength = "numeric",
##         abs = "numeric",
##         u.abs = "numeric",
##         ec = "numeric",
##         oc = "numeric",
##         aae = "numeric",
##         u.aae = "numeric",
##         scale = "numeric",
##         u.scale = "numeric",
##         chisq.aae.fit = "numeric",
##         a = "numeric",
##         u.a = "numeric",
##         b = "numeric",
##         u.b = "numeric",
##         alpha.brc = "numeric",
##         u.alpha.brc = "numeric",
##         chisq.type.fit = "numeric",
##         a.p = "numeric",
##         u.a.p = "numeric",
##         b.p = "numeric",
##         u.b.p = "numeric",
##         chisq.source.fit = "numeric",
##         brc = "numeric",
##         brc.fraction = "numeric",
##         bc.ff = "numeric",
##         bc.ff.fraction = "numeric",
##         bc.wb = "numeric",
##         bc.wb.fraction = "numeric",
##         ec.ff = "numeric",
##         ec.ff.fraction = "numeric",
##         ec.wb = "numeric",
##         ec.wb.fraction = "numeric",
##         oc.ff = "numeric",
##         oc.ff.fraction = "numeric",
##         oc.wb = "numeric",
##         oc.wb.fraction = "numeric",
##         oc.nc = "numeric",
##         oc.nc.fraction = "numeric"
##     )
## )

## sample <- setClass(
##     "sample",
##     slots = list(
##         name = "character",
##         properties = "property"
##     )
## )

## DATA READING SECTION

data.read.mw <- function(configuration.file) {
    ## Read the data from an input file.
    
    ## Parameters:
    ##   config.file : str,
    ##       the path to the JSON configuration file which contains, 
    ##       among other information, the input data file path under the 
    ##       key 'input file'.
    ##
    ## Returns:
    ##   data : sequence,
    ##       a list of sample objects, each one corresponding to a physical
    ##       sample for analysis.


    nn <- c("name", "data.type","wavelength","u.wavelength","abs","u.abs","ec","oc","aae","u.aae","scale","u.scale","chisq.aae.fit","a","u.a","b","u.b","alpha.brc","u.alpha.brc","chisq.type.fit","a.p","u.a.p","b.p","u.b.p","chisq.source.fit","brc","brc.frac","bc.ff","bc.ff.frac","bc.wb","bc.wb.frac","ec.ff","ec.ff.frac","ec.wb","ec.wb.frac","oc.ff","oc.ff.frac","oc.wb","oc.wb.frac","oc.nc","oc.nc.frac","add.meas")
    
    ## Open the input data file
    config <- jsonlite::fromJSON(config.file)
    names(config) <- gsub(" ", ".", names(config), fixed = TRUE) ## replace space with dot from the names in the config file
    tryCatch({
        input.file <- config$input.file
        add.meas.name <- config$additional
        data.tp <- config$data.type
        rawdata <- read.csv(input.file)
    }, error = function(e) {
        stop("Error: ", e$message, "\n")
    })
    
    ## Check if 'Name' column is present in rawdata
    if (!"Name" %in% colnames(rawdata)) {
        stop("INPUT FILE ERROR: the input file does not contain a column named 'Name'\n")
    }
    
    ## Read the sample ID
    names <- rawdata$Name
    
    ## Read the wavelength directly from the input file header
    wlength <- numeric(0)
    u.wlength <- numeric(0)
    for (col in colnames(rawdata)) {
        if (grepl("\\d+$", col)) {
            wlength <- c( wlength, as.numeric(stringr::str_extract(col, "\\d+$")) )
            u.wlength <- c(u.wlength, config$wavelength.error)
        }
    }
    
    ## Check if a mass apportionment is required
    mass.appo.requested <- "EC" %in% colnames(rawdata) && "OC" %in% colnames(rawdata)

    ## to be updated
    ## Create the list of sample objects. This will be the output
    data <- lapply(1:nrow(rawdata), function(x) {
        as.data.frame(matrix(ncol = length(nn), nrow = length(wlength), dimnames = list(NULL, nn)))} ## Every list element is a sample. Create a data.frame for each list element with the proper column names
    )

    ## update the values of the output
    data <- lapply(data, function(x) {
        transform(x,
                  "wavelength" = wlength,
                  "u.wavelength" = u.wlength,
                  "data.type" = data.tp)
    })
    names(data) <- rawdata$Name

    data <- lapply(seq_along(data), function(x) {
        transform(data[[x]],
                  "name" = names(data)[x],
                  "abs" = as.numeric(rawdata[x, match(paste0("X", wlength), names(rawdata))]),
                  "u.abs" = default.abs.unc * as.numeric(rawdata[x, match(paste0("X", wlength), names(rawdata))])
                  )
    })
    
    ## Fill with carbon concentrations if needed
    if (mass.appo.requested == TRUE) {
        data <- lapply(seq_along(data), function(x) {
            transform(data[[x]],
                      "ec" = rawdata[x, "EC"],
                      "oc" = rawdata[x, "OC"]
                      )
        })
    }
    
    ## add additional meas name, if needed
    if (length(add.meas.name) > 0) {
        data <- lapply(seq_along(data), function(x) {
            transform(data[[x]],
                      "add.meas" = rawdata[x, grep(add.meas.name, names(rawdata), fixed = TRUE)]
                      )
        })
    }

    names(data) <- rawdata$Name
    
    return(data)
}
