## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Main wrapper for the MWAA Toolkit  ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## specifiy the path to the configuration file
config.file <- "../../../R_debug/files/p_paper_config.json"

## load the script with the code for the apportionment
source("flow.R")

## execute the code for the apportionment
mwaa.mt(config.file)

