library(plumber)
library(servr)

#library(callr)
# 
# back <- r_bg(function(){
#   api_obj <- plumber::plumb('plumber.R')
#   api_obj$run(port = 8080, host = '0.0.0.0')
# })
# 
# front <- r_bg(function(){
#   servr::httw('', port = '80', daemon = FALSE)
# })

#setwd('~/Documents/projects/nag/sus-app/')

servr::httd('.', port = '80', host = '0.0.0.0', daemon = TRUE)

api_obj <- plumber::plumb('R/plumber.R')
api_obj$run(port = 8080, host = '0.0.0.0')
