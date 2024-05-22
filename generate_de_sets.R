
#parallel cores if on cluster
##########################
#library(foreach)
# if (!exists("n.cores")) {
#   
#   "initilizing cores..."
#   n.cores <- parallel::detectCores() - 1
#   my.cluster <- parallel::makeCluster(
#     n.cores, 
#     type = "PSOCK"
#   )
#   doParallel::registerDoParallel(cl = my.cluster)
#   
#   #check if it is registered (optional)
#   foreach::getDoParRegistered()
#   
#   "parallel cores initialized."
#   
# }




