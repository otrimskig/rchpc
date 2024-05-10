#install.packages("foreach")
library(foreach)

if (!exists("n.cores")) {

"initilizing cores..."
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
"parallel cores initialized."

}

your_vector<-c(1:10)


foreach(i=1:length(your_vector)) %dopar% {
  
#if you use a library you have to load it within the loop
#(or specify command with ::)

print(i)
  
}
