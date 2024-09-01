#########################################################################################################
#################################### GLOBAL SENSITIVITY ANALYSIS ########################################
######################################## CATLING ET AL. 2018 ############################################
###################################### VARIANCE-BASED ANALYSIS ##########################################
####################################### (SALTELLI ET AL. 2010) ##########################################
#########################################################################################################

library(sensitivity)
library(boot)

#-------------------------------------------------------------------------------------------------------#
# number of parameters combinations
n <- 1000000
# Matrix containing the n parameters combinations, dim = n x #parameters in the model (here 3)
X1 <- data.frame(matrix(runif(3 * n), nrow = n))
X2 <- data.frame(matrix(runif(3 * n), nrow = n))

#-------------------------------------------------------------------------------------------------------#
## Bayesian inference based on Catling et al. 2018 for 'number_bio' biosignatures 
## Returns a vector containing the n-th posterior probability calculated for each parameter combination
Catling_multipleTest <- function(M, number_test) {
  posterior <- data.frame(matrix(nrow=nrow(M), ncol=number_test+1))
  posterior[,1] <- M[,1]
    for (k in 2:(number_test+1)) {
      posterior[,k] <- (M[,2]*posterior[,k-1])/((M[,2]*posterior[,k-1])+(M[,3]*(1-posterior[,k-1])))
    }
  return(posterior[,number_test+1])
}

# Returns a matrix containing the first order indices for each parameter and for max_test number of biosignatures
Sobol_index <- function(max_test=15){
  Sobol_index <- c()
  for (i in 1:max_test) {
    y <- sobolSalt(model = Catling_multipleTest, X1 = X1, X2 = X2, number_test = i, scheme = 'B', nboot = 100)
    Sobol_index <- rbind(Sobol_index,cbind(t(y$S[,1]), t(y$S2[,1])))
  }
  colnames(Sobol_index) <- c('Prior_1st', 'Likelihood_1st', 'False_positive_1st', 'Prior_2nd', 'Likelihood_2nd', 'False_positive_2nd')
  return(Sobol_index)
}







