#########################################################################################################
#                                     GLOBAL SENSITIVITY ANALYSIS                                       #
#                                         CATLING ET AL. 2018                                           #
#-------------------------------------------------------------------------------------------------------#
#                                         Author: Julie Hartz                                           #
#                                      Last updated: 01.09.2024                                         #
#-------------------------------------------------------------------------------------------------------#
#                                         Full citation APA:                                            #
#                                 Hartz, J., & George, S. C. (2022)                                     #
#.        Quantitative Framework for Astrobiology Strategies and in situ Biogenic Assessments.          #
#                        Frontiers in Astronomy and Space Sciences, 9, 769607                           #
#########################################################################################################

#------------------------------------------------IMPORTS------------------------------------------------#
library(sensitivity)
library(boot)

#----------------------------------------------PARAMETERS----------------------------------------------#
# number of parameter combinations
n <- 1000000

# Build the matrices of variables randomly distributed between 0 and 1, can change min and max if necessary
X3 <- data.frame(matrix(runif(1 * n, min = 0, max = 1), nrow = n))
X4 <- data.frame(matrix(runif(1 * n, min = 0, max = 1), nrow = n))
X5 <- data.frame(matrix(runif(2 * n), nrow = n))
X6 <- data.frame(matrix(runif(2 * n), nrow = n))
X7 = cbind(X3, X5)
X8 = cbind(X4, X6)

#--------------------------------------------BAYESIAN MODEL---------------------------------------------#
Catling_multipleTest <- function(M, number_test) {
  # Bayesian inference based on Catling et al. 2018
  # Returns a vector containing the i-th posterior probability calculated for each parameter combination
  # M: matrix of ranging varibales
  # number_test: number of inferences (i.e. biosignatures)
  posterior <- data.frame(matrix(nrow=nrow(M), ncol=number_test+1))
  posterior[,1] <- M[,1]
    for (k in 2:(number_test+1)) {
      posterior[,k] <- (M[,2]*posterior[,k-1])/((M[,2]*posterior[,k-1])+(M[,3]*(1-posterior[,k-1])))
    }
  return(posterior[,number_test+1])
}

#-----------------------------------------SENSITIVITY ANALYSIS------------------------------------------#
Sobol_index <- function(max_test){
  # Performs Sobol sensitivity analysis on specified model
  # max_test: total number of inferences (i.e. biosignatures)
  # Returns a matrix containing the first order indices for each parameter and for max_test number of biosignatures
  Sobol_index <- c()
  for (i in 1:max_test) {
    y <- sobolSalt(model = Catling_multipleTest, X1 = X7, X2 = X8, number_test = i, scheme = 'B', nboot = 100)
    Sobol_index <- rbind(Sobol_index,cbind(t(y$S[,1]), t(y$S2[,1])))
  }
  colnames(Sobol_index) <- c('Prior_1st', 'Likelihood_1st', 'False_positive_1st', 'Prior_2nd', 'Likelihood_2nd', 'False_positive_2nd')
  return(Sobol_index)
}

# #----------------------------------------COMPUTATION & RESULTS------------------------------------------#
# max_test <- 15 # total number of inferences (i.e. biosignatures)
# Catling_sobol <- Sobol_index(max_test)
# print(Catling_sobol)
# 
# 
# #-----------------------------------------------PLOTTING------------------------------------------------#
# # First order Sobol indices
# par(mfrow=c(1,1))
# 
# plot_prior_1st          <- Catling_sobol[,1]
# plot_likelihood_1st     <- plot_prior_1st + Catling_sobol[,2]
# plot_false_positive_1st <- plot_likelihood_1st + Catling_sobol[,3]
#   
# #dev.new(width=5, height=4)
# plot(plot_prior_1st, xlim= c(1,max_test), ylim= c(0,1), type = "l", xlab = "Number of biosignatures", ylab = "1st order Sobol Index", xaxs="i", yaxs="i")
# lines(plot_likelihood_1st, col = "blue")
# lines(plot_false_positive_1st, col = "red")
# axis(1, at=1:max_test, lab=c(1:max_test))
# 
# xx_prior_1st = c(1, c(2:max_test), max_test, 1)
# yy_prior_1st = c(plot_prior_1st[1], plot_prior_1st[2:max_test], 0, 0)
# polygon(xx_prior_1st, yy_prior_1st, col= "#fdb863")
# 
# xx_likelihood_1st = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
# yy_likelihood_1st = c(plot_prior_1st[1], plot_likelihood_1st[1], plot_likelihood_1st[2:max_test], plot_prior_1st[max_test], plot_prior_1st[(max_test-1):1])
# polygon(xx_likelihood_1st, yy_likelihood_1st, col= "#b2abd2")
# 
# xx_false.positive_1st = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
# yy_false.positive_1st = c(plot_likelihood_1st[1], plot_false_positive_1st[1], plot_false_positive_1st[2:max_test], plot_likelihood_1st[max_test], plot_likelihood_1st[(max_test-1):1])
# polygon(xx_false.positive_1st, yy_false.positive_1st, col= "#abd9e9")
# 
# xx_interactions_1st = c(1, 1, max_test, max_test, (max_test-1):1)
# yy_interactions_1st = c(plot_false_positive_1st[1], 1, 1, plot_false_positive_1st[max_test], plot_false_positive_1st[(max_test-1):1])
# polygon(xx_interactions_1st, yy_interactions_1st, density = 60, angle = 350)
# 
# i = 2
# text(1.95, 0.06, "Prior", cex = 1)
# text((1.5+i), 0.37, "Likelihood", cex = 1)
# text((1.5+i), 0.75, "False positive", cex = 1)
# text((2.5+i), 0.92, "Interactions", cex = 1)
# 
# #------------------------------------------OPTIONAL PLOTTING-------------------------------------------#
# # Second order Sobol indices
# plot_prior_2nd          <- Catling_sobol[,4]
# plot_likelihood_2nd     <- plot_prior_2nd + Catling_sobol[,5]
# plot_false_positive_2nd <- plot_likelihood_2nd + Catling_sobol[,6]
# 
# plot(plot_prior_2nd, xlim= c(1,max_test), ylim= c(0,1), type = "l", xlab = "Number of biosignatures", ylab = "2nd order Sobol Index", xaxs="i", yaxs="i")
# lines(plot_likelihood_2nd, col = "blue")
# lines(plot_false_positive_2nd, col = "red")
# axis(1, at=1:max_test, lab=c(1:max_test))
# 
# #X1X2
# xx_prior_2nd = c(1, c(2:max_test), max_test, 1)
# yy_prior_2nd = c(plot_prior_2nd[1], plot_prior_2nd[2:max_test], 0, 0)
# polygon(xx_prior_2nd, yy_prior_2nd, col= "#1b9e77")
# 
# #X1X3
# xx_likelihood_2nd = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
# yy_likelihood_2nd = c(plot_prior_2nd[1], plot_likelihood_2nd[1], plot_likelihood_2nd[2:max_test], plot_prior_2nd[max_test], plot_prior_2nd[(max_test-1):1])
# polygon(xx_likelihood_2nd, yy_likelihood_2nd, col= "#d95f02")
# 
# #X2X3
# xx_false.positive_2nd = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
# yy_false.positive_2nd = c(plot_likelihood_2nd[1], plot_false_positive_2nd[1], plot_false_positive_2nd[2:max_test], plot_likelihood_2nd[max_test], plot_likelihood_2nd[(max_test-1):1])
# polygon(xx_false.positive_2nd, yy_false.positive_2nd, col= "#7570b3")





