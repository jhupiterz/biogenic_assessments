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

# ---------------------------------------------------------------------------- #
# -------------------------------- IMPORTS ----------------------------------- #

setwd(getwd())
source("GSA_Catling_public.R")
library(ggplot2)
library(plotly)

max_test <- 15
Catling_sobol <- Sobol_index(max_test)
#print(Catling_sobol)

# ---------------------------------------------------------------------------- #
# --------------------------------- PLOTS ------------------------------------ #

# ------------------------ FIRST ORDER SOBOL INDICES ------------------------- #
par(mfrow=c(1,1))

plot_prior_1st          <- Catling_sobol[,1]
plot_likelihood_1st     <- plot_prior_1st + Catling_sobol[,2]
plot_false_positive_1st <- plot_likelihood_1st + Catling_sobol[,3]

#dev.new(width=5, height=4)
plot(plot_prior_1st, xlim= c(1,max_test), ylim= c(0,1), type = "l", xlab = "Number of biosignatures", ylab = "1st order Sobol Index", xaxs="i", yaxs="i")
lines(plot_likelihood_1st, col = "blue")
lines(plot_false_positive_1st, col = "red")
axis(1, at=1:max_test, lab=c(1:max_test))

xx_prior_1st = c(1, c(2:max_test), max_test, 1)
yy_prior_1st = c(plot_prior_1st[1], plot_prior_1st[2:max_test], 0, 0)
polygon(xx_prior_1st, yy_prior_1st, col= "#fdb863")

xx_likelihood_1st = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
yy_likelihood_1st = c(plot_prior_1st[1], plot_likelihood_1st[1], plot_likelihood_1st[2:max_test], plot_prior_1st[max_test], plot_prior_1st[(max_test-1):1])
polygon(xx_likelihood_1st, yy_likelihood_1st, col= "#b2abd2")

xx_false.positive_1st = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
yy_false.positive_1st = c(plot_likelihood_1st[1], plot_false_positive_1st[1], plot_false_positive_1st[2:max_test], plot_likelihood_1st[max_test], plot_likelihood_1st[(max_test-1):1])
polygon(xx_false.positive_1st, yy_false.positive_1st, col= "#abd9e9")

xx_interactions_1st = c(1, 1, max_test, max_test, (max_test-1):1)
yy_interactions_1st = c(plot_false_positive_1st[1], 1, 1, plot_false_positive_1st[max_test], plot_false_positive_1st[(max_test-1):1])
polygon(xx_interactions_1st, yy_interactions_1st, density = 60, angle = 350)

i <- 2
text(2, 0.1, "Prior", cex = 1)
text((1.5+i), 0.35, "Likelihood", cex = 1)
text((1.5+i), 0.7, "False positive", cex = 1)
text((2.5+i), 0.92, "Interactions", cex = 1)

# ----------------------- SECOND ORDER SOBOL INDICES ------------------------- #
# ------- Uncomment this section to plot second order Sobol indices ---------- #

#plot_prior_2nd          <- Catling_sobol[,4]
#plot_likelihood_2nd     <- plot_prior_2nd + Catling_sobol[,5]
#plot_false_positive_2nd <- plot_likelihood_2nd + Catling_sobol[,6]

#dev.new(width=5, height=4)
#plot(plot_prior_2nd, xlim= c(1,max_test), ylim= c(0,1), type = "l", xlab = "Number of biosignatures", ylab = "2nd order Sobol Index", xaxs="i", yaxs="i")
#lines(plot_likelihood_2nd, col = "blue")
#lines(plot_false_positive_2nd, col = "red")
#axis(1, at=1:max_test, lab=c(1:max_test))

#X1X2
#xx_prior_2nd = c(1, c(2:max_test), max_test, 1)
#yy_prior_2nd = c(plot_prior_2nd[1], plot_prior_2nd[2:max_test], 0, 0)
#polygon(xx_prior_2nd, yy_prior_2nd, col= "#1b9e77")

#X1X3
#xx_likelihood_2nd = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
#yy_likelihood_2nd = c(plot_prior_2nd[1], plot_likelihood_2nd[1], plot_likelihood_2nd[2:max_test], plot_prior_2nd[max_test], plot_prior_2nd[(max_test-1):1])
#polygon(xx_likelihood_2nd, yy_likelihood_2nd, col= "#d95f02")

#X2X3
#xx_false.positive_2nd = c(1, 1, c(2:max_test), max_test, (max_test-1):1)
#yy_false.positive_2nd = c(plot_likelihood_2nd[1], plot_false_positive_2nd[1], plot_false_positive_2nd[2:max_test], plot_likelihood_2nd[max_test], plot_likelihood_2nd[(max_test-1):1])
#polygon(xx_false.positive_2nd, yy_false.positive_2nd, col= "#7570b3")

#xx_interactions_2nd = c(1, 1, max_test, max_test, (max_test-1):1)
#yy_interactions_2nd = c(plot_false_positive_2nd[1], 1, 1, plot_false_positive_2nd[max_test], plot_false_positive_2nd[(max_test-1):1])
#polygon(xx_interactions_2nd, yy_interactions_2nd, density = 60, angle = 350)
