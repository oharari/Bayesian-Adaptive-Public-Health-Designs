library(ggplot2)
library(gridExtra)
library(foreach)
library(tcltk)
library(rlecuyer)
library(tictoc)
library(lemon)
library(scales)
library(dplyr)
library(knitr)
library(kableExtra)
library(magick)
library(gtools)
library(data.table)
library(doParallel)
library(doSNOW)
library(snow)
library(foreach)
library(tidyverse)
library(reshape2)
library(pracma)
library(paletteer)
library(grDevices)


#Setting working directory to current one
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Sourcing functions
source('Vaccine_Trial_Functions.R')
source('Negative_Binomial_Functions.R')


#*****************************************************************************
#*****************************************************************************
# 1. Bayesian adaptive vaccine trials ----------------------------------------
#*****************************************************************************
#*****************************************************************************


#*****************************************************************************
# 1.0 Allocation ratio -------------------------------------------------------
#*****************************************************************************
D_opt = D_optimal_alloc_ratio(p_ctrl = .05, n_trts = 3, RRR = .375)
D_opt$rho
D_opt$plot

ratio = D_opt$rho
#*****************************************************************************
# 1.1 Prior selection --------------------------------------------------------
#*****************************************************************************

# 1.1.1 Inputs ---------------------------------------------------------------
q = 1/(1 + ratio) # the (normalized) allocation probability to vaccine
alpha_trt = c(1, 2, 5, 15) # Vaccine Dirichlet hyperparameters to explore


# 1.1.2 Summary table --------------------------------------------------------

# Selecting control Dirichlet hyperparameters to match prior probability of
# vaccine efficacy of 50%
alpha_ctrl = sapply(alpha_trt, alpha_ctrl_func, q = q)

alpha = cbind(alpha_trt, alpha_ctrl)

# Summarizing prior VE properties
priors = do.call(rbind,
                 apply(alpha, 1, VE_prior_summary, q = q)
)

priors %>% 
  kable(align = 'c', escape = F) %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped"))


# 1.1.3 Prior density plots --------------------------------------------------
LL = -1.5 # lower x-axis limit
UL = 1 # upper x-axis limit
delta = 5e-3 # grid resolution

VE_prior_1D_density_plot(LL, UL, delta, alpha, trt_prob = q)


#*****************************************************************************
# 1.2 Random adaptive vaccine trial simulation -------------------------------
#*****************************************************************************

# 1.2.1 Inputs ---------------------------------------------------------------
VE = c(.5, .375, .25) # vaccine effects
d = length(VE) # number of vaccines
trt_prob = 1/(d + ratio) # vaccine allocation probability
ctrl_prob = ratio*trt_prob # control allocation probability
alloc_probs = c(ctrl_prob, rep(trt_prob, d))
n_analyses = 4 # number of analyses (including final)
effic_thresh = rep(.99015, n_analyses) # efficacy threshold vector
futil_thresh = c(.25, .5, .75, 1) # futility threshold vector

# number of events per remaining arm to trigger an interim analysis 
N_events_per_arm_per_analysis = 10
# Dirichlet prior hyperparameters: (vaccine, control)
hyperpar = c(1, 1.078011)


# 1.2.2 Random trial simulation ----------------------------------------------
LL = -1 # lower x-axis limit
UL = 1 # upper x-axis limit
delta = 5e-3 # grid resolution

trial = vaccine_trial_plots_and_table(N_events_per_arm_per_analysis, VE, 
                                      n_analyses, alloc_probs, hyperpar, 
                                      futil_thresh, effic_thresh, LL, UL, 
                                      delta, seed = 39)


# 1.2.3 Posterior probability of efficacy plots ------------------------------
trial$p_post_effic


# 1.2.4 Summary table and posterior inference --------------------------------
trial$sum_tab %>%
  kable(align = 'c') %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped")) 


# 1.2.5 VE posterior density plots -------------------------------------------
trial$VE_post_plot


#*****************************************************************************
# 1.3 Calculating design operating characteristics ---------------------------
#*****************************************************************************
# change to TRUE if you wish to rerun the full MC simulation, otherwise an old
# one will be loaded to save time
Simulate_now = FALSE  

# 1.3.1 Full Bayesian adaptive Vaccine trial design simulation ---------------
# N_simulations = number of random trials to be simulated

if(!Simulate_now){
  load("Full_Vaccine_Trial_Design_Simulation.RDATA")
} else{
  full_sim = full_vaccine_trial_simulation(N_simulations = 1e5,
                                           N_events_per_arm_per_analysis, 
                                           VE, n_analyses, 
                                           alloc_probs, hyperpar, 
                                           futil_thresh, effic_thresh)
}

# 1.3.2 Operating characteristics summary table ------------------------------
full_sim$operating_chars %>%
  kable(align = 'c', escape = F) %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped"))







#*****************************************************************************
#*****************************************************************************
# 2. Bayesian adaptive cluster randomized trials -----------------------------
#*****************************************************************************
#*****************************************************************************

#*****************************************************************************
# 2.1 Prior selection --------------------------------------------------------
#*****************************************************************************

# 2.1.1 Inputs ---------------------------------------------------------------
# vector of intracluster correlation coefficient parameters to explore
ICC_vec = c(.15, .2) 

# vector of mean event rate parameters prior means to explore
mu_mean_vec = c(.1) 

# vector of beta prior first shape parameters to explore
beta_a_vec = c(2, 4, 8, 16)

#vector of mean cluster sizes to explore
n_vec = c(50)


# 2.1.2 Prior attributes on the mean event rate scale ------------------------
# summarizing the resultant beta priors and displaying their properties
priors = mu_prior_summary(beta_a_vec, n_vec, mu_mean_vec, ICC_vec)

priors %>% 
  kbl(align = 'c', format = 'html') %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped", 'bordered')) %>% 
  collapse_rows(columns = 1:3)


# 2.1.3 Prior mean event rate density plots  ---------------------------------
mu_mean = .1 # mean event rate expected value
beta_a_vec = c(2, 4, 8, 16) # first beta prior shape hyperparameter
ICC = .15 # intracluster correlation
mean_clust_size = 50 # mean cluster size
LL = 5e-4 # x-axis lower limit
UL = 1 # x-axis upper limit
delta = 5e-3 # grid resolution

mu_prior_plot(mu_mean, beta_a_vec, ICC, mean_clust_size, LL, UL, delta)


#*****************************************************************************
# 2.2 Random adaptive cluster trial simulation -------------------------------
#*****************************************************************************

# 2.2.1 Inputs ---------------------------------------------------------------
mean_clust_size = 50 # mean cluster size
ICC = .15 # intracluster correlation coefficient
mu_vec = c(.1, .075, .05) # mean event rate means
Trt_nms = c('A', 'B', 'C') # treatment arm names
N_clusts_max = 550 # maximum number of clusters
N_analyses = 3 # maximum number of analyses (including final)

# superiority threshold on the posterior probability of superiority
sup_thresh = .988  

# futility threshold on the posterior probability of superiority
fut_thresh = .012

#first beta prior shape parameter
beta_a = 2


# 2.2.2 Matching second prior shape parameter to mean even rate --------------
#second beta prior shape parameter
(beta_b = beta_prior_b_hyperpar(beta_a, mu_mean = .1, ICC, mean_clust_size))


# 2.2.3 Matching random effect shape parameter to ICC ------------------------
# Gamma random effect shape parameters
(true_phi_vec = ICC_to_shape(mu_vec, ICC))


# 2.2.4 Simulating a random Bayesian adaptive cluster trial ------------------
# N_MC = the Monte Carlo sample size of the beta posterior samples for the
# evaluation of posterior probability of superiority
set.seed(2)
simulation = single_cluster_trial_simulation(mu_vec, N_clusts_max,
                                             N_analyses, 
                                             mean_clust_size,
                                             ICC, beta_a, beta_b, 
                                             N_MC = 1e6, sup_thresh,
                                             fut_thresh)



# 2.2.5 Plotting the posterior probability of superiority over time ----------
cluster_prob_sup_plot(simulation, fut_thresh, sup_thresh)


# 2.2.6 Empirical Bayes estimation of the random effect shape parameters -----
phi_vec = cluster_empirical_bayes_shape(simulation$data, beta_a, beta_b)


# 2.2.7 Plotting the log-marginal likelihood curves --------------------------
# Plotting the log-marginal likelihood curves along with the estimates and
# the true gamma random effect shape values

LL = .2 # x-axis lower limit
UL = 1 # x-axis upper limit
delta = 1e-3 # grid resolution

cluster_marginal_likelihood_plot_all(simulation$data, 
                                     beta_a, beta_b, true_phi_vec,
                                     LL, UL, delta)



# 2.2.8 Posterior inference and summary table --------------------------------
# Posterior inference
single_trial_results = single_cluster_simulation_results(simulation, 
                                                         beta_a, beta_b)

single_trial_results %>%
  kable(align = c(rep('c', 5), 'r'), escape = FALSE) %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped")) %>%
  row_spec(0, align = 'c')


# 2.2.9 Posterior mean event rate density plots ------------------------------
LL = .035 # x-axis lower limit
UL = .14 # x-axis upper limit
delata = 5e-3 # grid resolution

cluster_post_density_plot(simulation, phi_vec, beta_a, beta_b, LL, UL, delta)


#*****************************************************************************
# 2.3 Calculating design operating characteristics ---------------------------
#*****************************************************************************

# 2.3.1 Full adaptive cluster trial design simulation ------------------------

# each row in this matrix is a vector of mean event rates: one "alternative" 
# scenario and two "null" scenarios
mu_matrix = matrix(c(.1, .075, .05,
                     .1, .1, .1,
                     .1, .05, .05), 
                   ncol = 3, byrow = TRUE)

# N_MC = the Monte Carlo sample size of the beta posterior samples for the
# evaluation of posterior probability of superiority

# N_simulations = number of random trials to be simulated

# change to TRUE if you wish to rerun the full MC simulation, otherwise an old
# one will be loaded to save time
Simulate_now = FALSE  

if(!Simulate_now){
  load("Full_Cluster_Trial_Design_Simulation.RDATA")
} else{
  oper_charac = full_cluster_trial_simulation(mu_vec, N_clusts_max,
                                              N_analyses, mean_clust_size,
                                              ICC, beta_a, beta_b, 
                                              N_MC = 3e3, sup_thresh, 
                                              fut_thresh,
                                              N_simulations = 1e5)
}

# 2.3.2 Operating characteristics summary table ------------------------------
oper_charac %>% 
  kable(align = 'c', escape = F) %>%
  kable_styling(full_width = F, 
                bootstrap_options = c("striped"))

