#*******************************************************************************
#*******************************************************************************
#* The functions used for the manuscript "Adaptive Trial Designs for Global 
#* Health: Vaccine and Cluster Randomized Trials go Bayesian" (2023).
#* 
#* All code was written by Ofir Harari
#*******************************************************************************
#*******************************************************************************




#*******************************************************************************
#* 1. Generating random data for the Gamma-Poisson mixture (Negative Binomial) 
#* model.
#* 
#* Inputs: 
#*   Trt_nms - a vector of treatment names
#*   N_clust_vec - number of clusters in each treatment arm
#*   mean_clust_size - mean cluster size (for the Poisson distribution)
#*   true_phi_vec - a vector of true shape parameters 
#*   mu_vec - a vector of mean event rates
#*   
#* Output: a data frame consisting of treatment names, cluster sizes and 
#*         numbers of cases.
#*******************************************************************************
sample_NB_data = function(Trt_nms, N_clust_vec, mean_clust_size,
                          true_phi_vec, mu_vec){
  sample_data_single_arm = function(N_clusts, 
                                    mean_clust_size,
                                    true_phi, mu){
    n = rpois(N_clusts, mean_clust_size)
    r = rgamma(N_clusts, shape = true_phi, scale = mu/true_phi)
    y = rbinom(N_clusts, n, pmin(r,1))
    
    return(as.data.frame(cbind(n, y)))
  }
  
  Trt = rep(Trt_nms, N_clust_vec)
  data = do.call(rbind,
                 lapply(as.list(1:length(mu_vec)),
                        function(i){
                          sample_data_single_arm(N_clust_vec[i], 
                                                 mean_clust_size,
                                                 true_phi_vec[i], 
                                                 mu_vec[i])
                        }
                 )
  )
  
  out = cbind(Trt, data)
  
  return(out)
}




#*******************************************************************************
#* 2. Drawing beta random samples from the beta posteriors of the mean event 
#* rate parameters.
#* 
#* Inputs: 
#*   data - the output of 'sample_NB_data'
#*   phi_vec - a vector of (empirical Bayes estimated) shape parameters 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   N_MC - size of Monte Carlo samples to be drawn
#*   
#* Output: a matrix whose columns are the posterior samples corresponding to the 
#*         different treatment arms.
#*******************************************************************************
cluster_posterior_samps = function(data, phi_vec, beta_a, 
                                   beta_b, N_MC){
  arms = unique(data$Trt)
  
  post_samp_single_arm = function(y, n, phi, beta_a, beta_b, N_MC){
    K = length(y)
    alpha = beta_a + sum(y)
    beta = beta_b + K*phi
    p_post = rbeta(N_MC, alpha, beta)
    mu_post = phi*p_post/(1 - p_post)/(mean(n))
    
    return(mu_post)
  }
  
  posts = sapply(1:length(arms),
                 function(i){
                   data_temp = data %>% filter(Trt == arms[i])
                   y = data_temp$y
                   n = data_temp$n
                   
                   return(post_samp_single_arm(y, n, phi_vec[i],
                                               beta_a, beta_b, 
                                               N_MC))
                 }
  )
  
  return(posts)
}




#*******************************************************************************
#* 3. Calculating the posterior probability of superiority of each arm.
#* 
#* Inputs: 
#*   post_samps - the output of 'cluster_posterior_samps'
#*   
#* Output: a vector of posterior probabilities of superiority corresponding to 
#*         the different treatment arms.
#*******************************************************************************
cluster_post_superiority = function(post_samps){
  arms = 1:ncol(post_samps)
  
  p = table(apply(post_samps, 1, which.min))/nrow(post_samps)
  miss = arms[!(arms %in% names(p))]
  if(length(miss) > 0){
    p = c(p, rep(0, length(miss)))
    names(p)[names(p) == ''] = miss
    p = p[order(names(p))]
  }
  
  return(p)
}




#*******************************************************************************
#* 4. Simulating a single random adaptive CRT.
#* 
#* Inputs: 
#*   mu_vec - a vector of mean event rates
#*   N_clusts_max - maximum overall number of clusters
#*   N_analyses - maximum number of analyses (interim + final)
#*   mean_clust_size - mean cluster size (for the Poisson distribution)
#*   ICC - the intracluster correlation coefficient 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   N_MC - size of Monte Carlo samples to be drawn
#*   sup_thresh - superiority threshold (on the posterior probability of 
#*                superiority)
#*   fut_thresh - futility threshold (on the posterior probability of 
#*                superiority)
#*   
#* Output: a list consisting of the following slots - 
#*   out - a table containing the numbers of clusters, patients and cases for 
#*         each arm in each analysis as well as the end outcome (Early 
#*         inferiority/Superiority/none)
#*   data - the cluster-level data (arm, patients and cases)
#*   prob_sup - a matrix whose columns contain the posterior probability of 
#*              superiority by analysis for each arm 
#*******************************************************************************
single_cluster_trial_simulation = function(mu_vec, N_clusts_max,
                                           N_analyses, 
                                           mean_clust_size,
                                           ICC, beta_a, beta_b, 
                                           N_MC, sup_thresh,
                                           fut_thresh){
  true_phi_vec = ICC_to_shape(mu_vec, ICC)
  
  num_treatments = num_treatments_old = length(mu_vec)
  Trts = 1:num_treatments
  N_Clusts_remaining = N_clusts_max
  mu_vec_new = mu_vec
  Trts_new = Trts
  true_phi_vec_new = true_phi_vec
  prob_superiority = prob_superiority_old = rep(1/num_treatments, 
                                                num_treatments)
  Outcome = rep('', num_treatments)
  Num_analyses = rep(0, length(mu_vec))
  
  i = 1
  Winner = NA
  Futile = integer()
  data_full = data.frame(Trt = integer(),
                         n = integer(),
                         y = integer(),
                         Analysis = integer())
  
  while(N_Clusts_remaining > 0 & i <= N_analyses & length(Trts_new) > 1){
    if(length(Futile) > 0){
      Num_analyses[-Futile] = Num_analyses[-Futile] + 1
    } else{
      Num_analyses = Num_analyses + 1
    }
    
    N_analyses_remaining = N_analyses - i + 1
    
    N_clusts = ceiling(N_Clusts_remaining/N_analyses_remaining/num_treatments)
    N_clust_new = ifelse(N_Clusts_remaining <= N_clusts*num_treatments,
                         floor(N_Clusts_remaining/num_treatments),
                         N_clusts)
    N_clust_vec = rep(N_clust_new, num_treatments)
    
    data_new = sample_NB_data(Trts_new, N_clust_vec, mean_clust_size,
                              true_phi_vec_new, mu_vec_new) %>% 
      mutate(Analysis = i)
    
    data_full = rbind(data_full, data_new)
    
    phi_vec = cluster_empirical_bayes_shape(data_full, beta_a, beta_b)
    
    post_samps = cluster_posterior_samps(data_full %>% 
                                           filter(Trt %in% Trts_new), 
                                         phi_vec, beta_a, beta_b, N_MC)
    prob_superiority[Trts_new] = cluster_post_superiority(post_samps)
    prob_superiority_old = rbind(prob_superiority_old, prob_superiority)
    
    N_Clusts_remaining = N_Clusts_remaining - sum(N_clust_vec)
    i = i + 1
    
    sup_candidate = which.max(prob_superiority)
    if(prob_superiority[sup_candidate] >= sup_thresh){
      Winner = sup_candidate
      Outcome[sup_candidate] = 'Superiority'
      break
    } else{ 
      futiles = Trts[which(prob_superiority <= fut_thresh)]
      if(length(futiles) > 0){
        Futile = unique(c(Futile, futiles))
        mu_vec_new = mu_vec[-Futile]
        Trts_new = Trts[-Futile]
        true_phi_vec_new = true_phi_vec[-Futile]
        phi_vec = phi_vec[-Futile]
        num_treatments = num_treatments_old - length(Futile)
        Outcome[Futile] = 'Early inferiority'
      }
    }
  } 
  
  prob_superiority_old = sapply(1:length(Num_analyses), 
                                function(i){
                                  l = Num_analyses[i] + 1
                                  m = nrow(prob_superiority_old)
                                  x = prob_superiority_old[,i]
                                  if(l < m){
                                    x[(l+1):m] = NA
                                  }
                                  
                                  return(x)
                                })
  
  out = data_full %>% 
    group_by(Trt, Analysis) %>% 
    summarise(Clusters = n(),
              Patients = sum(n),
              Cases = sum(y)) %>% 
    mutate(Clusters = cumsum(Clusters),
           Patients = cumsum(Patients),
           Cases = cumsum(Cases)) %>% 
    mutate(`Clusters, Cases/Patients` = paste0(Clusters, ', ', 
                                               Cases, 
                                               '/', Patients)) %>% 
    dplyr::select(-c(Clusters, Patients, Cases)) %>% 
    reshape2::dcast(Trt ~ Analysis, 
                    value.var = 'Clusters, Cases/Patients') %>% 
    as.data.frame() %>% 
    mutate(Outcome = Outcome)
  
  out[is.na(out)] = ''
  
  names(out)[-c(1, ncol(out))] = paste0('Analysis ', 
                                        1:(ncol(out) - 2),
                                        '<br>Clusters, Cases/Patients')
  
  return(list(out = out, 
              prob_sup = unname(prob_superiority_old),
              data = data_full))
}




#*******************************************************************************
#* 5. A lighter version of 'single_cluster_trial_simulation' (used for 
#* simulation rather than display).
#* 
#* Inputs: 
#*   mu_vec - a vector of mean event rates
#*   N_clusts_max - maximum overall number of clusters
#*   N_analyses - maximum number of analyses (interim + final)
#*   mean_clust_size - mean cluster size (for the Poisson distribution)
#*   ICC - the intracluster correlation coefficient 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   N_MC - size of Monte Carlo samples to be drawn
#*   sup_thresh - superiority threshold (on the posterior probability of 
#*                superiority)
#*   fut_thresh - futility threshold (on the posterior probability of 
#*                superiority)
#*   
#* Output: a data frame containing the numbers of clusters, patients, cases and 
#*         the end outcome (Early inferiority/Superiority/none) for each arm at
#*         the end of the trial
#*******************************************************************************
single_cluster_trial_simulation_fast = function(mu_vec, N_clusts_max,
                                                N_analyses, 
                                                mean_clust_size,
                                                ICC, beta_a, beta_b, 
                                                N_MC, sup_thresh,
                                                fut_thresh){
  oper_char_single_cluster_trial = function(simulation, mu_vec){
    best = which.min(mu_vec)
    
    prop_clust = simulation$Clusters/sum(simulation$Clusters)
    Power = max(simulation$Outcome == 'Superiority' & 
                  simulation$Trt == best)
    Type_I_err = max(simulation$Outcome == 'Superiority' & 
                       simulation$Trt != best)
    
    Clusters = sum(simulation$Clusters)
    
    out = data.frame(cbind(Clusters, rbind(prop_clust), Power, Type_I_err))
    names(out)[2:(length(mu_vec) + 1)] = paste("Prop", 1:length(mu_vec))
    row.names(out) = c()
    
    return(out)
  }
  
  true_phi_vec = ICC_to_shape(mu_vec, ICC)
  
  num_treatments = num_treatments_old = length(mu_vec)
  Trts = 1:num_treatments
  N_Clusts_remaining = N_clusts_max
  mu_vec_new = mu_vec
  Trts_new = Trts
  true_phi_vec_new = true_phi_vec
  prob_superiority = rep(1/num_treatments, num_treatments)
  Outcome = rep('', num_treatments)
  
  i = 1
  Winner = NA
  Futile = integer()
  data_full = data.frame(Trt = integer(),
                         n = integer(),
                         y = integer())
  
  while(N_Clusts_remaining > 0 & i <= N_analyses & length(Trts) > 1){
    N_clusts = ceiling(N_clusts_max/N_analyses/num_treatments)
    N_clust_new = ifelse(N_Clusts_remaining <= N_clusts*num_treatments,
                         floor(N_Clusts_remaining/num_treatments),
                         N_clusts)
    N_clust_vec = rep(N_clust_new, num_treatments)
    
    data_new = sample_NB_data(Trts_new, N_clust_vec, mean_clust_size,
                              true_phi_vec_new, mu_vec_new)
    data_full = rbind(data_full, data_new)
    
    phi_vec = cluster_empirical_bayes_shape(data_full, beta_a, beta_b)
    
    post_samps = cluster_posterior_samps(data_full %>% 
                                           filter(Trt %in% Trts_new), 
                                         phi_vec, beta_a, beta_b, N_MC)
    
    prob_superiority[Trts_new] = cluster_post_superiority(post_samps)
    
    N_Clusts_remaining = N_Clusts_remaining - sum(N_clust_vec)
    i = i + 1
    
    sup_candidate = which.max(prob_superiority)
    if(prob_superiority[sup_candidate] > sup_thresh){
      Winner = sup_candidate
      Outcome[sup_candidate] = 'Superiority'
      break
    } else{ 
      futiles = Trts[which(prob_superiority < fut_thresh)]
      if(length(futiles) > 0){
        Futile = unique(c(Futile, futiles))
        mu_vec_new = mu_vec[-Futile]
        Trts_new = Trts[-Futile]
        true_phi_vec_new = true_phi_vec[-Futile]
        phi_vec = phi_vec[-Futile]
        num_treatments = num_treatments_old - length(Futile)
        Outcome[Futile] = 'Early inferiority'
      }
    }
  } 
  
  out = data_full %>% 
    group_by(Trt) %>% 
    summarise(Clusters = n(),
              Patients = sum(n),
              Cases = sum(y)) %>% 
    as.data.frame() %>% 
    mutate(Outcome = Outcome)
  
  
  op_char = oper_char_single_cluster_trial(out, mu_vec)
  
  return(list(out  = out, op_char = op_char))
}





#*******************************************************************************
#* 6. Plotting the posterior probability of superiority by analysis for all 
#* arms.
#* 
#* Inputs: 
#*   simulation - the output of 'single_cluster_trial_simulation'
#*   sup_thresh - superiority threshold (on the posterior probability of 
#*                superiority)
#*   fut_thresh - futility threshold (on the posterior probability of 
#*                superiority)
#*   
#* Output: a ggplot2 object
#*******************************************************************************
cluster_prob_sup_plot = function(simulation, fut_thresh, sup_thresh){
  prob_sup = simulation$prob_sup
  
  df = prob_sup %>% 
    as.data.frame() %>%
    `colnames<-`(as.character(1:ncol(prob_sup))) %>% 
    reshape2::melt(value.name = "p_sup", 
                   variable.name = "Arm", id.vars = NULL) %>% 
    group_by(Arm) %>% 
    mutate(Analysis = 0:(n() - 1)) %>% 
    as.data.frame() %>% 
    na.omit()
  
  df_points = df %>% 
    group_by(Arm) %>% 
    filter(p_sup > sup_thresh | p_sup < fut_thresh) %>% 
    as.data.frame()
  
  df_other = setdiff(df, df_points)
  
  ggplot(df, aes(Analysis, p_sup, col = Arm)) + 
    geom_hline(yintercept = fut_thresh, linetype = 'dashed') + 
    geom_hline(yintercept = sup_thresh, linetype = 'dashed') +
    geom_line(linewidth = 1) + 
    geom_point(df_points,
               mapping = aes(Analysis, p_sup, col = Arm),
               shape = 21, size = 2.5, fill="white", stroke = 2)  +
    geom_point(df_other, 
               mapping = aes(Analysis, p_sup, col = Arm, fill = Arm),
               shape = 21, size = 1, stroke = 2) + 
    scale_x_continuous(expand = expansion(add = c(0.025, .075)),
                       breaks = 0:nrow(prob_sup)) + 
    scale_y_continuous(expand = expansion(add = c(0.035, .05)),
                       limits = c(0, 1)) + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic"),
          strip.text.x = element_text(size = 12, face = "bold")) + 
    scale_colour_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$Arm))]) + 
    scale_fill_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$Arm))]) + 
    ylab('Posterior Probability of Superiority')
}




#*******************************************************************************
#* 7. Plotting the posterior densities of the different mean event rates.
#* 
#* Inputs: 
#*   simulation - the output of 'single_cluster_trial_simulation'
#*   phi_vec - a vector of (empirical Bayes estimated) shape parameters 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   LL - x-axis lower limit
#*   UL - x-axis upper limit
#*   delta - grid resolution
#*   
#* Output: a ggplot2 object
#*******************************************************************************
cluster_post_density_plot = function(simulation, phi_vec, beta_a, 
                                     beta_b, LL, UL, delta){
  data = simulation$data
  
  tab = data %>% 
    group_by(Trt) %>% 
    summarise(clust_size = mean(n),
              sum_y = sum(y),
              n_clust = n()) %>% 
    dplyr::select(-Trt)
  
  tab = cbind(tab, phi_vec)
  
  mu_vec = seq(log(LL), log(UL), by = delta)
  
  dens = function(par_vec){
    n = as.numeric(par_vec[1])
    sum_y = as.numeric(par_vec[2])
    K = as.numeric(par_vec[3])
    phi = as.numeric(par_vec[4])
    
    out = (beta_a + sum_y)*log(n*exp(mu_vec)) + (beta_b + K*phi)*log(phi)
    out = out - lbeta(beta_a + sum_y, beta_b + K*phi)
    out = out - (beta_a + beta_b + sum_y + K*phi)*log(phi + n*exp(mu_vec))
    
    return(exp(out))
  }
  
  mu_dens = c(apply(tab, 1, dens))
  
  df = data.frame(mu = rep(mu_vec, length(phi_vec)),
                  Arm = factor(rep(1:length(phi_vec), each = length(mu_vec))),
                  mu_dens = mu_dens)
  
  ticks = seq(log(LL), log(UL), length = 4)
  labs = paste0(round(exp(ticks)*100), '%')
  
  ggplot(df, aes(mu, mu_dens, col = Arm)) + 
    geom_line() + 
    scale_x_continuous(labels = labs, breaks = ticks) + 
    geom_ribbon(aes(ymin = 0, ymax = mu_dens, fill = Arm), alpha = .5) + 
    scale_y_continuous(expand = expansion(c(0, .025))) + 
    scale_colour_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$Arm))]) + 
    scale_fill_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$Arm))]) + 
    ylab('Posterior density') + 
    xlab('Mean event rate') + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic")) 
}



#*******************************************************************************
#* 8. Summarizing the results of a single random CRT.
#* 
#* Inputs: 
#*   simulation - the output of 'single_cluster_trial_simulation'
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   
#* Output: a table containing the numbers of clusters, patients and cases for 
#*         each arm in each analysis, the end outcome (Early inferiority/
#*         Superiority/none) and posterior medians and 95% CrIs for the mean
#*         event rate parameters.
#*******************************************************************************
single_cluster_simulation_results = function(simulation, beta_a, beta_b){
  CrIs = mu_posterior_quantile(beta_a, beta_b, simulation)
  
  out = simulation$out %>% 
    mutate(`Mean event rate<br>[95% CrI]` = CrIs)
  
  return(out)
}



#*******************************************************************************
#* 9. Parallel Simulation of a large number of random CRTs.
#* 
#* Inputs: 
#*   mu_vec - a vector of mean event rates
#*   N_clusts_max - maximum overall number of clusters
#*   N_analyses - maximum number of analyses (interim + final)
#*   mean_clust_size - mean cluster size (for the Poisson distribution)
#*   true_phi_vec - a vector of true shape parameters 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   N_MC - size of Monte Carlo samples to be drawn
#*   sup_thresh - superiority threshold (on the posterior probability of 
#*                superiority)
#*   fut_thresh - futility threshold (on the posterior probability of 
#*                superiority)
#*   N_simulations - number of random trials to simulate
#*   
#* Output: a data frame containing the numbers of clusters, patients, cases and 
#*         the end outcome (Early inferiority/Superiority/none) for each arm at
#*         the end of each simulated trial
#*******************************************************************************
cluster_trial_simulation = function(mu_vec, N_clusts_max,
                                    N_analyses, 
                                    mean_clust_size,
                                    ICC, beta_a, beta_b, 
                                    N_MC, sup_thresh,
                                    fut_thresh,
                                    N_simulations){
  
  cl <- makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = N_simulations, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  l = foreach(i = 1:N_simulations,
              .options.snow = opts,
              .packages = c('dplyr'),
              .export = c('single_cluster_trial_simulation_fast',
                          'cluster_post_superiority',
                          'cluster_posterior_samps',
                          'sample_NB_ %>% ',
                          'cluster_log_marginal_likelihood',
                          'cluster_empirical_bayes_shape',
                          'ICC_to_shape')
  ) %dopar% {
    set.seed(i)
    out = single_cluster_trial_simulation_fast(mu_vec, N_clusts_max,
                                               N_analyses, 
                                               mean_clust_size,
                                               ICC, beta_a, beta_b, 
                                               N_MC, sup_thresh,
                                               fut_thresh)
    
    out$out$Simulation = out$op_char$Simulation = i
    
    out
  }
  
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  return(do.call(Map, c(f = rbind, l)))
}




#*******************************************************************************
#* 10. Evaluating the operating characteristics of a Bayesian adaptive CRT.
#* 
#* Inputs: 
#*   mu_vec - a vector of mean event rates
#*   N_clusts_max - maximum overall number of clusters
#*   N_analyses - maximum number of analyses (interim + final)
#*   mean_clust_size - mean cluster size (for the Poisson distribution)
#*   true_phi_vec - a vector of true shape parameters 
#*   beta_a - first shape parameter for the beta prior
#*   beta_b - second shape parameter for the beta prior
#*   N_MC - size of Monte Carlo samples to be drawn
#*   sup_thresh - superiority threshold (on the posterior probability of 
#*                superiority)
#*   fut_thresh - futility threshold (on the posterior probability of 
#*                superiority)
#*   N_simulations - number of random trials to simulate
#*   
#* Output: a table summarizing the mean no. of clusters, the proportion of 
#*         patients who received each treatment, and the detection and error 
#*         rates under both the null and the alternative scenario.
#*******************************************************************************
full_cluster_trial_simulation = function(mu_matrix, N_clusts_max,
                                         N_analyses, 
                                         mean_clust_size,
                                         ICC, beta_a, beta_b, 
                                         N_MC, sup_thresh,
                                         fut_thresh,
                                         N_simulations){
  
  
  tic("\nTime to complete")
  sim = apply(mu_matrix, 1, cluster_trial_simulation,
              N_clusts_max = N_clusts_max,
              N_analyses = N_analyses,
              mean_clust_size = mean_clust_size,
              ICC = ICC,
              beta_a = beta_a, beta_b = beta_b,
              N_MC = N_MC, sup_thresh = sup_thresh,
              fut_thresh = fut_thresh,
              N_simulations = N_simulations)
  
  
  operating_characs_func = function(i){
    oper_charac_df = sim[[i]]$op_char
    mu = mu_matrix[i,]
    
    apply(oper_charac_df, 2, mean) %>% 
      rbind() %>% 
      as.data.frame() %>% 
      dplyr::select(-Simulation) %>% 
      mutate(Type_I_err = ifelse(sum(mu == min(mu)) > 1,
                                 Type_I_err + Power, Type_I_err),
             Power = ifelse(sum(mu == min(mu)) > 1,
                            0, Power)
      ) %>% 
      mutate(Error_rate = Type_I_err,
             Clusters = sprintf('%.1f', Clusters)) %>% 
      rename('Detection_rate' = 'Power') %>% 
      dplyr::select(-Type_I_err) %>% 
      mutate_at(vars(!matches('Clusters')), 
                function(x){
                  paste0(sprintf('%.1f', x*100), '%')
                }) %>% 
      mutate(mu_vec = paste0('(',
                             paste(paste0(sprintf('%.1f', mu_matrix[i,]*100), '%'), 
                                   collapse = ', '),
                             ')'),
             Detection_rate = ifelse(sum(mu == min(mu)) > 1,
                                     '--', Detection_rate)) %>% 
      relocate(mu_vec, .before = Clusters)
  }
  
  results = do.call(rbind, lapply(1:nrow(mu_matrix), 
                                  operating_characs_func)) %>% 
    rename('Mean<br>event rates' = 'mu_vec',
           'Mean no.<br>of clusters' = 'Clusters',
           'TPR' = 'Detection_rate',
           'FPR' = 'Error_rate')
  
  names(results)[3:c(3 + ncol(mu_matrix) - 1)] = paste0('% of patients<br>in Arm ',
                                                        1:ncol(mu_matrix))
  row.names(results) = c()
  
  toc()
  
  return(results)
}




#*******************************************************************************
#* 11. Matching the shape parameter to a given intracluster correlation 
#* coefficient and a mean event rate.
#* 
#* Inputs: 
#*   mu_vec - a vector of mean event rates
#*   ICC - the assumed intracluster correlation coefficient
#*   
#* Output: a vector of shape parameters
#*******************************************************************************
ICC_to_shape = function(mu_vec, ICC){
  mu_vec/(1 - mu_vec)/ICC
}




#*******************************************************************************
#* 12. Matching the second beta prior shape parameter to a given intracluster 
#* correlation coefficient, the mean of the mean event rate prior, the first 
#* beta shape parameter, and the cluster size.
#* 
#* Inputs: 
#*   a - the first beta prior shape parameter
#*   mu_mean - a vector of mean event rates
#*   ICC - the assumed intracluster correlation coefficient
#*   n - cluster size
#*   
#* Output: the second beta prior shape parameter.
#*******************************************************************************
beta_prior_b_hyperpar = function(a, mu_mean, ICC, n){
  b = 1 + a/(n*ICC*(1 - mu_mean))
  
  return(b)
}



#*******************************************************************************
#* 13. Calculating quantiles of the mean event rate prior distribution.
#* 
#* Inputs: 
#*   u - the desired quantile
#*   beta_a - the first beta prior shape parameter
#*   beta_b - the second beta prior shape parameter
#*   n - cluster size
#*   mu_mean - a vector of mean event rates
#*   ICC - the assumed intracluster correlation coefficient
#*   
#* Output: the u-th quantile of the prior on the mean event rate scale.
#*******************************************************************************
mu_prior_quantile = function(u, beta_a, beta_b, n, mu_mean, ICC){
  phi = ICC_to_shape(mu_mean, ICC)
  p_quant = qbeta(u, beta_a, beta_b)
  mu_quant = phi/n*p_quant/(1 - p_quant)
  
  return(mu_quant)
}



#*******************************************************************************
#* 14. Potting the (log-) mean event rate prior distribution for various beta 
#* prior choices.
#* 
#* Inputs: 
#*   mu_mean - a vector of mean event rates
#*   beta_a - a vector of first beta prior shape parameters
#*   ICC - the assumed intracluster correlation coefficient
#*   n - cluster size
#*   LL - lowe x-axis limit
#*   UL - upper x-axis limit
#*   delta - grid resolution
#*   
#* Output: a ggplot2 object.
#*******************************************************************************
mu_prior_plot = function(mu_mean, beta_a, ICC, n, LL, UL, delta){
  beta_b = beta_prior_b_hyperpar(beta_a, mu_mean, ICC, n)
  hyperpars = cbind(beta_a, beta_b)
  phi = ICC_to_shape(mu_mean, ICC)
  pars = cbind(hyperpars, phi)
  
  mu_vec = seq(log(LL), log(UL), by = delta)
  
  dens = function(par_vec){
    a = par_vec[1]
    b = par_vec[2]
    phi = par_vec[3]
    
    (n*exp(mu_vec))^a*phi^b/beta(a, b)/(phi + n*exp(mu_vec))^(a + b)
  }
  
  a = rep(hyperpars[,1], each = length(mu_vec))
  b = rep(hyperpars[,2], each = length(mu_vec))
  levs = paste0('(a,b) = (', sprintf('%.2f', beta_a),
                ',', sprintf('%.2f', beta_b), ')')
  
  ticks = seq(log(LL), log(UL), length = 4)
  labs = paste0(round(exp(ticks)*100), '%')
  
  mu_dens = c(apply(pars, 1, dens))
  
  
  df = data.frame(cbind(a, b, mu = rep(mu_vec, nrow(pars)), mu_dens)) %>% 
    mutate(`(a,b)` = factor(paste0('(a,b) = (', sprintf('%.2f', a),
                                   ',', sprintf('%.2f', b), ')'),
                            levels = levs))
  
  ggplot(df, aes(mu, mu_dens, col = `(a,b)`)) + 
    geom_line() + 
    scale_x_continuous(labels = labs, breaks = ticks) + 
    geom_ribbon(aes(ymin = 0, ymax = mu_dens, fill = `(a,b)`), alpha = .5) + 
    scale_y_continuous(expand = expansion(c(0, .025))) + 
    scale_colour_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$`(a,b)`))]) + 
    scale_fill_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$`(a,b)`))]) + 
    ylab('Prior density') + 
    xlab('Mean event rate') + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          legend.title=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic")) 
}




#*******************************************************************************
#* 15. Summarizing the mean event rate prior distribution for various 
#* configurations of beta prior choices and other inputs in a table.
#* 
#* Inputs: 
#*   beta_a_vec - a vector of first beta prior shape parameters
#*   n_vec - a vector of cluster sizes
#*   mu_mean_vec - a vector of mean event rate parameter prior means
#*   ICC_vec - a vector of intracluster correlation coefficients
#*   
#* Output: a table summarizing the mean, median shape parameter and 95% CrI for 
#*         each configuration.
#*******************************************************************************
mu_prior_summary = function(beta_a_vec, n_vec, mu_mean_vec, ICC_vec){
  single_mu_prior_summary = function(beta_a, n, mu_mean, ICC){
    phi = ICC_to_shape(mu_mean, ICC)
    beta_b = beta_prior_b_hyperpar(beta_a, mu_mean, ICC, n)
    CrI = sapply(c(.025, .975), mu_prior_quantile,
                 beta_a = beta_a, beta_b = beta_b, n = n, 
                 mu_mean = mu_mean, ICC = ICC)
    
    out = data.frame(
      `$\\mathbb{E}[\\mu]$` = sprintf('%.2f', mu_mean),
      `$\\rho$` = sprintf('%.2f', ICC),
      `$\\phi$` = sprintf('%.2f', phi),
      `$(a,b)$` = paste0('$(', sprintf('%.2f', beta_a),
                         ', ', sprintf('%.2f', beta_b), ')$'),
      `95% CrI` = paste0('[', sprintf('%.2f', CrI[1]*100), '%; ',
                         sprintf('%.2f', CrI[2]*100), '%]'),
      check.names = FALSE
    )
    
    return(out)
  }
  
  mu_prior_summary_vec = function(v){
    beta_a = v[1]
    n = v[2] 
    mu_mean = v[3] 
    ICC = v[4]
    
    single_mu_prior_summary(beta_a, n, mu_mean, ICC)
  }
  
  df = expand.grid(beta_a_vec, n_vec, mu_mean_vec, ICC_vec)
  out = do.call(rbind, apply(df, 1, mu_prior_summary_vec))
  
  return(out)
}




#*******************************************************************************
#* 16. The log-marginal likelihood function.
#* 
#* Inputs: 
#*   data - the data emerging from a CRT, such as the output of 'sample_NB_data'
#*   beta_a - the first beta prior shape parameter
#*   beta_b - the second beta prior shape parameter
#*   phi - the shape parameter
#*   
#* Output: the log-marginal likelihood evaluated at phi.
#*******************************************************************************
cluster_log_marginal_likelihood = function(data, beta_a, beta_b, phi){
  K = nrow(data)
  
  l = lbeta(beta_a + sum(data$y), beta_b + K*phi) + 
    sum(lgamma(data$y + phi)) - K*lgamma(phi)
  
  return(-l)
}




#*******************************************************************************
#* 17. Empirical Bayes estimation of the shape parameters by maximizing the 
#* log-marginal likelihood.
#* 
#* Inputs: 
#*   data - the data emerging from a CRT, such as the output of 'sample_NB_data'
#*   beta_a - the first beta prior shape parameter
#*   beta_b - the second beta prior shape parameter
#*   
#* Output: empirical Bayes estimates of all shape parameters.
#*******************************************************************************
cluster_empirical_bayes_shape = function(data, beta_a, beta_b){
  one_group = function(i){
    phi_hat = suppressWarnings(nlminb(.5, cluster_log_marginal_likelihood,
                                      data = data %>% 
                                        filter(Trt == i), 
                                      beta_a = beta_a, 
                                      beta_b = beta_b)$par)
    
    return(phi_hat)
  }
  
  sapply(unique(data$Trt), one_group)
}




#*******************************************************************************
#* 18. Log-marginal likelihood curves with true shape parameters and their empirical 
#* Bayes estimates
#* 
#* Inputs: 
#*   data - the data emerging from a CRT, such as the output of 'sample_NB_data'
#*   beta_a - the first beta prior shape parameter
#*   beta_b - the second beta prior shape parameter
#*   true_phi - the true shape parameter that was used to generate the data
#*   LL - lower x-axis limit
#*   UL - upper x-axis limit
#*   delta - grid resolution
#*   
#* Output: a ggplot2 object.
#*******************************************************************************
cluster_marginal_likelihood_plot_all = function(data, beta_a, 
                                                beta_b, true_phi,
                                                LL, UL, delta){
  phi_hat = cluster_empirical_bayes_shape(data, beta_a, beta_b)
  
  phi_vec = seq(LL, UL, by = delta)
  
  Arms = unique(data$Trt)
  df = data.frame(Arm = rep(Arms, each = length(phi_vec)),
                  phi = rep(phi_vec, length(Arms)))
  
  aux_func = function(Arm){
    sapply(phi_vec, cluster_log_marginal_likelihood,
           data = data %>% filter(Trt == Arm), 
           beta_a = beta_a, 
           beta_b = beta_b)
  }
  
  
  aux_func2 = function(Arm){
    sapply(phi_hat[Arm], cluster_log_marginal_likelihood,
           data = data %>% filter(Trt == Arm), 
           beta_a = beta_a, 
           beta_b = beta_b)
  }
  
  
  aux_func3 = function(Arm){
    sapply(true_phi[Arm], cluster_log_marginal_likelihood,
           data = data %>% filter(Trt == Arm), 
           beta_a = beta_a, 
           beta_b = beta_b)
  }
  
  df$l = c(sapply(Arms, aux_func))
  df$Arm = paste0('Arm ', df$Arm)
  
  df_point_line = data.frame(
    Arm = rep(Arms, 2),
    phi = c(phi_hat, true_phi),
    
    lmin = rep(df %>% 
                 group_by(Arm) %>% 
                 summarise(lmin = -max(l, na.rm = T)) %>% 
                 as.data.frame() %>% 
                 dplyr::select(lmin) %>% 
                 unlist(), 2),
    
    l = c(sapply(Arms, aux_func2), sapply(Arms, aux_func3)),
    
    Shape = rep(c('Estimated', 'True'),
                each = length(Arms))
  ) %>% 
    mutate(Arm = paste0('Arm ', Arm))
  
  ggplot(df, aes(phi, -l)) + 
    geom_line(linewidth = 1, colour = '#CD534CFF') + 
    geom_segment(df_point_line,
                 mapping = aes(x = phi, xend = phi,
                               y = lmin, yend = -l,
                               col = Shape),
                 linetype = 'dashed', linewidth = .75) + 
    geom_point(df_point_line, mapping = aes(phi, -l, col = Shape),
               shape = 21, size = 2.5, fill="white", stroke = 2) + 
    facet_rep_wrap(~Arm, repeat.tick.labels = TRUE,
                   scales = 'free_y') + 
    scale_x_continuous(expand = expansion(c(0.01, 0.01))) + 
    scale_y_continuous(expand = expansion(c(0, .035))) + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          legend.text=element_text(size=10),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=18, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic"),
          strip.text.x = element_text(size = 12, face = "bold")) + 
    labs(x = expression(phi), y = 'Log-marginal Likelihood') + 
    scale_colour_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df$Arm))])
}




#*******************************************************************************
#* 19. Calculating point estimates ands 95% CrIs for the mean event rate 
#* parameters.
#* 
#* Inputs: 
#*   beta_a - the first beta prior shape parameter
#*   beta_b - the second beta prior shape parameter
#*   simulation - the output of 'single_cluster_trial_simulation'
#*   
#* Output: a vector of posterior medians and 95% CrIs.
#*******************************************************************************
mu_posterior_quantile = function(beta_a, beta_b, simulation){
  data = simulation$data
  phi = cluster_empirical_bayes_shape(data, beta_a, beta_b)
  
  tab = data %>% 
    group_by(Trt) %>% 
    summarise(clust_size = mean(n),
              sum_y = sum(y),
              n_clust = n()) %>% 
    dplyr::select(-Trt)
  
  tab = cbind(tab, phi) %>% 
    mutate(a_par = beta_a + sum_y,
           b_par = beta_b + n_clust*phi) 
  
  p_quant = qbeta(rep(c(.5, .025, .975), nrow(tab)), 
                  rep(tab$a_par, each = nrow(tab)), 
                  rep(tab$b_par, each = nrow(tab)))
  mu_quant = rep(tab$phi, each = nrow(tab))
  mu_quant = mu_quant/rep(tab$clust_size, each = nrow(tab))
  mu_quant = matrix(mu_quant*p_quant/(1 - p_quant), 
                    nrow = nrow(tab), byrow = TRUE) %>% 
    as.data.frame() %>% 
    mutate(`Mean event rate<br>[95% CrI]` = paste0(sprintf('%.1f', V1*100),
                                                   '% [', 
                                                   sprintf('%.1f', V2*100),
                                                   '%; ',
                                                   sprintf('%.1f', V3*100),
                                                   '%]')) %>% 
    dplyr::select(`Mean event rate<br>[95% CrI]`)
  
  return(mu_quant)
}