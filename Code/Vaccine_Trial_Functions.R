#*******************************************************************************
#*******************************************************************************
#* The functions used for the manuscript "Adaptive Trial Designs for Global 
#* Health: Vaccine and Cluster Randomized Trials go Bayesian" (2023).
#* 
#* All code was written by Ofir Harari
#*******************************************************************************
#*******************************************************************************




#*******************************************************************************
#* 1. Calculating the posterior probability of efficacy of each arm.
#* 
#* Inputs: 
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   n_events - a vector of event numbers (1st entry is control)
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   
#* Output: a vector of posterior probabilities of efficacy corresponding to 
#*         the different vaccine arms.
#*******************************************************************************
vaccine_post_efficacy = function(hyperpar, n_events, alloc_probs){
  alpha = hyperpar[1] + n_events
  alpha[1] = hyperpar[2] + n_events[1]
  
  post_effic = pbeta(alloc_probs[-1]/(alloc_probs[-1] + alloc_probs[1]),
                     alpha[-1], alpha[1])
  
  return(post_effic)
}




#*******************************************************************************
#* 2. Simulating a random Bayesian adaptive multi-arm, multi-stage vaccine 
#* trial.
#* 
#* Inputs: 
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   R - a vector of vaccine effects (not including control)
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a data frame containing the number of analyses, the stopping reason,
#*         and the eventual outcome (efficacious or not) of each vaccine arm.
#*******************************************************************************
single_vaccine_simulation = function(N_events_per_arm_per_analysis, 
                                     R, n_analyses, alloc_probs,
                                     hyperpar, futil_thresh, 
                                     effic_thresh, seed){
  set.seed(seed)
  
  r = c(0, R)
  K = length(r)
  N_event_per_analyses = ceiling(N_events_per_arm_per_analysis*K)
  Futile = integer()
  Efficacious = integer()
  active = rep(TRUE, K)
  events = matrix(0, ncol = K, nrow = 0)
  cum_events = rep(0, K)
  post_probs = matrix(.5, nrow = n_analyses + 1, ncol = K-1)
  alloc_probs_new = alloc_probs
  
  i = 0
  
  while(i < n_analyses){
    i = i + 1
    
    alloc_probs_new[active] = alloc_probs_new[active]/sum(alloc_probs_new[active])
    alloc_probs_new[!active] = 0
    theta = alloc_probs_new*(1 - r)/sum(alloc_probs_new*(1 - r))
    
    N_event_per_analyses = ceiling(N_events_per_arm_per_analysis*sum(active))
    events_new = c(rmultinom(1, N_event_per_analyses, theta))
    events = rbind(events, events_new)
    cum_events[active] = cum_events[active] + events_new[active]
    
    post_probs[i+1, active[-1]] = vaccine_post_efficacy(hyperpar, 
                                                        cum_events[active], 
                                                        alloc_probs_new[active])
    
    post_probs[i+1, !active[-1]] = post_probs[i, !active[-1]]
    
    Futile = which(post_probs[i+1,] <= futil_thresh[i]) + 1
    Efficacious = which(post_probs[i+1,] > effic_thresh[i]) + 1
    active[unique(c(Futile, Efficacious))] = FALSE
    
    if(sum(active) < 2) break
  }
  
  aux_func = function(j){
    futil_time = which(post_probs[2:n_analyses, j] < futil_thresh[1:(n_analyses - 1)])
    if(length(futil_time) > 0){
      final_analysis = min(futil_time)
      stopping_reason = 'Early futility'
    } else{
      effic_time = which(post_probs[2:n_analyses, j] > effic_thresh[1:(n_analyses - 1)])
      if(length(effic_time) > 0){
        final_analysis = min(effic_time)
        stopping_reason = 'Early efficacy'
      } else{
        final_analysis = n_analyses
        stopping_reason = 'End of trial'
      }
    }
    
    Efficacy = ifelse(any(post_probs[-1,j] > effic_thresh), 'Yes', 'No')
    
    data.frame(Seed = seed,
               Arm = j,
               final_analysis = final_analysis,
               stopping_reason = stopping_reason,
               Efficacy = Efficacy)
  }
  
  out = do.call(rbind, lapply(as.list(1:(K-1)), aux_func))
  
  return(out)
}




#*******************************************************************************
#* 3. Simulating a large number of random Bayesian adaptive multi-arm, 
#* multi-stage vaccine trials in parallel.
#* 
#* Inputs: 
#*   N_simulations - the numbers of trials to be simulated
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   R - a vector of vaccine effects (not including control)
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a data frame containing the number of analyses, the stopping reason,
#*         and the eventual outcome (efficacious or not) of each vaccine arm.
#*******************************************************************************
parallel_vaccine_simulation = function(N_simulations,
                                       N_events_per_arm_per_analysis, 
                                       R, n_analyses, 
                                       alloc_probs, hyperpar, 
                                       futil_thresh, effic_thresh){
  cl = makeSOCKcluster(parallel::detectCores() - 1)
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = N_simulations, style=3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress=progress)
  
  l = foreach(i = 1:N_simulations,
              .options.snow = opts,
              .export = c("single_vaccine_simulation",
                          "vaccine_post_efficacy")) %dopar% {
                            out = single_vaccine_simulation(N_events_per_arm_per_analysis, 
                                                            R, n_analyses, 
                                                            alloc_probs, 
                                                            hyperpar, 
                                                            futil_thresh, 
                                                            effic_thresh, 
                                                            seed = i)
                            
                            out
                          }
  close(pb)
  stopCluster(cl)
  registerDoSEQ()
  
  sim = do.call(rbind, l)
  
  Summary = sim %>%
    group_by(Arm) %>%
    summarise(
      power = mean(Efficacy == 'Yes')
    ) %>%
    as.data.frame
  
  
  return(list(Summary = Summary, sim = sim))
}




#*******************************************************************************
#* 4. Evaluating the operating characteristics of a Bayesian adaptive multi-arm,
#* multi-stage vaccine trial under the alternative scenario.
#* 
#* Inputs: 
#*   N_simulations - the numbers of trials to be simulated
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   R - a vector of vaccine effects (not including control)
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a data frame containing the number of analyses, the stopping reason,
#*         and the eventual outcome (efficacious or not) of each vaccine arm.
#*******************************************************************************
vaccine_trial_simulation = function(N_simulations,
                                    N_events_per_arm_per_analysis, 
                                    R, n_analyses, 
                                    alloc_probs, hyperpar, 
                                    futil_thresh, effic_thresh){
  R = sort(R, decreasing = T)
  Arm = 1:length(R)
  
  cat('\nSimulating under the alternative:\n')
  tic("\nTime to complete")
  power_sim = parallel_vaccine_simulation(N_simulations,
                                          N_events_per_arm_per_analysis, 
                                          R, n_analyses, 
                                          alloc_probs, hyperpar, 
                                          futil_thresh, effic_thresh)
  toc()
  cat("\n")
  
  operating_chars = power_sim$Summary %>%
    mutate(
      Arm = Arm,
      VE = paste0(100*R, "%"),
      Power = paste0(sprintf("%.1f", 100*power), "%")) %>%
    dplyr::select(Arm, VE, Power)
  
  
  Pr_stop = t(apply(table(power_sim$sim$Arm, 
                          power_sim$sim$final_analysis)/N_simulations*100, 
                    1, cumsum))
  Pr_stop = Pr_stop[,-ncol(Pr_stop)]
  Pr_stop = apply(Pr_stop, 2, 
                  function(x){
                    paste0(sprintf("%.1f", x), "%")
                  })
  
  colnames(Pr_stop) = paste0("Pr(Stop <= ", 1:ncol(Pr_stop), ")")
  colnames(Pr_stop) = gsub('Pr', '$\\\\mathbb{P}', colnames(Pr_stop))
  colnames(Pr_stop) = gsub('Stop <=', '\\\\mathrm{Stop}\\\\leq', 
                           colnames(Pr_stop))
  colnames(Pr_stop) = gsub(')', ')$', colnames(Pr_stop))
  colnames(Pr_stop) = gsub('\\\\leq 1', ' = 1', colnames(Pr_stop))
  
  operating_chars = cbind(operating_chars, Pr_stop)
  
  df = power_sim$sim %>% 
    mutate(Scenario = 'Alternative')
  
  
  return(list(operating_chars = operating_chars, df = df))
}



#*******************************************************************************
#* 5. Evaluating the operating characteristics of a Bayesian adaptive multi-arm,
#* multi-stage vaccine trial under the null scenario.
#* 
#* Inputs: 
#*   N_simulations - the numbers of trials to be simulated
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a data frame containing the number of analyses, the stopping reason,
#*         and the eventual outcome (efficacious or not) of each vaccine arm.
#*******************************************************************************
null_vaccine_trial_simulation = function(N_simulations,
                                         N_events_per_arm_per_analysis, 
                                         n_analyses, 
                                         alloc_probs, hyperpar, 
                                         futil_thresh, effic_thresh){
  R = rep(0, length(alloc_probs) - 1)
  Arm = 1:length(R)
  
  cat('\nSimulating under the null:\n')
  tic("\nTime to complete")
  power_sim = parallel_vaccine_simulation(N_simulations,
                                          N_events_per_arm_per_analysis, 
                                          R, n_analyses, 
                                          alloc_probs, hyperpar, 
                                          futil_thresh, effic_thresh)
  toc()
  cat("\n")
  
  power_sim$sim$VE = paste0("VE = ", 100*rep(R, N_simulations), "%")
  df = power_sim$sim %>% 
    mutate(Scenario = 'Null',
           Arm = 'All') %>% 
    dplyr::select(-VE)
  
  
  operating_chars = power_sim$Summary %>%
    mutate(
      Arm = 'All',
      VE = '0%'
    ) %>%
    group_by(Arm, VE) %>% 
    summarise(Power = paste0(sprintf("%.1f", 
                                     100*mean(power)), "%")) %>% 
    as.data.frame()
  
  
  Pr_stop = t(apply(table(df$Arm, df$final_analysis)/nrow(df)*100, 1, cumsum))
  Pr_stop = Pr_stop[,-ncol(Pr_stop)]
  Pr_stop = rbind(paste0(sprintf("%.1f", Pr_stop), "%"))
  colnames(Pr_stop) = paste0("Pr(Stop <= ", 1:ncol(Pr_stop), ")")
  
  colnames(Pr_stop) = paste0("Pr(Stop <= ", 1:ncol(Pr_stop), ")")
  colnames(Pr_stop) = gsub('Pr', '$\\\\mathbb{P}', colnames(Pr_stop))
  colnames(Pr_stop) = gsub('Stop <=', '\\\\mathrm{Stop}\\\\leq', 
                           colnames(Pr_stop))
  colnames(Pr_stop) = gsub(')', ')$', colnames(Pr_stop))
  colnames(Pr_stop) = gsub('\\\\leq 1', ' = 1', colnames(Pr_stop))
  
  operating_chars = cbind(operating_chars, Pr_stop)
  
  
  return(list(operating_chars = operating_chars, df = df))
}




#*******************************************************************************
#* 6. Evaluating the operating characteristics of a Bayesian adaptive multi-arm,
#* multi-stage vaccine trial under both the null and alternative scenarios.
#* 
#* Inputs: 
#*   N_simulations - the numbers of trials to be simulated
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   R - a vector of vaccine effects (not including control)
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a data frame containing the number of analyses, the stopping reason,
#*         and the eventual outcome (efficacious or not) of each vaccine arm.
#*******************************************************************************
full_vaccine_trial_simulation = function(N_simulations,
                                         N_events_per_arm_per_analysis, 
                                         R, n_analyses, 
                                         alloc_probs, hyperpar, 
                                         futil_thresh, effic_thresh){
  
  sim = vaccine_trial_simulation(N_simulations,
                                 N_events_per_arm_per_analysis, 
                                 R, n_analyses, 
                                 alloc_probs, hyperpar, 
                                 futil_thresh, effic_thresh)
  
  null_sim = null_vaccine_trial_simulation(N_simulations,
                                           N_events_per_arm_per_analysis, 
                                           n_analyses, 
                                           alloc_probs, hyperpar, 
                                           futil_thresh, effic_thresh)
  
  df = rbind(sim$df, null_sim$df)
  operating_chars = rbind(null_sim$operating_chars,
                          sim$operating_chars)
  
  return(list(operating_chars = operating_chars, df = df))
}




#*******************************************************************************
#* 7. Plotting the posterior probability of efficacy of each vaccine arm at each 
#* analysis.
#* 
#* Inputs: 
#*   df_post_probs - a data frame containing the arm number, the analysis 
#*                   number, and the posterior probability of efficacy for each.
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*                
#* Output: a ggplot2 object.
#*******************************************************************************
vaccine_post_efficacy_plot = function(df_post_probs, 
                                      effic_thresh,
                                      futil_thresh){
  
  vec = unique(df_post_probs$Analysis)
  vec[1] = 1
  if(length(vec) > 1){
    vec[length(vec)] = vec[length(vec) - 1]
  }
  
  df_futil = data.frame(
    Arm = factor(1),
    Analysis = unique(df_post_probs$Analysis)
  ) %>% 
    mutate(y = futil_thresh[vec])
  
  df_effic = data.frame(
    Arm = factor(1),
    Analysis = unique(df_post_probs$Analysis)
  ) %>% 
    mutate(y = effic_thresh[vec])
  
  df_Cross_up = suppressWarnings(df_post_probs %>%
                                   group_by(Arm) %>%
                                   summarise(
                                     Arm = Arm,
                                     Analysis = min(which(post_probs >= effic_thresh)) - 1
                                   )
  ) %>%
    unique() %>%
    merge(df_post_probs)
  
  n_arms = max(as.numeric(df_post_probs$Arm))
  df_futil_all = do.call(rbind, 
                         replicate(n_arms, 
                                   df_futil, 
                                   simplify = FALSE))
  df_futil_all$Arm = as.factor(rep(1:n_arms, 
                                   each = nrow(df_futil_all)/n_arms))
  
  df_Cross_down = suppressWarnings(df_post_probs %>% 
                                     left_join(df_futil_all %>%
                                                 mutate(Arm = as.integer(Arm)),
                                               by = c("Arm", 
                                                      "Analysis")) %>%
                                     filter(Analysis != 0) %>%
                                     group_by(Arm) %>%
                                     summarise(
                                       Arm = min(as.numeric(Arm)),
                                       Analysis = min(which(post_probs < y))
                                     )
  ) %>%
    merge(df_post_probs) %>% 
    mutate(Arm = as.factor(Arm),
           Analysis = as.integer(Analysis))
  
  df_Cross = rbind(df_Cross_up, df_Cross_down)
  df_Cross$Arm = as.factor(df_Cross$Arm)
  
  df_points = setdiff(df_post_probs %>% 
                        mutate(Arm = as.factor(Arm)), 
                      df_Cross)
  
  ggplot(df_post_probs %>% mutate(Arm = as.factor(Arm)), 
         aes(Analysis, post_probs, col = Arm)) + 
    geom_line(linewidth = 1) + 
    geom_point(df_points, 
               mapping = aes(Analysis, post_probs, col = Arm, fill = Arm),
               shape = 21, size = 1, stroke = 2) + 
    geom_step(data = df_futil, mapping = aes(Analysis, y), 
              linetype = 2, col = 1) + 
    geom_point(df_Cross, mapping = aes(Analysis, post_probs, col = Arm),
               shape = 21, size = 3, fill="white", stroke = 2) + 
    geom_hline(yintercept = effic_thresh, linetype = 2) +
    scale_x_continuous(expand = expansion(add = c(0.025, .075)),
                       breaks = unique(df_post_probs$Analysis)) + 
    scale_y_continuous(expand = expansion(add = c(0, .05)),
                       limits = c(0, 1)) +
    ylab("Posterior Probability of Efficacy") + 
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
          strip.text.x = element_text(size = 10, face = "bold")) + 
    scale_colour_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df_post_probs$Arm))]) + 
    scale_fill_manual(values = paletteer_d("ggsci::default_jco")[1:length(unique(df_post_probs$Arm))])
}



#*******************************************************************************
#* 8. Plotting the posterior VE density of each vaccine arm at the final 
#* analysis.
#* 
#* Inputs: 
#*   events - a matrix of the numbers of events observed in each arm between  
#*            every two interim analyses.
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   LL - lower x-axis limit for the plot
#*   UL - upper x-axis limit for the plot
#*   delta - x axis resolution (distance between any two points for plotting)
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*                                 
#* Output: a ggplot2 object.
#*******************************************************************************
VE_posterior_1D_density_plot = function(events, hyperpar,
                                        LL, UL, delta, 
                                        alloc_probs){
  tab = apply(events, 2, cumsum)
  tab_new = rbind(0, tab)
  a = hyperpar[1] + tab_new
  a[,1] = hyperpar[2] + tab_new[,1]
  K = ncol(a)
  
  temp = t(sapply(1:(K-1),
                  function(j){
                    v = a[,j+1]
                    last = max(which(!is.na(v)))
                    alpha = a[last, c(j+1, 1)]
                    
                    alloc_probs_new = alloc_probs[c(j+1, 1)]
                    alloc_probs_new = alloc_probs_new/sum(alloc_probs_new)
                    
                    return(c(alpha, alloc_probs_new))
                  }))
  
  alpha = temp[,1:2]
  trt_probs = temp[,3:4]
  
  colors = paletteer_d("ggsci::default_jco")[1:nrow(alpha)]
  labels = paste0('Arm ', 1:nrow(alpha))
  
  r_seq = seq(LL, UL, by = delta)
  r = rep(r_seq, nrow(alpha))
  
  df = cbind(r,
             matrix(rep(t(cbind(alpha, 1:nrow(alpha))), length(r_seq)), 
                    ncol = ncol(alpha) + 1, byrow = TRUE) %>% 
               as.data.frame() %>% 
               arrange(V3) %>% 
               dplyr::select(-V3) %>% 
               rename('alpha_trt' = 'V1',
                      'alpha_ctrl' = 'V2'),
             matrix(rep(t(cbind(trt_probs, 1:nrow(alpha))), length(r_seq)), 
                    ncol = ncol(alpha) + 1, byrow = TRUE) %>% 
               as.data.frame() %>% 
               arrange(V3) %>% 
               dplyr::select(-V3) %>% 
               rename('alloc_trt' = 'V1',
                      'alloc_ctrl' = 'V2')) 
  
  df = df %>% 
    mutate(density = apply(df, 1, VE_prior_1D),
           Arm = as.factor(rep(1:nrow(alpha), each = length(r_seq))))
  
  colors = paletteer_d("ggsci::default_jco")[1:nrow(alpha)]
  
  ggplot(df, aes(r, density, col = Arm, fill = Arm)) + 
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = density), alpha = c(.5)) + 
    scale_x_continuous(expand = expansion(c(0, 0.025)),
                       labels = scales::percent) + 
    scale_y_continuous(expand = expansion(c(0, 0.025))) + 
    scale_fill_manual(values = colors,
                      name = '',
                      labels = labels) +
    scale_color_manual(values = colors,
                       name = '',
                       labels = labels,
                       guide = "none") +
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic")) + 
    xlab('Vaccine Effect') + 
    ylab('Posterior Density')
}



#*******************************************************************************
#* 9. Calculating posterior VE medians and 95% CrIs for all vaccines.
#* 
#* Inputs: 
#*   events - a matrix of the numbers of events observed in each arm between  
#*            every two interim analyses.
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*                                 
#* Output: a vector of `VE [95% CrI]` elements.
#*******************************************************************************
VE_CrI_all = function(events, hyperpar, alloc_probs){
  tab = apply(events, 2, cumsum)
  tab_new = rbind(0, tab)
  a = hyperpar[1] + tab_new
  a[,1] = hyperpar[2] + tab_new[,1]
  
  
  VE_quantile = function(u, alpha, q){
    v = qbeta(1 - u, alpha[2], alpha[1])
    
    1 - q[1]*v/(q[2]*(1 - v))
  }
  
  RRR_CrI = function(i){
    last = max(which(!is.na(a[,i+1])))
    a_vec = a[last,c(1, i+1)]
    q = alloc_probs[c(1, i+1)]
    
    sapply(c(.025, .5, .975), VE_quantile,
           alpha = a_vec, q = q)
  }
  
  t(sapply(1:(ncol(events) - 1), RRR_CrI))
}




#*******************************************************************************
#* 10. Simulating a random Bayesian adaptive multi-arm, multi-stage vaccine 
#* trial for display purposes.
#* 
#* Inputs: 
#*   N_events_per_arm_per_analysis - the average number of events observed per  
#*                                   active (not yet stopped) arm for an interim
#*                                   analysis to take place
#*   R - a vector of vaccine effects (not including control)
#*   n_analyses - maximum overall (interim + final) number of analyses
#*   alloc_probs - a vector of of allocation probabilities (1st entry is 
#*                 control)
#*   hyperpar - a vector of length 2 of Dirichlet prior hyperparameters 
#*              (2nd entry is control)
#*   futil_thresh - vector of futility thresholds (on the posterior probability 
#*                  of efficacy)
#*   effic_thresh - vector of efficacy thresholds (on the posterior probability 
#*                  of superiority)
#*                
#* Output: a list consisting of the following slots - 
#*   sum_tab - a summary table that includes the number of cases per analysis,
#*             the stopping reason, the eventual outcome (efficacious or not), 
#*             posterior VE median and 95% CrI for each arm.
#*   VE_post_plot - the output of 'VE_posterior_1D_density_plot'.
#*   p_post_effic - the output of 'vaccine_post_efficacy_plot'.
#*******************************************************************************
vaccine_trial_plots_and_table = function(N_events_per_arm_per_analysis, 
                                         R, n_analyses, 
                                         alloc_probs, hyperpar, 
                                         futil_thresh, effic_thresh, 
                                         LL, UL, delta, seed){
  
  set.seed(seed)
  
  r = c(0, R)
  K = length(r)
  N_event_per_analyses = ceiling(N_events_per_arm_per_analysis*K)
  Futile = integer()
  Efficacious = integer()
  active = rep(TRUE, K)
  events = matrix(0, ncol = K, nrow = 0)
  cum_events = rep(0, K)
  post_probs = matrix(.5, nrow = n_analyses + 1, ncol = K-1)
  alloc_probs_new = alloc_probs
  
  i = 0
  
  while(i < n_analyses){
    i = i + 1
    
    alloc_probs_new[active] = alloc_probs_new[active]/sum(alloc_probs_new[active])
    alloc_probs_new[!active] = 0
    theta = alloc_probs_new*(1 - r)/sum(alloc_probs_new*(1 - r))
    
    N_event_per_analyses = ceiling(N_events_per_arm_per_analysis*sum(active))
    events_new = c(rmultinom(1, N_event_per_analyses, theta))
    events = rbind(events, events_new)
    cum_events[active] = cum_events[active] + events_new[active]
    
    post_probs[i+1, active[-1]] = vaccine_post_efficacy(hyperpar, 
                                                        cum_events[active], 
                                                        alloc_probs_new[active])
    
    post_probs[i+1, !active[-1]] = post_probs[i, !active[-1]]
    
    Futile = which(post_probs[i+1,] <= futil_thresh[i]) + 1
    Efficacious = which(post_probs[i+1,] > effic_thresh[i]) + 1
    active[unique(c(Futile, Efficacious))] = FALSE
    
    if(sum(active) < 2) break
  }
  
  aux_func = function(j){
    futil_time = which(post_probs[2:n_analyses, j] < futil_thresh[1:(n_analyses - 1)])
    if(length(futil_time) > 0){
      final_analysis = min(futil_time)
      stopping_reason = 'Early futility'
    } else{
      effic_time = which(post_probs[2:n_analyses, j] > effic_thresh[1:(n_analyses - 1)])
      if(length(effic_time) > 0){
        final_analysis = min(effic_time)
        stopping_reason = 'Early efficacy'
      } else{
        final_analysis = n_analyses
        stopping_reason = 'End of trial'
      }
    }
    
    Efficacy = ifelse(any(post_probs[-1,j] > effic_thresh), 'Yes', 'No')
    
    data.frame(Seed = seed,
               Arm = j,
               final_analysis = final_analysis,
               stopping_reason = stopping_reason,
               Efficacy = Efficacy)
  }
  
  out = do.call(rbind, lapply(as.list(1:(K-1)), aux_func))
  
  df_post_probs = do.call(rbind, 
                          lapply(1:ncol(post_probs),
                                 function(j){
                                   data.frame(
                                     Analysis = 0:out$final_analysis[j],
                                     Arm = j,
                                     post_probs = post_probs[1:(out$final_analysis[j] + 1),j]
                                   )
                                 }
                          )
  )
  
  p_post_effic = vaccine_post_efficacy_plot(df_post_probs, effic_thresh, 
                                            futil_thresh)
  
  events[events == 0] = NA 
  VE_post_plot = VE_posterior_1D_density_plot(events, hyperpar,
                                                LL, UL, delta, 
                                                alloc_probs)
  
  
  RRRs = VE_CrI_all(events, hyperpar, alloc_probs)
  RRR_CrI_txt = paste0(sprintf("%.1f", 100*RRRs[,2]), "% [",
                       sprintf("%.1f", 100*RRRs[,1]), "%; ",
                       sprintf("%.1f", 100*RRRs[,3]), "%]"
  )
  
  events_display = t(apply(events, 2, cumsum))
  events_display[is.na(events_display)] = ''
  row.names(events_display) = c()
  colnames(events_display) = paste0("Analysis ", 
                                    1:ncol(events_display), 
                                    "\ncases")
  
  events_display = events_display %>% 
    as.data.frame()
  
  
  post_effic = df_post_probs %>% 
    group_by(Arm) %>% 
    filter(row_number() == n()) %>% 
    as.data.frame() %>% 
    dplyr::select(post_probs) %>% 
    unlist()
  
  sum_tab = cbind(rbind(data.frame(Arm = 'Control', 
                                   stopping_reason = '',
                                   Efficacy = ''),
                        out %>% 
                          dplyr::select(-c(Seed, final_analysis))),
                  events_display
  ) %>% 
    mutate(`VE [95% CrI]` = c('', RRR_CrI_txt)) %>% 
    rename('Stopping\nreason' = 'stopping_reason') %>% 
    relocate(Efficacy, .before = `VE [95% CrI]`) %>% 
    mutate(`Pr(Efficacy|Data)` = c('', 
                                   paste0(sprintf('%.1f', 
                                                  post_effic*100), '%'))
    ) %>% 
    relocate(`Pr(Efficacy|Data)`, .before = `VE [95% CrI]`)
  
  return(list(p_post_effic = p_post_effic, 
              VE_post_plot = VE_post_plot,
              sum_tab = sum_tab))
}




#*******************************************************************************
#* 11. Calculating the prior VE density.
#* 
#* Inputs: 
#*   v - a vector of length 4, whose 1st entry is the VE, 2nd and 3rd entries 
#*       are the vaccine and control Dirichlet hyperparameters, respectively, 
#*       and 4th entry is the vaccine group (normalized) allocation probability.   
#*                                 
#* Output: the prior density evaluated at the 1st element of v.
#*******************************************************************************
VE_prior_1D = function(v){
  r = v[1]
  alpha = v[2:3]
  trt_prob = v[4]
  
  coeff = trt_prob^alpha[1]*(1 - trt_prob)^alpha[2]/beta(alpha[1], alpha[2])
  d = coeff*(1-r)^(alpha[1]-1)/(1-trt_prob*r)^sum(alpha)
  
  return(d)
}



#*******************************************************************************
#* 12. Deriving the control Dirichlet hyperparameter to ensure prior probability 
#* of efficacy of 0.5, based on the vaccine hyperparameter and the treatment 
#* (normalized) allocation probability.
#* 
#* Inputs: 
#*   alpha_trt - the vaccine Dirichlet hyperparameter. 
#*   q - the vaccine (normalized) allocation probability.  
#*                                 
#* Output: the control Dirichlet hyperparameter.
#*******************************************************************************
alpha_ctrl_func = function(alpha_trt, q){
  uniroot(
    function(x){
      vaccine_post_efficacy(c(alpha_trt, x), 
                            c(0, 0), 
                            c(1-q, q)) - .5
    },
    lower = 0, upper = 100
  )$root
}



#*******************************************************************************
#* 13. Plotting the prior VE density of each vaccine arm.
#* 
#* Inputs: 
#*   LL - lower x-axis limit for the plot
#*   UL - upper x-axis limit for the plot
#*   delta - x axis resolution (distance between any two points for plotting)
#*   alpha - a 2-column matrix of Dirichlet prior hyperparameters (1st entry is 
#*           vaccine)
#*   trt_prob - vaccine (normalized) allocation probability
#*                                 
#* Output: a ggplot2 object.
#*******************************************************************************
VE_prior_1D_density_plot = function(LL, UL, delta, alpha, trt_prob){
  r_seq = seq(LL, UL, by = delta)
  r = rep(r_seq, nrow(alpha))
  
  df = cbind(r,
             matrix(rep(t(cbind(alpha, 1:nrow(alpha))), length(r_seq)), 
                    ncol = ncol(alpha) + 1, byrow = TRUE) %>% 
               as.data.frame() %>% 
               arrange(V3) %>% 
               dplyr::select(-V3) %>% 
               rename('alpha_trt' = 'V1',
                      'alpha_ctrl' = 'V2'),
             trt_prob) 
  
  
  df = df %>% 
    mutate(density = apply(df, 1, VE_prior_1D),
           alpha_txt = paste0("list(alpha[vax] == ",
                              round(alpha_trt, 2), 
                              ", ",
                              "alpha[ctrl] == ", 
                              round(alpha_ctrl, 2), ")")
    ) %>% 
    mutate(alpha_txt = factor(alpha_txt, levels = unique(alpha_txt)))
  
  colors = paletteer_d("ggsci::default_jco")[1:nrow(alpha)]
  
  ggplot(df, aes(r, density, col = alpha_txt, fill = alpha_txt)) + 
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = density), alpha = c(.5)) + 
    scale_x_continuous(expand = expansion(c(0, 0.025)),
                       labels = scales::percent) + 
    scale_y_continuous(expand = expansion(c(0, 0.025))) + 
    scale_fill_manual(values = colors,
                      name = '',
                      labels = parse(text = levels(df$alpha_txt))) + 
    scale_color_manual(values = colors,
                       name = '',
                       labels = parse(text = df$alpha_txt),
                       guide = "none") +
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic")) + 
    xlab('Vaccine Effect') + 
    ylab('Prior Density')
}



#*******************************************************************************
#* 14. Calculating prior VE quantiles.
#* 
#* Inputs: 
#*   u - the quantile to be calculated
#*   alpha - a length 2 vector of Dirichlet prior hyperparameters (1st entry 
#*           is vaccine)
#*   q - vaccine (normalized) allocation probability
#*                                 
#* Output: the u-th prior VE quantile.
#*******************************************************************************
VE_quantile = function(u, alpha, q){
  v = qbeta(1 - u, alpha[1], alpha[2])
  
  1 - (1-q)*v/(q*(1 - v))
}



#*******************************************************************************
#* 15. Calculating prior VE mean.
#* 
#* Inputs: 
#*   alpha - a length 2 vector of Dirichlet prior hyperparameters (1st entry 
#*           is vaccine)
#*   q - vaccine (normalized) allocation probability
#*                                 
#* Output: the prior VE mean.
#*******************************************************************************
VE_mean = function(alpha, q){
  1 - (1 - q)/q*alpha[1]/(alpha[2] - 1)
}



#*******************************************************************************
#* 16. Summarizing prior VE distributions.
#* 
#* Inputs: 
#*   alpha - a length 2 vector of Dirichlet prior hyperparameters (1st entry 
#*           is vaccine)
#*   q - vaccine (normalized) allocation probability
#*                                 
#* Output: prior VE distribution hyperparameters, mean, median, and 95% CrI.
#*******************************************************************************
VE_prior_summary = function(alpha, q){
  Mean = sprintf('%.2f', round(VE_mean(alpha, q), 2))
  Mean = ifelse(Mean == '-0.00', '0.00', Mean)
  
  Median = sprintf('%.2f', round(VE_quantile(.5, alpha, q), 2))
  Median = ifelse(Median == '-0.00', '0.00', Median)
  
  CrI = sapply(c(.025, .975), VE_quantile, alpha = alpha, q = q)
  CrI[1] = ifelse(CrI[1] == '-0.00', '0.00', CrI[1])
  CrI[2] = ifelse(CrI[2] == '-0.00', '0.00', CrI[2])
  CrI = paste0('[', sprintf('%.2f', round(CrI[1],2)),
               '; ', sprintf('%.2f', round(CrI[2], 2)), ']')
  
  data.frame('$\\alpha_{\\mathrm{vax}}$' = sprintf('%.2f', 
                                                   round(alpha[1], 2)), 
             '$\\alpha_{\\mathrm{ctrl}}$' = sprintf('%.2f', 
                                                    round(alpha[2], 2)),
             Mean = Mean, Median = Median,
             `95% CrI` = CrI,
             check.names = FALSE)
}




#*******************************************************************************
#* 17. D-optimal control/treatment allocation ratio.
#* 
#* Inputs: 
#*   p_ctrl - control arm event rate
#*   n_trts - number of (active) vaccine arms
#*   RRR - relative rate reduction
#*                                 
#* Output: a list containing the following slots - 
#*   rho - the control/treatment allocation ratio.
#*   plot - a ggplot of the log-determinant of the estimated treatment - control 
#*          contrast vector vs. rho.
#*******************************************************************************
D_optimal_alloc_ratio = function(p_ctrl, n_trts, RRR){
  n_total = 1000
  
  D_crit = function(rho){
    p_trt = p_ctrl*(1 - RRR)
    n_trt = n_total/(rho + n_trts)
    n_ctrl = n_trt*rho
    
    M = diag(p_trt*(1 - p_trt)/n_trt, n_trts)
    M = M + p_ctrl*(1 - p_ctrl)/n_ctrl
    L = chol(M)
    
    return(2*sum(log(diag(L))))
  }
  
  rho_vec = seq(.75, 1.5, by = .001)
  d = sapply(rho_vec, D_crit)
  
  opt = nlminb(1, D_crit)
  
  df = data.frame(rho = rho_vec, d = d)
  df_point = data.frame(rho = opt$par, d = opt$objective)
  
  p = ggplot(df, aes(rho, d)) +
    geom_line(col = 'royal blue', linewidth = 1) + 
    xlab(expression(rho)) + 
    ylab(expression(log~abs(Var(hat(c))))) +
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=18, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic")) + 
    scale_x_continuous(expand = expansion(c(0, 0.015))) + 
    scale_y_continuous(expand = expansion(c(0.15, 0))) + 
    geom_segment(df_point, 
                 mapping = aes(x = rho, xend = rho, y = -Inf, yend = d), 
                 col = 'dark red', linetype = 'dashed', linewidth = .75) + 
    geom_point(df_point, mapping = aes(rho, d), col = 'dark red',
               shape = 21, size = 2.5, stroke = 2, bg = 'white')
  
  return(list(rho = opt$par, plot = p))
}
