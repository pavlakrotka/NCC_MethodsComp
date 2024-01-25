# devtools::install_github("pavlakrotka/NCC", force = TRUE)
library(NCC)
library(tidyverse)

n_sim <- 10000

##############################################################################################################################################################

# FUNCTIONS BY DOMINIC

## given expected and 'maximum' jump in the mean between periods/buckets,
## what are the corresponding tau_a and tau_b?

from_jump_size_to_tau_a_b <- function(e_jump_size = 0.01, max_jump_size = 0.1){
  
  fun_to_solve <- function(tau_b, max_jump_size, e_jump_size){
    pgamma(1 / max_jump_size ^ 2, 1 / (e_jump_size ^ 2) * tau_b, tau_b) - 0.01
  }
  
  tau_b <- uniroot(fun_to_solve,
                   c(1e-12, 1e12),
                   max_jump_size = max_jump_size,
                   e_jump_size = e_jump_size)$root
  
  tau_a <- 1 / e_jump_size ^ 2 * tau_b
  
  c(tau_a = tau_a,
    tau_b = tau_b)
}

##############################################################################################################################################################


# SCENARIO I - 3 arms, equal and different time trends, vary d and lambda ----

# N_peak for different d's (always middle of the trial):

# sum(get_ss_matrix(num_arms = 3, n_arm = 250, d = 0*(0:2)), na.rm=T)/2
# sum(get_ss_matrix(num_arms = 3, n_arm = 250, d = 125*(0:2)), na.rm=T)/2
# sum(get_ss_matrix(num_arms = 3, n_arm = 250, d = 250*(0:2)), na.rm=T)/2
# sum(get_ss_matrix(num_arms = 3, n_arm = 250, d = 375*(0:2)), na.rm=T)/2
# sum(get_ss_matrix(num_arms = 3, n_arm = 250, d = 500*(0:2)), na.rm=T)/2

## EQUAL TIME TREND ----

### TYPE I ERROR ----

set.seed(1)

scenario_i_eq_alpha <- data.frame(num_arms = 3, 
                                  n_arm = 250, 
                                  d1 = seq(0, 500, length.out=5)*0,
                                  d2 = seq(0, 500, length.out=5)*1,
                                  d3 = seq(0, 500, length.out=5)*2,
                                  period_blocks = 2, 
                                  mu0 = 0,
                                  sigma = 1,
                                  theta1 = 0,
                                  theta2 = 0,
                                  theta3 = 0,
                                  lambda0 = 0.15, 
                                  lambda1 = 0.15, 
                                  lambda2 = 0.15, 
                                  lambda3 = 0.15,
                                  trend = rep(c("linear", "stepwise_2", "inv_u"), each = 5),
                                  N_peak = c(rep(NA, 10), c(500, 553, 605, 672, 750)),
                                  alpha = 0.025,
                                  prec_theta = 0.001,
                                  prec_eta = 0.001,
                                  tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                  tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                  tau_case = "reasonable",
                                  prec_a = 0.001,
                                  prec_b = 0.001,
                                  bucket_size = 25,
                                  opt = 2, 
                                  prior_prec_tau = 0.002, 
                                  prior_prec_eta = 1,
                                  n_samples = 1000,
                                  n_chains = 4,
                                  n_iter = 4000,
                                  n_adapt = 1000,
                                  robustify = TRUE,
                                  weight = 0.1)


results_i_eq_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_i_eq_alpha, arms = c(2:3), models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPpriorNew"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_eq_alpha, "results/results_i_eq_alpha.csv")


set.seed(2)

scenario_i_eq_alpha_TM_step <- data.frame(num_arms = 3, 
                                          n_arm = 250, 
                                          d1 = 250*0,
                                          d2 = 250*1,
                                          d3 = 250*2,
                                          period_blocks = 2, 
                                          mu0 = 0,
                                          sigma = 1,
                                          theta1 = 0,
                                          theta2 = 0,
                                          theta3 = 0,
                                          lambda0 = seq(-0.15, 0.15, length.out = 9), 
                                          lambda1 = seq(-0.15, 0.15, length.out = 9),
                                          lambda2 = seq(-0.15, 0.15, length.out = 9),
                                          lambda3 = seq(-0.15, 0.15, length.out = 9),
                                          trend = "stepwise_2",
                                          N_peak = NA,
                                          alpha = 0.025,
                                          prec_theta = 0.001,
                                          prec_eta = 0.001,
                                          tau_a = rep(c(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1],
                                                        from_jump_size_to_tau_a_b(e_jump_size = 0.001, max_jump_size = 0.015)[1],
                                                        from_jump_size_to_tau_a_b(e_jump_size = 10, max_jump_size = 15)[1]), each=9),
                                          tau_b = rep(c(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2],
                                                        from_jump_size_to_tau_a_b(e_jump_size = 0.001, max_jump_size = 0.015)[2],
                                                        from_jump_size_to_tau_a_b(e_jump_size = 10, max_jump_size = 15)[2]), each=9),
                                          tau_case = rep(c("reasonable", "small", "large"), each=9),
                                          prec_a = 0.001,
                                          prec_b = 0.001,
                                          bucket_size = 25,
                                          opt = NA, 
                                          prior_prec_tau = NA, 
                                          prior_prec_eta = NA,
                                          n_samples  = NA, 
                                          n_chains = NA,
                                          n_iter = NA, 
                                          n_adapt = NA, 
                                          robustify = NA, 
                                          weight = NA)

results_i_eq_alpha_TM_step <- sim_study_par(nsim = n_sim, scenarios = scenario_i_eq_alpha_TM_step, arms = c(3), models = c("timemachine", "fixmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_eq_alpha_TM_step, "results/results_i_eq_alpha_TM_step.csv")



set.seed(3)

scenario_i_eq_alpha_MAP <- data.frame(num_arms = 3, 
                                      n_arm = 250, 
                                      d1 = 250*0,
                                      d2 = 250*1,
                                      d3 = 250*2,
                                      period_blocks = 2, 
                                      mu0 = 0,
                                      sigma = 1,
                                      theta1 = 0,
                                      theta2 = 0,
                                      theta3 = 0,
                                      lambda0 = seq(-0.15, 0.15, length.out = 9), 
                                      lambda1 = seq(-0.15, 0.15, length.out = 9),
                                      lambda2 = seq(-0.15, 0.15, length.out = 9),
                                      lambda3 = seq(-0.15, 0.15, length.out = 9),
                                      trend = "stepwise_2",
                                      N_peak = NA,
                                      alpha = 0.025,
                                      prec_theta = NA,
                                      prec_eta = NA,
                                      tau_a = NA,
                                      tau_b = NA,
                                      tau_case = NA,
                                      prec_a = NA,
                                      prec_b = NA,
                                      bucket_size = NA,
                                      opt = 2,
                                      prior_prec_tau = rep(c(2, 0.2, 0.002), each=9),
                                      prior_prec_eta = c(rep(0.001, 9*3), rep(1, 9*3)),
                                      n_samples = 1000,
                                      n_chains = 4,
                                      n_iter = 4000,
                                      n_adapt = 1000,
                                      robustify = TRUE,
                                      weight = 0.1)

results_i_eq_alpha_MAP <- sim_study_par(nsim = n_sim, scenarios = scenario_i_eq_alpha_MAP, arms = c(3), models = c("MAPpriorNew", "fixmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_eq_alpha_MAP, "results/results_i_eq_alpha_MAP.csv")







### POWER ----


set.seed(4)

scenario_i_eq_pow <- data.frame(num_arms = 3, 
                                n_arm = 250, 
                                d1 = seq(0, 500, length.out=5)*0,
                                d2 = seq(0, 500, length.out=5)*1,
                                d3 = seq(0, 500, length.out=5)*2,
                                period_blocks = 2, 
                                mu0 = 0,
                                sigma = 1,
                                theta1 = 0.25,
                                theta2 = 0.25,
                                theta3 = 0.25,
                                lambda0 = 0.15, 
                                lambda1 = 0.15, 
                                lambda2 = 0.15, 
                                lambda3 = 0.15,
                                trend = rep(c("linear", "stepwise_2", "inv_u"), each = 5),
                                N_peak = c(rep(NA, 10), c(500, 553, 605, 672, 750)),
                                alpha = 0.025,
                                prec_theta = 0.001,
                                prec_eta = 0.001,
                                tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                tau_case = "reasonable",
                                prec_a = 0.001,
                                prec_b = 0.001,
                                bucket_size = 25,
                                opt = 2, 
                                prior_prec_tau = 0.002, 
                                prior_prec_eta = 1,
                                n_samples = 1000,
                                n_chains = 4,
                                n_iter = 4000,
                                n_adapt = 1000,
                                robustify = TRUE,
                                weight = 0.1)

results_i_eq_pow <- sim_study_par(nsim = n_sim, scenarios = scenario_i_eq_pow, arms = c(2:3), models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPpriorNew"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_eq_pow, "results/results_i_eq_pow.csv")





## DIFFERENT TREND IN ARM 1 ----


### TYPE I ERROR ----

set.seed(5)

scenario_i_diff1_alpha <- data.frame(num_arms = 3, 
                                     n_arm = 250, 
                                     d1 = 250*0,
                                     d2 = 250*1,
                                     d3 = 250*2,
                                     period_blocks = 2, 
                                     mu0 = 0,
                                     sigma = 1,
                                     theta1 = 0,
                                     theta2 = 0,
                                     theta3 = 0,
                                     lambda0 = 0.1,
                                     lambda1 = seq(-0.05, 0.25, length.out = 9),
                                     lambda2 = 0.1,
                                     lambda3 = 0.1,
                                     trend = rep(c("linear", "stepwise_2", "inv_u"), each = 9),
                                     N_peak = c(rep(NA, 18), rep(605, 9)),
                                     alpha = 0.025,
                                     prec_theta = 0.001,
                                     prec_eta = 0.001,
                                     tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                     tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                     tau_case = "reasonable",
                                     prec_a = 0.001,
                                     prec_b = 0.001,
                                     bucket_size = 25,
                                     opt = 2, 
                                     prior_prec_tau = 0.002, 
                                     prior_prec_eta = 1,
                                     n_samples = 1000,
                                     n_chains = 4,
                                     n_iter = 4000,
                                     n_adapt = 1000,
                                     robustify = TRUE,
                                     weight = 0.1)

results_i_diff1_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_i_diff1_alpha, arms = c(3), models = c("fixmodel", "sepmodel", "timemachine", "MAPpriorNew"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_diff1_alpha, "results/results_i_diff1_alpha.csv")



## DIFFERENT TREND IN ARMS 1 & 2 ----


### TYPE I ERROR ----

set.seed(6)

scenario_i_diff12_alpha <- data.frame(num_arms = 3, 
                                      n_arm = 250, 
                                      d1 = 250*0,
                                      d2 = 250*1,
                                      d3 = 250*2,
                                      period_blocks = 2, 
                                      mu0 = 0,
                                      sigma = 1,
                                      theta1 = 0,
                                      theta2 = 0,
                                      theta3 = 0,
                                      lambda0 = 0.1,
                                      lambda1 = seq(-0.05, 0.25, length.out = 9),
                                      lambda2 = seq(-0.05, 0.25, length.out = 9),
                                      lambda3 = 0.1,
                                      trend = rep(c("linear", "stepwise_2", "inv_u"), each = 9),
                                      N_peak = c(rep(NA, 18), rep(605, 9)),
                                      alpha = 0.025,
                                      prec_theta = 0.001,
                                      prec_eta = 0.001,
                                      tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                      tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                      tau_case = "reasonable",
                                      prec_a = 0.001,
                                      prec_b = 0.001,
                                      bucket_size = 25,
                                      opt = 2, 
                                      prior_prec_tau = 0.002, 
                                      prior_prec_eta = 1,
                                      n_samples = 1000,
                                      n_chains = 4,
                                      n_iter = 4000,
                                      n_adapt = 1000,
                                      robustify = TRUE,
                                      weight = 0.1)

results_i_diff12_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_i_diff12_alpha, arms = c(3), models = c("fixmodel", "sepmodel", "timemachine", "MAPpriorNew"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_i_diff12_alpha, "results/results_i_diff12_alpha.csv")


##############################################################################################################################################################

# SCENARIO II - 4 arms with non-equidistant entry time (varying d3), only equal time trend ----

# N_peak for different d's (always middle of the trial):

# sum(get_ss_matrix(num_arms = 4, n_arm = 250, d = c(0, 300, 300, 800)), na.rm = T)/2
# sum(get_ss_matrix(num_arms = 4, n_arm = 250, d = c(0, 300, 425, 800)), na.rm = T)/2
# sum(get_ss_matrix(num_arms = 4, n_arm = 250, d = c(0, 300, 550, 800)), na.rm = T)/2
# sum(get_ss_matrix(num_arms = 4, n_arm = 250, d = c(0, 300, 675, 800)), na.rm = T)/2
# sum(get_ss_matrix(num_arms = 4, n_arm = 250, d = c(0, 300, 800, 800)), na.rm = T)/2


## TYPE I ERROR ----

set.seed(7)

scenario_ii_eq_alpha <- data.frame(num_arms = 4, 
                                   n_arm = 250, 
                                   d1 = 0,
                                   d2 = 300,
                                   d3 = seq(300, 800, length.out=5),
                                   d4 = 800,
                                   period_blocks = 2, 
                                   mu0 = 0,
                                   sigma = 1,
                                   theta1 = 0,
                                   theta2 = 0,
                                   theta3 = 0,
                                   theta4 = 0,
                                   lambda0 = 0.15,
                                   lambda1 = 0.15,
                                   lambda2 = 0.15,
                                   lambda3 = 0.15,
                                   lambda4 = 0.15,
                                   trend = rep(c("linear", "stepwise_2", "inv_u"), each = 5),
                                   N_peak = c(rep(NA, 10), c(767, 774, 781, 790, 800)),
                                   alpha = 0.025,
                                   prec_theta = 0.001,
                                   prec_eta = 0.001,
                                   tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                   tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                   tau_case = "reasonable",
                                   prec_a = 0.001,
                                   prec_b = 0.001,
                                   bucket_size = 25,
                                   opt = 2, 
                                   prior_prec_tau = 0.002, 
                                   prior_prec_eta = 1,
                                   n_samples = 1000,
                                   n_chains = 4,
                                   n_iter = 4000,
                                   n_adapt = 1000,
                                   robustify = TRUE,
                                   weight = 0.1)

results_ii_eq_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_ii_eq_alpha, arms = c(4), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_ii_eq_alpha, "results/results_ii_eq_alpha.csv")


## POWER ----

# power.t.test(n=250, sd=1, sig.level = 0.025, power=0.8, type = "two.sample", alternative = "one.sided")

set.seed(8)

scenario_ii_eq_pow <- data.frame(num_arms = 4, 
                                 n_arm = 250, 
                                 d1 = 0,
                                 d2 = 300,
                                 d3 = seq(300, 800, length.out=5),
                                 d4 = 800,
                                 period_blocks = 2, 
                                 mu0 = 0,
                                 sigma = 1,
                                 theta1 = 0.25,
                                 theta2 = 0.25,
                                 theta3 = 0.25,
                                 theta4 = 0.25,
                                 lambda0 = 0.15,
                                 lambda1 = 0.15,
                                 lambda2 = 0.15,
                                 lambda3 = 0.15,
                                 lambda4 = 0.15,
                                 trend = rep(c("linear", "stepwise_2", "inv_u"), each = 5),
                                 N_peak = c(rep(NA, 10), c(767, 774, 781, 790, 800)),
                                 alpha = 0.025,
                                 prec_theta = 0.001,
                                 prec_eta = 0.001,
                                 tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                 tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                 tau_case = "reasonable",
                                 prec_a = 0.001,
                                 prec_b = 0.001,
                                 bucket_size = 25,
                                 opt = 2, 
                                 prior_prec_tau = 0.002, 
                                 prior_prec_eta = 1,
                                 n_samples = 1000,
                                 n_chains = 4,
                                 n_iter = 4000,
                                 n_adapt = 1000,
                                 robustify = TRUE,
                                 weight = 0.1)

results_ii_eq_pow <- sim_study_par(nsim = n_sim, scenarios = scenario_ii_eq_pow, arms = c(4), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_ii_eq_pow, "results/results_ii_eq_pow.csv")


##############################################################################################################################################################

# SCENARIO III-I - 10 arms, equidistant entry times, different time trends in arms 1-9, inference on arm 10 ----

# sum(get_ss_matrix(num_arms = 10, n_arm = 250, d = c(300*(0:9))), na.rm=T)/2
# 
# plot_trial(treatments = datasim_bin(num_arms = 10, n_arm = 250, d =c(300*(0:9)), p0 = 0.7, OR = rep(1.8, 10), lambda = rep(0.15, 11), trend="stepwise")$treatment)

## TYPE I ERROR ----

set.seed(9)

scenario_iii_i_diff_alpha <- data.frame(num_arms = 10, 
                                        n_arm = 250,
                                        d1 = 300*0,
                                        d2 = 300*1,
                                        d3 = 300*2,
                                        d4 = 300*3,
                                        d5 = 300*4,
                                        d6 = 300*5,
                                        d7 = 300*6,
                                        d8 = 300*7,
                                        d9 = 300*8,
                                        d10 = 300*9,
                                        period_blocks = 2,
                                        mu0 = 0,
                                        sigma = 1,
                                        theta1 = 0,
                                        theta2 = 0,
                                        theta3 = 0,
                                        theta4 = 0,
                                        theta5 = 0,
                                        theta6 = 0,
                                        theta7 = 0,
                                        theta8 = 0,
                                        theta9 = 0,
                                        theta10 = 0,
                                        lambda0 = 0, 
                                        lambda1 = seq(-0.15, 0.15, length.out = 9),
                                        lambda2 = seq(-0.15, 0.15, length.out = 9),
                                        lambda3 = seq(-0.15, 0.15, length.out = 9),
                                        lambda4 = seq(-0.15, 0.15, length.out = 9),
                                        lambda5 = seq(-0.15, 0.15, length.out = 9),
                                        lambda6 = seq(-0.15, 0.15, length.out = 9),
                                        lambda7 = seq(-0.15, 0.15, length.out = 9),
                                        lambda8 = seq(-0.15, 0.15, length.out = 9),
                                        lambda9 = seq(-0.15, 0.15, length.out = 9),
                                        lambda10 = 0,
                                        trend = "linear",
                                        N_peak = NA,
                                        alpha = 0.025,
                                        prec_theta = 0.001,
                                        prec_eta = 0.001,
                                        tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                        tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                        tau_case = "reasonable",
                                        prec_a = 0.001,
                                        prec_b = 0.001,
                                        bucket_size = 25,
                                        opt = 2, 
                                        prior_prec_tau = 0.002, 
                                        prior_prec_eta = 1,
                                        n_samples = 1000,
                                        n_chains = 4,
                                        n_iter = 4000,
                                        n_adapt = 1000,
                                        robustify = TRUE,
                                        weight = 0.1)

results_iii_i_diff_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_iii_i_diff_alpha, arms = c(10), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_iii_i_diff_alpha, "results/results_iii_i_diff_alpha.csv")





## POWER ----

set.seed(10)

scenario_iii_i_diff_pow <- data.frame(num_arms = 10, 
                                      n_arm = 250,
                                      d1 = 300*0,
                                      d2 = 300*1,
                                      d3 = 300*2,
                                      d4 = 300*3,
                                      d5 = 300*4,
                                      d6 = 300*5,
                                      d7 = 300*6,
                                      d8 = 300*7,
                                      d9 = 300*8,
                                      d10 = 300*9,
                                      period_blocks = 2,
                                      mu0 = 0,
                                      sigma = 1,
                                      theta1 = 0.25,
                                      theta2 = 0.25,
                                      theta3 = 0.25,
                                      theta4 = 0.25,
                                      theta5 = 0.25,
                                      theta6 = 0.25,
                                      theta7 = 0.25,
                                      theta8 = 0.25,
                                      theta9 = 0.25,
                                      theta10 = 0.25,
                                      lambda0 = 0,
                                      lambda1 = seq(-0.15, 0.15, length.out = 9),
                                      lambda2 = seq(-0.15, 0.15, length.out = 9),
                                      lambda3 = seq(-0.15, 0.15, length.out = 9),
                                      lambda4 = seq(-0.15, 0.15, length.out = 9),
                                      lambda5 = seq(-0.15, 0.15, length.out = 9),
                                      lambda6 = seq(-0.15, 0.15, length.out = 9),
                                      lambda7 = seq(-0.15, 0.15, length.out = 9),
                                      lambda8 = seq(-0.15, 0.15, length.out = 9),
                                      lambda9 = seq(-0.15, 0.15, length.out = 9),
                                      lambda10 = 0,
                                      trend = "linear",
                                      N_peak = NA,
                                      alpha = 0.025,
                                      prec_theta = 0.001,
                                      prec_eta = 0.001,
                                      tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                      tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                      tau_case = "reasonable",
                                      prec_a = 0.001,
                                      prec_b = 0.001,
                                      bucket_size = 25,
                                      opt = 2, 
                                      prior_prec_tau = 0.002, 
                                      prior_prec_eta = 1,
                                      n_samples = 1000,
                                      n_chains = 4,
                                      n_iter = 4000,
                                      n_adapt = 1000,
                                      robustify = TRUE,
                                      weight = 0.1)


results_iii_i_diff_pow <- sim_study_par(nsim = n_sim, scenarios = scenario_iii_i_diff_pow, arms = c(10), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_iii_i_diff_pow, "results/results_iii_i_diff_pow.csv")





##############################################################################################################################################################


# SCENARIO III-II - 10 arms, random time trend in arms 1-9, inference on arm 10 ----


## TYPE I ERROR ----

set.seed(11)

scenario_iii_ii_rand_alpha <- data.frame(num_arms = 10, 
                                         n_arm = 250, 
                                         d1 = 300*0,
                                         d2 = 300*1,
                                         d3 = 300*2,
                                         d4 = 300*3,
                                         d5 = 300*4,
                                         d6 = 300*5,
                                         d7 = 300*6,
                                         d8 = 300*7,
                                         d9 = 300*8,
                                         d10 = 300*9,
                                         period_blocks = 2,
                                         mu0 = 0,
                                         sigma = 1,
                                         theta1 = 0,
                                         theta2 = 0,
                                         theta3 = 0,
                                         theta4 = 0,
                                         theta5 = 0,
                                         theta6 = 0,
                                         theta7 = 0,
                                         theta8 = 0,
                                         theta9 = 0,
                                         theta10 = 0,
                                         lambda0 = seq(-0.15, 0.15, length.out = 9), 
                                         lambda1 = "random",
                                         lambda2 = "random",
                                         lambda3 = "random",
                                         lambda4 = "random",
                                         lambda5 = "random",
                                         lambda6 = "random",
                                         lambda7 = "random",
                                         lambda8 = "random",
                                         lambda9 = "random",
                                         lambda10 = seq(-0.15, 0.15, length.out = 9),
                                         trend = "linear",
                                         N_peak = NA,
                                         trend_mean = seq(-0.15, 0.15, length.out = 9),
                                         trend_var = 0.5,
                                         alpha = 0.025,
                                         prec_theta = 0.001,
                                         prec_eta = 0.001,
                                         tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                         tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                         tau_case = "reasonable",
                                         prec_a = 0.001,
                                         prec_b = 0.001,
                                         bucket_size = 25,
                                         opt = 2, 
                                         prior_prec_tau = 0.002, 
                                         prior_prec_eta = 1,
                                         n_samples = 1000,
                                         n_chains = 4,
                                         n_iter = 4000,
                                         n_adapt = 1000,
                                         robustify = TRUE,
                                         weight = 0.1)

results_iii_ii_rand_alpha <- sim_study_par(nsim = n_sim, scenarios = scenario_iii_ii_rand_alpha, arms = c(10), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_iii_ii_rand_alpha, "results/results_iii_ii_rand_alpha.csv")


## POWER ----

# power.t.test(n=250, sd=1, sig.level = 0.025, power=0.8, type = "two.sample", alternative = "one.sided")

set.seed(12)

scenario_iii_ii_rand_pow <- data.frame(num_arms = 10, 
                                       n_arm = 250, 
                                       d1 = 300*0,
                                       d2 = 300*1,
                                       d3 = 300*2,
                                       d4 = 300*3,
                                       d5 = 300*4,
                                       d6 = 300*5,
                                       d7 = 300*6,
                                       d8 = 300*7,
                                       d9 = 300*8,
                                       d10 = 300*9,
                                       period_blocks = 2,
                                       mu0 = 0,
                                       sigma = 1,
                                       theta1 = 0.25,
                                       theta2 = 0.25,
                                       theta3 = 0.25,
                                       theta4 = 0.25,
                                       theta5 = 0.25,
                                       theta6 = 0.25,
                                       theta7 = 0.25,
                                       theta8 = 0.25,
                                       theta9 = 0.25,
                                       theta10 = 0.25,
                                       lambda0 = seq(-0.15, 0.15, length.out = 9), 
                                       lambda1 = "random",
                                       lambda2 = "random",
                                       lambda3 = "random",
                                       lambda4 = "random",
                                       lambda5 = "random",
                                       lambda6 = "random",
                                       lambda7 = "random",
                                       lambda8 = "random",
                                       lambda9 = "random",
                                       lambda10 = seq(-0.15, 0.15, length.out = 9),
                                       trend = "linear",
                                       N_peak = NA,
                                       trend_mean = seq(-0.15, 0.15, length.out = 9),
                                       trend_var = 0.5,
                                       alpha = 0.025,
                                       prec_theta = 0.001,
                                       prec_eta = 0.001,
                                       tau_a = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[1]),
                                       tau_b = unname(from_jump_size_to_tau_a_b(e_jump_size = 0.01, max_jump_size = 0.15)[2]),
                                       tau_case = "reasonable",
                                       prec_a = 0.001,
                                       prec_b = 0.001,
                                       bucket_size = 25,
                                       opt = 2, 
                                       prior_prec_tau = 0.002, 
                                       prior_prec_eta = 1,
                                       n_samples = 1000,
                                       n_chains = 4,
                                       n_iter = 4000,
                                       n_adapt = 1000,
                                       robustify = TRUE,
                                       weight = 0.1)

results_iii_ii_rand_pow <- sim_study_par(nsim = n_sim, scenarios = scenario_iii_ii_rand_pow, arms = c(10), models = c("fixmodel", "sepmodel", "poolmodel", "MAPpriorNew", "timemachine"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_iii_ii_rand_pow, "results/results_iii_ii_rand_pow.csv")


















