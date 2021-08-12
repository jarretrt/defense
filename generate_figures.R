library(tci)
library(ggplot2)
library(gganimate)
source("display_functions.R")
theme_set(theme_bw())


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## generate plot of open-loop control

## example patient
npatients = 10
set.seed(1)
ids <- sample(1:nrow(eleveld_pk), npatients)
patient_covariates <- eleveld_pk[eleveld_pk$ID %in% ids, 
                                 c("ID","AGE","PMA","WGT","HGT","TECH","BMI","M1F2","A1V2")]
pkpd_vars <- c("V1","V2","V3","CL","Q2","Q3","KE0","CE50","BIS0","BIS0",
               "GAMMA","GAMMA2","SIGMA","BIS_DELAY")

theta <- eleveld_poppk(patient_covariates)[,pkpd_vars]
set.seed(1)
theta0 <- eleveld_poppk(patient_covariates, rand = TRUE)[,pkpd_vars]


## Calculate infusion schedule to reach BIS = 50
fixed_tci <- vector("list",npatients)
for(i in 1:npatients){
  fixed_tci[[i]] <- tci_pd(pdresp = c(50,50), 
                           tms = c(0,5), 
                           pdinv = inv_emax_eleveld, 
                           pdmod = emax_eleveld,
                           pkmod = pkmod3cptm,
                           pars_pk = unlist(theta[i,c("V1","V2","V3","CL","Q2","Q3","KE0")]),
                           pars_pd = unlist(theta[i,c("CE50","GAMMA","GAMMA2","BIS0","BIS0")]))
}


## generate observed data 
fixed_tci_datasim <- vector("list", npatients)
set.seed(1)
for(i in 1:npatients){
  # generate data based on "true" parameter values
  pars_pki0 <- unlist(theta0[i,c("V1","V2","V3","CL","Q2","Q3","KE0")])
  pars_pdi0 <- unlist(theta0[i,c("CE50","GAMMA","GAMMA2","BIS0","BIS0")])
  
  # fixed bis target
  fixed_tci_datasim[[i]] <- gen_data(inf = fixed_tci[[i]], 
                                     pkmod = pkmod3cptm, 
                                     pdmod = emax_eleveld, 
                                     pars_pk0 = pars_pki0, 
                                     pars_pd0 = pars_pdi0, 
                                     sigma_add = theta0[i,"SIGMA"], 
                                     delay = theta0[i,"BIS_DELAY"], 
                                     init = c(0,0,0,0))
  fixed_tci_datasim[[i]]$sim <- cbind(id = i, fixed_tci_datasim[[i]]$sim)
  fixed_tci_datasim[[i]]$inf <- cbind(id = i, fixed_tci_datasim[[i]]$inf)
}

dat_plot_m <- process_datasim_list(fixed_tci_datasim)

p <- ggplot(dat_plot_m, aes(x = time, y = value, color = id)) + 
  geom_hline(aes(yintercept = target_val)) +
  facet_wrap(~ variable, nrow = 3, scales = "free_y", 
             strip.position = "left",
             labeller = as_labeller(c(`Predicted response` = "Bispectral Index", 
                                      `Observed response` = "Bispectral Index",
                                      `Infusion rate` = "Infusion rate (mg/min)"))) +
  ylab(NULL) +
  geom_line() +
  scale_color_viridis_d() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(x = "Minutes")
  
p + transition_reveal(plot_order) 
anim_save(filename = "images/open_loop_sim.gif")

## ----------------------------------------------

rand_vars <- c("V1","V2","V3","CL","Q2","Q3","KE0","CE50","SIGMA")
set.seed(1)

## set up prior and true parameters for CLC
prior_par_list <- true_par_list <- vector("list", npatients)
for(i in 1:npatients){
  theta_samples <- replicate(1e3, unlist(eleveld_poppk(patient_covariates[i,], 
                                                       rand = TRUE)[,rand_vars]))
  ltheta_vcov <- cov(t(log(theta_samples))) + diag(diag(cov(t(log(theta_samples))))*0.01)
  prior_par_list[[i]] <- list(pars_pkpd = unlist(theta[i,1:12]), sig = ltheta_vcov,
                     pk_ix = 1:7, pd_ix = c(8,11,12,9,10), fixed_ix = 9:12,
                     err = theta[i,"SIGMA"], delay = theta[i,"BIS_DELAY"]/60)
  
  true_par_list[[i]] <- list(pars_pkpd = unlist(theta0[i,c("V1","V2","V3","CL","Q2","Q3","KE0",
                                                        "CE50","GAMMA","GAMMA2","BIS0","BIS0")]),
                             pk_ix = 1:7, pd_ix = c(8,11,12,9,10), fixed_ix = 9:12,
                             err = theta0[i,"SIGMA"],
                             delay = theta0[i,"BIS_DELAY"]/60) 
}

## targets and update times
targets <- data.frame(time = c(0,10), 
                      target = c(50,50))
update_times <- data.frame(time = seq(1,10,1), 
                           full_data = rep(TRUE,10))

set.seed(1)
bayes_sim <- vector("list", npatients)
for(i in 1:npatients){
  print(i)
  tmp <- bayes_control(targets = targets, 
                                  updates = update_times,
                                  prior = prior_par_list[[i]], 
                                  true_pars = true_par_list[[i]])
  tmp$dat$sim <- cbind(id = i, tmp$dat$sim)
  tmp$dat$inf <- cbind(id = i, tmp$dat$inf)
  bayes_sim[[i]] <- tmp
}

bayes_datsim_list <- lapply(bayes_sim, `[[`, "dat")
dat_plot_m_bis50 <- process_datasim_list(bayes_datsim_list)

# subset to observed & infusions and create initial value
dat_plot_m_bis50 <- dat_plot_m_bis50[dat_plot_m_bis50$variable %in% 
                                       c("Observed response", "Infusion rate"),]
dat_plot_m_bis50$value[is.na(dat_plot_m_bis50$value)] <- 93

pb50 <- ggplot(dat_plot_m_bis50, aes(x = time, y = value, color = factor(id))) + 
  geom_hline(aes(yintercept = target_val), size = 1.2) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", 
             strip.position = "left",
             labeller = as_labeller(c(`Predicted response` = "Bispectral Index", 
                                      `Observed response` = "Bispectral Index",
                                      `Infusion rate` = "Infusion rate (mg/min)"))) +
  ylab(NULL) +
  geom_line() +
  scale_color_viridis_d() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(x = "Minutes", title = "Fixed reference point at BIS = 50")

pb50 + transition_reveal(plot_order) 
anim_save(filename = "images/bis50_sim.gif")


## ----------------------------------------------

## reference function parameters and values
tfin <- 6 # time at which Bf+ep is targeted
Hill <- 3.5
B0 <- 75
ep <- 1
Bf <- 50 # final target + ep
t50 <- exp(log(tfin) + 1/Hill * (log(ep) - log(B0 - Bf - ep)))
tm_eval <- seq(0,10,1/6)
bis_targets <- B0-(B0-Bf)*(tm_eval^Hill / (tm_eval^Hill + t50^Hill))
plot(tm_eval, bis_targets, type ="l", ylim = c(40,100))

## targets and update times
targets <- data.frame(time = tm_eval, 
                      target = bis_targets)

bayes_sim_rf <- vector("list", npatients)
set.seed(1)
for(i in 1:npatients){
  print(i)
  tmp <- bayes_control(targets = targets, 
                       updates = update_times,
                       prior = prior_par_list[[i]], 
                       true_pars = true_par_list[[i]])
  tmp$dat$sim <- cbind(id = i, tmp$dat$sim)
  tmp$dat$inf <- cbind(id = i, tmp$dat$inf)
  bayes_sim_rf[[i]] <- tmp
}

bayes_datsim_rf_list <- lapply(bayes_sim_rf, `[[`, "dat")
dat_plot_m_rf <- process_datasim_list(bayes_datsim_rf_list)

# subset to observed & infusions and create initial value
dat_plot_m_rf <- dat_plot_m_rf[dat_plot_m_rf$variable %in% 
                                       c("Observed response", "Infusion rate"),]
dat_plot_m_rf$value[is.na(dat_plot_m_rf$value)] <- 93

targets$time <- round(targets$time,3)
dat_plot_m_rf <- merge(dat_plot_m_rf, targets, all.x = TRUE)
dat_plot_m_rf$target[is.na(dat_plot_m_rf$target_val)] <- NA
target_df <- dat_plot_m_rf[,c("time","target","variable")]
target_df <- target_df[!duplicated(target_df),]
names(target_df)[names(target_df) %in% c("time","target")] <- c("rf_time","rf_target")

prf <- ggplot(dat_plot_m_rf, aes(x = time, y = value, color = factor(id))) + 
  geom_hline(aes(yintercept = target_val)) +
  geom_line(data = target_df, aes(x = rf_time, y = rf_target), color = "black", size = 1.2) + 
  facet_wrap(~ variable, nrow = 2, scales = "free_y", 
             strip.position = "left",
             labeller = as_labeller(c(`Predicted response` = "Bispectral Index", 
                                      `Observed response` = "Bispectral Index",
                                      `Infusion rate` = "Infusion rate (mg/min)"))) +
  ylab(NULL) +
  geom_line() +
  scale_color_viridis_d() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(x = "Minutes", title = "Sigmoidal reference function")

prf + transition_reveal(plot_order) 
anim_save(filename = "Presentation/images/rf_sim.gif")


ggsave(filename = "images/rf_static.jpeg", prf)
# ptitle <- paste0("Sigmoid reference function: BIS0=",B0,
#                  ", t50=",round(t50,2),
#                  ", Hill=",Hill)


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## Show distribution of points for BIS 50 and sigmoidal 
library(dplyr)
td_df <- dat_plot_m_rf %>% 
  group_by(id) %>% 
  filter(variable == "Infusion rate") %>%
  summarise(total_dose = sum(value)*1/6)

f <- function(){
  x <- sample(td_df$total_dose,1)
  x + rnorm(1, sd = x/20)
}
td_vec <- replicate(1000, f())
tdq <- quantile(td_vec, c(0.25,0.5,0.75))
dose_hist <- ggplot(data.frame(x = td_vec), aes(x = x)) + 
  geom_histogram(bins = 30, color = "grey") + 
  xlab("Total dosage propofol (mg)") +
  ylab("Number of patients") +
  ggtitle("Distribution of doses at BIS0 = 75, Hill = 3.5, t50 = 2.4") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  geom_vline(xintercept = mean(td_vec), color = "blue") +
  geom_text(x = mean(td_vec) + 40, y = 150, label = "Mean = 167", color = "blue") +
  geom_vline(xintercept = max(td_vec), color = "red") +
  geom_text(x = max(td_vec) - 40, y = 150, label = "Max = 433", color = "red") +
  annotate("segment", x = tdq[1], xend = tdq[3], y = 130, yend = 130, color = "green")
# +
  # geom_text(x = max(td_vec) - 40, y = 150, label = "Max = 433", color = "green")


dose_hist
ggsave(filename = "images/dose_hist.jpeg", dose_hist)


Hill <- 3.5
B0 <- 75
ep <- 1
Bf <- 50 # final target + ep
t50 <- exp(log(tfin) + 1/Hill * (log(ep) - log(B0 - Bf - ep)))



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## illustration of closed-loop control
theta <- c(v1 = 5, v2 = 17.297, v3 = 90.963, cl = 1.382, 
           q2 = 0.919, q3 = 0.609, ke0 = 1.289, c50 = 2.8, gamma = 1.47, 
           gamma2 = 1.89, e0 = 93, emx = 93, sigma = 8.03, bis_delay = 28.263)
theta0 <- c(v1 = 9, v2 = 19.857, v3 = 63.423, cl = 2.142, q2 = 1.237, 
            q3 = 0.258, ke0 = 2.088, c50 = 3.026, gamma = 1.47, gamma2 = 1.89,
            e0 = 93, emx = 93, 
            # sigma = 6.936, 
            sigma = 4,
            bis_delay = 28.263)

target = 50
inf1 <- tci_pd(pdresp = rep(target,2), tms = c(0,10), pdinv = inv_emax, 
                  pdmod = emax, pkmod = pkmod3cptm, pars_pk = theta[1:7], 
                  pars_pd = theta[8:12])

sim1 <- gen_data(inf = inf1, pkmod = pkmod3cptm, pdmod = emax,
                     pars_pk0 = theta0, pars_pd0 = theta0,
                     sigma_add = theta0["sigma"], delay = theta0["bis_delay"], 
                 init = c(0,0,0,0))

init_vals <- sim1$sim[nrow(sim1$sim), c("c1","c2","c3","c4")]
inf2 <- tci_pd(pdresp = rep(target,2), tms = c(10,20), pdinv = inv_emax, 
               pdmod = emax, pkmod = pkmod3cptm, pars_pk = theta0[1:7], 
               pars_pd = theta0[8:12], init = init_vals)

sim2 <- gen_data(inf = inf2, pkmod = pkmod3cptm, pdmod = emax,
                 pars_pk0 = theta0, pars_pd0 = theta0,
                 sigma_add = theta0["sigma"], delay = theta0["bis_delay"],
                 init = init_vals)

sim_all <- combine_sim(sim1, sim2)
inf_all <- as.data.frame(sim_all$inf[,c("infrt","begin")])
inf_all$begin <- round(inf_all$begin,3)

dat_clc <- as.data.frame(sim_all$sim[,c("time","pd0","pdobs")])
dat_clc$time <- round(dat_clc$time,3)

dat_clcm <- merge(dat_clc, inf_all, by.x = "time", by.y = "begin")
dat_clcm <- reshape2::melt(dat_clcm, id.vars = c("time","pdobs"))
dat_clcm$pdobs[dat_clcm$variable == "infrt"] <- NA
dat_clcm$target <- ifelse(dat_clcm$variable == "infrt", NA, 50)

pclc <- ggplot(dat_clcm, aes(x = time, y = value)) + 
  facet_wrap(~ variable, nrow = 2, scales = "free_y", 
             strip.position = "left",
             labeller = as_labeller(c(`pd0` = "Bispectral Index", 
                                      `infrt` = "Infusion rate (mg/min)"))) +
  geom_line(size = 1.2) +
  geom_point(aes(x= time, y = pdobs, group = seq_along(pdobs)), 
             color = "blue", alpha = 0.4) +
  labs(x = "Minutes", y = "Bispectral Index (BIS)", 
       title = "Targeting BIS=50, update at 10 min") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_hline(aes(yintercept = target))

pclc + transition_reveal(time)
anim_save(filename = "images/clc_sim.gif")



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## Results 
load("sims/bis50_opt_res_N122.rda")
load("sims/sig_opt_res_N122_rand.rda")

## stable entry plot
sim_list <- list(bis50_opt_res_N122, sig_opt_res_N122_rand[[3]])
crits <- c("mean","mean")
obj_fns <- c("phi_stable_entry","phi_stable_entry")
alphas <- c(0.5,0.5)
titles <- c("BIS50", "RFs-AvgSET")
p_set <- display_results_new2(sim_list, obj_fns, crits, titles, alphas, 
                          return_plot_list = TRUE, nrow_facet = 2, legend_titles = "SET")[[1]]
p_set <- p_set + xlab("Minutes") + ylab("Bispectral Index")
ggsave("images/results_set.jpeg", p_set, device = "jpeg")
# p_set + transition_reveal(time)


## total dose / SD
sim_list <- list(bis50_opt_res_N122, sig_opt_res_N122_rand[[5]])
obj_fns <- c("phi_min_dose","phi_min_dose")
titles <- c("BIS50", "RFs-AvgTD")
p_td <- display_results_new2(sim_list, obj_fns, crits, titles, alphas, 
                              return_plot_list = TRUE, nrow_facet = 2, legend_titles = "TD")[[1]]
p_td <- p_td + xlab("Minutes") + ylab("Bispectral Index")
ggsave("images/results_td.jpeg", p_td, device = "jpeg")




## kaplan meier curves
library(prodlim)
library(rootSolve)
library(numDeriv)

## in/out of sample
set.seed(123736491)
sample_size = 50
in_sample <- sample(1:122,sample_size, replace = FALSE)
out_sample <- setdiff(1:122, in_sample)

simres_list <- list(bis50_opt_res_N122[out_sample], 
                    sig_opt_res_N122_rand[[3]][out_sample])

group_names <- c("BIS50","RFs-SET")
entry_times <- data.frame(
  time = c(sapply(simres_list, function(x) 
    sapply(x, phi_stable_entry))),
  group = rep(group_names, each = length(out_sample))
)
entry_times$cens <- 0

km0 <- prodlim(Hist(time,cens)~group, 
               data=entry_times, 
               reverse = TRUE)

rise_times <- data.frame(
  time = c(sapply(simres_list, function(x) 
    sapply(x, phi_rise_time))),
  group = rep(group_names, each = length(out_sample))
)

rise_times$cens <- 0
km1 <- prodlim(Hist(time,cens)~group, 
               data=rise_times, 
               reverse = TRUE)
library(colorspace)
dcols <- darken(adjustcolor(pal100[c(55,95,40,10,80)], alpha.f = 1), amount = 0.3)


dev.off()
pdf("images/kaplan-meier.pdf",
    width = 10, height = 5)
par(mfrow = c(1,2))

plot(km0,  
     xlab = "Minutes",
     ylab = "Percent of patients in BIS = [40,60]", 
     col = dcols[2:3],
     atrisk = FALSE,
     legend = TRUE,
     legend.x = "bottomright",
     legend.title = "",
     axis1.cex.axis = 1.2,
     axis2.cex.axis = 1.2,
     plot.cex.lab = 1.2,
     plot.main = "Stable Entry")

plot(km1,  
     xlab = "Minutes",
     ylab = "Percent of patients with BIS < 60", 
     col = dcols[2:3],
     atrisk = FALSE,
     legend = TRUE,
     legend.x = "bottomright",
     legend.title = "",
     axis1.cex.axis = 1.2,
     axis2.cex.axis = 1.2,
     plot.cex.lab = 1.2,
     plot.main = "Rise Time")
dev.off()


pdf("images/dose-difference.pdf",
    width = 8, height = 8)
## change in dosages
BIS50_TD <- sapply(bis50_opt_res_N122[out_sample], phi_min_dose)
TFS_TD <- sapply(sig_opt_res_N122_rand[[5]][out_sample], phi_min_dose)
tdd_list <- list(BIS50_TD = BIS50_TD, TFS_TD =TFS_TD)
tdd <-  TFS_TD-BIS50_TD
breaks <- pretty(tdd,20)
col <- ifelse(breaks >= 0, dcols[1], dcols[4])

hist(tdd, xlab = "Difference in mg propofol (RFs-TD - BIS50)",
     main = "Difference in total dose among test set patients",
     breaks = breaks,
     col = col, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0, lwd = 5)
text(x = -120, y = 7, labels = paste0("Lower dose\nusing RFs-TD\n\nN=",sum(tdd<0)), cex = 1.5)
text(x = 40, y = 7, labels = paste0("Higher dose\nusing RFs-TD\n\nN=",sum(tdd>=0)), cex = 1.5)
dev.off()






## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## Results 
## separate testing set
set.seed(123736491)
sample_size = 50
in_sample <- sample(1:122,sample_size, replace = FALSE)
out_sample <- setdiff(1:122, in_sample)
eleveld_pkpd <- merge(eleveld_pk, eleveld_pd)

covariates_test <- eleveld_pkpd[out_sample, 
                               c("ID","AGE","PMA","WGT","HGT","TECH","BMI","M1F2","A1V2")]
pkpd_vars <- c("V1","V2","V3","CL","Q2","Q3","KE0","CE50","BIS0","BIS0",
               "GAMMA","GAMMA2","SIGMA","BIS_DELAY")


## test set parameters - predicted and true
## Note: PD variable names are different than output of eleveld_poppk()
theta_test <- eleveld_poppk(covariates_test)[,pkpd_vars]
theta0_test <- eleveld_pkpd[out_sample,c("V1","V2","V3","CL","Q2","Q3","KE0",
                                          "E50","EMAX","EMAX","GAM","GAM1","RESD","ALAG1")]
theta0_test$ALAG1 <- theta0_test$ALAG1*60
names(theta0_test) <- names(theta_test)


## Set up prior for patients
rand_vars <- c("V1","V2","V3","CL","Q2","Q3","KE0","CE50","SIGMA")
set.seed(1)

## set up prior and true parameters for CLC
npatients <- nrow(theta_test)
prior_par_list <- true_par_list <- vector("list", npatients)
for(i in 1:npatients){
  theta_samples <- replicate(1e2, unlist(eleveld_poppk(covariates_test[i,], 
                                                       rand = TRUE)[,rand_vars]))
  ltheta_vcov <- cov(t(log(theta_samples))) + diag(diag(cov(t(log(theta_samples))))*0.01)
  prior_par_list[[i]] <- list(pars_pkpd = unlist(theta_test[i,1:12]), sig = ltheta_vcov,
                              pk_ix = 1:7, pd_ix = c(8,11,12,9,10), fixed_ix = 9:12,
                              err = theta_test[i,"SIGMA"], delay = theta_test[i,"BIS_DELAY"]/60)
  
  true_par_list[[i]] <- list(pars_pkpd = unlist(theta0_test[i,c("V1","V2","V3","CL","Q2","Q3","KE0",
                                                           "CE50","GAMMA","GAMMA2","BIS0","BIS0")]),
                             pk_ix = 1:7, pd_ix = c(8,11,12,9,10), fixed_ix = 9:12,
                             err = theta0_test[i,"SIGMA"],
                             delay = theta0_test[i,"BIS_DELAY"]/60) 
}

## targets and update times
bis50_targets <- data.frame(time = c(0,20), 
                      target = c(50,50))
update_times <- data.frame(time = c(1,2,4,8,10,12,15,20), 
                           full_data = c(rep(TRUE,3),rep(FALSE,5)))


## simulate responses for BIS50
set.seed(1)
bayes_sim_bis50_test <- vector("list", npatients)
for(i in 1:npatients){
  print(i)
  tmp <- bayes_control(targets = bis50_targets, 
                       updates = update_times,
                       prior = prior_par_list[[i]], 
                       true_pars = true_par_list[[i]])
  tmp$dat$sim <- cbind(id = i, tmp$dat$sim)
  tmp$dat$inf <- cbind(id = i, tmp$dat$inf)
  bayes_sim_bis50_test[[i]] <- tmp
}

bayes_sim_bis50_test <- lapply(bayes_sim, `[[`, "dat")
dat_plot_m_bis50_test <- process_datasim_list(bayes_datsim_list)
setwd("./Presentation/")
save(dat_plot_m_bis50_test, file = "./sims/dat_plot_m_bis50_test.rda")








## Sigmoid time in range
Lam_tir <- c(55,5.7,5) # (BIS0,t50,gamma)
tm_eval <- seq(0,20,1/6)





## reference function parameters and values
tfin <- 6 # time at which Bf+ep is targeted
Hill <- 3.5
B0 <- 75
ep <- 1
Bf <- 50 # final target + ep
t50 <- exp(log(tfin) + 1/Hill * (log(ep) - log(B0 - Bf - ep)))
tm_eval <- seq(0,10,1/6)
bis_targets <- B0-(B0-Bf)*(tm_eval^Hill / (tm_eval^Hill + t50^Hill))
plot(tm_eval, bis_targets, type ="l", ylim = c(40,100))

## targets and update times
targets <- data.frame(time = tm_eval, 
                      target = bis_targets)

bayes_sim_rf <- vector("list", npatients)
set.seed(1)
for(i in 1:npatients){
  print(i)
  tmp <- bayes_control(targets = targets, 
                       updates = update_times,
                       prior = prior_par_list[[i]], 
                       true_pars = true_par_list[[i]])
  tmp$dat$sim <- cbind(id = i, tmp$dat$sim)
  tmp$dat$inf <- cbind(id = i, tmp$dat$inf)
  bayes_sim_rf[[i]] <- tmp
}

bayes_datsim_rf_list <- lapply(bayes_sim_rf, `[[`, "dat")
dat_plot_m_rf <- process_datasim_list(bayes_datsim_rf_list)

# subset to observed & infusions and create initial value
dat_plot_m_rf <- dat_plot_m_rf[dat_plot_m_rf$variable %in% 
                                 c("Observed response", "Infusion rate"),]
dat_plot_m_rf$value[is.na(dat_plot_m_rf$value)] <- 93

targets$time <- round(targets$time,3)
dat_plot_m_rf <- merge(dat_plot_m_rf, targets, all.x = TRUE)
dat_plot_m_rf$target[is.na(dat_plot_m_rf$target_val)] <- NA
target_df <- dat_plot_m_rf[,c("time","target","variable")]
target_df <- target_df[!duplicated(target_df),]
names(target_df)[names(target_df) %in% c("time","target")] <- c("rf_time","rf_target")

prf <- ggplot(dat_plot_m_rf, aes(x = time, y = value, color = factor(id))) + 
  geom_hline(aes(yintercept = target_val)) +
  geom_line(data = target_df, aes(x = rf_time, y = rf_target), color = "black", size = 1.2) + 
  facet_wrap(~ variable, nrow = 2, scales = "free_y", 
             strip.position = "left",
             labeller = as_labeller(c(`Predicted response` = "Bispectral Index", 
                                      `Observed response` = "Bispectral Index",
                                      `Infusion rate` = "Infusion rate (mg/min)"))) +
  ylab(NULL) +
  geom_line() +
  scale_color_viridis_d() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(x = "Minutes")

prf + transition_reveal(plot_order) 
anim_save(filename = "Presentation/images/rf_sim.gif")




### display region uncertainty
library(mvtnorm)
library(ggplot2)
theme_set(theme_bw())
N = 50
group <- rep(c("Group 1","Group 2"), each = N/2)
mu1 = c(1,3)
mu2 = c(1.5,-3)
X <- rbind(rmvnorm(N/2, mu1), rmvnorm(N/2, mu2))
df <- as.data.frame(X)
df$group <- group
names(df) <- c("X1","X2","group")
ggplot(df, aes(x = X1, y = X2, color = group)) + geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  scale_color_hue(l=50, c=80) +
  theme(legend.position = "top") + 
  geom_vline(xintercept = 0.9, linetype = "dotted", size = 1) +
  annotate("text", x = 1.2, y = 0.4, label = "(0.9,NA)") +
  geom_rug(outside = TRUE)
  


