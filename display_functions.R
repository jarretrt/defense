## functions for displaying results

process_datasim_list <- function(datasim_list, npatients = 10){
  # subset and format for plotting
  data_obs <- subset(
    as.data.frame(do.call("rbind", 
                          lapply(datasim_list, `[[`, "sim"))), 
    select = c(id,time,pd0,pdobs)
  )
  data_obs$time <- round(data_obs$time,3)
  
  inf_obs <- subset(
    as.data.frame(do.call("rbind", 
                          lapply(datasim_list, `[[`, "inf"))), 
    select = c(id,infrt,begin,pdresp_start)
  )
  inf_obs$begin <- round(inf_obs$begin, 3)
  names(inf_obs)[names(inf_obs) == "pdresp_start"] <- "pdpred"
  
  dat_plot <- merge(data_obs, inf_obs, by.x = c("id","time"), by.y = c("id","begin"), 
                    all.y = TRUE)
  dat_plot <- dat_plot[,-which(names(dat_plot) == "pdobs")]
  dat_plot$plot_order <- rep(1:(nrow(dat_plot)/npatients),npatients)
  
  # transform into long format
  dat_plot_m <- reshape2::melt(dat_plot, id.vars = c("id","time","plot_order"))
  dat_plot_m$variable <- factor(dat_plot_m$variable, levels = c("pdpred","pd0","infrt"))
  levels(dat_plot_m$variable) <- c("Predicted response","Observed response","Infusion rate")
  id_order <- with(patient_covariates, order(M1F2, AGE, BMI))
  dat_plot_m$id <- factor(dat_plot_m$id, levels = id_order) 
  
  levels(dat_plot_m$id) <- paste0(ifelse(patient_covariates$M1F2[id_order] == 1, "Male","Female"), 
                                  ",\nAge=", round(patient_covariates$AGE[id_order]), 
                                  ",\nBMI=", round(patient_covariates$BMI[id_order]),"\n")
  
  dat_plot_m$target_val <- NA
  dat_plot_m$target_val[dat_plot_m$variable %in% c("Predicted repsonse", "Observed response")] <- 50
  
  return(dat_plot_m)
}



display_results_new2 <- function(sim_list, obj_fns, crits, titles, alphas, 
                                 return_plot_list = FALSE, nrow_facet = 1,
                                 legend_titles = NULL)
{
  
  tmax = 20
  bist = 50
  tm_eval = seq(0,tmax,0.01)
  
  set.seed(123736491)
  sample_size = 50
  in_sample <- sample(1:122,sample_size, replace = FALSE)
  out_sample <- setdiff(1:122, in_sample)
  patient_ix <- out_sample
  
  target_df_list <- vector("list", length(sim_list))
  pobj_list <- vector("list", length(sim_list))
  for(i in 1:length(sim_list)){
    sim <- sim_list[[i]]
    pbis <- sapply(sim[patient_ix], function(x) x$dat$bis_t)
    pobj <- sapply(sim[patient_ix], obj_fns[i], alpha = alphas[i])
    if(obj_fns[i] == "phi_stable_entry" & crits[i] == "sd"){
      # use root squared deviation from mean for illustration
      pobj <- sqrt((pobj - mean(pobj))^2)
    }
    pobj_list[[i]] <- pobj
    datp <- data.frame(id = as.factor(rep(patient_ix, each = nrow(pbis))),
                       time = sim[[1]]$dat$time,
                       objective = rep(pobj, each = nrow(pbis)),
                       bis = c(pbis))
    datp$type <- ifelse(datp$id %in% in_sample, "training", "test")
    target_pars <- sim[[1]]$beta_est
    
    if(length(target_pars) == 2){
      bis0 = target_pars[1]
      lam = target_pars[2]
      bis_targets <- (bis0-bist) * exp(-lam * tm_eval) + bist
      tf <- "Exponential"
    } else{
      t50 <- target_pars[1]
      Hill <- target_pars[2]
      bis0 <- target_pars[3]
      bis_targets <- bis0-(bis0-bist)*(tm_eval^Hill / (tm_eval^Hill + t50^Hill))
      tf <- "Sigmoidal"
    }
    
    target_df <- rbind(datp, data.frame(id = NA, time = tm_eval, 
                                        objective = NA,
                                        bis = bis_targets, type = "target"))
    target_df$simulation <- titles[i]
    target_df$objective_function <- obj_fns[i]
    target_df$target_function <- tf
    target_df_list[[i]] <- target_df
  }
  
  ## combine results
  df_all <- do.call("rbind",target_df_list)
  df_all$simulation <- factor(df_all$simulation, levels = titles)
  
  xq <- vector("list", length(unique(df_all$objective_function)))
  for(i in 1:length(unique(df_all$objective_function))){
    obji <- unique(df_all$objective_function)[i]
    tmp <- df_all[df_all$objective_function == obji,"objective"]
    
    if(i == 1){
      xq[[i]] <- quantile(unique(tmp), probs = c(seq(0,0.49,length.out = 40),
                                                 seq(0.50,0.94,length.out = 30),
                                                 seq(0.95,1,length.out = 30)),
                          na.rm = TRUE)
    } else{
      xq[[i]] <- quantile(unique(tmp), probs = c(seq(0,0.15,length.out = 30),
                                                 seq(0.16,0.84,length.out = 40),
                                                 seq(0.85,1,length.out = 30)),
                          na.rm = TRUE) 
    }
  }
  
  pop_crits <- sapply(1:length(pobj_list), function(i) do.call(crits[i], pobj_list[i]))
  labels <- data.frame(target_function=factor(titles, levels = titles),
                       label = paste("Population", Hmisc::capitalize(crits), 
                                     round(pop_crits,2)))
  
  df_all$id_objfn <- paste(df_all$objective_function, as.character(df_all$id), sep = "-")
  
  plot_list <- vector("list",3)
  if(is.null(legend_titles)){
    legend_titles <- c(expression("WOS"~alpha==0.05),"SET","MD") 
  }
  xmax = 20
  for(i in 1:length(unique(df_all$objective_function))){
    objfni <- unique(df_all$objective_function)[i]
    tmp1 <- df_all[df_all$objective_function == objfni,]
    tmp1$target_function <- as.factor(tmp1$target_function)
    levels(tmp1$target_function) <- list(titles[1:2], titles[3:4], titles[5:6])[[i]]
    rf_df <- tmp1[tmp1$type == "target",]
    names(rf_df)[names(rf_df) %in% c("time","bis")] <- c("rf_time","rf_bis")
    plot_list[[i]] <- ggplot(tmp1, aes(x = time, y = bis, color = objective)) + 
      facet_wrap(~target_function, labeller = label_parsed, nrow= nrow_facet) +
      geom_line(aes(group = id), alpha = 0.4) + 
      scale_color_gradientn(colors = pal100, 
                            values = rescale(xq[[i]], limits = range(xq[[i]]))) +
      geom_line(data = rf_df, aes(x = rf_time, y = rf_bis), color = "black", size = 1) +
      geom_hline(yintercept = c(40,60), linetype = "dashed") + 
      # theme(axis.title = element_blank(), legend.title.align = 0) +
      labs(color = legend_titles[i]) + 
      xlim(0,xmax) +
      # ylim(20,100) +
      theme(legend.justification = "left") +
      geom_text(data = labels[list(c(1,2),c(3,4),c(5,6))[[i]],], 
                aes(label=label),
                x = 14, y = 80, size = 4.5,
                inherit.aes = FALSE)
  }
  
  if(return_plot_list) return(plot_list)
  
  ## else combine plots
  library(gridExtra)
  library(cowplot)
  
  p <- plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]], 
                 nrow = 3, align = "hv", 
                 labels = "AUTO")
  
  arrangeGrob(p,left = textGrob("Bispectral Index (BIS)",
                                rot = 90, vjust = 1), 
              bottom = textGrob("Minutes"))
}
