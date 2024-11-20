
# Cleaning data -------------------------------------------------

clean_data <- function(data, mice_filter = T) {
  cleaned_data <- data %>% 
    
    # Turn plasma total folate concentrations into a categorical variable
    dplyr::mutate(pl_fol3.t1 = case_when(tot_fol.t1 < quantile(tot_fol.t1, 0.1, na.rm = T) ~ 1,
                                         tot_fol.t1 > quantile(tot_fol.t1, 0.8, na.rm = T) ~ 3,
                                         is.na(tot_fol.t1) ~ NA,
                                         .default = 2)) %>% # Between 10th and 80th %ile is assigned 2
    
    # Not using this variable. Forgot to remove it.
    dplyr::mutate(pl_fol3.t3 = case_when(tot_fol.t3 < quantile(tot_fol.t3, 0.1, na.rm = T) ~ 1,
                                         tot_fol.t3 > quantile(tot_fol.t3, 0.8, na.rm = T) ~ 3,
                                         is.na(tot_fol.t3) ~ NA,
                                         .default = 2)) %>% # Between 10th and 80th %ile is assigned 2
    
    dplyr::mutate(umfa2.t1 = ifelse(umfa.t1.res > quantile(umfa.t1.res, 0.8, na.rm = T),
                                                            1, # ppl w/ umfa conc. above 80th percentile assigned 1
                                                            0)) %>%  # otherwise they are assigned 2
    
    # Turn year of enrollment into a categorical variable
    dplyr::mutate(year.enroll4 = year.enroll4 - 2007) %>% 
    
    # Fix the city variable, only cities where SRS was assessed will be included. 
    mutate(city6 = case_when(city10 == 1 ~ 1,
                             city10 == 4 ~ 2,
                             city10 == 5 ~ 3,
                             city10 == 7 ~ 4,
                             city10 == 9 ~ 5,
                             city10 == 10 ~ 6)) %>% 
    
    # Define reference level for categorical variables 
    dplyr::mutate(fol.intake3 =  factor(fol.intake3, levels = c(2, 1, 3))) %>% # Adequate FA is reference
    dplyr::mutate(pl_fol3.t1 =   factor(pl_fol3.t1, levels = c(2, 1, 3))) %>% # 10th - 80th %ile is reference
    dplyr::mutate(pl_fol3.t3 =   factor(pl_fol3.t3, levels = c(2, 1, 3))) %>% # 10th - 80th %ile is reference
    
    dplyr::mutate(city6 =        factor(city6, levels = c(5, 1, 2, 3, 4, 6))) %>% # Montreal is reference
    dplyr::mutate(income4 =      factor(income4, levels = c(4,1,2,3))) %>%  # > $100K/ year is reference
    dplyr::mutate(edu4 =         factor(edu4, levels = c(3,1,2,4))) %>%  # Undergraduate is reference
    dplyr::mutate(parity3 =      factor(parity3, levels = c(1,2,3))) %>% # nulliparous is reference
    dplyr::mutate(year.enroll4 = factor(year.enroll4, levels = c(3, 1, 2, 4))) %>%  # 2010 is reference
    
    dplyr::mutate(fish3 =        factor(fish3, levels = c(1,2,3))) %>%  # 0-2 times per month is reference
    
    dplyr::mutate(mom.birthplace5 = factor(mom.birthplace5, levels = c(1, 2, 3, 4, 5))) %>%  # Born in Canada is reference
    
    
    # Make mean = 0 for cnts confounding variables
    dplyr::mutate(home.score = home.score - mean(data$home.score, na.rm = T)) %>% 
    dplyr::mutate(srs.age = srs.age - mean(data$srs.age, na.rm = T)) %>% 
    dplyr::mutate(mom.age = mom.age - mean(data$mom.age))
    
    if(mice_filter == T){
      # When I run the analysis, the data will only include variables I need and variables that are necessary for the
      # analysis and variables that will aid with imputation
      
      # But to check that everything is working correctly, I need to look at all variables. To do this quickly, 
      # I will set mice_filter to FALSE, then reset it back to true when I run the analysis. 
      
      # Select variables used in the analysis
      cleaned_data <- cleaned_data %>%
        dplyr::select(c(srs, # outcome
                      arsenic.t1.flag:lead.t1.res, mercury.t1.flag:mercury.t1.res, # metals
                      bbhc.t1.flag:transnona.t1.res, # OC pesticides
                      pfhxs.t1.flag:pfoa.t1.res, # PFAS
                      pcb118.t1.flag:pcb180.t1.res, # PCBs
                      bde47.t1.flag:bde47.t1.res, # PBDE
                      fol.intake3, pl_fol3.t1, sex2, umfa2.t1, # modifiers
                      income4, edu4, living.status2, home.score, race.white2, # confounders
                      mom.age, parity3, year.enroll4, srs.age, city6, smoker2,
                      fish3,
                      
                      # Variables that aid with imputation 
                      fol.intake, tot_fol.t3, umfa.t3.res,
                      birth.wt, gest.age, alc2, prepreg.bmi, mom.birthplace5)) 
    } 

  return(cleaned_data)
}

pvalue_rounder <- function(p) 
{
  require(dplyr)
  formatted_p <- case_when(
    p > 0.99 ~ "p-int.>.99", # for P values greater than .99, report as "P>.99."
    p <= 0.99 & p > 0.05 ~ paste0("p-int.=", (format(round(p, 2), nsmall = 2))), # or P values greater than or equal to .01, report the value to the nearest hundredth
    p <= 0.05 & p >= 0.01 ~ paste0("p-int.=", format(round(p, 2), nsmall = 2), "*"),
    p < 0.01 & p >= 0.001 ~ paste0("p-int.=", format(round(p, 3), nsmall = 3), "*"), # for P values between .001 and .01, report the value to the nearest thousandth;
    p < 0.001 ~ "p.int<.001*" # for P values less than .001, report as "P<.001"
  )
  
  return(formatted_p)
}

pvalue_rounder_v2 <- function(p) 
{
  require(dplyr)
  formatted_p <- case_when(
    p > 0.99 ~ "p > .99", # for P values greater than .99, report as "P>.99."
    p <= 0.99 & p > 0.05 ~ paste0("p = ", (format(round(p, 2), nsmall = 2))), # or P values greater than or equal to .01, report the value to the nearest hundredth
    p <= 0.05 & p >= 0.01 ~ paste0("p = ", format(round(p, 2), nsmall = 2), "*"),
    p < 0.01 & p >= 0.001 ~ paste0("p = ", format(round(p, 3), nsmall = 3), "*"), # for P values between .001 and .01, report the value to the nearest thousandth;
    p < 0.001 ~ "p < .001*" # for P values less than .001, report as "P<.001"
  )
  
  return(formatted_p)
}


# quantile g-comp ---------------------------------------------
get_qgcomp_weights <- function(fit, expnms = NULL, mixture_name = NULL, 
                               mod_var_name = NULL, imputation_number = NULL) {
  require(dplyr)
  
  if(is.null(expnms) == T ){ # if no exposure names are provided (default)
    weights <- c(fit$pos.weights, -1 * fit$neg.weights) %>% 
      as.data.frame() %>% 
      rename(weight = ".") %>% 
      mutate(name = row.names(.)) %>% 
      select(c(name, weight)) %>% 
      arrange(name) ## arrange weights alphabetically 
    
  } else { # if exposure names (expnms) are provided
    weights <- c(fit$pos.weights, -1 * fit$neg.weights) %>% 
      as.data.frame() %>% 
      rename(weight = ".") %>% 
      mutate(name = row.names(.)) %>% 
      select(c(name, weight)) %>%     
      arrange(match(name, expnms)) ## arrange weights in the order given by "expnms"
  }
  
  # Add the name of the mixture, if provided
  if(is.null(mixture_name) == F){
    weights <- weights %>%
      mutate(mixture_name = mixture_name) %>%
      select(c(mixture_name, name, weight))
  }
  
  # If effect modification is being assessed, then the modifying variable, its level will also be included 
  if(class(fit)[1] == "qgcompemmweights"){
    weights <- weights %>% 
      mutate(mod_var_name = mod_var_name) %>% 
      mutate(mod_level = fit$emmval)
  }
  
   # Add the imputation number, if provided
  if(is.null(imputation_number) == F){
    weights <- weights %>% 
      mutate(imputation_number = imputation_number)
  }
  
  return(weights)
}


analysis_qgcomp <- function(chemical_names_list, mixture_names, y, confounder_names = c(), 
                            mice_data,
                            q = 4, round = 2, round_p = 3, alpha = 0.05){
  require(qgcomp)
  
  # Prep ----
  table_qgcomp <- setNames(data.frame(matrix(data = NA, nrow = length(chemical_names_list), 
                                             ncol = 6)), #make dataframe
                           c("mixture_name", "psi", "lb", "ub", "se", "p")) #set column names
  all_weights <- list()
  all_pooled_weights <- list()
  
  # Analysis ----
  m <- mice_data$m # Store number of imputations
  
  for (i in 1:length(chemical_names_list)) {
    this_mixture_name <- mixture_names[i]
    this_chemical_names <- chemical_names_list[[i]]
    
    this_formula <- as.formula(paste(y, " ~ ", 
                                     paste(c(this_chemical_names, confounder_names), collapse = "+")))
    
    # Loop through 'm' imputed datasets ----
    
    results_all_k_datasets <- list() # I need these 2 lists to reset for every new i-loop
    weights_all_k_datasets <- list()
    
    for (k in 1:m){
      # Extract the kth imputed dataset
      imputed_data_k <- complete(imputed_data, action = k)
      
      # Fit the model with the kth imputed dataset
      fit <- qgcomp(this_formula, q = q, expnms = this_chemical_names, 
                    data = imputed_data_k)
      
      # Save the fitted model
      results_all_k_datasets[[k]] <- fit
      
      # Save the weights
      weights_all_k_datasets[[k]] <- get_qgcomp_weights(fit, expnms = this_chemical_names, 
                                                        mixture_name = this_mixture_name,
                                                        imputation_number = k)
    }
    
    # Pool results ----
    pooled_results <- pool(results_all_k_datasets)
    qgcomp_res <- summary(pooled_results)
    
    # Store pooled results ----
    psi <- qgcomp_res[2, 2]
    se <- qgcomp_res[2, 3]
    df <- qgcomp_res[2, 5]
    p <- qgcomp_res[2, 6]
    
    t <- qt(1 - (alpha / 2), df) 
    lb <- psi - (t * se)
    ub <- psi + (t * se)
    
    
    table_qgcomp[i, 1] <- this_mixture_name   # Store mixture name
    table_qgcomp[i, 2] <- round(psi, round)   # psi
    table_qgcomp[i, 3] <- round(lb, round)    # lb
    table_qgcomp[i, 4] <- round(ub, round)    # ub
    table_qgcomp[i, 5] <- round(se, round)    # se
    table_qgcomp[i, 6] <- round(p, round_p)   # p value
    
    # Pool weights ----
    df_weights_list <- bind_rows(weights_all_k_datasets)
    
    pooled_weights <- df_weights_list %>% 
      group_by(name, mixture_name) %>% 
      summarize(ave_weight = mean(weight)) %>% 
      arrange(match(name, this_chemical_names)) # Put chemicals in the desired order
    
    # Store weights ----
    all_weights[[i]] <- df_weights_list
    all_pooled_weights[[i]] <- pooled_weights
  }
  
  df_all_weights <- as.data.frame(do.call(rbind, all_weights))
  df_all_pooled_weights <- as.data.frame(do.call(rbind, all_pooled_weights))
  
  # Store all results as a list 
  return(list(table_qgcomp, df_all_weights, df_all_pooled_weights, 
              all_weights, all_pooled_weights))
}


analysis_qgcomp.emm_binary <- function(chemical_names_list, mixture_names, y, emmvar, 
                                       confounder_names = c(), mice_data, alpha = 0.05,
                                       mod_var_name = NULL, level_names = c("1", "2"), 
                                       q = 4, round = 2, round_p = 3){
  
  # Prep ----
  if(is.null(mod_var_name) == T){ # If the character string 'mod_var_name' is not provided, 
    mod_var_name == emmvar # Then use what was written for the 'emmvar' arguement
  }
  
  # Make table to store results
  table_qgcomp_mod <- setNames(data.frame(matrix(data = NA, nrow = length(chemical_names_list) * 2, 
                                                 ncol = 7)),
                               c("mixture_name", "psi", "lb", "ub", 
                                 "p_int", "mod_var", "mod_level"))
  all_weights <- list()
  all_pooled_weights <- list()
  
  # Analysis ----
  m <- mice_data$m # Store number of imputations
  
  for (i in 1:length(chemical_names_list)) {
    
    # Select exposures for this iteration of the loop
    this_mixture_name <- mixture_names[i]
    this_chemical_names <- chemical_names_list[[i]]
    this_formula <- as.formula(paste(y, " ~ ", 
                                     paste(c(this_chemical_names, confounder_names), collapse = "+")))
    
    # Make lists to store data needed to calculate results
    results_all_k_datasets <- list()
    vcov_all_k_datasets <- list()
    weights_all_k_datasets <- list()
    
    # Loop through 'm' imputed datasets----
    for (k in 1:m){
      
      # Extract the kth imputed dataset
      imputed_data_k <- complete(imputed_data, action = k)
      
      # Add a new column for the modifying variable so the next part works
      modifying_variable <- imputed_data_k[,emmvar]
      
      imputed_data_k <- imputed_data_k %>% 
        mutate(modifying_variable = modifying_variable) 
      ## The qgcomp.emm.noboot() function requires that a character string be 
      ## passed to the 'emmvar' argument. So I had to come up with the above 
      ## work-around to allow this custom function from to work
      ## with any modifying variable. 
      
      # Fit the model with the kth imputed dataset
      fit_mice <- qgcomp.emm.noboot(this_formula, q = q, expnms = this_chemical_names, 
                                    data = imputed_data_k,
                                    emmvar = "modifying_variable")
      
      # Extract model, variance-covarian matrix for all k datasets
      results_all_k_datasets[[k]] <- fit_mice
      vcov_all_k_datasets[[k]] <- fit_mice$covmat.coef
      
      # Extract weights for all k datasets 
      weights_all_k_datasets[[k]] <- lapply(0:1, function(x) 
        # do the task below when emmval equals 0 and 1
        get_qgcomp_weights(getstratweights(fit_mice, emmval = x), 
                           expnms = this_chemical_names,
                           mixture_name = this_mixture_name, 
                           mod_var_name = mod_var_name, 
                           imputation_number = k) )
    }
    
    # Pool results----
    mice_pooled_results <- pool(results_all_k_datasets)
    summary <- summary(mice_pooled_results)
    
    # Extract model parameters
    psi_x <- summary$estimate[2]
    psi_xz <- rev(summary$estimate)[1] # selects last entry
    p_int <- rev(summary$p.value)[1] # selects last entry
    
    # calculate variance-covariance matrix (Using Rubin's rules)
    ubar <- Reduce("+", vcov_all_k_datasets) / (m)
    b <- mice_pooled_results$pooled$b 
    vcov <- ubar + (1 + 1 / (m)) * b
    
    # store values needed for calculating intervals
    var_x <- diag(vcov)[2]
    var_xz <- rev(diag(vcov))[1] # selects last entry
    cov_x_xz <- rev(vcov[,2])[1] # selects the last entry of the 2nd column
    
    df <- summary$df[2]
    t <- qt(1 - (alpha / 2), df)
    
    # Calculate psi and 95% intervals
    psi_lvl0 <- psi_x
    psi_lvl1 <- psi_x + psi_xz
    
    lb_lvl0 <- psi_lvl0 - ( t * sqrt(var_x ) )
    ub_lvl0 <- psi_lvl0 + ( t * sqrt(var_x ) )
    
    lb_lvl1 <- psi_lvl1 - ( t * sqrt(var_x + var_xz + 2 * cov_x_xz) )
    ub_lvl1 <- psi_lvl1 + ( t * sqrt(var_x + var_xz + 2 * cov_x_xz) )
    
    # Store pooled results ----
    ## Information for all levels
    table_qgcomp_mod[(1 + (i-1) * 2): (2 + (i-1) * 2), 1] <- this_mixture_name  # Mixture name
    table_qgcomp_mod[(1 + (i-1) * 2): (2 + (i-1) * 2), 6] <- mod_var_name       # modifying variable name
    
    ## level 0
    table_qgcomp_mod[1 + (i-1) * 2, 2] <- round(psi_lvl0, round)  # psi
    table_qgcomp_mod[1 + (i-1) * 2, 3] <- round(lb_lvl0, round)   # lb
    table_qgcomp_mod[1 + (i-1) * 2, 4] <- round(ub_lvl0, round)   # ub
    table_qgcomp_mod[1 + (i-1) * 2, 5] <- NA                      # p-interaction
    table_qgcomp_mod[1 + (i-1) * 2, 7] <- level_names[1]          # modifying variable level
    
    ## level 1
    table_qgcomp_mod[2 + (i-1) * 2, 2] <- round(psi_lvl1, round)  # psi
    table_qgcomp_mod[2 + (i-1) * 2, 3] <- round(lb_lvl1, round)   # lb
    table_qgcomp_mod[2 + (i-1) * 2, 4] <- round(ub_lvl1, round)   # ub
    table_qgcomp_mod[2 + (i-1) * 2, 5] <- round(p_int, round_p)   # p-interaction
    table_qgcomp_mod[2 + (i-1) * 2, 7] <- level_names[2]          # modifying variable level
    
    
    # Pool weights ----
    df_weights_all_k_datasets<- bind_rows(weights_all_k_datasets)
    
    pooled_weights <- df_weights_all_k_datasets %>%
      group_by(name, mixture_name, mod_var_name, mod_level) %>%
      summarize(ave_weight = mean(weight)) %>%
      arrange(match(name, this_chemical_names)) # Put chemicals in the desired order
    
    # Store weights ----
    # Store weights for each chemical mixture into a separate list
    all_weights[[i]] <- df_weights_all_k_datasets
    all_pooled_weights[[i]] <- pooled_weights
  }
  
  # Output 
  # Turn the lists made above into a dataframe (easier to interpret)
  df_all_weights <- as.data.frame(do.call(rbind, all_weights))
  df_all_pooled_weights <- as.data.frame(do.call(rbind, all_pooled_weights))
  
  # print("The output is a list with 5 items: \n1) A dataframe with the regression results
  #       \n2) A dataframe with all of the weights for each imputation
  #       \n3) A dataframe with the pooled weights
  #       \n4) A list for each chemical mixture with all of the weights for each imputation
  #       \n5) A list for each chemical mixture with the pooled weights")
  return(list(table_qgcomp_mod, df_all_weights, df_all_pooled_weights, 
              all_weights, all_pooled_weights))
}


analysis_qgcomp.emm_cat3 <- function(chemical_names_list, mixture_names, y, emmvar, 
                                     confounder_names = c(), mice_data, alpha = 0.05,
                                     mod_var_name = NULL, level_names = c("1", "2", "3"), 
                                     q = 4, round = 2, round_p = 3){
  
  # Prep ----
  if(is.null(mod_var_name) == T){ # If the character string 'mod_var_name' is not provided, 
    mod_var_name == emmvar # Then use what was written for the 'emmvar' arguement
  }
  
  # emmvar <- as.character(emmvar)
  
  # Make table to store results
  table_qgcomp_mod <- setNames(data.frame(matrix(data = NA, nrow = length(chemical_names_list) * 3, 
                                                 ncol = 7)),
                               c("mixture_name", "psi", "lb", "ub", 
                                 "p_int", "mod_var", "mod_level"))
  all_weights <- list()
  all_pooled_weights <- list()
  
  # Analysis ----
  m <- mice_data$m # Store number of imputations
  
  for (i in 1:length(chemical_names_list)) {
    
    # Select exposures for this iteration of the loop
    this_mixture_name <- mixture_names[i]
    this_chemical_names <- chemical_names_list[[i]]
    this_formula <- as.formula(paste(y, " ~ ", 
                                     paste(c(this_chemical_names, confounder_names), collapse = "+")))
    
    # Make lists to store data needed to calculate results
    results_all_k_datasets <- list()
    vcov_all_k_datasets <- list()
    weights_all_k_datasets <- list()
    
    # Loop through 'm' imputed datasets----
    for (k in 1:m){
      
      # Extract the kth imputed dataset
      imputed_data_k <- complete(imputed_data, action = k)
      
      # Add a new column for the modifying variable so the next part works
      modifying_variable <- imputed_data_k[,emmvar]
      
      imputed_data_k <- imputed_data_k %>% 
        mutate(modifying_variable = modifying_variable) 
      ## The qgcomp.emm.noboot() function requires that a character string be 
      ## passed to the 'emmvar' argument. So I had to come up with the above 
      ## work-around to allow this custom function from to work
      ## with any modifying variable. 
      
      # Fit the model with the kth imputed dataset
      fit_mice <- qgcomp.emm.noboot(this_formula, q = q, expnms = this_chemical_names, 
                                    data = imputed_data_k,
                                    emmvar = "modifying_variable")
      
      # Extract model, variance-covarian matrix for all k datasets
      results_all_k_datasets[[k]] <- fit_mice
      vcov_all_k_datasets[[k]] <- fit_mice$covmat.coef
      
      # Extract weights for all k datasets 
      weights_all_k_datasets[[k]] <- lapply(1:3, function(x) 
        # do the task below when emmval equals 1, 2, and 3
        get_qgcomp_weights(getstratweights(fit_mice, emmval = x), 
                           expnms = this_chemical_names,
                           mixture_name = this_mixture_name, 
                           mod_var_name = mod_var_name, 
                           imputation_number = k) )
    }
    
    # Pool results----
    mice_pooled_results <- pool(results_all_k_datasets)
    summary <- summary(mice_pooled_results)
    
    # Extract model parameters
    psi_x <- summary$estimate[2]
    psi_xz2 <- summary$estimate[4] 
    psi_xz3 <- summary$estimate[6] 
    p_int_xz2 <- summary$p.value[4] 
    p_int_xz3 <- summary$p.value[6] 
    
    # calculate variance-covariance matrix (Using Rubin's rules)
    ubar <- Reduce("+", vcov_all_k_datasets) / (m)
    b <- mice_pooled_results$pooled$b 
    vcov <- ubar + (1 + 1 / (m)) * b
    
    # store values needed for calculating intervals
    var_x <- diag(vcov)[2] 
    
    var_xz2 <- diag(vcov)[4] 
    var_xz3 <- diag(vcov)[6] 
    
    cov_x_xz2 <- vcov[4,2] 
    cov_x_xz3 <- vcov[6,2]
    
    df <- summary$df[2]
    t <- qt(1 - (alpha / 2), df)
    
    # Calculate psi and 95% intervals
    psi_lvl1 <- psi_x  # just equals psi_x when lvl1 = 0 (default)
    psi_lvl2 <- psi_x + psi_xz2
    psi_lvl3 <- psi_x + psi_xz3
    
    lb_lvl1 <- psi_lvl1 - ( t * sqrt(var_x))
    ub_lvl1 <- psi_lvl1 + ( t * sqrt(var_x))
    
    lb_lvl2 <- psi_lvl2 - ( t * sqrt(var_x + var_xz2 + 2*cov_x_xz2) )
    ub_lvl2 <- psi_lvl2 + ( t * sqrt(var_x + var_xz2 + 2*cov_x_xz2) )
    
    lb_lvl3 <- psi_lvl3 - ( t * sqrt(var_x + var_xz3 + 2*cov_x_xz3) )
    ub_lvl3 <- psi_lvl3 + ( t * sqrt(var_x + var_xz3 + 2*cov_x_xz3) )
    
    # Store pooled results ----
    ## Information for all levels
    table_qgcomp_mod[(1 + (i-1) * 3): (3 + (i-1) * 3), 1] <- this_mixture_name  # Mixture name
    table_qgcomp_mod[(1 + (i-1) * 3): (3 + (i-1) * 3), 6] <- mod_var_name       # modifying variable name
    
    ## level 1
    table_qgcomp_mod[1 + (i-1) * 3, 2] <- round(psi_lvl1, round)  # psi
    table_qgcomp_mod[1 + (i-1) * 3, 3] <- round(lb_lvl1, round)   # lb
    table_qgcomp_mod[1 + (i-1) * 3, 4] <- round(ub_lvl1, round)   # ub
    table_qgcomp_mod[1 + (i-1) * 3, 5] <- NA                      # p-interaction
    table_qgcomp_mod[1 + (i-1) * 3, 7] <- level_names[1]          # modifying variable level
    
    ## level 2
    table_qgcomp_mod[2 + (i-1) * 3, 2] <- round(psi_lvl2, round)    # psi
    table_qgcomp_mod[2 + (i-1) * 3, 3] <- round(lb_lvl2, round)     # lb
    table_qgcomp_mod[2 + (i-1) * 3, 4] <- round(ub_lvl2, round)     # ub
    table_qgcomp_mod[2 + (i-1) * 3, 5] <- round(p_int_xz2, round_p) # p-interaction
    table_qgcomp_mod[2 + (i-1) * 3, 7] <- level_names[2]            # modifying variable level
    
    ## level 3
    table_qgcomp_mod[3 + (i-1) * 3, 2] <- round(psi_lvl3, round)    # psi
    table_qgcomp_mod[3 + (i-1) * 3, 3] <- round(lb_lvl3, round)     # lb
    table_qgcomp_mod[3 + (i-1) * 3, 4] <- round(ub_lvl3, round)     # ub
    table_qgcomp_mod[3 + (i-1) * 3, 5] <- round(p_int_xz3, round_p) # p-interaction
    table_qgcomp_mod[3 + (i-1) * 3, 7] <- level_names[3]            # modifying variable level
    
    
    # Pool weights ----
    df_weights_all_k_datasets<- bind_rows(weights_all_k_datasets)
    
    pooled_weights <- df_weights_all_k_datasets %>%
      group_by(name, mixture_name, mod_var_name, mod_level) %>%
      summarize(ave_weight = mean(weight)) %>%
      arrange(match(name, this_chemical_names)) # Put chemicals in the desired order
    
    # Store weights ----
    # Store weights for each chemical mixture into a separate list
    all_weights[[i]] <- df_weights_all_k_datasets
    all_pooled_weights[[i]] <- pooled_weights
  }
  
  # Output
  # Turn the lists made above into a dataframe (easier to interpret)
  df_all_weights <- as.data.frame(do.call(rbind, all_weights))
  df_all_pooled_weights <- as.data.frame(do.call(rbind, all_pooled_weights))
  
  # print("The output is a list with 5 items: \n1) A dataframe with the regression results
  #       \n2) A dataframe with all of the weights for each imputation
  #       \n3) A dataframe with the pooled weights
  #       \n4) A list for each chemical mixture with all of the weights for each imputation
  #       \n5) A list for each chemical mixture with the pooled weights")
  return(list(table_qgcomp_mod, df_all_weights, df_all_pooled_weights, 
              all_weights, all_pooled_weights))
}

# Saving -----------------------------------------------------------------------
my_write.csv <- function(x, path, file_name, archive = T,
                         archive_path = "archive/", ...) {
  
  main_pathway <- paste0(path, file_name, ".csv", sep = "")
  
  # Save file
  write.csv(x, row.names = F, file = main_pathway, ...)
  
  # Save file in archive
  if(archive == T){
    
    archive_pathway <- paste0(path, archive_path, file_name, "_",  Sys.Date(), ".csv", sep = "")
    
    # save file
    write.csv(x, row.names = F, file = archive_pathway, ...)
    
  }
  
}


my_ggsave <- function(plot = NULL, path, file_name,  archive = T,
                      archive_path = "archive/", device = "png",
                      width = 7.5, height = 5, dpi =300, ...) {
  
  file_type <- paste0(".", device)
  main_pathway <- paste0(path, file_name, file_type, sep = "")
  
  # Save file
  ifelse(is.null(plot) == T,
         ggsave(file = main_pathway, width = width, height = height, dpi = dpi, device = device, ...),
         ggsave(plot, file = main_pathway, width = width, height = height, dpi = dpi, device = device, ...))
  # If no plot is provided, then do not include this argument and ggsave() will default to the last plot. 
  # Otherwise, the provided plot will be saved
  
  # Save file in archive
  if(archive == T){
    
    archive_pathway <- paste0(path, archive_path, file_name, "_",  Sys.Date(), file_type, sep = "")
    
    # save file
    ifelse(is.null(plot) == T,
           ggsave(file = archive_pathway, width = width, height = height, dpi = dpi, device = device, ...),
           ggsave(plot, file = archive_pathway, width = width, height = height, dpi = dpi, device = device, ...))
    
  }
  
}
