#functions to estimate applicant scores from partial reviews
library(tidyverse)
library(ibd)

########## coin flip test ###########
run_coin_flip_test <- function(score_frame){
  #function to run coin flip test of whether or not reviewers exhibit differential harshness
  #currently set up for 2 reviews per applicant
  #the null is that (under random assignment) reviewers have probability 0.5 of reviewing a candidate more harshly then their paired reviewer
  #input: 
    #score_frame: a dataframe of observed data, with column for reviewer, applicant, and score
  #output:
    #dataframe with reviewers and columns for summaries, test statistics, and p-value. 
  
  # Add column that indicates if reviewer score was top score for applicant
  data_top_score <- score_frame %>% 
    group_by(applicant) %>%
    mutate(top_score = ifelse(score == max(score), 1, 0))
  # Sum the number of top scores each reviewer gave and find the proportion of top score
  # and p-value from binomial test
  data_coin_flip <- data_top_score %>% group_by(reviewer) %>%
    summarise(num_top_score = sum(top_score), count = n(), prop_top_score = num_top_score/count,
              p.value = binom.test(num_top_score, count, p = 0.5, alternative = "two.sided")$p.value)
  # Return data with column of pvalues testing whether each reviewer is biased
  return(data_coin_flip)
}


########## assign missing values to worst case #########
#auxiliary function for assign_NA
get_objective <- function(imputed_frame){
  #given a frame return value of objective function describing difference between top scores between Males and Non-Males. 
  #this is the objective to be minimized by assign_NA
  #inputs:
    #imputed_frame: a complete dataframe (with missing values imputed)
  #outputs:
    #the imbalance between proportions of top scores within each group, assessed by absolute difference
  proportions_gender <- tapply(imputed_frame$top_score, imputed_frame$applicant_gender, mean)
  abs(proportions_gender["Male"] - proportions_gender['NonMale'])
}

#auxiliary function for assign_NA()
impute_NAs <- function(missing_frame, NA_vals){
  #given a data.frame with missing values and a vector of values to fill them, fill them.
  #input:
    #missing_frame: a dataframe with missing applicant_gender values
    #NA_vals: a vector of values to fill in the missing ones in the df
  #output:
    #a dataframe with no missing values
  filled_frame <- missing_frame %>%
    mutate(applicant_gender = ifelse(is.na(applicant_gender), NA_vals, applicant_gender))
  filled_frame
}


assign_NA <- function(score_sub_frame){
  #function to assign NA values such that 
  #input:
    #score_sub_frame: a score_frame for a single reviewer with columns applicant, applicant_gender, reviewer, reviewer_sex, score, top_score. applicant_gender may be missing some values
  #output:
    #a score frame where applicant gender has been assigned so as to minimize the imbalance between top_score proportions
  
  num_NA <- sum(is.na(score_sub_frame$applicant_gender))
  if(num_NA > 0){
      #objective to be minimized is the absolute difference between top score propotions
      
      #generate permutations of "Male" and "Non-Male"
      possible_imputations <- ifelse(gtools::permutations(n = 2, r = num_NA, repeats.allowed = TRUE) == 1, "Male", "NonMale")
      
      possible_imputed_frames <- apply(possible_imputations, 1, impute_NAs, missing_frame = score_sub_frame)
    
    
      objectives <- possible_imputed_frames %>%
        map(get_objective) %>%
        unlist()
      minimum <- which.min(objectives)
      imputed_frame <- impute_NAs(missing_frame = score_sub_frame, NA_vals = possible_imputations[minimum,])
      imputed_frame
  } else if(num_NA == 0){
      imputed_frame <- score_sub_frame
      imputed_frame
  }
}


########### gender difference test #########
run_gender_test <- function(score_frame, na_handling = "worst-case"){
  #function to test whether reviewers exhibit bias when scoring applicants
  #the null is that the proportion of higher reviews for men and non-men is the same for each reviewer
  
  data_top_score <- score_frame %>% 
    group_by(applicant) %>%
    mutate(top_score = ifelse(score == max(score), 1, 0)) %>%
    ungroup()

  if(na_handling == "worst-case"){
    complete_top_score <- data_top_score %>%
      split(data_top_score$reviewer) %>%
      map(assign_NA) %>%
      reduce(bind_rows)
  } else if(na_handling == "omit"){
    complete_top_score <- data_top_score %>%
      na.omit()
  }
  
  #check to ensure that all reviewers reviewed at least 1 applicant of each gender category
  distinct_gender_reviews <- complete_top_score %>%
    group_by(reviewer) %>%
    summarize(n_gender_reviews = n_distinct(applicant_gender))
  #stop if a reviewer reviewed only one gender
  if(any(distinct_gender_reviews$n_gender_reviews < 2)){
    stop(paste("One or more reviewers (", paste(distinct_gender_reviews$reviewer[distinct_gender_reviews$n_gender_reviews < 2], collapse = ","), ") reviewed only one gender of applicant. Drop these reviewers."))
  }
  
  reviewers_top_score <- complete_top_score %>%
    group_by(reviewer, applicant_gender) %>% 
    summarize(prop_top_score = mean(top_score), num_top_score = sum(top_score), num_lower_score = sum(top_score == 0)) %>% 
    pivot_wider(names_from = "applicant_gender", values_from = c("prop_top_score","num_top_score","num_lower_score")) %>%
    mutate(p_value = fisher.test(x = rbind(c(num_lower_score_Male, num_lower_score_NonMale), c(num_top_score_Male, num_top_score_NonMale)), alternative = "two.sided")$p.value)
  
  reviewers_top_score
}


########### estimate mean #############
estimate_observed_mean <- function(score_frame){
  #input: 
    #score_frame: a dataframe of observed data, with column for reviewer, applicant, and score
  #output:
    #a dataframe with a column for applicant ids and a column for their observed means
  score_frame %>%
    group_by(applicant) %>%
    summarize(observed_mean = mean(score))
}

estimate_weighted_mean <- function(score_frame){
  I <- n_distinct(score_frame$applicant)
  J <- n_distinct(score_frame$reviewer)
  reviews_per_applicant <- unique(table(score_frame$applicant))
  
  score_frame %>% 
    group_by(reviewer) %>%
    mutate(workload = n()) %>%
    ungroup() %>%
    group_by(applicant) %>% 
    summarize(weighted_mean = sum(I / (workload * J) * score)) 
}

######### estimate ANOVA ###########
estimate_anova <- function(score_frame){
  #input: 
    #score_frame: a dataframe of observed data, with column for reviewer, applicant, and score
  #output:
    #a dataframe with column for applicant ids and a column for their score as predidcted using ANOVA
  #make a dataframe where every reviewer scores every applicant, to generate predicted scores
  complete_frame <- expand.grid("reviewer" = unique(score_frame$reviewer), "applicant" = unique(score_frame$applicant))
  fitted_model <- lm(score ~ reviewer + applicant - 1, data = score_frame)
  score_predictions <- predict(fitted_model, newdata = complete_frame)
  ANOVA_scores <- complete_frame %>% 
    mutate(score = score_predictions) %>%
    group_by(applicant) %>%
    summarize(ANOVA_score = mean(score))
  ANOVA_scores
}

######## estimate LSE ##########
estimate_LSE <- function(score_frame){
  #input:
    #score_frame: a dataframe of observed data, with column for reviewer, applicant, and score
  #output:
    #a dataframe with a column for applicants and a column for scores estimated using Longford's shrinkage estimator (LSE)
  #the three sum-of-squares
  #within applicants
  S_E <- score_frame %>% 
    group_by(applicant) %>%
    summarize(squared_diff = (first(score) - last(score))^2) %>%
    ungroup() %>%
    summarize(SSQ = (1/2) * sum(squared_diff)) %>%
    pull(SSQ)
  
  #within reviewers
  S_R <- score_frame %>%
    group_by(reviewer) %>% 
    mutate(reviewer_mean = mean(score)) %>%
    ungroup() %>%
    mutate(squared_diff = (score - reviewer_mean)^2) %>%
    summarize(SSQ = sum(squared_diff)) %>%
    pull(SSQ)
  
  #total 
  S_T <- score_frame %>% 
    mutate(squared_diff = (score - mean(score))^2) %>%
    summarize(SSQ = sum(squared_diff)) %>%
    pull(SSQ)
  
  #Longford uses the method of moments, matching the expectations of these statistics with the parameter values to be estimated. 
  #this allows us to estimate variance terms using the three SSQs above
  
  I <- length(unique(score_frame$applicant))
  J <- length(unique(score_frame$reviewer))
  n2 <- score_frame %>%
    group_by(reviewer) %>%
    summarize(squared_workload = n()^2) %>%
    ungroup() %>%
    summarize(normalized_SSQ_workload = sum(squared_workload) / (2*I)) %>%
    pull(normalized_SSQ_workload)
  
  #see table 4 of longford for estimates from that study
  #variance of reviewer parameters
  sigma_b_sq <- (S_T - ((2*I - 1)/(2*I - J)) * S_R) / (2*I - n2)
  #variance of error (inconsitency)
  #sigma_e can be NA, because S_E/I - sigma_b^2 may be less than 0, in this case set it equal to 0
  sigma_e_sq <- S_E/I - sigma_b_sq
  #variance of applicant parameters (alpha)
  sigma_a_sq <- S_R / (2*I - J) - sigma_e_sq
  
  #variance ratios
  tau_b <- sigma_b_sq / sigma_a_sq
  tau_e <- sigma_e_sq / sigma_a_sq
  
  #correlations
  #within applicant correlation
  r_a <- 1 / (1 + tau_b + tau_e)
  #correlation of mean of observed score with true score 
  r_b <- 1 / sqrt(1 + (tau_b + tau_e) / 2)
  
  #adjustments
  applicant_means <- score_frame %>%
    group_by(reviewer) %>%
    mutate(reviewer_mean = mean(score)) %>%
    ungroup() %>%
    group_by(applicant) %>%
    summarize(applicant_mean = mean(score), reviewer_applicant_mean = mean(reviewer_mean))
  
  overall_mean <- mean(score_frame$score)
  
  #mean workload of reviewers who grade each essay
  n_plus <- score_frame %>%
    group_by(reviewer) %>%
    mutate(workload = n()) %>%
    ungroup() %>%
    group_by(applicant) %>%
    summarize(n_plus = mean(workload))
  
  #mean reciprocal workload
  n_minus <- score_frame %>%
    group_by(reviewer) %>%
    mutate(reciprocal_workload = 1/n()) %>%
    ungroup() %>%
    group_by(applicant) %>%
    summarize(n_minus = mean(reciprocal_workload))
  
  #adjust using LSE, with approximately optimal adjustment factor u_i^star
  LSE_score_frame <- n_plus %>%
    left_join(n_minus, by = "applicant") %>%
    left_join(applicant_means, by = "applicant") %>%
    #mutate(u_star = (sigma_b^2 * (I - n_plus) + sigma_e^2 * (I * n_minus - 1) / 2) / (sigma_a^2 * (I * n_minus - 2) + sigma_b^2 * (I - 2 * n_plus + n2) + sigma_e^2 * (I * n_minus -1))) %>%
    mutate(u_star = (sigma_b_sq * (I - n_plus) + sigma_e_sq * (I * n_minus - 1) / 2) / (sigma_a_sq * (I * n_minus - 2) + sigma_b_sq * (I - 2 * n_plus + n2) + sigma_e_sq * (I * n_minus -1))) %>%
    mutate(LSE_score = applicant_mean - u_star * (reviewer_applicant_mean - overall_mean)) %>%
    dplyr::select(applicant, LSE_score)
  LSE_score_frame
}



############## estimate all scores ##############
estimate_all_scores <- function(score_frame){
  #inputs:
    #score_frame: a dataframe of observed data, with column for reviewer, applicant, and score
  #outputs:
    #a dataframe with column for applicant ids and columns for scores estimated from observed means, ANOVA, and LSE
  observed_scores <- estimate_observed_mean(score_frame = score_frame)
  ANOVA_scores <- estimate_anova(score_frame = score_frame)
  LSE_scores <- estimate_LSE(score_frame = score_frame)
  all_scores <- observed_scores %>%
    left_join(ANOVA_scores, by = "applicant") %>%
    left_join(LSE_scores, by = "applicant")
  all_scores
}


############### compute ranks ##############
compute_ranks <- function(estimated_scores_frame){
  #inputs: 
    #estimated_scores_frame: dataframe of applicants and scores as output by estimate_all_scores
  #outputs:
    #a dataframe with column for applicant ID and columns for ranks according to each score
  rank_reverse <- function(x){rank(-x, ties.method = "first")}
  ranks_frame <- estimated_scores_frame %>%
    mutate_at(2:4, rank_reverse)
  ranks_frame
}

######## pairs plot results ##########
plot_results_pairs <- function(estimated_scores){
  estimated_scores <- estimated_scores %>%
    rename("Observed Mean" = observed_mean, "ANOVA" = ANOVA_score, "LSE" = LSE_scores)
  pairs(~ `Observed Mean` + ANOVA + LSE, data = estimated_scores, pch = 20)
}


########## run a single simulation #############
run_simulation <- function(n_applicants = 100, n_reviewers = 10, model_violation = "none"){
  #given some design parameters, run a simulation of estimating scores with subsetting
  #inputs:
    #n_applicants: the number of applicants
    #n_reviewers: the number of reviewers
    #model_violation: none is the default (model is correct); "constant bias" adds a differential bias depending on applicant covariates; "poisson noise" adds additional poisson noise to the observed estimates
  #output:
    #a dataframe with each estimate (observed mean, ANOVA, and LSE) and the true overall score (s_i)
  applicant_id <- 1:n_applicants
  reviewer_id <- 1:n_reviewers
 
  #fixed true applicant score (alpha)
  true_scores <- data.frame(applicant = as_factor(applicant_id), true_score = rep(24:42, length.out = n_applicants))
  #fixed true reviewer harshness (beta)
  true_reviewer_harshness <- data.frame(reviewer = as_factor(reviewer_id), harshness = rep(-3:3, length.out = n_reviewers))
  
  #complete scores
  complete_scores <- expand.grid(applicant = true_scores$applicant, reviewer = true_reviewer_harshness$reviewer) %>%
    left_join(true_scores, by = "applicant") %>%
    left_join(true_reviewer_harshness, by = "reviewer") %>%
    mutate(score = true_score + harshness)
  
  if(model_violation == "constant bias"){
    #add covariate that reviewers are biased on
    applicant_cov_bias <- data.frame(applicant = as_factor(applicant_id), covariate = sample(c("A", "B"), n_applicants, replace = TRUE)) %>%
      mutate(extra_bias = ifelse(covariate == "A", -1, 0))

    complete_scores <- complete_scores %>%
      left_join(applicant_cov_bias, by = "applicant") %>%
      mutate(score = score + extra_bias)
  } else if(model_violation == "minor bias"){
    #applicant_cov_bias <- data.frame(applicant = as_factor(applicant_id), app_covariate = rep(c("A","B"), each = n_applicants / 2))
    #reviewer_cov_bias <- data.frame(reviewer = as_factor(reviewer_id), rev_covariate = rep(c("A","B"), each = n_reviewers / 2))
    applicant_cov_bias <- data.frame(applicant = as_factor(applicant_id), app_covariate = sample(c("A", "B"), n_applicants, replace = TRUE))
    reviewer_cov_bias <- data.frame(reviewer = as_factor(reviewer_id), rev_covariate = sample(c("A", "B"), n_reviewers, replace = TRUE))
    # add -1 to score if applicant and reviewer both possess covariate value "A"
    complete_scores <- complete_scores %>%
      left_join(applicant_cov_bias, by = "applicant") %>%
      left_join(reviewer_cov_bias, by = "reviewer") %>%
      mutate(score = ifelse((app_covariate == "A" & rev_covariate == "A"), score - 1, score))
  } else if(model_violation == "poisson noise"){
    # add poisson error, lambda = 2
    complete_scores <- complete_scores %>% 
      mutate(score = score + rpois(nrow(sim_scores), 2))
  } else if(model_violation == "major bias"){
    #applicant_cov_bias <- data.frame(applicant = as_factor(applicant_id), app_covariate = rep(c("A","B"), each = n_applicants / 2))
    #reviewer_cov_bias <- data.frame(reviewer = as_factor(reviewer_id), rev_covariate = rep(c("A","B"), each = n_reviewers / 2))
    applicant_cov_bias <- data.frame(applicant = as_factor(applicant_id), app_covariate = sample(c("A", "B"), n_applicants, replace = TRUE))
    reviewer_cov_bias <- data.frame(reviewer = as_factor(reviewer_id), rev_covariate = sample(c("A", "B"), n_reviewers, replace = TRUE))
    # add -1 to score if applicant and reviewer both possess covariate value "A"
    complete_scores <- complete_scores %>%
      left_join(applicant_cov_bias, by = "applicant") %>%
      left_join(reviewer_cov_bias, by = "reviewer") %>%
      mutate(score = ifelse((app_covariate == "A" & rev_covariate == "A"), score - 3, score))
  }
  
  #simulate subsetting by random assignment
  sim_scores <- complete_scores %>%
    group_by(applicant) %>%
    sample_n(size = 2, replace = FALSE) %>%
    ungroup()
  
  
  #compute complete row means (s_i), the estimand
  true_panel_scores <- complete_scores %>%
    group_by(applicant) %>%
    summarize(true_panel_score = mean(score))
  
  estimates <- estimate_all_scores(sim_scores) %>%
    left_join(true_panel_scores, by = "applicant")
  estimates
}

########### compute summaries simulation #######
compute_summaries <- function(sim_estimates){
  #input: 
    #sim_estimates: a dataframe of applicants, true scores, and estimates, as output by run_simulation
  #output:
    #the mean squared error of the score estimates (computed across applicants)
  estimates_long <- sim_estimates %>%
    pivot_longer(cols = c("observed_mean", "ANOVA_score", "LSE_score"), names_to = "estimator", values_to = "estimate")
  summaries <- estimates_long %>%
    group_by(estimator) %>%
    mutate(squared_error = (true_panel_score - estimate)^2) %>%
    summarize(MSE = mean(squared_error))
  summaries
}

########### run a batch of simulations ###########
run_simulations_batch <- function(n_sims = 1000, n_applicants = 100, n_reviewers = 10){
  simulations <- replicate(n_sims, run_simulation(n_applicants = n_applicants, n_reviewers = n_reviewers, model_violation = "poisson noise"), simplify = FALSE)
  results <- simulations %>%
    map(compute_summaries) %>%
    reduce(bind_rows) %>%
    group_by(estimator) %>%
    summarize(avg_MSE = mean(MSE))
  results
}



############ check for existence of a BIBD #########
check_for_BIBD <- function(number_applicants, number_reviewers, number_reviews_per_applicant = 2){
  #function to check for a BIBD given number of assignments, returning true (if one exists) or false
  number_applicants_per_reviewer <- (number_reviews_per_applicant * number_applicants) / number_reviewers
  lambda <- number_reviews_per_applicant * (number_applicants_per_reviewer - 1) / (number_applicants - 1)
  if(number_applicants_per_reviewer %% 1 != 0){
    #"No BIBD exists. Number of needed reviews is not a multiple of the number of reviewers --> can't assign equal workload to all reviewers."
    FALSE
  } else if(lambda %% 1 != 0){
    #"No BIBD exists. There is no assignment that groups applicants within reviewers an equal number of times (lambda is not an integer). Change the size of the panel (number of reviewers), the number of reviews per applicant, or consider an alternative design."
    FALSE
  } else{
    #"A BIBD exists."
    TRUE
  }
}


######### check for BIBD along a grid of panel sizes #########
find_BIBD <- function(min_panel_size = 10, max_panel_size = 30, number_applicants, number_reviews_per_applicant = 2){
  panel_grid <- min_panel_size:max_panel_size
  BIBD_exists <- rep(FALSE, length(panel_grid))
  for(i in 1:length(panel_grid)){
    BIBD_exists[i] <- check_for_BIBD(number_applicants = number_applicants, number_reviewers = panel_grid[i], number_reviews_per_applicant = number_reviews_per_applicant)
  }
  
  if(all(!BIBD_exists)){
    "No BIBD exists for any size panel between min_panel_size and max_panel_size"
  } else{
    data.frame("Panel Size" = panel_grid[BIBD_exists], "Number of Applicants" = number_applicants, "Reviews per Applicant" = number_reviews_per_applicant)
  }
}



