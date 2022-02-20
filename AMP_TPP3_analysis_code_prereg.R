### Confirmatory data analysis code of the AMP-TPP2 project at preregistration

#############################################################
#                                                           #
#                        Packages                           #
#                                                           #
#############################################################

library(lme4) # for glmer()
library(tidyverse)


#############################################################
#                                                           #
#                   Custom functions                        #
#                                                           #
#############################################################


### Functions for Bayes factor caclulation using beta prior
# These functions are required to run the Bayes factor analysis
# The custom code is necessary because we use beta priors, and 
# the BayesFactor package by default does not have built in beta priors
# We thank Richard Morey for his help in developing these functions!


fullAlt_beta = Vectorize(function(p, y, N, alpha, beta){
  exp(dbinom(y, N, p, log = TRUE) + dbeta(p, alpha, beta, log = TRUE)) 
},"p")

normalize_beta = function(alpha, beta, interval){
  diff(pbeta(interval, alpha, beta))
}

restrictedAlt_beta = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  fullAlt_beta(p, y, N, alpha, beta) / normalize_beta(alpha, beta, interval) * (p>interval[1] & p<interval[2])
}

margLike_beta = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_beta, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")


### Function calculating the highest density interval using sampling
# We use hdi() from the library(HDInterval)
# this function is needed for the Bayesian parameter estimation robustness test

mode_HDI <- function(scale, density, crit_width = 0.95, n_samples = 1e5){
  samp <- sample(x = scale, size = n_samples, replace = TRUE, prob = density)
  hdi_result = hdi(samp, credMass=crit_width)
  result = c(scale[which(density == max(density))], # mode
             hdi_result[1], # lower bound
             hdi_result[2]) # upper bound
  
  # only needed for the names of the result
  Crit_lb = (1-crit_width)/2
  Crit_ub = crit_width + (1-crit_width)/2
  
  names(result) = c("mode", paste(Crit_lb*100, "%", sep = ""), paste(Crit_ub*100, "%", sep = ""))
  return(result)
}


# rule for hypothesis testing inference for Bayesian proportion tests
BF_inference_function = function(BF){
  if(Inference_threshold_BF_low >= BF) {return("M1")
  } else if(Inference_threshold_BF_high <= BF) {return("M0")
  } else {return("Inconclusive")}
}


### to convert logit to probability
### this is used for conversion of the results of the
### logistic regression to the probability scale

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


#############################################################
#                                                           #
#                Load and manage data                       #
#                                                           #
#############################################################


raw_data = read.csv("DATA ACCESS PATH")


raw_data[,"sides_match"] = as.factor(tolower(as.logical(raw_data[,"sides_match"])))
raw_data[,"participant_ID"] = as.factor(raw_data[,"participant_ID"])

# sides matching as a numerical variable
raw_data[,"sides_match_numeric"] = as.numeric(as.logical(raw_data[,"sides_match"]))

raw_data = raw_data %>% 
  mutate(sides_match_sign = recode(sides_match_numeric,
                                   "0" = -1,
                                   "1" = 1))


# only sessions conducted with the AMP-TPP3 account is included
data_nontest = raw_data %>% 
  filter(experimenter_ID_code == "9a406d975878bdb99a5654687abb56fdfeff48289ccb5292f2639266a75e016d")


######################################################################
#                                                                    #
#                    AMP-TPP 1 replication analysis                   #
#                                                                    #
######################################################################

##################################################################
#                      Set analysis parameters                   #
##################################################################

# number of erotic trials performed per participant
erotic_trial_size_per_participant_per_trial_type = 9

# probability of successful guess if M0 is true
M0_prob = 0.5

# interim analysis points (in total number of erotic trials performed)
when_to_check = c(37836, 62388, 86958)

# thresholds to infer support for M0 (high) or M1 (low)
Inference_threshold_BF_high = 25
Inference_threshold_BF_low = 1/Inference_threshold_BF_high

# this information is used both for calculating replication Bayes factor, and the Bayesian parameter estimation robustness test. 
# Here we use data from Bem's experiment 1, 828 successes within 1560 erotic trials
y_prior = 828 #number of successes in erotic trials in Bem's experiment 1
N_prior = 1560 # number of erotic trials in Bem's experiment 1


# smallest effect size of interest in the NHST equivalence tests 
# these are used in both the mixed model analysis in the primary analysis
# and the proportion test in the robustness analysis
minimum_effect_threshold_NHST = 0.01
# p threshold for the NHST tests
# these are used in both the mixed model analysis in the primary analysis
# and the proportion test in the robustness analysis
# although, in the primary analysis this is adjusted for multiple testing
# using Bonferroni's correction
Inference_threshold_NHST_AMP1rep = 0.005

##################################################################
#                      Data management                           #
##################################################################

# only pure sessions are included to replicate AMP-TPP2
data_nontest_AMP1rep = data_nontest[data_nontest$available_trial_type == 4,]

# add a row_counter, which will be useful to distinguish data coming in after the stopping rule was met.
data_nontest_AMP1rep[, "row_counter"] = 1:nrow(data_nontest_AMP1rep)

# extract trial data only
data_nontest_AMP1rep_trials = data_nontest_AMP1rep[!is.na(data_nontest_AMP1rep[, "trial_number"]),]

# extract data from erotic trials 
data_nontest_AMP1rep_trials_erotic = data_nontest_AMP1rep_trials[data_nontest_AMP1rep_trials[, "reward_type"] == "erotic", ]
# drop unused factor levels
data_nontest_AMP1rep_trials_erotic[,"participant_ID"] = droplevels(data_nontest_AMP1rep_trials_erotic[,"participant_ID"])

# extract data from nonerotic trials 
data_nontest_AMP1rep_trials_nonerotic = data_nontest_AMP1rep_trials[data_nontest_AMP1rep_trials[, "reward_type"] == "neutral", ]
# drop unused factor levels
data_nontest_AMP1rep_trials_nonerotic[,"participant_ID"] = droplevels(data_nontest_AMP1rep_trials_nonerotic[,"participant_ID"])

# split data by trial type
data_nontest_AMP1rep_trials_erotic_sham = data_nontest_AMP1rep_trials_erotic[data_nontest_AMP1rep_trials_erotic[,"trial_type"] == "sh",]
data_nontest_AMP1rep_trials_erotic_true = data_nontest_AMP1rep_trials_erotic[data_nontest_AMP1rep_trials_erotic[,"trial_type"] == "t",]

data_nontest_AMP1rep_trials_nonerotic_sham = data_nontest_AMP1rep_trials_nonerotic[data_nontest_AMP1rep_trials_nonerotic[,"trial_type"] == "sh",]
data_nontest_AMP1rep_trials_nonerotic_true = data_nontest_AMP1rep_trials_nonerotic[data_nontest_AMP1rep_trials_nonerotic[,"trial_type"] == "t",]

data_nontest_AMP1rep_trials_sham = data_nontest_AMP1rep_trials[data_nontest_AMP1rep_trials[,"trial_type"] == "sh",]
data_nontest_AMP1rep_trials_true = data_nontest_AMP1rep_trials[data_nontest_AMP1rep_trials[,"trial_type"] == "t",]

# cumulative sum of success sign

data_nontest_AMP1rep_trials_erotic_true$sides_match_sign_cumsum = cumsum(data_nontest_AMP1rep_trials_erotic_true$sides_match_sign)
data_nontest_AMP1rep_trials_erotic_sham$sides_match_sign_cumsum = cumsum(data_nontest_AMP1rep_trials_erotic_sham$sides_match_sign)
data_nontest_AMP1rep_trials_nonerotic_true$sides_match_sign_cumsum = cumsum(data_nontest_AMP1rep_trials_nonerotic_true$sides_match_sign)
data_nontest_AMP1rep_trials_nonerotic_sham$sides_match_sign_cumsum = cumsum(data_nontest_AMP1rep_trials_nonerotic_sham$sides_match_sign)



##################################################################
#                    Confirmatory analysis                       #
##################################################################


############################# Sham trials ############################

# This section conducts the primary confirmatory analysis at each stopping point.
# It also cuts the data at the point where one of the stopping rules has been met.

results_table = data.frame(matrix(NA, nrow = 1, ncol = 7))
names(results_table) = c("Mixed_mod_CIlb", "Mixed_mod_CIub", "mixed_CI_width","BF_replication", "BF_uniform", "BF_BUJ", "checked_at")

# this is a counter to count the number of tests conducted using the mixed model
# due to sequential testing. This is used to adjust the p-value threshold 
# for the number of comparions made
comparisons_Mixed_NHST = 0

for(i in 1:length(when_to_check)){
  
  # determin current stopping point and next stopping point
  current_stopping_point = when_to_check[i]
  if(i < length(when_to_check)){next_stopping_point = when_to_check[i+1]} else {next_stopping_point = "last"}
  print(paste("analyzing at reaching", current_stopping_point, "erotic trials"))
  
  # sampling starting from the beggining of the full simulated dataset (from the first trial of the first participant) 
  # until reaching the next interim analysis point
  data_BF = data_nontest_AMP1rep_trials_erotic_sham[1:current_stopping_point,]
  last_row = data_BF[nrow(data_BF), "row_counter"]
  # number of successes and total N of trials
  successes = sum(as.logical(data_BF[,"sides_match"]))
  total_N = current_stopping_point
  results_table[i, "checked_at"] = current_stopping_point
  
  #================================================================#
  #            Mixed effect logistic regression analysis           #
  #================================================================#
  
  # advance the counter to see how much adjustment needs to be made to the
  # NHST inference threshold due to multiple testing
  comparisons_Mixed_NHST = comparisons_Mixed_NHST + 2 # we add 2 at each sequential stopping point because we do two tests at each stop point, one for M0 and one for M1
  
  # build mixed logistic regression model and extract model coefficient and SE
  mod_mixed = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_BF, family = "binomial")
  estimate_mixed = summary(mod_mixed)$coefficients[1,1]
  se_mixed = summary(mod_mixed)$coefficients[1,2]
  
  # compute confidence interval on the probability scale, and save into results_table
  results_table[i,"mixed_CI_width"] = 1-(Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)
  wald_ci_mixed_logit <- c(estimate_mixed - se_mixed* qnorm(1-((Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)/2)),
                           estimate_mixed + se_mixed* qnorm(1-((Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)/2)))
  wald_ci_mixed = logit2prob(wald_ci_mixed_logit)
  
  SE_in_probability_from_mixed_model = logit2prob(estimate_mixed) - logit2prob(estimate_mixed - se_mixed)
  # just to verify the result, this sould provide very similar result, based on http://www.r-tutor.com/elementary-statistics/interval-estimation/interval-estimate-population-proportion
  # SE_in_probability_from_proportions = sqrt((successes/total_N)∗ (1 − (successes/total_N))/total_N)
  
  
  results_table[i, "Mixed_mod_CIlb"] = wald_ci_mixed[1]
  results_table[i, "Mixed_mod_CIub"] = wald_ci_mixed[2]
  
  
  # Statistical inference based on the results of the mixed model analysis  
  
  minimum_effect = M0_prob+minimum_effect_threshold_NHST
  if(results_table[i, "Mixed_mod_CIub"] < minimum_effect){Mixed_NHST_inference = "M0"
  } else if(results_table[i, "Mixed_mod_CIlb"] > M0_prob){Mixed_NHST_inference = "M1"
  } else {Mixed_NHST_inference = "Inconclusive"}
  
  #================================================================#
  #        Calculating Bayes factors using different priors        #
  #================================================================#
  
  # as determined in the analysis plan, three different prior distributions are used for M1
  # to ensure the robustness of the statistical inference to different analytical choices
  # the same 
  
  ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
  
  BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_replication"] = round(BF_replication, 3)
  BF_replication_inference = BF_inference_function(BF_replication)
  
  
  ### Bayes factor with uniform prior
  # using a non-informative flat prior distribution with alpha = 1 and beta = 1
  
  BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_uniform"] = round(BF_uniform, 3)
  BF_uniform_inference = BF_inference_function(BF_uniform)
  
  ### Bayes factor with BUJ prior
  # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
  # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
  # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
  # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
  # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
  # This distribution has 90% of its probability mass under p = 0.712, which we determined 
  # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
  # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
  # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
  # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
  # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
  
  BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_BUJ"] = round(BF_BUJ, 3)
  BF_BUJ_inference = BF_inference_function(BF_BUJ)
  
  
  #================================================================#
  #                    Main analysis inference                     #
  #================================================================#
  
  # determine final inference (supported model) based on the inferences drawn
  # from the mixed model and the Bayes factors 
  if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M1")) {
    primary_analysis_inference = "M1"
    which_threshold_passed = as.character(Inference_threshold_BF_low)
    print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
    break} else if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M0")) {
      primary_analysis_inference = "M0"
      which_threshold_passed = as.character(Inference_threshold_BF_high)
      print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
      break} else if((next_stopping_point != "last") & (nrow(data_nontest_AMP1rep_trials_erotic_sham) < next_stopping_point)){
        primary_analysis_inference = "Ongoing"
        which_threshold_passed = "Ongoing"
        print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
        break} else {
          primary_analysis_inference = "Inconclusive"
          which_threshold_passed = paste("neither ", Inference_threshold_BF_low, " or ", Inference_threshold_BF_high, sep = "")
          print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))}
  
}

data_BF_sham = data_BF
last_row_sham = last_row
results_table_sham = results_table
primary_analysis_inference_sham = primary_analysis_inference
total_N_sham = total_N
successes_sham = successes
SE_in_probability_from_mixed_model_sham = SE_in_probability_from_mixed_model
CIlb_in_probability_from_mixed_model_sham = wald_ci_mixed[1]
CIub_in_probability_from_mixed_model_sham = wald_ci_mixed[2]
which_threshold_passed_sham = which_threshold_passed
BF_replication_sham = BF_replication
BF_uniform_sham = BF_uniform
BF_BUJ_sham = BF_BUJ
Inference_threshold_NHST_AMP1rep_final_sham = Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST


############################# True trials ############################

# This section conducts the primary confirmatory analysis at each stopping point.
# It also cuts the data at the point where one of the stopping rules has been met.

results_table = data.frame(matrix(NA, nrow = 1, ncol = 7))
names(results_table) = c("Mixed_mod_CIlb", "Mixed_mod_CIub", "mixed_CI_width","BF_replication", "BF_uniform", "BF_BUJ", "checked_at")

# this is a counter to count the number of tests conducted using the mixed model
# due to sequential testing. This is used to adjust the p-value threshold 
# for the number of comparions made
comparisons_Mixed_NHST = 0

for(i in 1:length(when_to_check)){
  
  # determin current stopping point and next stopping point
  current_stopping_point = when_to_check[i]
  if(i < length(when_to_check)){next_stopping_point = when_to_check[i+1]} else {next_stopping_point = "last"}
  print(paste("analyzing at reaching", current_stopping_point, "erotic trials"))
  
  # sampling starting from the beggining of the full simulated dataset (from the first trial of the first participant) 
  # until reaching the next interim analysis point
  data_BF = data_nontest_AMP1rep_trials_erotic_true[1:current_stopping_point,]
  last_row = data_BF[nrow(data_BF), "row_counter"]
  # number of successes and total N of trials
  successes = sum(as.logical(data_BF[,"sides_match"]))
  total_N = current_stopping_point
  results_table[i, "checked_at"] = current_stopping_point
  
  #================================================================#
  #            Mixed effect logistic regression analysis           #
  #================================================================#
  
  # advance the counter to see how much adjustment needs to be made to the
  # NHST inference threshold due to multiple testing
  comparisons_Mixed_NHST = comparisons_Mixed_NHST + 2 # we add 2 at each sequential stopping point because we do two tests at each stop point, one for M0 and one for M1
  
  # build mixed logistic regression model and extract model coefficient and SE
  mod_mixed = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_BF, family = "binomial")
  estimate_mixed = summary(mod_mixed)$coefficients[1,1]
  se_mixed = summary(mod_mixed)$coefficients[1,2]
  
  # compute confidence interval on the probability scale, and save into results_table
  results_table[i,"mixed_CI_width"] = 1-(Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)
  wald_ci_mixed_logit <- c(estimate_mixed - se_mixed* qnorm(1-((Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)/2)),
                           estimate_mixed + se_mixed* qnorm(1-((Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST)/2)))
  wald_ci_mixed = logit2prob(wald_ci_mixed_logit)
  
  SE_in_probability_from_mixed_model = logit2prob(estimate_mixed) - logit2prob(estimate_mixed - se_mixed)
  # just to verify the result, this sould provide very similar result, based on http://www.r-tutor.com/elementary-statistics/interval-estimation/interval-estimate-population-proportion
  # SE_in_probability_from_proportions = sqrt((successes/total_N)∗ (1 − (successes/total_N))/total_N)
  
  
  results_table[i, "Mixed_mod_CIlb"] = wald_ci_mixed[1]
  results_table[i, "Mixed_mod_CIub"] = wald_ci_mixed[2]
  
  
  # Statistical inference based on the results of the mixed model analysis  
  
  minimum_effect = M0_prob+minimum_effect_threshold_NHST
  if(results_table[i, "Mixed_mod_CIub"] < minimum_effect){Mixed_NHST_inference = "M0"
  } else if(results_table[i, "Mixed_mod_CIlb"] > M0_prob){Mixed_NHST_inference = "M1"
  } else {Mixed_NHST_inference = "Inconclusive"}
  

  
  #================================================================#
  #        Calculating Bayes factors using different priors        #
  #================================================================#
  
  # as determined in the analysis plan, three different prior distributions are used for M1
  # to ensure the robustness of the statistical inference to different analytical choices
  # the same 
  
  ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
  
  BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_replication"] = round(BF_replication, 3)
  BF_replication_inference = BF_inference_function(BF_replication)
  
  
  ### Bayes factor with uniform prior
  # using a non-informative flat prior distribution with alpha = 1 and beta = 1
  
  BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_uniform"] = round(BF_uniform, 3)
  BF_uniform_inference = BF_inference_function(BF_uniform)
  
  ### Bayes factor with BUJ prior
  # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
  # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
  # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
  # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
  # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
  # This distribution has 90% of its probability mass under p = 0.712, which we determined 
  # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
  # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
  # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
  # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
  # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
  
  BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_BUJ"] = round(BF_BUJ, 3)
  BF_BUJ_inference = BF_inference_function(BF_BUJ)
  
  
  #================================================================#
  #                    Main analysis inference                     #
  #================================================================#
  
  # determine final inference (supported model) based on the inferences drawn
  # from the mixed model and the Bayes factors 
  if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M1")) {
    primary_analysis_inference = "M1"
    which_threshold_passed = as.character(Inference_threshold_BF_low)
    print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
    break} else if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M0")) {
      primary_analysis_inference = "M0"
      which_threshold_passed = as.character(Inference_threshold_BF_high)
      print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
      break} else if((next_stopping_point != "last") & (nrow(data_nontest_AMP1rep_trials_erotic_true) < next_stopping_point)){
        primary_analysis_inference = "Ongoing"
        which_threshold_passed = "Ongoing"
        print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
        break} else {
          primary_analysis_inference = "Inconclusive"
          which_threshold_passed = paste("neither ", Inference_threshold_BF_low, " or ", Inference_threshold_BF_high, sep = "")
          print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))}
  
}


data_BF_true = data_BF
last_row_true = last_row
results_table_true = results_table
primary_analysis_inference_true = primary_analysis_inference
total_N_true = total_N
successes_true = successes
SE_in_probability_from_mixed_model_true = SE_in_probability_from_mixed_model
CIlb_in_probability_from_mixed_model_true = wald_ci_mixed[1]
CIub_in_probability_from_mixed_model_true = wald_ci_mixed[2]
which_threshold_passed_true = which_threshold_passed
BF_replication_true = BF_replication
BF_uniform_true = BF_uniform
BF_BUJ_true = BF_BUJ
Inference_threshold_NHST_AMP1rep_final_true = Inference_threshold_NHST_AMP1rep/comparisons_Mixed_NHST

last_row_whole_study = c(last_row_sham, last_row_true)[which.max(c(last_row_sham, last_row_true))]













######################################################################
#                                                                    #
#                    AMP-TPP 2 replication analysis                  #
#                                                                    #
######################################################################

##################################################################
# Set analysis parameters for AMP-TPP 2 replication analysis     #
##################################################################

max_num_trials_AMP2rep = 127000
Inference_threshold_NHST_AMP2rep = 0.05

##################################################################
#                      Data management                           #
##################################################################

# only pure sessions are included to replicate AMP-TPP2
data_nontest_AMP2rep = data_nontest[data_nontest$available_trial_type == 1,]

# add a row_counter, which will be useful to distinguish data coming in after the stopping rule was met.
data_nontest_AMP2rep[, "row_counter"] = 1:nrow(data_nontest_AMP2rep)

data_nontest_AMP2rep_trials = data_nontest_AMP2rep[!is.na(data_nontest_AMP2rep[, "trial_number"]),]

## extract data from erotic trials 
data_nontest_AMP2rep_trials_erotic = data_nontest_AMP2rep_trials[data_nontest_AMP2rep_trials[, "reward_type"] == "erotic", ]

# drop unused factor levels
data_nontest_AMP2rep_trials_erotic[,"participant_ID"] = droplevels(data_nontest_AMP2rep_trials_erotic[,"participant_ID"])

# drop any data that is above the maximum trial size
if(nrow(data_nontest_AMP2rep_trials_erotic) > max_num_trials_AMP2rep){
  data_nontest_AMP2rep_trials_erotic_maxtrialnum = data_nontest_AMP2rep_trials_erotic[1:max_num_trials_AMP2rep,]
} else {data_nontest_AMP2rep_trials_erotic_maxtrialnum = data_nontest_AMP2rep_trials_erotic}


## extract data from non-erotic trials 
data_nontest_AMP2rep_trials_nonerotic = data_nontest_AMP2rep_trials[data_nontest_AMP2rep_trials[, "reward_type"] == "neutral", ]

# drop unused factor levels
data_nontest_AMP2rep_trials_nonerotic[,"participant_ID"] = droplevels(data_nontest_AMP2rep_trials_nonerotic[,"participant_ID"])

# drop any data that is above the maximum trial size
if(nrow(data_nontest_AMP2rep_trials_nonerotic) > max_num_trials_AMP2rep){
  data_nontest_AMP2rep_trials_nonerotic_maxtrialnum = data_nontest_AMP2rep_trials_nonerotic[1:max_num_trials_AMP2rep,]
} else {data_nontest_AMP2rep_trials_nonerotic_maxtrialnum = data_nontest_AMP2rep_trials_nonerotic}


# cumulative sum of success sign

data_nontest_AMP2rep_trials_erotic_maxtrialnum$sides_match_sign_cumsum = cumsum(data_nontest_AMP2rep_trials_erotic_maxtrialnum$sides_match_sign)
data_nontest_AMP2rep_trials_nonerotic_maxtrialnum$sides_match_sign_cumsum = cumsum(data_nontest_AMP2rep_trials_nonerotic_maxtrialnum$sides_match_sign)


##################################################################
#                    Confirmatory analysis                       #
##################################################################

#### Hypothesis 1

### Primary confirmatory analysis: mixed model binary logistic regression

mod_mixed_H1_AMP2rep = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_nontest_AMP2rep_trials_erotic_maxtrialnum, family = "binomial")

estimate_mixed_H1_AMP2rep = summary(mod_mixed_H1_AMP2rep)$coefficients[1,1]
se_mixed_H1_AMP2rep = summary(mod_mixed_H1_AMP2rep)$coefficients[1,2]

wald_ci_mixed_logit_H1_AMP2rep <- c(estimate_mixed_H1_AMP2rep - se_mixed_H1_AMP2rep* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)),
                            estimate_mixed_H1_AMP2rep + se_mixed_H1_AMP2rep* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)))
wald_ci_mixed_H1_AMP2rep = logit2prob(wald_ci_mixed_logit_H1_AMP2rep)

CI_lower_mixed_H1_AMP2rep = wald_ci_mixed_H1_AMP2rep[1]
CI_upper_mixed_H1_AMP2rep = wald_ci_mixed_H1_AMP2rep[2]

# results of the mixed model analysis
CI_lower_mixed_H1_AMP2rep
CI_upper_mixed_H1_AMP2rep


# final statistical inference based on the mixed model
if(CI_upper_mixed_H1_AMP2rep < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1_AMP2rep > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
conclusion

### Robustness analysis using binomial test

successes_H1_AMP2rep = sum(as.logical(data_nontest_AMP2rep_trials_erotic_maxtrialnum[,"sides_match"]))
total_n_of_trials_H1_AMP2rep = nrow(data_nontest_AMP2rep_trials_erotic_maxtrialnum)


CI_lower_binomtest_H1_AMP2rep = binom.test(x = successes_H1_AMP2rep, n = total_n_of_trials_H1_AMP2rep, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[1]
CI_upper_binomtest_H1_AMP2rep = binom.test(x = successes_H1_AMP2rep, n = total_n_of_trials_H1_AMP2rep, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[2]

## results of the binomial test
CI_lower_binomtest_H1_AMP2rep
CI_upper_binomtest_H1_AMP2rep



#### Hypothesis 2

### Primary confirmatory analysis: mixed model binary logistic regression

mod_mixed_H2_AMP2rep = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_nontest_AMP2rep_trials_nonerotic_maxtrialnum, family = "binomial")

estimate_mixed_H2_AMP2rep = summary(mod_mixed_H2_AMP2rep)$coefficients[1,1]
se_mixed_H2_AMP2rep = summary(mod_mixed_H2_AMP2rep)$coefficients[1,2]

wald_ci_mixed_logit_H2_AMP2rep <- c(estimate_mixed_H2_AMP2rep - se_mixed_H2_AMP2rep* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)),
                            estimate_mixed_H2_AMP2rep + se_mixed_H2_AMP2rep* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)))
wald_ci_mixed_H2_AMP2rep = logit2prob(wald_ci_mixed_logit_H2_AMP2rep)

CI_lower_mixed_H2_AMP2rep = wald_ci_mixed_H2_AMP2rep[1]
CI_upper_mixed_H2_AMP2rep = wald_ci_mixed_H2_AMP2rep[2]

# results of the mixed model analysis
CI_lower_mixed_H2_AMP2rep
CI_upper_mixed_H2_AMP2rep


# final statistical inference based on the mixed model
if(CI_upper_mixed_H2_AMP2rep < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H2_AMP2rep > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
conclusion

### Robustness analysis using binomial test

successes_H2_AMP2rep = sum(as.logical(data_nontest_AMP2rep_trials_nonerotic_maxtrialnum[,"sides_match"]))
total_n_of_trials_H2_AMP2rep = nrow(data_nontest_AMP2rep_trials_nonerotic_maxtrialnum)


CI_lower_binomtest_H2_AMP2rep = binom.test(x = successes_H2_AMP2rep, n = total_n_of_trials_H2_AMP2rep, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[1]
CI_upper_binomtest_H2_AMP2rep = binom.test(x = successes_H2_AMP2rep, n = total_n_of_trials_H2_AMP2rep, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[2]

## results of the binomial test
CI_lower_binomtest_H2_AMP2rep
CI_upper_binomtest_H2_AMP2rep







######################################################################
#                                                                    #
#       AMP-TPP1 replication analysis CI change over time plot       #
#                                                                    #
######################################################################

############################# True trials ############################

trial_cap_AMP1rep = c(seq(100, total_N_true, by = 100), total_N_true)

wide_plot_data_erotic_AMP1rep_true = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)



wide_plot_data_nonerotic_AMP1rep_true = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)


wide_plot_data_erotic_AMP1rep_sham = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)



wide_plot_data_nonerotic_AMP1rep_sham = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)



for(i in 1:length(trial_cap_AMP1rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP1rep[i]))
  
  data_segment_true = data_nontest_AMP1rep_trials_erotic_true[1:trial_cap_AMP1rep[i],]
  
  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1_true = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment_true, family = "binomial")
  
  estimate_mixed_H1_true = summary(mod_mixed_H1_true)$coefficients[1,1]
  se_mixed_H1_true = summary(mod_mixed_H1_true)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1_true <- c(estimate_mixed_H1_true - se_mixed_H1_true* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_true)/2)),
                              estimate_mixed_H1_true + se_mixed_H1_true* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_true)/2)))
  wald_ci_mixed_H1_true = logit2prob(wald_ci_mixed_logit_H1_true)
  
  CI_lower_mixed_H1_true = wald_ci_mixed_H1_true[1]
  CI_upper_mixed_H1_true = wald_ci_mixed_H1_true[2]
  
  # results of the mixed model analysis
  wide_plot_data_erotic_AMP1rep_true[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1_true)
  wide_plot_data_erotic_AMP1rep_true[i, "CI_lb_mixmod"] = CI_lower_mixed_H1_true
  wide_plot_data_erotic_AMP1rep_true[i, "CI_ub_mixmod"] = CI_upper_mixed_H1_true
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1_true < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1_true > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1_true = sum(as.logical(data_segment_true[,"sides_match"]))
  total_n_of_trials_H1_true = nrow(data_segment_true)
  
  
  CI_lower_binomtest_H1_true = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$conf.int[1]
  CI_upper_binomtest_H1_true = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_erotic_AMP1rep_true[i, "probability_of_success_binom"] = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$estimate
  wide_plot_data_erotic_AMP1rep_true[i, "CI_lb_binom"] = CI_lower_binomtest_H1_true
  wide_plot_data_erotic_AMP1rep_true[i, "CI_ub_binom"] = CI_upper_binomtest_H1_true
  
  wide_plot_data_erotic_AMP1rep_true[i, "trial_number"] = trial_cap_AMP1rep[i]
}


for(i in 1:length(trial_cap_AMP1rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP1rep[i]))
  
  data_segment_true = data_nontest_AMP1rep_trials_nonerotic_true[1:trial_cap_AMP1rep[i],]
  
  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1_true = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment_true, family = "binomial")
  
  estimate_mixed_H1_true = summary(mod_mixed_H1_true)$coefficients[1,1]
  se_mixed_H1_true = summary(mod_mixed_H1_true)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1_true <- c(estimate_mixed_H1_true - se_mixed_H1_true* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_true)/2)),
                                   estimate_mixed_H1_true + se_mixed_H1_true* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_true)/2)))
  wald_ci_mixed_H1_true = logit2prob(wald_ci_mixed_logit_H1_true)
  
  CI_lower_mixed_H1_true = wald_ci_mixed_H1_true[1]
  CI_upper_mixed_H1_true = wald_ci_mixed_H1_true[2]
  
  # results of the mixed model analysis
  wide_plot_data_nonerotic_AMP1rep_true[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1_true)
  wide_plot_data_nonerotic_AMP1rep_true[i, "CI_lb_mixmod"] = CI_lower_mixed_H1_true
  wide_plot_data_nonerotic_AMP1rep_true[i, "CI_ub_mixmod"] = CI_upper_mixed_H1_true
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1_true < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1_true > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1_true = sum(as.logical(data_segment_true[,"sides_match"]))
  total_n_of_trials_H1_true = nrow(data_segment_true)
  
  
  CI_lower_binomtest_H1_true = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$conf.int[1]
  CI_upper_binomtest_H1_true = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_nonerotic_AMP1rep_true[i, "probability_of_success_binom"] = binom.test(x = successes_H1_true, n = total_n_of_trials_H1_true, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_true))$estimate
  wide_plot_data_nonerotic_AMP1rep_true[i, "CI_lb_binom"] = CI_lower_binomtest_H1_true
  wide_plot_data_nonerotic_AMP1rep_true[i, "CI_ub_binom"] = CI_upper_binomtest_H1_true
  
  wide_plot_data_nonerotic_AMP1rep_true[i, "trial_number"] = trial_cap_AMP1rep[i]
}


for(i in 1:length(trial_cap_AMP1rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP1rep[i]))
  
  data_segment_sham = data_nontest_AMP1rep_trials_erotic_sham[1:trial_cap_AMP1rep[i],]
  
  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1_sham = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment_sham, family = "binomial")
  
  estimate_mixed_H1_sham = summary(mod_mixed_H1_sham)$coefficients[1,1]
  se_mixed_H1_sham = summary(mod_mixed_H1_sham)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1_sham <- c(estimate_mixed_H1_sham - se_mixed_H1_sham* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_sham)/2)),
                                   estimate_mixed_H1_sham + se_mixed_H1_sham* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_sham)/2)))
  wald_ci_mixed_H1_sham = logit2prob(wald_ci_mixed_logit_H1_sham)
  
  CI_lower_mixed_H1_sham = wald_ci_mixed_H1_sham[1]
  CI_upper_mixed_H1_sham = wald_ci_mixed_H1_sham[2]
  
  # results of the mixed model analysis
  wide_plot_data_erotic_AMP1rep_sham[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1_sham)
  wide_plot_data_erotic_AMP1rep_sham[i, "CI_lb_mixmod"] = CI_lower_mixed_H1_sham
  wide_plot_data_erotic_AMP1rep_sham[i, "CI_ub_mixmod"] = CI_upper_mixed_H1_sham
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1_sham < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1_sham > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1_sham = sum(as.logical(data_segment_sham[,"sides_match"]))
  total_n_of_trials_H1_sham = nrow(data_segment_sham)
  
  
  CI_lower_binomtest_H1_sham = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$conf.int[1]
  CI_upper_binomtest_H1_sham = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_erotic_AMP1rep_sham[i, "probability_of_success_binom"] = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$estimate
  wide_plot_data_erotic_AMP1rep_sham[i, "CI_lb_binom"] = CI_lower_binomtest_H1_sham
  wide_plot_data_erotic_AMP1rep_sham[i, "CI_ub_binom"] = CI_upper_binomtest_H1_sham
  
  wide_plot_data_erotic_AMP1rep_sham[i, "trial_number"] = trial_cap_AMP1rep[i]
}


for(i in 1:length(trial_cap_AMP1rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP1rep[i]))
  
  data_segment_sham = data_nontest_AMP1rep_trials_nonerotic_sham[1:trial_cap_AMP1rep[i],]
  
  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1_sham = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment_sham, family = "binomial")
  
  estimate_mixed_H1_sham = summary(mod_mixed_H1_sham)$coefficients[1,1]
  se_mixed_H1_sham = summary(mod_mixed_H1_sham)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1_sham <- c(estimate_mixed_H1_sham - se_mixed_H1_sham* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_sham)/2)),
                                   estimate_mixed_H1_sham + se_mixed_H1_sham* qnorm(1-((Inference_threshold_NHST_AMP1rep_final_sham)/2)))
  wald_ci_mixed_H1_sham = logit2prob(wald_ci_mixed_logit_H1_sham)
  
  CI_lower_mixed_H1_sham = wald_ci_mixed_H1_sham[1]
  CI_upper_mixed_H1_sham = wald_ci_mixed_H1_sham[2]
  
  # results of the mixed model analysis
  wide_plot_data_nonerotic_AMP1rep_sham[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1_sham)
  wide_plot_data_nonerotic_AMP1rep_sham[i, "CI_lb_mixmod"] = CI_lower_mixed_H1_sham
  wide_plot_data_nonerotic_AMP1rep_sham[i, "CI_ub_mixmod"] = CI_upper_mixed_H1_sham
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1_sham < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1_sham > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1_sham = sum(as.logical(data_segment_sham[,"sides_match"]))
  total_n_of_trials_H1_sham = nrow(data_segment_sham)
  
  
  CI_lower_binomtest_H1_sham = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$conf.int[1]
  CI_upper_binomtest_H1_sham = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_nonerotic_AMP1rep_sham[i, "probability_of_success_binom"] = binom.test(x = successes_H1_sham, n = total_n_of_trials_H1_sham, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP1rep_final_sham))$estimate
  wide_plot_data_nonerotic_AMP1rep_sham[i, "CI_lb_binom"] = CI_lower_binomtest_H1_sham
  wide_plot_data_nonerotic_AMP1rep_sham[i, "CI_ub_binom"] = CI_upper_binomtest_H1_sham
  
  wide_plot_data_nonerotic_AMP1rep_sham[i, "trial_number"] = trial_cap_AMP1rep[i]
}




wide_plot_data_AMP1rep_true = rbind(wide_plot_data_erotic_AMP1rep_true, wide_plot_data_nonerotic_AMP1rep_true)
wide_plot_data_AMP1rep_true$reward_type = c(rep("erotic", nrow(wide_plot_data_erotic_AMP1rep_true)), rep("nonerotic", nrow(wide_plot_data_nonerotic_AMP1rep_true)))

wide_plot_data_AMP1rep_sham = rbind(wide_plot_data_erotic_AMP1rep_sham, wide_plot_data_nonerotic_AMP1rep_sham)
wide_plot_data_AMP1rep_sham$reward_type = c(rep("erotic", nrow(wide_plot_data_erotic_AMP1rep_sham)), rep("nonerotic", nrow(wide_plot_data_nonerotic_AMP1rep_sham)))

wide_plot_data_AMP1rep = rbind(wide_plot_data_AMP1rep_true, wide_plot_data_AMP1rep_sham)
wide_plot_data_AMP1rep$trial_type = c(rep("true", nrow(wide_plot_data_AMP1rep_true)), rep("sham", nrow(wide_plot_data_AMP1rep_sham)))



wide_plot_data_AMP1rep %>% 
  ggplot() +
  aes(x = trial_number, y = probability_of_success_mixmod, fill = reward_type) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = CI_lb_mixmod, ymax = CI_ub_mixmod), alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(trial_type ~ .)


######################################################################
#                                                                    #
#       AMP-TPP1 replication analysis cumulative sum of sing plot    #
#                                                                    #
######################################################################


cumulative_sums_AMP1rep = c(
  data_nontest_AMP1rep_trials_erotic_true$sides_match_sign_cumsum,
  data_nontest_AMP1rep_trials_erotic_sham$sides_match_sign_cumsum,
  data_nontest_AMP1rep_trials_nonerotic_true$sides_match_sign_cumsum,
  data_nontest_AMP1rep_trials_nonerotic_sham$sides_match_sign_cumsum
)

cumulative_sums_origin_AMP1rep = c(rep("erotic_true", length(data_nontest_AMP1rep_trials_erotic_true$sides_match_sign_cumsum)),
                           rep("erotic_sham", length(data_nontest_AMP1rep_trials_erotic_sham$sides_match_sign_cumsum)),
                           rep("nonerotic_true", length(data_nontest_AMP1rep_trials_nonerotic_true$sides_match_sign_cumsum)),
                           rep("nonerotic_sham", length(data_nontest_AMP1rep_trials_nonerotic_sham$sides_match_sign_cumsum)))

trial_number_AMP1rep = c(1:length(data_nontest_AMP1rep_trials_erotic_true$sides_match_sign_cumsum), 
             1:length(data_nontest_AMP1rep_trials_erotic_sham$sides_match_sign_cumsum), 
             1:length(data_nontest_AMP1rep_trials_nonerotic_true$sides_match_sign_cumsum),
             1:length(data_nontest_AMP1rep_trials_nonerotic_sham$sides_match_sign_cumsum))



cumulative_sums_dataframe_AMP1rep = data.frame(cumulative_sums = cumulative_sums_AMP1rep, trial_number = trial_number_AMP1rep, origin = cumulative_sums_origin_AMP1rep)

cumulative_sums_dataframe_AMP1rep %>% 
  ggplot() +
    aes(y = cumulative_sums, x = trial_number, color = origin) +
      geom_line()





######################################################################
#                                                                    #
#       AMP-TPP2 replication analysis CI change over time plot       #
#                                                                    #
######################################################################

trial_cap_AMP2rep = c(seq(100, max_num_trials_AMP2rep, by = 100), max_num_trials_AMP2rep)

wide_plot_data_erotic_AMP2rep = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)



wide_plot_data_nonerotic_AMP2rep = data.frame(
  probability_of_success_mixmod = NA,
  CI_lb_mixmod = NA,
  CI_ub_mixmod = NA,
  probability_of_success_binom = NA,
  CI_lb_binom = NA,
  CI_ub_binom = NA,
  trial_number = NA)



for(i in 1:length(trial_cap_AMP2rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP2rep[i]))
  
  data_segment = data_nontest_AMP2rep_trials_erotic_maxtrialnum[1:trial_cap_AMP2rep[i],]

  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1 = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment, family = "binomial")
  
  estimate_mixed_H1 = summary(mod_mixed_H1)$coefficients[1,1]
  se_mixed_H1 = summary(mod_mixed_H1)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1 <- c(estimate_mixed_H1 - se_mixed_H1* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)),
                              estimate_mixed_H1 + se_mixed_H1* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)))
  wald_ci_mixed_H1 = logit2prob(wald_ci_mixed_logit_H1)
  
  CI_lower_mixed_H1 = wald_ci_mixed_H1[1]
  CI_upper_mixed_H1 = wald_ci_mixed_H1[2]
  
  # results of the mixed model analysis
  wide_plot_data_erotic_AMP2rep[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1)
  wide_plot_data_erotic_AMP2rep[i, "CI_lb_mixmod"] = CI_lower_mixed_H1
  wide_plot_data_erotic_AMP2rep[i, "CI_ub_mixmod"] = CI_upper_mixed_H1
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1 < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1 > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1 = sum(as.logical(data_segment[,"sides_match"]))
  total_n_of_trials_H1 = nrow(data_segment)
  
  
  CI_lower_binomtest_H1 = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[1]
  CI_upper_binomtest_H1 = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_erotic_AMP2rep[i, "probability_of_success_binom"] = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$estimate
  wide_plot_data_erotic_AMP2rep[i, "CI_lb_binom"] = CI_lower_binomtest_H1
  wide_plot_data_erotic_AMP2rep[i, "CI_ub_binom"] = CI_upper_binomtest_H1
  
  wide_plot_data_erotic_AMP2rep[i, "trial_number"] = trial_cap_AMP2rep[i]
}




for(i in 1:length(trial_cap_AMP2rep)){
  
  print(paste0("computing for trial = ", trial_cap_AMP2rep[i]))
  
  data_segment = data_nontest_AMP2rep_trials_nonerotic_maxtrialnum[1:trial_cap_AMP2rep[i],]

  ### Primary confirmatory analysis: mixed model binary logistic regression
  
  mod_mixed_H1 = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_segment, family = "binomial")
  
  estimate_mixed_H1 = summary(mod_mixed_H1)$coefficients[1,1]
  se_mixed_H1 = summary(mod_mixed_H1)$coefficients[1,2]
  
  wald_ci_mixed_logit_H1 <- c(estimate_mixed_H1 - se_mixed_H1* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)),
                              estimate_mixed_H1 + se_mixed_H1* qnorm(1-((Inference_threshold_NHST_AMP2rep)/2)))
  wald_ci_mixed_H1 = logit2prob(wald_ci_mixed_logit_H1)
  
  CI_lower_mixed_H1 = wald_ci_mixed_H1[1]
  CI_upper_mixed_H1 = wald_ci_mixed_H1[2]
  
  # results of the mixed model analysis
  wide_plot_data_nonerotic_AMP2rep[i, "probability_of_success_mixmod"] = logit2prob(estimate_mixed_H1)
  wide_plot_data_nonerotic_AMP2rep[i, "CI_lb_mixmod"] = CI_lower_mixed_H1
  wide_plot_data_nonerotic_AMP2rep[i, "CI_ub_mixmod"] = CI_upper_mixed_H1
  
  
  
  # final statistical inference based on the mixed model
  if(CI_upper_mixed_H1 < M0_prob){conclusion = "M1"} else if(CI_lower_mixed_H1 > M0_prob){conclusion = "M1"} else {conclusion = "Inconclusive"}
  conclusion
  
  ### Robustness analysis using binomial test
  
  successes_H1 = sum(as.logical(data_segment[,"sides_match"]))
  total_n_of_trials_H1 = nrow(data_segment)
  
  
  CI_lower_binomtest_H1 = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[1]
  CI_upper_binomtest_H1 = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$conf.int[2]
  
  ## results of the binomial test
  
  wide_plot_data_nonerotic_AMP2rep[i, "probability_of_success_binom"] = binom.test(x = successes_H1, n = total_n_of_trials_H1, p = 0.5, conf.level = (1-Inference_threshold_NHST_AMP2rep))$estimate
  wide_plot_data_nonerotic_AMP2rep[i, "CI_lb_binom"] = CI_lower_binomtest_H1
  wide_plot_data_nonerotic_AMP2rep[i, "CI_ub_binom"] = CI_upper_binomtest_H1
  
  wide_plot_data_nonerotic_AMP2rep[i, "trial_number"] = trial_cap_AMP2rep[i]
}


wide_plot_data_AMP2rep = rbind(wide_plot_data_erotic_AMP2rep, wide_plot_data_nonerotic_AMP2rep)
wide_plot_data_AMP2rep$reward_type = c(rep("erotic", nrow(wide_plot_data_erotic_AMP2rep)), rep("nonerotic", nrow(wide_plot_data_nonerotic_AMP2rep)))


wide_plot_data_AMP2rep %>% 
  ggplot() +
  aes(x = trial_number, y = probability_of_success_mixmod, fill = reward_type) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = CI_lb_mixmod, ymax = CI_ub_mixmod), alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0))


######################################################################
#                                                                    #
#       AMP-TPP2 replication analysis cumulative sum of sing plot    #
#                                                                    #
######################################################################



cumulative_sums_AMP2rep = c(
  data_nontest_AMP2rep_trials_erotic_maxtrialnum$sides_match_sign_cumsum,
  data_nontest_AMP2rep_trials_nonerotic_maxtrialnum$sides_match_sign_cumsum
)

cumulative_sums_origin_AMP2rep = c(rep("erotic", length(data_nontest_AMP2rep_trials_erotic_maxtrialnum$sides_match_sign_cumsum)),
                                   rep("nonerotic", length(data_nontest_AMP2rep_trials_nonerotic_maxtrialnum$sides_match_sign_cumsum)))

trial_number_AMP2rep = c(1:length(data_nontest_AMP2rep_trials_erotic_maxtrialnum$sides_match_sign_cumsum), 
                         1:length(data_nontest_AMP2rep_trials_nonerotic_maxtrialnum$sides_match_sign_cumsum))



cumulative_sums_dataframe_AMP2rep = data.frame(cumulative_sums = cumulative_sums_AMP2rep, trial_number = trial_number_AMP2rep, origin = cumulative_sums_origin_AMP2rep)

cumulative_sums_dataframe_AMP2rep %>% 
  ggplot() +
  aes(y = cumulative_sums, x = trial_number, color = origin) +
  geom_line()



