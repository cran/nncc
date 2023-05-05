## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4
)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(nncc)
library(survival)
library(dplyr)
library(ggplot2)
data(anifood)

## -----------------------------------------------------------------------------
dim(anifood)

## -----------------------------------------------------------------------------
# exposures of interest. In a real study, this list can be much longer
exp_interest <- c("exp01","exp09", "exp27")

# exposures to be controlled for any exposure of interest
exp_match <- setdiff(names(anifood), "case")

# variables to be excluded from matching for each exposure
# both exp_var and rm_vars are character variables
excl_vars %>% head

## -----------------------------------------------------------------------------
threshold_results <- get_threshold(anifood, exp_match, p_threshold = 0.50)

## ---- warning=FALSE-----------------------------------------------------------
distance_density_plot(threshold_results) + ggtitle("Example of distance_density_plot")

## -----------------------------------------------------------------------------
threshold_model_plot(threshold_results, p_threshold = 0.50) + ggtitle("Example of threshold_model_plot")

## -----------------------------------------------------------------------------
# create a variable (i.e., pair) indicating originally matched pairs
anifood_matched <- anifood %>% group_by(case) %>% mutate(pair = seq_along(case)) %>% ungroup

## -----------------------------------------------------------------------------
p <- original_compare_plot(anifood_matched, case, pair, threshold_results)


# the density plot of distance between originally matched cases and controls
p$plot_density + ggtitle("Example of original_compare_plot")

# proportion of originally matched cases and controls with a distance greater than the threshold
p$prop_distance_gt_threshold

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(furrr)
strata1 <- future_map(exp_interest, make_knn_strata,  rmvars = excl_vars, matchvars = exp_match, df = anifood) %>%
  setNames(exp_interest)

## ---- warning=FALSE-----------------------------------------------------------
length(strata1) == length(exp_interest)

# rows in a matched data set
all.equal(anifood %>% filter(case == 1) %>% NROW %>% `*`(250 + 1), NROW(strata1[[1]]))

## ---- warning = FALSE---------------------------------------------------------
strata2 <- future_map(exp_interest, make_analysis_set, stratified_data = strata1, data = anifood, maxdist = threshold_results$threshold) %>% setNames(exp_interest)

## ---- message=FALSE-----------------------------------------------------------
strata3 <- finalize_data(strata2)

## -----------------------------------------------------------------------------
# exposures to which neither cases nor controls were exposed
strata3 %>%
  lapply(function(dfm) dfm %>%
                            mutate(case = as.character(case), exp = as.character(exp)) %>%
                            filter(case != "" & exp != "") %>% 
                            with(table(case, exp))) %>%
  lapply(function(x) x==0) %>%
  lapply(function(x) {sum = sum(x); length = length(x); cbind(sum, length)}) %>%
  do.call(rbind,.) %>%
  as.data.frame %>%
  mutate(var = names(strata3)) %>%
  filter(sum >= 2 | length < 4) %>% select(var) %>%
  unclass %>%
  unlist -> expvars_invalid

expvars_invalid

# none exposed and odds ratio cannot be estimated
strata3[["exp09"]] %>% with(table(case, exp))

data_final <- strata3[setdiff(names(strata3), expvars_invalid)]

## -----------------------------------------------------------------------------
or_mh <- future_map(data_final, function(dfm) with(dfm, test_mh(case = case, exp = exp, strata = strata)))

or_mh[["exp01"]]

## -----------------------------------------------------------------------------
data_final[["exp27"]] %>% select(case, exp) %>% table

or_mh[["exp27"]]

## ---- warning=FALSE-----------------------------------------------------------
results_clogit <- future_map(data_final, function(dfm){clogit(case ~ exp + strata(strata) , data = dfm)}) 

results_clogit[["exp27"]] %>% summary() %>% `$`(conf.int)

## ---- warning=FALSE-----------------------------------------------------------
# regression coefficients
coef_logistf <- future_map(data_final, function(dfm){

  if(length(unique(dfm[["strata"]])) == 1) {
        x <- model.matrix(case ~ exp, data = dfm)
        plconf <- grep("exp", colnames(x))
        o <- logistf(case ~ exp, data = dfm, plconf = plconf) %>%
                 `[`(c("terms", "coefficients", "ci.lower", "ci.upper", "prob", "call", "loglik", "model"))
    } else {
        x <- model.matrix(case ~ exp + strata, data = dfm)
        plconf <- grep("exp", colnames(x))
        o <- logistf(case ~ exp + strata, data = dfm, plconf = plconf) %>%
                 `[`(c("terms", "coefficients", "ci.lower", "ci.upper", "prob", "call", "loglik", "model"))
    }
},
.progress = TRUE)


# odds ratios
or_logistf <- lapply(names(coef_logistf), function(var_name){
  coef_logistf[[var_name]] %>%
    `[`(c("terms", "coefficients", "ci.lower", "ci.upper", "prob")) %>%
    bind_cols %>%
    filter(grepl("exp", terms)) %>%
    mutate(variable = var_name, or = exp(coefficients), ci.lower = exp(ci.lower), ci.upper = exp(ci.upper)) %>%
  select(variable, terms, or, ci.lower, ci.upper, prob)}) %>%
  setNames(names(coef_logistf))

# odds ratio for exp27
or_logistf[["exp27"]] 

## -----------------------------------------------------------------------------

# prepare a data frame for calculating PAF
df_or <- bind_rows(or_logistf)

df_or

## ---- message=FALSE-----------------------------------------------------------
# point estimate of PAF
paf <- get_paf(df_or = df_or, which_or = or, exp_var = variable, exp_level = terms, df_matched = data_final)

paf

# lower confidence limit of PAF
paf_ci.lower <- get_paf(df_or = df_or, which_or = ci.lower, exp_var = variable, exp_level = terms, df_matched = data_final)

## ---- eval = FALSE------------------------------------------------------------
#  strata1 <- cacheit("abc",
#                     future_map(exp_interest, make_knn_strata,  rmvars = excl_vars, matchvars = exp_match, df = anifood) %>%
#    setNames(exp_interest),
#    clearcache = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  library(nncc)
#  library(dplyr)
#  library(furrr)
#  library(future.batchtools)
#  # the workers argument is used to define the number of cores for the analysis. By default, all cores will be used.
#  plan(multisession, workers = 3)
#  
#  strata1 <- future_map(exp_interest, make_knn_strata,  rmvars = excl_vars, matchvars = exp_match, df = anifood) %>%
#    setNames(exp_interest)
#  

## ---- eval=FALSE, comment = ""------------------------------------------------
#  #!/bin/bash -l
#  
#  # name of the job is helloR
#  #$ -N helloR
#  
#  #$ -cwd
#  
#  #$ -V
#  
#  #$ -pe smp 2-12
#  
#  #$ -q all.q
#  
#  module load R/3.6.2
#  
#  # to execute fugure-map.R
#  Rscript future-map.R
#  
#  exit 0

## ---- eval=FALSE--------------------------------------------------------------
#  library(furrr)
#  library(future.batchtools)
#  
#  plan(batchtools_sge)
#  
#  future_map(<a_list_or_vector>, function(x){...}, .progress = TRUE)

