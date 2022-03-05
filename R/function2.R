# modefiend new function 

# function to get boostraping "percentile" CI - only one random effect is allowed
# we use the case boostrap (fully non-parametric)


#' Title
#'
#' @param model: put a lmer model 
#' @param bootN: the number of boostrapping (defaul = 1000), but you could increase to 10,000 for your real results 
#'
#' @return
#' @export
#'
#' @examples
nonpara_boot2 <- function(model, bootN = 1000){
  
  # helper function for boostrap()
  my_get <- function(model) { 
    beta <- getME(model, "beta")
    sigma <- getME(model, "sigma")
    sigma2 <- sigma^2
    tau2 <- unname(sigma * getME(model, "theta"))^2
    rep <- tau2/(tau2 + sigma2)
    c(beta = beta, sigma2 = sigma2,  tau2 = tau2, rep = rep) 
  }
  
  # boostrapping - type = "case" fully
  # we are only respmpleing "individual" - resample = c(TRUE, FALSE)
  #  for why - see "Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel models. In J. de Leeuw and E. Meijer, editors, Handbook of Multilevel Analysis, pages 401–433. New York: Springer."
  # The book chapter says "The cases bootstrap, finally, requires minimal assumptions: Only the hierarchical dependency in the data is assumed to be specified correctly."
  boot_res <- lmeresampler::bootstrap(model = model, .f = my_get, type = "case", B = bootN, resample = c(TRUE, FALSE))
  
  # get the CI for everything
  paraN <- length(boot_res$observed)
  
  # get the list of CI results
  ci_low <- map_dbl(1:paraN, ~ quantile(boot_res$replicates[[.x]],0.025))
  ci_upr <- map_dbl(1:paraN, ~ quantile(boot_res$replicates[[.x]],0.975))
  est_mid <- map_dbl(1:paraN, ~ median(boot_res$replicates[[.x]]))
  
  # getting names of the parameters
  Name <- c(names(fixef(model)), names(boot_res$observed)[-(1:length(names(fixef(model))))])
  
  # putting all togther
  res <- tibble(Name = Name, Median_NP = as.numeric(boot_res$stats$rep.mean) , CI_L = ci_low, CI_U =  ci_upr)
  res
  
}

# test

#test <- nonpara_boot2(model = lmer1, bootN = 1000)
#test
