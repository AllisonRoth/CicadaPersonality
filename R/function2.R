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
  boot_res <- bootstrap.merMod(model = model, .f = my_get, type = "case", B = bootN, resample = c(TRUE, FALSE))
  
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

# from = https://github.com/aloy/lmeresampler/blob/master/R/bootstrap_lme4.R

#' @rdname bootstrap
#' @export
#' @method bootstrap merMod
#' @importFrom stats as.formula cov formula model.matrix na.exclude 
#' na.omit predict resid simulate sd confint quantile
bootstrap.merMod <- function(model, .f, type, B, resample, reb_type, hccme, aux.dist){
  switch(type,
         parametric = parametric_bootstrap.merMod(model, .f, B),
         residual = resid_bootstrap.merMod(model, .f, B),
         case = case_bootstrap.merMod(model, .f, B, resample),
         reb = reb_bootstrap.lmerMod(model, .f, B, reb_type),
         wild = wild_bootstrap.lmerMod(model, .f, B, hccme, aux.dist))
}


#' @rdname parametric_bootstrap
#' @export
#' @method parametric_bootstrap merMod
parametric_bootstrap.merMod <- function(model, .f, B){
  .f <- match.fun(.f)
  
  # model.fixef <- lme4::fixef(model) # Extract fixed effects
  ystar <- simulate(model, nsim = B, na.action = na.exclude)
  
  # refit here
  refits <- refit_merMod(ystar, model, .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "parametric", warnings = refits$warnings)
}


#' @rdname case_bootstrap
#' @export
#' @method case_bootstrap merMod
case_bootstrap.merMod <- function(model, .f, B, resample){
  
  data <- model@frame
  flist <- lme4::getME(model, "flist")
  re_names <- names(flist)
  clusters <- c(rev(re_names), ".id")
  
  if(length(clusters) != length(resample))
    stop("'resample' is not the same length as the number of grouping variables. 
         Please specify whether to resample the data at each level of grouping,
         including at the observation level.")
  
  if(!all(re_names %in% colnames(data))) {
    missing_re <- setdiff(re_names, colnames(data))
    data <- dplyr::bind_cols(data, flist[missing_re])
  }
  
  # rep.data <- purrr::map(integer(B), function(x) .cases.resamp(model = model, dat = data, cluster = clusters, resample = resample))
  refits <- purrr::map(integer(B), function(x) .resample_refit.cases(model = model, .f = .f, dat = data, cluster = clusters, resample = resample))
  
  tstar <- purrr::map(refits, ~.x$value)
  warnings <- list(
    warning = lapply(refits, function(.x) unlist(.x$warning)$message),
    message = lapply(refits, function(.x) unlist(.x$message)$message),
    error = lapply(refits, function(.x) unlist(.x$error)$message)
  )
  
  .bootstrap.completion(model, tstar = tstar, B, .f, type = "case", warnings = warnings)
}



#' @rdname resid_bootstrap
#' @export
#' @method resid_bootstrap merMod
resid_bootstrap.merMod <- function(model, .f, B){
  
  .f <- match.fun(.f)
  glmm <- lme4::isGLMM(model)
  
  setup <- .setup(model, type = "residual")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.cgr(
        glmm = glmm,
        b = setup$b, 
        e = setup$e, 
        level.num = setup$level.num, 
        Ztlist = setup$Ztlist, 
        Xbeta = setup$Xbeta,
        vclist = setup$vclist,
        sig0 = setup$sig0,
        invlink = ifelse(glmm, model@resp$family$linkinv, NULL)
      )
    )
  )
  
  if(glmm){
    fam <- stats::family(model)
    wts <- stats::weights(model)
    
    # simulate y
    simfun <- simfunList[[fam$family]]
    ystar <- purrr::map(ystar,
                        ~simfun(model, nsim = 1, ftd = .x, wts = wts)
    )
    
  }
  
  refits <- refit_merMod(ystar, model, .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "residual", warnings = refits$warnings)
}


#' @rdname wild_bootstrap
#' @export
#' @method wild_bootstrap lmerMod
wild_bootstrap.lmerMod <- function(model, .f, B, hccme = c("hc2", "hc3"), 
                                   aux.dist = c("f1", "f2")){
  
  .f <- match.fun(.f)
  hccme <- match.arg(hccme)
  aux.dist <- match.arg(aux.dist)
  
  setup <- .setup(model, type = "wild")
  
  ystar <- as.data.frame(
    replicate(
      n = B, 
      .resample.wild(
        Xbeta = setup$Xbeta, 
        mresid = setup$mresid, 
        .hatvalues = setup$.hatvalues, 
        hccme = hccme, 
        aux.dist = aux.dist,
        n.lev = setup$n.lev,
        flist = setup$flist
      )
    )
  )
  
  
  refits <- refit_merMod(ystar, model, .f)
  
  .bootstrap.completion(model, tstar = refits$tstar, B, .f, type = "wild", warnings = refits$warnings)
}




#' @rdname reb_bootstrap
#' @export
#' @method reb_bootstrap lmerMod
reb_bootstrap.lmerMod <- function(model, .f, B, reb_type){
  
  if(missing(reb_type)){
    reb_type <- 0
    warning("'reb_type' unspecified, performing REB 0 bootstrap")
  }
  
  if(lme4::getME(object = model, name = "n_rfacs") > 1) {
    stop("The REB bootstrap has not been adapted for 3+ level models.")
  }
  
  .f <- match.fun(.f)
  
  # Set up for bootstrapping
  setup <- .setup(model, type = "reb", reb_type = reb_type)
  
  # Generate bootstrap responses
  y.star <- replicate(
    n = B, 
    .resample.reb(
      Xbeta = setup$Xbeta, 
      Ztlist = setup$Ztlist, 
      Uhat = setup$b, 
      estar.vec = as.numeric(setup$e), 
      flist = setup$flist, 
      levs = setup$levs
    )
  )
  
  ystar <- as.data.frame(y.star)
  
  # Extract bootstrap statistics
  if(reb_type == 2) .f <- extract_parameters.merMod
  
  refits <- refit_merMod(ystar, model, .f)
  tstar <- refits$tstar
  # Extract original statistics
  t0 <- .f(model)
  
  # Postprocessing for REB/2
  if(reb_type == 2) 
    tstar <- .postprocess.reb2(t0, tstar, nbeta = length(lme4::getME(model, "beta")), B = B)
  
  # Format for return
  .bootstrap.completion(model, tstar = tstar, B, .f, 
                        type = paste("reb", reb_type, sep = ""), 
                        warnings = refits$warnings)
  
}


# addition 


#' Case resampler for mixed models
#' @keywords internal
#' @noRd
.resamp.cases <- function(dat, cluster, resample) {
  res <- dat
  
  for(i in seq_along(cluster)) {
    if(i==1 & resample[i]) {
      dots <- as.name(cluster[1])
      grouped <- dplyr::group_by(res, !!dots)
      g_rows <- dplyr::group_rows(grouped)
      cls <- sample(seq_along(g_rows), replace = resample[i])
      idx <- unlist(g_rows[cls], recursive = FALSE)
      res <- res[idx, ]
    } else{
      if(i == length(cluster) & resample[i]) {
        dots <- as.name(cluster[-i])
        grouped <- dplyr::group_by(res, !!dots) 
        res <- dplyr::sample_frac(grouped, size = 1, replace = TRUE)
      } else{
        if(resample[i]) {
          dots <- as.name(cluster[i])
          res <- split(res, res[, cluster[1:(i-1)]], drop = TRUE)
          res <- purrr::map_dfr(res, function(df) { # ldply to purrr map from list to df
            grouped <- dplyr::group_by(df, !!dots)
            g_rows <- dplyr::group_rows(grouped)
            cls <- sample(seq_along(g_rows), replace = resample[i])
            idx <- unlist(g_rows[cls], recursive = FALSE)
            grouped[idx, ]
          }, .id = NULL)
        }
      }
    }
  }
  
  res
}

#' @importFrom catchr catch_expr
.resample_refit.cases <- function(model, .f, dat, cluster, resample){
  resamp_data <- .resamp.cases(dat, cluster, resample)
  error <- NULL
  
  if(class(model) == "lmerMod"){
    # Refit the model and apply '.f' to it using map
    form <- model@call$formula
    reml <- lme4::isREML(model)
    
    tstar <- catchr::catch_expr(
      .f(lme4::lmer(formula = form, data = resamp_data, REML = reml)),
      warning, message, error
    )
    # tstar <- purrr::map(res, function(x) {
    #   .f(lme4::lmer(formula = form, data = as.data.frame(x), REML = reml)) 
    # })
  } else if(class(model) == "lme"){
    tstar <- updated.model(model = model, new.data = resamp_data)  
    tstar$value <- .f(tstar$value)
  } else if(class(model) == "glmerMod") {
    form <- update(model@call$formula, y ~ .)
    colnames(resamp_data)[1] <- "y"
    fam  <- family(model)
    tstar <- catchr::catch_expr(
      .f(lme4::glmer(formula = form, data = resamp_data, family = fam)),
      warning, message, error
    )
  } else{
    stop("model class must be one of 'lme', 'lmerMod', or 'glmerMod'")
  }
  tstar
}



#' CGR resampling from lmerMod objects
#' @keywords internal
#' @noRd
.resample.cgr <- function(glmm, b, e, level.num, Ztlist, Xbeta, vclist, sig0, invlink){
  
  # Resample Uhat
  ustar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  ustar <- purrr::map2(ustar, vclist, scale_center_ranef)
  
  # Structure u*
  if(level.num == 1){
    if(is.data.frame(ustar[[1]])){
      ustar <- purrr::map(ustar, .f = as.list)[[1]]
    }
    names(ustar) <- names(Ztlist)
  } else {
    ustar <- purrr::map(ustar, .f = as.data.frame)
    ustar <- do.call(c, ustar)
    names(ustar) <- names(Ztlist)
  }
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  if(glmm) {
    eta <- as.numeric(Xbeta + Zbstar.sum)
    ystar  <- invlink(eta) # not really ystar, need inv. link
  } else{
    # Resample resids
    estar <- sample(x = e, size = length(e), replace = TRUE)
    ehat <- scale_center_e(estar, sigma = sig0)
    
    # Calc. bootstrap y
    ystar <- as.numeric(Xbeta + Zbstar.sum + estar)
  }
  
  ystar
}



#' Wild bootstrap resampling from lmerMod and lme objects
#' @keywords internal
#' @noRd
.resample.wild <- function(Xbeta, mresid, .hatvalues, hccme, aux.dist, 
                           flist, n.lev){
  
  # Sample from auxillary distribution
  if(aux.dist == "f1") {
    prob <- (sqrt(5) + 1) / (2 * sqrt(5))
    w <- sample(
      c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2), 
      size = n.lev,
      prob = c(prob, 1 - prob),
      replace = TRUE
    )
  } 
  
  if(aux.dist == "f2") {
    w <- sample(c(1, -1), size = n.lev, replace = TRUE)
  }
  
  # Calc. bootstrap y
  if(hccme == "hc2") v <- (1 / sqrt(1 - .hatvalues)) * mresid
  if(hccme == "hc3") v <- (1 / (1 - .hatvalues)) * mresid
  
  as.numeric(Xbeta + v * w[flist])
}


#' REB resampling procedure for lmerMod objects
#' 
#' @param Xbeta marginal fitted values
#' @param Ztlist design matrix separated by variance
#' @param Uhat ranefs organized as Ztlist
#' @param estar.vec vector of level-1 residuals
#' @param flist a list of the grouping variables (factors) involved in the random effect terms
#' @param levs a list of levels of the grouping variables in flist
#' @inheritParams bootstrap
#' @import Matrix
#' @keywords internal
#' @noRd
.resample.reb <- function(Xbeta, Ztlist, Uhat, estar.vec, flist, levs){
  # For now, assume only a two-level model
  grps <- levs[[1]]
  units <- flist[[1]]
  resamp_u_ids <- sample(seq_along(grps), size = length(grps), replace = TRUE)
  resamp_e_grps <- sample(grps, size = length(grps), replace = TRUE)
  
  # resample uhats
  ustar <- purrr::map(Uhat, ~data.frame(.x[resamp_u_ids, ]))
  
  # Resample residuals, e
  estar <- numeric(length = length(units))
  for(i in seq_along(resamp_e_grps)) {
    target_units <- which(units == grps[i])
    donor_units <- which(units == resamp_e_grps[i])
    estar[target_units] <- sample(estar.vec, size = length(target_units), replace = TRUE)
  }
  
  # since only working with 2-levels models now
  ustar <- ustar[[1]]
  
  names(ustar) <- names(Ztlist) 
  ustar.df <- as.data.frame(ustar)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine(bstar = ustar.df, zstar = Ztlist)
  Zbstar.sum <- Reduce("+", Zbstar)
  
  # ystar
  as.numeric(Xbeta + Zbstar.sum + estar)
}


#' Resampling residuals from lme objects
#' @keywords internal
#' @noRd
.resample.resids.lme <- function(b, e, Xbeta, Zlist){
  bstar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  estar <- sample(e, size = length(e), replace = TRUE)
  
  Zbstar.sum <- .Zbstar.combine.lme(bstar, Zlist)
  
  # Combine function
  as.numeric(Xbeta + Zbstar.sum + estar)
}

#' CGR resampler for lme objects
#' @keywords internal
#' @noRd
.resample.cgr.lme <- function(b, e, Xbeta, Zlist, vclist, sig0){
  # Resample resids
  estar <- sample(x = e, size = length(e), replace = TRUE)
  ehat <- scale_center_e(estar, sigma = sig0)
  
  # Resample Uhat
  ustar <- purrr::map(b, .f = dplyr::slice_sample, prop = 1, replace = TRUE)
  ustar <- purrr::map2(ustar, vclist, scale_center_ranef)
  
  Zbstar.sum <- .Zbstar.combine.lme(ustar, Zlist)
  
  # Combine function
  as.numeric(Xbeta + Zbstar.sum + estar)
}


#' REB resampler for lme objects
#' @keywords internal
#' @noRd
.resample.reb.lme <- function(Xbeta, Zlist, Uhat, estar.vec, flist, levs){
  # For now, assume only a two-level model
  grps <- levs[[1]]
  units <- flist[[1]]
  resamp_u_ids <- sample(seq_along(grps), size = length(grps), replace = TRUE)
  resamp_e_grps <- sample(grps, size = length(grps), replace = TRUE)
  
  # resample uhats
  ustar <- purrr::map(Uhat, ~data.frame(.x[resamp_u_ids, ]))
  
  # Resample residuals, e
  estar <- numeric(length = length(units))
  for(i in seq_along(resamp_e_grps)) {
    target_units <- which(units == grps[i])
    donor_units <- which(units == resamp_e_grps[i])
    estar[target_units] <- sample(estar.vec, size = length(target_units), replace = TRUE)
  }
  
  # since only working with 2-levels models now
  # ustar <- ustar[[1]]
  # 
  # names(ustar) <- names(Zlist) 
  # ustar.df <- as.data.frame(ustar)
  
  # Get Zb*
  Zbstar <- .Zbstar.combine.lme(bstar = ustar, Zlist = Zlist)
  # Zbstar.sum <- Reduce("+", Zbstar)
  
  # ystar
  as.numeric(Xbeta + Zbstar + estar)
}

#' @title Bootstrap Completion
#'
#' @description
#' Finishes the bootstrap process and makes the output readable.
#'
#' @details
#' This function is given \code{model, tstar, B, .f} and uses them to complete
#' the bootstrap process. They are then structured into a list for output and returned.
#'
#' @param tstar The tstar being passed in
#' @inheritParams bootstrap
#'
#' @return list
#' @keywords internal
#' @noRd
.bootstrap.completion <- function(model, tstar, B, .f, type = type, warnings){
  t0 <- .f(model)
  
  nsim <- length(tstar)
  
  # tstar <- do.call("cbind", tstar) # Can these be nested?
  # row.names(tstar) <- names(t0)
  
  # if((nfail <- sum(bad.runs <- apply(is.na(tstar), 2, all))) > 0) {
  #   warning("some bootstrap runs failed (", nfail, "/", nsim, ")")
  #   fail.msgs <- purrr::map_chr(tstar[bad.runs], .f = attr, FUN.VALUE = character(1),
  #                               "fail.msgs")
  # } else fail.msgs <- NULL
  
  # prep for stats df
  
  observed <- t0
  
  if(is.numeric(t0)) {
    if(length(t0) == 1) {
      replicates <- unlist(tstar)
      rep.mean <- mean(replicates)
      se <- sd(replicates)
      bias <- rep.mean - observed
      stats <- dplyr::tibble(observed, rep.mean, se, bias)
    } else{
      replicates <- dplyr::bind_rows(tstar)
      rep.mean <- colMeans(replicates)
      se <- unlist(purrr::map(replicates, sd))
      bias <- rep.mean - observed
      stats <- dplyr::tibble(term = names(t0), observed, rep.mean, se, bias)
    }
    
  } else{
    if(is.data.frame(t0)) {
      .ids <- rep(seq_along(tstar), times = vapply(tstar, nrow, FUN.VALUE = 0L))
      replicates <- dplyr::bind_rows(tstar) %>% dplyr::mutate(.n = .ids)
    }
    stats <- NULL
  }
  
  
  if (class(model) == "lme") data <- model$data
  else data <- model@frame
  
  RES <- structure(list(observed = observed, model = model, .f = .f, replicates = replicates,
                        stats = stats, B = B, data = data,
                        seed = .Random.seed, type = type, call = match.call(),
                        message = warnings$message, warning = warnings$warning, error = warnings$error), 
                   class = "lmeresamp")
  
  # attr(RES,"bootFail") <- nfail
  # attr(RES,"boot.fail.msgs") <- fail.msgs
  return(RES)
}
