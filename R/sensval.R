

#' This is data to be included in my package.
#'
#' @references del Ninno et al. (2001), \url{https://doi.org/10.7910/DVN/UJA9N5}
"bangladesh_floods"


#' Classifies numeric variables as binary or continuous.
#'
#' This function is internal to other functions and checks if a particular
#' numeric variable should be treated as binary or continuous.
#'
#' @param x pair-matched or individual-level data.
#' @return Boolean indicator of whether the variable is binary.
is_binary <- function(x) {
  all( is.na(x) | x==0 | abs(x)==1 )
}


#' Performs sensitivity analysis for pair-matched binary data.
#'
#' This function is internal to other functions and performs sensitivity
#' analysis for pair-matched binary data using McNemar's test statistic.
#'
#' @param mpdifs vector of length nsets, where entry k contains the treated
#' minus control value corresponding to the kth matched pair.
#' @param gamma value of the sensitivity parameter.
#' @param alt character denoting direction of alternative hypothesis.
#' @return p-value yielded by sensitivity analysis.
sensitivity_analysis_mcnemar <- function(mpdifs,gamma=gamma,alt=alt) {
  if (alt=='greater'){
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
  }
  else {
    mpdifs <- mpdifs*(-1)
    D=sum(abs(mpdifs)==1)
    Tobs=sum(mpdifs==1)
  }
  p.positive=gamma/(1+gamma);
  upperbound=1-stats::pbinom(Tobs-1,D,p.positive);
  return(upperbound)
}


#' Computes sensitivity value for pair-matched data.
#'
#' This function is internal to other functions and calculates the sensitivity
#' value of pair-matched data using numerical methods.
#'
#' @param mpdifs vector of length nsets, where entry k contains the treated
#' minus control value corresponding to the kth matched pair.
#' @param a alpha-level of test.
#' @param alt character denoting direction of alternative hypothesis.
#' @return Sensitivity value.
#' @export
sensitivity_value <- function(mpdifs,a=a,alt=alt) {
  if (is_binary(mpdifs)){
    ## Use McNemar's statistic for binary data
    if (alt=='greater'){
      D=sum(abs(mpdifs)==1)
      Tobs=sum(mpdifs==1)
    }
    else {
      mpdifs <- mpdifs*(-1)
      D=sum(abs(mpdifs)==1)
      Tobs=sum(mpdifs==1)
    }
    senfx_bin <- function(g,a){
      1-stats::pbinom(Tobs-1,D,g/(1+g))-a
    }
    g<-1
    try(g<-stats::uniroot(senfx_bin,c(1,1e10),a=a)$root,silent=TRUE)
    kappa<-g/(1+g)
    return(kappa)
  }
  else{
    ## Use Wilcoxon's signed rank statistic for continuous data
    senfx_cont <- function(mpdifs,g,a,alt){
      DOS2::senWilcox(mpdifs, gamma=g, alpha=a, alternative=alt)$pval-a
    }
    g<-1
    try(g<-stats::uniroot(senfx_cont,c(1,1e10),mpdifs=mpdifs,a=a,alt=alt)$root,
        silent=TRUE)
    kappa<-g/(1+g)
    return(kappa)
  }
}


#' Sens-Val method for pair-matched data.
#'
#' This function runs the Sens-Val procedure, introduced by Bekerman et al.
#' (2024), on pair-matched data from an observational study. The method 
#' identifies and tests hypotheses from a split sample that are robust to
#' potential unmeasured biases. As defaults, binary data are analyzed using
#' McNemar's test statistic and continuous data are analyzed using Wilcoxon's
#' signed rank statistic.
#'
#' @param data_outcomes dataframe whose rows contain N individual subjects and
#' whose columns have L particular hypotheses.
#' @param control_idx matrix of size nsetsx1, where entry k contains the control
#' index of data_outcomes corresponding to the kth matched pair.
#' @param treated_idx vector of length nsets, where entry k contains the treated
#' index of data_outcomes corresponding to the kth matched pair.
#' @param alt_vec character vector of size L denoting direction of one-sided 
#' alternative hypotheses.
#' @param planning_sample_prop proportion of the data to allocate to the
#' planning sample.
#' @param alpha alpha-level of test.
#' @param gamma value of the sensitivity parameter.
#' @param nboot number of bootstraps to calculate variance of each sensitivity
#' value on the planning sample.
#' @param seed value of seed for reproducibility.
#' @return A list containing the names of hypotheses the Sens-Val procedure
#' chooses to test and those it rejects.
#' @export
sv <- function(data_outcomes,control_idx,treated_idx,
               alt_vec=NULL,planning_sample_prop=0.20,
               alpha=0.05,gamma=1,nboot=250,seed=0){
  
  ## Set seed for reproducibility
  set.seed(seed)
  num_outcomes <- ncol(data_outcomes)
  nsets <- length(treated_idx)
  
  ## Create objects to return
  outcomes_tested <- outcomes_rejected <- logical(length=num_outcomes)
  
  ## Divide data into planning and analysis samples
  planning_sample_size_sets <- floor(planning_sample_prop*nsets)
  planning_sample_size_subjects <- planning_sample_size_sets*2
  ix_subjects_treated_planning <- 
    sample(x=treated_idx,size=planning_sample_size_sets)
  ix_subjects_control_planning <- 
    control_idx[match(ix_subjects_treated_planning,treated_idx),]
  ix_subjects_treated_analysis <- 
    treated_idx[-match(ix_subjects_treated_planning,treated_idx)]
  ix_subjects_control_analysis <- 
    control_idx[match(ix_subjects_treated_analysis,treated_idx),]
  
  ## Some relevant parameters
  n1 <- planning_sample_size_sets
  n2 <- nsets-planning_sample_size_sets
  
  ## Planning stage
  if (is.null(alt_vec)){
    ## Get vector for alternatives to test using planning data
    alt_vec <- character(num_outcomes)
    for (outcome in 1:num_outcomes) {
      y_plan <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n1){
        if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],
                        y[ix_subjects_treated_planning[set]])))) {
          next
        }
        y_plan <- c(y_plan,y[as.matrix(ix_subjects_control_planning)[set,]],
                    y[ix_subjects_treated_planning[set]])
      }
      mpdifs <- diff(y_plan)[seq(1,length(y_plan),2)]
      alt_vec[outcome] <- ifelse(mean(mpdifs)>=0,"greater","less")
    }
  }
  ixsetlist <- list()
  for (boot in 1:nboot){
    ixset <- sample(x=1:planning_sample_size_sets,
                    size=planning_sample_size_sets,
                    replace=T)
    ixsetlist[[boot]] <- c(as.matrix(ix_subjects_control_planning)[ixset,],
                           ix_subjects_treated_planning[ixset])
  }
  sigma_hat_vec <- kappahatvec_whole <- numeric(length = num_outcomes)
  for (outcome in 1:num_outcomes){
    ## Compute sigmahat via bootstrapping
    y <- as.matrix(data_outcomes)[,outcome]
    kappahatvec <- numeric()
    for (boot in 1:nboot){
      y_boot <- y[ixsetlist[[boot]]]
      mpdifs <- diff(y_boot)[seq(1,length(y_boot),2)]
      if (!(all(is.na(mpdifs)))){
        kappahatvec <- c(kappahatvec,
                         sensitivity_value(mpdifs[!is.na(mpdifs)],
                                           a=alpha,
                                           alt=alt_vec[outcome]))
        }
    }
    sigma_hat_vec[outcome] <- sqrt(stats::var(kappahatvec)*n1)
    ## Compute kappahat on the whole planning sample
    y_boot <- numeric()
    for (set in 1:n1){
      if (any(is.na(c(y[as.matrix(ix_subjects_control_planning)[set,]],
                      y[ix_subjects_treated_planning[set]])))) { next }
      y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_planning)[set,]],
                  y[ix_subjects_treated_planning[set]])
    }
    kappahatvec_whole[outcome] <- sensitivity_value(
      diff(y_boot)[seq(1,length(y_boot),2)], a=alpha, alt=alt_vec[outcome])
  }
  aprimegrid <- alpha/(1:num_outcomes)
  tested_list <- vector("list", num_outcomes)
  difvec <- numeric(length = num_outcomes)
  for (aprime_idx in 1:num_outcomes){
    aprime <- aprimegrid[aprime_idx]
    our_tested_ix = logical(length = num_outcomes)
    for (outcome in 1:num_outcomes){
      sigma_hat <- sigma_hat_vec[outcome]
      kappahat <- kappahatvec_whole[outcome]
      sigma_r <- sigma_hat / sqrt(planning_sample_prop * 
                                    (1 - planning_sample_prop))
      c1_hat <- -sqrt(4/3) * stats::qnorm(1 - alpha) * 
        sqrt(kappahat * (1 - kappahat))
      c2_hat <- -sqrt(4/3) * stats::qnorm(1 - aprime) * 
        sqrt(kappahat * (1 - kappahat))
      mu_r <- c1_hat / sqrt(planning_sample_prop) - 
        c2_hat / sqrt(1 - planning_sample_prop)
      our_tested_ix[outcome] <- kappahat > 
        gamma / (1 + gamma) + 
        (mu_r - stats::qnorm(1 - alpha) * sigma_r) / (sqrt(n1 + n2))
    }
    tested_list[[aprime_idx]] <- our_tested_ix
    difvec[aprime_idx] <- abs(alpha-aprime*sum(our_tested_ix, na.rm=TRUE))
  }
  aprime <- aprimegrid[which.min(difvec)]
  outcomes_tested <- tested_list[[which.min(difvec)]]
  outcomes_tested[is.na(outcomes_tested)] <- F
  
  ## Analysis stage
  for (outcome in 1:num_outcomes){
    if (outcomes_tested[outcome]){
      y_boot <- numeric()
      y <- as.matrix(data_outcomes)[,outcome]
      for (set in 1:n2){
        if (any(is.na(c(y[as.matrix(ix_subjects_control_analysis)[set,]],
                        y[ix_subjects_treated_analysis[set]])))) { next }
        y_boot <- c(y_boot,y[as.matrix(ix_subjects_control_analysis)[set,]],
                    y[ix_subjects_treated_analysis[set]])
      }
      mpdifs <- diff(y_boot)[seq(1,length(y_boot),2)]
      if (is_binary(mpdifs)){
        reject <- sensitivity_analysis_mcnemar(mpdifs=mpdifs,gamma=gamma,
                                               alt=alt_vec[outcome]) < 
          alpha/sum(outcomes_tested)
      }
      else{
        reject <- DOS2::senWilcox(d=mpdifs,gamma=gamma,alpha=alpha,
                                  alternative=alt_vec[outcome])$pval < 
          alpha/sum(outcomes_tested)
      }
      if (!is.na(reject)&reject) { outcomes_rejected[outcome] <- T }
    }
  }
  
  return(list(outcomes_tested=names(data_outcomes)[which(outcomes_tested)],
              outcomes_rejected=names(data_outcomes)[which(outcomes_rejected)]))
}

