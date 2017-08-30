library(mvtnorm)
library(plyr)

### is.scalar ##################################################################
# Description: determine if object is scalar.
# Input:
#   x (object): object to test
# Output:
#   (logical): is object a scalar
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

### is.whole.number ############################################################
# Description: determine if numeric value is a whole number 
#   (see ?is.integer) for more information.
# Input:
#   x (numeric): numeric value
# Output:
#   (logical): is whole number to a given tolerance
is.whole.number <-
  function(x, tol=.Machine$double.eps^0.5) {
    stopifnot(is.numeric(x))
    abs(x - round(x)) < tol
  }

### in.open.interval ###########################################################
# Description: determine if numeric value is in an open interval (lb, ub)
# Input:
#   x (numeric): numeric value
# Output:
#   (logical): is whole number in the interval
in.open.interval <-
  function(x, lb=0, ub=1) {
    stopifnot(is.numeric(x), is.numeric(lb), is.numeric(ub), ub>lb)
    (x > lb) & (x < ub)
  }

### reals.to.probability #######################################################
# Description: map matrix of reals to a vector of probabilities by a logistic
#   transformation and normalization.
# Input:
#   x (numeric matrix): real values to be mapped
# Output:
#   (numeric matrix): non-negative values in (0, 1)
reals.to.probability <- function(x) plogis(x)/sum(plogis(x))

### squash #####################################################################
# Description: map matrix of reals to an interval (min, max) using a linear
#   threshold (method="linear"), a logistic function (method="logit"),
#   or other (tbd) methods.
# Input:
#   x (numeric matrix): real values to be mapped to (min, max)
# Output:
#   (numeric matrix): values in (min, max)
squash <- function(x, x.min=0, x.max=1,
                   method=c("linear", "logit") [1]) {
  stopifnot(is.numeric(x), 
            x.max>x.min)
  switch(method,
         linear=pmin(pmax(x.min, x), x.max), # Linear Interpolation 
         logit=x.min + (x.max-x.min)*plogis(x, location=(x.max+x.min)/2,
                                            scale=(x.max-x.min)/6)
         
  )
}

### get.sample.sizes.per.stage #################################################
# Description: Calculate the number of patients per arm in multi-stage clinical
#   trials with one or more treatments compared with controls in one population,
#   potentially divided into disjoint subgroups. Single subpopulations and 
#   single-stage trials are also supported.
#
# Input:
#   n.arms (numeric scalar - whole): number of total arms (L+1)
#   n.per.arm (numeric scalar - whole: number of patients per arm
#   subpopulation.sizes (numeric vector - vector of probabilities): the 
#     proportion of the population in each disjoint subgroup. Must sum to 1.
#     Using subpopulation.sizes = NULL or 1 gives a trial in a single
#     population.
#   interim.info.times (numeric vector - sorted proportions, ending in 1): The 
#     proportion of patients enrolled by each interim analysis. If no interim
#     analyses are desired, set to NULL or 1. 
#   accrual.rate (numeric scalar): number of patients enrolled per year
#   delay - numeric : time in years from randomization to observation of primary
#     outcome
#   warn.last.accrual.before.interim (logical scalar): warn if an interim 
#     analysis occurs after enrollment has stopped: such an interim analysis can
#     not reduce sample size, but may reduce costs of follow-up visits.
#
# Output: 
#   analytic.n.per.stage ([K x J(L+1)] matrix), where K is the number of stages,
#     J is the number of subpopulations, and L is the number of treatments
#     compared to the control arm. For continuous or binary outcomes, this 
#     represents the patients with primary outcome at each interim analysis; for
#     survival outcome, this represents number of patients enrolled. Each row is
#     one stage of the study, with each column indicating the number of outcomes
#     available in that (treatment x subpopulation) stratum. Columns are
#     arranged by subpopulation (S1, ..., SJ) within treatments (T0, ..., TL):
#
#     stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#     stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#       ...
#     stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#
#   total.n.per.stage ([K x J(L+1)] matrix): all enrolled patients at each
#     interim analysis, irrespective of when the primary outcome is observed.
#     stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#     stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#       ...
#     stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#
#   accrual.years (K x 1 numeric vector): time to accrue the patients in each
#     stage of a trial, assuming no interim stopping.
#
#   outcome.years (K x 1 numeric vector): time to observe the outcomes of the
#     patients accrued in each of K stages of a trial: may exceed the length
#     of enrollment, depending on enrollment rate and delay. Again, this assumes
#     no interim stopping.
#
# NOTE: When delay = 0, total.n.per.stage and analytic.n.per.stage will be
#   equal when sample sizes are integers: when not integers, the number of
#   observed outcomes is rounded down using floor() and the number of enrolled
#   participants is rounded up using ceiling() to be conservative. For survival
#   outcomes, delay should be set to 0.
get.sample.sizes.per.stage <-
  function(n.arms, n.per.arm, subpopulation.sizes,
           interim.info.times, accrual.rate, delay,
           warn.last.accrual.before.interim=TRUE){

    max.arms <- 8 # This can be increased, but not beyond 26 (=dim(letters))
    
    # Input Checks
    stopifnot(is.whole.number(n.arms) & n.arms > 0 & n.arms <= max.arms,
              is.finite(n.arms) & is.scalar(n.arms),
              is.whole.number(n.per.arm) & n.per.arm > 0,
              is.finite(n.per.arm) & is.scalar(n.per.arm),
              is.numeric(accrual.rate) & accrual.rate > 0,
              is.finite(accrual.rate) & is.scalar(accrual.rate),
              is.numeric(delay) & delay >=0,
              is.finite(delay) & is.scalar(delay),
              is.logical(warn.last.accrual.before.interim) &
                is.scalar(warn.last.accrual.before.interim))
    
    # Check interim.info.times and subpopulation.sizes for coherence
    if(!is.null(subpopulation.sizes)){
      stopifnot(is.numeric(subpopulation.sizes),
                length(subpopulation.sizes) >= 1,
                all.equal(sum(subpopulation.sizes), 1))
      if(length(subpopulation.sizes)>1) {
        stopifnot(in.open.interval(subpopulation.sizes, lb=0, ub=1))
      }
    } else {
      subpopulation.sizes <- 1
    }
    
    if(!is.null(interim.info.times)){
      # Numeric, Not empty, strictly increasing
      stopifnot(is.numeric(interim.info.times),
                length(interim.info.times) >= 1, 
                !is.unsorted(interim.info.times, strictly=TRUE))
      if(length(interim.info.times)>1) {
        # Positive, less than or equal to 1, last value == 1
        stopifnot(sum(interim.info.times<0)==0,
                  mean(interim.info.times<=1)==1,
                  tail(interim.info.times, 1)==1)
      } else stopifnot(interim.info.times==1)
    } else {
      interim.info.times <- 1
    }
    
    # Enforce n.per.arm limit
    if(length(subpopulation.sizes)>1){
      max.per.subpop <- round(head(subpopulation.sizes, -1)*n.per.arm)
      max.per.subpop <- c(max.per.subpop, n.per.arm-sum(max.per.subpop))
    } else max.per.subpop <- n.per.arm
    
    # See if a subpopulation is empty or near empty
    if(sum(round(subpopulation.sizes*n.per.arm)<1)){
      warning(paste0("Empty subgroup: Check n.per.arm (", n.per.arm,
                     ") and subpopulation.sizes (", 
                     paste(subpopulation.sizes, collapse=", "), ")."))
    }
    
    # Calculate time to accrue n.per.arm
    stages <- length(interim.info.times)
    subgroups <- length(subpopulation.sizes)
    accrual.years <- n.arms*interim.info.times*n.per.arm/accrual.rate
    outcome.years <- accrual.years + delay
    
    # Calculate total patients accrued at interim analyses
    accrued.patients <-
      floor(kronecker(outcome.years*accrual.rate/n.arms,
                        t(subpopulation.sizes)))
    
    # Cap max enrollment - can overshoot with long delay
    for(i in 1:nrow(accrued.patients)) {
      if(sum(accrued.patients[i,] > max.per.subpop)>0) {
        accrued.patients[i,] <- max.per.subpop
      }
    }
    
    # Count the number of outcomes observed at each information time
    accrued.outcomes <- floor(kronecker(interim.info.times,
                                        t(subpopulation.sizes))*n.per.arm)
    
    # Enforce n.per.arm limit
    accrued.outcomes[nrow(accrued.outcomes),] <- tail(accrued.patients, 1)
    
    # Apply names to each arm, subgroup; C for control is column 1 in each block
    arm.names <- c(LETTERS[3], LETTERS[1:max.arms][-3])[1:n.arms]
    analytic.n.per.stage <- kronecker(accrued.outcomes, t(rep(1, n.arms)))
    colnames(analytic.n.per.stage) <-
      paste0(arm.names, rep(1:subgroups, each=n.arms))
    
    total.n.per.stage <- kronecker(accrued.patients, t(rep(1, n.arms)))
    colnames(total.n.per.stage) <-
      paste0(arm.names, rep(1:subgroups, each=n.arms))

    # Optionally warn if accrual completes before interim analysis - 
    # If all patients are enrolled before an interim analysis, sample size
    # can't be reduced by efficacy/futility stopping.
    if(warn.last.accrual.before.interim &
       sum(duplicated(total.n.per.stage)>0)) {
      warning("Accrual completed before interim analyses.")
    }
    
    if(sum(analytic.n.per.stage==0)>0){
      warning("Some stage did not accrue a patient in one subpopulation.")
    }
    
    return(list(analytic.n.per.stage=analytic.n.per.stage,
                total.n.per.stage=total.n.per.stage,
                accrual.years=accrual.years,
                outcome.years=outcome.years))
  }

### get.sample.sizes.per.stage example #########################################
# get.sample.sizes.per.stage(n.arms=1, n.per.arm=100,
#                            subpopulation.sizes=NULL, interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=1, n.per.arm=100,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=2, n.per.arm=100,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=3, n.per.arm=500,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=c(.25, .5, .75, 1),
#                            accrual.rate=100, delay=0.1)


### list.potential.trials ######################################################
# Description: Create a look-up table that gives the sample size and duration
# of each potential outcome of a clinical trial. This can be used to calculate
# ESS and duration from a matrix of simulated trial results. If each row is a
# simulated trial, and each column specifies the stage at which each (arm x 
# subgroup) strata ended the trial, simply tabulate the number of outcomes of
# each type, merge with the look-up table, and then use weighted.mean() on
# sample size and duration.

list.potential.trials <-
  function(total.n.per.stage,
           accrual.rate,
           subpopulation.sizes){
    
    # Get number of arms, subgroups, and stages from total.n.per.stage
    n.arms <- 
      length(unique(substring(colnames(total.n.per.stage), first=1, last=1)))
    subgroups <- length(unique(substring(colnames(total.n.per.stage), first=2)))
    stages <- nrow(total.n.per.stage)
    
    # Get the increment in sample size at each stage
    incremental.n <- diff(rbind(0, total.n.per.stage))
    
    # Create a list of columns for each subgroup;
    # Vector of columns for each control arm
    subgroup.col <- list()
    for(j in 1:subgroups) subgroup.col[[j]] <-
      which(substring(colnames(total.n.per.stage), first=2)==j)
    control.col <- match(paste0("C", 1:subgroups), colnames(total.n.per.stage))
    
    ### Calculate possible sample sizes, durations
    # Create all potential outcomes: Arm x Subopulation x Stage combinations
    all.trials <- matrix(1:stages, nrow=stages, ncol=subgroups*n.arms)
    colnames(all.trials) <- colnames(total.n.per.stage)
    all.trials <- expand.grid(data.frame(all.trials))
    
    if(nrow(all.trials)>1){
      
      # Restrict such that active treatments can't enroll without controls:
      # Keep only trials where the last stage reached by any treatment is equal
      # to the last stage reached by control
      for(j in 1:subgroups){
        tx.cols <- setdiff(subgroup.col[[j]], control.col[j])
        all.trials <-
          all.trials[which(apply(all.trials[,tx.cols], 1, max) == 
                             all.trials[,control.col[j]]),]
      }
      
      stage.duration <- matrix(0, nrow=nrow(all.trials), ncol=subgroups)
      trial.duration <- matrix(0, nrow=nrow(all.trials))
      
      for (k in 1:stages){
        patients.stage.k <- 
          (all.trials >= k)* # If Treatment x Subpopulation strata in stage k
          kronecker(matrix(1, nrow=nrow(all.trials)), 
                    t(incremental.n[k,])) # Incremental sample size by strata
        
        # Count arms enrolling; Get accrual rate in subpopulation by adjusting
        # subpopulation accrual rate by arms enrolling in current stage
        for (j in 1:subgroups){
          arms.enrolling <- rowSums((all.trials >= k)[, subgroup.col[[j]]])
          randomization.rate <- 
            (accrual.rate)*subpopulation.sizes[j]*n.arms/arms.enrolling
          stage.duration[,j] <- 
            rowSums(patients.stage.k[, subgroup.col[[j]]])/randomization.rate
        }
        trial.duration <- trial.duration + apply(stage.duration, 1, max)
        
      }
      
      # Calculate sample size of each trial
      potential.trials <-
        data.frame(all.trials,
                   sample.size=apply(t(all.trials), 2, function(x) 
                     sum(total.n.per.stage[cbind(x, 1:ncol(all.trials))])),
                   duration=trial.duration)

    } else {
      potential.trials <- data.frame(all.trials,
                                     sample.size=all.trials,
                                     duration=all.trials/accrual.rate)
    }
    
    return(potential.trials)
}

per.stage.sample.sizes <-
  get.sample.sizes.per.stage(n.arms=3, n.per.arm=500,
                             subpopulation.sizes=c(0.1, 0.9),
                             interim.info.times=c(.5, .75, 1),
                             accrual.rate=100, delay=0.1)
# NOTE: this does not match because total.n.per.stage is rounded, resulting in
# over-estimation of duration.
list.potential.trials(total.n.per.stage=
                        per.stage.sample.sizes$total.n.per.stage,
                      accrual.rate=100,
                      subpopulation.sizes=c(0.1, 0.9))

### calculate.trial.criteria ###################################################
# Description: determine the maximum sample size, expected sample size, expected
#   duration, power, type I error, and familywise type I error rate (FWER) for
#   a clinical trial based on a trial.simulation object, a matrix of hypothesis
#   rejections (decisions.rejections), a matrix of stages at which decisions
#   were made (decisions.stages), the columns corresponding to null hypotheses
#   (null.hypotheses), a look-up table of potential trials with their sample
#   sizes and durations calculated (potential.trials)
#
#
#   if no hypotheses are null, specify null.hypotheses=NULL
calculate.trial.criteria <-
  function(decisions.rejections,
           decisions.stages,
           null.hypotheses=NULL,
           potential.trials){
    
    # Add warning if names of decisions.stages don't match potential.trials
    
    # Get number of arms, subgroups, and stages from total.n.per.stage
    n.arms <- 
      length(unique(substring(colnames(decisions.stages), first=1, last=1)))
    subgroups <- length(unique(substring(colnames(decisions.stages), first=2)))
    n.simulations <- nrow(decisions.stages)
    
    # Create a list of columns for each subgroup;
    # Vector of columns for each control arm
    subgroup.col <- list()
    for(j in 1:subgroups) subgroup.col[[j]] <-
      which(substring(colnames(decisions.stages), first=2)==j)
    
    conj.power.subgroups <-
      matrix(NA, nrow=n.simulations, ncol=subgroups)
    disj.power.subgroups <-
      matrix(NA, nrow=n.simulations, ncol=subgroups)
    
    # Copy decisions.rejections
    empirical.alpha <- empirical.power <- decisions.rejections
    if(!is.null(null.hypotheses)){
      empirical.power[, null.hypotheses] <- NA
      empirical.alpha[, -null.hypotheses] <- NA
    }
    
    for(j in 1:subgroups){
      non.null.sub.cols <- setdiff(subgroup.col[[j]], null.hypotheses)
      if(length(non.null.sub.cols)>0){
        sub.cols <- as.matrix(empirical.power[,non.null.sub.cols])
        conj.power.subgroups[,j] <- mean(rowMeans(sub.cols)==1)
        disj.power.subgroups[,j] <- mean(rowMeans(sub.cols)>0)
      }
    }
    
    
    conj.power.subgroups <- colMeans(conj.power.subgroups)
    disj.power.subgroups <- colMeans(disj.power.subgroups)
    names(disj.power.subgroups) <- paste0("DIS_Power_", 1:subgroups)
    names(conj.power.subgroups) <- paste0("CON_Power_", 1:subgroups)

    conj.power <- mean(rowMeans(empirical.power, na.rm=TRUE)==1)
    disj.power <- mean(rowMeans(empirical.power, na.rm=TRUE)>0)
    empirical.power <- colMeans(empirical.power==1)
    
    if(!is.null(null.hypotheses)){
      fwer <- mean(rowMeans(empirical.alpha, na.rm=TRUE)>0)
      type.1.error <- colMeans(empirical.alpha==1)
    } else {
      fwer <- NA
      type.1.error <- head(decisions.stages, 1)*NA
    }
    
    distribution.of.trials <- count(decisions.stages)
    names(distribution.of.trials) <- c(colnames(decisions.stages), "frequency")
    distribution.of.trials$proportion <- 
      distribution.of.trials$frequency/n.simulations
    
    distribution.of.trials <- merge(potential.trials,
                                    distribution.of.trials)
    distribution.of.trials <-
      distribution.of.trials[union(names(potential.trials),
                                   names(distribution.of.trials))]
    
    return(list(distribution.of.trials=distribution.of.trials,
                power=empirical.power,
                conj.power=conj.power,
                disj.power=disj.power,
                conj.power.subgroups=conj.power.subgroups,
                disj.power.subgroups=disj.power.subgroups,
                type.1.error=type.1.error,
                fwer=fwer))
  }
                   
### linear.threshold.penalty ###################################################
# Description: determine if a matrix of values meets or exceeds a matrix of
#   thresholds. A linear penalty is applied to values below the thresholds.
#   A loss matrix is returned of the same dimension as the power matrix and 
#   thresholds.
#
# If any scenario x hypothesis combinations do not contribute to the penalty,
#   they should be replaced with NA values before using power.penalty().
#
# Input:
#   x (numeric matrix): a (RxC) matrix of values
#   threshold.matrix (numeric matrix): a (RxC) matrix of thresholds
#   linear.penalty (numeric scalar - positive): a positive linear penalty to
#     scale the absolute difference between x and threshold.matrix
#
# Output: 
#   loss (numeric matrix): a (MxH) matrix of loss values.
linear.threshold.penalty <-
  function(x,
           threshold.matrix,
           linear.penalty) {
    # Check to see if x are probability scale and equal dimensions to
    # power thresholds, and thresholds are feasible.
    stopifnot(is.numeric(linear.penalty) & is.scalar(linear.penalty),
              linear.penalty>=0,
              identical(dim(x), dim(threshold.matrix)))
    
    return(linear.penalty*(x<threshold.matrix)*abs(x-threshold.matrix))
  }
                   
### power.penalized.weighted.ess ###############################################
# Description: an objective function that encorporates a weighted average of
#   expected sample size (ESS) across outcome scenarios, with an added penalty
#   for each violation of the power constraints. 
#
# Input:
#   trial.performance (trial performance object): trial performance to evaluate.
#     This should have the following named elements: mss (numeric scalar),
#     ess (numeric vector), duration (numeric vector), power (numeric matrix),
#     fwer (numeric vector), and type.1.error (numeric matrix).
#   scenario.weights (numeric vector): a numeric vector of weights for
#     calculating a weighted mean of ESS over outcome scenarios. If no weight is
#     supplied, the weights are assumed equal.
#   power.constraints (numeric matrix): a numeric matrix of power constraints
#     matching the dimensions of trial.performance$power. If no power constraint
#     is supplied, 80% power is taken as the constraint.
#   power.penalty (numeric scalar): a penalty applied to the difference between
#     power achieved and power desired.
#
#
# Output: an objective function value.
power.penalized.weighted.ess <- 
  function(trial.performance,
           scenario.weights=NULL,
           power.constraints=NULL,
           power.penalty=10,
           objective.scale=1) {
    
    stopifnot("ess" %in% names(trial.performance),
              "power" %in% names(trial.performance),
              is.numeric(objective.scale) & is.scalar(objective.scale),
              is.finite(objective.scale) & objective.scale > 0)
    
    trial.power <- as.matrix(trial.performance$power)
    
    if(is.null(power.constraints)){
      power.constraints <-
        matrix(0.8, nrow=nrow(trial.power), ncol=ncol(trial.power))
    }
    
    if(is.null(scenario.weights)){
      scenario.weights <- rep(1, length(trial.performance$ess))
    }
    
    objective.value <- 
      weighted.mean(trial.performance$ess, w=scenario.weights) +
      sum(linear.threshold.penalty(trial.performance$power,
                               threshold.matrix=power.constraints,
                               linear.penalty=power.penalty), na.rm=TRUE)
    
    stopifnot(is.scalar(objective.value))
    return(objective.value/objective.scale)
  }

# power.penalized.weighted.ess example #########################################
# # Trial performance - 6 scenarios: 1) global null 2-6) all hypotheses non-null
# trial.performance <-
#   list(mss=100,
#        ess=seq(50, 75, 5),
#        duration=seq(5, 10, 1),
#        power=
#          rbind(NA,
#                matrix(runif(4*5, min=0.801, max=0.802), nrow=5)),
#        fwer=c(0.025, rep(NA, 5)),
#        type.1.error=
#          rbind(rep(0.025/4, 4),
#                matrix(NA, nrow=5, ncol=4)))
# # No constraints violated - Just average of ESS
# mean(trial.performance$ess)
# power.penalized.weighted.ess(trial.performance)
# weighted.mean(trial.performance$ess, w=c(10, rep(1, 5)))
# power.penalized.weighted.ess(trial.performance,
#                              scenario.weights=c(10, rep(1, 5)))
# # 20 constraints violated by ~0.10 x penalty=10 => Increase objective by ~20
# power.penalized.weighted.ess(trial.performance,
#                              power.constraints=matrix(0.90, nrow=6, ncol=4),
#                              power.penalty=10)


### transform.parameters #######################################################
# Description: transform parameters to restricted parameters spaces, such as
#   non-negative reals, positive reals, probabilities, normalized vectors, etc.
#
#
# Input:
#   parameters (list): a named list of parameters which are assumed to take any
#     real value.
#
#   transforms (list): a named list of functions used to transform the 
#     parameters from real values to another set. This can be used to handle
#     parameters that are strictly positive, require normalization, and so on,
#     such as alpha (re)allocation and information times. The names of the list
#     should correspond to the parameters intended for the transform.
#
# Output:
#   tf.params (list): a named list of transformed parameters
#
transform.parameters <- function(parameters,
                                 transforms) {
  
  # 1. Make sure parameters and transforms are named
  if(sum(nchar(names(parameters))==0) > 0) {
    stop("parameters must be named to avoid ambiguity")
  }
  
  if(sum(nchar(names(transforms))==0) > 0) {
    stop("transforms must be named to avoid ambiguity")
  }
  
  if(!is.null(transforms)){
    
    # 2. Make sure transforms map uniquely to parameters
    which.tf <- match(names(transforms), names(parameters))
    
    if(length(unique(which.tf)) != length(transforms)){
      if(sum(duplicated(which.tf))>0) {
        stop("transforms did not uniquely match to parameters")
      }
      if(sum(is.na(which.tf))>0){
        stop(paste0("parameters did not contain: ",
                    names(transforms)[which(is.na(which.tf))]))
      }
    }
  }
  
  tf.params <- parameters
  for(i in 1:length(which.tf)) {
    tf.params[[which.tf[i]]] <- 
      transforms[[i]](parameters[[which.tf[i]]])
  }
  
  return(tf.params)
}

### transform.parameters example ###############################################
# search.parameters <- list(n.per.arm=100,
#                        total.alpha=0.025,
#                        alpha.allocation=c(0, 0, 0, 0))
# search.transforms <- list(alpha.allocation=reals.to.probability)
# 
# transform.parameters(search.parameters,
#                      search.transforms)


### get.optim.trajectory ######################################################
# Description: by default, optim() only returns the value of the objective
# function at convergence. By passing options to optim, the iteration history
# can be directed to console, and capture.output() can be used to turn this
# into a character vector. This function takes in such a character vector,
# extracts the iteration history, and places it in a data frame.
#
# Input: the output of a call to capture.output(optim(...))
#
# Output: a data frame with iteration (the iteration number) and obj.fx
# (the objective function at that iteration).

get.optim.trajectory <- function(optim.text) {
  optim.text <- optim.text[grep("iter\\s*\\d*\\s*value\\s*", optim.text)]
  optim.text <- data.frame(iteration=NA, obj.fx=NA, optim.text=optim.text)
  optim.text$iteration <-
    as.numeric(gsub("(iter)|(value.*)", "", optim.text$optim.text))
  optim.text$obj.fx <- 
    as.numeric(gsub("iter\\s*\\d*\\s*value", "", optim.text$optim.text))
  optim.text[c("iteration", "obj.fx")]
}

### sa.optimize ################################################################
# Description: this is a wrapper for optimizing objects, which requires two
#   functions: a function from the parameters to create an object, and a 
#   function that assigns a value to the object (a loss function that takes
#   in an object and produces a numeric scalar value). This allows optimization
#   of objects using simulated annealing by attempting to minimize the loss for
#   a given set of parameters.
#
#   create.object() is a function that takes the lists search.parameters and 
#   fixed.parameters, and produces an object. Optional transforms specified
#   by search.transforms can constrain this parameter space. The function 
#   evaluate.object() maps the object to a numeric scalar, using the optional
#   parameters specified by ellipsis (...).
#
#   search.parameters = {p1, p2, ..., pd}
#   search.transforms = {f1(), f2(), ... , fd()}
#   transformed.parameters := {x1=f1(p1), x2=f2(p2), ..., xd=fd(pd)
#   object.to.optimize <- create.object(transformed.parameters)
#   objective.value = evaluate.object(object.to.optimize, ...)
#
#   If G is create.object and F is evaluate.object, SA will try to find: 
#     xmin = argmin(F(G(x1, x2, ..., xd)))
#
# Input:
#   search.parameters (list): a named list of all parameters that are to be
#     searched over by simulated annealing, with starting values: passed to
#     create.object()
#
#   search.transforms (list): a named list of functions used to transform the
#     search parameters from real values to another set. This can be used to
#     handle parameters that are strictly positive, require normalization, and
#     so on. If NULL, parameters are assumed to be able to assume any real
#     value. See the documentation of transform.parameters() for more
#     information.
#
#   fixed.parameters (list): a named list of all parameters that are required to
#     fully specify the creation of the object, along with their values: passed
#     to create.object()
#
#   create.object (function): a function which accepts the arguments 
#     search.parameters (after transformation by search.transforms) and
#     fixed.parameters, and returns an object.
#
#   evaluate.object (function): a function that takes an object and assigns it
#     a scalar numeric value, with smaller values being preferred.
#
#   max.iterations (positive numeric scalar): number of iterations of simulated
#     annealing. See ?optim() for more details.
#
#   temperature (positive numeric scalar): 'temperature' parameter for simulated
#     annealing. See ?optim() for more details.
#
#   evals.per.temp (positive numeric scalar): number of evaluations per 
#     temperature. See ?optim() for more details.
#
#   report.iteration (positive numeric scalar): interval for reporting objective
#     function - report every x iterations. See ?optim() for more details.
#
# Output:
# parameters: the list of parameters for the optimized object, after any 
#   transformations were applied.
# optimization.trajectory: a data frame of the objective function values by
#   iteration. See get.optim.trajectory().
# optim.result: the object returned by the call to optim().
sa.optimize <-
  function(search.parameters,
           search.transforms=NULL,
           fixed.parameters=NULL,
           create.object,
           evaluate.object,
           max.iterations,
           temperature,
           evals.per.temp,
           report.iteration,
           function.scale=1,
           parameter.scale=1,
           ...){

    
    # optim() can only take vectors: need to unlist/relist to retain structure.
    initial.param <- as.relistable(search.parameters)
    
    optim.trajectory <-
      capture.output(
        sa.result <- 
          optim(par=unlist(initial.param),
                fn=function(varying,
                            transforms=search.transforms,
                            fixed=fixed.parameters,
                            objective.fx=evaluate.object){ 
                  # Re-structure parameters as list
                  varying <- relist(varying, skeleton=search.parameters)
                  
                  # Call user-specified create.object function: capture any text
                  discard.text <-
                    capture.output(
                      sa.object <-
                        do.call(what=create.object,
                                args=c(transform.parameters(varying,
                                                            transforms),
                                       fixed))
                    )
                      
                  # Call user-specified objective function: capture any text
                  discard.text <-
                    capture.output(
                      obj.fx <- do.call(what=evaluate.object,
                                        args=list(object=sa.object, ...))
                    )
                  
                  obj.fx
                },
                method="SANN",
                control=list(trace=6,
                             maxit=max.iterations,
                             temp=temperature,
                             tmax=evals.per.temp,
                             REPORT=report.iteration,
                             fnscale=function.scale,
                             parscale=parameter.scale))
      )
        
    # Re-structure result as list
    search.parameters <- relist(sa.result$par, skeleton=search.parameters)
    
    # Transform results and return
    parameters <- 
      c(transform.parameters(search.parameters, search.transforms),
        fixed.parameters)

    return(list(parameters=parameters,
                optimization.trajectory=get.optim.trajectory(optim.trajectory),
                optim.result=sa.result))
  }

### sa.optimize example ########################################################
do.nothing <- function(...) {}
my.result <-
  sa.optimize(search.parameters = list(mean=10),
              search.transforms = list(mean=function(x) pmax(abs(x), 1e-4)),
              fixed.parameters = list(n=1, sd=1),
              create.object = rnorm,
              evaluate.object = function(object, ...) {
                do.nothing(...)
                sum(abs(object))
              },
              max.iterations=10000,
              temperature=8,
              evals.per.temp=50,
              report.iteration=1,
              par.1=data.frame(A=1:4, B=letters[1:4]),
              par.2=list("something", "something else"))
my.result$optim.result
my.result$optimization.trajectory
my.result$parameters


### sa.temperature #############################################################
# Description: this is a function for viewing the default cooling schedule for
#   optim()'s implementation of simulated annealing. The behavior of simulated
#   annealing can strongly depend on the initial temperature.
sa.temperature <- 
  function(iteration=1,
           max.iterations=10000,
           temperature=10,
           evals.per.temp=10) {
  temperature / log(((iteration-1) %/% evals.per.temp)*evals.per.temp + exp(1))
}
### sa.temp example ############################################################
# sa.temperature(iteration=20, max.iterations=1e4,
#                temperature=10, evals.per.temp=10)
# plot(sa.temperature(iteration=seq(1, 1e4), max.iterations=1e4,
#                     temperature=10, evals.per.temp=500)~seq(1, 1e4), type='l',
#      ylab="SA Temperature", xlab="Iteration")
