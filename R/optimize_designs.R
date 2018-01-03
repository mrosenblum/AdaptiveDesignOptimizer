#' Adaptive Enrichment Design Optimization Using Simulated Annealing
#' Authors: Josh Betz (jbetz@jhu.edu), Tianchen Qian, Michael Rosenblum
#'
#' @param ui.n.arms (the prefix 'ui' abbreviates 'user-input') Number of Arms (including control arm), e.g., 2 arms means one treatment arm and one control arm; the design classes that come with this package can handle 2 or 3 arms
#' @param ui.type.of.outcome.data "continuous", "binary", or "time-to-event" (i.e., survival outcome)
#' @param ui.time.to.event.trial.type "superiority" or "non-inferiority" trial type; only implemented for time-to-event outcomes
#' @param ui.time.to.event.non.inferiority.trial.margin Non-inferiority margin; only relevant if outcome is time-to-event and trial type is non-inferiority. Represented as hazard ratio, required to be at least 1.
#' @param ui.subpopulation.1.size Proportion of overall population in subpopulation 1. Must be between 0 and 1.
#' @param ui.total.alpha Familywise Type I error rate (1-sided)
#' @param ui.max.size Maximum allowed sample size
#' @param ui.max.duration Maximum allowed trial duration in years
#' @param ui.accrual.yearly.rate Number of participants enrolled per year; assumed constant throughout trial
#' @param ui.followup.length Time from enrollment to measurement of primary outcome (only used for continuous or binary outcome types)
#' @param ui.optimization.target Quantity being optimized (objective function); "size" represents expected sample size
#' @param ui.time.to.event.censoring.rate probability that primary outcome is censored, assumed to be independent of the outcome and subpopulation; only implemented for time-to-event outcomes
#' @param ui.mcid Minimum, Clinically Important Treatment Effect (as difference of population means for binary/continuous outcomes; as hazard ratio for time-to-event outcomes)
#' @param ui.incorporate.precision.gain Incorporate into analysis a precision gain from adjustment for prognostic baseline variables; allowed values: TRUE or FALSE
#' @param ui.relative.efficiency If ui.incorporate.precision.gain==TRUE, this specifies relative efficiency (number > 1), representing the assumed precision gain from adjustment for prognostic baseline variables.
#' @param ui.max.stages Maximum number of stages allowed (currently not used by default classes of trial designs)
#' @param ui.include.designs.start.subpop.1 Search over designs that allow only subpopulation 1 to be enrolled during stage 1; TRUE or FALSE
#' @param ui.population.parameters Matrix encoding scenarios (data generating distributions) used to define power constraints and  objective function
#' @param ui.desired.power Matrix encoding power requirements for each scenario
#' @param ui.scenario.weights Matrix encoding weights used to define objective function
#' @param min.n.per.arm Minimum sample size per arm allowed
#' @param min.enrollment.period Minimum enrollment duration for trial
#' @param simulated.annealing.parameter.function.scale Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.n.scale Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.period.scale Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.max.iterations Maximum Number of Different Designs to Search Over using Simluated Annealing optimization
#' @param simulated.annealing.parameter.n.simulations Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.means.temperature Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.survival.temperature Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.evals.per.temp Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.report.iteration Used by Simulated Annealing Optimization algorithm
#' @param simulated.annealing.parameter.power.penalty Used in Objective Function to incorporate Power Constraints by Simulated Annealing Optimization algorithm
#' @return 4 element list containing optimized designs from four classes (with increasing complexity):
#' @section Designs:
#' Single.Stage.Equal.Alpha.Allocation.Design
#' Single.Stage.Optimized.Alpha.Allocation.Design
#' Two.Stage.Equal.Alpha.Allocation.Design
#' Two.Stage.Optimized.Alpha.Allocation.Design
#'
#' The first two designs are not adaptive, while the second two designs are adaptive.
#'
#' Each optimized design is a list containing design.parameters and design.performance
#' @section design.parameters:
#' design.parameters contains the following values:
#' cumulative.sample.sizes.and.calendar.time.per.stage=cumulative number enrolled (if no early stopping) per stage and calendar times of analyses just after each stage. In column names, ``A'' and ``C'' denote the treatment arm and control arm, respectively; numbers 1 and 2 indicate the corresponding subpopulation. Sample sizes represent the number enrolled at the time of the corresponding analysis (which may exceed the number of participants with outcomes observed, due to the time between enrollment and outcome measurement for each participant
#'
#' alpha.allocation=Alpha allocation using Error Spending Approach
#'
#' futility.boundaries=Boundaries for stopping subpopulation accrual, on the z-scale (or in designs with more than one treatment arm compared to control, this is gives for each treatment arm by subpopulation combination.
#'
#' @section design.performance:
#' design.performance contains the following values:
#' Power=Power to reject each null hypothesis under each scenario (NA indicates null hypothesis is true, so no power is presented)
#'
#' Type.1.Error=Type I error for each null hypothesis under each scenario (NA indicates null hypothesis is false)
#'
#' Expected.Sample.Size
#'
#' Expected.Duration (in years)
#'
#' Distribution.of.sample.size.and.duration.per.scenario=For each scenario, every possible combination of early stopping is considered. Columns C1, A1, etc. have the same meaning as described about for sample.sizes.and.calendar.time.per.stage. The value listed under each such column gives the analysis number at which accrual for that arm by subpopulation combination is stopped. E.g., C1=2,C2=1,A1=2,A2=1 corresponds to stopping the control and treatment A for subpopualtion 1 at the end of stage 1, while these continue to the end of stage 2 for subpopulation 2. The subsequent columns give the sample size, duration, and person-time when this pattern occurs. The columns frequency and proportion tell how often this pattern occurred under the corresponding scenario number (based on simulation).
#' @export
#' @examples
#' #For demonstration purposes, the examples below only execute 2 iterations of simulated annealing.
#' #In general, it is recommended to use at least 500 iterations.
#' #Example 1: Time-to-event outcome; 1 treatment arm versus control; non-inferiority design
#' optimized_designs <- optimize_designs(
#'   ui.n.arms=2,
#'   ui.type.of.outcome.data="time-to-event",
#'   ui.time.to.event.trial.type="non-inferiority",
#'   ui.time.to.event.non.inferiority.trial.margin=1.35,
#'   ui.subpopulation.1.size=0.33,
#'   ui.total.alpha=0.05,
#'   ui.max.size=10000,
#'   ui.max.duration=10,
#'   ui.accrual.yearly.rate=1000,
#'   ui.followup.length=0,
#'   ui.optimization.target="size",
#'   ui.time.to.event.censoring.rate=0,
#'   ui.mcid=0.1,
#'   ui.incorporate.precision.gain=FALSE,
#'   ui.relative.efficiency=1,
#'   ui.max.stages=2,
#'   ui.include.designs.start.subpop.1=FALSE,
#'   ui.population.parameters= 0.08*matrix(c(1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.3500001,
#'     1.00, 1.00, 1.00, 2.14, 1.00, 1.00, 1.3500001, 1.3500001), ncol=4, byrow=TRUE,
#'     dimnames=list(c(),c("lambda1_con","lambda2_con","lambda1_trt","lambda2_trt"))),
#'   ui.desired.power=0.8*matrix(c(1.00, 1.00, 0, 1.00, 0, 0, 1.00, 0, 0, 0, 0, 0), ncol=3,
#'     byrow=TRUE,dimnames=list(c(),c("Pow_H(0,1)","Pow_H(0,2)","Pow_Reject_H0,1_and_H0,2"))),
#'   ui.scenario.weights=matrix(rep(0.25,4),ncol=1,dimnames=list(c(),c("weight"))),
#'   simulated.annealing.parameter.max.iterations=2,
#'   )
#'
#'   #Example 2: continuous outcome; 1 treatment arm versus control; superiority design
#' optimized_designs <- optimize_designs(
#'   ui.n.arms=2,
#'   ui.type.of.outcome.data="continuous",
#'   ui.time.to.event.trial.type="",
#'   ui.time.to.event.non.inferiority.trial.margin=NULL,
#'   ui.subpopulation.1.size=0.5,
#'   ui.total.alpha=0.05,
#'   ui.max.size=1000,
#'   ui.max.duration=5,
#'   ui.accrual.yearly.rate=250,
#'   ui.followup.length=1,
#'   ui.optimization.target="size",
#'   ui.time.to.event.censoring.rate=0,
#'   ui.mcid=NULL,
#'   ui.incorporate.precision.gain=FALSE,
#'   ui.relative.efficiency=1,
#'   ui.max.stages=5,
#'   ui.include.designs.start.subpop.1=FALSE,
#'   ui.population.parameters=matrix(c(15,15,3600,3600,3600,3600,15,0,3600,3600,3600,3600,
#'     0,15,3600,3600,3600,3600,0,0,3600,3600,3600,3600),nrow=4, ncol=6, byrow=TRUE,dimnames=
#'     list(c(),c("delta1","delta2","sigma1_trt","sigma1_con","sigma2_trt","sigma2_con"))),
#'   ui.desired.power=matrix(c(0,0,0.8,0.8,0,0,0,0.8,0,0,0,0), nrow=4, ncol=3, byrow=TRUE,
#'     dimnames=list(c(),c("Pow_H(0,1)","Pow_H(0,2)","Pow_Reject_H0,1_and_H0,2"))),
#'   ui.scenario.weights=matrix(c(0.25,0.25,0.25,0.25),ncol=1,dimnames=list(c(),c("weight"))),
#'   simulated.annealing.parameter.max.iterations=2
#' )
#'
#'  #Example 3: binary outcome; 1 treatment arm versus control; superiority design
#' optimized_designs <- optimize_designs(
#'   ui.n.arms=2,
#'   ui.type.of.outcome.data="binary",
#'   ui.time.to.event.trial.type="",
#'   ui.time.to.event.non.inferiority.trial.margin=NULL,
#'   ui.subpopulation.1.size=0.4,
#'   ui.total.alpha=0.05,
#'   ui.max.size=2000,
#'   ui.max.duration=5,
#'   ui.accrual.yearly.rate=400,
#'   ui.followup.length=0,
#'   ui.optimization.target="size",
#'   ui.time.to.event.censoring.rate=0,
#'   ui.mcid=NULL,
#'   ui.incorporate.precision.gain=TRUE,
#'   ui.relative.efficiency=1.2,
#'   ui.max.stages=5,
#'   ui.include.designs.start.subpop.1=FALSE,
#'   ui.population.parameters=matrix(c(0.4,0.3,0.5,0.4,0.4,0.3,0.4,0.4,0.3,0.3,0.4,0.4),
#'     nrow=3, ncol=4, byrow=TRUE,dimnames=list(c(),c("p1_trt","p1_con","p2_trt","p2_con"))),
#'   ui.desired.power=matrix(c(0,0,0.8,0.8,0,0,0,0,0), nrow=3, ncol=3, byrow=TRUE,
#'     dimnames=list(c(),c("Pow_H(0,1)","Pow_H(0,2)","Pow_Reject_H0,1_and_H0,2"))),
#'   ui.scenario.weights=matrix(c(0.33,0.33,0.34),ncol=1,dimnames=list(c(),c("weight"))),
#'   simulated.annealing.parameter.max.iterations=2
#' )
#' @importFrom stats plogis
#' @importFrom mvtnorm pmvnorm GenzBretz
optimize_designs <- function(
  ui.n.arms,
  ui.type.of.outcome.data,
  ui.time.to.event.trial.type,
  ui.time.to.event.non.inferiority.trial.margin,
  ui.subpopulation.1.size,
  ui.total.alpha,
  ui.max.size,
  ui.max.duration,
  ui.accrual.yearly.rate,
  ui.followup.length,
  ui.optimization.target,
  ui.time.to.event.censoring.rate,
  ui.mcid,
  ui.incorporate.precision.gain,
  ui.relative.efficiency,
  ui.max.stages,
  ui.include.designs.start.subpop.1,
  ui.population.parameters,
  ui.desired.power,
  ui.scenario.weights,
  min.n.per.arm =25,       # For Continuous/Binary Outcomes
  min.enrollment.period =0.5,    # For Survival Outcomes
  simulated.annealing.parameter.function.scale =1,
  simulated.annealing.parameter.n.scale =100,
  simulated.annealing.parameter.period.scale =2,
  simulated.annealing.parameter.max.iterations =1000,
  simulated.annealing.parameter.n.simulations =1e4,
  simulated.annealing.parameter.means.temperature =100,
  simulated.annealing.parameter.survival.temperature =10,
  simulated.annealing.parameter.evals.per.temp =10,
  simulated.annealing.parameter.report.iteration =1,
  simulated.annealing.parameter.power.penalty =100000
){
  # Get start time
  isa.start.time <- proc.time()

  ### NOTE: restricted to two subpopulations ###
  n.subpopulations <- 2
  n.arms <- ui.n.arms
  ui.subpopulation.sizes <- c(ui.subpopulation.1.size, 1-ui.subpopulation.1.size)
  # If random seed is supplied, specify seeds. Otherwise pseudorandom seeds
  # are chosen based on the initial RNG state.
  if(!exists("initial.seed")){
    initial.seed <- sample(x=1:1e8, size=1)
  }

  # Set random seed
  set.seed(initial.seed)

  # Source design.evaluation code corresponding to number of arms in trial
  if(n.arms==2){
    # Computes distribution of test statistics in a given scenario,
    # using canonical joint distribution
    construct.joint.distribution.of.test.statistics <-
      function(...){
        construct.joint.distribution.of.test.statistics.OneTreatmentArm(...)
      }
    # Computes efficacy stopping boundaries
    generate.efficacy.boundaries <-
      function(...){
        get.eff.bound.OneTreatmentArm(...)
      }
    # Evaluates performance of simulated trials
    design.evaluate <-
      function(...){
        design.evaluate.OneTreatmentArm(...)
      }
    summarize.design.parameters.and.performance <-
      function(...){
        summarize.design.parameters.and.performance.OneTreatmentArm(...)
      }
  } else if(n.arms==3){
    # Computes distribution of test statistics in a given scenario,
    # using canonical joint distribution
    construct.joint.distribution.of.test.statistics <-
      function(...){
        construct.joint.distribution.of.test.statistics.TwoTreatmentArms(...)
      }
    # Computes efficacy stopping boundaries
    generate.efficacy.boundaries <-
      function(...){
        get.eff.bound.TwoTreatmentArms(...)
      }
    # Evaluates performance of simulated trials
    design.evaluate <-
      function(...){
        design.evaluate.TwoTreatmentArms(...)
      }
  }
  # Set functions for computing design features and design evaluation
  ##
  ## Format User Inputs from Graphical User Interface
  ##

  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    if(n.arms==2){
      if(ui.type.of.outcome.data=="binary") {
        ui.outcome.mean <- subset(ui.population.parameters,select=c(2,4,1,3))
        ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
      } else{
        ui.outcome.mean <- cbind(array(0,c(nrow(ui.population.parameters),2)),subset(ui.population.parameters,select=c(1,2)))
        ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(4,6,3,5)))
      }
    } else if(n.arms==3){
      if(ui.type.of.outcome.data=="binary") {
        ui.outcome.mean <- ui.population.parameters
        ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
      } else{
        ui.outcome.mean <- subset(ui.population.parameters,select=c(1,3,5,7,9,11))
        ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(2,4,6,8,10,12)))
      }
    }
    arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
    colnames(ui.outcome.sd) <- colnames(ui.outcome.mean) <-
      paste0(rep(arm.names, each=n.subpopulations),
             rep(1:n.subpopulations, n.arms))
  } else { # Survival Cases
    ui.hazard.rate <- ui.population.parameters
    if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}
    arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
    colnames(ui.hazard.rate) <-
      paste0(rep(arm.names, each=n.subpopulations),
             rep(1:n.subpopulations, n.arms))
  }
  ## Run optimizations, starting with 1 stage
  n.stages <- 1 # Single Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    # if(ui.optimization.target=="ESS") {
    #   Switch Objective Function and Parameters
    # }
    max.enrollment.period <- (ui.max.duration-ui.followup.length)
    max.possible.accrual <- ui.accrual.yearly.rate*max.enrollment.period
    #Binary search to minimize sample size over feasible designs
    feasible.n.per.arm <- osea.result <- NULL
    n.per.arm.upper.bound <- min(ui.max.size, max.possible.accrual)/n.arms
    n.per.arm.lower.bound <- min.n.per.arm
    #Increase upper bound if necessary to meet power requirements
    #repeat{
    #  osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type=ui.type.of.outcome.data,
    #                                                      n.per.arm=floor(n.per.arm.upper.bound),
    #                                                      n.arms=n.arms,
    #                                                      accrual.rate=ui.accrual.yearly.rate,
    #                                                      delay=ui.followup.length,
    #                                                      subpopulation.sizes=ui.subpopulation.sizes,
    #                                                      interim.info.times=NULL,
    #                                                      outcome.mean=ui.outcome.mean,
    #                                                      outcome.sd=ui.outcome.sd,
    #                                                      mcid=ui.mcid,
    #                                                      futility.boundaries=NULL,
    #                                                      relative.efficiency=ui.relative.efficiency,        #                                                      n.simulations=simulated.annealing.parameter.n.simulations,
    #                                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
    #                                                                           number.of.alpha.allocation.components),
    #                                                      total.alpha=ui.total.alpha)
    #  discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
    #  feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
    #                                  discrepancy.between.desired.power.empirical.power<=0)
    #  if(feasibility.indicator){
    #    break
    #  }
    #  n.per.arm.upper.bound <- 2*n.per.arm.upper.bound
    #}
    while(n.per.arm.upper.bound-n.per.arm.lower.bound>0.1){
      candidate.n.per.arm <- mean(c(n.per.arm.lower.bound,n.per.arm.upper.bound))
      osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type=ui.type.of.outcome.data,
                                                                         n.per.arm=floor(candidate.n.per.arm),
                                                                         n.arms=n.arms,
                                                                         accrual.rate=ui.accrual.yearly.rate,
                                                                         delay=ui.followup.length,
                                                                         subpopulation.sizes=ui.subpopulation.sizes,
                                                                         interim.info.times=NULL,
                                                                         outcome.mean=ui.outcome.mean,
                                                                         outcome.sd=ui.outcome.sd,
                                                                         mcid=ui.mcid,
                                                                         futility.boundaries=NULL,
                                                                         relative.efficiency=ui.relative.efficiency,
                                                                         n.simulations=simulated.annealing.parameter.n.simulations,
                                                                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                                                              number.of.alpha.allocation.components),
                                                                         total.alpha=ui.total.alpha,
                                                                         construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                                                         generate.efficacy.boundaries=generate.efficacy.boundaries,
                                                                         design.evaluate=design.evaluate)
      discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
      feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
                                      discrepancy.between.desired.power.empirical.power<=0)
      if(feasibility.indicator==TRUE){#Current sample size is feasible; store results and explore smaller sample sizes
        #osea.result <- osea.design.performance.evaluation;
        feasible.n.per.arm <- candidate.n.per.arm;
        n.per.arm.upper.bound <- candidate.n.per.arm} else{ #Current sample size is infeasible; explore larger sample sizes
          n.per.arm.lower.bound <- candidate.n.per.arm
        }
    }
    if(is.null(feasible.n.per.arm)){feasible.n.per.arm <- min(ui.max.size, max.possible.accrual)/n.arms} #If none feasible, consider max allowed sample size
    #Placeholder to force output into format expected by .Rnw file for report building
    osea.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=floor(feasible.n.per.arm)),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x){floor(feasible.n.per.arm)}),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        interim.info.times=NULL,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                             number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.n.scale,
                  max.iterations=2,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)

  } else {#Survival Outcome
    feasible.enrollment.period <- osea.design.performance.evaluation <- NULL
    enrollment.period.upper.bound <- min(ui.max.duration,
                                         ui.max.size/ui.accrual.yearly.rate)
    enrollment.period.lower.bound <- min.enrollment.period
    feasible.max.duration <- ui.max.duration
    #repeat{
    #  osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type='survival',
    #                                                        enrollment.period=enrollment.period.upper.bound,
    #                                                        n.arms=n.arms,
    #                                                        accrual.rate=ui.accrual.yearly.rate,
    #                                                        subpopulation.sizes=ui.subpopulation.sizes,
    #                                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
    #                                                        hazard.rate=ui.hazard.rate,
    #                                                        time=max(ui.max.duration,enrollment.period.upper.bound),
    #                                                        max.follow=Inf,
    #                                                        censoring.rate=ui.time.to.event.censoring.rate,
    #                                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
    #                                                        restrict.enrollment=FALSE,
    #                                                        mcid=ui.mcid,
    #                                                        futility.boundaries=NULL,
    #                                                        relative.efficiency=ui.relative.efficiency,
    #                                                         n.simulations=simulated.annealing.parameter.n.simulations,
    #                                                        alpha.allocation=
    #                                                          rep(1/number.of.alpha.allocation.components,
    #                                                              number.of.alpha.allocation.components),
    #                                                        total.alpha=ui.total.alpha)
    #  discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
    #  feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
    #                                  discrepancy.between.desired.power.empirical.power<=0)
    #  if(feasibility.indicator){
    #    break
    #  }
    #  enrollment.period.upper.bound <- 2*enrollment.period.upper.bound
    #}
    while(enrollment.period.upper.bound-enrollment.period.lower.bound>0.01){
      candidate.enrollment.period <- mean(c(enrollment.period.lower.bound,enrollment.period.upper.bound))
      osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type='survival',
                                                                         enrollment.period=candidate.enrollment.period,
                                                                         n.arms=n.arms,
                                                                         accrual.rate=ui.accrual.yearly.rate,
                                                                         subpopulation.sizes=ui.subpopulation.sizes,
                                                                         non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                                                         hazard.rate=ui.hazard.rate,
                                                                         time=max(ui.max.duration,candidate.enrollment.period),
                                                                         max.follow=Inf,
                                                                         censoring.rate=ui.time.to.event.censoring.rate,
                                                                         ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                                                         restrict.enrollment=FALSE,
                                                                         mcid=ui.mcid,
                                                                         futility.boundaries=NULL,
                                                                         relative.efficiency=ui.relative.efficiency,
                                                                         n.simulations=simulated.annealing.parameter.n.simulations,
                                                                         alpha.allocation=
                                                                           rep(1/number.of.alpha.allocation.components,
                                                                               number.of.alpha.allocation.components),
                                                                         total.alpha=ui.total.alpha,
                                                                         construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                                                         generate.efficacy.boundaries=generate.efficacy.boundaries,
                                                                         design.evaluate=design.evaluate)
      discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
      feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
                                      discrepancy.between.desired.power.empirical.power<=0)
      if(feasibility.indicator){#Current sample size is feasible; store current results and explore smaller sample sizes
        #osea.result <- osea.design.performance.evaluation;
        feasible.enrollment.period <- candidate.enrollment.period;
        feasible.max.duration <- max(ui.max.duration,feasible.enrollment.period)
        enrollment.period.upper.bound <- candidate.enrollment.period} else{
          #Current sample size is infeasible; explore larger sample sizes
          enrollment.period.lower.bound <- candidate.enrollment.period
        }
    }
    if(is.null(feasible.enrollment.period)){feasible.enrollment.period <- min(ui.max.duration,
                                                                              ui.max.size/ui.accrual.yearly.rate)} #if no feasible solution, use maximum duration
    #Placeholder to force output into format expected by .Rnw file for report building
    osea.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=feasible.enrollment.period),
                  search.transforms=
                    list(enrollment.period=function(x){feasible.enrollment.period}),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        time=max(feasible.enrollment.period,ui.max.duration),
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=
                                          rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.period.scale,
                  max.iterations=2,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ## 1SOA 1 stage optimized alpha
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    osoa.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms)),
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        interim.info.times=NULL,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)

  } else { # Survival Cases
    osoa.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(c(feasible.max.duration,
                                                        ui.max.size/ui.accrual.yearly.rate))),
                         alpha.allocation=
                           rep(1/number.of.alpha.allocation.components,
                               number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,feasible.enrollment.period),
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        time=feasible.max.duration,
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.period.scale,rep(1,number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ## Two stage design
  n.stages <- 2 # Two Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations

  ## 2SEA 2 stage equal alpha
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    # if(ui.optimization.target=="ESS") {
    #   Switch Objective Function and Parameters
    # }
    tsea.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual)),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms))
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        interim.info.times=c(1/2,1),
                                        outcome.type=ui.type.of.outcome.data,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                             number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.n.scale,
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  } else { # Survival Cases
    if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}


    tsea.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(feasible.max.duration,
                                                      ui.max.size/ui.accrual.yearly.rate))),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,
                             min(ui.max.duration,
                                 ui.max.size/ui.accrual.yearly.rate))),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        time=c(feasible.max.duration/2,feasible.max.duration),
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=
                                          rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.period.scale,
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ## 2SOA 2 stage optimized alpha
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    tsoa.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                         interim.info.times=c(1/2,1),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms)),
                      interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate
                  ),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)

  } else { # Survival Cases
    tsoa.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(feasible.max.duration,
                                                      ui.max.size/ui.accrual.yearly.rate)),
                         time=c(feasible.max.duration/2,feasible.max.duration),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=
                           rep(1/number.of.alpha.allocation.components,
                               number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,
                             min(ui.max.duration,
                                 ui.max.size/ui.accrual.yearly.rate)),
                      time=function(t){t1 <- squash(t[1],0.01,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ####
  # Optimize 2 Stage, Group Sequential Design, for 2 arm trials
  ####
  if(n.arms==2){
    # Computes distribution of test statistics in a given scenario,
    # using canonical joint distribution
    construct.joint.distribution.of.test.statistics <-
      function(...){
        construct.joint.distribution.of.test.statistics.GroupSequential.OneTreatmentArm(...)
      }
    # Computes efficacy stopping boundaries
    generate.efficacy.boundaries <-
      function(...){
        get.eff.bound.GroupSequential.OneTreatmentArm(...)
      }
    # Evaluates performance of simulated trials
    design.evaluate <-
      function(...){
        design.evaluate.GroupSequential.OneTreatmentArm(...)
      }
    ## Optimize:
    if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
      group.sequential.tsoa.result <-
        sa.optimize(search.parameters=
                      list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                           interim.info.times=c(1/2,1),
                           futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                           alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                number.of.alpha.allocation.components)
                      ),
                    search.transforms=
                      # Cap sample size at minimum of the maximum specified size
                      # and the accrual rate x maximum allowable duration
                      list(n.per.arm=function(x)
                        ceiling(
                          squash(x,
                                 min.n.per.arm,
                                 min(ui.max.size, max.possible.accrual)/n.arms)),
                        interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                        alpha.allocation=reals.to.probability
                      ),
                    fixed.parameters=list(n.arms=n.arms,
                                          accrual.rate=ui.accrual.yearly.rate,
                                          delay=ui.followup.length,
                                          subpopulation.sizes=ui.subpopulation.sizes,
                                          outcome.type=ui.type.of.outcome.data,
                                          outcome.mean=ui.outcome.mean,
                                          outcome.sd=ui.outcome.sd,
                                          mcid=ui.mcid,
                                          relative.efficiency=ui.relative.efficiency,
                                          n.simulations=simulated.annealing.parameter.n.simulations,
                                          total.alpha=ui.total.alpha,
                                          construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                          generate.efficacy.boundaries=generate.efficacy.boundaries,
                                          design.evaluate=design.evaluate
                    ),
                    create.object=triage.based.on.outcome.type,
                    evaluate.object=power.penalized.weighted,
                    function.scale=simulated.annealing.parameter.function.scale,
                    parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                    max.iterations=simulated.annealing.parameter.max.iterations,
                    temperature=simulated.annealing.parameter.means.temperature,
                    evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                    report.iteration=simulated.annealing.parameter.report.iteration,
                    scenario.weights=ui.scenario.weights,
                    power.penalty=simulated.annealing.parameter.power.penalty,
                    power.constraints=ui.desired.power,
                    optimization.target=ui.optimization.target)

    } else { # Survival Cases
      group.sequential.tsoa.result <-
        sa.optimize(search.parameters=
                      list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                    min(feasible.max.duration,
                                                        ui.max.size/ui.accrual.yearly.rate)),
                           time=c(feasible.max.duration/2,feasible.max.duration),
                           futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                           alpha.allocation=
                             rep(1/number.of.alpha.allocation.components,
                                 number.of.alpha.allocation.components)
                      ),
                    search.transforms=
                      list(enrollment.period=function(x)
                        squash(x, min.enrollment.period,
                               min(ui.max.duration,
                                   ui.max.size/ui.accrual.yearly.rate)),
                        time=function(t){t1 <- squash(t[1],0.01,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))},
                        alpha.allocation=reals.to.probability
                      ),
                    fixed.parameters=list(n.arms=n.arms,
                                          accrual.rate=ui.accrual.yearly.rate,
                                          subpopulation.sizes=ui.subpopulation.sizes,
                                          outcome.type='survival',
                                          non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                          hazard.rate=ui.hazard.rate,
                                          max.follow=Inf,
                                          censoring.rate=ui.time.to.event.censoring.rate,
                                          ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                          restrict.enrollment=FALSE,
                                          mcid=ui.mcid,
                                          relative.efficiency=ui.relative.efficiency,
                                          n.simulations=simulated.annealing.parameter.n.simulations,
                                          total.alpha=ui.total.alpha,
                                          construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                          generate.efficacy.boundaries=generate.efficacy.boundaries,
                                          design.evaluate=design.evaluate),
                    create.object=triage.based.on.outcome.type,
                    evaluate.object=power.penalized.weighted,
                    function.scale=simulated.annealing.parameter.function.scale,
                    parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                    max.iterations=simulated.annealing.parameter.max.iterations,
                    temperature=simulated.annealing.parameter.survival.temperature,
                    evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                    report.iteration=simulated.annealing.parameter.report.iteration,
                    scenario.weights=ui.scenario.weights,
                    power.penalty=simulated.annealing.parameter.power.penalty,
                    power.constraints=ui.desired.power,
                    optimization.target=ui.optimization.target)
    }
  }

  return(lapply(list(Single.Stage.Equal.Alpha.Allocation.Design=osea.result,Single.Stage.Optimized.Alpha.Allocation.Design=osoa.result,Two.Stage.Group.Sequential.Design=group.sequential.tsoa.result,Two.Stage.Equal.Alpha.Allocation.Design=tsea.result,Two.Stage.Optimized.Alpha.Allocation.Design=tsoa.result),summarize.design.parameters.and.performance))
}
