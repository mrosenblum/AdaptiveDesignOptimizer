### sa.temperature #############################################################
# Description: this is a function for viewing the default cooling schedule for
#   optim()'s implementation of simulated annealing. The behavior of simulated
#   annealing can strongly depend on the initial temperature.
#' Tempature of Cooling Schedule
#'
#' Thhis is a function for viewing the default cooling schedule for
#   optim()'s implementation of simulated annealing. The behavior of simulated
#   annealing can strongly depend on the initial temperature.
#' @param iteration
#' @param max.iterations
#' @param temperature
#' @param evals.per.temp
#'
#' @return
#' @export
#'
#' @examples
#' sa.temperature(iteration=20, max.iterations=1e4,
#' temperature=10, evals.per.temp=10)
#' plot(sa.temperature(iteration=seq(1, 1e4), max.iterations=1e4,
#' temperature=10, evals.per.temp=500)~seq(1, 1e4), type='l',
#' ylab="SA Temperature", xlab="Iteration")
sa.temperature <-
  function(iteration=1,
           max.iterations=10000,
           temperature=10,
           evals.per.temp=10) {
    temperature / log(((iteration-1) %/% evals.per.temp)*evals.per.temp + exp(1))
  }
### sa.temp example ############################################################

