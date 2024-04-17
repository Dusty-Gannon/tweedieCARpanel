#' Fit spatial Conditional Autoregressive (CAR) Model to panel data with Tweedie error distribution
#'
#'
#' @param form A formula specifying the model for the fixed effects.
#' @param speff A character string specifying the name of the spatial grouping variable (i.e., spatial random effect).
#' @param B A numeric matrix specifying the adjacency matrix (should be a matrix of 0's and 1's indicating which spatial groups are neighbors).
#' @param data A data frame containing the variables specified in the formula.
#' @param fixed_speff Logical indicating whether the spatial diffusion is fixed in time. If \code{TRUE}, this fits a model of the form
#' \deqn{y_{it} = x_{it}^T \beta + s_i + \rho \sum_{j\ne i}w_{ij}s_j + \epsilon_{it}}
#' where the spatial diffusion is constant in time.
#' @param sim Logical indicating whether to simulate data to be used with \pkg{DHARMa} to check quantile residuals.
#' @return A fitted model object containing the estimated fixed effects, random effects, dispersion parameter, Tweedie power parameter, and simulated data (if sim = TRUE).
#' @examples
#' fit_panel_CAR_tw(form = response ~ predictor1 + predictor2,
#'                   speff = "group",
#'                   B = adjacency_matrix,
#'                   data = my_data,
#'                   fixed_speff = TRUE)
#' @export
fit_panel_CAR_tw <- function(form, speff, B, data, fixed_speff = TRUE, sim = FALSE) {

  # Convert data to data frame
  data <- as.data.frame(data)
  root <- system.file(package = "tweedieCARpanel")

  # Check existence of model files
  if (!file.exists(paste0(root, "/src/tweedieCAR_fixspdiff.o")) &
      !file.exists(paste0(root, "src/tweedieCAR_fixspdiff.cpp"))) {
    stop("Looking for the model C++ code in ../src/, but could not find it.")
  }

  # Compile model code if necessary
  if (!file.exists(paste0(root, "/src/tweedieCAR_fixspdiff.o")) &
      file.exists(paste0(root, "/src/tweedieCAR_fixspdiff.cpp"))) {
    print("Compiling model code...\n")
    TMB::compile(paste0(root, "/src/tweedieCAR_fixspdiff.cpp"))
    print("Done.\n")
  }

  # Load compiled model
  dyn.load(TMB::dynlib(paste0(root, "/src/tweedieCAR_fixspdiff")))

  # Prepare data for fitting
  X <- model.matrix(lm(form, data = data))
  varnames <- all.vars(form)
  datlist <- list(
    y = data[, varnames[1]],
    X = X,
    group = factor(data[, speff]),
    B = B
  )

  # Set initial parameter values for likelihood evaluation
  parslist <- list(
    s = rep(0, length(unique(data[, speff]))),
    log_sigma_s = 0,
    logit_rho = 0,
    beta = rep(0, ncol(X)),
    log_phi = 0,
    logit_xi = 0
  )

  # Create objective function
  obj <- TMB::MakeADFun(datlist, parslist, DLL = "tweedieCAR_fixspdiff", random = "s", silent = T)

  # Optimize parameters
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  # Generate report
  report <- TMB::sdreport(obj)

  # Create matrix of simulations for DHARMA
  if (sim) {
    y_rep <- replicate(250, obj$simulate()$y)
  } else {
    y_rep <- NULL
  }

  # Construct list of results
  betaids <- grep("beta", names(report$par.fixed))
  fit <- list(
    fixef = list(
      estims = report$par.fixed[betaids],
      se = sqrt(diag(report$cov.fixed))[betaids],
      V = report$cov.fixed[betaids, betaids]
    ),
    ranef = list(
      s = report$par.random,
      sigma_s = report$value["sigma_s"],
      spatial_corr = report$value["rho"],
      adjMat = B
    ),
    dispersion = report$value["phi"],
    tweedie_power = report$value["xi"],
    y_rep = y_rep
  )

  return(fit)
}








