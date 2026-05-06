#' @title indent
#' @export
indent <- function(n = 1) {
  paste(rep("  ", n), collapse = "")
}


#' @title rtmb_template
#'
#' @export
rtmb_template <- function(
    name,
    wd = NULL,
    random = NULL,
    control = NULL,
    newton_steps = FALSE,
    jitter_test = FALSE,
    likeprof = FALSE,
    mcmc = FALSE) {
  if(is.null(wd)) {
    ofile <- paste0(name, ".R")
  } else if (!is.null(wd)) {
    ofile <- file.path(wd, paste0(name, ".R"))
  }

  # clean environment
  olin <- c(
    "# clean environment",
    "rm(list = ls())",
    "gc()"
  )

  # load RTMB
  olin <- c(
    olin, "",
    "# packages",
    "library(RTMB)"
  )

  # data and parameters lists
  olin <- c(
    olin, "",
    "# data and parameters lists",
    "dat <- list()",
    "pars <- list()"
  )

  # function
  olin <- c(
    olin, "",
    "# function",
    "f <- function(dat, pars) {",
    paste0(indent(1), "getAll(dat, pars)"),
    "",
    paste0(indent(1), "return(jnll)"),
    "}"
  )

  # optimizer
  # random effect
  if (!is.null(random)) {
    olin <- c(
      olin, "",
      "# run optimizer",
      "cmb <- function(f, d) function(p) f(d, p)",
      paste0("obj <- MakeADFun(cmb(f, dat), pars, random = ", deparse(random), ")"),
      "# obj$fn()",
      "# obj$gr()"
    )
  } else {
    olin <- c(
      olin, "",
      paste0("obj <- MakeADFun(f, pars)")
    )
  }
  if (!is.null(control)) {
    control_text <- paste0(
      "control = list(",
      paste(
        sprintf(
          "%s = %s",
          names(control),
          vapply(control, toString, character(1))
        ),
        collapse = ", "
      ),
      ")"
    )
    olin <- c(
      olin,
      "opt <- nlminb(obj$par, obj$fn, obj$gr,",
      paste0(indent(1), control_text),
      ")",
      "# opt"
    )
  } else {
    olin <- c(
      olin,
      "opt <- nlminb(obj$par, obj$fn, obj$gr)"
    )
  }

  # newton steps
  if (newton_steps == TRUE) {
    olin <- c(olin, newton_steps(print_console = FALSE))
  }

  # jitter test
  if (jitter_test == TRUE) {
    olin <- c(olin, jitter_test(print_console = FALSE))
  }

  # likelihood profile
  if (likeprof == TRUE) {
    olin <- c(olin, likelihood_profile(print_console = FALSE))
  }

  # MCMC
  if (mcmc == TRUE) {
    olin <- c(olin, mcmc_stan(print_console = FALSE, random = random))
  }

  # standard errors
  olin <- c(
    olin,
    "",
    "# standard errors",
    "sd_rep <- sdreport(obj)",
    "print(summary(sd_rep))"
  )

  # outputs
  olin <- c(
    olin,
    "",
    "# outputs",
    "# pl <- as.list(sd_rep, 'Est')",
    "# plsd <- as.list(sd_rep, 'Std')",
    "# plr <- as.list(sd_rep, 'Est', report = TRUE)",
    "# plrsd <- as.list(sd_rep, 'Std', report = TRUE)"
  )
  writeLines(olin, ofile)
  }


#' @title newton_steps
#'
#' @export
newton_steps <- function(print_console = TRUE) {
  olin <- c(
    "",
    "# newton steps",
    "for (n in 1:3) {",
    paste0(indent(1), "g <- as.numeric(obj$gr(opt$par))"),
    paste0(indent(1), "h <- numDeriv::jacobian(obj$gr, opt$par)"),
    paste0(indent(1), "new_par <- opt$par - solve(h, g)"),
    paste0(indent(1), "opt <- nlminb(new_par, obj$fn, obj$gr,"),
    paste0(indent(2), "control = list(eval.max = 1e4, iter.max = 1e4)"),
    paste0(indent(1), ")"),
    "}",
    "# opt"
  )
  if (print_console == TRUE) {
    writeLines(olin)
  } else {
    return(olin)
  }
}


#' @title jitter_test
#'
#' @export
jitter_test <- function(print_console = TRUE) {
  olin <- c(
    "",
    "# jitter test",
    "doone <- function() {",
    paste0(indent(1), "fit <- nlminb(opt$par + rnorm(length(opt$par), sd = 0.1),"),
    paste0(indent(2), "obj$fn, obj$gr,"),
    paste0(indent(2), "control = list(eval.max = 5e3, iter.max = 5e3)"),
    paste0(indent(1), ")"),
    paste0(indent(1), "c(fit$par, \"convergence\" = fit$convergence)"),
    "}",
    "set.seed(123456)",
    "jit <- replicate(100, doone())",
    "boxplot(t(jit))"
  )
  if (print_console == TRUE) {
    writeLines(olin)
  } else {
    return(olin)
  }
}


#' @title likelihood_profile
#'
#' @export
likelihood_profile <- function(print_console = TRUE) {
  olin <- c(
    "",
    "# likelihood profile",
    "names(obj$par)",
    "idx <- 1",
    "pro <- TMB:::tmbprofile(obj, name = idx)",
    "plot(pro, ylab = \"NLL\", xlab = \"\")",
    "abline(v = opt$par[idx], col = \"red\", lty = 2, lwd = 2.5)",
    "confint(pro)"
  )
  if (print_console == TRUE) {
    writeLines(olin)
  } else {
    return(olin)
  }
}

#' @title mcmc_stan
#'
#' @export
mcmc_stan <- function(print_console = TRUE, random = NULL) {
  olin <- c(
    "",
    "# MCMC",
    "library(tmbstan)",
    "fitmcmc <- tmbstan(obj, chains = 1,",
    paste0(indent(1), "iter = 1e4,")
  )
  if (!is.null(random)) {
    olin <- c(
      olin,
      paste0(indent(1), "init = list(c(opt$par, sd_rep$par.random))"),
      ")",
      "mc <- extract(fitmcmc, par = names(c(opt$par, sd_rep$par.random)),",
      paste0(indent(1), "inc_warmup = TRUE, permuted = FALSE)")
    )
  } else if (is.null(random)) {
    olin <- c(
      olin,
      "  init = list(opt$par)",
      ")",
      "mc <- extract(fitmcmc, par = names(c(opt$par)),",
      paste0(indent(1), "inc_warmup = TRUE, permuted = FALSE)")
    )
  }
  if (print_console == TRUE) {
    writeLines(olin)
  } else {
    return(olin)
  }
}
