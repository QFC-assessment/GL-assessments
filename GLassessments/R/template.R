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
    random = NULL,
    control = NULL,
    newton_steps = FALSE,
    jitter_test = FALSE,
    likeprof = FALSE,
    mcmc = FALSE) {
  ofile <- paste0(name, ".R")

  # clean environment
  olin <- c(
    "rm(list = ls())",
    "gc()"
  )

  # load RTMB
  olin <- c(
    olin, "# packages",
    "library(RTMB)"
  )

  # function
  olin <- c(
    olin, "",
    "# function",
    "f <- function(par) {",
    paste0(indent(1), "getAll(data, par)"),
    "",
    paste0(indent(1), "return(jnll)"),
    "}", ""
  )

  # optimizer
  # random effect
  if (!is.null(random)) {
    olin <- c(
      olin,
      paste0("obj <- MakeADFun(f, par, random = ", deparse(random), ")"),
      "# obj$fn()",
      "# obj$gr()"
    )
  } else {
    olin <- c(
      olin,
      paste0("obj <- MakeADFun(f, par)")
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
  if(newton_steps == TRUE) {
    olin <- c(
      olin,
      "",
      "# newton steps",
      "for (n in 1:3) {",
      "  g <- as.numeric(obj$gr(opt$par))",
      "  h <- stats::optimHess(opt$par, obj$fn, obj$gr)",
      "  new_par <- opt$par - solve(h, g)",
      "  opt <- nlminb(new_par, obj$fn, obj$gr,",
      "    control = list(eval.max = 1e4, iter.max = 1e4)",
      "  )",
      "}",
      "# opt"
    )
  }

  # jitter test
  if(jitter_test == TRUE) {
    olin <- c(
      olin,
      "",
      "# jitter test",
      "doone <- function() {",
      "  fit <- nlminb(opt$par + rnorm(length(opt$par), sd = 0.1),",
      "    obj$fn, obj$gr,",
      "    control = list(eval.max = 5e3, iter.max = 5e3)",
      "  )",
      "  c(fit$par, \"convergence\" = fit$convergence)",
      "}",
      "set.seed(123456)",
      "jit <- replicate(100, doone())",
      "boxplot(t(jit))"
    )
  }
  
  # likelihood profile
  if(likeprof == TRUE) {
    olin <- c(
      olin,
      "",
      "# likelihood profile",
      "names(obj$par)",
      "pro <- TMB:::tmbprofile(obj, name = 1)",
      "plot(pro, ylab = \"NLL\", xlab = \"\")",
      "abline(v = opt$par[1], col = \"red\", lty = 2, lwd = 2.5)",
      "confint(pro)"
    )
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
