#' @title read_table_nowarn
#' function to supress incomplete final lines warning
read_table_nowarn <- function(...) {
  try_catch <- function(expr) {
    W <- NULL
    w.handler <- function(w) { # warning handler
      if (!grepl("incomplete final line", w)) W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
  }
  lis <- try_catch(read.table(...))
  if (!is.null(lis$warning)) warning(lis$warning)

  return(lis$value)
}


#' @title read_data()
#' function to read in data files in .dat format
#' .dat file
#' first line are details
#' second line are years
#' third lines are ages
#' fourth line is data type
#'   - 1: matrix
#'   - 2: vector
read_data <- function(file) {
  type <- scan(file, skip = 1, n = 1, quiet = TRUE)
  if (type == 1) { # matrix
    header <- scan(file, skip = 2, n = 4, quiet = TRUE)
    min_y <- header[1]
    max_y <- header[2]
    min_a <- header[3]
    max_a <- header[4]

    out <- as.matrix(read_table_nowarn(file, skip = 4, header = FALSE))
    rownames(out) <- min_y:max_y
    colnames(out) <- min_a:max_a
  } else if (type == 2) { # vector
    header <- scan(file, skip = 2, n = 3, quiet = TRUE)
    min_x <- header[1]
    max_x <- header[2]

    out <- as.matrix(read_table_nowarn(file, skip = 4, header = FALSE))
    rownames(out) <- min_x:max_x
  }

  return(out)
}


#' @title setup_data()
#' function to set up data for model
setup_data <- function(data_list,
                       fishery_type,
                       maturity_at_age = NULL,
                       weight_at_age = NULL,
                       natural_mortality = NULL) {
  obs <- list()
  n_data <- length(data_list)
  for (i in 1:n_data) {
    data_dim <- dim(data_list[[i]])
    if(data_dim[2] > 1) { # matrices
      d <- data.frame(data_list[[i]], check.names = FALSE)
      obs[[i]] <- reshape(
        data = d,
        direction = "long",
        varying = c(1:ncol(d)),
        v.name = "obs",
        times = as.numeric(colnames(d)),
        timevar = "age",
        ids = as.numeric(rownames(d)),
        idvar = "year"
      )
      rownames(obs[[i]]) <- NULL
      obs[[i]]$type <- fishery_type$type[i]
      obs[[i]]$fishery <- fishery_type$fishery[i]
      obs[[i]]$logobs <- log(obs[[i]]$obs)
      obs[[i]] <- obs[[i]][, c("obs", "logobs", "year", "age", "type", "fishery")]
    } else if(data_dim[2] == 1) { # vectors
      obs[[i]] <- data.frame(data_list[[i]])
      names(obs[[i]]) <- "obs"
      obs[[i]]$year <- rownames(data_list[[i]])
      rownames(obs[[i]]) <- NULL
      obs[[i]]$age <- NA
      obs[[i]]$type <- fishery_type$type[i]
      obs[[i]]$fishery <- fishery_type$fishery[i]
      obs[[i]]$logobs <- log(obs[[i]]$obs)
    }
  }

  # combine fishery and survey datasets
  obs_all <- do.call(rbind, obs)
  for(i in 1:ncol(obs_all)) {
    obs_all[, i] <- as.numeric(obs_all[, i])
  }

  # put all into one list
  data_all <- list()
  data_all$year <- obs_all$year
  data_all$age <- obs_all$age
  data_all$min_age <- min(data_all$age)
  data_all$max_age <- max(data_all$age)
  data_all$type <- obs_all$type
  data_all$fishery <- obs_all$fishery
  data_all$obs <- obs_all$obs
  data_all$logobs <- obs_all$logobs
  data_all$logobs[is.infinite(data_all$logobs)] <- NA

  return(data_all)
}
