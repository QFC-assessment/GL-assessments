#' @title .melt_df
#' @export
.melt_df <- function(
    df,
    years,
    ages,
    type = "default") {
  df <- as.data.frame(df)
  if (type == "jitter") {
    melt_df <- data.frame(
      variables = rep(row.names(df), each = ncol(df)),
      iterations = rep(1:nrow(df), nrow(df)),
      values = as.vector(t(df))
    )
  } else if (type == "default") {
    melt_df <- reshape(
      data = df, # data frame to manipulate
      direction = "long", # how to manipulate (switch to long form)
      # Columns you are combining (melting together)
      varying = c(1:ncol(df)), # ages to combine (melt)
      v.name = "v1", # name of new combined column
      # The new subsetting column
      times = ages, # values for the grouping column
      timevar = "age", # name of the groups
      # The new individual column
      ids = years, # values for the individuals
      idvar = "year", # name of the individuals ages
      # Makes the data frame easier to read
      new.row.names = 1:(ncol(df) * nrow(df))
    )
  }
  return(melt_df)
}


#' @title age_year_plot
#' @export
age_year_plot <- function(
    df,
    years,
    ages,
    title = NULL,
    y_lab = NULL) {
  # save and restore par on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # change margins
  par(mar = c(4, 4, 3, 1), oma = c(0.5, 0.5, 0.5, 6))

  melt_df <- .melt_df(as.data.frame(df), years = years, ages = ages)

  # sort data by year to ensure consistent coloring
  years_unique <- sort(unique(melt_df$year))

  # define color scale from orange to blue
  ncols <- length(years_unique)
  cols <- colorRampPalette(c("orange", "blue"))(ncols)
  year_cols <- setNames(cols, years_unique)

  plot(
    x = range(melt_df$age),
    y = range(as.numeric(melt_df$v1), na.rm = TRUE),
    type = "n",
    xlab = "Age",
    ylab = if (is.null(y_lab)) "" else y_lab,
    main = if (is.null(title)) "" else title
  )

  # draw lines for each year
  for (yr in years_unique) {
    sub <- melt_df[melt_df$year == yr, ]
    lines(
      sub$age,
      as.numeric(sub$v1),
      col = year_cols[as.character(yr)],
      lty = 1
    )
  }

  grid()

  # add legend
  usr <- par("usr")
  par(xpd = NA)
  x_left <- usr[2] + 1.2 # push to the right of plot
  x_right <- usr[2] + 1.8
  plot_height <- usr[4] - usr[3]
  y_bottom <- usr[3] + 0.2 * plot_height # start 10% above bottom
  y_top <- usr[4] - 0.2 * plot_height # end 10% below
  # draw gradient strip
  rects <- seq(y_bottom, y_top, length.out = ncols + 1)
  for (i in 1:ncols) {
    rect(x_left, rects[i], x_right, rects[i + 1], col = cols[i], border = NA)
  }
  # add axis labels (first and last year)
  axis(
    side = 4, at = c(y_bottom - 0.03, y_top + 0.03),
    labels = c(min(years_unique), max(years_unique)),
    las = 1, tick = FALSE
  )
}


#' @title year_plot
#' @export
year_plot <- function(
    vec,
    years,
    y_min,
    y_max,
    title = NULL,
    y_lab = NULL) {
  # Setup the plot
  plot(
    x = years,
    y = vec,
    type = "l",
    xlab = "Years",
    ylab = if (is.null(y_lab)) "" else y_lab,
    main = if (is.null(title)) "" else title,
    ylim = c(y_min, y_max),
    col = "black",
    lwd = 2
  )

  grid()
}


#' @title year_age_b_plot
#' @export
year_age_b_plot <- function(
    df,
    years,
    ages,
    title = NULL,
    y_lab = NULL,
    legend = FALSE,
    legend_title = NULL) {
  # save and restore par on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # change margins
  par(mar = c(4, 4, 3, 1), oma = c(0.5, 0.5, 0.5, 6))

  plot(
    x = range(years),
    y = range(df, na.rm = TRUE),
    type = "n",
    xlab = "Years",
    ylab = if (is.null(y_lab)) "" else y_lab,
    main = if (is.null(title)) "" else title
  )

  # choose a color for this series
  cols <- seq_len(ncol(df))

  # Loop through columns
  for (i in 1:ncol(df)) {
    yvals <- df[, i]
    lines(years, yvals, col = cols[i], type = "b", pch = "")
    text(years, yvals, labels = ages[i], col = cols[i])
  }

  grid()

  # Add legend in right margin
  par(xpd = NA)
  x_inner <- grconvertX(1, from = "nic", to = "user") # right edge of the inner region (inside oma)
  x_dev <- grconvertX(1, from = "ndc", to = "user") # right edge of the device
  x_leg <- (x_inner + x_dev) / 2 # horizontal center of the right outer margin
  y_leg <- grconvertY(0.5, from = "npc", to = "user") # vertical center of the plot
  legend(
    x_leg, y_leg,
    xjust = 0.5, yjust = 0.5,
    legend = ages,
    col = cols,
    lty = 1,
    title = if (is.null(legend_title)) "" else legend_title,
    bty = "n" # no box
  )
}


#' @title bubble_plot
#' @export
bubble_plot <- function(
    pa,
    years,
    ages,
    bubble_size = 3,
    mlab = NULL) {
  n_year <- length(years)
  n_age <- length(ages)

  # Ensure pa is a data frame
  pa <- as.data.frame(pa, check.names = FALSE)
  # Use .melt_df to get long format: age, year, obs
  comp_mat <- .melt_df(pa, years = years, ages = ages)

  x <- as.numeric(as.character(comp_mat$year))
  y <- as.numeric(as.character(comp_mat$age))
  z <- as.numeric(as.character(comp_mat$v1))

  # settings
  n <- length(x)
  col_res <- rep("grey40", n)
  bg_open <- gray(0.95, 0.3)
  cex_res <- sqrt(abs(z)) * bubble_size
  pch_res <- rep(16, n)
  bg_res <- rep(bg_open, n)
  bg_res[z < 0] <- "grey80"

  xlim_res <- range(x)
  ylim_res <- range(y)

  # plot
  plot(x, y,
    type = "n",
    xlim = xlim_res, ylim = ylim_res, xlab = "", ylab = "",
    axes = FALSE, main = mlab
  )
  axis(1, at = years, labels = years)
  axis(2, at = seq(ages[1], ages[n_age], 1))
  box()
  mtext("", side = 1, line = 2.8)
  mtext("Ages", side = 2, line = 2.8)
  points(x, y, pch = pch_res, col = "grey30", cex = cex_res, bg = bg_res)

  cohort_years <- (years[1] - n_age):years[n_year]
  for (i in 1:length(cohort_years)) {
    abline(-cohort_years[i], 1, lty = 2, col = "grey80")
  }
}


#' @title effort_F_plot
#' @export
effort_F_plot <- function(
    E_t,
    F_ta,
    vul_age,
    x_min,
    x_max,
    y_min,
    y_max,
    years,
    title = NULL,
    x_lab = NULL,
    y_lab = NULL) {
  # Compute average F over vulnerable ages
  F_t <- rowMeans(F_ta[, vul_age:ncol(F_ta)])

  # Prepare data
  ef_df <- data.frame(
    eff = E_t,
    F = F_t,
    label = as.character(years)
  )

  # Set up empty plot
  plot(ef_df$eff, ef_df$F,
    type = "n", # no points, just set axes
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    xlab = if (is.null(x_lab)) "" else x_lab,
    ylab = if (is.null(y_lab)) "" else y_lab,
    main = if (is.null(title)) "" else title
  )
  # Add text labels at each point
  text(ef_df$eff, ef_df$F, labels = ef_df$label)
  abline(a = 0, b = 1, lty = 2, col = "grey50")  
  
  grid()
}


#' @title uncertainty_plot
#' @export
uncertainty_plot <- function(
    vec,
    vec_sd,
    years,
    quantile = 0.975,
    y_lab = NULL,
    title = NULL) {
  # Compute confidence bounds
  z_score <- qnorm(quantile)
  vec_sd_low <- vec - (vec_sd * z_score)
  vec_sd_high <- vec + (vec_sd * z_score)

  # Set up empty plot
  plot(years, vec,
    type = "n",
    xlab = "Years",
    ylab = if (is.null(y_lab)) "" else y_lab,
    main = if (is.null(title)) "" else title,
    ylim = range(c(vec_sd_low, vec_sd_high, vec), na.rm = TRUE)
  )

  # Draw shaded ribbon (polygon)
  polygon(
    x = c(years, rev(years)),
    y = c(vec_sd_low, rev(vec_sd_high)),
    col = rgb(1, 0, 0, alpha = 0.15), # red with transparency
    border = NA
  )

  # Draw main line
  lines(years, vec, col = "black", lwd = 2)

  # Draw points
  points(years, vec, pch = 16, col = "black")

  # Draw dashed lines for lower and upper bounds
  lines(years, vec_sd_low, col = "red", lty = 2)
  lines(years, vec_sd_high, col = "red", lty = 2)

  grid()
}


#' @title fit_vs_data_year
#' @export
fit_vs_data_year <- function(
    data,
    fit,
    years,
    y_min,
    y_max,
    title = NULL,
    y_lab = NULL) {
  # Set up empty plot
  plot(years, fit,
    type = "n",
    xlab = "Years",
    ylab = if (is.null(y_lab)) "" else y_lab,
    ylim = c(y_min, y_max),
    main = if (is.null(title)) "" else title
  )

  # Add fit line
  lines(years, fit, col = "blue", lty = 1, lwd = 2)

  # Add data points
  points(years, data, pch = 16, col = "black")

  grid()

  # Add legend
  legend("topright",
    legend = c("Data", "Fit"),
    col = c("black", "blue"),
    pch = c(16, NA),
    lty = c(NA, 1),
    lwd = c(NA, 2),
    bty = "n"
  )
}


#' @title fit_vs_data_pa
#' @export
fit_vs_data_pa <- function(
    data,
    fit,
    sample_size = NULL,
    years,
    ages,
    title = NULL) {
  # save and restore par on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Convert matrices to data frames
  data <- as.data.frame(data)
  fit <- as.data.frame(fit)

  # Melt using .melt_df
  data_df <- .melt_df(data, years = years, ages = ages)
  fit_df <- .melt_df(fit, years = years, ages = ages)

  # Determine plotting limits
  max_X <- max(as.numeric(ages))
  max_Y <- max(c(as.numeric(data_df$v1), as.numeric(fit_df$v1)))

  # Setup plotting layout: one row per year, or adjust
  n_years <- length(years)
  
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  # Loop over each year to create “facets”
  for (y in 1:length(years)) {
    # Subset data for the current year
    data_sub <- data_df[data_df$year == years[y], ]
    fit_sub <- fit_df[fit_df$year == years[y], ]

    # Plot setup
    plot(as.numeric(data_sub$age), as.numeric(data_sub$v1),
      type = "p",
      ylim = c(0, max_Y),
      xlim = c(min(ages), max(ages)),
      xlab = "Age",
      ylab = "",
      pch = 16,
      col = "black"
    )
    text(x = max(ages) - 2, y = max_Y * 0.95, labels = paste(years[y], "\nn =", sample_size[y]))
    # Add fit line
    lines(as.numeric(fit_sub$age), as.numeric(fit_sub$v1), col = "blue", lwd = 2)
    # title
    if (y == 2) {
      if (!is.null(title)) {
        mtext(title, side = 3, line = 1.5, cex = 1.25, font = 2)
      }
    }
    grid()
  }
}
