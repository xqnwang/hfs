#--------------------------------------------------------------------
# Combine results of different methods
#--------------------------------------------------------------------
combine_table <- function(data_label, methods, measure, scenario = NULL, horizons) {
  for (method_label in methods) {
    data_method <- paste0(data_label, "_", method_label)
    for (h in horizons) {
      assign(
        paste0(data_method, "_", measure, "_h", h),
        readRDS(file = paste0(
          "data_new/", data_method, "_reconsf",
          ifelse(is.null(scenario), "", paste0("_", scenario)), "_", measure, "_", h, ".rds"
        ))
      )
    }
    example <- get(paste0(data_method, "_", measure, "_h", horizons[1]))
    levels <- colnames(example)[-1]
    method <- sub("_", "-", example$Method)
    out <- lapply(levels,
      function(level) {
        mget(paste0(data_method, "_", measure, "_h", horizons), inherits = TRUE) |>
        sapply(function(lentry) as.numeric(lentry[[level]]))
      }
    )
    out <- data.frame(method, do.call(cbind, out))
    colnames(out) <- c("Method", c(ifelse(horizons == 1, paste0("h=", 1), paste0("1--", horizons))) |> rep(length(levels)))
    assign(paste0(data_method, "_", measure), out)
  }

  out_com <- mget(paste0(data_label, "_", methods, "_", measure), inherits = TRUE) %>% do.call(rbind, .)
  rownames(out_com) <- NULL
  out_com <- out_com[!duplicated(out_com), ]
  if (grepl("tourism", data_label) | grepl("labour", data_label)) {
    candidates <- c("OLS", "WLSs", "WLSv", "MinTs")
  } else if (data_label == "simulation") {
    candidates <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  }
  target <- c(
    "Base", "BU",
    sapply(
      candidates,
      function(len) c(len, paste0(len, "-", c("subset", "intuitive", "lasso")))
    ) |> as.vector(),
    "EMinT", "Elasso"
  )
  out_com <- out_com[match(target, out_com$Method), ]
  colnames_out <- colnames(out_com)
  table_out <- sapply(1:NROW(out_com), function(i) {
    r1 <- out_com[1, -1] |> as.numeric()
    rt <- out_com[i, -1] |> as.numeric()
    if (i == 1) {
      return(r1)
    } else {
      return(100 * (rt - r1) / r1)
    }
  }) |> t()
  table_out <- data.frame(out_com[, 1], table_out)
  colnames(table_out) <- colnames_out
  return(list(table_out = table_out, levels = levels, candidates = candidates))
}

#--------------------------------------------------------------------
# Output table latex
#--------------------------------------------------------------------
latex_table <- function(out_all) {
  out <- out_all$table_out
  levels <- sub("_", " x ", out_all$levels)
  horizons <- (ncol(out) - 1) / length(levels)
  header <- c("", rep(as.character(horizons), length(levels)))
  names(header) <- c("", levels)
  candidates <- out_all$candidates

  # Blue entries identify the best performing approaches
  # Bold entries identify methods that perform better than the corresponding benchmark method
  out[, -1] <- lapply(out[, -1], function(x) {
    x <- round(x, 1)
    origin_x <- x
    min_x <- min(origin_x)
    comp_x <- c(
      origin_x[1:2],
      rep(origin_x[3 + 4 * (0:(length(candidates) - 1))], each = 4),
      rep(origin_x[length(origin_x) - 1], each = 2)
    )
    out_x <- ifelse(x == min(x),
      cell_spec(format(x, nsmall = 1), bold = TRUE, color = "blue"),
      ifelse(x < comp_x, cell_spec(format(x, nsmall = 1), bold = TRUE), format(x, nsmall = 1))
    )
    out_x <- sub("-", "--", out_x)
    out_x
  })

  out |>
    kable(
      format = "latex",
      booktabs = TRUE,
      digits = 1,
      align = c("l", rep("r", ncol(out) - 1)),
      escape = FALSE,
      linesep = ""
    ) |>
    row_spec(2 + 4 * (0:length(candidates)), hline_after = TRUE) |>
    kable_paper(full_width = FALSE) |>
    kable_styling(
      latex_options = c("repeat_header", "scale_down"),
      font_size = 10
    ) |>
    row_spec(
      which(grepl("-", out$Method) | grepl("Elasso", out$Method)),
      background = "#e6e3e3"
    ) |>
    add_header_above(header, align = "c") |>
    # footnote(
    #   general_title = "",
    #   general = "Note: The Base row shows the average RMSE of the base forecasts. Entries below this row indicate the percentage decrease (negative) or increase (positive) in the average RMSE of the reconciled forecasts compared to the base forecasts. The lowest values in each column are highlighted in blue. Proposed methods are shaded in gray, and those outperforming the benchmark method are marked in bold.",
    #   footnote_as_chunk = TRUE,
    #   threeparttable = TRUE,
    #   fixed_small_size = FALSE
    # ) |>
    print()
}

latex_table_h <- function(out_all) {
  out <- out_all$table_out
  colnames(out) <- sub("_", " x ", colnames(out))
  levels <- sub("_", " x ", out_all$levels)
  horizons <- (ncol(out) - 1) / length(levels)
  header <- c("", rep(as.character(horizons), length(levels)))
  names(header) <- c("", levels)
  candidates <- out_all$candidates
  
  # Blue entries identify the best performing approaches
  # Bold entries identify methods that perform better than the corresponding benchmark method
  out[, -1] <- lapply(out[, -1], function(x) {
    x <- round(x, 1)
    origin_x <- x
    min_x <- min(origin_x)
    comp_x <- c(
      origin_x[1:2],
      rep(origin_x[3 + 4 * (0:(length(candidates) - 1))], each = 4),
      rep(origin_x[length(origin_x) - 1], each = 2)
    )
    out_x <- ifelse(x == min(x),
                    cell_spec(format(x, nsmall = 1), bold = TRUE, color = "blue"),
                    ifelse(x < comp_x, cell_spec(format(x, nsmall = 1), bold = TRUE), format(x, nsmall = 1))
    )
    out_x <- sub("-", "--", out_x)
    out_x
  })
  
  out |>
    kable(
      format = "latex",
      booktabs = TRUE,
      digits = 1,
      align = c("l", rep("r", ncol(out) - 1)),
      escape = FALSE,
      linesep = ""
    ) |>
    row_spec(2 + 4 * (0:length(candidates)), hline_after = TRUE) |>
    # kable_paper(full_width = FALSE) |>
    kable_styling(latex_options = c("repeat_header"), font_size = 10) |>
    # kable_styling(
    #   latex_options = c("repeat_header", "scale_down"),
    #   font_size = 8
    # ) |>
    row_spec(
      which(grepl("-", out$Method) | grepl("Elasso", out$Method)),
      background = "#e6e3e3"
    ) |>
    column_spec(1:ncol(out), width = paste0(5.5/ncol(out), "in")) %>%
    print()
}

#--------------------------------------------------------------------
# Extract element from forecast results
#--------------------------------------------------------------------
extract_element <- function(data, index, method, element) {
  out <- data[[index]][[method]][[element]]
  if (length(out) == 1) {
    if (is.na(out)) {
      out <- NULL
    }
  } else {
    if (is.vector(out)) {
      out <- cbind(matrix(out, nrow = 1), Index = index)
    } else {
      out <- cbind(out, Index = index)
    }
  }
  out
}

#--------------------------------------------------------------------
# Calculate error
#--------------------------------------------------------------------
calc_error <- function(fc, test, h) {
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  err <- cbind(err, Index = subset(test, select = Index))
  err_h <- err |>
    as_tibble() |>
    group_by(Index) |>
    mutate(Horizon = row_number()) |>
    filter(Horizon <= h) |>
    ungroup() |>
    select(!c("Index", "Horizon"))
  return(err_h)
}

#--------------------------------------------------------------------
# Calculate RMSE
#--------------------------------------------------------------------
calc_rmse <- function(fc, test, h) {
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  err <- cbind(err, Index = subset(test, select = Index))
  rmse <- err |>
    as_tibble() |>
    group_by(Index) |>
    mutate(Horizon = row_number()) |>
    filter(Horizon <= h) |>
    summarise_at(1:NCOL(err), function(x) sqrt(mean(x^2))) |>
    ungroup() |>
    select(!c("Index", "Horizon")) |>
    summarise_all(mean)
  return(rmse)
}

calc_rmse_series <- function(fc, test, h){
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  err <- cbind(err, Index = subset(test, select = Index))
  rmse <- err |> 
    as_tibble() |> 
    group_by(Index) |> 
    mutate(Horizon = row_number()) |> 
    filter(Horizon <= h) |>
    summarise_at(1:NCOL(err), function(x) sqrt(mean(x^2))) |>
    ungroup() |>
    select(!c("Index", "Horizon")) |>
    rowwise() |>
    mutate(Average = mean(c_across(1:(NCOL(err)-1)))) |>
    pull(Average)
  return(rmse)
}

#--------------------------------------------------------------------
# Calculate MASE
#--------------------------------------------------------------------
calc_mase <- function(fc, train, test, freq, h) {
  x <- subset(train, select = -Index)
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  scaling <- apply(x, 2, function(s) mean(abs(diff(as.vector(s), freq))))
  q <- cbind(t(t(err) / scaling), Index = subset(test, select = Index))

  mase <- q |>
    as_tibble() |>
    group_by(Index) |>
    mutate(Horizon = row_number()) |>
    filter(Horizon <= h) |>
    summarise_at(1:NCOL(q), function(x) mean(abs(x))) |>
    ungroup() |>
    select(!c("Index", "Horizon")) |>
    summarise_all(mean)

  return(mase)
}

#--------------------------------------------------------------------
# Calculate RMSSE
#--------------------------------------------------------------------
calc_rmsse <- function(fc, train, test, freq, h) {
  x <- subset(train, select = -Index)
  err <- subset(test, select = -Index) - subset(fc, select = -Index)
  scaling <- apply(x, 2, function(s) mean(abs(diff(as.vector(s), freq))))
  q <- cbind(t(t(err) / scaling), Index = subset(test, select = Index))
  
  rmsse <- q |>
    as_tibble() |>
    group_by(Index) |>
    mutate(Horizon = row_number()) |>
    filter(Horizon <= h) |>
    summarise_at(1:NCOL(q), function(x) sqrt(mean(x^2))) |>
    ungroup() |>
    select(!c("Index", "Horizon")) |>
    summarise_all(mean)
  
  return(rmsse)
}

#--------------------------------------------------------------------
# Combine z outputs from results of different methods
#--------------------------------------------------------------------
combine_z <- function(data_label, methods, scenarios, series_name) {
  for (scenario in scenarios) {
    for (method_label in methods) {
      target <- paste0(
        data_label, "_", method_label, "_reconsf",
        ifelse(scenario == "s0", "", paste0("_", scenario))
      )

      assign(
        target,
        readRDS(file = paste0("data_new/", target, ".rds"))
      )
      z_method <- lapply(get(target), function(lentry) {
        select_methods <- names(lentry)[grepl("_", names(lentry)) | grepl("Elasso", names(lentry))]
        z_index <- sapply(lentry[select_methods], function(len) c(len$z, sum(len$z))) |>
          t() |>
          round(0)
        colnames(z_index) <- c(series_name, "Total")
        return(z_index)
      }) %>% do.call(rbind, .)
      assign(paste0("z_", method_label), rowsum(subset(z_method, select = -Total), row.names(z_method)))
      assign(paste0("n_", method_label), data.frame(Method = row.names(z_method), Total = subset(z_method, select = Total)))
    }
    z_scenario <- mget(paste0("z_", methods)) %>% do.call(rbind, .)
    rownames(z_scenario) <- sub("_", "-", rownames(z_scenario))
    n_scenario <- mget(paste0("n_", methods)) %>% do.call(rbind, .)
    n_scenario$Method <- sub("_", "-", n_scenario$Method)
    rownames(n_scenario) <- NULL
    assign(
      paste0("out_", scenario),
      list(
        z = assign(paste0("z_", scenario), z_scenario),
        n = assign(paste0("n_", scenario), n_scenario)
      )
    )
  }
  mget(paste0("out_", scenarios))
}

#--------------------------------------------------------------------
# Output table latex for number of time series retained after subset selection for the simulation data
#--------------------------------------------------------------------
latex_sim_nos_table <- function(z_out, n_out, rmsse_out, label_out) {
  candidates <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  target <- sapply(candidates, function(len) {
    c(paste0(len, "-", c("subset", "parsim", "lasso")))
  }) |> as.vector()
  target <- c(target, "Elasso")
  z_out <- z_out[match(target, row.names(z_out)), ] / 500
  z_out <- data.frame(row.names(z_out), z_out, Summary = "")
  if ("AA" %in% colnames(z_out)) colnames(z_out) <- sub("Top", "Total", colnames(z_out))
  header <- c("", colnames(z_out)[-1])
  colnames(z_out) <- c("(RMSSE)", paste0("(", format(round(as.numeric(rmsse_out), digits = 2), 2), ")"), "")
  row.names(z_out) <- NULL
  n_img <- split(n_out$Total, n_out$Method)
  n_img <- n_img[match(target, names(n_img))]
  n_img <- lapply(n_img, function(N) data.frame(table(N), Method = "method"))

  colors_used <- hcl.colors(10, "Purples")[3:6]
  names(colors_used) <- 7:4

  inline_bars <-
    n_img %>%
    map(~ ggplot(.x, aes(x = Method, y = Freq, fill = N)) +
      geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
      scale_fill_manual(values = colors_used) +
      coord_cartesian(ylim = c(0, 500), clip = "off") +
      coord_flip() +
      theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_blank()
      ))

  map(1:length(n_img), function(i) {
    ggsave(
      filename = paste0(label_out, "_", names(n_img)[i], ".png"),
      path = "_figs/",
      plot = inline_bars[[i]], height = 1.5, width = 7, dpi = 300,
      create.dir = TRUE
    )
  })
  ls_inline_plots <- file.path(getwd(), paste0("_figs/", label_out, "_", names(n_img), ".png"))

  z_out |>
    kable(
      format = "latex",
      booktabs = TRUE,
      digits = 2,
      align = c("l", rep("r", NCOL(z_out)-1)),
      escape = FALSE,
      linesep = ""
    ) |>
    add_header_above(header, line = FALSE) |>
    kable_styling(font_size = 10) |>
    kable_paper(full_width = FALSE) |>
    row_spec(3 * (1:length(candidates)), hline_after = TRUE) |>
    # footnote(
    #   general_title = "",
    #   general = "Note: The last column displays a stacked barplot for each method, based on the total number of selected series from 500 simulation instances. A darker sub-bar indicates a higher count.",
    #   escape = TRUE,
    #   footnote_as_chunk = TRUE,
    #   threeparttable = TRUE,
    #   fixed_small_size = FALSE
    # ) |>
    column_spec(NCOL(z_out),
                     image = spec_image(ls_inline_plots, width = 140, height = 30)
    )
}

#--------------------------------------------------------------------
# Combine results of different methods and correlation coefficients for corr sumulation data
#--------------------------------------------------------------------
combine_corr_table <- function(data_label, methods, corr, index, measure) {
  for (method_label in methods) {
    for (i in index) {
      data_method <- paste0(data_label, "_", i, "_", method_label)
      assign(
        data_method,
        readRDS(file = paste0(
          "data_new/", data_method, "_reconsf",
          "_", measure, "_1.rds"
        ))
      )
    }
    example <- get(paste0(data_label, "_", index[1], "_", method_label))
    levels <- colnames(example)[-1]
    method <- sub("_", "-", example$Method)
    out <- lapply(levels, function(level) {
      mget(paste0(data_label, "_", index, "_", method_label), inherits = TRUE) |>
        sapply(function(lentry) as.numeric(lentry[[level]]))
    })
    out <- data.frame(method, do.call(cbind, out))
    colnames(out) <- c("Method", c(ifelse(corr[index] == -0.8, "$\\rho$=--0.8", sub("-", "--", corr[index]))) |> rep(length(levels)))
    assign(paste0(data_label, "_", method_label), out)
  }

  out_com <- mget(paste0(data_label, "_", methods), inherits = TRUE) %>% do.call(rbind, .)
  rownames(out_com) <- NULL
  out_com <- out_com[!duplicated(out_com), ]
  if (data_label %in% c("tourism", "labour")) {
    candidates <- c("OLS", "WLSs", "WLSv", "MinTs")
  } else if (data_label %in% c("simulation", "corr")) {
    candidates <- c("OLS", "WLSs", "WLSv", "MinT", "MinTs")
  }
  target <- c(
    "Base", "BU",
    sapply(
      candidates,
      function(len) c(len, paste0(len, "-", c("subset", "intuitive", "lasso")))
    ) |> as.vector(),
    "EMinT", "Elasso"
  )
  out_com <- out_com[match(target, out_com$Method), ]
  colnames_out <- colnames(out_com)
  table_out <- sapply(1:NROW(out_com), function(i) {
    r1 <- out_com[1, -1] |> as.numeric()
    rt <- out_com[i, -1] |> as.numeric()
    if (i == 1) {
      return(r1)
    } else {
      return(100 * (rt - r1) / r1)
    }
  }) |> t()
  table_out <- data.frame(out_com[, 1], table_out)
  colnames(table_out) <- colnames_out
  return(list(table_out = table_out, levels = levels, candidates = candidates))
}

#--------------------------------------------------------------------
# Output table latex for corr simulation data
#--------------------------------------------------------------------
latex_corr_table <- function(out_all) {
  out <- out_all$table_out
  correlations <- (ncol(out) - 1) / length(levels)
  levels <- sub("_", " x ", out_all$levels)
  header <- c("", rep(as.character((ncol(out) - 1) / length(levels)), length(levels)))
  names(header) <- c("", levels)
  candidates <- out_all$candidates

  # Blue entries identify the best performing approaches
  # Bold entries identify methods that perform better than the corresponding benchmark method
  out[, -1] <- lapply(out[, -1], function(x) {
    x <- round(x, 1)
    origin_x <- x
    min_x <- min(origin_x)
    comp_x <- c(
      origin_x[1:2],
      rep(origin_x[3 + 4 * (0:(length(candidates) - 1))], each = 4),
      rep(origin_x[length(origin_x) - 1], each = 2)
    )
    out_x <- ifelse(x == min(x), cell_spec(format(x, nsmall = 1), bold = TRUE, color = "blue"), ifelse(x < comp_x, cell_spec(format(x, nsmall = 1), bold = TRUE), format(x, nsmall = 1)))
    out_x <- sub("-", "--", out_x)
    out_x
  })

  out |>
    kable(
      format = "latex",
      booktabs = TRUE,
      digits = 1,
      align = c("l", rep("r", ncol(out) - 1)),
      escape = FALSE,
      linesep = ""
    ) |>
    row_spec(2 + 4 * (0:length(candidates)), hline_after = TRUE) |>
    kable_paper(full_width = FALSE) |>
    kable_styling(
      latex_options = c("repeat_header", "scale_down"),
      font_size = 10
    ) |>
    row_spec((grepl("-", out$Method) | grepl("Elasso", out$Method)) |> which(),
      background = "#e6e3e3"
    ) |>
    add_header_above(header, align = "c") |>
    # footnote(
    #   general_title = "",
    #   general = "Note: The Base row shows the average RMSE of the base forecasts. Entries below this row indicate the percentage decrease (negative) or increase (positive) in the average RMSE of the reconciled forecasts compared to the base forecasts. The lowest values in each column are highlighted in blue. Proposed methods are shaded in gray, and those outperforming the benchmark method are marked in bold.",
    #   footnote_as_chunk = TRUE,
    #   threeparttable = TRUE,
    #   fixed_small_size = FALSE
    # ) |>
    print()
}

#--------------------------------------------------------------------
# paste element-wise
#--------------------------------------------------------------------
paste_matrices <- function(mat1, mat2, sep = ",") {
  result <- matrix(paste(mat1, mat2, sep = sep), nrow = nrow(mat1), ncol = ncol(mat1))
  return(result)
}
