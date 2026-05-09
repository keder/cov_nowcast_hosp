# Load Packages -----------------------------------------------------------
pacman::p_load(
  aTSA, stats, forecast, prophet, urca, tidyverse, xtable, xgboost,
  uroot
)

## Faster than auto.arima way to get optimal model, first CSS then ML
extract_loglik <- function(p, d, q, P, D, Q) {
  # print(paste("Order is", p, ",",d, ",", q, " | Season is", P, ",",D, ",", Q))
  arima(ts_dat,
    order = c(p, d, q), method = "CSS",
    seasonal = c(P, D, Q)
  )$loglik
}
extract_aic <- function(p, d, q, P, D, Q) {
  print(paste("Order is", p, ",", d, ",", q, " | Season is", P, ",", D, ",", Q))
  arima(ts_dat,
    order = c(p, d, q), method = "ML", transform.pars = FALSE,
    seasonal = c(P, D, Q)
  )$aic
}

add_perturbation <- function(x) {
  len <- length(x)
  scale_factor <- 1 / 1000
  if (sum(x, na.rm = T) == 0) {
    rand_vector <- runif(len)
    return((log(rand_vector) + 1) * scale_factor)
  } else {
    return(x)
  }
}

# Fit Arima for each state/multivariate model
fit_arima_regressors <- function(cut_date, data, state, freq, xreg = NULL) {
  train_data <- data %>% filter(date < cut_date)
  ts_dat <- ts(train_data$hosp_transformed, frequency = freq)

  training_xreg <- NULL
  test_xreg <- NULL
  if (!is.null(xreg)) {
    training_xreg <- xreg %>%
      filter(date < cut_date) %>%
      select(!date) %>%
      as.matrix()

    test_xreg <- xreg %>%
      filter(date >= cut_date) %>%
      select(!date) %>%
      as.matrix()
  }

  tryCatch(
    {
      if (is.null(training_xreg)) {
        model <- auto.arima(ts_dat,
          max.order = 30,
          # Max values come from arima fits of TS alone
          max.p = 7, max.q = 8,
          max.P = 3, max.Q = 3,
          max.d = 2, max.D = 2,
          stepwise = FALSE, approximation = FALSE,
          trace = FALSE, ic = "aic", seasonal.test = "ch"
        )
      } else {
        model <- auto.arima(ts_dat,
          max.order = 30,
          # Max values come from arima fits of TS alone
          max.p = 7, max.q = 8,
          max.P = 3, max.Q = 3,
          max.d = 2, max.D = 2,
          stepwise = FALSE, approximation = FALSE,
          trace = FALSE, ic = "aic", seasonal.test = "ch",
          xreg = training_xreg
        )
      }

      train_data$fitted <- model$fitted
      if (is.null(training_xreg)) {
        n <- n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
        forecasted <- as_tibble(forecast(model, h = n)) %>%
          bind_cols(data %>% filter(date >= cut_date)) %>%
          select(date, hosp, everything())
      } else {
        n <- n_distinct(data %>% filter(date >= cut_date) %>% pull(date))
        forecasted <- as_tibble(forecast(model, h = n, xreg = test_xreg)) %>%
          bind_cols(data %>% filter(date >= cut_date)) %>%
          select(date, hosp, everything())
      }

      train_data <- train_data %>%
        bind_rows(forecasted) %>%
        mutate(
          state = state, cut_date = cut_date,
          freq = freq,
          fitted = exp(fitted) - 1,
          `Point Forecast` = exp(`Point Forecast`) - 1
        ) %>%
        select(state, cut_date, date, hosp, everything())
      return(train_data)
    },
    error = function(err) {
      print(paste0("ARIMA for state ", state, " on date ", cut_date, " with xreg=", paste(names(xreg), collapse = ", "), ", failed with error: ", err$message))
      return(tibble())
    }
  )
}


# Set Prediction Horizons -------------------------------------------------
joined_df <- read.csv("data/joined_df.csv", check.names = FALSE)
full_data <- joined_df
full_data$date <- as.Date(full_data$date)
max_date <- max(full_data$date)
prediction_horizons <- c(
  max_date %m-% months(3), max_date %m-% months(6),
  max_date %m-% months(9), max_date %m-% months(12)
) %>%
  as.character()

# ARIMAX ------------------------------------------------------------------
freq <- 26
df_list <- tibble()

# For each state
for (state_nm in full_data$state) {
  print(state_nm)
  data <- full_data %>%
    filter(
      state == state_nm,
      date <= "2024-05-01"
    ) %>%
    janitor::clean_names() %>%
    dplyr::select(date,
      hosp = total_admissions_all_covid_confirmed, count_pos,
      percent_pos, percent_visits_covid,
      weighted_raw_ww = avg_detect_prop_weighted_raw_wastewater,
      weighted_postgrit_ww = avg_detect_prop_weighted_post_grit_removal,
      weighted_sludge_ww = avg_detect_prop_weighted_primary_sludge
    ) %>%
    # use box-cox transform on hospital data to avoid negatives in predictions
    # we assume box-cox order 0 (so log() transform, but add 1 to avoid fitting
    # log(0) in our model)
    mutate(
      count_pos = replace_na(count_pos, 0),
      count_transformed = log(count_pos + 1),
      percent_pos = replace_na(percent_pos, 0),
      percent_transformed = log(percent_pos + 1),
      percent_visits_covid = replace_na(percent_visits_covid, 0),
      ed_transformed = log(percent_visits_covid + 1),
      weighted_raw_ww = replace_na(weighted_raw_ww, 0),
      raw_pct_transformed = log(weighted_raw_ww + 1),
      weighted_postgrit_ww = replace_na(weighted_postgrit_ww, 0),
      postgrit_pct_transformed = log(weighted_postgrit_ww + 1),
      weighted_sludge_ww = replace_na(weighted_sludge_ww, 0),
      sludge_pct_transformed = log(weighted_sludge_ww + 1),
      hosp_transformed = log(hosp + 1)
    ) %>%
    filter(date >= "2022-10-01") %>%
    drop_na() %>%
    mutate(across(
      sludge_pct_transformed,
      ~ add_perturbation(.x)
    ))
  for (predict_cut in prediction_horizons) {
    # Run models for each prediction horizon and set of covariates
    print(paste0(state_nm, predict_cut))
    no_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = NULL
    ) %>%
      mutate(model_type = "No exogenous regressors")
    ed_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = data %>%
        select(date, ed_transformed)
    ) %>%
      mutate(model_type = "% ED Visits")
    ed_count_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = data %>%
        select(
          date, ed_transformed,
          count_transformed,
          percent_transformed
        )
    ) %>%
      mutate(model_type = "% ED Visits and Test Positivity")
    ww_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = data %>%
        select(date, raw_pct_transformed)
    ) %>%
      mutate(model_type = "Wastewater % Detection (Raw)")
    ww_all_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = data %>%
        select(
          date, raw_pct_transformed,
          postgrit_pct_transformed,
          sludge_pct_transformed
        )
    ) %>%
      mutate(model_type = "Wastewater % Detection (Raw, Post-Grit, Sludge)")
    all_reg <- fit_arima_regressors(predict_cut, data, state_nm, freq,
      xreg = data %>%
        select(
          date, ed_transformed,
          count_transformed,
          percent_transformed,
          raw_pct_transformed,
          postgrit_pct_transformed,
          sludge_pct_transformed
        )
    ) %>%
      mutate(model_type = "All Regressors")
    df_list <- bind_rows(no_reg, ed_reg) %>%
      bind_rows(ed_count_reg) %>%
      bind_rows(ww_reg) %>%
      bind_rows(ww_all_reg) %>%
      bind_rows(all_reg)
  }
}

# Write out results
dir.create("results/Multivariate", recursive = TRUE)
write_rds(df_list, paste0("results/Multivariate/arimaX_", freq, "final.rds"))

# Plot
for (predict_cut in prediction_horizons) {
  for (type in unique(df_list$model_type)) {
    p <- df_list %>%
      select(
        state, cut_date, date, hosp, hosp_transformed,
        fitted, `Point Forecast`, model_type
      ) %>%
      arrange(state) %>%
      rename(`Fitted Period` = fitted, `Forecast Period` = `Point Forecast`) %>%
      #      filter(is.na(`Forecast Period`) | `Forecast Period` <= max(hosp) * 100 ) %>%
      pivot_longer(`Fitted Period`:`Forecast Period`) %>%
      # mutate(value = exp(value) - 1,
      #        hosp = exp(hosp) - 1) %>%
      drop_na() %>%
      filter(
        cut_date == predict_cut,
        date >= "2022-07-01",
        model_type == type
      ) %>%
      ggplot(aes(x = date)) +
      geom_line(aes(y = hosp, color = "Reported \nHospitalizations")) +
      geom_point(aes(y = value, color = name), size = 0.5) +
      theme_minimal() +
      scale_color_manual(values = c("#1B9E77", "#7570B3", "#000000")) +
      labs(y = "Hospitalizations", x = "Date", color = "") +
      facet_wrap(state ~ ., ncol = 4, scales = "free_y")
    ggsave(
      paste0(
        "results/Multivariate/Multivariate Plots/ArimaX seasonal results for ", gsub("%", "PCT", type, fixed = T),
        " freqency time series at ", predict_cut, " horizon", ".pdf"
      ),
      plot = p, width = 16, height = 10, units = "in", bg = "white",
      dpi = "retina"
    )
  }
}
