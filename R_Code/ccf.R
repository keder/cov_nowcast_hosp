# Load Packages -----------------------------------------------------------
pacman::p_load(plotly, tidyverse, gridExtra, RColorBrewer)

joined_df <- read.csv("data/joined_df.csv", check.names = FALSE)

# Helper: tidy CCF output -------------------------------------------------
tidy_ccf <- function(ccf_obj, label) {
  tibble(
    acf = as.vector(ccf_obj$acf),
    lag = as.vector(ccf_obj$lag),
    type = label
  )
}

# Function ---------------------------------------------------------------
get_ccf <- function(df, state_nm) {
  message("Processing: ", state_nm)

  filtered_df <- df %>%
    filter(state == state_nm)

  # --- NVSS ---
  nvss <- filtered_df %>%
    select(state, date, percent_pos, count_pos,
           total_admissions_all_covid_confirmed) %>%
    filter(date >= "2022-10-01") %>%
    mutate(
      percent_pos = replace_na(percent_pos, 0),
      count_pos = replace_na(count_pos, 0),
      total_admissions_all_covid_confirmed =
        replace_na(total_admissions_all_covid_confirmed, 0)
    )

  # --- ED ---
  ed <- filtered_df %>%
    select(state, date, percent_visits_covid,
           total_admissions_all_covid_confirmed) %>%
    filter(date >= "2022-10-01") %>%
    mutate(
      percent_visits_covid = replace_na(percent_visits_covid, 0),
      total_admissions_all_covid_confirmed =
        replace_na(total_admissions_all_covid_confirmed, 0)
    )

  # --- NWSS ---
  nwss <- filtered_df %>%
    select(
      state, date,
      `avg_detect_prop_weighted_post grit removal`,
      `avg_detect_prop_weighted_raw wastewater`,
      `avg_detect_prop_weighted_primary sludge`,
      `avg_detect_prop_unweighted_raw wastewater`,
      total_admissions_all_covid_confirmed
    ) %>%
    filter(date >= "2022-10-01") %>%
    drop_na()

  # Save intermediate datasets (optional but useful)
  write_csv(nvss, paste0("results/debug_nvss_", state_nm, ".csv"))
  write_csv(ed, paste0("results/debug_ed_", state_nm, ".csv"))
  write_csv(nwss, paste0("results/debug_nwss_", state_nm, ".csv"))

  # --- Run CCFs ---
  ccf_list <- list(
    "% Positive" = ccf(nvss$percent_pos,
                       nvss$total_admissions_all_covid_confirmed,
                       lag.max = 5, plot = FALSE),

    "Count Positive" = ccf(nvss$count_pos,
                           nvss$total_admissions_all_covid_confirmed,
                           lag.max = 5, na.action = na.pass, plot = FALSE),

    "ED Visit" = ccf(ed$percent_visits_covid,
                     ed$total_admissions_all_covid_confirmed,
                     lag.max = 5, na.action = na.pass, plot = FALSE),

    "Weighted WW Post-Grit Removal % Detect" =
      ccf(nwss$`avg_detect_prop_weighted_post grit removal`,
          nwss$total_admissions_all_covid_confirmed,
          lag.max = 5, na.action = na.pass, plot = FALSE),

    "Weighted WW Raw % Detect" =
      ccf(nwss$`avg_detect_prop_unweighted_raw wastewater`,
          nwss$total_admissions_all_covid_confirmed,
          lag.max = 5, na.action = na.pass, plot = FALSE),

    "Weighted WW Primary Sludge % Detect" =
      ccf(nwss$`avg_detect_prop_weighted_primary sludge`,
          nwss$total_admissions_all_covid_confirmed,
          lag.max = 5, na.action = na.pass, plot = FALSE)
  )

  # --- Combine results cleanly ---
  result <- purrr::imap_dfr(ccf_list, tidy_ccf) %>%
    mutate(
      acf = replace_na(acf, 0),
      abs_acf = abs(acf)
    ) %>%
    filter(lag <= 0) %>%
    arrange(desc(abs_acf)) %>%
    slice_head(n = 1, by = type) %>%
    mutate(
      state = state_nm,
      lag = if_else(acf == 0, 0, lag)
    )

  # Save per-state CCF summary
  write_csv(result, paste0("results/debug_ccf_", state_nm, ".csv"))

  return(result)
}

# Iterate through states --------------------------------------------------
state_ccf <- purrr::map_dfr(states$state, function(st) {
  if (st == "USA") return(NULL)
  get_ccf(joined_df, st)
})

# Save final results
write_csv(state_ccf, "results/ccf.csv")


# Plot Heatmap ------------------------------------------------------------
ccf_df <- state_ccf %>%
  select(state, type, acf)

write_csv(ccf_df, "results/ccf_heatmap_input.csv")

p <- ggplot(ccf_df, aes(x = state, y = type, fill = acf)) +
  geom_tile() +
  geom_text(aes(label = if_else(acf == 0, "-", as.character(round(acf, 2)))),
            size = 2.7) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 18)) +
  scale_fill_distiller("Correlation", palette = "RdBu", limits = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0.45)
  )

ggsave("results/CCF_Heatmap.pdf",
       plot = p, width = 16, height = 3,
       units = "in", bg = "white", dpi = "retina")