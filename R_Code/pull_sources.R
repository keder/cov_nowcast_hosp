# Load Packages -----------------------------------------------------------
pacman::p_load(tidyverse, httr, lubridate, MMWRweek)

# Constants ---------------------------------------------------------------
START_DATE <- as.Date("2022-10-01")
END_DATE <- as.Date("2024-05-01")

# Helpers -----------------------------------------------------------------
pull_api <- function(url, limit = 5000000) {
  limit <- as.character(as.integer(limit))
  read_csv(paste0(url, "?$limit=", limit), show_col_types = FALSE)
}

add_mmwr <- function(df, date_col = "date") {
  df %>%
    mutate(
      date = as.Date(.data[[date_col]]),
      mmwr = MMWRweek::MMWRweek(date),
      mmwr_date = sprintf("%d-%02d", mmwr$MMWRyear, mmwr$MMWRweek)
    ) %>%
    select(-mmwr)
}

add_hhs_region <- function(df) {
  df %>%
    mutate(
      level = case_when(
        grepl("CT|ME|MA|NH|RI|VT", state) ~ "Region 1",
        grepl("NJ|NY|PR|VI", state) ~ "Region 2",
        grepl("DE|DC|MD|PA|VA|WV", state) ~ "Region 3",
        grepl("AL|FL|GA|KY|MS|NC|SC|TN", state) ~ "Region 4",
        grepl("IL|IN|MI|MN|OH|WI", state) ~ "Region 5",
        grepl("AR|LA|NM|OK|TX", state) ~ "Region 6",
        grepl("IA|KS|MO|NE", state) ~ "Region 7",
        grepl("CO|MT|ND|SD|UT|WY", state) ~ "Region 8",
        grepl("AZ|CA|HI|NV|AS|GU|MP", state) ~ "Region 9",
        grepl("AK|ID|OR|WA", state) ~ "Region 10",
        state == "USA" ~ "National",
        TRUE ~ NA_character_
      )
    )
}

interpolate_or_zero <- function(x, min_points = 2) {
  non_na <- !is.na(x)

  if (sum(non_na) >= min_points) {
    return(
      approx(
        x = which(non_na),
        y = x[non_na],
        xout = seq_along(x),
        rule = 2
      )$y
    )
  }

  replace(x, non_na, 0)
}

# State crosswalk ---------------------------------------------------------
states <- tibble(
  state = c(state.abb, "DC", "USA"),
  geography = c(state.name, "District of Columbia", "United States")
)

state_to_abbr <- setNames(states$state, states$geography)

# Hospitalizations --------------------------------------------------------
hospitalizations <- pull_api("https://data.cdc.gov/resource/aemt-mg7g.csv") %>%
  transmute(
    date = as.Date(week_end_date),
    state = jurisdiction,
    total_admissions_all_covid_confirmed = total_admissions_all_covid_confirmed
  ) %>%
  add_mmwr("date") %>%
  add_hhs_region() %>%
  filter(date <= END_DATE)

# NVSS --------------------------------------------------------------------
nvss <- pull_api("https://data.cdc.gov/resource/gvsb-yw6g.csv") %>%
  mutate(date = as.Date(mmwrweek_end)) %>%
  add_mmwr("date") %>%
  transmute(
    date,
    level,
    mmwr_date,
    percent_pos,
    count_pos = percent_pos * number_tested,
    posted
  ) %>%
  group_by(date, level) %>%
  filter(posted == max(posted), date <= END_DATE) %>%
  ungroup() %>%
  distinct()

# NSSP --------------------------------------------------------------------
nssp <- pull_api("https://data.cdc.gov/resource/rdmq-nq56.csv") %>%
  filter(county == "All", week_end <= END_DATE) %>%
  mutate(date = as.Date(week_end)) %>%
  add_mmwr("date") %>%
  transmute(
    state = state_to_abbr[geography],
    date,
    mmwr_date,
    percent_visits_covid = replace_na(percent_visits_covid, 0)
  )

# ----------------------------- NWSS METRIC -----------------------------
nwss_metric <- pull_api("https://data.cdc.gov/resource/2ew6-ywp6.csv") %>%
  select(
    wwtp_jurisdiction:key_plot_id,
    population_served,
    date = date_end,
    detect_prop_15d
  ) %>%
  distinct() %>%
  mutate(
    date = as.Date(date),
    mmwr = MMWRweek::MMWRweek(date),
    mmwr_date = sprintf("%d-%02d", mmwr$MMWRyear, mmwr$MMWRweek),
    state = state.abb[match(reporting_jurisdiction, state.name)]
  ) %>%
  filter(!is.nan(detect_prop_15d)) %>%
  group_by(state, key_plot_id, sample_location, mmwr_date) %>%
  summarise(
    population_served = mean(population_served, na.rm = TRUE),
    detect_prop_15d = mean(na.omit(detect_prop_15d), na.rm = TRUE),
    .groups = "drop"
  ) %>% # Average data across weekdays for each location
  group_by(state, sample_location, key_plot_id) %>%
  arrange(mmwr_date, .by_group = TRUE) %>% # Sort by mmwr_date
  mutate(
    detect_prop_15d = interpolate_or_zero(detect_prop_15d)
  )

# ----------------------------- NWSS FILTER -----------------------------
nwss <- nwss_metric %>%
  filter(sample_location == "Treatment plant") %>% # Leave only last point of water treatment
  mutate(
    location = gsub("90_|89_", "", str_split_i(key_plot_id, "Treatment plant_", 2))
  ) %>% # Extract only type of sampling (sample location as in where in the facility), e.g. raw wastewater, post grit removal, primary sludge
  ungroup() %>%
  transmute(
    state,
    key_plot_id,
    mmwr_date,
    location,
    population_served,
    detect_prop_15d
  )

# ----------------------------- NWSS GRID COMPONENTS --------------------
nwss_dates_state <- nwss %>% distinct(state)
nwss_dates_location <- nwss %>% distinct(location)
nwss_dates_time <- nwss %>% distinct(mmwr_date)

# ----------------------------- NWSS AGGREGATION ------------------------
nwss_metric_wide <- nwss %>%
  filter(!is.nan(detect_prop_15d), !is.na(population_served)) %>%
  group_by(state, mmwr_date, location) %>%
  summarise(
    total_pop = sum(population_served, na.rm = TRUE),
    avg_detect_prop_weighted = weighted.mean(detect_prop_15d, population_served, na.rm = TRUE), # Percentage detected weighted by population
    avg_detect_prop_unweighted = mean(detect_prop_15d, na.rm = TRUE), # Percentage detected unweighted average
    .groups = "drop"
  ) %>%
  select(-total_pop)

# ----------------------------- NWSS WIDE -------------------------------
nwss_wide <- nwss_metric_wide %>%
  full_join(
    tidyr::crossing(
      state = nwss_dates_state$state,
      location = nwss_dates_location$location,
      mmwr_date = nwss_dates_time$mmwr_date
    ), # Project data on fully covering grid
    by = c("state", "location", "mmwr_date")
  ) %>%
  group_by(state, location) %>%
  arrange(mmwr_date, .by_group = TRUE) %>%
  mutate(
    across(
      c(
        avg_detect_prop_weighted,
        avg_detect_prop_unweighted
      ),
      interpolate_or_zero
    )
  ) %>% # Fill data missing when we projected existing data on all combinations of "state", "location", "mmwr_date"
  ungroup() %>%
  filter(!is.na(state)) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0) / 100)) %>%
  pivot_wider(
    names_from = location,
    values_from = c(avg_detect_prop_weighted, avg_detect_prop_unweighted)
  )

# FINAL JOIN --------------------------------------------------------------
joined_df <- hospitalizations %>%
  left_join(nvss, by = c("level", "date", "mmwr_date")) %>%
  left_join(nssp, by = c("state", "date", "mmwr_date")) %>%
  left_join(nwss_wide, by = c("state", "mmwr_date")) %>%
  filter(
    date <= as.Date(END_DATE),
    date >= as.Date(START_DATE)
  ) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

write_csv(joined_df, "data/joined_df.csv", na = "NA")
