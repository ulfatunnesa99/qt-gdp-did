# =====================================================================
# Quantitative Tightening and U.S. GDP
# Event-study / interrupted time series on public FRED data.
# =====================================================================

## ---- 0. packages ---------------------------------------------------
pkgs <- c("dplyr", "lubridate", "lmtest", "sandwich", "ggplot2", "readr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create("results", showWarnings = FALSE)

## ---- 1. load data --------------------------------------------------
# Primary: a combined FRED CSV at data/fredgraph.csv with columns
#   observation_date, GDP_PC1, UNRATE, FEDFUNDS, WALCL, DGS10
# Fallback: pull the same series directly from FRED with quantmod.
load_data <- function(csv_path = "data/fredgraph.csv") {
  if (file.exists(csv_path)) {
    raw <- readr::read_csv(csv_path, show_col_types = FALSE)
    names(raw)[names(raw) == "observation_date"] <- "date"
    raw$date <- as.Date(raw$date)
    return(raw)
  }
  message("data/fredgraph.csv not found; pulling from FRED via quantmod ...")
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    install.packages("quantmod", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(quantmod))
  syms <- c("A191RL1Q225SBEA", "UNRATE", "FEDFUNDS", "WALCL", "DGS10")
  quantmod::getSymbols(syms, src = "FRED")
  g <- data.frame(date = zoo::index(A191RL1Q225SBEA),
                  GDP_PC1 = as.numeric(A191RL1Q225SBEA))
  ctrls <- list(UNRATE = UNRATE, FEDFUNDS = FEDFUNDS, WALCL = WALCL, DGS10 = DGS10)
  for (nm in names(ctrls)) {
    x <- ctrls[[nm]]
    df <- data.frame(date = zoo::index(x), val = as.numeric(x))
    df$q <- lubridate::floor_date(df$date, "quarter")
    agg <- aggregate(val ~ q, df, mean, na.rm = TRUE)
    names(agg) <- c("q", nm)
    g$q <- lubridate::floor_date(g$date, "quarter")
    g <- merge(g, agg, by = "q", all.x = TRUE)
    g$q <- NULL
  }
  g
}

## ---- 2. build quarterly panel -------------------------------------
build_panel <- function(raw) {
  raw %>%
    mutate(q = lubridate::floor_date(date, "quarter")) %>%
    group_by(q) %>%
    summarise(
      growth   = mean(GDP_PC1,  na.rm = TRUE),
      UNRATE   = mean(UNRATE,   na.rm = TRUE),
      FEDFUNDS = mean(FEDFUNDS, na.rm = TRUE),
      WALCL    = mean(WALCL,    na.rm = TRUE),
      DGS10    = mean(DGS10,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(date = q) %>%
    filter(!is.nan(growth)) %>%
    arrange(date) %>%
    mutate(
      qt1   = as.integer(date >= as.Date("2017-10-01") & date <= as.Date("2019-09-30")),
      qt2   = as.integer(date >= as.Date("2022-04-01")),
      covid = as.integer(date >= as.Date("2020-01-01") & date <= as.Date("2020-06-30")),
      t     = as.numeric(date - min(date)) / 365.25
    )
}

## ---- 3. interrupted time series -----------------------------------
# growth_t = b0 + b1*QT1_t + b2*QT2_t + trend + COVID + cyclical control
# b1, b2 are the average growth differentials during each QT window versus
# the non-QT baseline, with HC1 robust standard errors. This is a
# quasi-experimental event-study design, not a fully identified causal
# effect (see README).
main <- function() {
  panel <- build_panel(load_data())
  readr::write_csv(panel, "results/panel_quarterly.csv")
  cat("Quarters:", nrow(panel),
      " range:", format(min(panel$date)), "to", format(max(panel$date)), "\n")
  
  f_primary <- growth ~ qt1 + qt2 + t + covid + UNRATE
  f_full    <- growth ~ qt1 + qt2 + t + covid + UNRATE + FEDFUNDS + WALCL + DGS10
  
  m1 <- lm(f_primary, data = panel)
  r1 <- coeftest(m1, vcov = vcovHC(m1, type = "HC1"))
  cat("\n==== Primary ITS (trend + unemployment) ====\n"); print(r1)
  
  m2 <- lm(f_full, data = panel)
  r2 <- coeftest(m2, vcov = vcovHC(m2, type = "HC1"))
  cat("\n==== Robustness ITS (full controls; note: FEDFUNDS and WALCL",
      "proxy the policy itself and partly absorb the QT channel) ====\n")
  print(r2)
  
  b <- coef(m1)
  cat("\n--- Headline estimates (primary spec) ---\n")
  cat(sprintf("QT1 growth differential vs baseline: %+.2f pp\n", b[["qt1"]]))
  cat(sprintf("QT2 growth differential vs baseline: %+.2f pp\n", b[["qt2"]]))
  cat(sprintf("QT2 minus QT1 contrast:              %+.2f pp\n", b[["qt2"]] - b[["qt1"]]))
  cat("\nUse THESE printed numbers in the resume bullet. Nothing is hardcoded.\n")
  
  out <- data.frame(term = rownames(r1), estimate = r1[, 1], std_error = r1[, 2],
                    t = r1[, 3], p_value = r1[, 4], row.names = NULL)
  readr::write_csv(out, "results/its_estimates.csv")
  
  ## ---- 4. event-study figure (data-driven) ------------------------
  make_event <- function(panel, start, win = 6, tag) {
    panel <- dplyr::arrange(panel, date)
    k <- which.min(abs(as.numeric(panel$date - as.Date(start))))
    lo <- max(1, k - win); hi <- min(nrow(panel), k + win)
    seg <- panel[lo:hi, ]
    data.frame(event_time = seq(lo, hi) - k, growth = seg$growth, episode = tag)
  }
  es <- rbind(make_event(panel, "2017-10-01", 6, "QT1 (2017)"),
              make_event(panel, "2022-04-01", 6, "QT2 (2022)"))
  gg <- ggplot(es, aes(event_time, growth, color = episode, group = episode)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 1.1) + geom_point(size = 2.4) +
    labs(title = "Real GDP growth around Fed balance-sheet runoff episodes",
         subtitle = "Event study; quarter 0 is the QT start",
         x = "Quarters relative to QT start",
         y = "Real GDP growth (%)", color = NULL) +
    theme_minimal(base_size = 13)
  ggsave("results/event_study.png", gg, width = 8, height = 4.5, dpi = 150)
  cat("\nSaved: results/its_estimates.csv, results/panel_quarterly.csv,",
      "results/event_study.png\n")
}

main()