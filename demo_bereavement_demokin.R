# =============================================================================
# Demo: Period Bereavement Estimation with DemoKin
# =============================================================================
#
# PURPOSE
# -------
# This script demonstrates how to use the bereavement() function to estimate
# how many people in a population lose a close relative in a given year.
# We call this a "period bereavement estimate": a snapshot of kinship loss
# for everyone alive in a specific calendar year.
#
# METHOD (the "product" or "Kike" method)
# ----------------------------------------
# For a focal woman of age a in year y, the probability of losing at least
# one relative of type k (e.g. a parent) in that year is:
#
#   P(lose ≥1 kin k | age a, year y)
#     = 1 - prod_i [ (1 - qx_i)^living_i ]
#
# where the product runs over all ages i of kin type k:
#   - living_i = expected number of living relatives of type k aged i
#   - qx_i     = probability that a person aged i dies this year
#
# Intuitively: each living relative faces a small independent risk of dying.
# We compute the probability that at least one of them dies.
#
# Multiplying P(lose ≥1 kin k | age a) by the number of women aged a in year y
# gives the expected count of bereaved women of age a. Summing over all ages
# gives the total bereaved count for kin type k in year y.
#
# INPUTS
# ------
# 1. kin_full  — kinship table from DemoKin::kin2sex() or DemoKin::kin()
#                Contains, for each (year, focal age, kin type, kin age, kin sex),
#                the expected number of living and dead relatives.
#
# 2. qx        — probability of dying by age, sex, and year.
#                Can be all-cause (as here) or cause-specific.
#
# 3. pop       — population counts by age and year, for the focal individual.
#                Used to convert per-capita probabilities to head counts.
#
# DATA
# ----
# All data come from the DemoKin package (no external files needed):
#   swe_px   — Swedish female survival probabilities  (ages 0-100 x 1900-2018)
#   swe_asfr — Swedish female age-specific fertility  (ages 0-100 x 1900-2018)
#   swe_pop  — Swedish female population counts       (ages 0-100 x 1900-2018)
#
# Male rates are synthetic (following the DemoKin two-sex vignette):
#   pm = pf^1.5           — lower male survival
#   fm = ff shifted 5 yrs — later peak age at fathering
#
# REFERENCE
# ---------
# Caswell, H. (2022). Formal demography of kinship IV: Two-sex models.
# Demographic Research, 47, 359-396.
# =============================================================================

rm(list = ls())

library(DemoKin)
library(tidyverse)

source("function/bereavement.R")

# =============================================================================
# 1.  LOAD DEMOGRAPHIC RATES
# =============================================================================
# swe_px and swe_asfr are matrices bundled with DemoKin.
# Rows = ages 0-100, Columns = years 1900-2018.

pf <- swe_px    # female annual survival probability: P(alive at age a+1 | alive at age a)
ff <- swe_asfr  # female age-specific fertility rate (live births per woman per year)

n_ages <- nrow(pf)   # 101 (ages 0 to 100)
n_yrs  <- ncol(pf)   # 119


# --- Synthetic male rates ----------------------------------------------------
# DemoKin's built-in data are female-only. We create plausible male rates.
# In a real study, use observed male rates if available.

# Male survival: raised to power 1.5 → lower survival at all ages (male excess mortality)
pm <- pf ^ 1.5
dimnames(pm) <- dimnames(pf)

# Male fertility: shift the female ASFR 5 years to the right → later peak age at fathering.
# The top 5 rows of ff are dropped and replaced with zeros at the bottom.
fm <- rbind(matrix(0, 5, n_yrs), ff[-((n_ages - 4):n_ages), ])
dimnames(fm) <- dimnames(ff)

# =============================================================================
# 2.  COMPUTE KINSHIP SURFACE  (with caching)
# =============================================================================
# kin2sex() runs the two-sex time-varying kinship model.
# It returns, for every (year, focal age, kin type, kin age, kin sex),
# the expected number of living and dead relatives — the "kinship surface".
#
# time_invariant = FALSE: uses demographic rates that change over time
#   (important for capturing mortality decline, fertility transition, etc.)
# sex_focal = "f": we track kinship from the perspective of a woman
# output_kin: which relationship types to compute (saves time & memory)
#
# CACHE: kin2sex() is computationally intensive. We save the result to disk
# after the first run and reload it on subsequent runs.

cache_file <- "data_int/kin_full_swe_demo.rds"

if (file.exists(cache_file)) {
  kin_full <- readRDS(cache_file)
} else {
  kin_out <- kin2sex(
    pf = pf, pm = pm,
    ff = ff, fm = fm,
    sex_focal      = "f",
    time_invariant = FALSE,
    birth_female   = 1 / 2.04,   # ~49% of births are female (standard assumption)
    output_kin     = c("d", "gd", "gm", "m", "s")
  )
  kin_full <- kin_out$kin_full
  dir.create("data_int", showWarnings = FALSE)
  saveRDS(kin_full, cache_file)
}

# Quick check: dimensions and coverage
nrow(kin_full)
range(kin_full$year)
sort(unique(kin_full$kin))

# =============================================================================
# 3.  COMPUTE PROBABILITY OF DYING  (qx)
# =============================================================================
# The bereavement model requires qx: the probability that a person aged x
# dies within the year. We derive this from the survival probabilities:
#   qx = 1 - px
#
# TERMINAL AGE FIX:
# swe_px[age = 100] = 0 because everyone in the open age group (100+) is
# assumed to die within the year in this life table convention. That gives
# qx = 1, which collapses the entire product to zero and triggers an error
# in bereavement(). We avoid this by replacing px = 0 with a tiny positive
# value (1e-6), so qx stays just below 1.
#
# In a formal analysis, derive qx from mx using DemoTools::lt_single_mx()
# for a proper life-table treatment.

qx_f <- 1 - pmax(pf, 1e-6)   # female qx, terminal age capped
qx_m <- 1 - pmax(pm, 1e-6)   # male   qx, terminal age capped

dimnames(qx_f) <- dimnames(pf)
dimnames(qx_m) <- dimnames(pm)

# =============================================================================
# 4.  POPULATION DATA  (focal counts)
# =============================================================================
# swe_pop contains observed female population counts for Sweden by age and year.
# This is the denominator: we need to know how many women of each age are alive
# in each year to convert per-capita probabilities into head counts.
#
# bereaved[age a, year y, kin k] = P(lose ≥1 kin k | age a, year y) * pop[a, y]
#
# bereaved_prop = bereaved / total women in year y
#   → interpretable as "share of Swedish women who lost this kin type that year"

pop_mat <- swe_pop
dimnames(pop_mat) <- list(as.character(0:100), colnames(swe_pop))

# =============================================================================
# 5.  ESTIMATE BEREAVEMENT
# =============================================================================
# bereavement() applies the product formula to every combination of
# (year, focal age, kin type) and returns a tidy data frame.
#
# Arguments:
#   kin_full    — the kinship surface from kin2sex()
#   qx          — named list of qx matrices (f = female, m = male) because
#                 kin_full has sex-specific kin (sex_kin column)
#   pop         — single population matrix (female counts); passed without a
#                 sex column so it applies to all focal women regardless of
#                 which kin they might lose
#   output_year — NULL means "compute for every year in kin_full"
#
# Returns one row per (year, age_focal, kin):
#   bereaved      = expected number of bereaved women of this age
#   bereaved_prop = bereaved / total women in that year

result <- bereavement(
  kin_full    = kin_full,
  qx          = list(f = qx_f, m = qx_m),
  pop         = pop_mat,
  output_year = NULL
)

# =============================================================================
# 6.  AGGREGATE: total bereaved women by kin type and year
# =============================================================================
# result has one row per (year, age_focal, kin). To get the total for each
# kin type in a given year, sum bereaved over all focal ages.
# bereaved_prop is already scaled by total population, so summing it gives
# the total share of women bereaved for that kin type in that year.

summary_by_year <- result |>
  summarise(
    bereaved      = sum(bereaved),
    bereaved_prop = sum(bereaved_prop),
    .by = c(year, kin)
  )

# Add human-readable kin labels using DemoKin's built-in helper.
# rename_kin() maps codes like "m", "d", "gm" to "Mothers", "Children", etc.
summary_by_year <- rename_kin(summary_by_year, sex = "f")

# Bereaved share (% of Swedish women) in selected years
summary_by_year |>
  filter(year %in% c(1950, 1970, 1990, 2010)) |>
  mutate(pct = round(bereaved_prop * 100, 2)) |>
  select(year, kin_label, pct) |>
  pivot_wider(names_from = year, values_from = pct) |>
  arrange(kin_label) |>
  print()

# =============================================================================
# 7.  PLOT: bereaved women by kin type, faceted by year
# =============================================================================
# Style mirrors script 8 (8_sensitivity_methods_compare.R):
#   - Horizontal bar chart (coord_flip)
#   - One panel per selected year
#   - Bars show total bereaved women (summed over all focal ages)

plot_years <- c(1950, 1970, 1990, 2010)

p <- summary_by_year |>
  filter(year %in% plot_years) |>
  ggplot(aes(y = bereaved, x = kin_label, fill = kin_label)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25,
           show.legend = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
  coord_flip() +
  facet_wrap(~ year, nrow = 1) +
  labs(
    title    = "Swedish Women Who Lost a Relative — All-Cause Mortality",
    subtitle = paste0(
      "Bereaved = P(\u22651 kin loss this year) \u00d7 female population, summed over focal ages.\n",
      "Two-sex time-varying model. Female focal. Synthetic male rates. Source: DemoKin."
    ),
    x = NULL,
    y = "Women who lost this relative in the given year"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 7.5, colour = "grey40"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dir.create("output", showWarnings = FALSE)
ggsave("output/bereavement_demo_demokin.png", p, width = 14, height = 4, dpi = 300)
cat("\nPlot saved to output/bereavement_demo_demokin.png\n")
