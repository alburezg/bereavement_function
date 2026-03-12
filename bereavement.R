# bereavement() — period bereavement estimator for DemoKin
#
# For each year (and each focal age), estimates how many people in the
# population lose at least one kin of each type *in that year*, using the
# product (kike) method:
#
#   P(lose ≥1 kin k in year y | age a) = 1 - prod_i (1 - qx_i)^living_i
#
# where the product runs over all ages i of kin type k present in kin_full.
# This is a pure period estimate — no cohort accumulation.

library(tidyverse)

# Internal helper: convert matrix, data frame, or named list to canonical long form
# value_col: name for the value column in the output ("qx" or "pop")
# sex_label: if non-NULL, adds a constant 'sex' column (used for single matrix/df case)
.to_long <- function(x, value_col, sex_label = NULL) {
  if (is.matrix(x)) {
    if (is.null(rownames(x)) || is.null(colnames(x))) {
      stop(value_col, " matrix must have dimnames: rownames = ages, colnames = years.")
    }
    out <- as.data.frame(x) |>
      tibble::rownames_to_column("age") |>
      pivot_longer(-age, names_to = "year", values_to = value_col) |>
      mutate(age = as.integer(age), year = as.integer(year))
  } else if (is.data.frame(x)) {
    req <- c("year", "age", value_col)
    missing <- setdiff(req, names(x))
    if (length(missing) > 0) {
      stop(value_col, " data frame must have columns: ", paste(req, collapse = ", "),
           ". Missing: ", paste(missing, collapse = ", "))
    }
    out <- x |>
      select(year, age, any_of("sex"), all_of(value_col)) |>
      mutate(year = as.integer(year), age = as.integer(age))
  } else {
    stop(value_col, " must be a matrix, data frame, or named list(f=..., m=...).")
  }
  if (!is.null(sex_label) && !"sex" %in% names(out)) {
    out$sex <- sex_label
  }
  out
}

# Internal: parse qx or pop argument into a single canonical long data frame.
# Returns a data frame with columns: year, age, [sex,] value_col
# require_sex: if TRUE, a sex column is mandatory when two_sex = TRUE (for qx).
#              if FALSE, sex column is optional (for pop, which may be total).
.parse_input <- function(x, value_col, two_sex, require_sex = TRUE) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!all(c("f", "m") %in% names(x))) {
      stop(value_col, " list must be named list(f = ..., m = ...) with both sexes.")
    }
    f_long <- .to_long(x$f, value_col) |> mutate(sex = "f")
    m_long <- .to_long(x$m, value_col) |> mutate(sex = "m")
    out <- bind_rows(f_long, m_long)
  } else {
    out <- .to_long(x, value_col)
  }
  has_sex_col <- "sex" %in% names(out)
  if (require_sex && two_sex && !has_sex_col) {
    stop(value_col, ": kin_full has two sexes but ", value_col,
         " has no 'sex' column. Pass a named list(f=..., m=...) or",
         " a data frame with a 'sex' column.")
  }
  if (!two_sex && has_sex_col && n_distinct(out$sex) > 1) {
    stop(value_col, ": kin_full is single-sex but ", value_col,
         " contains multiple sexes. Pass a single matrix or data frame.")
  }
  out
}


bereavement <- function(
  kin_full,           # data.frame: output of DemoKin::kin() or DemoKin::kin2sex()
  qx,                 # probability of dying (all-cause or cause-specific)
  pop,                # focal population counts (by age/year of focal individual)
  output_year = NULL  # integer vector or NULL (NULL = all years in kin_full)
) {

  # ── Step 0: Validate inputs & detect mode ──────────────────────────────────

  req_cols <- c("kin", "age_kin", "age_focal", "year", "living")
  missing_cols <- setdiff(req_cols, names(kin_full))
  if (length(missing_cols) > 0) {
    stop("kin_full is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  two_sex <- "sex_kin" %in% names(kin_full) && dplyr::n_distinct(kin_full$sex_kin) > 1

  # Filter to requested years
  if (!is.null(output_year)) {
    kin_full <- filter(kin_full, year %in% output_year)
    if (nrow(kin_full) == 0) {
      stop("output_year = ", paste(output_year, collapse = ", "),
           " not found in kin_full$year. Extend kin_full to include these years first.")
    }
  }

  # ── Step 1: Normalise qx and pop to long form ──────────────────────────────

  qx_long  <- .parse_input(qx,  "qx",  two_sex, require_sex = TRUE)
  pop_long <- .parse_input(pop, "pop", two_sex, require_sex = FALSE)

  # ── Step 2: Join qx onto kin_full ──────────────────────────────────────────

  # Build rename map and join keys
  qx_join <- qx_long |> rename(age_kin = age, qx_kin = qx)
  if ("sex" %in% names(qx_join)) qx_join <- rename(qx_join, sex_kin = sex)

  join_keys_qx <- if (two_sex) c("year", "age_kin", "sex_kin") else c("year", "age_kin")

  kin_work <- left_join(kin_full, qx_join, by = join_keys_qx)

  # Error checks
  n_na <- sum(is.na(kin_work$qx_kin))
  if (n_na > 0) {
    stop(
      n_na, " NA values in qx after joining onto kin_full. ",
      "Ensure qx covers all ages 0-100 and all years present in kin_full."
    )
  }

  bad_ages <- sort(unique(kin_work$age_kin[kin_work$qx_kin == 1 & kin_work$living > 0]))
  if (length(bad_ages) > 0) {
    stop(
      "qx = 1 detected at age(s) ", paste(bad_ages, collapse = ", "),
      " with living kin > 0. This collapses the product to 0 for any focal person ",
      "with kin at those ages, producing unreliable bereavement estimates.\n",
      "Fix before calling bereavement() — e.g., replace qx[age == 100] with ",
      "1 - exp(-mx[age == 100])."
    )
  }

  # ── Step 3: Per-cell survival probability ──────────────────────────────────

  kin_work <- mutate(kin_work, p0_cell = (1 - qx_kin)^living)

  # ── Step 4: Product over all kin ages per (year, age_focal, kin) ───────────
  # P(none of kin type k die in year y for a focal person of age a)

  group_vars <- c("year", "age_focal", "kin")
  if (two_sex && "sex_focal" %in% names(kin_work)) {
    group_vars <- c(group_vars, "sex_focal")
  }

  p0 <- summarise(kin_work, p0 = prod(p0_cell), .by = all_of(group_vars))

  # ── Step 5: Join focal population ─────────────────────────────────────────

  pop_join <- pop_long |> rename(age_focal = age, pop_focal = pop)
  if ("sex" %in% names(pop_join)) pop_join <- rename(pop_join, sex_focal = sex)

  pop_join_keys <- if ("sex_focal" %in% names(pop_join)) {
    c("year", "age_focal", "sex_focal")
  } else {
    c("year", "age_focal")
  }

  p0 <- left_join(p0, pop_join, by = pop_join_keys)

  # ── Step 6: Compute bereaved and bereaved_prop ────────────────────────────

  # Total population per year (denominator for proportion)
  pop_denom <- summarise(pop_long, pop_total = sum(pop, na.rm = TRUE), .by = year)

  result <- p0 |>
    left_join(pop_denom, by = "year") |>
    mutate(
      bereaved      = (1 - p0) * pop_focal,
      bereaved_prop = bereaved / pop_total
    ) |>
    select(year, age_focal, kin, bereaved, bereaved_prop) |>
    arrange(year, kin, age_focal)

  result
}
