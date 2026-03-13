# bereavement() — period bereavement estimator for DemoKin
#
# Estimates how many people in a population lose at least one relative of each
# type in a given year, using the product (Kike) method:
#
#   b(a, y, k) = P(lose ≥1 kin k in year y | focal age a)
#              = 1 - prod_{a', g} (1 - qx(a', y, g))^living(a, a', k, y, g)
#
# where the product runs over all ages a' and sexes g of kin type k.
# This is a pure period estimate — no cohort accumulation.
#
# Works with output from DemoKin::kin()  (one-sex, time-invariant or -varying)
#                   and DemoKin::kin2sex() (two-sex, time-invariant or -varying).
# Time-invariant models have year = NA in kin_full; the function detects this
# and drops year from all join keys automatically.

library(tidyverse)


# ── Internal helpers ───────────────────────────────────────────────────────────

# Convert matrix, numeric vector, or data frame to canonical long form.
# Numeric vectors are treated as age-indexed (age 0, 1, ..., n-1), no year.
.to_long <- function(x, value_col, sex_label = NULL) {

  # Numeric vector → single-column matrix with dummy year "1"
  if (is.numeric(x) && !is.matrix(x) && !is.data.frame(x)) {
    rn <- if (!is.null(names(x))) names(x) else as.character(seq_along(x) - 1L)
    x  <- matrix(x, ncol = 1, dimnames = list(rn, "1"))
  }

  if (is.matrix(x)) {
    if (is.null(rownames(x)) || is.null(colnames(x))) {
      stop(value_col, " matrix must have dimnames (rownames = ages, colnames = years).")
    }
    out <- as.data.frame(x) |>
      tibble::rownames_to_column("age") |>
      pivot_longer(-age, names_to = "year", values_to = value_col) |>
      mutate(age = as.integer(age), year = as.integer(year))

  } else if (is.data.frame(x)) {
    req <- c("age", value_col)
    missing_cols <- setdiff(req, names(x))
    if (length(missing_cols) > 0) {
      stop(value_col, " data frame must have columns 'age' and '", value_col,
           "'. Missing: ", paste(missing_cols, collapse = ", "))
    }
    out <- x |>
      select(any_of("year"), age, any_of("sex"), all_of(value_col)) |>
      mutate(age = as.integer(age))
    if ("year" %in% names(out)) out <- mutate(out, year = as.integer(year))

  } else {
    stop(value_col, " must be a numeric vector, matrix, data frame, or ",
         "named list(f = ..., m = ...).")
  }

  if (!is.null(sex_label) && !"sex" %in% names(out)) out$sex <- sex_label
  out
}

# Parse qx or pop argument into a single canonical long data frame.
# In two-sex mode, require_sex = TRUE forces a sex column to be present.
.parse_input <- function(x, value_col, two_sex, require_sex) {

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

  has_sex <- "sex" %in% names(out)

  if (require_sex && two_sex && !has_sex) {
    stop(value_col, ": kin_full has two sexes (sex_kin column present) but ",
         value_col, " has no sex column. Pass a named list(f = ..., m = ...).")
  }
  if (!two_sex && has_sex && n_distinct(out$sex) > 1) {
    stop(value_col, ": kin_full is one-sex but ", value_col,
         " supplies two sexes. Pass a single matrix or data frame.")
  }
  out
}


# ── Main function ──────────────────────────────────────────────────────────────

bereavement <- function(
    kin_full,           # data.frame from DemoKin::kin() or DemoKin::kin2sex()
    qx,                 # probability of dying: named list(f,m) in two-sex mode
    pop,                # population counts:    named list(f,m) in two-sex mode
    output_year = NULL, # integer vector of years to compute, or NULL (= all)
    group_kin   = NULL  # named list of kin groupings, e.g. list(nuclear=c("m","d","s"))
) {

  # ── Step 0: Validate inputs & detect model mode ─────────────────────────────

  req_cols <- c("kin", "age_kin", "age_focal", "year", "living")
  missing_cols <- setdiff(req_cols, names(kin_full))
  if (length(missing_cols) > 0) {
    stop("kin_full is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Two-sex mode: sex_kin column present with ≥2 distinct values
  two_sex <- "sex_kin" %in% names(kin_full) && n_distinct(kin_full$sex_kin) > 1

  # Time-invariant mode: kin_full has year = NA (DemoKin convention)
  time_inv <- all(is.na(kin_full$year))

  # Filter to requested years (time-varying only; ignored in time-invariant mode)
  if (!is.null(output_year) && !time_inv) {
    kin_full <- filter(kin_full, year %in% output_year)
    if (nrow(kin_full) == 0) {
      stop("output_year = {", paste(output_year, collapse = ", "),
           "} not found in kin_full$year.")
    }
  }

  # ── Step 1: Normalise qx and pop to long form ───────────────────────────────
  # Both must be sex-specific (list(f,m)) in two-sex mode.
  # In two-sex mode, pop's f and m entries are summed to give total population.

  qx_long  <- .parse_input(qx,  "qx",  two_sex, require_sex = TRUE)
  pop_long <- .parse_input(pop, "pop", two_sex, require_sex = TRUE)

  if (two_sex && "sex" %in% names(pop_long)) {
    pop_long <- summarise(pop_long, pop = sum(pop, na.rm = TRUE),
                          .by = any_of(c("year", "age")))
  }

  # ── Step 2: Build join keys (time-invariant drops "year") ───────────────────

  if (time_inv) {
    join_keys_qx  <- if (two_sex) c("age_kin", "sex_kin") else "age_kin"
    pop_join_keys <- "age_focal"
    denom_key     <- character(0)          # denominator = scalar (sum over all ages)
  } else {
    join_keys_qx  <- if (two_sex) c("year", "age_kin", "sex_kin") else c("year", "age_kin")
    pop_join_keys <- "year"                # extended below after pop_join is built
    denom_key     <- "year"
  }

  # ── Step 3: Join qx onto kin_full ───────────────────────────────────────────

  qx_join <- qx_long |> rename(age_kin = age, qx_kin = qx)
  if ("sex" %in% names(qx_join)) qx_join <- rename(qx_join, sex_kin = sex)
  # In time-invariant mode the dummy year column must be dropped before joining,
  # otherwise it collides with kin_full$year (= NA).
  if (time_inv && "year" %in% names(qx_join)) qx_join <- select(qx_join, -year)

  kin_work <- left_join(kin_full, qx_join, by = join_keys_qx)

  n_na <- sum(is.na(kin_work$qx_kin))
  if (n_na > 0) {
    stop(n_na, " NA values in qx after joining. ",
         "Ensure qx covers ages 0-100 and all years in kin_full.")
  }
  bad_ages <- sort(unique(kin_work$age_kin[kin_work$qx_kin == 1 & kin_work$living > 0]))
  if (length(bad_ages) > 0) {
    stop("qx = 1 at age(s) ", paste(bad_ages, collapse = ", "),
         " with living kin > 0. Fix terminal age before calling bereavement():\n",
         "  e.g. qx[age == 100] <- 1 - exp(-mx[age == 99])")
  }

  # ── Step 4: Per-cell survival probability → product per (year, age_focal, kin) ─

  kin_work <- mutate(kin_work, p0_cell = (1 - qx_kin)^living)

  gvars <- unique(c("year", "age_focal", "kin",
                    if ("sex_focal" %in% names(kin_work)) "sex_focal"))

  p0 <- summarise(kin_work, p0 = prod(p0_cell), .by = all_of(gvars))

  # ── Step 4b: Grouped kin (optional) ─────────────────────────────────────────
  # For each user-defined group, multiply p0 across the constituent kin types.
  # This uses the same independence assumption as the individual calculation:
  #   P(lose ≥1 in group) = 1 - prod_{k in group} p0_k

  if (!is.null(group_kin)) {
    if (!is.list(group_kin) || is.null(names(group_kin))) {
      stop("group_kin must be a named list, e.g. list(nuclear = c('m', 'd', 's')).")
    }
    non_kin <- setdiff(gvars, "kin")

    group_p0_list <- lapply(names(group_kin), function(gname) {
      codes   <- group_kin[[gname]]
      present <- codes[codes %in% unique(p0$kin)]
      if (length(present) == 0) {
        warning("group_kin '", gname, "': none of {",
                paste(codes, collapse = ", "), "} found in kin_full. Skipping.")
        return(NULL)
      }
      if (length(present) < length(codes)) {
        warning("group_kin '", gname, "': only {",
                paste(present, collapse = ", "), "} present in kin_full.")
      }
      p0 |>
        filter(kin %in% present) |>
        summarise(p0 = prod(p0), .by = all_of(non_kin)) |>
        mutate(kin = gname)
    })

    group_p0 <- bind_rows(Filter(Negate(is.null), group_p0_list))
    if (nrow(group_p0) > 0) p0 <- bind_rows(p0, group_p0)
  }

  # ── Step 5: Join focal population ───────────────────────────────────────────

  pop_join <- pop_long |> rename(age_focal = age, pop_focal = pop)
  if ("sex" %in% names(pop_join)) pop_join <- rename(pop_join, sex_focal = sex)

  if (!time_inv) {
    pop_join_keys <- c("year", "age_focal",
                       if ("sex_focal" %in% names(pop_join)) "sex_focal")
  }
  # In time-invariant mode drop the dummy year column to avoid collision with p0$year (NA).
  if (time_inv && "year" %in% names(pop_join)) pop_join <- select(pop_join, -year)

  p0 <- left_join(p0, pop_join, by = pop_join_keys)

  # ── Step 6: Total-population denominator ─────────────────────────────────────

  if (time_inv) {
    pop_total <- sum(pop_long$pop, na.rm = TRUE)
    result <- p0 |>
      mutate(
        bereaved      = (1 - p0) * pop_focal,
        bereaved_prop = bereaved / pop_total
      )
  } else {
    pop_denom <- summarise(pop_long, pop_total = sum(pop, na.rm = TRUE),
                           .by = all_of(denom_key))
    result <- p0 |>
      left_join(pop_denom, by = denom_key) |>
      mutate(
        bereaved      = (1 - p0) * pop_focal,
        bereaved_prop = bereaved / pop_total
      )
  }

  result |>
    select(year, age_focal, kin, bereaved, bereaved_prop) |>
    arrange(year, kin, age_focal)
}
