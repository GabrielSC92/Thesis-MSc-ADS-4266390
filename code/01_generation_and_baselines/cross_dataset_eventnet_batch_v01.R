# Cross-dataset batch run: for each (seed x MR level), simulate data, fit REM
# and RHEM, and save tidy coefficient tables for cross_dataset_averaging_v02.qmd.

# run from the repo root or this folder:
#   Rscript code/thesis/01_generation_and_baselines/cross_dataset_eventnet_batch_v01.R

library(survival)
library(remify)
library(remstats)
library(remstimate)
library(remulate)

N_SEEDS       <- 30
MAIN_SEED     <- 20260423
SEEDS         <- MAIN_SEED + 0:(N_SEEDS - 1)
MR_LEVELS_PCT <- c(0, 2, 8, 16, 24)

n_actors   <- 20
n_events   <- 2000
eps_factor <- 0.10   # tie-break noise scale (see noise_reshuffle_rem_v02.qmd)

# ground-truth simulator parameters, used later to compute bias
true_params <- c(
  baseline = -5.00, inertia = 0.02, reciprocity = 0.10, otp = 0.03,
  difference_active = 0.15, same_group = 0.40
)

# locate project folders; needs to be run via Rscript so --file= is set
script_path  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir   <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

data_dir   <- file.path(project_root, "code", "data", "synthetic")
out_root   <- file.path(data_dir, "cross_dataset")
inputs_dir <- file.path(out_root, "inputs")      # eventnet input CSVs
enet_dir   <- file.path(out_root, "eventnet_out") # eventnet design matrices
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enet_dir,   recursive = TRUE, showWarnings = FALSE)

rem_csv  <- file.path(out_root, "rem_coefs_by_seed.csv")
rhem_csv <- file.path(out_root, "rhem_coefs_by_seed.csv")

# Same eventnet setup as the GUI runs in rhem_eventnet_prep_v01.qmd.
# We only swap the input/output paths per cell; everything else stays fixed.
eventnet_config <- file.path(data_dir, "eventnet_gui_statistics", "eventnet.configuration.txt")
eventnet_jar    <- normalizePath(file.path(project_root, "..", "eventnet-1.3.jar"))
obs_name        <- "MR24_COND_SENDER"  # must match the observation block in the config file


# Step 1: draw one clean baseline (0% MR) for this seed.
# same as baseline_and_mr_generation_v02.qmd
generate_baseline <- function(seed) {
  set.seed(seed)

  actor_attrs <- data.frame(
    id = 1:n_actors, time = 0,
    group = sample(c(0, 1), n_actors, replace = TRUE),
    active = rnorm(n_actors)
  )

  # bquote() pastes the actual actor_attrs table into the formula, otherwise
  # remulate cannot find it because it lives inside this function.
  effects <- eval(bquote(
    ~ baseline(-5) + inertia(0.02) + reciprocity(0.10) + otp(0.03) +
      difference(0.15, variable = "active", attr_actors = .(actor_attrs)) +
      same(0.4, variable = "group", attr_actors = .(actor_attrs))
  ))

  reh <- remulateTie(
    effects = effects, actors = 1:n_actors,
    endTime = Inf, events = n_events, initial = 0
  )
  reh$event_id <- seq_len(nrow(reh))

  long0 <- data.frame(
    event_id = reh$event_id, sender = reh$sender,
    receiver = reh$receiver, time = reh$time, mr_size = 1L
  )
  list(long0 = long0, attrs = actor_attrs)
}


# Step 2: randomly tag a fraction of events and add one extra receiver (MR v2 design)
make_mr_v02 <- function(reh_df, p_mr_event, seed) {
  set.seed(seed)
  n_events <- nrow(reh_df)
  n_to_tag <- round(p_mr_event * n_events)

  # which events become multi-receiver (none when p = 0)
  tagged_ids <- sample(reh_df$event_id, n_to_tag)

  mr_long <- data.frame()
  for (i in 1:n_events) {
    row <- reh_df[i, ]
    is_tagged <- row$event_id %in% tagged_ids
    this_size <- if (is_tagged) 2 else 1

    # original sender -> receiver row (kept for every event)
    mr_long <- rbind(mr_long, data.frame(
      event_id = row$event_id, sender = row$sender,
      receiver = row$receiver, time = row$time, mr_size = this_size
    ))

    # tagged events also get one extra receiver at the same time
    if (is_tagged) {
      extra_receiver <- sample(setdiff(1:n_actors, c(row$sender, row$receiver)), 1)
      mr_long <- rbind(mr_long, data.frame(
        event_id = row$event_id, sender = row$sender,
        receiver = extra_receiver, time = row$time, mr_size = this_size
      ))
    }
  }

  mr_long
}


# Step 3 (REM only): break timestamp ties with small random noise so remify gets a strict order
inject_noise <- function(long_df, seed) {
  in_mr <- long_df$mr_size > 1
  shared_time <- duplicated(long_df$time) | duplicated(long_df$time, fromLast = TRUE)
  needs_noise <- in_mr | shared_time

  unique_times <- sort(unique(long_df$time))
  eps_max <- min(diff(unique_times)) * eps_factor

  set.seed(seed)
  noise <- numeric(nrow(long_df))
  noise[needs_noise] <- runif(sum(needs_noise), 0, eps_max)

  long_df$time_noisy <- long_df$time + noise
  long_df <- long_df[order(long_df$time_noisy, long_df$event_id), ]

  stopifnot(sum(duplicated(long_df$time_noisy)) == 0)
  list(df = long_df, eps_max = eps_max)
}


# fit dyadic REM on the noise-adjusted edgelist (same effects as parameter recovery notebook)
fit_rem <- function(long_df, actor_attrs) {
  attrs_rs <- data.frame(
    name = as.character(actor_attrs$id), time = actor_attrs$time,
    group = actor_attrs$group, active = actor_attrs$active
  )
  edgelist <- data.frame(
    time = long_df$time_noisy,
    actor1 = long_df$sender, actor2 = long_df$receiver
  )

  reh <- remify(edgelist, directed = TRUE, ordinal = FALSE, model = "tie",
                actors = 1:n_actors)
  stats <- remstats(reh,
    tie_effects = ~ 1 + inertia() + reciprocity() + otp() +
      difference(variable = "active") + same(variable = "group"),
    attr_actors = attrs_rs)
  fit <- remstimate(reh, stats, method = "MLE")

  list(coef = coef(fit), se = sqrt(diag(vcov(fit))))
}


# Build the long-format CSV that eventnet expects
build_eventnet_csv <- function(long_df, actor_attrs, t_dummy) {
  add_actors <- data.frame(
    message.id = 0, sender.id = actor_attrs$id, receiver.id = actor_attrs$id,
    time = t_dummy, type = "add.actor", weight = 1
  )
  attr_group <- data.frame(
    message.id = 0,
    sender.id = actor_attrs$id[actor_attrs$group == 1],
    receiver.id = actor_attrs$id[actor_attrs$group == 1],
    time = t_dummy, type = "is.group", weight = 1
  )
  attr_active <- data.frame(
    message.id = 0, sender.id = actor_attrs$id, receiver.id = actor_attrs$id,
    time = t_dummy, type = "active", weight = actor_attrs$active
  )
  events <- data.frame(
    message.id = long_df$event_id, sender.id = long_df$sender,
    receiver.id = long_df$receiver, time = long_df$time,
    type = "email", weight = 1
  )
  out <- rbind(add_actors, attr_group, attr_active, events)
  out[order(out$time, out$message.id), ]
}


# Copy the GUI config and replace only the three paths that change per run.
write_eventnet_config <- function(input_dir, input_file, output_dir) {
  txt <- readLines(eventnet_config)
  txt <- sub('(<input\\.directory name=")[^"]*(")',
             paste0("\\1", normalizePath(input_dir), "\\2"), txt)
  txt <- sub('(<file name=")[^"]*(")', paste0("\\1", input_file, "\\2"), txt)
  txt <- sub('(<output\\.directory name=")[^"]*(")',
             paste0("\\1", normalizePath(output_dir), "\\2"), txt)
  cfg <- file.path(out_root, paste0(input_file, ".config.txt"))
  writeLines(txt, cfg)
  cfg
}


# Fit RHEM by running coxph on the eventnet design matrix (same as Section 7 in rhem_eventnet_prep_v01.qmd).
fit_rhem <- function(output_csv) {
  ev <- read.csv(output_csv)
  ev <- ev[ev$TYPE == "email" | ev$IS_OBSERVED == 0, ]

  fit <- coxph(
    Surv(time = rep(1, nrow(ev)), event = ev$IS_OBSERVED) ~
      s.r.sub.rep.1 + reciprocation + otp + absdiff.active + catdiff.group +
      strata(EVENT_INTERVAL),
    data = ev
  )

  cf <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  effect_map <- c(
    "s.r.sub.rep.1" = "inertia", reciprocation = "reciprocity", otp = "otp",
    absdiff.active = "difference_active", catdiff.group = "same_group"
  )
  data.frame(
    scale = "raw",
    statistic = names(cf),
    effect = effect_map[names(cf)],
    estimate = as.numeric(cf),
    se = as.numeric(se),
    # eventnet CATDIFF has opposite sign to remstats::same(); flip for comparison.
    estimate_aligned = ifelse(names(cf) == "catdiff.group", -as.numeric(cf), as.numeric(cf))
  )
}


# main loop: one pass over all (seed, MR level) cells

rem_all  <- data.frame()
rhem_all <- data.frame()

for (seed in SEEDS) {

  # Baseline is shared across MR levels within a seed.
  bl <- generate_baseline(seed)
  long0 <- bl$long0
  attrs <- bl$attrs
  t_dummy <- min(long0$time) - 1   # dummy time for attribute rows in eventnet CSV

  # Skip seeds where remulate already produced tied timestamps (rare).
  if (any(duplicated(long0$time))) next

  for (pct in MR_LEVELS_PCT) {

    # --- REM pipeline: MR inject -> tie-break noise -> MLE ---
    mr_long <- make_mr_v02(long0, pct / 100, seed + 100000 + pct)
    noisy   <- inject_noise(mr_long, seed + 200000 + pct)

    rem_fit <- fit_rem(noisy$df, attrs)
    rem_all <- rbind(rem_all, data.frame(
      seed = seed, mr_level_pct = pct,
      statistic = names(rem_fit$coef),
      estimate = as.numeric(rem_fit$coef),
      se = as.numeric(rem_fit$se),
      n_rows = nrow(noisy$df), eps_max = noisy$eps_max,
      true_value = true_params[names(rem_fit$coef)],
      bias = as.numeric(rem_fit$coef) - true_params[names(rem_fit$coef)]
    ))

    # RHEM pipeline: export CSV -> eventnet -> coxph
    # RHEM uses original timestamps (no tie-break noise).
    tag        <- paste0("s", seed, "_mr", sprintf("%02d", pct))
    input_file <- paste0(tag, "_eventnet.csv")
    write.csv(build_eventnet_csv(mr_long, attrs, t_dummy),
              file.path(inputs_dir, input_file), row.names = FALSE)

    cfg <- write_eventnet_config(inputs_dir, input_file, enet_dir)
    system2("java", c("-jar", shQuote(eventnet_jar), shQuote(cfg)),
            stdout = FALSE, stderr = FALSE)

    out_file <- file.path(enet_dir, paste0(tag, "_eventnet_", obs_name, ".csv"))
    rhem_fit <- fit_rhem(out_file)
    rhem_fit$seed <- seed
    rhem_fit$mr_level_pct <- pct
    rhem_fit$true_value <- true_params[rhem_fit$effect]
    rhem_all <- rbind(rhem_all, rhem_fit[, c(
      "seed", "mr_level_pct", "scale", "statistic", "effect",
      "estimate", "estimate_aligned", "se", "true_value"
    )])
  }
}

write.csv(rem_all,  rem_csv,  row.names = FALSE)
write.csv(rhem_all, rhem_csv, row.names = FALSE)
