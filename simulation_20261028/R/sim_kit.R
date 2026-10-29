# SimDesign 模擬：四種估法在已知真值下，能不能找回 KIT 語言量尺（到世代內部連結為止）。
# 從 simulation_20261028/ 執行（先跑 R/truth.R）：
#   SIM_REPS=1   Rscript R/sim_kit.R   # 試跑：每個情境 1 次
#   SIM_REPS=100 Rscript R/sim_kit.R   # 正式
# 結果：results/sim_kit_<reps>reps.rds（runSimulation() 的結果，含每次複製的結果）、results/sessionInfo.txt
#
# 四種估法（同一份模擬資料都跑，所以是成對的比較）：
#   1 CML + WLE：kits_item()、kits_subj()
#   2 MML 一個常態 + WLE：TAM::tam.mml(irtmodel = "PCM")、tam.wle()
#   3 MML 每個波次一個常態 + WLE：tam.mml(..., group = wave)、tam.wle()
#   4 MML 每個波次一個常態 + EAP：和 3 同一個校準，用 $person$EAP、SD.EAP
# 每次複製：兩個世代分開校準 → 原點對齊（M03）→ 同一套連結（M36 → M03）→ 和真值比。

suppressPackageStartupMessages({
  library(SimDesign)
  library(kits)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

REPS <- as.integer(Sys.getenv("SIM_REPS", "1"))
NCORES <- as.integer(Sys.getenv("SIM_CORES", "8"))

# ---- 設定 ----

# TAM 的積分點：比較教學（第 2.1 節）用 -12 到 12、間距 0.3（81 點）。
# 每個波次一個常態時，TAM 把第一個波次（3 個月）的平均定成 0，72 個月的平均在 +9.3；
# 情境 C 的 SD 放大 2 倍，上尾（+3 SD）到 +14。所以上界放到 15，間距不變：-12 到 15，91 點。
NODES <- seq(-12, 15, by = 0.3)
MAXITER <- 5000

SCENARIOS <- tribble(
  ~scenario, ~dist,     ~sd_factor, ~model, ~label,
  "A",       "normal",  1,          "PCM",  "每個波次常態，PCM",
  "B1",      "skew",    1,          "PCM",  "左偏（偏態 -1，下尾長）",
  "B2",      "mixture", 1,          "PCM",  "混合：90% 常態 + 10% 遲緩（-2 SD）",
  "C",       "normal",  2,          "PCM",  "SD x 2：地板和天花板",
  "D",       "normal",  1,          "GPCM", "資料是 GPCM（斜率 lognormal(0, 0.3)），配 PCM"
)

CALIBS <- c(CML = "CML", MML1 = "MML 一個常態", MMLw = "MML 每波一個常態")
METHODS <- c(m1 = "1 CML + WLE", m2 = "2 MML 一個常態 + WLE",
             m3 = "3 MML 每波常態 + WLE", m4 = "4 MML 每波常態 + EAP")

# ---- 產生資料 ----

# 一個孩子各波次的相關：相鄰波次的相關（truth.R），不相鄰的是中間相鄰相關的乘積
wave_cor <- function(cors, coh, ws) {
  adj <- cors$r[cors$cohort == coh]
  R <- diag(length(ws))
  for (j in seq_along(ws)[-1]) {
    for (i in seq_len(j - 1)) R[i, j] <- R[j, i] <- prod(adj[i:(j - 1)])
  }
  R
}

# 標準化的分佈（平均 0、SD 1），經過常態的 copula，所以孩子在波次之間的名次相關不變
#   skew：反轉的標準化 gamma（shape 4，偏態 -1，下尾長）
#   mixture：10% 的孩子「遲緩」，每個波次都低 2 SD，再縮放回平均 0、SD 1
shape_z <- function(Z, dist) {
  if (dist == "skew") {
    return(-(qgamma(pnorm(Z), shape = 4) - 4) / 2)
  }
  if (dist == "mixture") {
    delayed <- runif(nrow(Z)) < 0.1
    return((Z - 2 * delayed + 0.2) / sqrt(1.36))
  }
  Z
}

# GPCM（PCM 是 a = 1）：P(X = k) 正比於 exp(a * sum_{j <= k} (theta - b_j))
draw_item <- function(theta, b, a) {
  m <- length(b)
  eta <- a * (outer(theta, 0:m) - matrix(c(0, cumsum(b)), length(theta), m + 1, byrow = TRUE))
  p <- exp(eta - apply(eta, 1, max))
  cum <- (p / rowSums(p)) %*% upper.tri(diag(m + 1), diag = TRUE)
  as.integer(rowSums(runif(length(theta)) > cum))
}

sim_cohort <- function(coh, sc, slopes, tr) {
  ws <- tr$WAVES[[coh]]
  n <- tr$N_CHILDREN
  wv <- tr$waves[tr$waves$cohort == coh, ]
  Z <- matrix(rnorm(n * length(ws)), n) %*% chol(wave_cor(tr$cors, coh, ws))
  Z <- shape_z(Z, sc$dist)
  kids <- tibble(subj_id = rep(sprintf("%s_%03d", coh, seq_len(n)), length(ws)),
                 wave_m = rep(ws, each = n),
                 theta = rep(wv$mean, each = n) + sc$sd_factor * rep(wv$sd, each = n) * as.vector(Z),
                 p_wave = rep(wv$p_wave, each = n)) |>
    filter(runif(n()) < p_wave) |>
    mutate(visit_id = paste0(subj_id, "_m", wave_m), cohort = coh, .before = 1) |>
    select(-p_wave)
  form <- tr$forms[tr$forms$cohort == coh, ]
  items <- sort(unique(form$item_id))
  resp <- matrix(NA_integer_, nrow(kids), length(items), dimnames = list(kids$visit_id, items))
  for (w in ws) {
    rows <- which(kids$wave_m == w)
    fw <- form[form$wave_m == w, ]
    for (k in seq_len(nrow(fw))) {
      i <- fw$item_id[k]
      x <- draw_item(kids$theta[rows], tr$betas$beta[tr$betas$item_id == i], slopes[[i]])
      x[runif(length(rows)) > fw$p_item[k]] <- NA
      resp[rows, i] <- x
    }
  }
  list(resp = resp, persons = kids)
}

# CML 只用原始分不是 0、也不是滿分的人次；在這些人次裡，每題的分數要從 0 到最高分都有人
# （kits 的條件是每一步都有人達到、有人沒達到；分數中間有缺時步驟參數也估不出來）。
# 不合的題在那一次複製、那個世代裡不放進校準（四種估法都一樣），記下來。
# 拿掉一題會改變誰是 0 分或滿分，所以重複檢查到沒有為止。
# 例如共同題 lanc01 用 M03 的值時對 M36 太容易，M36 幾乎沒有人得 0、1 分。
gap_items <- function(resp, betas) {
  steps <- table(betas$item_id)
  top <- as.vector((!is.na(resp)) %*% as.vector(steps[colnames(resp)]))
  raw <- rowSums(resp, na.rm = TRUE)
  inner <- resp[raw > 0 & raw < top, , drop = FALSE]
  colnames(resp)[map_lgl(colnames(resp), \(i) !all(0:steps[[i]] %in% inner[, i]))]
}

drop_gaps <- function(d, betas) {
  d$dropped <- character()
  repeat {
    gap <- gap_items(d$resp, betas)
    if (length(gap) == 0) break
    d$dropped <- c(d$dropped, gap)
    d$resp <- d$resp[, !colnames(d$resp) %in% gap, drop = FALSE]
  }
  d
}

Generate <- function(condition, fixed_objects) {
  tr <- fixed_objects$truth
  sc <- fixed_objects$scenarios[fixed_objects$scenarios$scenario == condition$scenario, ]
  items <- unique(tr$betas$item_id)
  a <- if (sc$model == "GPCM") rlnorm(length(items), 0, 0.3) else rep(1, length(items))
  slopes <- set_names(a, items)
  map(c(M03 = "M03", M36 = "M36"), \(coh) drop_gaps(sim_cohort(coh, sc, slopes, tr), tr$betas))
}

# ---- 估計 ----

fit_cml <- function(d) {
  item <- suppressMessages(kits_item(d$resp, scale = paste("lan", d$persons$cohort[1])))
  subj <- as_tibble(suppressMessages(kits_subj(d$resp, item)))
  o <- item$fit$optim
  list(steps = select(as_tibble(item$thresholds), item_id, item_score, beta),
       persons = tibble(visit_id = subj$visit_id, wle = subj$theta, wle_se = subj$se,
                        eap = NA_real_, eap_sd = NA_real_),
       iter = NA_integer_,
       converged = o$convergence == 0 || (o$max_gradient <= 1e-4 && isTRUE(o$hessian_pd)))
}

fit_tam <- function(d, group) {
  wave <- d$persons$wave_m[match(rownames(d$resp), d$persons$visit_id)]
  fit <- TAM::tam.mml(d$resp, irtmodel = "PCM", pid = rownames(d$resp),
                      group = if (group) wave else NULL, verbose = FALSE,
                      control = list(nodes = NODES, acceleration = "Ramsay", maxiter = MAXITER))
  w <- TAM::tam.wle(fit, progress = FALSE)
  steps <- tibble(name = rownames(fit$xsi), beta = fit$xsi$xsi) |>
    separate_wider_delim(name, delim = "_Cat", names = c("item_id", "item_score")) |>
    mutate(item_score = as.integer(item_score))
  persons <- tibble(visit_id = as.character(w$pid), wle = w$theta, wle_se = w$error) |>
    left_join(tibble(visit_id = as.character(fit$person$pid), eap = fit$person$EAP,
                     eap_sd = fit$person$SD.EAP), by = "visit_id")
  list(steps = steps, persons = persons, iter = fit$iter, converged = fit$iter < MAXITER)
}

# 一個估法、兩個世代分開校準，記下時間
calibrate <- function(dat, fun, ...) {
  sec <- system.time(fits <- map(c(M03 = "M03", M36 = "M36"), \(coh) fun(dat[[coh]], ...)))
  fits$seconds <- sec[["elapsed"]]
  fits
}

# Thurstonian 門檻（得到 k 分以上的機率 = .5 的 theta）：每個估法都用 kits 同一個函數算
thurstone <- function(steps) {
  steps |>
    mutate(thr = map_dbl(item_score, \(k) kits:::step_threshold(beta, k)[["thr"]]), .by = item_id)
}

# 連結 M36 → M03：每個估法都一樣，共同題每一步的門檻當點，kits_link_screen()、kits_link_line()
link_cohorts <- function(steps36, steps03, common) {
  from <- thurstone(filter(steps36, item_id %in% common))
  to <- thurstone(filter(steps03, item_id %in% common))
  points <- inner_join(select(from, item_id, item_score, from = thr),
                       select(to, item_id, item_score, to = thr), by = c("item_id", "item_score")) |>
    transmute(id = paste0(item_id, ":", item_score), item = item_id, from, to) |>
    as.data.frame()
  line <- points |>
    kits_link_screen(from_scale = "lan M36", to_scale = "lan M03") |>
    kits_link_line()
  list(intercept = line$intercept, slope = line$slope, dropped = length(line$dropped))
}

# 放到 M03 的量尺上
#   M03：平移一次，讓 M03 步驟參數的平均等於真值的平均（每個估法同一個規則）
#   M36：先置中（步驟參數平均 0），再只用連結換到 M03；不直接對齊真值
true_mean <- function(steps, betas) {
  mean(semi_join(betas, steps, by = c("item_id", "item_score"))$beta)
}

place <- function(fit, tr) {
  off03 <- true_mean(fit$M03$steps, tr$betas) - mean(fit$M03$steps$beta)
  off36 <- -mean(fit$M36$steps$beta)
  steps03 <- mutate(fit$M03$steps, beta = beta + off03)
  steps36 <- mutate(fit$M36$steps, beta = beta + off36)
  ln <- link_cohorts(steps36, steps03, tr$common)
  ln$true_intercept <- true_mean(fit$M36$steps, tr$betas)   # 置中的 M36 量尺 -> M03：真的截距，斜率 1
  to03 <- function(x) ln$intercept + ln$slope * (x + off36)
  persons <- bind_rows(
    mutate(fit$M03$persons, cohort = "M03", across(c(wle, eap), \(x) x + off03)),
    mutate(fit$M36$persons, cohort = "M36", across(c(wle, eap), to03),
           across(c(wle_se, eap_sd), \(x) ln$slope * x))
  )
  items <- bind_rows(
    mutate(steps03, cohort = "M03"),
    mutate(steps36, cohort = "M36", beta = ln$intercept + ln$slope * beta)
  )
  list(items = items, persons = persons, link = ln)
}

# ---- 每次複製的結果 ----

q <- function(x, p) unname(quantile(x, p, type = 7))
group_stats <- function(x) {
  c(mean = mean(x), sd = sd(x), p10 = q(x, 0.1), p50 = q(x, 0.5), p90 = q(x, 0.9))
}

person_outcomes <- function(est, edges) {
  by_method <- est |>
    summarise(n = n(), bias = bias(est, theta), rmse = RMSE(est, theta),
              coverage = ECR(cbind(est - 1.96 * se, est + 1.96 * se), theta),
              .by = c(method, cohort))
  bins <- est |>
    mutate(bin = findInterval(theta, edges, all.inside = TRUE)) |>
    summarise(n = n(), sum_err = sum(est - theta), sum_true = sum(theta),
              .by = c(method, cohort, bin))
  groups <- est |>
    reframe(stat = names(group_stats(est)), est = unname(group_stats(est)),
            true = unname(group_stats(theta)),
            .by = c(method, cohort, wave_m))
  list(persons = by_method, bins = bins, groups = groups)
}

floor_ceiling <- function(dat, betas) {
  steps <- table(betas$item_id)
  map(c("M03", "M36"), \(coh) {
    r <- dat[[coh]]$resp
    top <- (!is.na(r)) %*% as.vector(steps[colnames(r)])
    tibble(visit_id = rownames(r), raw = rowSums(r, na.rm = TRUE), top = as.vector(top)) |>
      left_join(dat[[coh]]$persons, by = "visit_id") |>
      summarise(n = n(), floor = mean(raw == 0), ceiling = mean(raw == top), .by = c(cohort, wave_m))
  }) |>
    list_rbind()
}

Analyse <- function(condition, dat, fixed_objects) {
  tr <- fixed_objects$truth
  fits <- list(CML = calibrate(dat, fit_cml),
               MML1 = calibrate(dat, fit_tam, group = FALSE),
               MMLw = calibrate(dat, fit_tam, group = TRUE))
  placed <- map(fits, \(f) place(f, tr))
  truth_p <- bind_rows(dat$M03$persons, dat$M36$persons) |>
    select(visit_id, wave_m, theta)
  pick <- function(calib, method, value, se) {
    placed[[calib]]$persons |>
      transmute(method = METHODS[[method]], visit_id, cohort, est = .data[[value]], se = .data[[se]])
  }
  est <- bind_rows(pick("CML", "m1", "wle", "wle_se"), pick("MML1", "m2", "wle", "wle_se"),
                   pick("MMLw", "m3", "wle", "wle_se"), pick("MMLw", "m4", "eap", "eap_sd")) |>
    left_join(truth_p, by = "visit_id")
  out <- person_outcomes(est, fixed_objects$edges[[condition$scenario]])
  out$items <- imap(placed, \(p, k) mutate(p$items, calib = k)) |>
    list_rbind() |>
    rename(est = beta)
  out$link <- imap(placed, \(p, k) as_tibble(c(list(calib = k), p$link))) |>
    list_rbind()
  out$run <- imap(fits, \(f, k) tibble(calib = k, seconds = f$seconds)) |>
    list_rbind()
  out$conv <- imap(fits, \(f, k) tibble(calib = k, cohort = c("M03", "M36"),
                                         iter = c(f$M03$iter, f$M36$iter),
                                         converged = c(f$M03$converged, f$M36$converged))) |>
    list_rbind()
  out$floor <- floor_ceiling(dat, tr$betas)
  out$dropped <- tibble(cohort = c("M03", "M36"), items = map_chr(dat, \(d) paste(d$dropped, collapse = " ")))
  out
}

# ---- 彙整：bias()、RMSE()、ECR() ----

stack <- function(results, part) {
  map(seq_along(results), \(r) mutate(as_tibble(results[[r]][[part]]), rep = r)) |>
    list_rbind()
}

Summarise <- function(condition, results, fixed_objects) {
  tr <- fixed_objects$truth
  steps <- stack(results, "items") |>
    left_join(select(tr$betas, item_id, item_score, true = beta), by = c("item_id", "item_score")) |>
    summarise(true = first(true), bias = bias(est, true), rmse = RMSE(est, true),
              .by = c(calib, cohort, item_id, item_score))
  items <- steps |>
    summarise(steps = n(), mean_bias = mean(bias), mean_abs_bias = mean(abs(bias)),
              mean_rmse = mean(rmse), .by = c(calib, cohort))
  link <- stack(results, "link") |>
    summarise(intercept_bias = bias(intercept, true_intercept),
              intercept_rmse = RMSE(intercept, true_intercept),
              slope_bias = bias(slope, 1), slope_rmse = RMSE(slope, 1),
              mean_dropped = mean(dropped), true_intercept = mean(true_intercept), .by = calib)
  persons <- stack(results, "persons") |>
    summarise(n = mean(n), bias = mean(bias), rmse = sqrt(mean(rmse^2)),
              coverage = mean(coverage), .by = c(method, cohort))
  bins <- stack(results, "bins") |>
    summarise(n = sum(n), bias = sum(sum_err) / sum(n), theta = sum(sum_true) / sum(n),
              .by = c(method, cohort, bin))
  groups <- stack(results, "groups") |>
    summarise(bias = bias(est, true), rmse = RMSE(est, true), true = mean(true),
              .by = c(method, cohort, wave_m, stat))
  floor <- stack(results, "floor") |>
    summarise(n = mean(n), floor = mean(floor), ceiling = mean(ceiling), .by = c(cohort, wave_m))
  run <- stack(results, "run") |>
    summarise(max_seconds = max(seconds), seconds = mean(seconds), .by = calib)
  conv <- stack(results, "conv") |>
    summarise(fits = n(), not_converged = sum(!converged), mean_iter = mean(iter),
              max_iter = max(iter), .by = c(calib, cohort))
  dropped <- stack(results, "dropped") |>
    separate_longer_delim(items, " ") |>
    filter(items != "") |>
    count(cohort, item_id = items, name = "reps_dropped") |>
    mutate(reps = length(results))
  list(items = items, steps = steps, link = link, persons = persons, bins = bins,
       groups = groups, floor = floor, run = run, conv = conv, dropped = dropped)
}

# ---- 執行 ----

truth <- readRDS("results/truth.rds")

# 依真 theta 分 10 組：每個情境固定的等寬區間，涵蓋兩個世代所有波次的 平均 ± 2.5 SD（超出的歸到兩端）
edges_for <- function(k) {
  w <- truth$waves
  seq(min(w$mean - 2.5 * k * w$sd), max(w$mean + 2.5 * k * w$sd), length.out = 11)
}
edges <- map(set_names(SCENARIOS$sd_factor, SCENARIOS$scenario), edges_for)

if (sys.nframe() == 0) {
  Design <- createDesign(scenario = SCENARIOS$scenario)
  fixed <- list(truth = truth, scenarios = SCENARIOS, edges = edges)
  started <- Sys.time()
  res <- runSimulation(Design, replications = REPS, generate = Generate, analyse = Analyse,
                       summarise = Summarise, fixed_objects = fixed,
                       packages = c("kits", "TAM", "dplyr", "tidyr", "purrr"),
                       parallel = REPS > 1, ncores = NCORES,
                       seed = 20261028 + seq_len(nrow(Design)),
                       store_results = TRUE, save = REPS > 1,
                       save_details = list(out_rootdir = "results"), verbose = FALSE)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  attr(res, "run_minutes") <- minutes
  attr(res, "settings") <- list(nodes = range(NODES), n_nodes = length(NODES), maxiter = MAXITER,
                                ncores = NCORES, reps = REPS)
  saveRDS(res, sprintf("results/sim_kit_%dreps.rds", REPS))
  writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
  cat("minutes:", round(minutes, 1), "\n")
}
