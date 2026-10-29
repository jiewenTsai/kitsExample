# 真實模型：全部從 kits 的示範資料（kitm03、kitm36）的語言領域建出來，不用真的資料。
# 從 simulation_20261028/ 執行：Rscript R/truth.R
# 結果：results/truth.rds（sim_kit.R 讀它）。
#
# 1 題目參數：kits 分別校準 kitm03、kitm36 的語言題；kits_linking() 連結 M36 → M03，
#   只留平移（斜率固定 1）：shift = 篩選後留下的點的 mean(to - from)。
#   真值是一個 PCM：M03 有的題（含 18 題共同題）用 M03 的值；只有 M36 有的題用 M36 的值 + shift。
# 2 題本：每個波次的標準題本（示範資料裡那個波次最常見的作答組合），每題沒作答的比例照示範資料。
# 3 到訪：每個孩子在每個波次出現的機率照示範資料（該波次的人次 / 500）。
# 4 每個波次的平均、真分數 SD：題目參數固定在 kits 的值，TAM 每個波次一個常態估出來（M36 再加 shift）。
#   古典的修正 sqrt(var(WLE) - mean(SE^2)) 在 3、6 個月是負的，所以不用。
# 5 相鄰波次的相關：示範資料 WLE 的相關，除以 sqrt(信度的乘積) 修正誤差，最多 0.95；
#   信度 = 真分數變異數（第 4 步）/ WLE 的變異數。不相鄰的波次是相鄰相關的乘積（一階馬可夫）。

library(kits)
library(dplyr)
library(tidyr)
library(purrr)

WAVES <- list(M03 = c(3, 6, 12, 18, 24, 36, 48, 60, 72), M36 = c(36, 48, 60, 72, 84, 96))
N_CHILDREN <- 500

lan <- list(M03 = kitm03$responses, M36 = kitm36$responses) |>
  map(\(r) filter(r, domain == "lan"))

to_wide <- function(r) {
  r |>
    select(visit_id, item_id, item_score) |>
    pivot_wider(names_from = item_id, values_from = item_score, names_sort = TRUE)
}
X <- map(lan, to_wide)
visits <- map(lan, \(r) distinct(r, visit_id, subj_id, cohort, wave_m))

# ---- 1 題目參數 ----

item <- imap(X, \(x, coh) suppressMessages(kits_item(x, scale = paste("lan", coh))))
link <- kits_linking(item$M36, item$M03)
kept <- filter(as_tibble(link$points), !dropped)
shift <- mean(kept$to - kept$from)

steps_of <- function(it) select(as_tibble(it$thresholds), item_id, item_score, beta)
only36 <- setdiff(unique(item$M36$thresholds$item_id), unique(item$M03$thresholds$item_id))
betas <- bind_rows(
  mutate(steps_of(item$M03), source = "M03"),
  steps_of(item$M36) |> filter(item_id %in% only36) |> mutate(beta = beta + shift, source = "M36")
)
common <- intersect(unique(item$M03$thresholds$item_id), unique(item$M36$thresholds$item_id))

# ---- 2 題本 ----

forms <- imap(lan, \(r, coh) {
  sets <- r |>
    summarise(set = paste(sort(item_id), collapse = " "), .by = c(visit_id, wave_m))
  standard <- sets |>
    count(wave_m, set) |>
    slice_max(n, n = 1, by = wave_m, with_ties = FALSE)
  answered <- r |>
    count(wave_m, item_id, name = "n_item") |>
    left_join(count(sets, wave_m, name = "n_wave"), by = "wave_m")
  standard |>
    mutate(item_id = strsplit(set, " ")) |>
    unnest(item_id) |>
    left_join(answered, by = c("wave_m", "item_id")) |>
    transmute(cohort = coh, wave_m, item_id, p_item = n_item / n_wave)
}) |>
  list_rbind()

# ---- 3 到訪 ----

attend <- imap(visits, \(v, coh) {
  n_kids <- n_distinct(v$subj_id)
  v |>
    count(wave_m) |>
    transmute(cohort = coh, wave_m, p_wave = n / n_kids, n_kids = n_kids)
}) |>
  list_rbind()

# ---- 4 每個波次的平均、真分數 SD ----

GRID <- seq(-12, 12, length.out = 81)
wave_normals <- function(x, it, v) {
  R <- x |>
    tibble::column_to_rownames("visit_id") |>
    as.matrix()
  wave <- v$wave_m[match(rownames(R), v$visit_id)]
  stopifnot(identical(colnames(R), unique(it$thresholds$item_id)))
  fixed <- cbind(seq_len(nrow(it$thresholds)), it$thresholds$beta)
  fit <- TAM::tam.mml(R, irtmodel = "PCM", group = wave, xsi.fixed = fixed, verbose = FALSE,
                      control = list(nodes = GRID, acceleration = "Ramsay", maxiter = 5000))
  ws <- sort(unique(wave))
  tibble(wave_m = ws, mean = fit$beta[, 1], sd = sqrt(fit$variance[match(ws, wave)]))
}
normals <- pmap(list(X, item, visits), wave_normals)

# 對照：WLE 的平均、變異數、平均誤差變異數
subj <- map2(X, item, \(x, it) as_tibble(kits_subj(x, it)))
wle_moments <- imap(subj, \(s, coh) {
  s |>
    left_join(visits[[coh]], by = "visit_id") |>
    summarise(wle_mean = mean(theta), wle_var = var(theta), err_var = mean(se^2), .by = wave_m) |>
    mutate(cohort = coh)
}) |>
  list_rbind()

waves <- imap(normals, \(nw, coh) mutate(nw, cohort = coh)) |>
  list_rbind() |>
  mutate(mean = if_else(cohort == "M36", mean + shift, mean)) |>
  left_join(wle_moments, by = c("cohort", "wave_m")) |>
  mutate(rel = pmin(sd^2 / wle_var, 1)) |>
  left_join(attend, by = c("cohort", "wave_m")) |>
  select(cohort, wave_m, mean, sd, rel, p_wave, wle_mean, wle_var, err_var)

# ---- 5 相鄰波次的相關 ----

adjacent_cor <- function(coh) {
  ws <- WAVES[[coh]]
  w <- subj[[coh]] |>
    left_join(visits[[coh]], by = "visit_id") |>
    select(subj_id, wave_m, theta) |>
    pivot_wider(names_from = wave_m, values_from = theta, names_prefix = "w")
  rel <- waves |> filter(cohort == coh)
  tibble(cohort = coh, w1 = head(ws, -1), w2 = tail(ws, -1)) |>
    mutate(r_obs = map2_dbl(w1, w2, \(a, b) {
      cor(w[[paste0("w", a)]], w[[paste0("w", b)]], use = "complete.obs")
    }),
    rel1 = rel$rel[match(w1, rel$wave_m)],
    rel2 = rel$rel[match(w2, rel$wave_m)],
    r = pmin(r_obs / sqrt(rel1 * rel2), 0.95))
}
cors <- map(names(WAVES), adjacent_cor) |>
  list_rbind()

truth <- list(betas = betas, common = common, shift = shift,
              link_demo = c(intercept = link$intercept, slope = link$slope, n = link$n,
                            dropped = length(link$dropped)),
              forms = forms, waves = waves, cors = cors, WAVES = WAVES, N_CHILDREN = N_CHILDREN,
              made = list(kits = as.character(packageVersion("kits")),
                          TAM = as.character(packageVersion("TAM")), date = Sys.time()))
saveRDS(truth, "results/truth.rds")

print(waves, width = 120)
print(cors)
cat("shift =", round(shift, 4), " link:", round(link$intercept, 4), round(link$slope, 4), "\n")
