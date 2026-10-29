# 模擬：四種估法在已知真值下找不找得回 KIT 語言量尺（2026-10-28）

## 做什麼

真實資料上 kits（CML + WLE）和 TAM 的數字不一樣，沒有真值就不知道誰對。
這裡用 SimDesign 在已知的真實模型下比較四種估法，一路看到兩個世代（M03、M36）的內部連結：

1. CML + WLE（kits）
2. MML 一個常態 + WLE（TAM）
3. MML 每個波次一個常態 + WLE（TAM，`group = wave`）
4. MML 每個波次一個常態 + EAP（同 3 的校準）

真實模型全部從 kits 的示範資料（`kitm03`、`kitm36`）的語言領域建出來，不用真實資料。
五個情境：A 常態、B1 左偏、B2 混合（10% 遲緩）、C SD × 2（地板和天花板）、D 資料是 GPCM。
接到 D 的外部連結不在範圍內。

## 檔案

* `R/truth.R`：真實模型（題目參數、題本、到訪、波次的平均和 SD、相鄰波次的相關）→ `results/truth.rds`
* `R/sim_kit.R`：SimDesign 的 Design、Generate、Analyse、Summarise、runSimulation
* `results/sim_kit_100reps.rds`：runSimulation() 的結果（`SUMMARISE` 欄是彙整，`SimResults()` 是每次複製的結果）
* `results/sim_kit_1reps.rds`：試跑（每個情境 1 次）
* `results/sessionInfo.txt`、`results/run_100reps.log`
* `sim_report_zh.Rmd`、`sim_report_zh.html`：報告

## 怎麼跑

從這個資料夾執行：

```sh
export OPENBLAS_NUM_THREADS=1 OMP_THREAD_LIMIT=1
Rscript R/truth.R
SIM_REPS=1   Rscript R/sim_kit.R                 # 試跑
SIM_REPS=100 SIM_CORES=8 Rscript R/sim_kit.R     # 正式
Rscript -e 'rmarkdown::render("sim_report_zh.Rmd")'
```

需要 kits 0.23.0、TAM、SimDesign 2.27。

## 時間

* 正式：每個情境 100 次複製、5 個情境、8 個核心：73 分鐘（2026-10-28；同時有另一個工作用 4 個核心）。
  每次複製約 1 分鐘，大部分花在 TAM 每個波次一個常態的校準。
* 試跑：每個情境 1 次，不平行：約 3 分鐘。
* 跑完以後修了 `Summarise()` 裡 `max_seconds` 的順序，用存下來的每次複製結果重新彙整一次（同一個檔案）；其他數字不變。
* TAM 的 2000 次校準都收斂（沒有到 `maxiter = 5000`）；情境 D 有 1 次 TAM（一個常態、M03）在 Ramsay 加速時出錯，SimDesign 丟掉那一次、重抽。
