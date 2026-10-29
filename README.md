# kitsExample

[kits](https://github.com/jiewenTsai/kits) 和 [kitsPreview](https://github.com/jiewenTsai/kitsPreview) 的教學、例子和模擬。
全部只用套件的示範資料（`kitm03`、`kitm36`）和模擬資料，沒有任何真實孩子的資料。

網頁版（GitHub Pages）：<https://jiewentsai.github.io/kitsExample/>

## 教學

| 資料夾 | 內容 | 網頁 |
|---|---|---|
| `tutorial_beginner/` | 入門：從作答到常模（中文、英文） | [中文](https://jiewentsai.github.io/kitsExample/tutorial_beginner/kits_beginner_zh.html)、[English](https://jiewentsai.github.io/kitsExample/tutorial_beginner/kits_beginner.html) |
| `tutorial_algorithms/` | 七份循序漸進的演算法教學，公式可以複製（LaTeX），每一步都和 kits 對答案 | 見下表 |
| `tutorial_models/` | kits 的模型：寫一次，導數交給 RTMB | [網頁](https://jiewentsai.github.io/kitsExample/tutorial_models/kits_models.html) |
| `tutorial_preview/` | kitsPreview 導覽：dexter、immer、cNORM 原本怎麼做 | [網頁](https://jiewentsai.github.io/kitsExample/tutorial_preview/kitsPreview_tutorial_zh.html) |
| `simulation_20261028/` | 已知真值下比較四種估法（CML + WLE、TAM 的 MML + WLE/EAP） | [報告](https://jiewentsai.github.io/kitsExample/simulation_20261028/sim_report_zh.html) |

演算法教學：

| 份 | 主題 |
|---|---|
| [01](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/01_calibration.html) | 校準：從作答到門檻（PCM、CML、ESF） |
| [02](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/02_scoring.html) | 計分：每個孩子的 WLE（Jeffreys 眾數估計）和 SE |
| [03](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/03_linking.html) | 連結與換算：世代連結、接到 D |
| [04](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/04_norms.html) | 常模：shash 模型 |
| [05](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/05_checks.html) | 檢查題目：適配、gamma、DIF、偏移、每組的通過率圖 |
| [06](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/06_scale.html) | 量尺的性質：訊息、精準度、原始分對照、垂直量尺 |
| [07](https://jiewentsai.github.io/kitsExample/tutorial_algorithms/07_subscales.html) | 分量表：共同量尺上的分量表分數、SE、成長模型和模擬 |

## 自己重新產生

每份 `.qmd` 都可以用 Quarto 產生，需要 kits 和 kitsPreview：

```r
remotes::install_github("jiewenTsai/kits")
remotes::install_github("jiewenTsai/kitsPreview")
```

```sh
quarto render tutorial_algorithms/01_calibration.qmd --to html
```

模擬要先跑 `simulation_20261028/R/truth.R`，再跑 `R/sim_kit.R`（見該資料夾的 README）。

## 之後會加上的

- 用 KIT 資料的兩份教學（kits 用在 KIT 的完整流程、kits 和 TAM 的比較）會改用示範資料重寫後再放上來。
  `tutorial_models` 裡連到 TAM 比較的連結，在那之前打不開。

## 授權

GPL-3，和 kits、kitsPreview 相同。
