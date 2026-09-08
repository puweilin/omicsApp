# 当前版本审查证据

对应 [CURRENT-VERSION-REVIEW.md](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/CURRENT-VERSION-REVIEW.md)。审查基线为 `a5ae7a310297985969f1f2c81498880d6c042644`，版本及验证边界见 `metadata.json`。

本目录属于框架审计证据，不是科研分析产物。没有复制真实研究矩阵、私人项目、安装包库或测试 RDS；真实数据补测只保留测试脚本、汇总及日志。

## 文件索引

| 文件 | 用途 |
|---|---|
| `review_probes.R` / `review-probes.log` | C03–C05、C08、C09、C12：批量参数/错误、任务负载、迟到回调、导出 |
| `review_probes_2.R` / `review-probes-2.log` | C01、C06：脚本层选择和尺度 provenance；导入暂存状态 |
| `review_replay_process.R` / `replay-process.log` | C06：归档原文件存在时，在全新 R 进程执行实际导出的脚本 |
| `review_probes_3.R` / `review-probes-3.log` | C07、C11：完整 app_server 新建后重导入；无关更新触发 QC |
| `review_probes_4.R` / `review-probes-4.log` | C02、C11：协变量整合签名失配；完整 app_server 重复整合 |
| `review_probes_5.R` / `review-probes-5.log` | C01：真实整合模块选择同一层两次 |
| `probe_fixtures.R` | 上述探针的确定性小型数据和真实分析 bundle |
| `progress_browser.R`、`progress-app/app.R`、`progress-browser-2.log` | C10：Chrome + Shiny + multisession(2) 的实际进度轮询复现 |
| `*-suite.log`、`*-totals.log`、`*-tests.csv` | 两套主测试的原日志及恢复的统计 |
| `extra-tests.log`、`*-extra.csv`、`run_extra.R` | 真实 RNA/蛋白组、可选性能预算补测 |
| `fuzz-tests.log`、`core-fuzz.csv`、`run_fuzz.R` | 可选重型 API 检查；包含一部分主套件的重复检查 |
| `package-checks.log`、`*-00check.log`、`check_packages.R` | 本机隔离安装和 R CMD check |
| `deploy-current/`、`demo-study.log` | 审查期间新增的部署脚本；本次运行了生成示例和整合核对，计时脚本仅读码 |
| `live-delta/`、`live-exports.log` | 工作区绘图改动的捕获副本、patch 和安装后导出接口验证 |
| `msigdbr-preflight.log` | 多库测试的预检查实际因错误参数而提前跳过 |
| `*-replay*.R` | 导出函数当时生成的脚本文本；输入文件未随证据复制 |
| `metadata.json`、`sha256.json` | 环境、版本、范围和证据文件校验清单 |

## 如何复现

依赖需与包 DESCRIPTION 相符，浏览器探针还需要 Chrome、shinytest2 和 callr。建议使用独立临时目录，避免在开发者正在修改的工作区运行，也避免把探针输出误写入产品仓库。

```sh
repo=/Users/puweilin/SciProject/CHISSS/omicsApp
audit_dir=$(mktemp -d /tmp/omicsapp-review-replay.XXXXXX)
git -C "$repo" archive a5ae7a310297985969f1f2c81498880d6c042644 | tar -x -C "$audit_dir"
mkdir -p "$audit_dir/audit-results"
cp -R "$repo/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/." "$audit_dir/audit-results/"
cd "$audit_dir"
mkdir -p audit-results/library
R CMD INSTALL --library=audit-results/library packages/omicsCore
Rscript audit-results/review_probes.R
Rscript audit-results/review_probes_2.R
Rscript audit-results/review_probes_3.R
Rscript audit-results/review_probes_4.R
Rscript audit-results/review_probes_5.R
Rscript audit-results/review_replay_process.R
Rscript audit-results/progress_browser.R
```

每个探针调用真实产品函数。需要保持可控的地方才模拟调度或后端异常：C08 控制迟到任务的返回顺序；C11 控制完整 app_server 的任务释放顺序，但实际执行 Core 计算；C04 只让底层 enricher 抛错。C10 使用真实后台 worker 和浏览器，不使用 testServer 同步分支。C02 的 RNA 次级差异实际使用 DESeq2。

探针是诊断脚本，会打印条件和 `PROBE_ERROR`，并不是现成的 CI pass/fail 门禁。修复时应将相关条件转成有明确数值断言和失败退出码的回归测试。

`run_suite.R` 是修正 CSV 采集后的版本。原始主测试日志末尾保留了修正前的 CSV list 列错误：测试已经 DONE、测试结果 RDS 已写完，随后写 CSV 失败。现在的 CSV 与 totals 来自那些保存的 testthat 结果，未通过重复测试替换原记录。

`run_extra.R` 中两个真实数据位置是本机测试资源，其他机器需提供对应环境数据后才能复现。不要为了消除 skip 把私人研究数据提交到仓库。

L01 的证据保留了当时的未提交源码和 NAMESPACE；它是一个时间点的工作区状态。若 CC 随后修正导出声明，应在最终提交重新安装检查，而不能继续把旧的 L01 状态当成最终版本结论。
