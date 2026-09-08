# omicsApp 当前版本独立审查

审查日期：2026-09-08。结论：**可以继续受控部署测试，但当前证据不足以放行为正式结果交付版本。先修复本文 P1 问题，再冻结代码、构建镜像并验收。**

框架已经有扎实的可用基础：差异分析的既有测试、真实数据对照、较大矩阵性能测试均通过。本轮仍复现了取错数据层、结果签名失配、富集失败被当作空结果、导出不能复现等缺陷。这些问题不能靠“页面能打开、分析能跑完”排除。

本次只新增审查文档和证据，没有修改业务代码，没有部署或重启服务器。

## 1. 版本和验证范围

**完整测试基线：`a5ae7a310297985969f1f2c81498880d6c042644`。** 开始时通过 `git archive` 建立隔离快照，避免 CC 同时修改源码影响结果。

审查期间 HEAD 前进到 `bba74b5592333be52e580a8da1c7984335ebdb11`。已提交差异只增加三个部署辅助脚本、删除意外提交的输出文件，`packages/` 业务代码与上述基线一致。新增示例生成/核对脚本另行运行；计时脚本另行读码检查。

审查末尾工作区还有尚未提交的绘图抽稀改动。本文 **L01** 单独记录该改动的即时安装验证，不能将基线测试结果视为这批未提交代码已通过测试。交付修复时应以新的固定提交重新验收。完整版本、状态及证据摘要见 [metadata.json](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/metadata.json)。

验证环境为本机 macOS arm64、R 4.4.3，Shiny 1.13.0、future 1.69.0、promises 1.5.0、DESeq2 1.46.0、limma 3.62.2、edgeR 4.4.2。**未验证测试服务器的实际镜像、CPU/内存峰值、OIDC 登录或宿主机卷权限**；本机 Docker daemon 不可用。下文源码行号对应审查基线，后续修改可能使行号移动。

### 已完成的验证

| 验证 | 结果 | 解释 |
|---|---|---|
| omicsCore 主套件 | 2,677 项断言通过，0 failure / 0 error；11 skip、1 warning | 小型 DESeq2 数据的拟合 warning；大部分可选 skip 已在补测中覆盖 |
| omicsApp 主套件 | 1,658 项断言通过，0 failure / 0 error；5 skip、5 warning | 包含真实 Chrome 的 journey、smoke、two-sessions；warning 为 C locale 中文路径和局部 NA 系数 |
| 真实 RNA-seq 对照及性能补测 | 33 项通过，无失败、警告或跳过 | 真实 follicle RNA 数据 63,241 × 129，edgeR 对照效应量和统计列在既有 1e-8 容差内一致；另含 63,241 × 258 矩阵与 20,000 × 258 limma 性能预算 |
| 真实蛋白组补测 | 23 项通过，无失败、警告或跳过 | 原始工作簿、归一化和既有分析流程对照 |
| 重型 API 合同检查 | 43 项通过，无失败、警告或跳过 | 包含已在主套件运行的轻量检查，不能作为 43 项全新覆盖重复累加 |
| Core `R CMD check` | 0 error / 0 warning / 0 note | 使用 `--no-manual --no-tests --no-vignettes`，测试由上面的独立运行覆盖 |
| App `R CMD check` | 0 error / **1 warning / 1 note** | 非 ASCII R 字符串；报告模块闭包的 `bundle` 静态绑定提示，详见第 5 节 |
| 新增双组学示例 | 60 人、RNA 3,900 基因、蛋白 1,200 特征；60 对 donor 配对 | 进入两层显著集合的 154 个预设变化对，方向误判 0；并非全部 240 个预设变化对都达到显著 |
| 定向边界复现 | 见 C01–C12、L01 | 使用真实模块/Core 函数；部分探针控制任务返回次序或注入底层错误，逐项注明 |

主套件合计 **4,335 项通过断言**，真实数据/性能补测另 **56 项**。这不是 4,391 个独立使用场景，更不是目标 Docker 环境的压测结论。

审计采集脚本在两套主测试显示 `DONE` 并保存结果 RDS 后，写 CSV 时发生 list 列处理错误。随后已从保存的测试结果恢复 CSV 和统计；原日志保留。这是审计输出脚本的问题，不是产品测试失败，也没有将日志尾部错误删去。

## 2. 审查期间的工作区即时回归

### L01 · P1 · 新绘图改动使 `plot_volcano()` 失去导出

**范围：未提交工作区，单独验证；不属于 a5ae7a3 基线。**

新增 `thin_volcano_points()` 被插在原来属于 `plot_volcano()` 的 roxygen 文档与函数之间。生成后的 NAMESPACE 删除了 `export(plot_volcano)`，新增 `export(thin_volcano_points)`；App 仍然通过 [mod_diff_view.R:608](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_diff_view.R:608) 调用 `omicsCore::plot_volcano()`。

我把相关工作区改动复制到隔离包并安装，实际得到：

```text
VOLCANO_EXPORTED=FALSE
QUALIFIED_ACCESS='plot_volcano' is not an exported object from 'namespace:omicsCore'
```

因此，若按捕获到的状态构建镜像，火山图调用会失败。开发态测试能访问包内函数，不能替代安装后的导出接口检查。

**修复及验收：** 为两个函数分别明确文档和导出归属；辅助函数若不作为公共 API，应设为内部函数。重新生成 NAMESPACE/Rd 后安装包，直接调用 `omicsCore::plot_volcano()`，并在安装版 App 验证火山图。证据：[live-exports.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/live-exports.log)。

## 3. 应在正式结果交付前修复的问题

### C01 · P1 · 同一组学的多个数据层可能被混淆，整合可以变成自身比较

**位置：** [mod_integration_view.R:72](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:72)、[export-script.R:232](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-script.R:232)。

这两条路径通过 `omics_type` 找第一个匹配层，忽略结果已有的 `layer_tag`。两个蛋白组层 `skin`、`follicle` 同时存在，主差异结果来自 follicle 时，整合把主层名认成 skin，第二层选成 follicle。主结果实际上已经是 follicle，第二层再次计算 follicle，最终可以输出一张看似正常的“skin 与 follicle”整合表。

真实模块探针得到 30 行结果、两侧效应量最大差异为 0；同类项目导出的脚本把 follicle 的分析写成 `run_diff(input_skin, ...)`，构造数据上重放效应量最大差异为 8。

**适用边界：** 普通 UI 导入通常每种组学只有一个槽，此问题集中于从 Core/其他流程建立并恢复的同类型多层项目。不能据此说所有 RNA–蛋白整合都错误，但框架既然允许这些项目，就不能静默取第一个层。

**修复及验收：** 统一依据明确的层标签和来源绑定定位；旧项目只有类型信息且匹配不唯一时明确拒绝或要求选择。测试“同类型双层、主结果来自第二层、交换层顺序、保存恢复、脚本重放”，核对真实层标签和数值。

证据：[review-probes-5.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-5.log)、[review-probes-2.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-2.log)。

### C02 · P1 · 带协变量的整合结果在签名后又被修改，随即被判为无效

**位置：** [mod_integration_view.R:382](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:382)。

回调先 `stamp_from_launch(result, launch)`，再写入 `secondary_covariates`、`secondary_paired_col`、`unadjusted_secondary`。签名对应修改前的参数，项目校验对应修改后的参数，二者冲突。

使用蛋白组 limma + age 协变量、第二层 raw-count RNA + DESeq2 的真实计算路径复现：整合计算成功，但 `params_signature` 不等于最终参数摘要，注册校验提示 `params signature does not match the saved parameters`，`retire_invalid()` 后整合被移除。未带这些额外参数的默认路径不一定触发。

**修复及验收：** 先形成完整最终参数，再生成签名及 analysis ID；签名后不再修改参与签名的字段。分别覆盖协变量、配对列、允许第二层不调整三条路径，验证计算→注册→失效扫描→保存→恢复全链条仍保留有效结果。

证据：[review-probes-4.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-4.log)。

### C03 · P1 · 批量富集未继承页面阈值，并丢失实际运行参数

**位置：** [mod_report_view.R:137](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_report_view.R:137)、[run-enrichment-grid.R:136](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/run-enrichment-grid.R:136)。

单次富集已经接入差异页的 p 值口径、阈值和效应量阈值；批量按钮只传入 diff bundle 与进度函数，使用默认 adjusted p < 0.05、无效应量过滤。它也未像单次路径一样设置 `output_p_cutoff = 1`，ORA/GSEA 的存储输出仍可能在 0.05 截断。

因此，用户在差异页调整筛选后，“一次跑完”并不等于按当前设置把各库逐一跑完。GSEA 应继续使用完整排序，问题在输出截断和记录口径，不应通过将差异基因先筛一遍来修复 GSEA。

即使直接调用 Core grid 并明确传入 `p_cutoff=0.001, effect_cutoff=2, output_p_cutoff=1`，子调用确实收到 0.001，最终 bundle 却没有记录这些参数，也没有合并各子任务的基因集来源信息。

**修复及验收：** App 形成一份明确的批量请求参数快照；Core 保存解析后的全部统计参数和每库来源/版本。对同一 diff、同一阈值，逐库比较单次与批量的 ORA 输入基因及结果；验证 GSEA 排序不受 ORA 阈值影响，且保存输出满足既定完整性约定。

### C04 · P1 · 富集后端失败会被批量任务显示成正常“空结果”

**位置：** [enrich-ora.R:23](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/enrich-ora.R:23)、[enrich-gsea.R:31](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/enrich-gsea.R:31)、[run-enrichment-grid.R:95](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/run-enrichment-grid.R:95)。

内部后端先把 error 转成 warning + NULL，外层 grid 又 `suppressWarnings()`。于是外层的错误捕获见不到 error，最后记为 `status="empty"`、`message=NA`，结果 warnings 也为空。

在真实 Core 调用链中，仅将底层 `clusterProfiler::enricher()` 替换为抛错函数，复现了上述结果。这验证的是错误处理机制，**不是宣称真实数据库此次发生了相同故障**。现有测试若只 mock 外层 `run_enrichment()` 抛错，会漏掉它。

**修复及验收：** 明确区分成功、有效空集合、后端失败；失败可以不终止其他库，但必须保存库名、阶段和错误消息。分别注入 ORA/GSEA 底层异常，要求日志为 failed、页面与 Excel 可见；正常无富集结果仍应为 empty。

C03/C04 证据：[review-probes.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes.log)。

### C05 · P1 · 整合任务传输仍携带分析登记表、QC 矩阵和历史结果

**位置：** [async_helpers.R:301](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/async_helpers.R:301)、[mod_integration_view.R:371](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:371)。

先前的差异模型裁剪及 `bundle_for_worker()` 已改进单次富集的大对象问题。但 `project_for_worker(..., assays=FALSE)` 只删除 `$bundles` 和 experiments 下矩阵，保留了 `$analyses`、`$preprocessing`、`$retired`。这些位置仍可以保存旧项目的 DESeq2 模型、实际 QC 矩阵和历史负载。

探针在真实项目登记表放入 4 MiB 模型并应用 QC：处理后的 experiment 矩阵为 0×0，但登记表仍有全部 4 MiB 模型，preprocessing 仍有 30×8 矩阵，退休条目仍存在。

这是对遗漏路径的有界复现，**没有把用户的 1.85 GiB 错误在当前版本完整重现，也不能据此断定当前每个项目都会超限**。但大项目/旧项目可以重新暴露同类 future globals 限制和序列化开销。

**修复及验收：** concordance 任务按实际需要白名单构造轻量项目描述，不再复制整个 project 后逐项删除。测试必须包含应用过 QC、有 analyses 和 retired、恢复过保留模型的旧项目；增加任务请求字节数检查，同时比较裁剪前后结果一致性。不要把单纯提高 `future.globals.maxSize` 当作完成修复。

证据：[review-probes.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes.log)。

### C06 · P1 · 手动修正输入尺度后，复现脚本仍使用旧尺度

**位置：** [mod_import_view.R:185](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:185)、[export-script.R:117](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-script.R:117)。

尺度控件修改实际输入的 `assay_type` 和 fingerprint，但未同步记录到导入 provenance。脚本导出则重放旧 `read_omics()` 参数。

使用取整的对数尺度数据：自动识别 raw_count，用户手动改成 logcpm 后 App 的 limma 可以运行；生成的脚本仍读成 raw_count，随后被模型入口的尺度检查拒绝。补充验证中保留了实际归档的原文件，并在全新 R 进程直接执行导出的脚本，仍得到同一尺度错误。当前分析结果不因此变错，但下载的“复现脚本”不能复现。

**修复及验收：** 记录最终确认的尺度决策，脚本重放最终设置，并与文件角色、归一化处理一起保持一致。覆盖自动识别→人工改尺度→分析→保存恢复→全新 R 进程执行脚本；比较重放的有效输入摘要及差异结果。

证据：[review-probes-2.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-2.log)、[replay-process.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/replay-process.log)。

## 4. 状态、性能和导出问题

### C07 · P2 · 新建项目后，原文件再次确认仍不能导入

**位置：** [mod_project_view.R:372](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_project_view.R:372)、[mod_import_view.R:618](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:618)。

新建清空 current_project，但导入模块仍保留 parsed/confirmed_input。再次确认相同对象时 reactiveVal 没变化，主服务收不到新提交事件。

真实 `app_server` 探针得到：首次导入有项目→新建后为空→再次确认仍为空，导入模块却仍为 confirmed。

**修复及验收：** 将新建/激活显式传入导入模块以重置暂存状态，或采用带提交编号的导入事件。用实际浏览器覆盖“导入文件→新建→重新选择或复用相同文件→确认”。证据：[review-probes-3.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-3.log)。

### C08 · P2 · 新批量富集模块没有接入已有的过期回调合同

**位置：** [mod_report_view.R:100](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_report_view.R:100)。

该模块使用 `activation=NA_integer_`，调用 run_async 时没有 admissible/request_seq/on_discard。项目为空、无已存 grid 或 grid 失效时，恢复 observer 只 return，没有清理 grid_bundle/grid_error。

人为控制任务返回次序的真实模块探针显示：A 任务的晚返回仍进入 B 的模块状态，晚错误可覆盖错误提示；新建项目后旧 grid 状态仍保留。**不能据此宣称 A 的结果已经写入 B 的持久项目**：app_server 的 project_id 与 input-current 检查仍拒绝该结果，这层防护确实有效。另外，当前普通 promise 路径会阻塞发起会话的反应式处理，控制回调的测试不等于已经在生产浏览器实测切换中的并发竞争。

**修复及验收：** report 模块接入实际 activation、每视图请求序号与完整来源检查；成功/失败回调同受限制。切项目、新建、上游失效都主动清空旧状态。修复 C10 的真正非阻塞任务时，必须同时保留此合同。

### C09 · P2 · 已导出的批量富集结果没有对应复现代码

**位置：** [export-script.R:358](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-script.R:358)。

脚本导出只处理 qc/diff/gsva/enrich/integration，没有 enrich_all 分支。包含批量结果的项目生成脚本后，既没有 `run_enrichment_grid()`，也没有遗漏提示。Excel 支持 grid，脚本却不支持，是新结果类型未接完所有出口。

**修复及验收：** 在 C03 完整参数记录基础上增加 grid 重放。最好从统一分析类型定义生成存储/导出分发，不再让各出口维护各自的类型清单。带 grid 的 `.omp` 在新进程执行导出脚本，应得到相同任务组合、参数和结果；无法重放的项明确说明。

C08/C09 证据：[review-probes.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes.log)。

### C10 · P2 · 当前真实异步运行中，新增进度轮询未能更新进度条

**位置：** [async_helpers.R:91](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/async_helpers.R:91)、[async_helpers.R:123](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/async_helpers.R:123)。

进度轮询是普通 Shiny observer，而调用链把 promise 返回到发起任务的 observer。在真实 Chrome + Shiny + multisession(2) 最小应用中，使用项目实际的 run_async/progress_reader/detached_call，worker 每秒写一步、持续六秒：浏览器六次采样都为 30%，服务器轮询仅发生一次；任务结束才显示 done。此次没有走 testServer 的同步捷径。

普通反应式环境里的 promise 可以让其他会话继续服务，但仍会阻塞发起会话的常规处理；这一机制与本次现象一致。[Shiny ExtendedTask 官方说明](https://shiny.posit.co/r/reference/shiny/1.8.1/extendedtask.html)。

**修复及验收：** 评估将长期计算统一迁移到 ExtendedTask + future_promise，确保启动事件及时返回，进度轮询能独立推进，并保留任务准入、关闭会话清理和有界并发。使用安装版 App、实际 future plan 做浏览器测试：任务未结束时至少收到两次中间进度，取消/关闭后没有遗留轮询；在 Docker 再验证 multicore 路径。[官方实现模式](https://shiny.posit.co/r/articles/improve/nonblocking)。

证据：[progress-browser-2.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/progress-browser-2.log)。

### C11 · P2 · 对整个项目的依赖导致无关更新触发 QC 和重复整合

**位置：** [mod_integration_view.R:66](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:66)、[mod_integration_view.R:419](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:419)、[mod_qc_view.R:350](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_qc_view.R:350)。

can_run 依赖整个 current_project，observer 在其重新失效且 ok 时自动计算。整合自己的结果写回项目也会触发这条链。

对完整 app_server 控制调度、实际执行计算并释放结果：一次主差异结果触发 **两次整合任务**，第二次之后稳定；没有复现无限循环。独立 QC 探针中，仅修改无关项目 metadata 就额外调用一次 run_qc。

**修复及验收：** 昂贵计算的自动触发应依赖明确任务键：层/有效输入版本、模型参数、上游 analysis ID；不依赖整个项目对象。相同任务键复用结果，手动重跑仍可显式发起。增加调用计数断言：整合自身保存不重算、改项目名/写富集结果不重跑 QC、真正改变输入才重新计算。

证据：[review-probes-3.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-3.log)、[review-probes-4.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes-4.log)。

### C12 · P2 · Excel 仍以单个显示槽导出，多层项目的分析说明不完整

**位置：** [export-workbook.R:69](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-workbook.R:69)、[export-workbook.R:115](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-workbook.R:115)。

工作簿只读 `bundles$qc` 和 `bundles$diff`，未遍历按层 preprocessing，也未包括整合所依赖且已经登记的第二层差异结果。参数页的 features/samples 取原始 experiments，未区分基线与实际分析输入数量。

双层均应用 QC 的项目，探针只得到一组 QC 表及一张 Differential 表。该导出可以作为“当前显示槽表格”，但不足以作为“完整多层分析交付物”；当前说明没有明确列出遗漏范围。

**修复及验收：** 明确报告范围，至少包含当前结果的所有有效上游、逐层 QC 参数、实际输入数量，并显示 layer_tag/analysis_id。无需默认导出所有退休历史，但省略哪些数据必须明确。双层 QC + 整合项目导出后，应能逐层追溯两侧统计结果及预处理。

证据：[review-probes.log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/review-probes.log)。

## 5. 为什么已有大量测试仍会漏掉这些问题

问题集中在测试边界，而不只是测试数量不够。

1. **开发态与安装态不同。** L01 中包内调用可能通过，`omicsCore::plot_volcano()` 的安装接口却不可用。App `R CMD check` 还出现新的非 ASCII warning 和闭包 `bundle` note；note 本身不是已证实的运行时异常，但应通过明确闭包接口消除，不宜用全局变量声明掩盖未定义依赖。证据见 [App check log](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/omicsApp-00check.log)。
2. **同步测试无法证明异步会话正确。** run_async 在 testServer 中会走同步分支；现有 Chrome 测试直接启动 `inst/app/app.R`，没有像 `launch()` 一样设置生产 future plan。它们仍有价值，但不能覆盖 C08/C10 的实际调度机制。
3. **模拟错误的位置过高。** mock run_enrichment 抛错能验证 grid 外层捕获，却漏掉 C04 的底层 error→warning→被吞掉链。
4. **输入夹具过简单。** 只包含 bundles 和 experiments 的轻量项目不能发现 analyses/preprocessing/retired 的负载；每种类型只有一层不能发现 C01；默认无协变量不能发现 C02。
5. **只验证对象形状，没有验证完整使用路径。** 返回 analysis_bundle、生成 XLSX/R 文件，都不等于分析来源正确、脚本能在新进程执行、工作簿包含所有约定内容。
6. **一个多库测试错误地提前 skip。** [test-enrich-edge.R:67](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/tests/testthat/test-enrich-edge.R:67) 调用 `msigdbr::msigdbr("Homo sapiens", db)`。本机 msigdbr 的第一参数是 db_species，实际立即报 `nchar(db_species) not equal to 2`，却被解释为“No kegg gene sets available”。这不是证实数据库不可用。应通过 Core 的正确数据库映射检查资源，并分别运行各库测试。

新增 `check_demo_integration.R` 的数据方向核对很有帮助；本次确实得到 0/154 方向错误。但它目前只是打印错配计数，没有在方向错误、覆盖过低、配对方式错误时返回失败退出码，不能直接作为部署门禁。

新增 `time_pipeline.R` 也应修订后再用于诊断用户的 DESeq2 慢：它从单一 bundles$diff 取参数而未核对指定层，重跑时未传协变量/配对参数；单独拟合第二次 DESeq2 后用两次总时间相减，不能精确隔离第一次的 extract 耗时；serialize 时间是序列化指标，不是实测 worker 完整回传时间。应在同一次真实调用内分段计时。

## 6. 建议的修复顺序与部署验收出口

**第一批：固定结果正确性和可复现性。** 修 L01、C01、C02、C03、C04、C06。每项都以完整触发链回归测试验收；尤其检查安装接口、真实层标签、最终签名、底层失败状态和全新进程脚本重放。C09 与 C03 的参数合同一起完成。

**第二批：完成任务及状态合同。** 修 C05、C07、C08、C10、C11。重点是轻量任务输入、真实 activation/request_seq、真正异步浏览器测试和避免无关项目写回触发昂贵计算。不要先通过增加 worker 或提高 globals 上限掩盖任务负载问题。

**第三批：统一交付出口。** 修 C12，清理 check warning/note、误 skip 和示例验收脚本。将分析类型、参数、来源绑定、导出能力集中定义，减少“新模块能跑，但忘接清理/恢复/脚本/Excel”的重复缺陷。不需要推倒现有架构，也不需要改写成熟的 DESeq2/edgeR/limma 统计后端。

**随后在冻结提交的 Docker 镜像上完成：**

| 场景 | 必须记录或核对的内容 |
|---|---|
| 用户约 200 例 RNA-seq 原场景 | DESeq2 输入维度、设计公式、拟合/整理结果/任务传输/浏览器绘图/自动保存分段耗时；实际 CPU 限额、内存峰值及镜像版本 |
| 两层真实数据端到端 | 两层应用 QC、主层分析、整合、单次/批量富集、Excel 和脚本；来源及参数一致 |
| 状态切换 | 新建并重用文件、A→B→A、恢复含非默认 QC 参数的旧项目、任务晚成功和晚失败 |
| 真正异步与并发 | 任务进行中收到中间进度；新请求正确淘汰旧回调；队列有边界；会话结束可清理；多用户资源总量不超预算 |
| 持久化故障 | 自动保存失败提示、磁盘满/配额边界、容器中断后恢复最近一次已成功保存状态、备份实际还原 |
| 部署配置 | OIDC 登录、不同用户目录隔离、容器实际 UID 的读写权限、TLS/WebSocket、部署镜像 digest 与冻结提交一致 |

仓库部署模板已有非 root 运行、OIDC sub 作为稳定用户目录、网络分离、线程限制、资源上限、依赖/镜像验证等基础。这些配置值得保留；源码模板存在不等于测试服务器已经正确应用。模板的 2 CPU / 6 GiB 也不能替代 200 例样本的容量测量。

自动保存当前每次项目更新就同步写入。它保住了较短恢复窗口，不能未经测量就改成长 debounce；先消除重复计算和冗余模型、测量实际保存耗时，再决定是否需要带持久化保证的优化。

## 7. 对当前框架的总体评价

**核心统计计算的已有证据较好；App 的跨模块合同尚未全部落实。** 新增批量富集和导出代码也位于 omicsCore 包内，因此本轮修复不能全部放在 App。这里应区分“成熟统计后端”与“Core 中新增加的编排、错误处理、导出封装”，不能把后者的问题归为 DESeq2 本身不稳健。

之前修复的数据版本、来源绑定、会话 token、按层 preprocessing 等设计是有价值的。C08 中持久写入检查拒绝错误项目的结果，正说明这些保护已经在工作。当前主要问题是：一些新模块没有统一采用这些合同，或在签名后又改变对象，或仍从旧的显示槽/类型推断信息。

下一步最有价值的投入是把这些合同变成统一入口及跨模块验收，而不是继续增加功能数量。修复本文问题并完成目标镜像上的验收后，才有充分依据把它从“功能基本成型”提升为“可以稳定交付研究结果”。

证据目录：[current-version-a5ae7a3](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-08/evidence/current-version-a5ae7a3/README.md)。所有探针、关键原始日志、结果 CSV、版本信息和 SHA-256 清单均保留；未复制真实研究原始矩阵、生成的示例工作簿或安装包库。
