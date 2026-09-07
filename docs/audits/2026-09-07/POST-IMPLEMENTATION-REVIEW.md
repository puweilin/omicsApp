# CC 实施后复审：部署测试准入

日期：2026-09-07。审查提交：`2b64bed0e7d849c13e20c1307c594206698c0295`。
对比基线：`f756751c24c12b01cf99fca399a6a5a17eeea4f4`，共 11 个提交。
审查开始及验证结束时业务工作区均干净。本轮只增加审查文档与证据，没有修改业务源码、部署配置或原整改方案。

## 1. 结论

**不能确认“完整整改已完成”，暂不通过批次 2–5 完成后进入部署验收的门槛。**

可以在隔离环境继续做镜像构建、启动、登录、反向代理与挂载的基础演练；这些演练不能替代本报告列出的业务正确性修复，也不应使用真实研究结果作为可交付产物。先修复 N01–N06 及 N10 涉及的正确性/合同缺口，再进行完整业务链的部署测试。N07–N09 和后文列出的保存状态、CI 等要求也应纳入本轮交付。

此次不是推翻前几批修复。安全模型设计、RDS 入口、有效 QC 输入、完整富集输出等已经有实际改进。主要问题仍是：**底层合同已经实现，但部分页面没有调用；局部测试通过，没有覆盖项目激活、跨模块写回与脚本重放的完整路径。**

特别是上一轮 C1 和 C3，当前实现仍没有在真实 App 路径闭合：完整 `callback_admissible()` 没有接入页面，整合内部第二层差异分析没有进入来源登记与导出流程。

## 2. 验证范围与证据边界

从被审提交 `git archive` 建立独立副本 `/tmp/omicsapp-review-20260907.SsorGW`，在副本中运行，未使用本机旧安装包代替新源码。

| 验证 | 结果 |
|---|---|
| omicsCore 主测试套件，`NOT_CRAN=true`，包含本机 RNA-seq golden 数据 | testthat 完成，无测试断言失败；1 条警告，6 项跳过 |
| omicsApp 主测试套件，`NOT_CRAN=true`，指定本机蛋白组工作簿 | testthat 完成，无测试断言失败；5 条警告，无跳过 |
| App 的真实浏览器测试 | 现有 app journey、七页面 smoke、双会话测试实际运行，未报告失败；它们不等于覆盖本报告新增的状态场景 |
| 另行开启 Core heavy API sweep 和性能预算 | 52 项断言通过，0 失败、0 错误、0 警告、0 跳过 |
| 独立 library 安装 omicsCore | 成功 |
| 两个包的 R CMD check | 均 0 error、0 warning、0 note；使用 `--no-manual --no-tests --no-vignettes`，测试已另行运行 |
| 新增针对真实模块的复现 | 见 N01–N08，探针脚本最终执行退出码为 0，无 `PROBE_ERROR` |
| 大体积正常备份输入 | 复现误拒绝，见 N09；没有访问真实 Keycloak/PostgreSQL |
| Docker 镜像、目标服务器 | 未验证：本机 Docker daemon 未启动，socket 不存在；没有目标服务器放行证据 |

主测试日志末尾存在**审查记录脚本**的 CSV 写出错误：`as.data.frame(testthat_results)` 包含 list 类型的 `result` 列，我在 testthat 输出 DONE 后直接 `write.csv()` 导致记录阶段退出码为 1。它不是包测试断言失败；原日志完整保留，也没有将这两个进程伪报为退出码 0。本报告不声称没有成功取得的主套件断言总数。额外测试的汇总脚本去除了 list 列，52 项结果有完整 CSV。

Core 主套件的 6 个跳过中，heavy sweep 与 4 项性能测试已由额外运行补齐；仍有一项 KEGG 多数据库测试因当地缓存未提供相应基因集而跳过。警告主要来自小样本 DESeq2、C locale 中文路径测试和一个含 NA 系数的 QC/整合用例。

证据：

- [复现脚本](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/review_probes.R)、[复现输出](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/review-probes.log)。
- [Core 主套件日志](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/core-test.log)、[App 主套件日志](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/app-test.log)。
- [额外测试](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/optional-results.csv)、[包检查](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/package-check.log)、[审查元数据](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/summary.json)。

回调复现使用真实模块构造的工作函数、接收检查和回调，只用可控调度替身延后返回，便于确定性地改变状态。**这证明现有回调合同存在缺陷；不声称已经在当前仍阻塞发起会话的真实浏览器中完成所有“计算中切项目”的操作。**

## 3. 必须处理的发现

### N01｜P1：完整回调检查没有接入，旧结果还能被重新标成新项目的结果

位置：[mod_diff_view.R:301](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_diff_view.R:301)、[mod_diff_view.R:334](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_diff_view.R:334)、[async_helpers.R:158](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/async_helpers.R:158)。

Core 的 [callback_admissible](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/analysis-binding.R:499) 已实现来源与令牌的组合检查，但 App 三个分析模块全部调用 `newest_request_only()`。`activation` 没有传入这三个模块，项目/数据失效也没有使其请求序号失效。该简化检查甚至在无法读取当前请求序号时默认放行。

更严重的是，差异分析成功后用 `current_project()` 调用 `bind_result()`：输入版本来自启动时，项目身份却来自结束时。如果 A、B 使用同一份数据，晚到的 A 结果会得到 B 的身份，并通过写回的输入版本检查。

真实 diff 模块的可控晚到回调复现：

```text
VIEW_EMPTY_AFTER_SWITCH=TRUE
LATE_SUCCESS_ADMITTED=TRUE
A_RESULT_STAMPED_AS_B=TRUE
WRONGLY_STAMPED_RESULT_VALID_IN_B=TRUE
LATE_ERROR_ADMITTED=TRUE
CURRENT_ERROR=late error from A
```

整合回调还在完成时重新读取两层 `effective_input(current_project(), ...)` 的版本，能够把旧计算描述成新输入版本。这需要与差异、富集一起修。

最小修复：发起时捕获真实来源、项目身份与完整令牌；所有成功/失败回调统一调用完整门禁；检查与提交连续完成；不在回调时从新项目补写来源。会话结束与无法验证的上下文都必须拒绝。验收必须驱动真实模块创建的门禁，而非仅给 `run_async()` 注入一个永远返回 FALSE 的测试函数。

### N02｜P1：QC 恢复仍丢失未暴露参数，并会在普通项目更新时覆盖草稿

位置：[mod_qc_view.R:223](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_qc_view.R:223)、[mod_qc_view.R:262](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_qc_view.R:262)、[mod_qc_view.R:559](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_qc_view.R:559)。

预览仍只传三个可见参数，没有合并已保存的 `outlier_sd_threshold`、`sample_missing_threshold` 或填补方法额外参数。恢复项目后点击 Apply，会以默认值覆盖它们：

```text
SAVED_OUTLIER_SD=7
REAPPLIED_OUTLIER_SD=3
SAVED_SAMPLE_MISSING=0.9
REAPPLIED_SAMPLE_MISSING=<NULL>
```

这直接违反已确定的 D6。已有清洗矩阵可以正确恢复，并不能证明下一次应用时仍保留原参数。

此外，恢复 observer 监听 `list(current_project(), active()$tag)`，并非仅监听激活/层变化。普通项目对象更新也发恢复消息。实测当前草稿为 `none` 时，仅改变项目名称，就向浏览器发送把方法恢复成 `min` 的更新；分析结果写回也会更新项目对象。

最小修复：将已保存的完整参数作为基底，可见控件只覆盖对应字段；以激活批次和层身份触发一次恢复。无已应用状态的新层使用明确默认值；当前环境缺少已保存方法时要明确提示，不能静默换选项。补“恢复→编辑→无关分析完成→应用”的浏览器用例。

### N03｜P2：QC 刚应用就清掉自己的输出，项目没有保存对应 QC 结果包

位置：[app_server.R:130](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/app_server.R:130)、[mod_qc_view.R:136](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_qc_view.R:136)。

应用 observer 写入 `preprocessing` 后增加 `layer_generation`。QC 自己监听这个信号并清空 `applied()`；优先级较低的 bundle 写回 observer 随后只能读到空 QC 输出。

真实 `app_server` 复现：

```text
PREPROCESSING_PRESENT=TRUE
QC_BUNDLE_PRESENT=FALSE
QC_APPLIED_MODULE_PRESENT=FALSE
```

因此，清洗数据确实已经到达下游，这是 A03 修复的有效部分；但 QC 结果包、诊断图和该结果的来源登记没有一起保存。报告仍能从预处理 provenance 列出步骤，**不能将此缺陷夸大成完全丢失 QC 参数或清洗矩阵**，但“完整保存并恢复 QC 分析”尚不成立。

最小修复：QC 应用在同一次项目更新中登记清洗层、QC bundle 与来源，只使依赖旧有效输入的下游分析失效。不要用清空自身输出的全局信号完成提交。

### N04｜P1：整合的第二层来源仍未登记，导出脚本不能执行完整整合

位置：[mod_integration_view.R:274](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:274)、[mod_integration_view.R:306](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:306)、[export-script.R:393](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/export-script.R:393)。

`sec_diff` 仍在 worker 的局部作用域中，完成后只返回 integration。整合的来源绑定只记录主 diff；第二层没有被登记到 `project$analyses`。新增了 `analysis_id_b` 列，但普通 App 计算路径中的值为 NA。

本轮使用真实双层工作簿、真实 diff/integration 模块、真实脚本生成与执行得到：

```text
SECONDARY_ANALYSIS_IDS=NA
INTEGRATION_UPSTREAM_COUNT=1
INTEGRATION_REPLAY_ERROR=`diff_bundles` is required for method = 'concordance'.
```

[实际生成的脚本](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/integration-replay-generated.R) 只计算一份 diff，再调用没有 `diff_bundles` 的 `run_integration()`。附加 NOTE 不能使它成为可执行复现脚本。

单层 log2/vsn/QC 的 provenance 重放已有实质修复；双层重放和上一轮 C3 仍未落实。脚本也尚未按来源图选择各下游实际依赖的历史结果，仍按固定 `bundles$diff` / `bundles$enrich` / `bundles$integration` 生成。

最小修复：第二层 diff 与 integration 作为一个提交单元返回并登记；整合绑定两个真实上游 ID；导出按依赖关系重建两层模型及其输入。验收必须在独立目录从归档原文件执行整份脚本，并对齐键后比较两层统计量和整合结果。

### N05｜P1：limma 尺度入口仍能放行明确标为原始强度的小数值矩阵

位置：[diff-utils.R:151](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/diff-utils.R:151)。

把入口校验放在 `auto` 解析之后是正确的，也能阻止大数值原始强度。但现在 `max(expr_mat)` 不超过“合理对数值”上限就放行，包括明确的 `raw_intensity`。数值小并不能证明已经做过对数转换。

构造明确为线性尺度、最大值约 17 的正值原始强度矩阵，通过真实 `run_diff(method="auto")` 得到：

```text
RAW_INTENSITY_ACCEPTED=TRUE
EFFECT_TYPE=log2FC
MAX_EFFECT_VS_LOG2_DIFFERENCE=3.65358853183165
```

构造器会给尺度警告，但随后仍返回标为 log2FC 的线性均值差。对显式 raw 标签的校验应遵守尺度合同，不能用量级猜测覆盖它。旧项目的模糊 `intensity` 标签兼容，应作为单独迁移/确认路径处理；保留正确的 lm raw-count 转换不要求保留此绕过。

### N06｜P1：改工作表角色重导入仍被旧指纹去重提前拦截

位置：[mod_import_view.R:555](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:555)、[mod_import_view.R:627](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:627)、[mod_import_view.R:662](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:662)。

`read_omics()` 已保存 `sheet_roles`，但 App 的 `input_fingerprint()` 仍没有这个参数，调用点也没有传入。新的内容版本机制到不了这个提前返回的路径。

在真实导入模块中，用同一工作簿的两个表达矩阵，第一次选 A，随后把角色改为选择 B。候选矩阵每个值变化 5，点击 Confirm 后仍认为是同一个文件：

```text
ROLE_CHANGE_CANDIDATE_MATRIX_DIFF=5.0000000000001
ROLE_CHANGE_NEW_INPUT_COMMITTED=FALSE
ROLE_CHANGE_PROJECT_MATRIX_DIFF=0
ROLE_CHANGE_FINGERPRINT_STILL_MATCHES=TRUE
```

最小修复：指纹纳入完整解析配置/角色的稳定摘要，并统一预览、确认、去重的计算入口。验收要点击真实 Confirm，不能只比较两个 `input_version()` 的值。

### N07｜P2：打开项目保住了磁盘结果，但没有恢复分析页面

位置：[mod_project_view.R:329](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_project_view.R:329)、[mod_diff_view.R:253](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_diff_view.R:253)。

打开/恢复会调用激活处理并清空模块，避免旧视图覆盖存档，这部分修复有效。但 diff/enrich/integration 没有把验证后的项目结果装回各自 reactiveVal；diff 的方法、对比与协变量也没有从保存结果完整恢复。

在真实 `app_server` 打开带有合法来源和登记结果的 `.omp`：

```text
SAVED_PROJECT_DIFF_PRESENT=TRUE
RESTORED_DIFF_VIEW_PRESENT=FALSE
```

这里测试禁止自动计算，目的是验证“恢复已有结果”，不能靠重新算一次充当恢复。最小修复是打开时验证来源与依赖后恢复页面与参数，避免先显示默认设置、再误导用户重跑。

### N08｜P2：示例隔离只修了成功分支，真实项目空状态/失败时仍显示示例

位置：[mod_enrich_view.R:183](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_enrich_view.R:183)、[mod_integration_view.R:515](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:515)。

ActivePathways 新增的空状态依据是模块 `is_demo()`，但这个变量在真实项目未计算、缺上游或失败时仍为 TRUE。富集和整合其他图表也延续同一回退路径。

真实非空项目、没有上游分析时实测：

```text
REAL_ENRICH_IS_DEMO=TRUE
REAL_ENRICH_DEMO_ROWS=15
REAL_INTEGRATION_IS_DEMO=TRUE
REAL_INTEGRATION_AP_ROWS=6
```

最小修复：演示模式由明确项目模式决定，与运行状态分开。真实项目的未计算、失败、失效均显示对应空状态。不能仅断言“`is_demo(FALSE)` 时 AP 返回 NULL”。

### N09｜P2：新的备份脚本会把正常大转储误判为失败

位置：[backup-keycloak.sh:77](/Users/puweilin/SciProject/CHISSS/omicsApp/deploy/scripts/backup-keycloak.sh:77)。

临时文件、权限、pipefail、失败保留旧轮换文件都是正确改进。但 `gunzip -c ... | grep -q ...` 在找到文件开头的标识后提前关闭管道；大文件的 gunzip 随后遇到 SIGPIPE，在 pipefail 下整个管道非零，进入“不是 pg_dump 输出”的失败分支。

使用只输出标准头和约 10 MB 正文的本地 docker 替身，不访问数据库，复现：

```text
LARGE_VALID_DUMP_EXIT=1
PREVIOUS_BACKUP_PRESERVED=True
FAILED: the dump does not look like pg_dump output
```

最小修复：使用会完整消费输入的检查，或先解压到受控临时文件后校验。补充“大正常转储成功更新”“上游失败不覆盖旧摘要”“压缩失败不覆盖旧摘要”以及真实恢复演练。已有测试的小输入不足以暴露管道提前结束。[复现脚本](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/post-implementation-2b64bed/backup_probe.py)可从仓库根目录运行。

### N10｜P1：D3 未确认前暂停合并显著性结论的要求没有落实

位置：[integration-concordance.R:77](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/integration-concordance.R:77)、[plot-integration.R:188](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/plot-integration.R:188)。

方案 T4.1/D3 明确要求：统计口径确认前，只给描述性一致性和分层显著性。当前仍对按 `p_preference` 选出的两层 P 值计算 Fisher 合并值、再次校正，并在 dual volcano 中展示 `-log10(combined p)`；`is_significant` 也没有按方案拆分贯通消费者。

这是源码与已约定交付口径不一致的确认，**本轮没有另行证明某一具体数据集的合并 P 值错误**。应执行已约定的范围收缩，或者先完成 D3 的明确统计确认与对应实现/验证，不能把字段补齐视为该要求已经完成。

## 4. 原 A01–A17 的当前状态

| 原项 | 本轮判断 |
|---|---|
| A01 公式/对比代码执行 | 已改为生成的内部名称与数值对比；危险名字、中文、冲突及后端回归测试运行通过 |
| A02 网页 RDS | 网页传 `trusted=FALSE`，服务端固定拒绝；可信 Core 本地用法分开，入口测试通过 |
| A03 QC 到下游 | 有效清洗输入已接入 diff 与整合第二层；但 N02/N03 表明恢复与提交尚未完整 |
| A04 归一化/尺度 | 归一化失败停止提交、auto 后校验已修；显式小值 raw 绕过见 N05 |
| A05 脚本重放 | 单层 log2/vsn/QC 重放通过；双层依赖仍失败，见 N04；归档摘要验证、sample_link、基因集版本等仍需完整验收 |
| A06 归属/项目切换 | 激活、来源结构、登记和递归失效已有实现；真实门禁、视图恢复及导入角色路径见 N01/N06/N07 |
| A07 两层散点 | 已保存真实 effect_a/effect_b 并移除反推；整合来源和统计口径仍见 N04/N10 |
| A08 模型设置 | 拒绝与裁剪并存；第二层传入协变量/配对，缺列有提示与用户选择；相应用例运行通过 |
| A09 富集截断 | ORA/GSEA 两条 UI 调用均传 `output_p_cutoff=1`；参数保存与阈值用例运行通过 |
| A10 autosave | 失败可见及启动写权限检查已加入；未保存/保存中/最后成功时间等完整状态合同还未完成 |
| A11 项目目录 | 读取/删除解析为存储目录列出的真实项目；越界与中文目录测试通过 |
| A12 备份 | 失败保留旧文件已改；正常大转储失败，见 N09；真实恢复尚未演练 |
| A13 中文 PDF | 已加 XeLaTeX、CJK 字体配置、工具链检查和构建自检；本轮没有完成正式镜像内中文字形/图形验收 |
| A14 launch(project) | 显式启动项目已接入；启动加载相关用例运行通过，恢复视图仍受 N07 影响 |
| A15 示例/报告 | 报告新增完整结果、参数与主要图形；示例失败路径仍有 N08，QC 结果章节受 N03 影响 |
| A16 CI/发布镜像 | 未完成：workflow paths 仍不包含 deploy/**，检查仍使用 R release，未构建并测试正式发布 digest |
| A17 异步 | 注释已准确说明发起会话阻塞；只有请求序号门禁接入，完整来源、激活、会话结束与队列验收未闭合 |

## 5. 其他尚未完成的合同与目标环境验收

以下不是新增复现数量，但会影响“全部完成”的判断：

- **自动保存状态**：`wire_autosave()` 只在写入返回后发布成功/失败状态；失败对象把 `time` 设为 NULL，没有保留上次成功时间，也没有完整的 dirty/saving 状态。stderr 不等于容器删除后日志仍可追溯。
- **另存副本**：`store_save_project()` 仍直接保存同一个项目身份，没有实现方案约定的新 `project_id`、`derived_from` 与来源迁移。
- **来源验证完整性**：结构、项目归属和输入版本已补；但 `valid_result_ids()` 没有验证登记键与 analysis_id 一致、完整多节点无环性、绑定参数摘要对应实际保存参数等上一轮 C2 要求。不要把 helper 的部分通过扩大为完整合同验收。
- **CI 与发布**：需要覆盖 deploy/** 变更，并在最终发布镜像 digest 上跑关键业务链。目标服务器要验证用户 A/B 的文件和网络隔离、目录权限、重启恢复、日志持久化、磁盘满和备份恢复。现有本机双会话测试不能代替容器隔离。
- **范围决策**：PDF 已实际选择实现，报告也已选择扩充；方案 D1/D7 仍停留在旧决策文字，应与最终交付同步。D3 必须按 N10 处理。

## 6. 建议的下一步

1. **保留此次已通过的 Core 修复，不再重写整份计划。** 先把 N01–N06 的复现场景变成真实模块回归；N10 按既定口径收缩或完成统计确认。
2. 统一项目激活、QC 提交和分析完成的写入入口。来源在任务发起时确定；QC 与第二层 diff 的登记随其产物一起提交；页面、项目对象、自动保存和导出读同一份有效状态。
3. 修复 N07–N09 及保存/另存/CI 等未完合同，验证正常路径也能成功。只有“错误被拒”而没有“正常可完成”的测试不够。
4. 在新的固定提交上取得这些新增用例和完整回归的证据后，进入**正式镜像的部署测试**。至少覆盖“导入→QC 预览/应用→两层模型→富集/整合→保存→切项目→重开→独立脚本重放→报告”。

**当前提交尚不能签署为整改完成。** 已复现的参数重置、来源错配和脚本失败应先在代码层解决，部署环境的安装和连通性演练可以独立推进。
