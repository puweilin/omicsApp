# 整改方案 v3 复审

日期：2026-09-07。源码仍为 `f756751c24c12b01cf99fca399a6a5a17eeea4f4`，业务源码没有实施整改。
本轮完整阅读 v3 的 845 行，并运行候选函数与针对性验证。
方案 SHA-256：`c662c1df082194636b520d97b4d4a501e669298b8703602659a2b7b1c2525ce8`。

**结论：v3 有实质进步，架构方向可以认可；§1 数据合同还不能签字通过。建议补齐下述六项后进入批次 2–5 的实施。无需推翻按层预处理和统一有效输入入口的设计。**

批次 1 可以开始独立实施与测试，但 T1.4 需要补上现有归一化开关的处理，不能仅修改异常分支。候选方案通过设计评审与产品通过部署验收是两件事，后者仍需要真实修复、最终镜像和目标服务器证据。

为避免后续更新导致行号漂移，本报告引用 [本轮方案快照](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md)。另附 [复现输出](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-v3-probes.txt) 和 [验证脚本快照](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-v3-probes-source.txt)。

## 1. 已经改善、应当保留的部分

| 项目 | 本轮判断与证据 |
|---|---|
| R01 的参数裁剪回归 | 候选逻辑重新包含“拒绝实际不支持值 + 裁剪不支持的键”。按当前 t-test/lm 签名补齐示意能力表后，真实 `run_diff()` 默认调用均成功，实际 covariates 值被 t-test 拒绝 |
| R02 的系数混淆 | 新的内部项名、内部因子水平和数值 contrast 修正了上一轮已复现的组名 `v2` 冲突；与手工参考模型相比，effect 和 P 值最大差异均为 0 |
| 多模型项和恶意列名 | 抽测 ≥10 个模型项能够构造 contrast；恶意列名进入候选设计及真实 limma 拟合后，marker 未产生，数值与参考模型一致 |
| 按层保存 QC | `preprocessing[[tag]]` 能表达各层独立的已应用状态，解决了上一版单一全局槽位的结构限制 |
| 统一有效输入 | 明确将整合第二层纳入调用清单，修正了原先“整合只吃 diff bundle”的错误认识 |
| 失效归属与恢复 | 已认识到项目状态必须主动失效，以及 project_id 与一次打开/恢复操作不能混用；方向正确 |
| 导出与尺度范围 | 明确了不同步骤的返回值、offset、原文件与中间结果数值校验；取消临时新增 allow_raw 是合理收缩 |
| 发布与工期 | 删除“安全入口修完即可发布”的表述，撤销未建立接口依据的工期承诺，是正确调整 |

这些测试只覆盖所列场景，不能替代五个后端、连续/样条/ANOVA、配对模型及正式依赖环境的回归。T1.1 的 DESeq2 full/reduced 迁移仍需保持现有检验类型和数值语义；不能在修公式安全的同时无意改变分析方法。

## 2. §1 合同仍需补齐的六项

### V01｜P1：结果有效性没有完整校验请求，也没有定义传递失效

位置：[§1.4，第 137 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:137)、[T3.3，第 570 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:570)。

状态树声明了 `request_seq`，但 `binding_valid()` 根本不读取它。一个绑定到正确项目、激活版本和输入版本的旧请求，即使当前请求序号已从 1 变成 2，仍返回 TRUE。本轮直接执行候选函数得到了 `V3_OLD_REQUEST_ACCEPTED=TRUE`。参数选择发生变化但输入矩阵相同，也是这个问题的实际场景。

同时，`upstream_still_present()` 只有名称，没有定义。存在不等于仍然有效。候选 `invalidate_dependents()` 在整个循环中检查的是未更新的原项目：若按“上游 ID 仍存在”的语义补齐这个助手，B 层 diff 会被移入 retired，但依赖该 diff 的 integration 仍被保留。本轮以这种明确标注的助手实现验证，得到：

```text
V3_INVALIDATION_DIFF_B_RETIRED=TRUE
V3_INVALIDATION_INTEGRATION_RETAINED=TRUE
```

这不是对所有可能的助手实现的结论，而是说明：**只有将上游仍有效及传递失效写进合同，当前方案才能保证正确性。**

还存在两个相关缺口：binding 只记录一个 `tag/input_version`，整合却依赖两层；方案使用多个 `analysis_id`，但没有定义其生成、保存、查询和替换规则。当前 bundle 构造器也没有这种 ID，不能假设现成可用。

**必须补定：**

- 每个结果有不可变 analysis_id；引用稳定的上游结果 ID，并能找到这些结果及其数据来源。
- 多层分析保存按 tag 命名的全部输入版本，或者使用能证明覆盖全部输入的完整依赖图。
- 任务回调校验最新 request_seq、相关参数签名、激活令牌与全部输入版本。成功和失败回调同样校验。
- 上游失效后，所有依赖结果传递失效；采用拓扑校验、递归有效性校验或迭代至稳定的算法，不能只检查 ID 是否存在。
- `invalidate_dependents()` 的 `tag/new_version` 目前没有使用，`activation` 也不是显式形参。将调用合同写清楚，避免 Core 从不明确的外部作用域读取会话变量。

**新增验收：** 同输入连续两个请求的旧回调拒收；只更新第二层时整合失效；diff→enrich→integration 多级依赖失效；打乱 bundle 顺序结论不变。

### V02｜P1：会话激活校验与磁盘结果有效性需要分开

位置：[状态树第 64 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:64)、[恢复规则第 605 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:605)。

activation 被定义为不持久化，但它又进入会随 `.omp` 保存的 binding。重新打开时 activation 递增，所以保存结果的旧 binding 会被同一个 `binding_valid()` 判为无效。候选函数在 activation 从 1 变成 2 时返回 FALSE，这对旧的计算回调是正确的，对已经完成且随项目保存的结果却不能直接如此处置。

v3 提到了“明确文件恢复流程迁移绑定”，但还没有定义普通新格式文件如何恢复。只覆盖“旧格式迁移”不够。

**建议分成两类合同：**

1. **持久化结果来源**：项目/输入版本、完整参数、上游 analysis_id、处理记录和结果 ID。用于证明一个已保存结果对应哪份数据和分析。
2. **会话计算令牌**：session 标识、activation、request_seq。只用于决定一个正在进行的任务是否还允许提交结果。

打开文件时，先验证持久化来源与依赖链，再装载/恢复结果视图；为当前会话建立新的上下文。不要原样使用历史 activation，也不要不经验证就给所有历史结果贴上新令牌。另存副本还需明确：如何生成新 project_id、如何处理已完成结果的来源，以及如何拒绝旧会话任务写入副本。

**新增验收：** 正常保存的结果重开后仍能查看和导出；同时旧 activation 的晚到回调被拒绝；同 project_id 的旧保存真正生效；新格式、旧格式与另存副本分别验证。

### V03｜P1：输入版本少了特征注释，清洗结果也没有绑定其导入来源

位置：[effective_input，第 88 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:88)、[input_version，第 116 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:116)。

`input_version()` 包含矩阵、样本元数据、assay_type 和 provenance，却没有 `feature_df` 与 `omics_type`。本轮只改变一个 feature_symbol，摘要完全不变。但 feature_symbol 决定富集与整合的基因映射，这属于实际分析输入。

`preprocessing[[tag]]` 只保存处理后的版本，没有保存它基于哪个导入版本。`effective_input()` 只要看到 pp 存在就返回它。按当前合同构造状态：替换 `experiments$A` 后留着旧 pp，函数仍返回旧清洗数据。若实现希望靠所有写入者维护这一不变量，必须明确规定替换/删除层时原子清理 pp；更稳妥的是 pp 同时保存并校验 `base_input_version`。不能让任何一次恢复或替换遗漏成为静默读错数据。

还需补齐三处：

- **未应用 QC 时的缓存位置。** pp 不存在，`effective_input()` 每次调用都重新 hash。实测连续取三次输入就 hash 三次，与第 128–130 行“不在反应式反复计算”相矛盾。导入层也需缓存版本，不应为了缓存伪造一个已应用 QC 槽位。
- **provenance 的权威位置。** pp 有 `provenance`，其 `cleaned` 还可能有 `input$provenance`；hash 读取后者，导出读取前者。必须规定其中一份是唯一权威，另一份不保存或严格由它派生并校验。
- **导入去重的旧入口。** 当前 [mod_import_view.R:612](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:612) 仍先用 source_fingerprint 判断“无需重新导入”；[第 647 行](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:647) 的计算不含 sheet_roles。v3 虽要求记录 sheet_roles，却未明确改造这个早退判断。把它加入调用点清单，否则用户修改角色配置可能根本到不了新的版本机制。

**新增验收：** 单独改变 feature_symbol、样本注释、组学类型或解析角色，相关身份/去重判断相应变化；替换一层后不能返回旧 pp；重复查询输入不反复序列化大矩阵；导出使用的 provenance 与被 hash 的记录一致。

### V04｜P1：QC 应该从哪个版本重算尚未明确，示例调用还会报错

位置：[统一入口规则第 82 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:82)、[预览签名第 477 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:477)。

`preview_signature()` 绑定的是当前 effective input，但 `run_qc(qc_input, …)` 没有定义 qc_input 从哪里来。

若 QC 从导入层重算，签名必须对应导入层，而不是上次应用后的清洗层；若把“全部从 effective_input 取数”也用于 QC，则会在已经过滤/填补后的数据上再次 QC。此时放宽过滤阈值找不回删掉的特征，换填补方法也无法重新处理已经填好的缺失值。

本轮用真实 `run_qc()` 验证：严格过滤后再放宽，基于导入层重算得到 30 个特征，基于上次 cleaned 重算只有 29 个。这证明两种选择的语义不同，必须在合同中确定。

**建议本版固定为：每次 QC 预览从不可变导入层重新计算，应用替换这一层的 QC 结果，不隐式累计。** 下游分析仍统一使用 effective input。将“QC 从导入基线重算”作为统一入口规则中的显式例外，提供明确的基线读取函数；签名记录实际基线版本、解析后参数和当前激活令牌。

此外，第 488 行的调用：

```r
omicsCore::run_qc(qc_input, !!!params)
```

在这个普通 R 函数调用中没有参数注入上下文。直接运行得到 `invalid argument type`。可以使用：

```r
do.call(omicsCore::run_qc, c(list(input = qc_input), params))
```

本轮同样参数用后者调用成功。这是执行示例错误，不能以“只是示意”掩盖；如果保留伪代码，需要清楚标注与实际 API 的转换方式。

**新增验收：** 应用后放宽阈值能恢复应保留的特征；切换填补方法等于从同一导入基线直接运行对应方法；同参数重复应用不会累计过滤/重复堆叠 provenance；预览签名对应实际计算输入。

### V05｜P1：T1.4 只处理中途失败，现有开关仍允许原始强度进入 limma

位置：[缩范围第 405 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:405)、[延期范围第 798 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:798)。

取消新增 allow_raw 合理，但“没有新增继续原始分析的复选框”不代表不存在这条路径。当前 [mod_import_view.R:348](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:348) 已有 `Normalize on import` 复选框。取消勾选后，[第 536 行](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:536) 返回 method=none，[第 563 行](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_import_view.R:563) 跳过整个归一化分支，因此也不会进入新增的失败处理。

原始强度之后仍可交给 limma。用当前真实公共入口再次验证，得到 `effect_type=log2FC`、最大绝对 effect 约 869,311。数值具体大小依赖 fixture，问题是原始强度均值差被标成 log2 fold change。

**必须补一个最小的服务端限制，但不必本次设计完整的原始尺度分析功能：**

- 明确 raw intensity 是否只允许导入后做 QC；若允许，差异分析入口必须拒绝未经支持转换、却会被标为 log2FC 的模型输入。
- 如果产品规定蛋白组导入必须完成归一化，则提交端必须处理已有的 normalize=false，并同时考虑旧项目恢复与直接 Core 调用的边界。
- 完整的 method × analysis_type × assay scale 效果量扩展可以延期；不能把已知的“错误尺度被标为 log2FC”路径一并延期后宣布 A04 完成。TPM/FPKM 等仍开放的入口也要有明确支持或拒绝范围。

**新增验收：** 归一化失败、主动关闭归一化、旧项目保存着 raw intensity、直接调用 Core 四条路径分别有明确结果。不能仅凭 UI 隐藏或默认勾选证明服务端正确。

### V06｜P1：T4.3 漏写 GSEA 的完整结果保存要求

位置：[T4.3，第 690 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v3-reviewed.md:690)。

v3 只明确给 ORA 调用补 `output_p_cutoff=1`。GSEA 的“不要预筛特征”是正确要求，但与“保存全部通路结果”不同。

当前 [mod_enrich_view.R:110](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_enrich_view.R:110) 的 GSEA 调用没有传该参数，[enrich-gsea.R:81](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/enrich-gsea.R:81) 仍将默认 0.05 用作输出截断。通过真实 bundle 适配函数捕获送往底层的阈值，本轮确认：默认为 0.05，显式完整保存时为 1。此验证检查的是参数传递，不是重新执行完整 GSEA 通路计算。

**修订要求：** ORA 与 GSEA 两个分支都明确设置完整保存范围并记录 params；两者分别有真实计算路径的测试。GSEA 输入完整排名列表和输出完整通路表必须同时满足，不能相互替代。

## 3. 避免 v3 重写时遗失的原任务要求

v3 说“替换 v2”，因此以下要求需要留在新版任务中，不能依赖读者去拼接旧版：

| 项目 | 需要明确保留的内容 |
|---|---|
| T4.2 整合继承模型（A08） | 新版写了遗留 input 和“有意未校正”例外，但应继续明确正常路径：第二层继承协变量/配对设置，检查列与有效数据是否支持；不支持时停止，例外需明确选择并持久记录。只修 prune 不会自动修复第二层模型 |
| T5.5 报告范围（A15） | 目前主要写示例隔离。需明确报告是补齐参数、结果和主要图形，还是限定为“分析概要”并同步所有声明；两者应选定，而不是在任务标题里保留“报告”就算覆盖 |
| T5.6 async（A17） | 保留对现有同会话仍阻塞的注释/产品说明修正；队列与并发边界、会话结束后的回调规则也应列为实现/验收项。不能只增加版本字段 |
| RDS 部署开关 | 第 336 行仍写环境开关控制，随后又说多用户固定关闭不给开关。明确区分可信本地脚本与多用户网页入口；网页路径固定拒绝是最清楚的 v1 合同 |
| 失效历史 | §1.5 仍把“删除或保留历史”列为待定，T3.3 已选择累积 retired。若保留，明确唯一 ID、可见范围和容量/保留策略，避免每次重新分析使 autosave 无限积累大型模型 |

小的文字错误可以直接修：第 51 行说 §8 有初估，但实际没有；第 150 行说 §2 之后实施先过合同，与第 162 行批次 1 独立开工容易产生歧义；第 642 行把原函数描述为已经承担拒绝职责，基线实际是静默裁剪；第 622 行声称新增字段值都已在 joined 中，实际上 effect_type 和 analysis_id 仍需明确来源。它们不应被抬高为新的架构阻塞项。

## 4. 这次可以定下来的事项

**估工作废：同意。** 先明确上述生命周期和依赖合同，再按任务估算。可以先估独立入口修复的小任务，但不要沿用旧总工期作为上线承诺。

**D6：建议不新增 outlier_sd_threshold 控件。** 新建项目采用当前 `run_qc()` 的默认值 3，保存实际解析后的值并在适用位置展示。但是，加载一个已经保存了非默认值的项目时，应保留并展示原值，不能因没有 UI 控件就悄悄重置为 3。其他未暴露的 QC 参数也遵循这一原则。

**D5：暂缓通过具体合同，不否定架构方向。** 修改不需要再写一轮长篇复盘；建议直接补一页明确字段表和一张状态转移表：

| 操作 | 必须明确的状态变化 |
|---|---|
| 导入/替换/删除层 | 导入版本如何变化、pp 如何清理、哪些依赖结果失效、导入角色去重如何判断 |
| QC 预览 | 从导入基线计算；捕获实际基线版本与完整解析参数，不改已应用状态 |
| QC 应用 | 原子替换该层 pp；记录基线和有效版本；传递失效；同步视图与 autosave |
| 启动分析 | 生成唯一结果/请求标识；捕获全部输入、上游依赖和参数签名 |
| 回调完成/失败 | 校验最新会话令牌与请求序号；过期不写入、不改变当前错误状态 |
| 保存/打开/恢复 | 持久化来源与会话令牌分离；保存结果验证后恢复，旧在途任务仍拒收 |
| 另存副本 | 项目身份、已完成结果来源与当前会话关系明确更新 |
| 报告/导出 | 只读取来源可验证的有效结果及对应输入版本；不将当前数据强配给历史结果 |

完成 V01–V06、保留第 3 节任务要求，并把上述状态转换落实为可执行测试后，可以批准数据合同进入实施。**目前最合适的评价是“方向通过，合同补充后通过，部署尚未放行”。**
