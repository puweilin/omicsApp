# 整改方案 v4 可行性确认

日期：2026-09-07。源码：`f756751c24c12b01cf99fca399a6a5a17eeea4f4`，业务源码未实施整改。
评审文件共 719 行，SHA-256：`049bbd20584e840d910b8c9a293ae3816ecbda9003cbcfc0b4a8ce74c4849d1c`。

**结论：v4 的实施路线可行，架构可以保留；具体合同仍有三处必须闭合的接口问题，不能无条件签署 D5。建议开始批次 1，并将这三处作为 T2.0 的首个契约验证提交，不必再重写整份方案。**

另有两项实施约束：尺度拒绝的调用顺序/兼容范围，以及历史结果裁剪。它们可以在对应任务中落实，但验收时必须检查实际行为。正式部署仍需完成修复、正式镜像和目标服务器验收，文档通过不等于部署通过。

证据：[v4 方案快照](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md)、[执行输出](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-v4-probes.txt)、[验证脚本](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-v4-probes-source.txt)。本轮直接执行候选函数及少量真实 Core 调用；状态场景使用构造的项目对象，未将其描述为已实施后的完整 Shiny/浏览器验收。

## 已经验证有效的改进

| 改进 | 本轮结果 |
|---|---|
| feature_df 纳入输入摘要 | 单独改变 feature_symbol，版本发生变化 |
| pp 绑定导入基线 | 正常 pp 被选中；基线变化后不再返回旧 cleaned，返回新基线并标记 stale |
| 输入版本缓存 | 缓存已有时连续三次查询，hash 调用次数为 0 |
| 保存来源独立于 activation | 正常持久化来源判断不再读取会话 activation |
| request_seq 校验 | 序号较旧的请求被候选回调函数拒绝 |
| 不动点失效 | 完整依赖图中，更新 B 层会移除 B 的 diff、enrich 和相关 integration，仅保留不受影响的 A 层 diff；打乱存储顺序结论相同 |
| QC 基线、参数调用、GSEA 保存 | v4 明确采用基线重算与 do.call，并为 ORA/GSEA 两支设置完整输出范围，回应了上一轮对应缺口；这些仍需实施后的真实路径测试 |

T1.1 保留 v3 方案合理，前轮验证的数值结果继续作为所列场景的证据。v4 第 301–302 行明确保留完整后端回归要求，避免把抽测扩大成全框架认证，这一点正确。

## 必须闭合的三个接口

### C1｜P1：回调必须同时验证会话令牌和实际分析来源

位置：[callback_admissible，第 177 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:177)、[回调状态转换，第 272 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:272)。

两类绑定拆开是正确的，但写入时不能只验证其中一类。当前 callback_admissible 没有检查输入版本或上游结果；替换层/QC 应用也没有按状态表递增 activation 或 request_seq。所以同一项目、同一会话、同一请求，在输入已经变化时仍可通过。

直接执行候选函数：

```text
V4_INPUT_CHANGED_CALLBACK_ACCEPTED=TRUE
V4_OLD_RESULT_SOURCE_CURRENT=FALSE
```

这说明“令牌仍匹配”和“产物仍属于当前数据”可以同时得出相反结论。它是候选合同的可复现场景，不代表已经证明当前仍阻塞同会话的 Shiny 实现允许用户在所有在途任务中切数据。

还有一个确定的字段错配：第 164–170 行的 session_token 没有 params_signature，但函数第 182 行检查它。ctx 有正常参数摘要时，按字段表生成的令牌会被拒；两侧都缺这个字段时，`identical(NULL, NULL)` 又会让此项校验通过。摘要还应按 view/request 管理，不能将所有模块共用一个不明确的 ctx$params_signature。

**可直接补入 T2.0/T5.6 的合同：**

> 每个任务在启动时捕获完整来源绑定与会话令牌。成功/失败回调写入前，必须同时满足：会话未结束、项目与激活版本一致、该视图请求序号最新、解析参数摘要一致、全部输入版本当前有效、全部上游依赖仍有效。检查和提交在一个受控操作中完成。

具体实现可以合成一个接收 token、来源草稿、ctx、project 的 admission 函数，也可以显式组合现有函数。两个函数概念上的分离应保留，写入检查必须取交集。缺必需字段应拒绝，不能让 NULL 与 NULL 相等代替有效性证明。

**验收：** 在没有启动新请求的情况下改变主层/QC/第二层/上游结果，原任务成功和失败回调均不写入；正常未变任务成功；缺参数摘要拒绝；diff 与 enrich 的请求序号互不串用。

### C2｜P1：持久化来源校验还缺项目归属、结构校验及旧文件处理

位置：[result_inputs_current，第 186 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:186)、[valid_result_ids，第 211 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:211)。

当前函数只比较输入版本，不检查 pb$project_id。两个项目使用相同数据时，另一个项目的来源绑定会通过。空绑定也会通过：`all(logical(0))` 返回 TRUE。

```text
V4_FOREIGN_PROJECT_SOURCE_ACCEPTED=TRUE
V4_EMPTY_SOURCE_ACCEPTED=TRUE
```

而将一个没有 provenance_binding 的旧 bundle 放进 valid_result_ids 时，vapply 提取 analysis_id 会抛出长度错误。这个行为与“旧 .omp 加载不报错”和“未知归属不当作当前结果”的要求均未接通。

**可直接补入 T2.0/T3.4 的合同：**

- 先验证绑定结构：唯一且非空的 analysis_id、非空 project_id、具名且有效的输入版本映射、合法上游 ID、参数摘要；按分析类型检查应有的输入及依赖，不能把空映射当作已验证。
- 必须满足 `pb$project_id == project$metadata$project_id`；另存副本仅在明确迁移操作中更新归属。
- 验证记录的参数摘要对应结果实际保存的解析参数；恢复历史结果不要求它等于当前页面尚未应用的草稿参数。
- 旧文件在调用通用有效性算法前迁移。来源不够证明时作为带原因的历史/待验证记录展示，不能通过直接填入当前来源“使其有效”。
- 校验 analysis_id 唯一、依赖图无自环/循环；不动点算法负责传递失效，不代替结构验证。

**验收：** 相同数据但错误 project_id 被拒；空/缺字段绑定有结构化结果；旧 .omp 能打开且不伪造来源；正常同项目保存结果仍可恢复；另存副本归属迁移后有效。

### C3｜P1：依赖图需要能查到整合内部生成的第二层差异结果

位置：[依赖算法，第 214 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:214)，与 [mod_integration_view.R:224](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/mod_integration_view.R:224)。

不动点算法本身在完整图中有效，但它只从 project$bundles 取上游 ID。当前实际链路是：第二层 sec_diff 在 worker 中计算，传给 run_integration 后丢在局部作用域；回调只返回 integration。[app_server.R:174](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsApp/R/app_server.R:174) 也只回写固定视图槽位。

如果 integration 的 upstream_ids 正确记录了两层差异 ID，而项目只持有当前主 diff 和 integration，新算法反而会把刚完成的整合判为无效，因为找不到 secondary diff。构造与当前保存形状相同的对象，得到：

```text
V4_UNREGISTERED_SECONDARY_INTEGRATION_VALID=FALSE
V4_REGISTERED_SECONDARY_INTEGRATION_VALID=TRUE
```

**可直接补入 T2.0/T3.2/T4.2 的合同：**

> 所有被 upstream_ids 引用的分析，必须有按 analysis_id 可查询、可持久化的来源记录。整合内部的 secondary diff 与最终整合结果应作为一个提交单元登记；视图当前选中结果与分析记录集合分开管理，切换视图不能删除仍被有效结果引用的节点。

可以保留现有视图槽位作为索引，再新增来源记录索引；也可以调整 bundles 的索引方式。无需强制保存全部大型模型对象，但至少保存可验证的来源、解析参数及脚本重放所需信息。不能靠删掉 upstream_ids 来躲过依赖校验。

此外，Core 的 [run-integration.R:188](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/run-integration.R:188) 当前仍从 project$experiments 填入样本数、特征数、尺度等 input_info。接通有效输入时应同步这条元数据路径，保证记录反映实际参与分析的数据，不能只换 UI 的第二层输入。

**验收：** 整合刚完成即在有效集合中；保存重开后两条上游来源都可查；切换 diff 视图不误删被引用的节点；任一真实上游变化后整合失效；input_info 数量与实际有效输入一致。

## 两项必须写进实施验收的约束

### I1：尺度检查必须在 auto 解析后调用，lm 的范围变化要准确说明

位置：[T1.4，第 343 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:343)。

候选检查能够拒绝显式 limma + raw_intensity，这项有效。但 method=auto 会直接被放行，而蛋白组 raw intensity 的 auto 随后会选择 limma。

这不证明一个尚未实施的 v4 已经存在绕过；它说明**调用位置必须写死**：在 [run-diff.R:213](/Users/puweilin/SciProject/CHISSS/omicsApp/packages/omicsCore/R/run-diff.R:213) 的 auto/fallback 最终解析后、后端执行前校验。验收加入 raw intensity + 默认 auto，不能只测显式 limma。

将所有 lm 线性输入一概拒绝，是额外的功能范围收缩。当前 lm 对 RNA raw_count 已自动 `log2(x+1)`，真实公共调用可以成功；它的组间 effect_type 是 beta，连续分析是 correlation，并不都报告 log2FC。因此错误文案“method='lm' reports a log2 fold change”不准确。对 RNA 建议 normalize_omics() 也不可执行，因为该函数只支持 proteomics。

建议保留已有正确的 lm raw_count 转换路径，按解析后的后端和实际变换检查尺度；若产品明确选择删除该支持，则将其作为兼容性变更更新 API 文档、方法列表和测试，并提供实际可行的替代路径。修复 limma 错误尺度不要求顺带删除已经正确处理 counts 的功能。

### I2：历史裁剪策略合理，但展示的 retire_invalid 尚未执行它

位置：[retire_invalid，第 226 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:226)、[容量策略第 258 行](/Users/puweilin/SciProject/CHISSS/omicsApp/docs/audits/2026-09-07/evidence/remediation-plan-v4-reviewed.md:258)。

“只保留来源与摘要，最多 20 条”的设计可以采用。但展示的函数只是 `c(retired, bundles[gone])`；运行 21 次会得到 21 个完整旧 bundle，并且重复使用 diff 作为键。实际差异模型常在 `results$diff_object`，不能只删除一个叫 model_object 的字段。

实现时需要显式的压缩/限额步骤：按允许字段创建历史摘要，以 analysis_id 唯一索引，记录退休时间和原因，按实际时间淘汰至 20 条。用户可见的历史应与有效依赖记录区分，不能因为历史显示上限而误删仍被有效分析引用的来源节点。

**验收：** 连续退休 21 次后不超过 20 条，键唯一，大型模型字段不在保存对象中；保存/恢复后限额仍成立。该项属于落实已写明策略，不要求重新选择架构。

## 建议的推进方式

1. **可以让 Claude Code 开始批次 1。** T1.1 继续执行完整模型/安全回归；T1.4 明确 auto 之后校验及 lm 兼容范围。
2. **将 C1–C3 作为 T2.0 的第一组可执行契约测试和实现。** 先用本轮场景证明缺口，再实现通过。D5 在这组测试通过、接口一致后签字，不以再次口头接受意见替代验证。
3. **通过后继续批次 2–5，按功能链验证。** I2 并入历史处理测试；D1、D4、D7 在对应功能交付前定案。D2 固定拒绝网页 RDS、D6 保留隐藏参数原值的决定可以保留。
4. **最后在正式发布镜像和目标服务器取得放行证据。** 代码与文档的任何局部通过都不能替代这一阶段。

**本轮不建议再做一次全篇 v5 重写。将三个接口缺口和两个实施约束直接并入对应任务，转向可执行测试与代码验证，会比继续改长篇方案更有效。** 本轮仅增加评审文档、快照及证据，没有修改业务代码或用户方案。
