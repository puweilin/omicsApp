# omicsApp 整改实施方案（v4）

源码基线：`f756751`
上游：[审查报告](README.md)（17 项）→ [方案复审](REMEDIATION-REVIEW.md)（R01–R10）→ [v3 复审](REMEDIATION-V3-REVIEW.md)（V01–V06）
本版变更：按 V01–V06 补齐 §1 数据合同；架构方向不变。

---

## 0. 对 v3 复审的处置

V01–V06 全部接受。V03/V05/V06 的源码缺口我已独立确认（§9），
V01/V02/V04 是对我自己 §1 的阅读，逐条核对成立——其中 V02 是 v3 的**自相矛盾**：
`activation` 被声明为不持久化，却进入了随 `.omp` 保存的 `binding`，
而 `binding_valid()` 要求它相等。按 v3 原样实现，**任何被恢复的项目，
它自己保存的分析结果都会被判为无效**——恰好摧毁 T3.3 要保住的东西。

| 编号 | 缺口 | 本版处置 |
|---|---|---|
| V01 | `binding_valid()` 不读 `request_seq`；`upstream_still_present()` 只有名字；传递失效未定义 | §1.4 拆分绑定 + §1.5 不动点算法 |
| V02 | 会话令牌与持久化来源混为一谈 | §1.4 拆成两类；§1.7 定义打开流程 |
| V03 | `input_version` 缺 `feature_df`/`omics_type`；pp 未绑定导入基线；未应用 QC 时反复 hash；provenance 两处；导入去重早退绕过 | §1.3 重写 + §1.2 调用点补 |
| V04 | QC 从哪个版本重算未定；`!!!` 语法错误 | §1.2 显式例外 + T3.1 修正 |
| V05 | `Normalize on import` 取消勾选可绕过 T1.4 | T1.4 增加模型入口拒绝 |
| V06 | GSEA 输出仍按 0.05 截断 | T4.3 两个分支都补 |

v3 复审 §3 点出的"重写时遗失的原任务要求"（第二层模型继承、报告范围决策、
异步队列边界、RDS 开关口径、retired 容量策略）已全部写回对应任务。
§3 末尾的 4 处文字错误已改正——包括我把基线 `prune_backend_args()` 描述为
"承担拒绝职责"，实际基线是**静默裁剪**，拒绝是本次新增的。

### 0.1 里程碑口径

| 节点 | 含义 | 不含义 |
|---|---|---|
| 批次 1 完成 | 安全入口整改完成，可进入后续验证 | 不是可发布，不是生产就绪 |
| §1 合同评审通过 | 批次 2–5 可分派 | 批次 1 **不等待**此项，已可开工 |
| 批次 2–5 完成 | 功能修复完成，可进入镜像验收 | 未取得目标服务器证据 |
| 批次 6 完成 | 取得放行证据 | —— |

**任何时点都不得据本方案放行部署。** 估工待 §1 评审通过后按任务重估；
批次 1 的四个任务可先单独估。v2 的"24–26 人日"作废。

---

## 1. 数据合同（完整定义）

### 1.1 状态树

```text
project
├─ metadata
│   ├─ project_id          持久。另存副本时重新生成，并记 derived_from
│   └─ schema_note         旧 .omp 迁移标注
├─ experiments[[tag]]      导入层，不可变基线
│   └─ version             ← 导入/替换时算一次并缓存（V03：不再每次 hash）
├─ preprocessing[[tag]]    按层的已应用预处理；未应用则不存在
│   ├─ base_input_version  ← 基于哪一版导入层（V03：不靠写入者维护不变量）
│   ├─ input_version       该 pp 产出的有效输入版本
│   ├─ cleaned             omics_input；provenance 存在这里（唯一权威）
│   └─ params              解析后的完整 QC 参数（含经 ... 传入的）
├─ bundles[[name]]         当前有效结果
├─ retired[[analysis_id]]  已失效结果（见 §1.6 容量策略）
└─ （会话内，不入 .omp）
    ├─ activation          每次 open / restore / 首次 import / 另存 递增
    └─ request_seq[[view]] 每次发起分析递增
```

### 1.2 唯一入口与它的一个例外

**下游分析统一读有效输入**：

```r
# 差异分析、整合的第二层差异、报告、导出脚本全部经此获取。
# 任何模块都不得再直接读 proj$experiments[[tag]] 做分析。
effective_input <- function(project, tag) {
  raw <- project$experiments[[tag]]
  pp  <- project$preprocessing[[tag]]
  # pp 必须证明自己基于当前的导入层。层被替换而 pp 未清理时,
  # 静默返回旧 cleaned 就是读错数据 —— 宁可退回基线并告警。
  if (!is.null(pp) && identical(pp$base_input_version, raw$version)) {
    return(list(input = pp$cleaned, version = pp$input_version,
                source = "qc_cleaned", provenance = pp$cleaned$provenance,
                stale = FALSE))
  }
  list(input = raw, version = raw$version, source = "imported",
       provenance = raw$provenance, stale = !is.null(pp))
}
```

**例外：QC 自己从导入基线重算**（V04）。

```r
# QC 是唯一不读 effective_input() 的消费者。
# 在已过滤/已填补的数据上再跑一次 QC,放宽阈值找不回删掉的特征,
# 换填补方法也无法重新处理已经填好的缺失值 —— 复审实测：
# 严格过滤后放宽,从导入层重算得 30 个特征,从上次 cleaned 重算只有 29 个。
baseline_input <- function(project, tag) {
  raw <- project$experiments[[tag]]
  list(input = raw, version = raw$version)
}
```

**每次应用 QC 都是替换该层的 pp，不叠加。**

#### 调用点清单（缺一个就是 R03/V03 复发）

| 位置 | 现状 | 改后 |
|---|---|---|
| `mod_diff_view.R:75` `active()` | `exps[[tag]]` | `effective_input()` |
| `mod_integration_view.R:82` `sec <- exps[[secondary_idx]]` | 导入层 | `effective_input()` |
| `export-script.R:216` 导出循环 | `experiments[[tag]]$provenance` | `effective_input()$provenance` |
| `default-report.Rmd` | 项目概要 | `effective_input()` |
| `mod_qc_view.R` 预览与应用 | `exps[[tag]]` | **`baseline_input()`**（例外） |
| `mod_import_view.R:612` `fingerprints_match()` 早退 | 指纹不含 `sheet_roles` | 见 §1.3 |

### 1.3 输入版本（V03）

```r
# 补 feature_df 与 omics_type：feature_symbol 决定富集的基因映射和整合的
# 连接键,改了它分析输入就变了 —— v3 的摘要对此完全不变。
input_version <- function(input) {
  digest::digest(list(
    matrix     = input$expr_mat,
    samples    = colnames(input$expr_mat),
    features   = rownames(input$expr_mat),
    meta       = input$meta_df,
    feature_df = input$feature_df,     # ← 补
    omics_type = input$omics_type,     # ← 补
    assay_type = input$assay_type,
    provenance = input$provenance
  ), algo = "xxhash64")
}
```

**缓存位置**（V03）：算出的版本存回 `experiments[[tag]]$version` 与
`preprocessing[[tag]]$input_version`，只在导入 / 替换 / 应用 QC / 打开项目时计算。
`effective_input()` 只读字段，不 hash。**不为了缓存伪造一个已应用 QC 槽位。**

**provenance 唯一权威**（V03）：存在 `omics_input$provenance` 上（导入层是
`raw$provenance`，清洗层是 `pp$cleaned$provenance`）。`preprocessing` **不再单独
持有 `provenance` 字段**——v3 两处并存会漂移。hash 与导出读同一处。

**导入去重**（V03）：`input_fingerprint()` 补 `sheet_roles`，
否则用户改了角色配置、重新导入同一文件，会在 `mod_import_view.R:612`
被"那是已经加载的文件"挡下，**根本到不了新的版本机制**。

### 1.4 两类绑定（V01 + V02）

v3 把"这个结果对应哪份数据"和"这个回调还允许写入吗"混成一个 `binding`，
于是不持久的 `activation` 进了持久化结构。拆开：

```r
# 一、持久化来源。随 .omp 走,用于证明一个已保存结果对应哪份数据和分析。
provenance_binding <- list(
  analysis_id = "<不可变 ID,本结果的身份>",
  project_id  = "<所属项目>",
  inputs      = list(proteomics = "<input_version>",     # 按 tag 的映射
                     rnaseq     = "<input_version>"),    # 整合有两个
  upstream_ids = c("<diff 的 analysis_id>", ...),
  params_signature = "<解析后全部参数的摘要>"
)

# 二、会话计算令牌。不入 .omp,只决定一个在途任务是否还允许提交。
session_token <- list(
  session_id  = "<Shiny session>",
  project_id  = "<发起时的项目>",
  activation  = 3L,
  request_seq = 7L,
  view        = "diff"
)
```

两个不同的判断，用两个函数：

```r
# 回调是否允许写入（会话内）。成功与失败回调都要过（V01/T5.6）。
callback_admissible <- function(token, ctx) {
  identical(token$session_id, ctx$session_id) &&
  identical(token$project_id, ctx$project_id) &&
  identical(token$activation,  ctx$activation) &&
  identical(token$request_seq, ctx$request_seq[[token$view]]) &&   # ← v3 漏了
  identical(token$params_signature, ctx$params_signature)
}

# 已保存结果是否仍有效（打开文件时）。只看来源,不看会话令牌。
result_inputs_current <- function(pb, project) {
  all(vapply(names(pb$inputs), function(tg) {
    identical(pb$inputs[[tg]], effective_input(project, tg)$version)
  }, logical(1)))
}
```

`analysis_id` 现有 bundle 构造器**没有**（V01），需在 `new_analysis_bundle()` 增加，
生成方式与 `project_id` 同（不用 `sample()`，见 T2.0）。

### 1.5 传递失效：迭代至不动点（V01）

v3 的 `invalidate_dependents()` 单遍扫描且检查的是未更新的原项目。
复审用一个按"上游 ID 仍存在"语义补齐的助手实测：

```text
V3_INVALIDATION_DIFF_B_RETIRED=TRUE
V3_INVALIDATION_INTEGRATION_RETAINED=TRUE
```

——B 层 diff 被退休了，依赖它的 integration 却留着。**存在 ≠ 仍然有效。**

```r
# 有效性是递归的。单遍扫描会让 integration 活在已经退休的 diff 上,
# 所以迭代到不动点；bundle 的存储顺序不得影响结论。
valid_result_ids <- function(project) {
  ids <- names(project$bundles)
  repeat {
    alive <- vapply(project$bundles[ids],
                    function(b) b$provenance_binding$analysis_id, character(1))
    keep <- Filter(function(nm) {
      pb <- project$bundles[[nm]]$provenance_binding
      result_inputs_current(pb, project) && all(pb$upstream_ids %in% alive)
    }, ids)
    if (identical(keep, ids)) return(ids)
    ids <- keep
  }
}

# 调用合同显式化：不从外部作用域读会话变量（V01）
retire_invalid <- function(project) {
  keep <- valid_result_ids(project)
  gone <- setdiff(names(project$bundles), keep)
  project$retired  <- c(project$retired, project$bundles[gone])
  project$bundles  <- project$bundles[keep]
  project
}
```

**失效必须同时作用于四处**：内存项目 `bundles`、视图 `reactiveVal`、
autosave 快照、导出/报告的读取集合。

### 1.6 字段表

| 字段 | 类型 | 入 `.omp` | 权威 | 写入者 |
|---|---|:---:|---|---|
| `metadata$project_id` | chr | ✓ | 自身 | 创建 / 另存 |
| `metadata$schema_note` | chr | ✓ | 自身 | 旧文件迁移 |
| `experiments[[tag]]` | omics_input | ✓ | 自身 | 导入 |
| `experiments[[tag]]$version` | chr | ✓ | `input_version()` | 导入 / 替换 |
| `experiments[[tag]]$provenance` | list | ✓ | 自身 | `read_omics` / `normalize_omics` |
| `preprocessing[[tag]]$base_input_version` | chr | ✓ | 复制自 `experiments$version` | 应用 QC |
| `preprocessing[[tag]]$input_version` | chr | ✓ | `input_version(cleaned)` | 应用 QC |
| `preprocessing[[tag]]$cleaned` | omics_input | ✓ | 自身 | 应用 QC |
| `preprocessing[[tag]]$cleaned$provenance` | list | ✓ | **唯一权威** | `run_qc` |
| `preprocessing[[tag]]$params` | list | ✓ | 由 provenance 末步派生 | 应用 QC |
| `bundles[[n]]$provenance_binding` | list | ✓ | 自身 | 分析完成 |
| `bundles[[n]]$session_token` | list | ✗ | 自身 | 发起分析 |
| `retired[[id]]` | list | ✓ | 自身 | `retire_invalid()` |
| `activation` | int | ✗ | 会话 | open/restore/首次 import/另存 |
| `request_seq[[view]]` | int | ✗ | 会话 | 发起分析 |

**`retired` 容量策略**（v3 复审 §3）：只保留 `provenance_binding` +
结果摘要表，**丢弃 `model_object` 等大对象**；每项目上限 20 条，超出按时间淘汰。
否则反复重算会让 autosave 无限膨胀。

### 1.7 状态转移表

| 操作 | 状态变化 |
|---|---|
| **导入新层** | 计算 `version` 并缓存；`provenance` 起头；`activation++`（若是本会话首个真实项目） |
| **替换层** | 新 `version`；**原子清理该 tag 的 `pp`**；`retire_invalid()` |
| **删除层** | 移除 `experiments[[tag]]` 与 `pp[[tag]]`；`retire_invalid()` |
| **QC 预览** | 从 `baseline_input()` 计算；签名记 `base_input_version` + 解析后参数 + `activation`；**不改已应用状态** |
| **QC 应用** | 验证预览成功且签名仍匹配 → 原子替换 `pp[[tag]]`（记 `base_input_version` 与 `input_version`）→ `retire_invalid()` → 同步视图与 autosave |
| **发起分析** | 生成 `analysis_id`；`request_seq[[view]]++`；捕获 `session_token` 与 `provenance_binding` 草稿 |
| **回调（成功/失败）** | `callback_admissible()` 全中才写入；过期则丢弃，**且不改变当前错误状态** |
| **保存** | 写 `bundles`（含 `provenance_binding`，剥掉 `session_token`）、`preprocessing`、`retired` |
| **打开 / 恢复** | ①读文件 ②对每个 bundle 跑 `result_inputs_current` + 递归上游 ③有效者装载、无效者移入 `retired` 并告知 ④`activation++` ⑤**不给历史结果贴新令牌**（它们本就没有）⑥在途旧回调因 `activation` 不匹配被拒 |
| **另存副本** | 新 `project_id`，记 `derived_from`；已完成结果的 `provenance_binding$project_id` 一并迁移；旧会话任务因 `token$project_id` 不匹配被拒 |
| **报告 / 导出** | 只读 `valid_result_ids()` 之内的结果及其 `provenance_binding$inputs` 对应的输入版本；**不把当前数据强配给历史结果** |

### 1.8 合同评审出口

评审需确认：`preprocessing` / `retired` 的 `.omp` 兼容策略、
`input_version()` 性能预算（见 T2.0 验收）、`activation` 的语义边界、
`is_significant` 拆列后旧 `.omp` 的读取语义。

---

## 2. 批次 1 —— 独立入口修复（不等合同，可立即开工）

### T1.1 设计矩阵与 contrast（A01 / R02）

**v3 版本已通过复审验证**（组名冲突时 effect 与 P 值和参考模型最大差异均为 0，
≥10 个模型项能正确构造 contrast，恶意列名无副作用且数值正确），**本版不改**。

要点保留：内部项名 `t001…`、内部因子水平 `L001…`、设计列名全部程序生成且唯一、
**contrast 按系数位置构造成数值向量**（弃用解析字符串的 `makeContrasts()`）、
用户标签只进显示映射。

逐分支处理表（limma 两组/ANOVA/连续/样条、lm、edgeR、DESeq2、`run-diff.R:399`）
不变。**复审特别提示**：DESeq2 的 full/reduced 迁移必须保持现有检验类型与
数值语义，`lm()` 周边逻辑不得直接换 `lm.fit()`——修公式安全的同时不能改变分析方法。

复审的测试只覆盖所列场景，**不能替代五个后端 × 连续/样条/ANOVA/配对模型
在正式依赖环境下的完整回归**。

### T1.2 网页端关闭 RDS（A02 / R10 / v3 复审 §3）

**口径统一**（v3 措辞矛盾，本版定死）：

- **网页上传路径：固定拒绝，不提供开关。** 这是 v1 的合同
- `read_omics(allow_r_objects = TRUE)` 只保留给**可信的本地脚本**调用
- 环境变量 `OMICSAPP_ALLOW_RDS` **取消**

文案：加载并使用不可信 R 对象可以执行代码（触发点是字段访问与方法派发）；
不把所有 `readRDS()` 一概说成立即执行任意代码。

测试须在**真实上传解析入口**证明拒绝发生在读取对象之前。

### T1.3 项目路径服务端校验（A11 / R10）

不变。要点：白名单用 `list_saved_projects()`（**不是** `store_list_projects()`）、
`normalizePath` 后校验在 store 内、越界与不存在统一回 `"Unknown project."`、
校验后文件消失单列结构化错误。测试补：目录外真实 `.omp`、绝对路径、软链接、
**正常中文项目名**、校验后文件消失。

### T1.4 归一化绕过路径（A04 / **V05**）

#### v3 的缺口

v3 只改了归一化**失败**分支。但 `mod_import_view.R:348` 有
`Normalize on import` 复选框（默认勾选）；取消勾选后 `:536` 的
`pending_normalize_method()` 返回 `"none"`，`:563` 的 `do_normalize` 为 FALSE，
**整个归一化分支被跳过，新增的失败处理根本不会执行**。
原始强度照样进 limma，`effect_type` 仍标 `log2FC`（复审实测最大 |effect| ≈ 869,311）。

#### 改法：一个服务端检查收敛四条路径

导入侧保留 v3 的失败中止。**新增模型入口的尺度拒绝**——这是四条绕过路径
（归一化失败 / 主动关闭 / 旧项目存着 raw / 直接调 Core）的共同收口：

```r
# 名字要与实际范围相符。v3 的 assert_log_scale() 只挡 raw_intensity,
# 而 LOG_SCALE_ASSAY_TYPES 只有 5 个 —— tpm / fpkm / raw_count 一个没挡,
# 名字却宣称在保证 log scale。
assert_effect_scale_supported <- function(input, method, analysis_type) {
  if (!method %in% c("limma", "lm")) return(invisible(TRUE))     # 各后端各自的合同
  at <- input$assay_type
  if (at %in% LOG_SCALE_ASSAY_TYPES) return(invisible(TRUE))
  stop(sprintf(
    "method = '%s' reports a log2 fold change, but assay_type '%s' is on a linear scale. Normalize the layer first (normalize_omics()), or choose a backend whose effect definition matches this scale.",
    method, at), call. = FALSE)
}
```

**范围声明（必须写进 README 与页面）**：v1 只支持**已在对数尺度**的输入做
limma/lm 差异分析。`tpm` / `fpkm` / `raw_count` / `raw_intensity` 进 limma
一律拒绝并给出可执行的下一步。完整的 method × analysis_type × assay scale
效果量扩展延期，**但已知的"错误尺度被标为 log2FC"路径不随之延期**。

#### 验收（V05）

四条路径分别有明确结果，**不能仅凭 UI 隐藏或默认勾选证明服务端正确**：

| 路径 | 预期 |
|---|---|
| 归一化失败 | 导入中止，可重试 |
| 主动取消 `Normalize on import` | 可导入（QC 仍可做），但 limma 差异分析被拒 |
| 旧项目保存着 `raw_intensity` | 打开可查看，差异分析被拒并说明 |
| 直接 `omicsCore::run_diff()` | 同样被拒 |

---

## 3. 批次 2 —— 数据合同实现

### T2.0 按层预处理、版本与绑定（R03/R04/R06/R07 + V01–V03 的共同地基）

1. `omics_project` 增加 `preprocessing`、`retired`、`metadata$project_id`；
   `load_project()` 对旧 `.omp` 补齐并标注 `schema_note`
2. 实现 §1.2–§1.5 全部函数：`effective_input()`、`baseline_input()`、
   `input_version()`、`callback_admissible()`、`result_inputs_current()`、
   `valid_result_ids()`、`retire_invalid()`
3. `new_analysis_bundle()` 增加 `analysis_id`（V01：现有构造器没有）
4. `project_id` / `analysis_id` **不用 `sample()`**（R07：消耗分析 RNG 流，
   且在相同 `set.seed()` 下重复）；用 `uuid::UUIDgenerate()` 或
   时间戳 + `basename(tempfile(""))`
5. `activation` / `request_seq` 会话计数器
6. `subset_omics()` 保留 `provenance` / `source_path` / `source_fingerprint`
   （实测现在全丢，§9）
7. `input_fingerprint()` 补 `sheet_roles`；改造 `mod_import_view.R:612` 早退

**验收**：

- 两层各自持有独立 `pp`，互不覆盖
- 只改 `feature_symbol` → `input_version` 变化（v3 不变）
- 替换某层后不返回旧 `pp`（`stale = TRUE` 退回基线）
- 连续三次取同一输入 → **只序列化一次**大矩阵
- `.omp` round-trip：`preprocessing` / `provenance_binding` / `retired` 完好
- 旧 `.omp` 加载不报错且 `schema_note` 可见
- **性能预算**：63k × 258 的 `input_version()` < 0.5 s（`test-perf-budget.R`）
- 打乱 `bundles` 顺序，`valid_result_ids()` 结论不变（V01）

---

## 4. 批次 3 —— QC、失效与项目激活

### T3.1 QC 预览/应用（A03 前置 / R05 / **V04**）

**基线固定**：预览与应用都从 `baseline_input()` 重算，**不在已清洗数据上叠加**。
复审实测两种语义确实不同（严格过滤后放宽：从导入层重算得 30 个特征，
从上次 cleaned 重算只有 29 个）。

**预览自带签名，应用只验证**（R05：`mod_qc_view.R:238` 失败时保留旧
`last_bundle`，v3 的应用按钮会给它贴当前层身份）：

```r
    preview_signature <- function(tag, params) {
      list(project_id = current_project()$metadata$project_id,
           activation = activation(),
           tag        = tag,
           base_input_version = baseline_input(current_project(), tag)$version,
           params     = params)
    }

    shiny::observe({
      a <- active(); params <- collect_qc_params(input, a)
      sig <- preview_signature(a$tag, params)
      base <- baseline_input(current_project(), a$tag)$input
      # v3 写的 run_qc(qc_input, !!!params) 是 rlang 注入语法,
      # 在普通调用里无效 —— 实测 invalid argument type。
      b <- tryCatch(do.call(omicsCore::run_qc, c(list(input = base), params)),
                    error = function(e) e)
      preview(if (inherits(b, "error"))
                list(bundle = NULL, sig = sig, ok = FALSE, error = conditionMessage(b))
              else list(bundle = b, sig = sig, ok = TRUE))
    })

    shiny::observeEvent(input$apply_to_downstream, {
      p <- preview(); a <- active()
      if (!isTRUE(p$ok)) { notify_error("当前设置没有跑通，修正后再应用。"); return() }
      if (!identical(p$sig, preview_signature(a$tag, collect_qc_params(input, a)))) {
        notify_warn("设置在预览之后变过了，请重新运行再应用。"); return()
      }
      apply_preprocessing(a$tag, p$bundle, p$sig)   # 原子替换 pp + retire_invalid
    })
```

**控件恢复**：`activation` 或层变化时执行一次明确恢复（已应用状态优先于当前草稿），
之后才允许编辑。页面三态分列：已保存参数 / 当前草稿 / 已应用状态。

**D6 已定**：不新增 `outlier_sd_threshold` 控件（QC UI 实际只有 3 个输入控件）。
新建项目用 `run_qc()` 默认值 3；**加载已保存非默认值的项目时必须保留并展示原值，
不得因没有控件就重置为 3**（v3 复审补充）。其他未暴露的 QC 参数同此原则。

**验收**（R05 + V04）：
"成功后失败再应用"拒绝、"切层后预览失败再应用"拒绝、
"A 项目 knn → B 项目 min"显示 min、"已保存方法在当前环境不可用"明确提示不静默改选、
**应用后放宽阈值能恢复应保留的特征**、换填补方法等价于从同一基线直接运行、
同参数重复应用不累计过滤也不重复堆叠 provenance。

### T3.2 有效输入接通全部消费者（A03 / R03）

按 §1.2 调用点清单逐一改造。**关键是 `mod_integration_view.R:82`**——
第二层是整合模块自己从 `proj$experiments` 取的，不改就会继续出现
"第一层清洗、第二层未清洗"。

**验收**：两层制造可识别的待过滤特征/样本，分别应用不同 QC，反复切层、保存、
重启、恢复，再跑两层差异与整合；逐层断言实际输入矩阵、ID 集、样本数、
特征数与已提交状态一致。

### T3.3 失效同时作用于四处（R06 / **V01**）

用 §1.5 的 `retire_invalid()` 替换 v3 的单遍 `invalidate_dependents()`。
调用合同显式化：`activation` 作形参传入，Core 不从外部作用域读会话变量。

`downstream_generation` 保留，但只负责**通知视图同步**，不再是失效机制本身。

**验收**（R06 + V01）：
"运行中的项目"与"刚恢复、未打开分析页面的项目"两种状态下改 QC，检查
**内存项目 / 视图 / autosave / 重启恢复 / 导出**五处；
只更新第二层时整合失效；diff→enrich→integration 多级传递失效；
打乱 bundle 顺序结论不变。

### T3.4 项目激活与严格归属（A06 / R07 / **V02**）

按 §1.4 拆分为持久化来源与会话令牌，按 §1.7 的"打开/恢复"与"另存副本"实现。

- 归属检查**不允许 `%||%` 填空**（R07）：无 `provenance_binding` 的结果
  一律视为不属于当前激活
- 打开文件时**先验证来源与依赖链，再装载视图**；不原样使用历史 `activation`，
  也不给历史结果贴新令牌
- 新格式、旧格式迁移、另存副本三条路径分别定义（V02：v3 只覆盖了旧格式）

**验收**（R07 + V02）：
正常保存的结果重开后仍能查看和导出（v3 会全判无效）；
同时旧 `activation` 的晚到回调被拒；同 `project_id` 重开旧存档确实生效；
另存副本、无归属结果、旧项目迁移各有用例。

---

## 5. 批次 4 —— 结果完整性

### T4.1 整合保存两层效应（A07 / R10）

`integration-concordance.R` 的 `out` 增加
`effect_a/effect_b/p_a/p_b/effect_type_a/effect_type_b/analysis_id_a/analysis_id_b`。
**来源要写清楚**（v3 说"值在 joined 里已有"不准确）：`effect_*`/`p_*` 来自 `joined`；
`effect_type_*` 与 `analysis_id_*` 来自对应的 `diff_bundles[[tag]]`。

删除 `plot-integration.R:232` 的 `integration_fill_effects()`；缺列时显示
"无法绘制"，不由差值反推。

**测试修正**（R10）：默认连接键是 `feature_symbol`，v3 的测试用主层
`feature_id` 的 `match()` 索引去取第二层——两层行序与 feature_id 都可能不同。
按真实连接键分别匹配；打乱两层行序、使用不同 feature_id、覆盖重复 symbol 的既定规则；
补空结果模板兼容。

**`is_significant` 拆列贯通消费者**：同步 schema、空结果模板、图表、过滤、报告、
旧 `.omp` 读取；明确 `p_a`/`p_b` 存 raw 还是 adjusted，不靠名字猜。

**统计口径见 D3**：在有经验证的合并方案之前，只出描述性一致性与分层显著性，
关闭未经验证的合并显著性结论。

### T4.2 模型设置：拒绝 + 裁剪 + 第二层继承（A08 / R01 / v3 复审 §3）

#### 拒绝与裁剪都要做

**措辞更正**：基线 `prune_backend_args()` 是**静默裁剪**，并不承担拒绝职责；
拒绝是本次新增的。v3 写反了。

```r
# 两件事都要做,不能用前者替代后者：
#   1. 用户显式传了后端不支持的设置 → 报错（本次新增,修 A08）
#   2. 后端签名里没有的形参键 → 移除（基线的必要职责,漏掉就 unused argument）
prune_backend_args <- function(method, analysis_type, args) {
  drop <- backend_unsupported_args(method, analysis_type)   # 按 method × analysis_type
  supplied <- names(args)[!vapply(args, is.null, logical(1))]
  hit <- intersect(drop, supplied)
  if (length(hit)) stop(sprintf(
    "method = '%s' with analysis_type = '%s' does not support %s. %s",
    method, analysis_type, paste(hit, collapse = ", "),
    suggest_alternative_backend(method, hit)), call. = FALSE)
  args[setdiff(names(args), drop)]        # ← 基线的裁剪,必须保留
}
```

**测试四项齐全**（R01）：不传成功 / 显式 NULL 成功 / 实际值拒绝 /
支持的后端确实应用（断言结果与未校正不同）。

#### 第二层模型继承（v3 复审 §3 —— v3 压缩掉了正常路径）

只修 `prune_backend_args()` **不会**自动修好第二层模型。
`mod_integration_view.R:224` 必须：

1. 从主分析继承 `covariates` / `paired_col`
2. 检查第二层 `meta_df` 是否有这些列、有效数据是否支持该模型
3. **不支持时停止**，说明缺什么（如"第二层缺少 age 列"）
4. 用户可显式选择"接受未校正的第二层"——该例外**持久记录**，
   并在页面与导出中**持续可见**，不得仍声称"两层模型一致"

UI 隐藏协变量控件后，**服务端仍要处理遗留 input**（Shiny 保留隐藏控件的旧值）。

### T4.3 富集：ORA 与 GSEA 都保存完整通路表（A09 / **V06**）

#### v3 的遗漏

v3 只给 ORA 补了 `output_p_cutoff = 1`。GSEA 的"不要预筛特征"是另一件事——
`mod_enrich_view.R:110` 的 GSEA 调用**不传任何阈值**，
`run_enrichment()` 默认 `p_cutoff = 0.05`，
`enrich-gsea.R:81` 的 `output_p_cutoff %||% p_cutoff` 于是把输出截断在 0.05。

```r
# 两个分支都要。输入完整排名列表和输出完整通路表是两个要求,不能相互替代。
if (identical(type, "ora")) {
  omicsCore::run_enrichment(bundle, type = type, database = db_arg, direction = dir_,
    p_cutoff = thr$p_cutoff,          # 哪些基因进 ORA
    output_p_cutoff = 1,              # 保存完整通路表
    p_preference = thr$p_preference, effect_cutoff = thr$effect_cutoff)
} else {
  omicsCore::run_enrichment(bundle, type = type, database = db_arg, direction = dir_,
    output_p_cutoff = 1)              # ← v3 漏了；GSEA 仍不预筛特征
}
```

三个阈值分别记入 params：`feature_p_cutoff` / `output_p_cutoff` / `display_p_cutoff`。
**Core 的 `run_enrichment()` 也须记录 `output_p_cutoff`**（任务范围不只 App）。
放宽显示不得改变 ORA 输入基因集；GSEA 不得按特征显著性预筛。

**测试**：ORA 与 GSEA **各有真实计算路径的用例**。
v3 对 `session$output$hits_table` 直接 `nrow()` 是错的——那是渲染后的 widget；
改为服务端检查过滤后的数据，浏览器用例检查可见行与阈值交互。

### T4.4 provenance 与脚本重放（A05 / R04 / V03）

provenance 唯一权威在 `omics_input$provenance`（§1.3），导出经
`effective_input()$provenance` 读取——与被 hash 的记录是同一份。

其余修正保留：`normalize_omics(input, method, offset)` **无 `...`**，记录 `offset`；
各步骤返回值取法不同（`read_omics()` 取 `$input`，`run_qc()` 取
`$results$cleaned_input`），不能无差别串成同一种赋值；
补齐归档文件重定位与摘要验证、`sheet_roles`、实际默认参数、随机种子与实现版本、
两个差异 bundle、sample link、基因集版本。

**旧项目缺记录时不得声称"复现分析但不复现导入"**——已知遗漏归一化会改变数值。
应禁用完整复现导出，或生成明确停止执行、要求补充信息的草稿。

**验收**：独立目录从归档原文件执行脚本，**先按稳定键对齐**再比较中间矩阵、
特征/样本 ID、参数、effect、statistic、P 值、校正 P 值与整合结果。

---

## 6. 批次 5 —— 失败可见与交付诚实

| ID | 任务 | 要求 |
|---|---|---|
| T5.1 | autosave 状态（A10） | 返回值 logical → list，**须审查所有调用者与测试**。状态区分：未保存修改 / 保存中 / 成功 / 失败 / 最后成功时间。`message()` 不代表容器删除后日志仍在，持久日志证据由 T6.2 取得 |
| T5.2 | 备份（A12） | pipefail + 临时文件 + 校验后原子替换 + 独立恢复演练。**> 1KB 只是粗筛**；故障注入必须验证**旧轮换文件摘要不变** |
| T5.3 | PDF（A13） | 若保留，装字体之外须实际配置使用 CJK 字体，覆盖正文/表格/图形/失败时临时文件清理；若不提供，**所有入口同步隐藏或降级** |
| T5.4 | `launch(project=)`（A14） | 与项目激活/恢复走**同一入口**（§1.7 的"打开"），不新增状态旁路 |
| T5.5 | 示例隔离 **+ 报告范围决策**（A15） | ①真实项目的"未计算/失效/失败"**全部**为空状态，不因 `is_demo(TRUE)` 显示示例，**不止 ActivePathways 一个面板**；②**报告范围必须二选一并同步所有声明**：补齐参数值清单/完整结果表/主要图形/模型说明，**或**明确限定为"分析概要"。保留"报告"这个标题不等于覆盖了这个决策 |
| T5.6 | async（A17） | ①注释与产品说明修正为"发起会话仍阻塞"；②`callback_admissible()` 校验成功与失败两类回调；③**队列与并发上限**、**会话结束后的回调丢弃规则**同为实现与验收项——不能只增加版本字段 |

---

## 7. 批次 6 —— 放行

- **CI**（A16）：加 `deploy/**` 触发；**构建并测试最终发布 digest**；
  CVE 扫描要有明确处置规则；修复后先跑完整包回归，再在正式镜像验证关键路径
- **部署**（§5）：容器东西向 A/B 负面测试、宿主机 loopback 路径、
  `.dockerignore` 补 `deploy/keycloak`/`deploy/cron`/**实际生成的 `host.env` 与 `.env`**
  （构建上下文包含 ≠ 镜像内包含，分别验证）、`render.sh` 密钥来源、
  文档步骤顺序、nginx 规范 hostname、基因集缓存权限与 seed 验证、
  软配额说明、多会话 autosave 覆盖策略、24h 回收 / 登出 / OOM 实测
- **验收单**：每项有测试 ID、对应问题、fixture、执行环境、预期数值/状态、
  结果日志、发布 digest

**放行标准**沿用审查报告 §8 第三批 9 项，并加两条：
对确定的代码缺陷，先验证测试在基线上暴露问题、再验证修复后通过；
对新增能力或关闭的功能，用适用的正反向契约测试，
不机械要求运维流程套用"先红后绿"。

场景矩阵（v3 复审 §3 + V01–V06 新增）：

| 场景 | 关键断言 |
|---|---|
| 恶意与正常列名 × 五后端 × 多 analysis_type | 无副作用；合法输入有有效结果；contrast/P 值与参考模型吻合 |
| 不支持参数：缺省 / NULL / 实际值 | 前两者成功，第三者拒绝；支持的后端确实应用 |
| 两层分别提交 QC，往返切层 | 两层已应用数据都保留；主/次层差异与整合用正确版本 |
| 预览失败后应用、切层后失败 | 拒绝，不给旧预览重新绑定身份 |
| **应用后放宽阈值** | 恢复应保留的特征（证明从基线重算） |
| 已保存项目在已有会话中恢复 | 三个 QC 控件 + 未暴露参数原值、已应用数据、结果归属一致 |
| 恢复后改变 QC | 项目/视图/autosave/重启/报告/脚本一致失效 |
| **同输入连续两个请求** | 旧回调（成功与失败）均被拒收 |
| **只更新第二层** | 整合失效；diff→enrich→integration 多级传递 |
| **打乱 bundle 顺序** | 有效集合结论不变 |
| 同数据不同项目；同项目重开旧存档 | 无串写；**保存的有效结果重开后仍可查看导出** |
| 另存副本 | 新身份、来源迁移、旧会话任务被拒 |
| 四条 raw scale 路径 | 各有明确结果（见 T1.4） |
| **ORA 与 GSEA 各自** | 保存完整通路表；显示阈值可逆 |
| 导出脚本独立重放 | 中间矩阵与关键数值吻合；缺记录时明确停止或降级 |
| 正式镜像与目标服务器 | 多账号隔离、并发资源、备份恢复、日志留存、升级回滚 |

---

## 8. 决策

| ID | 问题 | 状态 |
|---|---|---|
| D1 | v1 是否放弃 PDF 导出？ | 建议放弃，只承诺 HTML |
| D2 | 网页 RDS | **已定**：固定拒绝，不给开关；`OMICSAPP_ALLOW_RDS` 取消 |
| D3 | concordance 合并检验口径 | 待统计口径确认；在此之前只出描述性一致性 + 分层显著性 |
| D4 | 用户文件误删是否需可恢复？ | 待产品决策 |
| D5 | §1 数据合同 | 本版按 V01–V06 补齐，**待评审签字** |
| D6 | `outlier_sd_threshold` 控件 | **已定**：不新增；新项目用默认 3，**恢复项目保留已保存的非默认值** |
| D7 | T5.5 报告范围 | **新增待决**：补齐为完整报告，还是限定为"分析概要"并同步 README |

**记入待办、不在本次范围**：完整的 method × analysis_type × assay scale
效果量合同。但 T1.4 的模型入口拒绝**不延期**。

估工：待 §1 评审通过后按任务估算。批次 1 四项可先单独估。

---

## 9. 本版新增的复现证据

对 V03/V05/V06 的源码确认（V01/V02/V04 是对 §1 自身的核对）：

```text
V05  mod_import_view.R:348  checkboxInput(ns("normalize"), "Normalize on import", TRUE)
     mod_import_view.R:536  未勾选 → pending_normalize_method() 返回 "none"
     mod_import_view.R:563  do_normalize = FALSE → 整个归一化分支被跳过
     ⇒ v3 新增的失败处理不会执行,原始强度照常提交

V06  mod_enrich_view.R:110  GSEA 分支只传 bundle/type/database/direction
     enrich-gsea.R:81       p_cutoff = output_p_cutoff %||% p_cutoff
     ⇒ 默认 0.05 成为输出截断,与 ORA 是同一个缺陷的两个分支

V03  input_fingerprint(path, omics_type, assay_type, normalize)   ← 无 sheet_roles
     mod_import_view.R:612  fingerprints_match() 早退
     ⇒ 改角色重导入同一文件被判"nothing to re-import",到不了新版本机制
```

v3 中我自己写错、本版已改的：

```text
§1.1/§1.4  activation 声明不持久化,却进入随 .omp 保存的 binding,
           binding_valid() 又要求它相等 ⇒ 任何被恢复的项目,
           它自己保存的结果都会被判无效（与 T3.3 的目标相反）
§1.4       binding 声明了 request_seq,binding_valid() 不读
§1.4       upstream_still_present() 只有名字没有定义
§1.5       invalidate_dependents(project, tag, new_version) 未用 tag/new_version,
           且从外部作用域读 activation
§1.3       input_version() 缺 feature_df / omics_type；未应用 QC 时每次重新 hash
§1.3       provenance 在 pp 与 pp$cleaned 两处并存
T3.1       run_qc(qc_input, !!!params) 是 rlang 注入语法,普通调用无效
T4.2       把基线 prune_backend_args() 说成"承担拒绝职责",基线实为静默裁剪
T4.1       称新增字段值"都已在 joined 中",实际 effect_type/analysis_id 另有来源
§0.2       称 §8 给出初估,§8 实际没有
```

A01–A12 的原始复现、QC 参数链路实测、"覆盖法 vs 双槽位"实测、
R01/R02 的补丁缺陷复现，均保持有效，见上游三份文档的证据索引。
