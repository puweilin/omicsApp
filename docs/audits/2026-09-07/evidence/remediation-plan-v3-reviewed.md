# omicsApp 整改实施方案（v3）

源码基线：`f756751`
上游文档：[审查报告](README.md)（17 项，已全部复核属实）、[方案复审](REMEDIATION-REVIEW.md)（R01–R10）
本版变更：按复审 R01–R10 重写。**R01/R02 的补丁缺陷已独立复现**，见 §9。

---

## 0. 对复审意见的处置

R01–R10 全部接受。我复现了 R01、R02、R03、R04、R08、R10 的关键断言（§9），
其余由源码判定成立。v2 方案的两个补丁会引入新回归，两处数据流互不闭合，
**v2 不可执行，本版替换它。**

| 编号 | 处置 | 落到本版何处 |
|---|---|---|
| R01 | 接受。我的 `prune_backend_args()` 删掉了原函数必需的裁剪职责 | T4.2 重写 |
| R02 | 接受。前缀剥离 + 按名字建 contrast 会指向错误系数 | T1.1 重写 |
| R03 | 接受。单一 `bundles$qc` 槽位撑不起多层；整合第二层绕过 QC | **T2.0 新增数据合同** |
| R04 | 接受。写在 `cleaned_input`、读在 `experiments` 两处不闭合 | T2.0 + T4.4 |
| R05 | 接受。应用按钮会给失败后残留的旧预览贴当前层身份 | T3.1 重写 |
| R06 | 接受。清视图不等于清项目里已保存的结果 | T3.3 重写 |
| R07 | 接受。`project_id` 认不出"打开同一项目的旧存档" | T3.4 重写 |
| R08 | 接受。`allow_raw` 扩了 API 却没有完整合同 | T1.4 缩范围 |
| R09 | 接受。批次里程碑措辞与放行标准冲突 | §0.2 + §7 |
| R10 | 接受。7 处示例函数名/测试接口错误 | 各任务内逐条改 |

### 0.1 v2 的三个根本错误

1. **把"能编译"当成"能执行"。** T1.1、T3.2 的补丁没有在真实调用路径上跑过。
   T3.2 会让每一次正常的 t-test / lm 调用失败。
2. **把"两个槽位够用"当成"设计完备"。** 上一轮讨论里我确认了
   `experiments` + `bundles$qc` 两个槽位都持久化，就据此宣称"不需要新结构"。
   这对单层成立，对多层不成立——`bundles$qc` 只有一个槽。
3. **把接口一致性交给了叙述而不是定义。** provenance 写在 `cleaned_input`、
   导出读 `experiments`，两个任务各自看都对，合起来不闭合。

本版的对策是 §1 的数据合同——**先定读写契约，再分派实现**。

### 0.2 里程碑口径（R09）

| 节点 | 含义 | 不含义 |
|---|---|---|
| 批次 1 完成 | 安全入口整改完成，可进入后续验证 | **不是**可发布，**不是**科研生产就绪 |
| 批次 2–5 完成 | 功能修复完成，可进入镜像验收 | 仍未取得目标服务器证据 |
| 批次 6 完成 | 取得放行证据 | —— |

**任何时点都不得据本方案放行部署。** 放行需要 §7 的全部证据。

估工：**重新估算。** 在 §1 合同评审通过前，任何总工时数字都不可信；
v2 的"24–26 人日"作废。§8 给出的是按合同拆分后的初估，仍需评审确认。

---

## 1. 数据合同（先定义，后实现）

R03/R04/R06/R07 是同一个缺陷的四个面：**没有一个明确的"哪份数据、哪一版、
属于谁"的契约**。四个任务各自打补丁必然互不闭合。

### 1.1 状态树

```text
project
├─ project_id            持久化，随 .omp 走，另存副本时重新生成
├─ activation            每次 open / restore / import-first 递增（不持久化）
├─ experiments[[tag]]    导入层（可能已归一化），不可变
├─ preprocessing[[tag]]  ← 新增：按层的已应用预处理
│    ├─ input_version    该层有效输入的不可变版本 ID
│    ├─ cleaned          omics_input，下游必须使用的有效输入
│    ├─ params           解析后的完整 QC 参数（含经 ... 传入的）
│    └─ provenance       read_omics → normalize_omics → run_qc 的可重放步骤
└─ bundles[[name]]       分析结果
     └─ input_info$binding = list(project_id, activation, tag, input_version,
                                  upstream_analysis_ids, request_seq)
```

`preprocessing[[tag]]` 与 `experiments[[tag]]` 同为按层索引，两者一一对应。
QC 未应用时 `preprocessing[[tag]]` 不存在，有效输入即 `experiments[[tag]]`。

### 1.2 唯一入口

所有需要"这一层当前该用哪份数据"的地方，**必须**经过同一个函数——
不允许任何模块再直接读 `proj$experiments[[tag]]` 做分析：

```r
# 有效输入 + 它的版本。差异分析、整合的第二层差异、报告、导出脚本
# 全部经此获取,不再各自决定读哪个槽位。
effective_input <- function(project, tag) {
  pp <- project$preprocessing[[tag]]
  if (is.null(pp)) {
    raw <- project$experiments[[tag]]
    return(list(input = raw, version = input_version(raw),
                source = "imported", provenance = raw$provenance))
  }
  list(input = pp$cleaned, version = pp$input_version,
       source = "qc_cleaned", provenance = pp$provenance)
}
```

**调用点清单**（实现时逐一改造，缺一个就是 R03 复发）：

| 位置 | 现状 | 改后 |
|---|---|---|
| `mod_diff_view.R:75` `active()` | `exps[[tag]]` | `effective_input(proj, tag)` |
| `mod_integration_view.R:82` `sec <- exps[[secondary_idx]]` | 导入层 | `effective_input(proj, secondary_tag)` |
| `export-script.R:216` 导出循环 | `experiments[[tag]]$provenance` | `effective_input()$provenance` |
| `default-report.Rmd` | 项目概要 | 同上 |

### 1.3 输入版本

R06 指出：无指纹时用矩阵形状当身份，两个不同的 20×10 会同身份。

```r
# 不可变版本 ID。上传指纹只能标识"文件",标识不了"文件经过哪些处理之后"。
# 形状不能当身份 —— 两个不同的 20x10 矩阵是两份不同的数据。
input_version <- function(input) {
  digest::digest(list(
    matrix      = input$expr_mat,          # 内容,不是形状
    samples     = colnames(input$expr_mat),
    features    = rownames(input$expr_mat),
    meta        = input$meta_df,
    assay_type  = input$assay_type,
    provenance  = input$provenance
  ), algo = "xxhash64")
}
```

大矩阵的开销：63k × 258 的 `digest()` 约 0.2–0.4 s。只在提交 QC、
导入、打开项目时计算并缓存进 `preprocessing[[tag]]$input_version`，
不在反应式里反复算。**实现前需实测确认这个开销可接受**（性能预算见 T2.0 验收）。

### 1.4 结果归属与失效

```r
# 结果是否仍然有效：五个条件全中才算。
# 只比对 analysis_id + fingerprint 不够（复审 R06/T5.4）。
binding_valid <- function(binding, project, activation) {
  identical(binding$project_id, project$metadata$project_id) &&
  identical(binding$activation, activation) &&
  identical(binding$input_version, effective_input(project, binding$tag)$version) &&
  all(vapply(binding$upstream_analysis_ids, upstream_still_present, logical(1), project))
}
```

**失效必须同时作用于三处**（R06）：内存项目的 `bundles`、视图的 `reactiveVal`、
autosave 写出的快照。只清视图不够——恢复的项目可能带着 saved diff 而视图为空。

### 1.5 合同评审出口

§1 需要**先单独评审通过**，再分派 §2 之后的实现。评审要回答：

1. `preprocessing` 进 `.omp` 的兼容策略（旧文件无此字段）
2. `input_version` 的算法与性能预算
3. `activation` 的语义边界：open / restore / 首次 import / 另存副本 / 关闭到 demo
4. 失效是"删除"还是"移出有效集合但保留为历史"
5. `is_significant` 拆列后旧 `.omp` 的读取语义（复审补充 2）

---

## 2. 批次 1 —— 独立入口修复

四个任务互不依赖，不依赖 §1 合同，可立即并行开工。

### T1.1 设计矩阵与 contrast（A01 / 复审 R02）

#### v2 补丁为什么错

v2 用 `v1, v2, …` 做内部标识，再 `startsWith(cn, "v1")` 剥前缀还原用户名字。
两个缺陷已复现（§9）：

```text
组名恰为 "v2" + 一个协变量  →  设计列名 = A | v2 | v2
                              实际 contrast = -1 | 0 | 1   ← 指向协变量
                              期望 contrast = -1 | 1 | 0
九个协变量                  →  v10 被 startsWith("v1") 吞成 "0" → makeContrasts 报错
```

第一种**不报错**，出完整结果表和正常火山图，但比较的是错误系数。
这正是整份审查在讲的那类错误，v2 的"修复"制造了一个新的。

#### 正确做法：三样东西分开管理

| 层面 | 取值来源 | 是否接触用户字符串 |
|---|---|---|
| 公式项 | `t001, t002, …` 程序生成 | 否 |
| 因子水平 | `L001, L002, …` 程序生成 | 否 |
| 设计矩阵列名 | `t001L002` 等，程序生成且唯一 | 否 |
| 用户显示标签 | 仅在标准化输出时映射回去 | 是，但只进显示 |

**关键：contrast 由系数位置构造成数值向量，全程不按名字匹配。**
`limma::makeContrasts()` 解析字符串，直接弃用。

```r
# packages/omicsCore/R/diff-design.R（新增）
#
# 用户上传的列名和因子水平都可能与任何内部约定相撞 —— 组名可以恰好
# 叫 "v2",也可以恰好叫 "t001"。所以设计矩阵的列名全部由序号生成,
# 用户字符串一次都不进入 model.matrix 看得到的地方,contrast 按列位置构造。
build_design <- function(meta, terms, primary = NULL, intercept = TRUE) {
  missing <- setdiff(terms, colnames(meta))
  if (length(missing)) {
    stop("Columns not found in metadata: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ids <- sprintf("t%03d", seq_along(terms))
  d <- meta[, terms, drop = FALSE]
  names(d) <- ids

  level_map <- list()
  for (i in seq_along(ids)) {
    if (is.factor(d[[i]]) || is.character(d[[i]])) {
      f <- factor(d[[i]])
      level_map[[ids[i]]] <- levels(f)                 # 用户水平,只留作映射
      levels(f) <- sprintf("L%03d", seq_along(levels(f)))
      d[[i]] <- f
    }
  }

  rhs <- paste(c(if (!intercept) "0", ids), collapse = " + ")
  mm <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = d)
  if (anyDuplicated(colnames(mm))) {
    stop("Internal design column names collided; this is a bug.", call. = FALSE)
  }
  structure(mm,
    term_ids   = stats::setNames(as.list(terms), ids),
    level_map  = level_map,
    primary_id = if (is.null(primary)) NULL else ids[match(primary, terms)])
}

# 数值 contrast。makeContrasts() 解析字符串,不用。
design_contrast <- function(mm, control_level, case_level) {
  pid <- attr(mm, "primary_id")
  lv  <- attr(mm, "level_map")[[pid]]
  i_c <- match(control_level, lv); i_k <- match(case_level, lv)
  if (anyNA(c(i_c, i_k))) {
    stop("Contrast level not present in the group column.", call. = FALSE)
  }
  k <- match(paste0(pid, sprintf("L%03d", c(i_c, i_k))), colnames(mm))
  if (anyNA(k)) stop("Contrast column missing from the design.", call. = FALSE)
  v <- numeric(ncol(mm)); v[k[1L]] <- -1; v[k[2L]] <- 1
  matrix(v, ncol = 1L,
         dimnames = list(colnames(mm), paste(case_level, "-", control_level)))
}

# 用户标签只在这里出现,只进显示
design_labels <- function(mm) { ... }
```

`diff-limma.R:66-90` 改为：

```r
design <- build_design(target_meta, c(group_col, covariates),
                       primary = group_col, intercept = FALSE)
fit  <- limma::lmFit(expr_sub, design)                       # 或带 block 的分支
fit2 <- limma::eBayes(limma::contrasts.fit(
          fit, design_contrast(design, control_group, case_group)))
```

#### 逐分支处理（不用"同理"带过，复审 R02 特别点名）

| 位置 | 处理 | 需明确 |
|---|---|---|
| `diff-limma.R:66,74` 两组 | 上述 | —— |
| `diff-limma.R:308,316` ANOVA | `build_design` + 数值 contrast 矩阵（多列） | 检验对象、自由度 |
| `diff-limma.R:205-234` 连续/样条 | 响应固定 `y`/`cont`；样条项 `splines::ns(cont, df)` 由我们构造 | 矩阵秩、缺失值、df |
| `diff-lm.R:52,169,219,223` | 同上；`lm()` 周边逻辑（`summary`、系数提取）需逐项对照 | **不可直接换 `lm.fit()`** |
| `diff-edger.R:84,85` | `build_design` + `glmQLFTest` 的 coef 位置 | coef 索引来源 |
| `diff-deseq2.R:174,260` | DESeq2 用 `design =` 公式；改传 `full`/`reduced` 设计矩阵 | `results(contrast=)` 的列表形式 |
| `run-diff.R:399` | 删掉反引号方案，走 `build_design` | —— |

`diff-limma.R:177` 的 `model.matrix(~ ., data = design_df)` 无用户字符串，不改。

#### 测试

新增 `packages/omicsCore/tests/testthat/test-design-contrast.R`：

```r
# 一、安全：必须证明真的走到了目标建模路径,不能只 try() 后看 marker（复审 R10）
test_that("a payload column name reaches the model path and has no effect", {
  for (m in SUPPORTED_DIFF_METHODS) {
    marker <- withr::local_tempfile()
    inp <- fixture_for_method(m, group_col_name = sprintf('file.create("%s")', marker))
    res <- run_diff(inp, method = m, analysis_type = "group",
                    group_col = attr(inp, "payload_name"),
                    control_group = "A", case_group = "B")
    expect_s3_class(res, "analysis_bundle")          # 到达了建模路径并成功
    expect_false(file.exists(marker))                # 且没有副作用
  }
})

# 二、正确性：与手工构造的安全设计矩阵逐位比对（复审 R02 要求）
test_that("group levels that collide with internal ids still contrast correctly", {
  meta <- data.frame(grp = rep(c("A", "v2", "t001", "L001"), each = 4),
                     age = rnorm(16))
  for (case in c("v2", "t001", "L001")) {
    got <- run_diff(inp, method = "limma", analysis_type = "group",
                    group_col = "grp", control_group = "A", case_group = case,
                    covariates = "age")
    ref <- reference_limma_fit(meta, case)           # 手工安全矩阵 + 数值 contrast
    expect_equal(got$results$diff_result_df$effect,  ref$effect,  tolerance = 1e-10)
    expect_equal(got$results$diff_result_df$p_value, ref$p_value, tolerance = 1e-10)
  }
})

test_that("ten or more terms do not truncate any design column", { ... })  # v10 → "0"

# 三、合法的真实列名与水平
test_that("spaces, CJK and punctuation model correctly", {
  for (col in c("年龄 (岁)", "batch-1", "treatment group", "a.b", "a|b")) { ... }
  for (lv in c("对照组", "G 1", "G-2")) { ... }                    # 水平也来自上传
})
```

`test-backend-agreement.R`、`test-degenerate-designs.R`、`test-metamorphic-diff.R`
全跑，数值零变化。

#### 验收

- 5 个后端 × 3 类 payload：到达建模路径且**成功返回**，marker 不存在
- 组名为 `v2` / `t001` / `L001` 时，effect 与 p_value 和手工安全模型
  在 `1e-10` 内一致
- ≥10 个模型项时无列名截断
- 现有数值回归零变化

---

### T1.2 网页端关闭 RDS 上传（A02 / 复审 R10）

**现状**：`read-omics.R:299` `readRDS()` 后靠 `inherits()` 判断，再访问
`obj$expr_mat`。实测活动绑定在解析阶段触发。

**文案修正**（复审要求）：不把所有 `readRDS()` 一概说成"立即执行任意代码"。
准确表述是——**加载并使用不可信 R 对象可以执行代码**：触发点是字段访问、
`print`/`format` 等方法派发，而"反序列化后检查 class"不构成边界。

**改法**：`read_omics()` 增加 `allow_r_objects`，rds 分支加门；
`mod_import_view.R` 固定传 `allow_rds()`，默认由
`OMICSAPP_ALLOW_RDS` 控制且默认 `false`；**多用户部署固定关闭**（不给开关）。
`fileInput` 的 `accept` 去掉 `.rds` 只是前端提示，服务端拒绝才是防线。

**测试**（`test-hostile-inputs.R`）：必须在**真实上传解析入口**证明拒绝
发生在读取对象之前——不是只测 `read_omics()` 直调。

```r
test_that("the upload path refuses .rds before the object is read", {
  # 走 mod_import_view 的解析入口,不是直接调 read_omics()
  session$setInputs(file = list(datapath = evil_rds, name = "x.rds"))
  expect_match(session$output$notices, "RDS input is disabled")
  expect_false(file.exists(marker))
})
```

---

### T1.3 项目路径服务端校验（A11 / 复审 R10）

**现状**：`mod_project_view.R:329,336` 把 `input$saved_pick` 直接交给
`store_load_project()` / `store_delete_project()`。实测 `../outside/victim`
读、删均成功。

**函数名修正**：实际入口是 **`list_saved_projects()`**（`project_store.R:282`），
v2 写的 `store_list_projects()` 不存在。

```r
resolve_project_path <- function(slug, dir = omicsapp_data_dir()) {
  if (!is.character(slug) || length(slug) != 1L || is.na(slug) || !nzchar(slug)) {
    return(list(ok = FALSE, message = "No project selected."))
  }
  if (!slug %in% list_saved_projects(dir)$slug) {          # 服务端白名单
    return(list(ok = FALSE, message = "Unknown project."))
  }
  path <- tryCatch(normalizePath(project_path(slug, dir), mustWork = TRUE),
                   error = function(e) NA_character_)      # 校验后文件消失
  if (is.na(path)) return(list(ok = FALSE, message = "Project file is no longer available."))
  root <- normalizePath(dir, mustWork = TRUE)
  if (!startsWith(path, paste0(root, .Platform$file.sep))) {
    return(list(ok = FALSE, message = "Unknown project."))
  }
  list(ok = TRUE, path = path)
}
```

越界与不存在统一回 `"Unknown project."`，不泄露存储布局；
文件消失单列一个结构化错误（复审要求）。

**测试**（`test-store-adversarial.R`）：现有用例用 `../../etc/passwd`，
实际访问 `passwd.omp`（不存在），只证明了"不崩溃"。替换为在 store 外
**真实创建** `.omp`，断言读删后仍存在；补绝对路径、软链接、
**正常中文项目名**、校验后文件消失四个变体。

---

### T1.4 归一化失败中止提交（A04 / 复审 R08）

#### 缩范围

复审 R08 指出 v2 的 `allow_raw` + `assert_log_scale()` 有三个问题，
我核对属实：

| v2 的写法 | 实际情况 |
|---|---|
| `assert_log_scale()` 只挡 `raw_intensity` | `LOG_SCALE_ASSAY_TYPES` 只有 5 个；`tpm`/`fpkm`/`raw_count` 都是线性，一个没挡。函数名与错误文案名不副实 |
| 一律改标 `mean_difference` | 现有词表是 `mean_diff`（t-test）、`beta`（lm）、`correlation`（limma 连续）、`F_statistic`（ANOVA）。新名与 `mean_diff` 重复，且连续/ANOVA 语义不是均值差 |
| `allow_raw` 贯通分派与后端 | limma 后端不接受该参数，标准化函数也不收 `assay_type`；这是一整套合同，不是一个开关 |

**本任务只做导入侧**，不引入 `allow_raw`：

```r
# mod_import_view.R:568
if (inherits(normalized, "error")) {
  import_error(sprintf(
    "Normalization (%s) failed: %s. Nothing was imported — change the method and retry.",
    method, conditionMessage(normalized)))
  return(invisible(NULL))          # 不调用 confirmed_input()
}
```

文件、角色、方法选择不清空，可直接重试。**不提供"以原始强度继续"的勾选**——
那需要 R08 说的整套效果量合同，不在本次范围。

模型侧的尺度校验作为**独立议题**记入 §8 待办，与 method × analysis_type ×
assay scale 的输入/输出合同一起设计，不在本次塞进来。

**测试**：`test-import-normalization.R` 断言提交次数为 0、项目层数不变、
错误状态可读、可重试成功。**模块返回合同要先写清楚**——v2 用的
`session$returned()$error` 与现有接口不符（复审 R10）。

---

## 3. 批次 2 —— 数据合同实现

### T2.0 按层预处理存储与有效输入解析（新增，R03/R04/R06/R07 的共同地基）

**这是 §1 合同的落地，批次 3 之后的全部任务依赖它。**

1. `omics_project` 增加 `preprocessing`（命名列表，按层索引）与
   `metadata$project_id`；`load_project()` 对旧 `.omp` 补齐并标注
   `metadata$schema_note`
2. 实现 `effective_input()`、`input_version()`、`binding_valid()`
3. `project_id` **不用 `sample()` 生成**（复审 R07：会消耗分析 RNG 流，
   且在相同 `set.seed()` 下重复）。改用
   `paste0("prj_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", basename(tempfile("")))`
   或 `uuid::UUIDgenerate()`
4. `activation` 计数器：`open` / `restore` / 首次 import / 另存副本各自递增
5. `subset_omics()` 保留 `provenance` / `source_path` / `source_fingerprint`
   （实测现在都丢，§9）

**验收**：

- 两个层各自持有独立的 `preprocessing`，互不覆盖
- `.omp` round-trip 后 `preprocessing`、`input_version`、`binding` 完好
- 旧 `.omp` 加载不报错，且 `schema_note` 可见
- **性能预算**：63k × 258 矩阵的 `input_version()` < 0.5 s，
  且不在反应式里重复计算（`test-perf-budget.R` 追加）

---

## 4. 批次 3 —— QC、失效与项目激活

依赖 T2.0。

### T3.1 QC 预览/应用分离 + 预览签名（A03 前置 / 复审 R05）

#### v2 补丁为什么错

`mod_qc_view.R:238` 在 QC 失败时**保留上一次成功的 `last_bundle`**
（只写 `last_error`）。v2 的应用按钮只检查 `b` 与当前层存在，就把**当前**层的
标签和指纹写进那个旧 `b`。可达场景：A 层预览成功 → 切到 B 层 → 新预览因阈值
过严失败 → 点应用 → 提交结果是 `B | fp-B | 数据来自 A`。

#### 改法：预览自带签名，应用只做验证

```r
    # 预览启动时就把身份钉进结果,应用阶段只验证、不补发身份。
    # 失败时保留旧图供参考是可以的,但旧预览不能成为当前请求的成功产物。
    preview <- shiny::reactiveVal(NULL)   # list(bundle, sig, ok)

    preview_signature <- function(a, params) {
      list(project_id = current_project()$metadata$project_id,
           activation = activation(),
           tag        = a$tag,
           input_version = effective_input(current_project(), a$tag)$version,
           params     = params)
    }

    shiny::observe({
      a <- active(); params <- collect_qc_params(input, a)
      sig <- preview_signature(a, params)
      b <- tryCatch(omicsCore::run_qc(qc_input, !!!params), error = function(e) e)
      if (inherits(b, "error")) {
        preview(list(bundle = NULL, sig = sig, ok = FALSE, error = conditionMessage(b)))
      } else {
        preview(list(bundle = b, sig = sig, ok = TRUE))
      }
    })

    shiny::observeEvent(input$apply_to_downstream, {
      p <- preview(); a <- active()
      if (!isTRUE(p$ok)) {
        shiny::showNotification("The current QC settings did not run. Fix them before applying.",
                                type = "error"); return()
      }
      if (!identical(p$sig, preview_signature(a, collect_qc_params(input, a)))) {
        shiny::showNotification("Settings changed since this preview ran. Re-run, then apply.",
                                type = "warning"); return()
      }
      apply_preprocessing(p$sig$tag, p$bundle, p$sig)   # 写 preprocessing[[tag]]
    })
```

#### 控件恢复（复审 R05 第二段）

v2 的三级回退 `isolate(input) → applied params → 默认` 优先级错了：
用户已选 `knn` 时打开保存为 `min` 的项目，仍显示 `knn`——正是要消除的错配。

**改为：项目/层激活时执行一次明确的恢复，之后才允许编辑草稿。**

```r
    # activation 变化 = 换了项目或换了存档。此时以已应用状态为准,
    # 覆盖当前草稿；之后用户的编辑才是草稿。
    shiny::observeEvent(list(activation(), active()$tag), {
      p <- current_project()$preprocessing[[active()$tag]]
      restore_qc_controls(session, p$params %||% qc_defaults(active()$input))
    }, priority = 25L)
```

页面**三态分列显示**：已保存参数 / 当前草稿 / 已应用状态。

**控件数量修正**：QC UI 只有 3 个输入控件——`missing_threshold`（:457）、
`outlier_method`（:465）、`impute_method`（:204）。**没有** `outlier_sd_threshold`
控件（v2 说"四个控件"是错的）。是否新增该控件需决策 D6；若新增，
必须贯通到 `run_qc()`，不能只写"同样恢复"。

**测试**：四个场景全部必须覆盖（复审 R05 验收）——
"成功后失败再应用"（拒绝）、"切层后预览失败再应用"（拒绝）、
"A 项目 knn → B 项目 min"（显示 min）、"已保存方法在当前环境不可用"
（明确提示，不静默改选）。

---

### T3.2 有效输入接通全部消费者（A03 / 复审 R03）

按 §1.2 的调用点清单逐一改造。**关键是 `mod_integration_view.R:82`**——
v2 说"整合模块不需要改，它们吃的是 diff bundle"是错的：第二层是整合模块
自己从 `proj$experiments` 取的（`sec <- exps[[secondary_idx]]`），
再在 `:224` 用它跑差异分析。不改就会继续出现
"第一层用清洗数据、第二层用导入数据"。

**验收**（复审 R03）：两个层都制造可识别的待过滤特征/样本，分别应用不同 QC，
反复切层、保存、重启、恢复，再跑两层差异与整合；
**逐层断言实际输入矩阵、ID 集、样本数、特征数与已提交状态一致**。

---

### T3.3 失效同时作用于项目、视图与持久化（复审 R06）

#### v2 补丁为什么不够

`app_server.R:180` 的回写 observer 只在某视图**此前持有非空 bundle**
且 generation 匹配时，才从项目里删除对应 bundle。恢复一个带 saved diff、
视图为空的项目，此后改 QC——视图继续返回 NULL，**项目里的 diff 不会被删**。
报告、保存、导出仍读得到旧结果。

#### 改法

失效是**项目事实状态**的变更，不是视图状态的副产物：

```r
# 提交 QC 时,直接在项目里把依赖它的结果移出有效集合。
# 不依赖 observer 顺序,也不依赖某个视图"曾经显示过什么"。
invalidate_dependents <- function(project, tag, new_version) {
  keep <- list(); retired <- list()
  for (nm in names(project$bundles)) {
    b <- project$bundles[[nm]]
    if (binding_valid(b$input_info$binding, project, activation)) keep[[nm]] <- b
    else retired[[nm]] <- b
  }
  project$bundles <- keep
  project$retired <- c(project$retired, retired)   # 保留为历史,但不是有效结果
  project
}
```

`downstream_generation` 保留（v2 的这部分方向对），但它只负责**通知视图同步**，
不再是失效的唯一机制。

`qc_identity` 改用 §1.3 的 `input_version`（内容摘要），不再用矩阵形状。

**验收**（复审 R06）：在"运行中的项目"和"刚从磁盘恢复、尚未打开分析页面的
项目"两种状态下改 QC，同时检查**内存项目、视图、autosave 快照、重启恢复、
导出脚本**五处。只断言页面 bundle 为 NULL 不够。

---

### T3.4 项目激活版本与严格归属（A06 / 复审 R07）

#### v2 补丁的两个缺口

1. 打开**同一项目的较早存档**：`project_id` 与数据都没变，v2 的 observer
   不清空视图，内存里的新结果仍可能回写覆盖刚打开的旧结果
2. `v$input_info$project_id %||% pid` 把**无归属结果直接当作当前项目的**——
   这条检查实际允许未知归属，起不到兜底作用

#### 改法

- `project_id`（持久）与 `activation`（每次 open/restore/import-first/另存递增）
  **分开管理**。切换判断用 `activation`，不用 `project_id`
- 归属检查**不允许 `%||%` 填空**：无 `binding` 的结果一律视为不属于当前激活，
  只能经明确的文件恢复流程迁移绑定
- `project_id` 生成不用 `sample()`（见 T2.0）

**验收**（复审 R07）：除 A→B 相同数据外，必须覆盖
"同一 project_id 重开旧存档"、"另存副本"、"无归属视图结果"、
"有归属的旧项目迁移"四类。

---

## 5. 批次 4 —— 结果完整性

### T4.1 整合保存两层效应（A07 / 复审 R10、补充 2）

`integration-concordance.R` 的 `out` 增加 `effect_a/effect_b/p_a/p_b/
effect_type_a/effect_type_b/analysis_id_a/analysis_id_b`（值在 `joined` 里已有）；
删除 `plot-integration.R:232` 的 `integration_fill_effects()`，缺列时显示
"无法绘制"而不是由差值反推。

**测试修正**（复审 R10）：v2 的测试用主层 `feature_id` 的 `match()` 索引去取
第二层的值——默认连接键是 `feature_symbol`，两层行序和 feature_id 都可能不同。
改为按**真实连接键**分别匹配，并打乱两层行序、使用不同 feature_id、
覆盖重复 symbol 的既定规则。补空结果模板兼容。

**`is_significant` 拆列必须贯通消费者**（复审补充 2）：改成两种标记时，
同步 schema、空结果模板、图表、过滤、报告、旧 `.omp` 读取；
明确 `p_a`/`p_b` 存的是 raw 还是 adjusted，不靠名字猜。

**统计口径见决策 D3**——在有经验证的合并方案之前，
可先提供描述性一致性与分层显著性，关闭未经验证的合并显著性结论。

### T4.2 不支持的设置改为拒绝（A08 / 复审 R01）

#### v2 补丁为什么错（已复现，§9）

原 `prune_backend_args()` 承担**两件事**：拒绝 + 裁剪。v2 只做了拒绝，
最后 `return(args)` 原样返回。而 `run-diff.R:217-225` 构造的 `backend_args`
**恒含** `covariates` 与 `paired_col` 键（值可为 NULL），
`run_ttest_group()` 无 `covariates` 形参、`run_lm_group()` 无 `paired_col` 形参，
于是每次正常调用都 `unused argument (covariates = NULL)`。

#### 改法

```r
# 两件事都要做,不能用前者替代后者：
#   1. 用户显式传了后端不支持的设置 → 报错（原来是静默丢弃 = A08）
#   2. 后端签名里没有的形参键 → 移除（原函数的必要职责,漏掉就 unused argument）
prune_backend_args <- function(method, analysis_type, args) {
  drop <- backend_unsupported_args(method, analysis_type)
  supplied <- names(args)[!vapply(args, is.null, logical(1))]
  hit <- intersect(drop, supplied)
  if (length(hit)) {
    stop(sprintf("method = '%s' with analysis_type = '%s' does not support %s. %s",
                 method, analysis_type, paste(hit, collapse = ", "),
                 suggest_alternative_backend(method, hit)), call. = FALSE)
  }
  args[setdiff(names(args), drop)]        # ← 保留原有裁剪
}

# 按 method × analysis_type 查表（复审要求区分二者），不是只看 method
backend_unsupported_args <- function(method, analysis_type) { ... }
```

**测试必须四项齐全**（复审 R01 验收）——v2 只有负向断言，
证明不了正常调用不回归：

```r
test_that("unsupported settings are refused without breaking normal calls", {
  expect_s3_class(run_diff(inp, method = "ttest", ...), "analysis_bundle")          # 不传
  expect_s3_class(run_diff(inp, method = "ttest", ..., covariates = NULL), "analysis_bundle")  # 显式 NULL
  expect_error(run_diff(inp, method = "ttest", ..., covariates = "age"), "does not support")   # 实际值
  b <- run_diff(inp, method = "limma", ..., covariates = "age")                     # 支持的后端
  expect_equal(b$params$covariates, "age")
  expect_false(isTRUE(all.equal(b$results$diff_result_df$p_value, unadjusted_p)))   # 确实生效
})
```

UI 侧隐藏协变量控件后，**服务端仍要处理遗留 input**（Shiny 保留隐藏控件的
上次取值）。整合的"有意未校正"例外必须在**页面与导出中持续可见**，
不能仍声称"两层模型一致"。

### T4.3 富集三阈值分离（A09 / 复审 R10）

`mod_enrich_view.R:100` 的 ORA 调用补 `output_p_cutoff = 1`；
**Core 的 `run_enrichment()` 也须把 `output_p_cutoff` 记入 params**
（复审：任务范围不应只标 App）。

特征筛选与通路显示各自明确 P 值口径；**放宽显示不得改变 ORA 输入基因集**，
**GSEA 不得按特征显著性预筛**。

**测试修正**：v2 对 `session$output$hits_table` 直接 `nrow()`——那是渲染后的
widget，不是 data.frame。改为在服务端检查过滤后的数据，浏览器用例检查可见行
与阈值交互。

### T4.4 provenance 与脚本重放（A05 / 复审 R04）

#### v2 的不闭合

写在 `cleaned_input$provenance`，导出读 `experiments[[tag]]$provenance`。
T2 又明确不把清洗数据写回 `experiments`，所以导出侧永远看不到 `run_qc` 步骤。

v2 用来论证"provenance 权威"的理由——"`bundles$qc` 会被清掉但数据还在"——
**在 v2 自己的存储约定下不成立**：清洗数据本来就只存在于那个 bundle 里。

T2.0 引入 `preprocessing[[tag]]` 之后这个矛盾消失：
数据、参数、provenance 三者同处一个按层槽位，导出经 `effective_input()` 读取。

#### 其余修正

- `normalize_omics(input, method, offset)` **没有 `...`**（v2 写错）。
  应明确记录 `offset`
- 步骤返回值取法各不相同，不能无差别串成同一种赋值：
  `read_omics()` 取 `$input`，`run_qc()` 取 `$results$cleaned_input`
- 补齐：归档文件重定位与摘要验证、`sheet_roles`、实际默认参数、
  随机种子与相关实现版本、两个差异 bundle、sample link、基因集版本
- **旧项目缺记录时不得声称"复现分析但不复现导入"**——已知遗漏归一化会改变
  数值。应禁用完整复现导出，或生成明确停止执行、要求补充信息的草稿

#### 验收（复审 R04）

在独立目录从归档原文件执行导出脚本，**先按稳定键对齐**再比较
中间矩阵、特征/样本 ID、参数、最终 effect、statistic、P 值、校正 P 值
与整合结果。不只比 effect 向量。无 provenance 的旧项目要有明确的降级测试。

---

## 6. 批次 5 —— 失败可见与交付诚实

| ID | 任务 | 复审要求的修订 |
|---|---|---|
| T5.1 | autosave 状态（A10） | 返回值 logical → list，**须审查所有调用者与测试**。状态区分：未保存修改 / 保存中 / 成功 / 失败 / 最后成功时间。`message()` 不代表容器删除后日志仍在，持久日志证据由 T6.3 取得 |
| T5.2 | 备份（A12） | pipefail + 临时文件 + 校验后原子替换 + 独立恢复演练保留。**> 1KB 只是粗筛**，不能证明可恢复；故障注入必须验证**旧轮换文件摘要不变** |
| T5.3 | PDF（A13） | 若保留 PDF，装字体之外还须实际配置使用 CJK 字体，覆盖正文/表格/图形/失败时临时文件清理；若不提供，**所有入口同步隐藏或降级** |
| T5.4 | `launch(project=)`（A14） | 与项目激活/恢复走**同一入口**，不新增状态旁路 |
| T5.5 | 示例与报告（A15） | 真实项目的"未计算/失效/失败"**全部**应为空状态，不能因 `is_demo(TRUE)` 显示示例；**不能只修 ActivePathways 一个面板** |
| T5.6 | async 版本校验（A17） | 校验须含项目激活版本、输入/QC 版本、请求序号、相关分析参数；**成功与失败回调都要检查**。仅 `analysis_id + layer_fingerprint` 不足 |

---

## 7. 批次 6 —— 放行

### T6.1–T6.3

- **CI**（A16）：加 `deploy/**` 触发；**构建并测试最终发布 digest**；
  CVE 扫描要有明确处置规则，不能只生成报告。修复后先跑完整包回归，
  再在正式镜像验证关键路径
- **部署**（§5）：容器东西向隔离 A/B 负面测试、宿主机账号 loopback 路径、
  `.dockerignore` 补 `deploy/keycloak`/`deploy/cron`/**实际生成的 `host.env` 与 `.env`**
  （复审补充 3：构建上下文包含 ≠ 镜像内包含，两者分别验证）、
  `render.sh` 密钥来源、部署文档步骤顺序、nginx 规范 hostname、
  基因集缓存挂载权限与 seed 验证、软配额说明、多会话 autosave 覆盖策略、
  24h 回收 / 登出 / OOM 实测
- **验收单**：每项有测试 ID、对应问题、fixture、执行环境、预期数值/状态、
  结果日志、发布 digest

### 7.2 放行标准

沿用审查报告 §8 第三批 9 项清单，并加两条（复审）：

- **对确定的代码缺陷**：先验证测试在基线上暴露问题，再验证修复后通过
- **对新增能力或关闭的功能**：用适用的正反向契约测试，
  不机械要求运维流程套用单元测试"先红后绿"

必须覆盖的场景矩阵（复审 §3，直接采用）：

| 场景 | 关键断言 |
|---|---|
| 正常与恶意列名，多后端 × 多 analysis_type | 无副作用；合法输入有有效结果；contrast/P 值与安全参考模型吻合 |
| 不支持参数：缺省 / NULL / 实际值 | 前两者成功，第三者拒绝；支持的后端确实应用设置 |
| 两层分别提交 QC，往返切层 | 两层已应用数据都保留；主/次层差异与整合使用正确版本 |
| 预览失败后应用、切层后失败 | 不给旧预览重新绑定新身份 |
| 已保存项目在已有会话中恢复 | 三个 QC 控件、已应用数据、结果归属一致，不被当前草稿覆盖 |
| 恢复后改变 QC | 项目、视图、autosave、重启后状态、报告、脚本一致失效 |
| 同数据不同项目；同项目重开旧存档 | 无串写；打开的存档确实生效 |
| 导出脚本独立重放 | 中间矩阵与关键数值吻合；缺记录时明确停止或降级 |
| 异步回调到达时版本已过期 | 成功与失败回调均不覆盖当前状态 |
| 正式镜像与目标服务器 | 多账号隔离、并发资源、备份恢复、日志留存、升级回滚均有证据 |

---

## 8. 决策与待办

| ID | 问题 | 堵住 | 建议 |
|---|---|---|---|
| **D1** | v1 是否放弃 PDF 导出？ | T5.3 | 放弃，只承诺 HTML |
| **D2** | 网页 RDS 是否固定关闭？ | T1.2 | 多用户部署固定关闭，不给开关 |
| **D3** | concordance 的合并检验口径 | T4.1 | 在有经验证方案前，只出描述性一致性 + 分层显著性，关闭合并显著性结论 |
| **D4** | 用户文件误删是否需可恢复？ | T5.2 | 需产品决策 |
| **D5** | §1 数据合同评审 | 批次 2–5 全部 | **先评审，后分派** |
| **D6** | 是否新增 `outlier_sd_threshold` 控件？ | T3.1 | 不新增；参数保持 `run_qc()` 默认并在页面标注 |

**记入待办、不在本次范围**：method × analysis_type × assay scale 的
输入/输出效果量合同（R08）。本次只做导入侧中止，不引入 `allow_raw`。

估工待 §1 评审后重估。v2 的 24–26 人日作废——按层存储、激活/恢复、
provenance 三个接口未统一前，任何数字都不构成上线承诺。

---

## 9. 本版的复现证据

对复审 R01/R02 的独立复现（我自己构造，未使用复审方脚本）：

```text
R01  run-diff.R:217-225  backend_args 恒含 covariates / paired_col 键
     run_ttest_group(input, group_col, control_group, case_group, var_equal, paired_col)
       → 无 covariates 形参
     run_lm_group(input, group_col, control_group, case_group, covariates)
       → 无 paired_col 形参
     v2 补丁 return(args) 原样返回 ⇒ 每次正常调用 unused argument

R02  场景1：组名恰为 "v2" + 一个数值协变量
       设计列名      = A | v2 | v2
       实际 contrast = -1 | 0 | 1        ← 指向协变量系数
       期望 contrast = -1 | 1 | 0
       不报错,出完整结果表和正常火山图,但比较的是错误系数
     场景2：九个协变量
       设计列名 = A | B | v2 | ... | v9 | 0     ← v10 被 startsWith("v1") 吞掉
       makeContrasts: Non-valid names: 0
```

源码核对确认的其余各条：

```text
R03  mod_integration_view.R:82  sec <- exps[[secondary_idx]]     ← 直接读导入层
R04  normalize_omics(input, method = c("vsn","log2"), offset = 1) ← 无 ...
R08  LOG_SCALE_ASSAY_TYPES = normalized_intensity, imputed_intensity,
       filtered_intensity, vst, logcpm                            ← tpm/fpkm 未挡
     effect_type 现有词表 = F_statistic, beta, correlation, effect_diff,
       log2FC, log2FC_per_unit, mean_diff, neg_log10_padj         ← 无 mean_difference
R05  mod_qc_view.R:238  失败时只写 last_error,last_bundle 保留旧值
R10  project_store.R:282  list_saved_projects()                   ← 非 store_list_projects()
     mod_qc_view.R  输入控件只有 missing_threshold(:457) /
       outlier_method(:465) / impute_method(:204)                 ← 无 outlier_sd_threshold
```

上一版记录的 A01–A12 复现证据、QC 参数链路实测、
"覆盖法 vs 双槽位"实测，均保持有效，见 [审查报告](README.md) 与
[复审报告](REMEDIATION-REVIEW.md) 的证据索引。
