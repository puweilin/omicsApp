# 下一版本更新建议（v0.2 → v0.4）

依据：三轮系统验证（[review-2026-09.md](./review-2026-09.md)、[review-2026-10.md](./review-2026-10.md)）
中确认但尚未修复的问题，以及审查中发现的设计短板。按"先保证结论正确，再扩展方法，最后做体验与运维"排序。

优先级：**P0** 可能给出错误结论或丢数据/安全问题，下个版本必须做；**P1** 明显影响可用性或可复现性；
**P2** 增强与长期维护。工作量为一名熟悉代码的开发者的粗估。

---

## 1. 版本规划总览

| 版本 | 主题 | 主要内容 | 粗估 |
|---|---|---|---|
| **v0.2.0** | 稳定与可复现 | 第 2 节全部 P0；整合分析脚本导出；恢复项目后视图重建；上传加固；CI 与生产镜像一致；版本号/README 更新 | 3–4 周 |
| **v0.3.0** | 统计方法扩展 | voom、lfcShrink、交互/因子设计与公式输入、批次校正、时间序列、蛋白组专用方法、GSVA 差异分析、>2 层整合 | 5–7 周 |
| **v0.4.0** | 协作与体验 | 中文界面、无障碍、移动端、项目版本/未保存提醒、多用户隔离、部署加固与备份、性能优化 | 4–6 周 |

---

## 2. 已知问题（确认但未修复）

### 2.1 导入
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 | 样本在行、特征在列的表（如 Olink NPX 宽表）方向识别置信度默认 0.4，会被转置读入且无提示 | 置信度 < 0.6 或首列像样本名时在导入页要求确认方向 |
| P0 | 全数值的元数据 sheet 可能被识别成表达矩阵，压过真正的矩阵 sheet | 对行/列很少的 sheet 降权；多个候选时让用户选择 |
| P1 | CP1252/GBK 编码的 CSV 直接崩溃 | 依次尝试 UTF-8 → 本地编码 → latin1，报告使用的编码 |
| P1 | 表头带引号且含分隔符的 CSV 崩溃（自写的 `strsplit` 解析） | 改用 `utils::count.fields()` / `data.table::fread()` 解析表头 |
| P1 | 单样本的 Salmon / RSEM / kallisto 文件被读成 4 个"样本"（各列） | 识别这些格式，支持多文件合并（tximport） |
| P1 | `Total`/`Mean` 等汇总**列**被当作样本（汇总行已处理） | 与汇总行同样识别并剔除 |
| P1 | `infer_assay_type()` 对 RNA-seq 永远返回 `raw_count`；非整数（TPM、log 值）也被当作计数 | 依据整数性、取值范围、负值判断，并在导入页确认 |
| P1 | RDS 导入的 `omics_input` 不重新校验 | 走 `validate_omics_input()` |
| P2 | 不支持 `.gz`、`.xlsm`、`SummarizedExperiment`、`DESeqDataSet`、小鼠 ENSMUSG、`_PAR_Y` 后缀 | 逐项补充读取器与 ID 映射 |
| P2 | `winsorize_counts()` 会把"开/关"型基因压成 0，并产生非整数 | 仅对表达基因、按整数截断 |
| P2 | MinProb 在特征很少的数据上报错或无效 | 回退到 MinDet 并提示 |

### 2.2 QC
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 | 样本 ≤ 10 时 z 分数离群检测不可能触发（本轮已加说明，但没有替代方法） | 小样本用稳健方法：`rrcov::PcaHubert` 或基于留一法的 Grubbs 型检验 |
| P1 | 缺失率过滤是全局的，不按组 | 增加"至少一组中 ≥ x% 有值"的过滤（蛋白组标准做法） |
| P1 | QC 视图每次切换图层都用默认参数重算，覆盖恢复项目里保存的 QC 结果 | 有保存结果时先展示保存的结果，参数改变后再重算 |
| P2 | log2 归一化没有中位数对齐选项 | `normalize_omics(method = "log2", center = "median")` |

### 2.3 差异分析
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 | edgeR 结果的 `statistic` 是无符号 F，与 limma 的 t、DESeq2 的 Wald 含义不同 | 结果表增加带符号统计量列（GSEA 已在内部处理） |
| P1 | 全局检验（edgeR/DESeq2）未使用 tximport offset | 与两组比较一致处理 |
| P1 | 多个处理组时配对 t 检验要求过严 | 每个对比单独配对 |
| P2 | DESeq2 连续变量分析没有测试覆盖 | 补测试 |

### 2.4 富集与整合
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 | 整合分析导出脚本不完整：`cor_method`、`min_samples`、阈值、数据库等方法参数丢失；`sample_link` 未写出；伴随层的差异分析只在会话缓存里，脚本不可重现；作图代码总是写 dual_volcano | 整合 bundle 记录全部方法参数与伴随差异 bundle；脚本按方法生成对应图 |
| P1 | ORA 默认"both"把上调下调合并 | 默认分别做上/下调，合并作为选项 |
| P1 | 多数据库同时富集时没有跨库多重校正 | 提供跨库 BH 选项 |
| P1 | ActivePathways 的结果没有方向（只说"显著"） | 使用 directional ActivePathways（`merge_method = "DPM"`） |
| P1 | 非人物种依赖基因符号大小写，没有直系同源映射 | 引入 `msigdbr` 物种映射或 `babelgene` |
| P2 | 蛋白—基因匹配只按 symbol，同一基因的多个蛋白异构体只保留一个 | 引入 UniProt ↔ Gene 映射表与 feature_link |
| P2 | GSEA 未暴露 `eps` / `nPermSimple` | 透传参数 |

### 2.5 Shiny 应用
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 | 上传文件同步解析、无大小上限（xlsx/gzip 炸弹可拖垮会话） | `shiny.maxRequestSize` 分级上限；解析放到 future；解压后大小检查 |
| P0 | 同一用户两个浏览器标签共用一个自动保存文件，互相覆盖 | 自动保存按会话/标签区分，恢复时让用户选 |
| P1 | 打开或恢复项目后，各视图不从保存的结果重建（需要重新运行） | 视图启动时从 `project$bundles` 回填（rehydrate），并显示"已恢复"横幅 |
| P1 | 一个组学类型只能有一个图层（图层名 = omics_type） | 允许多个同类型图层（如两批蛋白组），图层名由用户指定 |
| P1 | 结果过期没有状态提示（改了参数但未重跑） | 控件变化后标记"结果对应旧参数" |
| P1 | DESeq2 等长任务没有进度 | 分阶段进度（离散度/拟合/检验） |
| P2 | QC 参数修改无防抖，大数据集上每动一次滑块就重算 | `debounce()` |

### 2.6 性能（实测 8000×60 蛋白 / 30000×60 RNA）
| 级别 | 热点 | 建议（已验证的加速） |
|---|---|---|
| P1 | 项目文件 122 MB，保存 2 秒：DESeqDataSet（57 MB）与 DGEGLM（40 MB）模型对象随项目保存 | 保存时不存模型对象或只存必要部分 → 31 MB，快 6 倍 |
| P1 | limma 连续变量分析逐特征 `lm` + `cor.test`：8.5 秒 | 闭式解，约 40 倍 |
| P1 | t 检验 / lm 逐特征循环：0.7 秒 / 8.4 秒 | 向量化（t 检验 0.014 秒） |
| P1 | QC 与 PCA 图每次都做完整 `prcomp`：RNA 5–7 秒 | Gram 矩阵特征分解（17 倍）或 `irlba`，并把得分缓存进 QC bundle |
| P2 | DESeq2 24 秒、edgeR 8.8 秒 | 预过滤（DESeq2 → 13 秒）、`glmGamPoi` |
| P2 | QC bundle 复制了一份 `cleaned_input`（25 MB） | 只存过滤/插补记录，按需重建 |

### 2.7 部署与 CI
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 | 备份与生产在同一主机；`rsync --delete` 镜像没有历史版本；`pg_dump` 无 `pipefail`，失败时会覆盖好的备份；无告警；基因集卷、`app.yml`、TLS、Keycloak 配置未备份 | 异地、带版本保留的备份；`set -o pipefail`；失败告警；补齐备份清单并定期恢复演练 |
| P0 | CI 用最新 R-release 测试，而生产镜像是 R 4.4.2 / Bioconductor 3.20 快照；CI 从不构建 Dockerfile | CI 在生产镜像中跑测试，PR 上构建镜像 |
| P1 | `/auth/admin` 在局域网可访问；容器无 `no-new-privileges`、pids 限制、出站限制；日志未持久化 | 反向代理屏蔽管理端点；容器加固；日志落盘/集中 |
| P1 | 无回滚：镜像标签可变、基础镜像按标签引用、`.omp` 格式无迁移钩子 | 不可变标签 + digest 固定；`load_project()` 版本迁移 |
| P1 | `load_project()` 反序列化不可信的 qs2 文件 | 只加载本服务写出的文件（签名/哈希），上传的项目先校验结构 |
| P2 | 缺 HSTS/安全头/限流；CI 只在 main 触发、单一 OS/R 版本、42 分钟；无 lint/覆盖率/shellcheck/hadolint；actions 未固定 SHA；19 个部署契约测试与 8 个 hygiene 测试在 check 下被跳过；性能/模糊测试从不运行 | 逐项补齐；性能与 fuzz 放到 nightly |

### 2.8 打包与可维护性
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 | limma、clusterProfiler（约 120 个依赖）、msigdbr 在 Imports，与 README 所说"可选"不符 | 移到 Suggests，通过 `ensure_*()` 按需提示安装 |
| P1 | 未使用的 Suggests（here、ggpubr、tximport、GenomicFeatures；app 中 GSVA、fgsea、enrichplot、ComplexHeatmap 等） | 清理 |
| P1 | 版本号仍是 0.0.0.9000；README 停留在 "Phase 0"；`omicsApp-package.Rd` 过时；`inst/app/app.R` 使用 `:::` | 发布 v0.2.0 时统一更新 |
| P2 | 模块服务器过大（diff 约 1100 行、import、integration 各 650+ 行） | 按卡片拆分子模块 |
| P2 | `%||%` 定义了 8 次；各后端的设计矩阵准备与 8 个 standardize 函数大量重复 | 抽公共工具函数 |
| P2 | 测试 `setup.R` 依赖 `load_all`，与安装后测试耦合 | 测试只依赖已安装包 |

---

## 3. 方法扩展建议（v0.3.0）

### 3.1 差异分析
1. **voom / voomWithQualityWeights**：RNA-seq 的 limma 路线，样本质量不一时显著优于 log-CPM + limma。
2. **lfcShrink（apeglm/ashr）**：DESeq2 结果的效应量收缩，用于排序、火山图与 GSEA。
3. **公式输入与因子/交互设计**：例如 `~ genotype * treatment`，支持"处理效应在两种基因型中是否不同"。
   界面上提供"高级：模型公式"折叠区，后端统一走 `model.matrix` + 对比字符串（已有的自定义对比可以复用）。
4. **批次校正视图**：`limma::removeBatchEffect` 仅用于可视化，`ComBat-seq` 可选；明确提示"统计检验应把批次放进模型而不是先校正"。
5. **时间序列**：样条 × 组的交互（limma），以及 DESeq2 的 LRT 时间效应。
6. **混合模型**：`dream`（variancePartition）处理重复测量/多层随机效应。
7. **蛋白组专用**：DEqMS（按肽段数校正方差）、proDA（概率性缺失模型，不需插补）、msqrob2。
8. **独立过滤**：DESeq2 的 independent filtering 参数暴露到界面。

### 3.2 富集与整合
1. **GSVA/ssGSEA 进入应用**，并支持对 GSVA 分数做差异分析（limma）。
2. **多层整合**：MOFA2 / DIABLO（mixOmics），支持 > 2 个图层。
3. **directional ActivePathways** 与方向一致性的整合 p 值。
4. 富集结果的网络图/emap 图（enrichplot）。

### 3.3 新组学类型
代谢组（log 强度、缺失多为 MNAR，可复用蛋白组路径），以及 Olink NPX（已是 log2，样本在行）。

---

## 4. 体验与协作建议（v0.4.0）

1. **中文界面（i18n）**：用 `shiny.i18n` 或自建字典，所有面向用户的字符串集中管理；图形使用支持 CJK 的字体。
2. **无障碍**：颜色之外的标记（形状/文字）、plotly 图提供文字摘要、`aria-live` 通知、键盘焦点顺序。
3. **移动端/窄屏布局**：侧栏折叠、表格横向滚动。
4. **项目管理**：当前项目指示、未保存更改提醒、项目版本历史与回退、导出/导入 `.omp`。
5. **结果状态**：每张结果卡片显示"基于哪个图层、什么参数、何时运行"，参数变化后标记过期。
6. **教程**：在示例项目上加入"多组设计""配对设计""连续变量"三条引导路线。

---

## 5. 可复现性与质量保障

1. **脚本自检**：`export_script()` 生成后，在隔离进程中运行，对比主要结果（差异表 p 值、富集通路）与项目中保存的结果，
   不一致时在报告中标红。本轮已修好归一化、研究设计、sheet 角色的导出，整合分析是剩下的主要缺口。
2. **结果溯源**：每个 bundle 记录 omicsCore 版本、输入指纹、随机种子；报告页脚统一显示。
3. **金标准数据集回归**：选 2–3 个公开数据集（如 airway RNA-seq、一个 DIA 蛋白组），固定期望结果，作为 nightly 测试。
4. **统计正确性测试**：本轮新增的"零假设下假阳性率""已知效应可检出"类测试推广到每个后端。

---

## 6. v0.2.0 发布前检查清单

- [ ] 第 2 节所有 P0 关闭
- [ ] `R CMD check --as-cran` 两个包 0 ERROR / 0 WARNING
- [ ] omicsCore、omicsApp 全量测试（含浏览器测试）在生产镜像中通过
- [ ] `export_script()` 自检在示例项目和两个真实项目上通过
- [ ] 备份恢复演练一次
- [ ] README、tutorial、CHANGELOG（NEWS.md）、版本号更新
- [ ] 旧版本 `.omp` 项目可加载（迁移测试）
