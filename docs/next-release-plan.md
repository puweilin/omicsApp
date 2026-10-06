# 下一版本更新建议（v0.2 → v0.4）

依据：三轮系统验证（[review-2026-09.md](./review-2026-09.md)、[review-2026-10.md](./review-2026-10.md)）
中确认但尚未修复的问题，以及审查中发现的设计短板。按"先保证结论正确，再扩展方法，最后做体验与运维"排序。

优先级：**P0** 可能给出错误结论或丢数据/安全问题，下个版本必须做；**P1** 明显影响可用性或可复现性；
**P2** 增强与长期维护。工作量为一名熟悉代码的开发者的粗估。

---

## 1. 版本规划总览

| 版本 | 主题 | 主要内容 | 粗估 |
|---|---|---|---|
| **v0.2.0** | 稳定与可复现 | 第 2 节全部 P0（✅ 已完成，见 1.1）；整合分析脚本导出；恢复项目后视图重建；上传加固；CI 与生产镜像一致；版本号/README 更新 | 3–4 周 |
| **v0.3.0** | 统计方法扩展 | voom、lfcShrink、交互/因子设计与公式输入、批次校正、时间序列、蛋白组专用方法、GSVA 差异分析、>2 层整合 | 5–7 周 |
| **v0.4.0** | 协作与体验 | 中文界面、无障碍、移动端、项目版本/未保存提醒、多用户隔离、部署加固与备份、性能优化 | 4–6 周 |

---

## 1.1 v0.2 P0 完成情况（2026-10）

第 2 节中标为 **P0 ✅** 的条目已全部完成，每项都有回归测试：

| 条目 | 做法 | 测试 |
|---|---|---|
| 样本在行的表（Olink）方向 | 新增"样本名共享前缀"判断（S01…、Patient_3…，基因名没有共同前缀），Olink 宽表能正确识别；置信度 < 0.6 时导入报告标明"是猜的"；`read_omics(orientation = )` 可显式指定并写进可复现脚本；导入页新增方向选择（猜测时以警告突出显示） | `test-audit-import.R` |
| 全数值样本表压过矩阵 | 若一张"像矩阵"的表，其取值是另一张表的样本名，就判为样本表；带样本表表头的全数值表降低置信度；多个候选时取更大、更可信的表，并在报告中列出 | `test-audit-import.R` |
| 整合分析脚本导出 | 整合 bundle 记录每层差异分析的参数（含应用内计算的伴随层）；脚本逐层重算 `run_diff()`（共享拟合后 `select_comparison()`）、写出 `sample_link`、方法参数（`cor_method`、`min_samples`、阈值、数据库等）全部写出、按方法生成对应图；导出脚本端到端运行，结果与项目中保存的一致（concordance 与 correlation 两种方法）。顺带修复：多对比结果的火山图在脚本中报错 | `test-audit-enrich.R` |
| 上传加固 | `read_omics()` 读取前检查压缩包声明的解压后大小（默认 2 GB，`options(omicsCore.max_unpacked_mb)`）与条目数，拦截 zip/gzip 炸弹；导入解析移到后台 future，不再冻结会话；新文件覆盖旧解析（过期结果丢弃）；RDS 中的 `omics_input` 重新校验 | `test-audit-import.R`、`test-audit-app.R` |
| 多标签页共用自动保存 | 每个浏览器会话写自己的 `_autosave-<id>.omp`（附带名称/图层数的 json 说明），保留当前会话外最近 3 份；"恢复"在有多份时弹出选择框；到达时自动恢复最新一份 | `test-audit-app.R`，原有 store/acceptance 测试同步更新 |
| 备份 | `deploy/scripts/backup.sh`：`set -Eeuo pipefail`；按日期的快照（`--link-dest` 硬链接去重），14 天日备 + 8 周周备；`pg_dump` 先写临时文件并校验后才替换；配置/证书/cron/基因集缓存一并备份；`MANIFEST.sha256` 与镜像 digest；`BACKUP_REMOTE` 异地复制（保留历史）；任何失败经 webhook/邮件/syslog 告警并非零退出；`restore_check.sh` 每周演练（备份时效、校验和、在临时 Postgres 中恢复账户库并计数、用镜像打开一个项目、核对异地副本） | `test-deploy-contract.R`（含端到端运行） |
| CI 与生产镜像一致 | 新增 `.github/workflows/production-image.yaml`：PR/主干/每周构建 `deploy/docker/Dockerfile`（GHA 缓存），启动容器做服务冒烟测试，并在镜像内运行两个包的全部测试 | `test-deploy-contract.R` |

## 1.2 第四轮系统测试后的进展（2026-10）

第四轮从准确性、效率、UI、用户体验四个视角做了完整测试，详见 [review-2026-10-round4.md](./review-2026-10-round4.md)。
本计划中的以下条目已一并完成：恢复项目后视图重建结果（2.5 P1）、QC 默认同时运行三种离群检测与按图层排除样本（2.2）、
RDS 输入重新校验（2.1 P1）、t 检验/lm/连续变量的向量化与项目文件瘦身（2.6 P1 中 4 项）、配对设计与连续变量进入界面、
单独样本表上传、物种选择、完整结果表下载、面向生物学用户的报告、导出脚本连同数据下载。

## 1.3 第四轮报告"下一版待办"完成情况（2026-10）

第四轮报告中列为"未修复 / 未采用"的 9 项已全部完成，详见
[review-2026-10-round4.md 的"后续：下一版待办"](./review-2026-10-round4.md)：
长标签与手机屏幕上的图表可读性、favicon 404、运行中按钮禁用与分步进度（DESeq2 / edgeR / QC / 报告下载）、
术语统一、基因集表缓存、`fread` 读取 CSV（同时修复带引号表头的 P1 问题，见 2.1）、DESeq2 可选并行、
自动保存合并写入。本节以下列表中 2.1 的"表头带引号且含分隔符的 CSV 崩溃"与 2.5 的"DESeq2 等长任务没有进度"随之完成。

## 1.4 发布计划中全部 P1 完成情况（2026-10）

第 2 节中所有 P1 条目均已完成（含此前已在第四轮完成、本次补标的 RDS 重新校验、恢复项目后视图回填、四项性能优化）。
每项都有回归测试。合并后的全量测试：omicsCore 与 omicsApp（含无头浏览器测试）均 0 失败；跳过的仅为需要真实数据集
（SkinProteomics、卵泡 RNA-seq）或需手动开启（性能预算、全量参数模糊测试）的测试。要点如下：

**导入（2.1）**
- 非 UTF-8 文件：按字节识别编码（UTF-16 BOM → UTF-8 → GB18030/GBK → Windows-1252/Latin-1），转为 UTF-8 后读取，矩阵与样本表都适用，导入页提示所用编码。
- Salmon / RSEM / kallisto：按列特征识别；新函数 `read_quant_files()` 合并多个样本文件（计数取 NumReads/expected_count/est_counts，有效长度存为 tximport 元数据，DESeq2/edgeR 自动使用长度校正），支持 tx2gene 汇总到基因（与 tximport 算法一致，不依赖 tximport 包）；导入页可一次选择多个文件；导出脚本用 `read_quant_files()` 重读归档的原文件。
- 汇总列（Total、Mean、Sum、SD 等）按列名或"等于其他列的行和/行均值"识别并剔除。
- RNA-seq 数值尺度：整数 → 计数；列和约 1e6 → TPM；有负值或最大值 < 30 → log 尺度；其余 → FPKM。导入页显示判断依据，可修改；非计数数据不提供 DESeq2/edgeR。
- 上传的 `.rds` 先做与不可信项目相同的"只含数据"结构检查。

**QC（2.2）**
- 新增留一法（leave-one-out）离群检测：比较每个样本与最近邻样本的距离和其余样本的同一分布（中位数/MAD），4 个样本即可使用；模拟中干净数据误报 0–1%，明显异常样本检出 ≥ 86–100%。加入默认的"全部方法"，其他三种方法结果不变。
- 按组缺失率过滤：`missing_filter = "any_group"`（至少一组满足）/ `"all_groups"`（每组都满足），默认仍为全局；界面可选，导出脚本可复现。
- 打开/恢复项目时显示保存的 QC 结果并恢复其参数，只有参数改变时才重算。

**差异分析（2.3）**
- 新增统一的带符号统计量列 `signed_stat`（limma t、DESeq2 Wald、edgeR sign(logFC)·√F、t 检验/线性模型 t），GSEA 排序优先使用；旧结果兼容。
- edgeR 全局检验使用与两组比较相同的 tximport 长度校正（DESeq2 本已使用）。
- 配对 t 检验按每个对比单独配对（单个对比也一样），并说明每个对比使用/剔除了哪些配对。

**富集与整合（2.4）**
- ORA 默认上调、下调分别富集（各自校正，结果中注明来自哪一组），合并作为可选项；旧结果与脚本含义不变。
- 多数据库时可选跨库 BH 校正（`p_adjust_scope = "all"`）；界面目前只选单个数据库，故只在核心函数与脚本中提供。
- ActivePathways 使用方向性方法（DPM，要求两层方向一致），结果给出每条通路在两层中的方向及是否一致。
- 物种：新增 `enrichment_species()`（人、小鼠、大鼠、斑马鱼、果蝇、酵母、线虫，及 msigdbr 支持的其他物种），使用 msigdbr 的直系同源映射；基因名大小写不匹配但忽略大小写匹配良好时自动忽略大小写并提示。

**Shiny 应用（2.5）**
- 同一组学类型可以有多个图层：导入时填写图层名（默认组学类型）；同名时可"替换"或"保留两者"（新图层命名为 `<名称>_2`）。所有结果记录所属图层（新函数 `bundle_layer()`），各视图、导出脚本按图层而非组学类型匹配；替换/删除图层只清除该图层的结果。
- 结果过期提示：差异页、富集页在控件（或差异阈值）改变后提示"结果对应之前的设置，请重新运行"；整合页随输入自动重算、QC 页实时计算，不会出现过期结果。

**部署（2.7）**
- nginx 只允许 `ADMIN_ALLOW_CIDR` 访问 Keycloak 管理端点、ShinyProxy `/admin`、`/actuator` 等；应用端口只绑定本机。
- 容器加固：ShinyProxy 支持的内存/CPU 限制、非特权、专用网络；镜像去除所有 setuid/setgid 与文件能力，进程数上限；Keycloak/Postgres 服务 no-new-privileges、去除能力、pids/内存/CPU 限制。ShinyProxy 的 Docker 后端不支持只读根文件系统、cap-drop 等选项，已在 README 说明替代方案（主机级 `daemon.json`）。
- 出站限制：`egress.sh` + systemd 单元阻止应用容器网络新建出站连接。
- 日志持久化与轮转，并纳入备份。
- 回滚：镜像标签为 `版本-提交号`，不再使用可变标签；`rollback.sh <tag>` 一键切换；基础镜像与 Postgres 按 digest 固定（`pin_base_digests.sh` 维护锁文件；Keycloak 因网络限制暂未固定，README 有 TODO）。
- `.omp` 项目文件记录格式版本，`load_project()` 按版本链执行迁移。
- 项目文件签名：保存时附加 HMAC-SHA256（密钥来自 `OMICSAPP_SIGNING_KEY` 或数据目录中自动生成的密钥），打开前验证；未签名或被篡改的文件被拒绝；首次运行时为通过结构检查的旧文件补签名。

**打包（2.8）**
- limma、clusterProfiler、msigdbr 移到 Suggests（调用处均已有安装检查，缺 limma 时 `auto` 退回 t 检验）；删除未使用的 Suggests（here、ggpubr、tximport、GenomicFeatures；应用中的 GSVA、fgsea、enrichplot、ComplexHeatmap、circlize、ggrepel、knitr）。
- 两个包版本号 0.2.0，新增 NEWS.md，README 更新，新增导出函数 `shiny_app()`，`inst/app/app.R` 不再使用 `:::`；删除误提交的文件。

**已知限制 / 后续**
- 多数据库跨库校正尚无界面入口；KEGG 在线刷新只支持人/小鼠；部署预热只包含人/小鼠基因集。
- 留一法离群检测不能发现"两组之间互换"的样本；单样本组可能被标为离群（仍保留）。
- 项目只保存一个差异结果槽位：在两个图层上分别做差异分析时，项目中保留最近一次。
- 部署改动无法在本环境中对真实 nginx/Docker/ShinyProxy 验证，需要在服务器上执行 README 中的检查步骤。

---

## 2. 已知问题清单（✅ = 已完成）

### 2.1 导入
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 ✅ | 样本在行、特征在列的表（如 Olink NPX 宽表）方向识别置信度默认 0.4，会被转置读入且无提示 | 置信度 < 0.6 或首列像样本名时在导入页要求确认方向 |
| P0 ✅ | 全数值的元数据 sheet 可能被识别成表达矩阵，压过真正的矩阵 sheet | 对行/列很少的 sheet 降权；多个候选时让用户选择 |
| P1 ✅ | CP1252/GBK 编码的 CSV 直接崩溃 | 依次尝试 UTF-8 → 本地编码 → latin1，报告使用的编码 |
| P1 ✅ | 表头带引号且含分隔符的 CSV 崩溃（自写的 `strsplit` 解析） | 改用 `utils::count.fields()` / `data.table::fread()` 解析表头 |
| P1 ✅ | 单样本的 Salmon / RSEM / kallisto 文件被读成 4 个"样本"（各列） | 识别这些格式，支持多文件合并（tximport） |
| P1 ✅ | `Total`/`Mean` 等汇总**列**被当作样本（汇总行已处理） | 与汇总行同样识别并剔除 |
| P1 ✅ | `infer_assay_type()` 对 RNA-seq 永远返回 `raw_count`；非整数（TPM、log 值）也被当作计数 | 依据整数性、取值范围、负值判断，并在导入页确认 |
| P1 ✅ | RDS 导入的 `omics_input` 不重新校验 | 走 `validate_omics_input()` |
| P2 | 不支持 `.gz`、`.xlsm`、`SummarizedExperiment`、`DESeqDataSet`、小鼠 ENSMUSG、`_PAR_Y` 后缀 | 逐项补充读取器与 ID 映射 |
| P2 | `winsorize_counts()` 会把"开/关"型基因压成 0，并产生非整数 | 仅对表达基因、按整数截断 |
| P2 | MinProb 在特征很少的数据上报错或无效 | 回退到 MinDet 并提示 |

### 2.2 QC
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 ✅ | 样本 ≤ 10 时 z 分数离群检测不可能触发（本轮已加说明，但没有替代方法） | 小样本用稳健方法：`rrcov::PcaHubert` 或基于留一法的 Grubbs 型检验 |
| P1 ✅ | 缺失率过滤是全局的，不按组 | 增加"至少一组中 ≥ x% 有值"的过滤（蛋白组标准做法） |
| P1 ✅ | QC 视图每次切换图层都用默认参数重算，覆盖恢复项目里保存的 QC 结果 | 有保存结果时先展示保存的结果，参数改变后再重算 |
| P2 | log2 归一化没有中位数对齐选项 | `normalize_omics(method = "log2", center = "median")` |

### 2.3 差异分析
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 ✅ | edgeR 结果的 `statistic` 是无符号 F，与 limma 的 t、DESeq2 的 Wald 含义不同 | 结果表增加带符号统计量列（GSEA 已在内部处理） |
| P1 ✅ | 全局检验（edgeR/DESeq2）未使用 tximport offset | 与两组比较一致处理 |
| P1 ✅ | 多个处理组时配对 t 检验要求过严 | 每个对比单独配对 |
| P2 | DESeq2 连续变量分析没有测试覆盖 | 补测试 |

### 2.4 富集与整合
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 ✅ | 整合分析导出脚本不完整：`cor_method`、`min_samples`、阈值、数据库等方法参数丢失；`sample_link` 未写出；伴随层的差异分析只在会话缓存里，脚本不可重现；作图代码总是写 dual_volcano | 整合 bundle 记录全部方法参数与伴随差异 bundle；脚本按方法生成对应图 |
| P1 ✅ | ORA 默认"both"把上调下调合并 | 默认分别做上/下调，合并作为选项 |
| P1 ✅ | 多数据库同时富集时没有跨库多重校正 | 提供跨库 BH 选项 |
| P1 ✅ | ActivePathways 的结果没有方向（只说"显著"） | 使用 directional ActivePathways（`merge_method = "DPM"`） |
| P1 ✅ | 非人物种依赖基因符号大小写，没有直系同源映射 | 引入 `msigdbr` 物种映射或 `babelgene` |
| P2 | 蛋白—基因匹配只按 symbol，同一基因的多个蛋白异构体只保留一个 | 引入 UniProt ↔ Gene 映射表与 feature_link |
| P2 | GSEA 未暴露 `eps` / `nPermSimple` | 透传参数 |

### 2.5 Shiny 应用
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 ✅ | 上传文件同步解析、无大小上限（xlsx/gzip 炸弹可拖垮会话） | `shiny.maxRequestSize` 分级上限；解析放到 future；解压后大小检查 |
| P0 ✅ | 同一用户两个浏览器标签共用一个自动保存文件，互相覆盖 | 自动保存按会话/标签区分，恢复时让用户选 |
| P1 ✅ | 打开或恢复项目后，各视图不从保存的结果重建（需要重新运行） | 视图启动时从 `project$bundles` 回填（rehydrate），并显示"已恢复"横幅 |
| P1 ✅ | 一个组学类型只能有一个图层（图层名 = omics_type） | 允许多个同类型图层（如两批蛋白组），图层名由用户指定 |
| P1 ✅ | 结果过期没有状态提示（改了参数但未重跑） | 控件变化后标记"结果对应旧参数" |
| P1 ✅ | DESeq2 等长任务没有进度 | 分阶段进度（离散度/拟合/检验） |
| P2 | QC 参数修改无防抖，大数据集上每动一次滑块就重算 | `debounce()` |

### 2.6 性能（实测 8000×60 蛋白 / 30000×60 RNA）
| 级别 | 热点 | 建议（已验证的加速） |
|---|---|---|
| P1 ✅ | 项目文件 122 MB，保存 2 秒：DESeqDataSet（57 MB）与 DGEGLM（40 MB）模型对象随项目保存 | 保存时不存模型对象或只存必要部分 → 31 MB，快 6 倍 |
| P1 ✅ | limma 连续变量分析逐特征 `lm` + `cor.test`：8.5 秒 | 闭式解，约 40 倍 |
| P1 ✅ | t 检验 / lm 逐特征循环：0.7 秒 / 8.4 秒 | 向量化（t 检验 0.014 秒） |
| P1 ✅ | QC 与 PCA 图每次都做完整 `prcomp`：RNA 5–7 秒 | Gram 矩阵特征分解（17 倍）或 `irlba`，并把得分缓存进 QC bundle |
| P2 | DESeq2 24 秒、edgeR 8.8 秒 | 预过滤（DESeq2 → 13 秒）、`glmGamPoi` |
| P2 | QC bundle 复制了一份 `cleaned_input`（25 MB） | 只存过滤/插补记录，按需重建 |

### 2.7 部署与 CI
| 级别 | 问题 | 建议 |
|---|---|---|
| P0 ✅ | 备份与生产在同一主机；`rsync --delete` 镜像没有历史版本；`pg_dump` 无 `pipefail`，失败时会覆盖好的备份；无告警；基因集卷、`app.yml`、TLS、Keycloak 配置未备份 | 异地、带版本保留的备份；`set -o pipefail`；失败告警；补齐备份清单并定期恢复演练 |
| P0 ✅ | CI 用最新 R-release 测试，而生产镜像是 R 4.4.2 / Bioconductor 3.20 快照；CI 从不构建 Dockerfile | CI 在生产镜像中跑测试，PR 上构建镜像 |
| P1 ✅ | `/auth/admin` 在局域网可访问；容器无 `no-new-privileges`、pids 限制、出站限制；日志未持久化 | 反向代理屏蔽管理端点；容器加固；日志落盘/集中 |
| P1 ✅ | 无回滚：镜像标签可变、基础镜像按标签引用、`.omp` 格式无迁移钩子 | 不可变标签 + digest 固定；`load_project()` 版本迁移 |
| P1 ✅ | `load_project()` 反序列化不可信的 qs2 文件 | 只加载本服务写出的文件（签名/哈希），上传的项目先校验结构 |
| P2 | 缺 HSTS/安全头/限流；CI 只在 main 触发、单一 OS/R 版本、42 分钟；无 lint/覆盖率/shellcheck/hadolint；actions 未固定 SHA；19 个部署契约测试与 8 个 hygiene 测试在 check 下被跳过；性能/模糊测试从不运行 | 逐项补齐；性能与 fuzz 放到 nightly |

### 2.8 打包与可维护性
| 级别 | 问题 | 建议 |
|---|---|---|
| P1 ✅ | limma、clusterProfiler（约 120 个依赖）、msigdbr 在 Imports，与 README 所说"可选"不符 | 移到 Suggests，通过 `ensure_*()` 按需提示安装 |
| P1 ✅ | 未使用的 Suggests（here、ggpubr、tximport、GenomicFeatures；app 中 GSVA、fgsea、enrichplot、ComplexHeatmap 等） | 清理 |
| P1 ✅ | 版本号仍是 0.0.0.9000；README 停留在 "Phase 0"；`omicsApp-package.Rd` 过时；`inst/app/app.R` 使用 `:::` | 发布 v0.2.0 时统一更新 |
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

- [x] 第 2 节所有 P0 关闭
- [ ] `R CMD check --as-cran` 两个包 0 ERROR / 0 WARNING
- [ ] omicsCore、omicsApp 全量测试（含浏览器测试）在生产镜像中通过
- [ ] `export_script()` 自检在示例项目和两个真实项目上通过
- [ ] 备份恢复演练一次
- [ ] README、tutorial、CHANGELOG（NEWS.md）、版本号更新
- [ ] 旧版本 `.omp` 项目可加载（迁移测试）
