**RNA-seq 差异分析耗时与富集任务体积核查 · 2026-09-08**

结论：富集失败的直接原因是任务输入过大；其来源是将完整 DESeq2 拟合对象作为差异结果包的一部分继续传递，模型的设计公式还保留了计算时的局部环境。已在独立模拟数据上复现相同的六个 globals 和体积膨胀机制。差异分析慢还涉及模型计算与容器资源，需要服务器分阶段计时，不能仅凭该错误确定耗时占比。

本次为诊断，未修改产品代码、未重启容器。用户确认运行环境为测试服务器 Docker，方法为 DESeq2。源码基线 HEAD 为 `2a7d660fd50bf1e6f1c915cda1e1dad6b6d8f3ba`，同时存在正在编辑的工作区修改；验证使用独立临时副本。未取得服务器实际镜像摘要、会话加载的包版本、原始数据或服务器计时日志。本结论不是全应用部署验收。

核查结束时，工作区公式隔离与 QC 改动已提交为 `4d69fa7`。已逐一比对本次测量涉及的五个关键源文件（run-diff、diff-deseq2、diff-design、async_helpers、mod_enrich_view），其内容与独立测量快照完全一致。以下“工作区”描述指复制快照时的状态；新增 QC 的其他改动仍未作完整验收。

**已确认的故障链**

1. `mod_enrich_view.R:104` 获取完整 `diff_bundle()`，在约 140–185 行将整个 bundle 放入 `detached_call()`。ORA/GSEA 读取标准化差异表及少量元信息，不读取 `diff_object`。ORA 在 `enrich-ora.R:63` 由全表形成背景集；GSEA 使用全表生成排序。因此富集需要全部受检基因的结果，通常不需要基因 × 样本的计数矩阵或拟合对象。
2. `async_helpers.R:132` 的 `detached_call()` 已隔离 Shiny 模块环境，但仍完整携带显式传入的 bundle。`run_async()` 在约 88 行创建 future 时扫描任务。模拟扫描得到 `func, type, bundle, db_arg, dir_, thr`，与用户错误的六个 globals 对应。
3. 用户错误中的 `func ≈ 949.65 MiB` 是函数连同环境的体积，绝非函数正文有这么大；`bundle ≈ 949.65 MiB` 是同份大对象又作为单独依赖计入。两者解释了约 1.85 GiB 的总额。这是 future 的导出体积统计，不可直接解释为实际内存一定复制了两份。
4. 旧实现始终将后端模型写入 `results$diff_object`。DESeq2 模型含计数、拟合量等数据；`as.formula()` 默认使用调用环境，而 DESeqDataSet 保存设计公式，进而保留后端局部变量。工作区修改使用 `as.formula(..., env = baseenv())`，已针对其实际生成的简单列名公式验证。
5. 错误发生在任务提交阶段，富集算法尚未开始。用户看到的 500 MiB 是 future 的全局对象导出限制，不是 Docker 的内存配额。[future 官方说明](https://future.futureverse.org/reference/zzz-future.options.html)

**独立复现结果**

本机 R 4.4.3、DESeq2 1.46.0、future 1.69.0；模拟 1,200 个基因 × 200 个样本，两组各 100 例。使用实际 `run_diff()` 与等价于 App 的富集闭包，通过 `future::getGlobalsAndPackages()` 测量。下表不是用户数据的测量值，也不是服务器速度基准。

| 结果包构成 | bundle 的 future 体积 | 整个任务 globals 总体积 |
| --- | ---: | ---: |
| 完整模型，原有公式环境 | 23.69 MiB | 47.37 MiB |
| 完整模型，公式环境改为 baseenv | 8.13 MiB | 16.26 MiB |
| 不携带完整模型，保留两份结果表 | 0.29 MiB | 0.58 MiB |
| 仅保留富集所需表和元信息 | 0.22 MiB | 0.44 MiB |

实际执行了保留模型/不保留模型及两种公式环境下的 DESeq2 调用，标准化结果表与原始结果表均 `identical()`。这个验证覆盖本次两组模拟设计，不等于所有设计和后端均已验收。轻量富集输入测量只验证其导出体积，未执行 ORA/GSEA 数值一致性测试。

原完整模型 bundle 的 `object.size()` 只有 8.68 MiB，而 future 测量为 23.69 MiB；仅看 `object.size()` 会漏掉环境携带的内容。两次保留模型调用的本机 elapsed 为 2.722 和 1.547 秒，含加载/缓存次序影响，不能据此评价公式修改的拟合加速效果。

**当前补丁完成到哪里**

- HEAD 已将 `run_diff(..., keep_model = FALSE)` 设为默认，在 `run-diff.R:295` 附近省略模型，仍可显式选择保留。这直接减少新结果的跨进程传递和保存负担；不减少 DESeq2 内部拟合步骤。
- 工作区已有 DESeq2 公式环境隔离修改，独立模拟结果支持其方向。
- App 富集入口仍原样传完整 bundle；项目恢复路径未发现针对 `diff_object` 的裁剪。因此恢复旧 `.omp`、已有会话中的大结果或显式 `keep_model = TRUE` 的结果仍可能重现问题。新默认不会追溯改变旧对象。
- 本机源码改变不代表运行中的容器已更新。需核实构建所用提交、镜像摘要及会话实际加载的 Core 函数；持久化恢复的数据也须纳入回归测试。

**差异分析为什么可能慢**

`diff-deseq2.R:92` 调用 `DESeq2::DESeq(dds, quiet = TRUE)`，未开启 `parallel = TRUE`；本机安装版本该参数默认 FALSE。App 的两个 future worker 表示任务并发容量，不会自动将一次 DESeq2 拟合拆成两核执行。仓库 ShinyProxy 模板目前限制每容器 2 CPU、6g 内存；Dockerfile 还将 OMP/BLAS/BiocParallel 的默认线程/进程数限制为 1。这些是源码配置，尚未验证服务器最终生效值。

DESeq2 拟合耗时随基因数、样本数和设计复杂度变化。官方建议可根据数据设计进行低计数预过滤，并谨慎评估并行收益；增加核心存在任务传递开销。预过滤与结果阶段的 independent filtering 是不同步骤，后者不能省掉已经发生的拟合计算。[DESeq2 官方说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#speed-up-and-parallelization-thoughts)

本次核查时工作区正在增加 raw-count QC 低表达过滤，默认 `min_count = 10, min_samples = 3`。这一规则能去掉部分低表达特征，但不能把 3 当作所有 200 例设计的通用最优值；阈值应结合实际分组/研究问题确定并记录。该 QC 改动的完整恢复、版本绑定和端到端效果不在本次验证范围内。[DESeq2 预过滤说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering)

此外，`wire_autosave()` 在项目状态改变后同步调用项目保存，保存内容含项目内 bundle。完整拟合对象会增加保存负担。富集模块还监听每次 diff bundle 改变并触发任务（源码注释说只自动一次，但 observer 没有 once 限制）。这些可能增加差异完成后的等待；未测服务器实际占比，不应认定它们就是主要耗时。

**建议交给 CC 的修复与验收顺序**

1. 在富集任务提交边界显式构造轻量输入：标准化的全基因结果表、必要的 input_info/comparison 等元信息和分析标识。保持 analysis_bundle 必需结构；上游来源与 session 准入检查继续使用原有绑定。不要把计数矩阵、模型、重复的原始结果表带入 worker。不能仅留下显著基因，以免改变 ORA 背景或 GSEA 排序。
2. 保留 Core 的可选完整模型能力，并完成公式环境隔离。App 即使收到旧结果或显式保留模型的结果，也应在提交边界得到同样的小输入；无需强迫用户重算已经有效的差异结果。
3. 扩展测试为真实 DESeq2 模型结果及旧项目恢复路径，检查 future 实际发现的 globals 和明确的体积预算。现有 `test-async-payload-sites.R` 的富集用例是 30 × 8 的 t-test，比较“闭包体积不超过输入 bundle + 常量”；它只能检测意外捕获更多环境，不能检测输入 bundle 本身不合理地大。新增测试还应比较轻量化前后的 ORA/GSEA 数值、背景集和完整排序，保留失效回调拒收的覆盖。
4. 为服务器加阶段时间：提交前扫描、构建 DDS、DESeq 拟合、results 提取、结果返回、保存、富集；同时记录有效基因数、样本数、设计列数和容器 CPU/内存。按测量结果决定预过滤和并行配置。将同时运行任务数与每任务核心数统一预算，避免直接把外层 future 与内层多进程都调大。
5. 在明确提交构建的新测试镜像中，分别验证新项目计算及恢复旧项目后的富集，并记录生效版本和实际 globals 体积。不要仅通过提高 `future.globals.maxSize` 或改用其他统计方法来宣布问题解决。

复现脚本与结果位于同目录 `evidence/rnaseq-payload-probe/`。`working.patch` 是复制源码时的工作区补丁快照，`snapshot-sha256.txt` 标识关键输入源码。脚本运行时需将工作目录设为包含 `packages/` 的对应源码副本；其恢复原公式环境的 mock 依赖快照中的 `detached_formula()` 修改。所有模拟数据均由脚本生成，证据目录不含用户样本数据。
