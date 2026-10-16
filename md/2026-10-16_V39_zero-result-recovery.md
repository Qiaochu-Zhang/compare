# 2026-10-16：V39 零结果识别补检

- 状态：代码、说明、完整回归及 ZIP 交付校验已完成
- 基于版本：V38 R2（73667d1，`SEM_0830_0831_compare_V38_bundle`）
- 目标版本：V39（`SEM_0830_0831_compare_V39_bundle`）
- 目的：增强整张图无有效测量时的定位与结构识别，保留 V38 slit/via 融合能力。

## 开始状态与原因

- main → origin/main 工作区干净；fetch 和快进检查成功，源码目录及 ZIP 的最大整数版本均为 V38。
- 原零结果补检仅包含两级对比度增强/定位阈值调整，不自动尝试相反极性。
- 已用真实 V38 处理入口复现：反向极性的 via/slit/trench 均为 0；清晰 2×2 via 也为 0。后者进入阵列选择的最低支撑数为 5，与下游允许 4 个支撑的模型校验不一致。
- 原 V37 局部搜索的整图 99.5 分位峰值估计对极小占比目标可能失效；新增候选搜索改在超过噪声下限的正响应内估计阈值。
- 本批真实 SEM 原图/TXT 尚不可用，继续使用代码分析、已知几何合成图及实际 CLI/worker 回归。

## 具体改动

- 从 V38 R2 的明确四文件清单建立 V39 独立目录，保留历史目录与 ZIP。
- 更新 V39 主脚本入口、选项类、输出名与版本标识；历史方法/字段按实际来源保留。
- 新增 `_v39_zero_candidates`：原图多尺度定位，slit 横向核/via 方形核，最大 1025 px，加入显式尺寸提示尺度，基于噪声以上响应估计阈值。最多检查 240 个去重定位轮廓、保留 120 个原图初审有效候选。响应掩膜不作为最终测量轮廓。
- 候选先接受现有原图边缘/完整性和形貌初审，再进入原主测量与最终用户参数审核；避免大量反极性背景碎片干扰孔径/阵列拟合。原型试验中未做初审的反向 via 仍为 0，加入初审后恢复 20 个。
- 原流程、L1、L2 均为零时才执行 L3 原极性多尺度/结构试验；仍为零再 L4 相反极性。trench 不重复已穷尽的同极性搜索，直接尝试 L4。新增试验保持严格模型上下文，所有临时极性、门槛和四孔启用标记均用 finally 恢复。
- `_choose_via_lattice_v21` 的正常五点要求保留；仅在新试验中接受完整 2×2 四孔，须夹角 80°–100°、四格唯一、残差不超过 min(2.5 px,0.08×间距)。最终原图、形貌、用户范围或跨度审核少一个孔，整组四孔启动结果拒测。
- 在真实单图处理函数中接入新试验，保持首个有效结果优先、原始图像标注及既有总时限。新增方法/极性/候选数审计、对象复核标记、L3/L4 状态与补检横幅；settings 记录规则和上限。原 measurement_confidence_flag 保留补检复核兼容值，新增 confidence_flag 区分 V39 搜索恢复。
- `tests/test_v39_recovery.py`：14 项直接加载 V39 的零结果、误检边界及 V38 保留能力验证。
- `tests/test_v39_workflows.py`：14 项实际加载 V39 的 single/compare、统计/导出、auto、CLI/worker 及超时测试。
- 本版 README/ALGORITHM、根 README 及 md 索引同步新入口、机制、对照结果、审计字段和已知限制；依赖与 V38 完全相同。

## 受控对照

| 场景 | V38 R2 | V39 |
| --- | --- | --- |
| 反向极性 via，含弱孔及弱下半部 | 0 | 20 |
| 反向极性 slit，三列四排 | 0 | 3 个最下排 |
| 反向极性 trench，弱下条 | 0 | 9 对 |
| 清晰 2×2 via，直径约 40 px | 0 | 4 |

- 对照运行真实同名 PNG/TXT、正常及 L1/L2 流程，结果见 `/tmp/compare-v39-probe-refined.log` 及自动回归日志。PixelSize=2 的输出按真实边缘约为 via 80×80 nm、slit 34×42 nm、trench 宽 40 nm/间隙 120 nm。
- 800×800 画布内四个约 14 px 直径小孔可从新候选搜索建立四点支撑；显式提示保持有效。三孔、单排、非矩形、星形异常、用户过滤后仅三孔、空白/斜坡/噪声均未强行产生结果。
- V38 弯曲/密集 via、旋转平滑椭圆、弯曲/密集与局部缺失 slit 数量继续保持；已有正常结果时，Mock 断言新增恢复函数未被调用。

## 验证

- 统一环境：`PATH=/tmp/compare-v38-venv/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl`，Python 3.12，未增加依赖或写缓存到发布目录。
- `python3 -m unittest discover -s tests -p 'test_v39*.py' -v`：初次 26 项全部通过，163.128 秒，日志 `/tmp/compare-v39-targeted.log`。
- 后补“仅一孔被用户筛掉”和“显式 bright 默认恢复”等验证，并加强异常恢复断言；`PYTHONPATH=tests python3 -m unittest test_v39_recovery.ZeroResultRecoveryTests.test_filtering_one_of_four_valid_contours_invalidates_sparse_bootstrap test_v39_recovery.ZeroResultRecoveryTests.test_opposite_trial_respects_an_explicit_bright_default_and_restores_it test_v39_recovery.ZeroResultRecoveryTests.test_failure_restores_polarity_sparse_gate_and_threshold_context -v`：3 项通过，2.522 秒，日志 `/tmp/compare-v39-context-targeted.log`。与初次重叠一项，共 28 项 V39 直接验证通过。
- 语法 AST、V39 `--help`、尾部空白检查通过；requirements 与 V38 逐字节相同。交付前 Markdown 本地链接全部有效，`git diff --check` 和已暂存内容空白检查通过。
- AST 逐定义对比（包含同名历史重定义）确认 12 项关键原图/形貌/trench 函数完全保持：V38 slit/via 局部背景、连续性、失败种子修复、原图精细追踪、via 形貌与整体候选、trench 搜索和主测量。slit 主测量仅更新终端版本文案。
- 额外真实 worker auto 检查：传入 `v25_pattern_key=None`，共享 120 秒时限，同一四孔图的 via40/via60 均以 L3 测到 4 个，复用 via60 胜出结果和标注，PixelSize=2 的宽高均值约 79.35 nm，父进程极性/四孔开关保持默认；日志 `/tmp/compare-v39-auto-probe.log`。
- 首轮全仓库回归在审查中发现 V39 的几条终端显示仍写 V38 后主动中止（exit 130），当时未观察到测试失败；日志保留为 `/tmp/compare-v39-regression-interrupted.log`。修正当前版本的批处理横幅、进度、via 异常及 slit 进度文案，并在真实 CLI 测试增加版本显示/settings 断言；算法行为未改。
- 修正后最终全仓库回归：`PATH=/tmp/compare-v38-venv/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl python3 -m unittest discover -s tests -v`，**288 项全部通过，844.127 秒，其中 28 项直接验证 V39**；日志 `/tmp/compare-v39-regression.log`。最终版本显示和实际 CLI/settings 断言均包含在此次回归。
- 未验证事项：真实 SEM 召回率/误检率及 Windows/macOS 原生运行。

## 交付

- ZIP：`SEM_0830_0831_compare_V39_bundle.zip`，209496 字节，同名顶层目录。
- 明确发布四文件：`sem_before_after_compare_0831_V39.py`、`requirements.txt`、`README_V39.txt`、`ALGORITHM_V39.md`；不包含缓存、虚拟环境、临时输出、测试数据或无关文件。
- ZIP CRC 完整性检查通过；归档条目与发布清单一致；实际解压后四个文件与源码目录逐字节一致。校验报告：`/tmp/compare-v39-package-check.json`。
- SHA256：`4c60c5557feea27e3cb4d02f9da3ee56b406175866c66378589570072f1b6fe9`。
- GitHub 同步目标：git@github.com:Qiaochu-Zhang/compare.git，main → origin/main；发布提交包含本记录、V39 源码/说明/测试、根 README、索引和 ZIP。交付回复提供提交 ID 及推送后远端一致性核对结果。
