# 2026-10-16：V41 slit 四排长度递增与 100 nm 下限

- 状态：已完成
- 基于版本：V40（6411236，`SEM_0830_0831_compare_V40_bundle`）
- 目标版本：V41（`SEM_0830_0831_compare_V41_bundle`）
- 目的：slit 按图像从上到下的四排，Y 向实测长度须递增；用户已确认仅最下排 Y 长度至少 100 nm，上三排可以短于 100 nm。

## 具体改动

- 从 V40 四个发布文件建立 V41 独立目录，保留历史源码和 ZIP；更新入口 `main_v41`、`V41DetectionOptions`、默认输出/工作簿名、script_version/source_algorithm 和当前版本显示。
- 用户已确认：100 nm 下限只作用于第 4 排，上三排可以短于 100 nm。四排按图像上下位置确定，实际 Y 长度仍须向下严格递增。
- `_v41_slot_length_audit`：使用原图细化后轮廓沿实际 Y 轴投影，乘每图 PixelSize 得到 nm；检查实际中心行序和可见长度递增，只有最下排执行 100 nm 下限。部分列不补造缺失行长度。
- `_v41_filter_slot_length_columns`：违规整列在计算阵列支撑前剔除，重拟合测量轴并重新检查直到稳定；四个实际行级的长度中位数也须递增。切换到仅完整列时再次审核，最后仍执行用户尺寸提示与容差。
- 发现相反极性补检可将原本约 98 nm 孔槽外的较长亮边条带当成 slit。新增 `_v41_exclude_slit_rim_bands`，将已验证亮圈附近属于同一对象的普通候选关联到实际内沿；该内沿仍须通过原图与长度审核，避免通过测量亮边材料绕过孔内尺寸下限。
- 所有正常/L1/L2/L3/L4 slit 路径统一应用新规则。不会按长度重排行号、拉长轮廓、固定输出 100 nm 或为缺失对象填零。
- 新增 `v41_slit_length_rule` settings、逐列实际长度 JSON、规则通过标记、下限/倒置/全图行序/亮边条带剔除审计。原 V40 亮边与 V39 恢复字段保留来源含义。
- `tests/test_v41_slit_rules.py`：19 项直接加载 V41 的长度规则、亮边/暗槽、缺失、单位换算、明确下限、补检拒测、独立 worker 和真实 compare 测试；短孔的完整补检分别覆盖 3 px / 7 px 亮边。
- `tests/test_v41_workflows.py`：14 项实际加载 V41 的两模式、统计/导出/作图、auto、CLI/worker 和超时测试。真实 CLI 的 slit 输入改为 21/31/43/55 px，PixelSize=2 对应约 42/62/86/110 nm，专门验证“上三排可短于 100 nm”。
- V41 README/算法说明、根 README、md 索引与记录同步；requirements 与 V40 逐字节相同，无新增依赖。via/trench 的主测量、原有原图边缘与亮圈提取保持。

## 验证

- 开始时 main → origin/main 干净，fetch 与 ff-only 成功，最大源码和 ZIP 版本均为 V40。
- 原 V40 的 55/43/31/21 px 递减 slit 保留为历史样例；新版本使用递增正例，并将递减作为拒测反例，历史源码/测试不改。
- 环境沿用 `/tmp/compare-v40-venv`，Python 3.12；统一 `PATH=/tmp/compare-v40-venv/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v40-mpl`。
- 首轮 V41 测试 33 项中 29 通过、3 失败、1 错误，126.138 秒：两项直接比较浮点数的断言改为数值精度比较；compare 无错误时继承的错误 CSV 可以仅含 BOM，因此改检查实际 Excel 的 errors 表；一项实际发现亮边材料可被反极性补检误测，已通过亮圈内沿关联修复。日志 `/tmp/compare-v41-targeted.log`。
- 修复后相关 4 项通过，49.260 秒，日志 `/tmp/compare-v41-fixes.log`；随后扩大边缘带关联范围并补充 7 px 亮边反例，纳入以下最终专项验证。
- 最终 V41 专项：`python3 -m unittest discover -s tests -p 'test_v41*.py' -v`，**33 项全部通过，155.039 秒**；日志 `/tmp/compare-v41-targeted-final.log`。
- AST/语法与 requirements 一致性检查通过；按所有同名重定义比较确认 11 项关键算法未改：普通/完整原图追踪、V40 亮圈追踪/定位/格位/补入、全部 via/trench 主测量及 via 形貌。V41 `--help` 已运行，明确注明递增和仅最下排至少 100 nm。
- 最终全仓库回归：`MPLBACKEND=Agg python3 -m unittest discover -s tests -v`（上述统一环境），**369 项全部通过，1348.641 秒，退出码 0**；日志 `/tmp/compare-v41-regression.log`。其中 33 项直接加载 V41，覆盖正常测量、全部补检、真实 worker/CLI、single/compare 和导出。
- 发布包使用明确四文件清单，CRC 完整性检查通过；实际解压后四个文件均与版本目录逐字节一致，报告 `/tmp/compare-v41-package-check.json`。未包含缓存、虚拟环境、临时输出或历史数据。
- `git diff --check`、任务文件空白/换行、Python AST、文档本地链接检查均通过；最大同系列源码及 ZIP 版本为 V41。
- 提交前再次 `git fetch origin` 成功，`git rev-list --left-right --count HEAD...@{upstream}` 为 `0 0`，没有待整合的远端变更。
- 未验证：本批真实 SEM 原图/TXT、实际召回率/误检率和 Windows/macOS 原生环境。

## 交付

- ZIP：[SEM_0830_0831_compare_V41_bundle.zip](../SEM_0830_0831_compare_V41_bundle.zip)，215202 字节。
- SHA-256：`1f705dd3ca478ee28f0dc876a955fff1e2e88d9e93c2c4ba773c74942a297430`。
- ZIP 内保留 `SEM_0830_0831_compare_V41_bundle/` 顶层目录，发布文件为 `sem_before_after_compare_0831_V41.py`、`requirements.txt`、`README_V41.txt`、`ALGORITHM_V41.md`。
- GitHub 同步目标：git@github.com:Qiaochu-Zhang/compare.git，main → origin/main；本记录与源码、测试、说明和 ZIP 一并提交，实际提交号与远端核对结果见本次交付回复。
