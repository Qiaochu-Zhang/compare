# 2026-10-16：V38 slit 识别回退分析与分图案改进

- 状态：代码、完整说明、验证与 ZIP 已完成
- 基于版本：V37（SEM_0830_0831_compare_V37_bundle），对照 V35 slit 流程
- 目标版本：V38（SEM_0830_0831_compare_V38_bundle）
- 目的：分析 V37 相比 V35 的 slit 漏检原因，取长补短并保留 V37 的 via 改进。

## 开始状态与范围

- `git status --short --branch` 为干净的 main → origin/main；远端为既有 compare 仓库，沿用 SSH 配置。
- `git fetch origin` 成功，`git merge --ff-only @{upstream}` 返回 Already up to date；数值扫描源码目录和 ZIP，最大版本均为 V37，建立独立 V38。
- 用户确认暂无可用真实图片目录，先做代码分析及受控回归。未修改任何历史目录、ZIP 或历史测试。

## 已确认的原因

1. **V37 外侧前景误判**：大 ROI 只拟合背景平面，无法完全去除弯曲照明；低阈值单点探测还会碰到邻位或模糊尾部。局部原图边缘正确的 slit 被标为 `UNRESOLVED_OUTER_BOUNDARY`，或经不必要重定位扩大尺寸。
2. **V35 的失败轮廓局部修复被删除**：V35 验证失败后在同格位从原始灰度重新定位并再次验证，V36/V37 新逐对象循环缺少该分支。
3. **完整列被新增整图门槛误删**：V36/V37 为允许部分列设置的重复阵列门槛也作用于完整列。仅一列四排齐全、其他列缺失较多时，V35 可保留该列，V37 全部拒绝。
4. **候选影响定位的风险**：V37 将更多候选在四排定位前合并，候选替换/非阵列候选可能扰动行距、相位。V38 将首轮与补检独立并保护首轮可靠对象。本次不把该静态风险当成每张实际漏检的已确认原因。

受控合成对照（有效最下排 slit 数量）：

| 场景 | V35 | V37 | V38 |
| --- | --- | --- | --- |
| 弯曲背景 A=4 / A=8，实际约 17×21 px | 3 | 0 | 3 |
| 密集三列，宽约 65 px、列距 80 px | 3 | 0 | 3 |
| 一列完整、其他列缺失较多 | 1 | 0 | 1 |
| 模糊 sigma=3、真实宽约 17 px | 3，宽约 16.3–17.6 | 3，宽约 21.8–22.0 | 3，宽约 16.3–17.6 |

A=4 场景的额外原因隔离：仅给 V37 换用 V38 slit 完整性审核即可恢复 3 个；仅去掉补充候选但仍用旧完整性审核为 0 个。因此这张合成图的主要问题确为完整性误判，不能笼统归因于候选池变大。对照参数和真实跨度断言见 `tests/test_v38_slit_fusion.py`。

## 具体改动

- `SEM_0830_0831_compare_V38_bundle/sem_before_after_compare_0831_V38.py`：从 V37 建立独立入口 `main_v38` / `V38DetectionOptions`，更新脚本版本、默认输出、XLSX、CLI 和 source_algorithm；保留历史算法函数名、兼容字段及来源标记。
- 新增 `_v38_slot_boundary_response` / `_v38_slot_outer_contact`：slit 用逐行横向背景与连续外侧前景审核，跨越背景间隔的邻位不单独触发误判，探测越过常见模糊尾部；只用于审核/定位，最终尺寸仍由原图 50% 灰度法线追踪得到。
- 共用完整性函数仅 `pattern='slot'` 进入新分支。via 保持 V37 原有背景、补检、形貌、格位/去重、尺度和最终原图审核，trench 核心检测保持。
- 新增 `_v38_refine_slot_observations`：首轮有效对象优先；失败对象恢复同格位 V30 原图局部修复，限制中心偏移，重新验证真实边缘/形状/最终图像边界，补检前后均去重。
- slit 首轮和补检分开定位。首轮四排行级均已观测时用实际行级约束补检；偏离行级的候选记录为 `slot_supplement_outside_primary_rows`。
- 组列保留 V35 的有效完整列；部分列继续要求四排行级及重复阵列。证据不足时仅保留完整列，并记录部分列拒绝原因；无最下排不测、整排消失不提升第三排、不生成虚构轮廓或零尺寸。
- 新增 `slit_locator_path`、`slit_boundary_rule` 和 `measure_v38_slit_*` 计数、行级锁定及仅完整列标记；settings 保存新规则。审核计数不等于最终测量数。
- `README_V38.txt`：完整安装、目录/排序、single/compare、参数、六项指标、统计、输出及排障说明，更新完整列/部分列规则与诊断字段。
- `ALGORITHM_V38.md`：原因、V35/V37/V38 对照、逐步骤改动、保留 via 的范围和实际限制。
- `requirements.txt` 与 V37 完全一致，未增加依赖；发布包不导入历史版本。
- `tests/test_v38_slit_fusion.py`：9 项直接加载 V35/V37/V38 的回退、尺寸、明亮极性、比例尺/提示、弱下半部、实际短对象、坏种子原图修复、空白原图拒绝及 via 一致性测试。
- `tests/test_v38_edges.py` / `tests/test_v38_morphology.py` / `tests/test_v38_workflows.py`：将 V37 的 38 项场景实际应用于 V38，验证完整轮廓、形貌、污染缺失、trench、统计、导出、auto、真实 CLI/worker 和超时。
- 根 `README.md` 与 `md/README.md`：更新当前版本、使用方式、下载和记录索引；保留历史版本链接。

## 验证

- 环境：新建 `/tmp/compare-v38-venv`，Python 3.12，安装完整 requirements；numpy 2.5.4、pandas 3.0.6、scipy 1.18.1、opencv-python 5.0.0.93、matplotlib 3.11.2、openpyxl 3.1.5。沙箱联网失败后通过获准的依赖安装完成。
- 运行统一使用 `OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl`，避免缓存写入发布目录。
- 第一组针对性验证：`python -m unittest discover -s tests -p 'test_v38*.py' -v`（当时已创建 edges/morphology/workflows 三个文件），38 项全部通过，244.091 秒。日志 `/tmp/compare-v38-targeted.log`。
- 首次回退测试中，坏种子 fixture 只有四个点，而描述器要求至少五点，因而无法构造对象；修正测试为同一几何的六点轮廓，没有修改生产描述器去绕过输入要求。
- 修正后：`python -m unittest discover -s tests -p 'test_v38_slit_fusion.py' -v`，9 项全部通过，92.595 秒。日志 `/tmp/compare-v38-fusion-final.log`。后续仅完善原图审核日志和完整列场景拒绝记录，进入最终全仓库回归。
- AST 对比：V37/V38 的 via 主测量、整体候选、形貌指标、形貌 QC、可选质量筛选、局部多尺度候选、原始精细边缘追踪，以及 trench 主测量/候选函数实现一致。共用完整性函数的 slit 分支明确隔离；V37/V38 弱孔阵列 20 个、旋转椭圆阵列 18 个的中心和尺寸逐对象一致。
- 弱下半部 slit：6/12 灰度的下半部在明暗两种极性下测回约 17×71 px；实际短对象维持约 35 px；窄坏种子在合成原始灰度中恢复约 17×21 px，测量原图换成空白时拒绝。
- 所有新 Python 文件语法解析、尾部空白、V38 `--help` 检查通过。生成最终 ZIP 后，根 README、记录索引及本次记录的本地 Markdown 链接均已检查有效。
- 最终全回归：`PATH=/tmp/compare-v38-venv/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl python3 -m unittest discover -s tests -v`。**252 项全部通过，耗时 592.748 秒，其中 47 项直接验证 V38**。日志 `/tmp/compare-v38-regression.log`。最终源码中的审核日志、完整列保留及部分列拒绝记录均已纳入此次回归。
- 未验证事项：本批真实 SEM 原图及 TXT 不可用，不能量化真实样品召回率/误检率或宣称漏检全部解决。超大模糊、极短弱尾部、连续污染、极低对比及未知构图仍可能影响可靠性。单独一个几乎无阵列参照的完整列可能在上游定位时失败。Windows/macOS 原生运行未验证。

## 交付

- ZIP：`SEM_0830_0831_compare_V38_bundle.zip`，202629 字节。
- 发布清单：独立主脚本、requirements、README、ALGORITHM，4 文件；同名顶层目录，排除缓存、虚拟环境、密钥、临时输出及测试数据。
- ZIP CRC 完整性检查通过；实际解压后，4 个发布文件与源码目录逐字节一致，归档没有其他条目。
- SHA256：`fff63437c687903476776557ea7fc37589c2a937445bfbedd8a111439e87906d`。
- 发布后本地 Markdown 链接、所有新增文件尾部空白及 `git diff --check` 检查通过。
- GitHub 目标：`git@github.com:Qiaochu-Zhang/compare.git`，main → origin/main；本次源码、测试、说明、修改记录与 ZIP 一并提交推送，最终远端核对结果在交付回复中说明。
