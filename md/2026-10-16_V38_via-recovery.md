# 2026-10-16：V38 补充 V35/V37 via 优点融合

- 状态：代码、完整说明、验证与修订 ZIP 已完成
- 基于版本：初次 V38（e8436bf），对照 V35 与 V37 via
- 目标版本：按用户指定修订 V38，不另建 V39；内部版本以 R2 结尾
- 目的：补足 V35 能识别、V37 仍漏检的 via，同时保留 V38 slit 修复及 V37 平滑椭圆/弱对象的优点。

## 开始状态与范围

- main → origin/main 工作区干净；`git fetch origin` 成功，`git merge --ff-only @{upstream}` 返回 Already up to date，数值扫描源码目录与 ZIP 最大版本均为 V38。
- 用户明确要求一起改进到 V38，故修订当前目录，完成后重建同名 ZIP。
- 使用前轮依赖齐备的 `/tmp/compare-v38-venv`（Python 3.12，requirements 不变）。用户真实 SEM 图片/TXT 仍不可用，本次证据来自受控合成图与真实程序流程。

## 复现和原因

构造 4×5 个半径 20 px 的圆孔，边缘 sigma=0.8，固定随机噪声 0.2 灰度；背景为倾斜平面叠加 A×sin(x/29)×cos(y/37)。

| 场景 | V35 | V37 | 修订 V38 |
| --- | --- | --- | --- |
| A=8，孔距 75 px | 20 | 6 | 20 |
| A=0，孔距 50 px | 20 | 0 | 20 |
| 含缺孔与星形污染的旋转长椭圆 | 少于 18 | 18 | 18 |

- **背景误判仍存在于 via**：V37 的大窗口背景平面与单点外探不能充分去除弯曲照明，出现 `UNRESOLVED_OUTER_BOUNDARY`，局部原灰度边缘正确的孔也被拒测。初次 V38 保留相同 via 实现，也继承该回退。
- **补检污染搜索尺度**：密集孔的局部对比候选包含较多孔间小碎片。V37 在合并后估计尺度，将真实约 39 px 的孔径变成约 13 px，后续边缘/格位/形状审核因错误尺度拒绝真实孔。
- 试验表明只把 slit 的横向背景响应直接用于 via 不足：A=8 只恢复 12/20。孔两侧的局部背景插值可恢复 20/20；密集孔还必须保护首轮尺度。两种问题分别处理，未直接放宽形貌门槛。

## 具体改动

- `SEM_0830_0831_compare_V38_bundle/sem_before_after_compare_0831_V38.py`：在同一 V38 入口中增加 via 修订，`V32_SCRIPT_VERSION=2026-10-16-V38-SLIT-VIA-FUSION-R2`；更新 source_algorithm 以表示 V35 首轮定位与 V37 整体候选融合，运行入口和输出文件名继续使用 V38。
- `_v38_via_boundary_response`：原始灰度经 sigma=0.8 平滑后，从孔 X 包络左右两侧灰度窄条逐行插值背景，允许弯曲照明；只在左右采样范围内使用响应，不以多项式外推制造前景。单侧样本不足时退回既有局部背景估计。该响应只用于定位/完整性审核，不是最终孔轮廓。
- via 复用 V38 的连续外侧前景采样，并保留 V37 有界区域重定位。最终仍经原图灰度 50% 法线追踪、支持率/连续缺边、形貌和格位审核，未改精细边缘追踪或形貌门槛。
- `_v32_via_search_scale` 增加仅用于尺度选择的 `prefer_round` 分支，沿用 V35 首轮完整近圆对象规则。显式尺寸提示优先；至少四个有效首轮对象可支撑时保留该尺度；近圆不足再使用 V37 首轮整体椭圆，仍不足才参考补检候选。近圆规则不限制最终椭圆测量。
- 径向定位使用首轮候选，保留 V35 的定位顺序；V37 弱对象/整体轮廓补检继续使用合并候选。补充碎片不替换已获支撑的首轮尺度。
- 新增 `via_boundary_rule`、`via_scale_locator`、`via_locator_path`，以及 `measure_v38_via_primary_scale_locked` / `primary_scale_candidate_count` / `primary_radial_count`（均带 measure_v38_via_ 前缀）；settings 保存尺度与背景规则。定位数量不等于最终统计数。
- `tests/test_v38_via_recovery.py`：8 项直接对照 V35/V37/V38 的 via 测试，核验孔数、真实中心、尺寸、尺度、补检、明亮极性、用户提示/筛选、弱下半部、空白测量原图拒绝、椭圆/异常形貌和独立 worker。
- `tests/test_v38_slit_fusion.py`：将原“via 逐浮点值相同”要求改为按图像格位匹配、核验已知真实中心及尺寸接近程度。新修订授权改变 via 边缘流程，不再要求复刻旧算法每个浮点数；保留对象数、真实位置和几何断言。
- `README_V38.txt` / `ALGORITHM_V38.md`：更新 V35/V37 via 融合、复现、参数/输出/诊断字段及限制说明；完整模式和统计说明保留，无新增依赖或 CLI 参数。
- 根 `README.md` / `md/README.md`：同步当前修订说明、对照结果与记录索引；在初次 V38 记录开头标明它是历史交付及旧 SHA，避免与当前同名 ZIP 混淆。
- AST 对比确认：前轮 slit 主测量、局部修复、背景与连续性采样，via 形貌/整体候选及原图精细追踪、trench 主测量实现均保持；共用包装只增加 via 背景分支。

## 验证

- 统一环境：`OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl`，使用 `/tmp/compare-v38-venv/bin/python`，未写缓存到发布目录。
- 初次继承验证 `python -m unittest discover -s tests -p 'test_v38*.py' -v`，当时发现 47 项：46 通过、1 失败，耗时 176.219 秒。失败项要求 V37/V38 via 浮点完全相同，且按浮点 X 首先排序，中心极小移动导致不同物理行错配。改为按实际格位排序，并增加已知中心断言；不是移除形状或数量检查。日志 `/tmp/compare-v38-via-inherited.log`。
- 修改后的保留能力测试：`PYTHONPATH=tests python -m unittest test_v38_slit_fusion.SlitFusionTests.test_whole_via_pipeline_preserves_v37_weak_and_elongated_objects -v`，1 项通过，34.073 秒。日志 `/tmp/compare-v38-via-preservation.log`。
- 新增 via 针对性测试：`python -m unittest discover -s tests -p 'test_v38_via_recovery.py' -v`，8 项全部通过，97.411 秒。日志 `/tmp/compare-v38-via-targeted.log`。
- 独立 worker：含中文/空格路径、PixelSize=2 的密集 via 实际测出 20 个，约 80×80 nm，首轮尺度约 39 px，标注存在，尺度/来源字段保留；真实图像读取及候选池流程已覆盖，不仅测试分割候选的直调。
- 弱下半部 via：6/12 灰度、明暗极性均测回约 45×77 px，采用原图边缘且不强制椭圆；原图替换为空白仍拒测；显式圆形度 1 会拒绝对应轮廓；空首轮候选仍可从补检恢复 20 个。
- 新修改 Python 文件语法、尾部空白、AST 保留范围、V38 `--help` 和本地 Markdown 链接检查通过；主脚本无历史版本导入，requirements 与 V37/初次 V38 一致。
- 最终回归：`PATH=/tmp/compare-v38-venv/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v38-mpl python3 -m unittest discover -s tests -v`，**260 项全部通过，668.085 秒，其中 55 项直接验证修订 V38**；日志 `/tmp/compare-v38-r2-regression.log`。最终 R2 版本标识、来源字段、尺度与背景分支均包含在此次回归中。
- 未验证事项：未取得本批实际 SEM 原图/TXT，无法量化真实样品召回率或误检率；不能保证 V35 的每个检测在任意输入中都保留。两侧连续污染、背景留白不足、极短弱延伸、错误首轮尺度及低信噪比仍可能影响判断；Windows/macOS 原生运行未验证。

## 交付

- ZIP：已重新生成 `SEM_0830_0831_compare_V38_bundle.zip`，205910 字节；按用户指定替换同名 V38 包。
- 明确发布清单保持四文件：主脚本、requirements、README、ALGORITHM，使用同名顶层目录，排除缓存、虚拟环境和临时数据。
- ZIP CRC 检查通过；实际解压后，4 个发布文件与版本目录逐字节一致，归档条目仅为明确发布清单。
- SHA256：`70a47e4778c8634e3e8d2f6b865c6f783cacde846cf8f6e0123e22601d80db71`。
- 最终 Markdown 本地链接、所有修改文本尾部空白与 `git diff --check` 检查通过。
- GitHub：git@github.com:Qiaochu-Zhang/compare.git，main → origin/main，源码、测试、说明、修改记录和修订 ZIP 一并提交推送；最终远端一致性结果见交付回复。
