# 2026-10-15：V35 单文件夹样品分析

- 状态：代码、完整说明、验证与 ZIP 已完成
- 基于版本：V34（SEM_0830_0831_compare_V34_bundle）
- 目标版本：V35（SEM_0830_0831_compare_V35_bundle）
- 目的：用户只提供一个根目录，每个数字子目录作为独立编号样品，沿用现有测量及统计作图规则；保留 before/after 模式。

## 具体改动

- `SEM_0830_0831_compare_V35_bundle/sem_before_after_compare_0831_V35.py`：独立 V35 脚本；增加 `--input`、`--mode single|compare` 与 `main_v35`。单输入自动选择 single，参数冲突在创建输出前报错。数字目录按数值排序并保留前导零；每个样品的图片数独立，沿用自然排序、trench/slot/via 及 Region 2 起始规则。空样品记警告，比例尺错误保留排序位置并跳过对应图片。
- 单模式使用中性 `stage=single`，共用 V34 原图边缘、四排 slot 最下排测量、via 圆形度/同图尺寸筛选、尺寸提示与 300 秒硬时限。auto 对每张 via 独立选模型，复用最优测量/标注；比较模式继续由 before 选模型并固定给对应 after。单模式采用通用检测分支，不触发历史 AFTER 专用增强。
- 新增 `sample_metric_summary.csv`、单模式 XLSX 和三类图：单样品各区域、全部样品区域总览、合并区域后的样品汇总。六项指标、nm、均值柱、样本 SD（ddof=1）误差棒与 220 dpi 沿用原规则；汇总按对象计算，n=1 的 SD/SEM 留空，无有效测量不补 0。单模式不构造前后配对表和差值图。
- CSV/Excel 共用表清单，单模式空结果仍保留可读取表头；settings 记录模式、输入、数量/排序/聚合/作图规则。输出目录须不存在或为空，避免跨模式或重跑时混入旧图。
- 同目录 `README_V35.txt`：完整中文说明，包含安装、Windows/macOS/Linux 运行、目录树、排序与比例尺、两种模式、全部 CLI 参数、SD/SEM 与汇总权重、输出字段、故障处理及验证限制。`ALGORITHM_V35.md`：模式分支、数据流、统计和兼容说明。`requirements.txt` 保持 V34 依赖清单。
- `tests/test_v35_single.py`：14 项直接加载 V35 的测试，覆盖单模式目录、参数、统计/图像、错误恢复、导出、真实 CLI/worker 及 compare 兼容。
- `README.md`、`md/README.md`：更新当前版本、下载链接、两种模式的快速使用方式及修改记录索引。
- 历史版本目录、ZIP 与历史测试保持不变。

## 验证

- 开始时 `git status --short --branch` 确认工作区干净，`git fetch origin` 后快进检查为 Already up to date；扫描数字版本确定最高源码/ZIP 为 V34。
- 测试环境：`/tmp/compare-v35-venv`，Python 3.12；安装包内 requirements。numpy 2.5.4、pandas 3.0.6、scipy 1.18.1、opencv-python 5.0.0.93、matplotlib 3.11.2、openpyxl 3.1.5。
- 初次 V35 针对性测试中 13 项通过、1 项测试预期错误：同主文件名不同扩展名的大小写次序预期不符。改为无重复主文件名的混合格式测试集，明确验证 `1,2,3,10,11` 自然排序；没有为适配测试改变原排序算法。
- 完整回归：`PATH=/tmp/compare-v35-venv/bin:$PATH PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v35-mpl OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s tests -v`。
  **140 项全部通过，耗时 215.182 秒**，包含原有 126 项和新增 V35 14 项。日志位于本机 `/tmp/compare-v35-regression.log`，不入发布包。
- 真实 CLI/子进程：中文及带空格目录，样品 007 有 3 张 trench/slot/via，PixelSize=2；样品 10 有 1 张 trench，PixelSize=1。实际 auto 选模；41 个有效对象、4 张图的标注路径存在、20 张 PNG、8 行样品指标及完整 XLSX/CSV。
- 尺寸验证：007 的 trench 宽度/间隙接近 40/120 nm，slot 最下排 34/42 nm，via 80/80 nm；10 的 trench 为 20/60 nm。全部对象记录独立宽度/长度容差 0.25/0.8。
- 统计验证：`[10,14]` 和 `[30]` 两个区域合并后均值为 18（而非各图等权的 21）；实际 matplotlib 柱高、误差线、图轴数值排序和 220 dpi 参数通过断言；n=1 的 SD 保持 NaN。
- 故障验证：缺 TXT/无效 PixelSize、空样品、超时、损坏图片不影响后续有效图片；真实空白图不产生 0 尺寸；真实 0.001 秒硬时限终止 worker 且不发布部分标注。
- compare 回归：数量不等条件被跳过，三种图案共 6 项指标保留 after−before=1 的测试差值，比较工作表保留，auto 模型冻结行为通过。
- `--help`、Python 语法解析、文档本地链接及新增文件尾部空白检查通过；`git diff --check` 通过。全回归后仅完善说明及两条错误提示的 V35 文案，未修改测量/执行逻辑。
- 未验证事项：本批真实 SEM 原图及 TXT 不在仓库；目前为合成图与流程验证，不能声称实际样品识别准确率已验证。Linux/Python 3.12 已运行；Windows/macOS 的原生运行未在本机验证。

## 交付

- ZIP：`SEM_0830_0831_compare_V35_bundle.zip`，190672 字节。
- 明确发布清单：主脚本、requirements、README、ALGORITHM，共 4 个文件，位于同名顶层目录。无虚拟环境、缓存、临时输出、测试数据或凭据。
- ZIP 完整性（CRC）检查通过；实际解压后每个文件与源码目录逐字节一致。
- SHA256：`aede64225386df65788c34836ff7f7b108dd6a791098eee68fdf8101b0fbdae7`。
- GitHub 同步目标：`git@github.com:Qiaochu-Zhang/compare.git`，当前 `main` → `origin/main`；源码、测试、说明、根 README、本记录与 ZIP 一并提交，沿用已有 SSH 配置推送。
