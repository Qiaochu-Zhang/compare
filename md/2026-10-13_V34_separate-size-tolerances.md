# 2026-10-13：宽度与长度 / 间隙默认容差分离

- 状态：代码、文档、验证与 ZIP 已完成
- 基于版本：V33（SEM_0830_0831_compare_V33_bundle）
- 目标版本：V34
- 目的：宽度与 via 直径默认 ±50%，slit 长度与 trench 间隙默认 ±100%。

## 具体改动

- `SEM_0830_0831_compare_V34_bundle/sem_before_after_compare_0831_V34.py`：新增 `width_tolerance` / `length_tolerance`，应用于候选搜索及最终尺寸检查；旧 `size_tolerance` 可统一覆盖默认，独立参数优先。宽度 / 直径仍使用 50%，长度 / 间隙使用 100%；实测尺寸必须为有限正数，trench 两端点仍需正确排列。配置和结果记录实际生效的两类容差。
- 同目录 `README_V34.txt` / `ALGORITHM_V34.md`：参数、兼容行为、正尺寸边界、使用示例与验证范围；运行入口、输出文件及标识升级至 V34。历史算法 helper 名称及兼容字段保留。
- `tests/test_v34_tolerances.py`：实际加载 V34，验证默认值、自定义值、旧参数优先级、序列化、参数校验、有效尺寸边界、slit / trench 实测、via / slit 宽度独立筛选、空白拒测、真实子进程及配置导出。
- `README.md`、`md/README.md`：当前版本入口、参数表、下载链接及修改记录索引。
- 历史版本目录和 ZIP 保留原样。

## 验证

- 针对性测试：`MPLBACKEND=Agg python -m unittest discover -s tests -p test_v34_tolerances.py -v`，10 项通过（6.875 秒）。
- 仓库完整回归：`MPLBACKEND=Agg python -m unittest discover -s tests -v`，126 项通过（192.523 秒）。使用现有 `/tmp/compare-v33-venv/bin/python`，设置 `PYTHONDONTWRITEBYTECODE=1`、`OPENBLAS_NUM_THREADS=1`、`MPLCONFIGDIR=/tmp/compare-v33-mpl`；无依赖变更。
- slit：实测长度约 21 nm，分别提示 12 / 50 nm；默认 100% 容差均保留 3 个最下排对象（共 12 个四排轮廓），收紧到 50% 时均剔除。
- trench：实测间隙约 60 nm，分别提示 35 / 140 nm；默认 100% 均保留 9 组，收紧到 50% 时均剔除。
- 真实 worker：PixelSize=2 nm/pixel、提示间隙 70 nm、宽度容差 0.25 / 长度容差 0.8，保留 9 组，实测间隙约 120 nm；结果记录两类有效容差，标注文件存在。
- 新脚本 `--help`、Python 语法、文档链接及 `git diff --check` 检查通过。
- 未验证事项：仓库无本批真实 SEM 原图和比例尺文件，不能声称实际图片识别准确率已验证。

## 交付

- ZIP：`SEM_0830_0831_compare_V34_bundle.zip`，180811 字节。
- 明确清单包含主脚本、requirements、README、ALGORITHM 共 4 文件，位于同名顶层目录；CRC 检查与解压后逐字节对比均通过，无虚拟环境、缓存、密钥或临时输出。
- ZIP SHA256：`fa551184d7e4ac8e076611d248a1c39d22bc40898f6703655864737e0ad874cd`。
- GitHub 同步目标：`git@github.com:Qiaochu-Zhang/compare.git`，分支 `main`；源码、测试、文档、本记录及 ZIP 一并提交，使用现有 SSH 配置同步。
