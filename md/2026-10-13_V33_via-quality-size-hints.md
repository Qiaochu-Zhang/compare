# 2026-10-13：via 质量筛选与尺寸辅助识别

- 状态：已完成
- 基于版本：V32（`SEM_0830_0831_compare_V32_bundle`）
- 目标版本：V33（`SEM_0830_0831_compare_V33_bundle`）
- 目的：提供 via 圆形度阈值及同图尺寸偏差阈值，并支持用户预先指定 slit、via、trench 的大致尺寸，改善模糊图片的搜索与补检。
- 用户确认：尺寸偏差相对同图 via，非 before/after 差异；尺寸单位 nm，可选填 via 直径、slit 宽度/长度、trench 宽度/间隙。

## 具体改动

- `SEM_0830_0831_compare_V33_bundle/sem_before_after_compare_0831_V33.py`：在 V32 独立脚本基础上新增 `V33DetectionOptions` 及 CLI 参数，默认 via 最低圆形度 0.75、同图等效直径最大偏差 0.30。先通过原图边缘、圆形度与所填尺寸范围检查，再从至少 3 个同图有效孔计算一次中位直径；筛掉异常后重新拟合坐标轴、编号并输出，前后阶段独立筛选。
- 同一脚本：新增可选的 `--via-diameter-nm`、`--slit-width-nm`、`--slit-length-nm`、`--trench-width-nm`、`--trench-gap-nm` 和 `--size-tolerance`；默认尺寸容差 ±50%。尺寸提示参与 via 多尺度搜索、slit 去背景候选和缺失位置搜索、trench 边缘窗口和成对端点选择，并按每张图的 PixelSize 换算。slit/slot 参数互为别名，长度范围只检查最下排。
- 同一脚本：参数随独立 worker 请求传递，并保存到 `settings.json`；对象和剔除日志补充 `input_*`、via 圆形度/等效直径/中位直径/偏差及原因码。补检不降低用户圆形度阈值，尺寸提示不替代原图测量；保留 300 秒硬超时、四排完整列和历史兼容字段。
- `tests/test_v33.py`：将继承回归测试指向 V33，覆盖灰度边缘精度、物理比例尺、四排 slit、真实 worker、超时和继续处理；版本及输出文件名预期更新到 V33。
- `tests/test_v33_controls.py`：新增 14 项用户参数、形状/尺寸筛选、模糊 slit 补检、空白拒测、参数隔离及进程传递测试。
- `SEM_0830_0831_compare_V33_bundle/README_V33.txt`、`ALGORITHM_V33.md`：新增参数说明、公式、使用示例、尺寸适用范围、输出字段和验证限制。运行依赖沿用 `requirements.txt`。
- `README.md` 与 `md/README.md`：更新当前版本、下载入口、参数表和本次修改记录；保留 V32 及更早版本目录与 ZIP。

## 验证

- 全量回归：`MPLBACKEND=Agg MPLCONFIGDIR=/tmp/compare-v33-mpl PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 /tmp/compare-v33-venv/bin/python -m unittest discover -s tests -v`，116 项全部通过（205.276 秒）。测试依赖安装在 `/tmp` 隔离环境，未加入版本包。
- 合成 via 阵列含 19 个直径 40 px 的孔和 1 个 56 px 的孔：关闭同图尺寸筛选保留 20 个，阈值 0.20 时保留 19 个，剔除日志、有效标注和统计数量一致。
- 带梯度和噪声的 slit 图额外高斯模糊 3.5 px：不提供尺寸时无有效列，提供 17 × 21 nm 的最下排尺寸后恢复 2 列，仍通过四排完整性及原图边缘支持验证。
- ZIP 解压后的独立 CLI 处理 6 张合成图：退出码 0，6 张状态均 OK，64 条有效对象记录，0 条处理错误，Excel、CSV 和标注均生成。每阶段得到 20 个 via、3 个 slit 最下排对象及 9 组 trench。
- 上述独立运行统一设置 via 提示直径 90 nm；两阶段使用不同 PixelSize，实测平均宽度分别约 79.98 nm、88.00 nm，提示未覆盖实际差异。
- 文档本地链接检查、最新版本数值识别及 Git 空白格式检查通过。
- 本仓库没有本批真实 SEM 图片及配套 PixelSize TXT，需明确合成图验证的范围。

## 交付

- ZIP：`SEM_0830_0831_compare_V33_bundle.zip`，179116 字节；只包含同名顶层目录下的主脚本、requirements、README 与算法文档。CRC 检查通过，4 个发布文件与源码目录逐字节一致。
- ZIP SHA-256：`ccb0654d7ea10e09577a863fc9a4ab9ce7f8a42e2a8e2284c95b3af77d7998de`。
- GitHub 同步目标：`Qiaochu-Zhang/compare` 的 `main` 分支。
