# SEM 样品分析与 before/after 对比

当前版本：[V39](SEM_0830_0831_compare_V39_bundle/README_V39.txt)，增强**整张图零结果时的识别**，继承 V38 R2 的 slit/via 融合和两种运行模式。

下载 [V39 完整压缩包](SEM_0830_0831_compare_V39_bundle.zip)。解压进入同名目录，安装依赖：

```bash
python -m pip install -r requirements.txt
```

## V39 零结果补检

- 原流程及既有 L1/L2 补检均没有有效测量后，再做原极性的多尺度原图搜索；已有结果时沿用原流程。
- 局部响应阈值由噪声以上的有效信号估计，减少小图案占比低时被整图分位数忽略的问题。候选先通过原图边缘与形貌初审，再参与阵列定位。
- 仍为零时自动尝试相反明暗极性，覆盖明亮目标；每次试验后恢复配置，单图审计记录实际极性与所用方法。
- via 可从完整、近直角的 2×2 四孔阵列启动；四个孔须全部通过最终原图、形貌和用户尺寸审核。三孔、单排或缺边无法凑成完整二维支撑。
- 新增补检仍执行原图测量及显式参数，沿用单图 300 秒硬时限。补检结果有复核标记，缺失尺寸不补零。

受控合成图的完整处理流程对照：

| 场景 | V38 R2 有效数 | V39 有效数 |
| --- | --- | --- |
| 反向极性 via，20 个孔 | 0 | 20 |
| 反向极性 slit，3 列四排 | 0 | 3 个最下排 |
| 反向极性 trench，9 对 | 0 | 9 |
| 清晰 2×2 via | 0 | 4 |

V38 对弯曲背景、密集邻位、弱孔、椭圆及部分缺失 slit 的识别继续验证。没有本批真实 SEM 原图/TXT，以上不能当作真实召回率。实现和限制见 [V39 算法说明](SEM_0830_0831_compare_V39_bundle/ALGORITHM_V39.md)。

## 单文件夹模式

只提供样品根目录，其下每个纯数字子目录代表一个样品：

```text
samples/
  001/
    image_1.tif + image_1.txt   # Region 2 trench
    image_2.tif + image_2.txt   # Region 2 slot/slit
    image_3.tif + image_3.txt   # Region 2 via
    image_4.tif + image_4.txt   # Region 3 trench，继续同样顺序
  002/
    image_1.tif + image_1.txt
    image_2.tif + image_2.txt
    image_3.tif + image_3.txt
```

```bash
python sem_before_after_compare_0831_V39.py --input "D:/SEM/samples" --output "D:/SEM/single_V39"
```

`--input` 自动选择单文件夹模式，也可显式加 `--mode single`。各样品图片数可以不同；编号保留前导零，按数值排序。

图片仍按文件名自然排序，每三张依次分配 **trench、slot、via**，Region 从 2 开始；类型由排序位置决定。尾部可以不满三张，中间缺图会导致后续类型错位。支持 TIF/TIFF/PNG/JPG/JPEG，每张图片须有同名 TXT，包含 `PixelSize=正数`，单位 nm/pixel。只扫描直接数字子目录内的图片。

统计作图继续使用六项尺寸指标、nm 单位、均值柱和样本标准差误差棒；slot 只统计边缘有效的完整列最下排，或满足阵列支持要求的部分列最下排。生成单样品各区域图、全部样品区域总览图、合并区域后的样品汇总图。汇总按有效对象计算，n=1 的标准差为空，无有效值不填 0。

主要结果：`SEM_0830_0831_single_V39_results.xlsx`、`sample_metric_summary.csv`、单图统计与对象明细、`annotated/`、`plots/`。单文件夹模式不输出前后配对和差值文件。

## before/after 对比（保留）

```bash
python sem_before_after_compare_0831_V39.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V39"
```

可显式加 `--mode compare`。从 after 发现数字目录，与 before 完全同名目录匹配；对应目录图片数量须相等，按自然排序配对，文件名和格式可以不同。输出 `SEM_0830_0831_before_after_V39_results.xlsx`，保留图片级、总体及逐对象配对结果与差值图。

两种模式的参数不能混用。输出目录必须是新的或空的，重跑请换一个目录名，避免混入旧结果。未指定输出时，single 默认使用输入根目录同级的 `SEM_0830_0831_single_V39`，compare 使用 after 同级的 `SEM_0830_0831_compare_V39`。

## 共同检测参数

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `--via-pattern` | `via60` / `via40` / `auto`，实际尺寸来自原图 | `via60` |
| `--via-min-circularity` | 可选圆形度下限；显式提高可能排除椭圆 | `0`（关闭） |
| `--via-max-size-deviation` | 可选同图等效直径偏差；非零时至少 3 个孔才执行 | `0`（关闭） |
| `--slit-min-column-rows` | 每个测量列的实际支持行数，须包含最下排；可选 2–4 | `2` |
| `--via-diameter-nm` | 预计 via 直径 | 不指定 |
| `--slit-width-nm` | 预计 slit 宽度，别名 `--slot-width-nm` | 不指定 |
| `--slit-length-nm` | 预计最下排 slit 长度，别名 `--slot-length-nm` | 不指定 |
| `--trench-width-nm` | 预计 trench 宽度 | 不指定 |
| `--trench-gap-nm` | 预计上下 trench 端点间隙 | 不指定 |
| `--width-tolerance` | via 直径及 slit/trench 宽度相对容差 | `0.50`（±50%） |
| `--length-tolerance` | slit 长度及 trench 间隙容差，别名 `--gap-tolerance` | `1.00`（±100%） |
| `--size-tolerance` | 旧版统一容差，独立宽度/长度参数优先 | 不指定 |
| `--image-timeout-seconds` | 单图硬时限，含补检及 auto 模型尝试 | `300` 秒，范围 `(0,300]` |

只填写已知尺寸，单位 nm，必须为有限正数。容差范围 `(0,1]`；同图 via 偏差筛选需显式启用且至少有 3 个有效孔。所有参数统一应用于本次样品/阶段，并保存到 `settings.json`；剔除原因保存到 `rejected_objects.csv`。auto 在 single 中每张 via 图独立选模型，在 compare 中由 before 选择并用于对应 after。

尺寸提示参与搜索与最终范围检查，测量仍须通过原图边缘验证。每张图按自己的 PixelSize 换算；超时或处理失败会记录并继续，缺失测量不补零。

## 使用说明、验证与历史

安装、目录示例、参数、统计含义、输出字段及故障处理见 [V39 完整使用说明](SEM_0830_0831_compare_V39_bundle/README_V39.txt)，算法见 [ALGORITHM_V39.md](SEM_0830_0831_compare_V39_bundle/ALGORITHM_V39.md)，修改、回归及 ZIP 校验见 [V39 修改记录](md/2026-10-16_V39_zero-result-recovery.md)。

全仓库 **288 项测试全部通过**（844.127 秒），其中 **28 项直接验证 V39**，覆盖零结果补检、原有识别能力与完整运行流程。本批真实 SEM 输入及 Windows/macOS 原生运行尚未验证。

```bash
MPLBACKEND=Agg python3 -m unittest discover -s tests -v
```

历史 [V38 使用说明](SEM_0830_0831_compare_V38_bundle/README_V38.txt)和 [V38 ZIP](SEM_0830_0831_compare_V38_bundle.zip)、[V37 使用说明](SEM_0830_0831_compare_V37_bundle/README_V37.txt)和 [V37 ZIP](SEM_0830_0831_compare_V37_bundle.zip)、[V36 使用说明](SEM_0830_0831_compare_V36_bundle/README_V36.txt)和 [V36 ZIP](SEM_0830_0831_compare_V36_bundle.zip)、[V35 使用说明](SEM_0830_0831_compare_V35_bundle/README_V35.txt)和 [V35 ZIP](SEM_0830_0831_compare_V35_bundle.zip)、[V34 使用说明](SEM_0830_0831_compare_V34_bundle/README_V34.txt)和 [V34 ZIP](SEM_0830_0831_compare_V34_bundle.zip)，以及更早版本目录、ZIP 和上传的 V27/Yakun 包均保留。

仓库工作规则见 [AGENTS.md](AGENTS.md)，所有修改记录见 [md/](md/README.md)。
