# SEM 样品分析与 before/after 对比

当前版本：[V38](SEM_0830_0831_compare_V38_bundle/README_V38.txt)，改善 **V37 相对 V35 的 slit 识别回退**，保留 V37 的 via / trench 改进，以及单文件夹、before/after 两种模式。

下载 [V38 完整压缩包](SEM_0830_0831_compare_V38_bundle.zip)。解压并进入同名目录，安装依赖：

```bash
python -m pip install -r requirements.txt
```

## V38 识别改进

- **slit 完整性按局部横向背景判断**：V37 的大窗口背景平面和单点外探，会把弯曲照明、邻位或模糊尾部误认为轮廓未完整。V38 检查连续外侧前景，减少误拒绝和不必要的尺寸扩张。
- **首轮可靠定位优先**：沿用 V35 首轮定位顺序；补检独立运行，已有有效轮廓和四排行级优先保留。失败种子恢复同格位原图局部修复，仍须通过原图边缘审核。
- **完整列可保留，部分列按阵列证据处理**：四排齐全且边缘有效的列不会因其他列缺失被一并删除。部分列继续要求重复阵列支持，缺最下排不测，不提升第三排、不补零。
- **保留 V37 via 改进**：弱孔补检、完整轮廓、平滑椭圆、异常形貌筛查和去重继续生效；trench 保留局部条体及端点识别。

受控合成对照（表中为有效最下排 slit 数量）：

| 场景 | V35 | V37 | V38 |
| --- | --- | --- | --- |
| 轻微弯曲背景 | 3 | 0 | 3 |
| 密集邻位 | 3 | 0 | 3 |
| 一列完整、其余列缺失较多 | 1 | 0 | 1 |

模糊场景中，真实宽约 17 px，V37 测到约 22 px，V38 恢复到约 16.3–17.6 px；弱下半部的完整轮廓恢复也继续验证。原因、复现参数及边界见 [V38 算法说明](SEM_0830_0831_compare_V38_bundle/ALGORITHM_V38.md)。这些是合成验证，尚无本批真实 SEM 原图/TXT，不能当作实际准确率。

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
python sem_before_after_compare_0831_V38.py --input "D:/SEM/samples" --output "D:/SEM/single_V38"
```

`--input` 自动选择单文件夹模式，也可显式加 `--mode single`。各样品图片数可以不同；编号保留前导零，按数值排序。

图片仍按文件名自然排序，每三张依次分配 **trench、slot、via**，Region 从 2 开始；类型由排序位置决定。尾部可以不满三张，中间缺图会导致后续类型错位。支持 TIF/TIFF/PNG/JPG/JPEG，每张图片须有同名 TXT，包含 `PixelSize=正数`，单位 nm/pixel。只扫描直接数字子目录内的图片。

统计作图继续使用六项尺寸指标、nm 单位、均值柱和样本标准差误差棒；slot 只统计边缘有效的完整列最下排，或满足阵列支持要求的部分列最下排。生成单样品各区域图、全部样品区域总览图、合并区域后的样品汇总图。汇总按有效对象计算，n=1 的标准差为空，无有效值不填 0。

主要结果：`SEM_0830_0831_single_V38_results.xlsx`、`sample_metric_summary.csv`、单图统计与对象明细、`annotated/`、`plots/`。单文件夹模式不输出前后配对和差值文件。

## before/after 对比（保留）

```bash
python sem_before_after_compare_0831_V38.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V38"
```

可显式加 `--mode compare`。从 after 发现数字目录，与 before 完全同名目录匹配；对应目录图片数量须相等，按自然排序配对，文件名和格式可以不同。输出 `SEM_0830_0831_before_after_V38_results.xlsx`，保留图片级、总体及逐对象配对结果与差值图。

两种模式的参数不能混用。输出目录必须是新的或空的，重跑请换一个目录名，避免混入旧结果。未指定输出时，single 默认使用输入根目录同级的 `SEM_0830_0831_single_V38`，compare 使用 after 同级的 `SEM_0830_0831_compare_V38`。

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

完整的安装步骤、目录示例、两种模式、全部参数、统计含义、输出字段和故障处理见 [V38 完整使用说明](SEM_0830_0831_compare_V38_bundle/README_V38.txt)，实现说明见 [ALGORITHM_V38.md](SEM_0830_0831_compare_V38_bundle/ALGORITHM_V38.md)，本次修改与交付验证见 [V38 修改记录](md/2026-10-16_V38_slit-via-fusion.md)。

全仓库 **252 项测试全部通过**，其中 **47 项直接验证 V38**：9 项 V35/V37/V38 回退对照与分图案验证，加上 38 项完整轮廓、形貌、污染缺失及工作流场景。实际命令、结果和交付校验见修改记录；本批真实 SEM 原图/TXT 尚未验证。

```bash
MPLBACKEND=Agg python3 -m unittest discover -s tests -v
```

历史 [V37 使用说明](SEM_0830_0831_compare_V37_bundle/README_V37.txt)和 [V37 ZIP](SEM_0830_0831_compare_V37_bundle.zip)、[V36 使用说明](SEM_0830_0831_compare_V36_bundle/README_V36.txt)和 [V36 ZIP](SEM_0830_0831_compare_V36_bundle.zip)、[V35 使用说明](SEM_0830_0831_compare_V35_bundle/README_V35.txt)和 [V35 ZIP](SEM_0830_0831_compare_V35_bundle.zip)、[V34 使用说明](SEM_0830_0831_compare_V34_bundle/README_V34.txt)和 [V34 ZIP](SEM_0830_0831_compare_V34_bundle.zip)，以及更早版本目录、ZIP 和上传的 V27/Yakun 包均保留。

仓库工作规则见 [AGENTS.md](AGENTS.md)，所有修改记录见 [md/](md/README.md)。
