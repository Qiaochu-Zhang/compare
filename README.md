# SEM 样品分析与 before/after 对比

当前版本：[V36](SEM_0830_0831_compare_V36_bundle/README_V36.txt)，新增**平滑非圆 via 识别与 slit 局部污染缺失处理**，保留单文件夹及 before/after 两种模式。

下载 [V36 完整压缩包](SEM_0830_0831_compare_V36_bundle.zip)。解压并进入同名目录，安装依赖：

```bash
python -m pip install -r requirements.txt
```

## V36 识别改进

- **via 接纳拉长的圆/椭圆**：默认关闭圆形度和同图等效直径偏差的硬筛选；通过原图边缘、凸性、椭圆归一化形状偏差和边缘粗糙度排除明显锯齿、星形、多瓣异常。拟合椭圆只用于检查，输出保留实际轮廓。
- **避免长孔重复计数**：补充整体轮廓定位，按实际阵列格位合并候选，并对细化后的真实轮廓去重。
- **slit 允许局部缺失**：全图仍须有四排阵列证据，单个测量列默认有“最下排 + 至少一个其他排”的真实边缘即可。上排缺失不再自动删掉清晰的最下排；最下排缺失时不补造尺寸。
- `--slit-min-column-rows 4` 可恢复每列必须四排完整的要求。缺失位置、形状拒绝和实际支持行数都写入结果，统计/作图规则保持一致。

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
python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36"
```

`--input` 自动选择单文件夹模式，也可显式加 `--mode single`。各样品图片数可以不同；编号保留前导零，按数值排序。

图片仍按文件名自然排序，每三张依次分配 **trench、slot、via**，Region 从 2 开始；类型由排序位置决定。尾部可以不满三张，中间缺图会导致后续类型错位。支持 TIF/TIFF/PNG/JPG/JPEG，每张图片须有同名 TXT，包含 `PixelSize=正数`，单位 nm/pixel。只扫描直接数字子目录内的图片。

统计作图继续使用六项尺寸指标、nm 单位、均值柱和样本标准差误差棒；slot 只统计满足阵列支持要求且边缘可确认的最下排。生成单样品各区域图、全部样品区域总览图、合并区域后的样品汇总图。汇总按有效对象计算，n=1 的标准差为空，无有效值不填 0。

主要结果：`SEM_0830_0831_single_V36_results.xlsx`、`sample_metric_summary.csv`、单图统计与对象明细、`annotated/`、`plots/`。单文件夹模式不输出前后配对和差值文件。

## before/after 对比（保留）

```bash
python sem_before_after_compare_0831_V36.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V36"
```

可显式加 `--mode compare`。从 after 发现数字目录，与 before 完全同名目录匹配；对应目录图片数量须相等，按自然排序配对，文件名和格式可以不同。输出 `SEM_0830_0831_before_after_V36_results.xlsx`，保留图片级、总体及逐对象配对结果与差值图。

两种模式的参数不能混用。输出目录必须是新的或空的，重跑请换一个目录名，避免混入旧结果。未指定输出时，single 默认使用输入根目录同级的 `SEM_0830_0831_single_V36`，compare 使用 after 同级的 `SEM_0830_0831_compare_V36`。

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

完整的安装步骤、目录示例、两种模式、全部参数、统计含义、输出字段和故障处理见 [V36 完整使用说明](SEM_0830_0831_compare_V36_bundle/README_V36.txt)，实现说明见 [ALGORITHM_V36.md](SEM_0830_0831_compare_V36_bundle/ALGORITHM_V36.md)，本次修改与交付验证见 [V36 修改记录](md/2026-10-15_V36_morphology-occlusion.md)。

全仓库 **167 项测试通过**，其中 **27 项直接验证 V36**，覆盖椭圆真实边缘、异常轮廓拒绝、缺失 slit、统计与真实 CLI/worker；完整结果见修改记录。仓库暂无本批真实 SEM 原图与比例尺 TXT，实际识别结果仍需结合标注复核。

```bash
MPLBACKEND=Agg python3 -m unittest discover -s tests -v
```

历史 [V35 使用说明](SEM_0830_0831_compare_V35_bundle/README_V35.txt)和 [V35 ZIP](SEM_0830_0831_compare_V35_bundle.zip)、[V34 使用说明](SEM_0830_0831_compare_V34_bundle/README_V34.txt)和 [V34 ZIP](SEM_0830_0831_compare_V34_bundle.zip)，以及更早版本目录、ZIP 和上传的 V27/Yakun 包均保留。

仓库工作规则见 [AGENTS.md](AGENTS.md)，所有修改记录见 [md/](md/README.md)。
