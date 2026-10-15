V35 完整使用说明 — 单文件夹样品分析 + before/after 对比
====================================================

V35 基于 V34，新增 single 模式：只提供一个根目录，其下每个纯数字子目录代表一个
独立样品。各样品分别识别、测量、统计并作图，不需要 before/after 文件夹。
原有 compare 模式保留；两种模式共用 V34 测量、质量筛选和尺寸提示规则。
本包可独立运行，不需要下载或导入历史版本。

一、安装与首次运行
------------------
1. 解压 SEM_0830_0831_compare_V35_bundle.zip，进入同名目录。
   应能看到以下四个文件：
     sem_before_after_compare_0831_V35.py
     requirements.txt
     README_V35.txt
     ALGORITHM_V35.md
2. 准备 Python 3.9 或更新版本（本版测试环境为 Python 3.12）。
   Windows 可在该目录打开 PowerShell / 命令提示符；macOS/Linux 打开终端并 cd 到目录。
   如果系统使用 python3 命令，请把下面的 python 换成 python3；Windows 也可使用 py。
3. 安装依赖：
     python -m pip install -r requirements.txt
   包括 numpy、pandas、scipy、opencv-python、matplotlib、openpyxl。
   可选：先 python -m venv .venv，再激活虚拟环境后安装。
     Windows PowerShell：.\.venv\Scripts\Activate.ps1
     Windows cmd：.venv\Scripts\activate.bat
     macOS/Linux：source .venv/bin/activate
4. 按下文组织输入，然后运行单文件夹模式：
     python sem_before_after_compare_0831_V35.py --input "D:/SEM/samples" --output "D:/SEM/single_V35"
   macOS/Linux 示例：
     python3 sem_before_after_compare_0831_V35.py --input "/data/SEM/samples" --output "/data/SEM/single_V35"
   路径含空格或中文时保留引号。程序通过命令行运行，无需修改脚本中的旧默认路径。
5. 运行后先核对 inventory_mapping.csv，再检查 image_status.csv、annotated/ 和统计图。
   不要仅依据命令退出成功就判断所有图片均有有效测量。

二、单文件夹模式：目录组织和图片顺序
----------------------------------
--input 指向所有样品的父目录，不是某一个编号目录，也不是单张图片。

  samples/
    001/                    样品 001
      image_1.tif           Region 2，trench
      image_1.txt
      image_2.png           Region 2，slot（slit，四排结构）
      image_2.txt
      image_3.jpg           Region 2，via
      image_3.txt
      image_4.tif           Region 3，trench
      image_4.txt
      image_5.tif           Region 3，slot
      image_5.txt
      image_6.tif           Region 3，via
      image_6.txt
    002/                    样品 002，可以比 001 少几个区域
      image_1.tif           Region 2，trench
      image_1.txt
      image_2.tif           Region 2，slot
      image_2.txt
      image_3.tif           Region 2，via
      image_3.txt
    notes/                  非数字目录，不处理

目录发现规则：
- 只读取根目录的直接子目录，且名称必须全部由 0–9 组成，例如 0、2、007、10、103。
- 按数字大小排序；相同数值按完整名称排序。007 和 7 是两个不同样品，保留前导零。
- 不递归扫描。样品图片须直接放在数字目录内；notes/、sample1/、1_test/ 等不会处理。
- 单文件夹模式各样品图片数独立，不要求相等，不要求固定为 9 张或为 3 的倍数。
- 空样品目录记录 EMPTY_SAMPLE_FOLDER，继续其他样品；完全没有可处理图片则退出报错。

图案/区域分配规则（沿用旧版）：
- 支持 .tif / .tiff / .png / .jpg / .jpeg，扩展名大小写不限，也可混用。
- 每个样品内部按文件名自然排序，例如 image_1、image_2、image_10。
- 排序后的第 1、2、3 张分别是 trench、slot、via，记为 Region 2；
  第 4、5、6 张同样依次为 trench、slot、via，记为 Region 3，以此类推。
- 程序按顺序分配类型，不根据文件名中的 via/slot/trench 文字自动切换类型。
  建议用 001_trench、002_slot、003_via、004_trench 等序号命名。
- 末尾不满三张可以处理，例如第 10 张为 Region 5 trench、第 11 张为 Region 5 slot。
  中途缺失图片会使后续类型错位；必须按完整拍摄顺序整理，不要随意删掉中间图片。
  仅有一张图片时会按 trench 处理，不能把单张 via 当成第 1 张使用。
- 缺 TXT 或比例尺无效的现存图片仍保留其位置，只跳过这张，后续图片不重新编号。
- 不要把标注图或历史结果放入样品目录，否则也会作为输入图片参与排序。
  推荐输出到输入根目录之外，并先检查 inventory_mapping.csv 的 image/pattern/region。

每张图片的比例尺：
- 每张图必须有同目录、同主文件名的 TXT，例如 image_1.tif 对应 image_1.txt。
- TXT 至少包含一行（可带其他仪器信息）：
    PixelSize=0.8268229
- 数值单位是 nm/pixel，必须为有限正数。不同图片可以使用各自不同的比例尺。
  不能把像素数量、标尺总长度或 µm/pixel 数值不换算就写入。
- TXT 扩展名大小写兼容；应确保 Windows 没有把文件保存为 image_1.txt.txt。
- 无 TXT、无 PixelSize、非正数等会写入 processing_errors.csv；不猜测比例尺。

三、选择模式与常用命令
----------------------
最简单文件夹分析（提供 --input 即自动选择 single）：
  python sem_before_after_compare_0831_V35.py --input "D:/SEM/samples" --output "D:/SEM/single_V35"

显式选择 single，效果相同：
  python sem_before_after_compare_0831_V35.py --mode single --input "D:/SEM/samples" --output "D:/SEM/single_V35_run2"

原有 before/after 对比（不加 --mode 仍可使用）：
  python sem_before_after_compare_0831_V35.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V35"
也可加 --mode compare。

模式参数不能混用：single 必须有 --input，不能再给 --before / --after；
compare 不能给 --input。未传 --input 时默认 compare，未传 before/after 路径时仍使用
脚本保留的历史 Windows 默认路径；实际使用请明确提供自己的路径。

compare 的附加规则：
- 从 after 的直接数字子目录发现条件，匹配 before 下完全同名的目录。
- 只存在于 before 的目录不处理；007 不能匹配 7。
- 同条件 before/after 图片数量必须相同，文件名及格式可不同，按各自自然排序配对。
- 缺匹配目录或数量不等时跳过该条件；某张比例尺预检失败时跳过对应前后图片组。
- 单文件夹模式完全不进行前后或样品间的对象配对。

输出路径：
- --output 可指定新的或空的目录；V35 拒绝非空目录，避免旧图或另一模式的结果混入。
  重新运行请换目录名，或自行整理旧输出后使用空目录。
- 不指定时，single 在输入根目录同级创建 SEM_0830_0831_single_V35；
  compare 在 after 根目录同级创建 SEM_0830_0831_compare_V35。

四、测量、统计与作图规则
----------------------
两种模式保留相同的六个尺寸指标，单位为 nm：
  trench：trench_width_nm（上下条平均宽度）、tip_to_tip_y_nm（上下端点 Y 间隙）
  slot：slot_x_width_nm（最下排 X 宽度）、slot_y_length_nm（最下排 Y 长度）
  via：via_x_width_nm（X 宽度）、via_y_height_nm（Y 高度）
对象明细中的源字段分别为 x_width_nm、y_height_nm、tip_gap_y_nm。
slit 和 slot 表示同类结构，输出沿用 slot 名称。

- slot 仍要求四排组成完整有效列，只统计完整列的最下排；不是把全部四排混在一起。
- via 仍使用原图轮廓、圆形度与同图尺寸筛选；trench 仍复核原图端点及宽度。
- 尺寸由原始灰度边缘测量；增强/补检与尺寸提示只辅助定位，不代替真实边缘。
- 单文件夹图按均值画柱，误差棒为样本标准差 SD（ddof=1），与原对比图一致。
  表中还保存 n、SEM=SD/sqrt(n)、中位数、最小/最大值；误差棒不使用 SEM。
- n 表示有效对象数，不是图片数。n=1 时保留均值，SD/SEM 在表中为空，图上不画误差范围。
  没有有效对象时不生成该图该指标的统计行/柱，不以 0 nm 填补缺失结果。
- 样品级汇总把同编号、同图案各区域的有效对象合并后统计，延续原 overall 规则。
  因此是按对象数量加权，不是把各张图的均值再次等权平均。
  sample_metric_summary.csv 同时给出 image_count、region_count（实际有有效值的数量）。
- 样品/条件编号按数值排序，区域按 Region 编号排序；保留 007 等原始编号。
- 图为 220 dpi PNG；三个图案各两项尺寸，存在有效数据时才生成相应图。

single 生成三类图：
  C001__trench__trench_width_nm__single.png
    一个样品的一项指标，横轴为各 Region，纵轴 nm，均值 ± SD。
  ALL__trench__trench_width_nm__single.png
    所有样品的各 Region 总览，横轴如 C001-R2、C001-R3、C002-R2。
  SAMPLES__trench__trench_width_nm__single.png
    合并区域后的样品汇总图，横轴为样品编号，各柱为该样品有效对象的均值 ± SD。
  slot/via 与另一尺寸指标按同样命名规则生成。

compare 继续生成每条件/指标的 before/after 并排柱（各自均值 ± SD），
以及 ALL__...__delta.png（After − Before，单位 nm）。差值图本身不加 SD 误差棒。

五、所有命令行选项
------------------
  -h, --help                  显示完整帮助后退出
  --mode single|compare       不指定时根据有无 --input 自动选择
  --input PATH                单文件夹模式样品根目录
  --before PATH               对比模式处理前根目录
  --after PATH                对比模式处理后根目录
  --output PATH               新的/空的结果目录，默认位置见第三节
  --image-timeout-seconds N    每张图片的总处理时限，默认 300 秒；必须 >0 且 <=300
  --via-pattern via60|via40|auto
                              默认 via60；标称模型名称不代表把结果固定为 60/40 nm
  --via-min-circularity N      via 圆形度下限，默认 0.75，范围 [0,1]；0 关闭该阈值
  --via-max-size-deviation N   via 等效直径相对同图中位数的最大偏差，默认 0.30
                              范围 [0,1)，0 关闭；至少 3 个有效 via 才执行
  --via-diameter-nm N          预计 via 直径，nm，默认不指定
  --slit-width-nm N            预计 slit 宽度，nm；别名 --slot-width-nm
  --slit-length-nm N           预计最下排 slit 长度，nm；别名 --slot-length-nm
  --trench-width-nm N          预计 trench 宽度，nm
  --trench-gap-nm N            预计上下 trench 端点间隙，nm，不是横向阵列间距
  --width-tolerance N         via 直径和 slit/trench 宽度相对容差，默认 0.50（±50%）
  --length-tolerance N        最下排 slit 长度和 trench 间隙容差，默认 1.00（±100%）
                              别名 --gap-tolerance
  --size-tolerance N          兼容旧版，统一设置两类容差，默认不指定

via 筛选详解：
- 圆形度 C=4π×面积/周长²，使用最终原图轮廓，越接近 1 越圆。
  提高阈值更严格；降低阈值不取消边缘、完整性、solidity、阵列等其他要求。
  该阈值只作用于 via，不删除细长的 slit/trench；补检也不能绕过阈值。
- 等效直径 D=2√(面积/π)。从同图通过边缘、圆形度及尺寸提示范围检查的 via 中
  计算一次中位数 M，|D/M−1| 大于设定值的对象剔除；不足 3 个时记录跳过筛选。
  single 按每张图独立计算，compare 的 before/after 也独立计算，不按样品间或前后
  差异删除真实尺寸变化。
- auto 会试 via40/via60，选通过验证的较可靠结果并复用标注。
  single 每张 via 图片独立选模型；compare 由 BEFORE 选并固定给对应 AFTER。
  单文件夹的中性 single 阶段沿用通用识别分支，不启用历史 AFTER 专用增强。

尺寸提示详解：
- 仅填写已知项；尺寸必须有限且 >0，不知道时省略，不填 0。
- 每张图按自己的 PixelSize 把 nm 转成像素，用于搜索尺度、候选范围和最终范围检查。
- 两类容差必须 >0 且 <=1。±50% 为预计值的 0.5–1.5 倍；±100% 为正值至预计值的 2 倍。
  0、负尺寸、没有边缘证据的对象不因容差宽而通过。
- slit 四排可以各自有不同长度，长度提示只检查最下排；四排完整性仍保留。
- 独立的 width/length-tolerance 优先于旧 size-tolerance，和命令参数顺序无关。
  例如 --size-tolerance 0.3 --length-tolerance 1 表示宽度 ±30%、长度/间隙 ±100%。
- 提示及阈值统一作用于本次所有样品/阶段。样品差异很大时可省略尺寸提示或分批运行。
  提示过窄可能排除真实变化；最终输出不会被强制替换为输入尺寸。

示例（数值仅演示写法，请替换成自己的尺寸）：
  python sem_before_after_compare_0831_V35.py --input "D:/SEM/samples" --output "D:/SEM/single_V35_hints" --via-diameter-nm 80 --slit-width-nm 80 --slit-length-nm 210 --trench-width-nm 160 --trench-gap-nm 120

调整 via 质量阈值并缩短单图等待时间：
  python sem_before_after_compare_0831_V35.py --input "D:/SEM/samples" --output "D:/SEM/single_V35_strict" --via-min-circularity 0.85 --via-max-size-deviation 0.20 --image-timeout-seconds 120

六、输出文件和复核顺序
----------------------
两种模式都有：
  inventory_mapping.csv       图片清单：编号、文件路径、顺序、图案、Region、比例尺预检
  image_status.csv            每张已处理/已记录跳过图片的状态、数量、比例尺、耗时、标注路径
  all_object_measurements.csv 每个有效对象的尺寸、质量及原图边缘检查、输入参数
  image_metric_summary.csv    每张有效图片的每个尺寸指标统计
  rejected_objects.csv        剔除对象及原因，可能包括不同检测/补检尝试的记录
  via_model_choice.csv        via 模型选择记录；single 不使用 before_image/after_image 字段
  processing_errors.csv      目录、比例尺、处理失败、超时等错误与警告
  settings.json              版本、mode、输入路径、编号、实际参数、统计/作图规则及数量
  annotated/                 原图标注，single 位于 single/编号/Region_N/ 下
                              compare 位于 before/编号/Region_N/ 和 after/编号/Region_N/ 下
  plots/                     统计图，命名见第四节

single 额外输出：
  sample_metric_summary.csv
  SEM_0830_0831_single_V35_results.xlsx
    工作表：inventory、image_status、all_objects、image_metrics、sample_summary、
            rejected_objects、via_model_choice、errors、settings。
  不生成前后配对表、变化百分比或差值图。

compare 额外输出：
  image_pair_comparison.csv                 图片级 before/after 统计与差值
  overall_comparison.csv                    同条件各区域合并后统计与差值
  object_before_after_comparison.csv        逐对象空间配对结果；不改变实际检测值
  SEM_0830_0831_before_after_V35_results.xlsx
    含通用工作表，以及 image_pair_compare、overall_compare、object_compare。
    前后差值是 after − before；change_percent 是差值 / before × 100%，
    缺失或 before 接近 0 时为空。逐对象匹配只是位置对应推断，需结合标注复核。

字段约定：
- condition 始终保存目录原始编号，single 中即样品编号；stage 为 single、before 或 after。
- pair_key 在 single 中仅为编号/Region/图案的分组标识，不意味着建立了前后配对。
- selected_v25_pattern_key、内部 _v32_* 等名称为继承兼容字段/函数，仍实际运行 V35。
- circularity、via_equivalent_diameter_nm、via_image_median_diameter_nm、
  via_size_deviation_fraction、via_size_reference_count、via_size_filter_applied 用于复核 via 筛选。
- input_width_tolerance / input_length_tolerance 是最终生效值；input_size_tolerance
  是是否显式指定旧统一参数的记录。edge_* 用于复核原图边缘证据。
- CSV 使用 UTF-8 BOM。Excel 直接打开 CSV 可能把 007 变成 7，建议打开本程序的 XLSX，
  或导入 CSV 时把 condition 列设为文本；程序导出的 XLSX 保留编号字符串。
- single 无有效测量时仍输出有表头的空明细/统计表，结合状态和错误表判断原因。

建议复核顺序：inventory → image_status/errors → annotated → rejected_objects → 统计表/plots。
标注与有效统计消费同一组通过筛选的对象；异常补检结果也应核对原图。

七、常见问题、失败处理和验证范围
------------------------------
1. 找不到数字子目录：确认 --input 指向样品父目录，目录名称为纯数字，图片没有额外嵌套。
2. 类型不对或结果异常：先核对自然排序；重命名图片时，同名 TXT 也要对应改名。
3. PixelSize 错误：检查 TXT 是否同名、是否有 PixelSize=正数、单位是否为 nm/pixel。
   single 只跳过对应图片，不影响其他样品。compare 跳过对应前后组。
4. 结果为 CHECK_NO_VALID_MEASUREMENT...：图片已处理，但原图证据或质量要求不足，
   不会编造尺寸。查看标注、剔除日志；核实类型、比例尺及尺寸提示，必要时换清晰原图。
   低圆形度/尺寸异常常见原因：via_circularity_below_threshold、via_size_outlier、
   via_diameter_outside_hint_range、slit_size_outside_hint_range、trench_size_outside_hint_range。
5. IMAGE_TIMEOUT / SKIPPED_TIMEOUT：单图超过设定时限即结束该 worker，记录并继续。
   时限包括补检及 auto 的两个模型尝试。超时图不会发布未完成标注。
6. 输出目录非空：换一个目录名；不要混用同一目录保存 single/compare 或不同参数结果。
7. 报 ModuleNotFoundError：确认在运行脚本的同一个 Python 环境执行 python -m pip install。
8. Excel 文件无法写入：关闭正在打开的同名 Excel 文件并使用新的可写输出目录。
9. 运行较久：终端会报告图片处理进度；300 秒是单图上限，不是整批时限。
   Ctrl+C 中断后不会自动续跑，已完成状态日志可查看，重跑需新的输出目录。
10. 非零退出码表示参数、目录、导出等致命错误或用户中断；退出码 0 表示批处理/导出完成，
    仍可能有图片失败、被跳过或没有有效对象，请检查 image_status 和 processing_errors。
    只有空白图、或所有比例尺失败时，统计图可能为空，这是没有有效数据，不是 0 nm。

V35 自动测试直接加载 V35：覆盖目录/排序/前导零、样品数量独立、比例尺错误、
超时继续、参数冲突、统计权重、SD 误差棒、两种模式导出和对比兼容。
真实 CLI 与独立 worker 使用合成 trench/slot/via 图片，验证不同比例尺、auto 模型、
标注、CSV/Excel 和 PNG；空白图拒测及硬超时也纳入验证。
仓库没有本批真实 SEM 原图及对应 TXT，合成图验证不代表真实样品识别准确率已确认。
请结合实际原图复核，严重模糊时尺寸提示不能保证恢复可靠测量。

开发测试（需完整仓库中的 tests/，发布 ZIP 不包含测试数据）：
  MPLBACKEND=Agg python3 -m unittest discover -s tests -v
算法、数据流和兼容说明见同目录 ALGORITHM_V35.md。
