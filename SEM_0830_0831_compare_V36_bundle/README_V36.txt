V36 完整使用说明 — 平滑非圆 via 与允许局部缺失的 slit 阵列
====================================================

V36 基于 V35：via 默认不按圆形度或同图等效直径偏差硬筛选，允许平滑拉长的孔；
通过原图边缘、轮廓凹陷/锯齿/不规则度及阵列位置排除异常对象。slit 仍按四排阵列
识别，但允许部分位置因污染而缺失，不再因一个上排成员缺失就删除整个测量列。
保留 single 单文件夹样品模式和 compare 前后对比模式，统计作图规则保持一致。
本包可独立运行，不需要下载或导入历史版本。

一、安装与首次运行
------------------
1. 解压 SEM_0830_0831_compare_V36_bundle.zip，进入同名目录。
   应能看到以下四个文件：
     sem_before_after_compare_0831_V36.py
     requirements.txt
     README_V36.txt
     ALGORITHM_V36.md
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
     python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36"
   macOS/Linux 示例：
     python3 sem_before_after_compare_0831_V36.py --input "/data/SEM/samples" --output "/data/SEM/single_V36"
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
  python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36"

显式选择 single，效果相同：
  python sem_before_after_compare_0831_V36.py --mode single --input "D:/SEM/samples" --output "D:/SEM/single_V36_run2"

原有 before/after 对比（不加 --mode 仍可使用）：
  python sem_before_after_compare_0831_V36.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V36"
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
- --output 可指定新的或空的目录；V36 拒绝非空目录，避免旧图或另一模式的结果混入。
  重新运行请换目录名，或自行整理旧输出后使用空目录。
- 不指定时，single 在输入根目录同级创建 SEM_0830_0831_single_V36；
  compare 在 after 根目录同级创建 SEM_0830_0831_compare_V36。

四、测量、统计与作图规则
----------------------
两种模式保留相同的六个尺寸指标，单位为 nm：
  trench：trench_width_nm（上下条平均宽度）、tip_to_tip_y_nm（上下端点 Y 间隙）
  slot：slot_x_width_nm（最下排 X 宽度）、slot_y_length_nm（最下排 Y 长度）
  via：via_x_width_nm（X 宽度）、via_y_height_nm（Y 高度）
对象明细中的源字段分别为 x_width_nm、y_height_nm、tip_gap_y_nm。
slit 和 slot 表示同类结构，输出沿用 slot 名称。

- slot 仍只统计最下排，允许上排部分成员缺失。全图必须有四排的实际观测证据，
  至少两个受支持的列，且至少三排各有两个真实对象；单个测量列默认至少有两排
  边缘支持，其中必须包括最下排。没有最下排轮廓的列不测量，也不把第三排当第四排。
- via 以原图轮廓为准，默认圆形度只记录，不作阈值淘汰；以形貌规则排除明显异常。
  trench 继续复核原图端点及宽度。
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
  --via-min-circularity N      可选 via 圆形度下限，默认 0（关闭），范围 [0,1]
  --via-max-size-deviation N   可选同图等效直径偏差限制，默认 0（关闭），范围 [0,1)
                              非零时至少 3 个有效 via 才执行；可能剔除真实拉长形变
  --slit-min-column-rows N     测量一个 slit 列所需的真实行数，默认 2；可选 2、3、4
                              别名 --slot-min-column-rows；4 恢复每列必须完整的要求
  --via-diameter-nm N          预计 via 直径，nm，默认不指定
  --slit-width-nm N            预计 slit 宽度，nm；别名 --slot-width-nm
  --slit-length-nm N           预计最下排 slit 长度，nm；别名 --slot-length-nm
  --trench-width-nm N          预计 trench 宽度，nm
  --trench-gap-nm N            预计上下 trench 端点间隙，nm，不是横向阵列间距
  --width-tolerance N         via 直径和 slit/trench 宽度相对容差，默认 0.50（±50%）
  --length-tolerance N        最下排 slit 长度和 trench 间隙容差，默认 1.00（±100%）
                              别名 --gap-tolerance
  --size-tolerance N          兼容旧版，统一设置两类容差，默认不指定

via 形貌与边缘筛选详解：
- 圆形度 C=4π×面积/周长² 仍记录在结果中，但默认 --via-min-circularity 0。
  3:1 或更长的椭圆不因低圆形度、旧版 1.62 长宽比限制或圆心半径变化大而被删除。
  圆形度不再是默认识别目标，原始像素的真实边缘始终是测量依据。
- 为避免拉长造成面积/等效直径变化后被另一项默认筛选删除，V36 也将同图尺寸
  偏差筛选默认设为关闭。仍保留图内搜索尺度、图像边界和阵列邻位范围检查。
  明确提供预计直径或非零 --via-max-size-deviation 时，仍按用户参数做尺寸筛选。
- 在轮廓上拟合椭圆仅用于归一化形状检查，消除正常拉长/旋转对不规则度的影响。
  输出轮廓不替换成拟合椭圆，不为提高圆形度移动边缘。
- 原图边缘验证后，检查 solidity（面积/凸包面积）>=0.85、轮廓凸性>=0.88，
  椭圆归一化半径偏差的 95 分位 <=0.18、高频边缘粗糙度 <=0.06。
  这些组合检查用于拒绝明显锯齿、星形、大凹口和多瓣异常形状；补检不放宽它们。
  门槛为合成图验证下的默认值，尚未经本批真实样品标定。
- 新增整体连通轮廓定位，避免旧圆形搜索漏掉椭圆长轴端点；按实际阵列格位合并
  多个定位结果，最终轮廓再去重，避免把一个拉长孔统计成两个或多个孔。
- 每个有效轮廓仍须满足原灰度边缘支持率至少 85%、连续缺边段不超过约 8% 的
  原有检查。只有少量缺边时可沿真实支持端点连接；大面积遮挡不会凭模型补出尺寸。
- 如确需恢复圆形度门槛，可显式设 --via-min-circularity 0.75；该设置可能删除椭圆。
  同图偏差筛选显式启用后，从通过形状/提示检查的 via 计算一次中位等效直径 M，
  剔除 |D/M−1| 超限对象；D=2√(面积/π)。不足 3 个有效孔时跳过同图筛选。
  single 每张图独立处理，compare 的 before/after 也各自独立，不按前后差异剔除。
- auto 会试 via40/via60，并复用较可靠的测量与标注。single 每张 via 图片独立选模型；
  compare 由 BEFORE 选并固定给 AFTER。single 不启用历史 AFTER 专用增强。

slit 污染缺失处理详解：
- 四排是全图形貌约束，单列不必每个位置都能识别。默认一个列有“最下排 + 任意
  一个以上其他排”的真实轮廓支持即可；全图还必须满足第四节的重复阵列证据。
- 例如三个列分别缺少第一、第二、第三排中的一个对象，仍可测量三个清晰的最下排。
  如果某列最下排本身被污染挡住，则该列没有尺寸，不会填零或使用邻列尺寸代替。
- 若整排消失、只剩零散对象、无法确认哪一排是第四排，程序会拒绝统计该阵列；
  不能把“四排先验”当成补造被污染对象的依据。
- --slit-min-column-rows 3 要求更强的列支持；设 4 要求每列四个对象都真实有效。
  这些设置都不改变逐对象原图边缘检查，也不自动恢复被遮挡的边缘。
- 标注只画可确认的实际轮廓；绿色最下排对象参与统计，辅助轮廓可能仅用于阵列
  定位，不参与统计。缺失位置记录到 rejected_objects.csv，不画虚构轮廓。

尺寸提示详解：
- 仅填写已知项；尺寸必须有限且 >0，不知道时省略，不填 0。
- 每张图按自己的 PixelSize 把 nm 转成像素，用于搜索尺度、候选范围和最终范围检查。
- 两类容差必须 >0 且 <=1。±50% 为预计值的 0.5–1.5 倍；±100% 为正值至预计值的 2 倍。
  0、负尺寸、没有边缘证据的对象不因容差宽而通过。
- slit 四排可以各自有不同长度，长度提示只检查最下排；全图四排形貌和每列最少真实行数要求仍保留。
- 独立的 width/length-tolerance 优先于旧 size-tolerance，和命令参数顺序无关。
  例如 --size-tolerance 0.3 --length-tolerance 1 表示宽度 ±30%、长度/间隙 ±100%。
- 提示及阈值统一作用于本次所有样品/阶段。样品差异很大时可省略尺寸提示或分批运行。
  提示过窄可能排除真实变化；最终输出不会被强制替换为输入尺寸。

示例（数值仅演示写法，请替换成自己的尺寸）：
  python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36_hints" --via-diameter-nm 80 --slit-width-nm 80 --slit-length-nm 210 --trench-width-nm 160 --trench-gap-nm 120

需要主动收紧圆形度/尺寸筛选时（会排除部分椭圆），可用：
  python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36_strict" --via-min-circularity 0.85 --via-max-size-deviation 0.20 --image-timeout-seconds 120

有污染缺失的 slit 默认无需额外参数；要切回每列完整检查：
  python sem_before_after_compare_0831_V36.py --input "D:/SEM/samples" --output "D:/SEM/single_V36_full_columns" --slit-min-column-rows 4

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
  SEM_0830_0831_single_V36_results.xlsx
    工作表：inventory、image_status、all_objects、image_metrics、sample_summary、
            rejected_objects、via_model_choice、errors、settings。
  不生成前后配对表、变化百分比或差值图。

compare 额外输出：
  image_pair_comparison.csv                 图片级 before/after 统计与差值
  overall_comparison.csv                    同条件各区域合并后统计与差值
  object_before_after_comparison.csv        逐对象空间配对结果；不改变实际检测值
  SEM_0830_0831_before_after_V36_results.xlsx
    含通用工作表，以及 image_pair_compare、overall_compare、object_compare。
    前后差值是 after − before；change_percent 是差值 / before × 100%，
    缺失或 before 接近 0 时为空。逐对象匹配只是位置对应推断，需结合标注复核。

字段约定：
- condition 始终保存目录原始编号，single 中即样品编号；stage 为 single、before 或 after。
- pair_key 在 single 中仅为编号/Region/图案的分组标识，不意味着建立了前后配对。
- selected_v25_pattern_key、内部 _v32_* 等名称为继承兼容字段/函数，仍实际运行 V36。
- circularity、via_equivalent_diameter_nm、via_image_median_diameter_nm、
  via_size_deviation_fraction、via_size_reference_count、via_size_filter_applied 用于复核 via 筛选。
- input_width_tolerance / input_length_tolerance 是最终生效值；input_size_tolerance
  是是否显式指定旧统一参数的记录。edge_* 用于复核原图边缘证据。
- via 新增 via_convexity、via_ellipse_residual_p95、via_edge_roughness、via_shape_rule；
  via_radial_cv 在 V36 中为椭圆归一化后的半径变异系数，不与旧版原始半径 CV 直接比较。
- slit 新增 slot_observed_rows / slot_missing_rows、slot_is_partial_column、
  slot_supported_column_count。slot_complete_column_count 表示实际四排齐全的列数，
  可以为 0 而仍有可测的部分列。object_id / slot_column_id_v23 可有间隔，表示阵列列号。
- image_status 中 measure_v36_slot_* 记录全图支持情况、各排行数、缺失格位、缺失最下排、
  部分列测量数量；measure_v36_observed_component_count 记录 via 整体轮廓补充数量。
  全阵列不足时返回已有的拒测状态；查看 slot_insufficient_array_support 原因。
- settings.json 保存 via_shape_limits、via_shape_rule、slot_array_support_rule、
  slot_missing_site_policy 及 slit_min_column_rows 等实际规则。
- CSV 使用 UTF-8 BOM。Excel 直接打开 CSV 可能把 007 变成 7，建议打开本程序的 XLSX，
  或导入 CSV 时把 condition 列设为文本；程序导出的 XLSX 保留编号字符串。
- single 无有效测量时仍输出有表头的空明细/统计表，结合状态和错误表判断原因。

建议复核顺序：inventory → image_status/errors → annotated → rejected_objects → 统计表/plots。
标注中的测量对象与统计表一致；辅助阵列轮廓不一定参与统计。异常补检结果也应核对原图。

七、常见问题、失败处理和验证范围
------------------------------
1. 找不到数字子目录：确认 --input 指向样品父目录，目录名称为纯数字，图片没有额外嵌套。
2. 类型不对或结果异常：先核对自然排序；重命名图片时，同名 TXT 也要对应改名。
3. PixelSize 错误：检查 TXT 是否同名、是否有 PixelSize=正数、单位是否为 nm/pixel。
   single 只跳过对应图片，不影响其他样品。compare 跳过对应前后组。
4. 结果为 CHECK_NO_VALID_MEASUREMENT...：图片已处理，但原图证据或质量要求不足，
   不会编造尺寸。查看标注、剔除日志；核实类型、比例尺及尺寸提示，必要时换清晰原图。
   V36 形貌/缺失相关原因：via_irregular_boundary、via_crosses_lattice_cell、
   slot_missing_or_obscured_site、slot_insufficient_column_support、slot_insufficient_array_support。
   显式启用旧筛选时还可能出现 via_circularity_below_threshold、via_size_outlier、
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

V36 自动测试直接加载 V36，覆盖旋转椭圆原图测量、缺孔/星形干扰、长孔去重、
平滑大椭圆与圆孔混合、可选严格阈值、锯齿/多瓣排除、空白原图拒测、slit 局部缺失、
缺少最下排/整排时不补造、非阵列污渍及完整列模式。真实 CLI 验证污染合成图的
CSV/Excel/标注/统计一致性；同时继续验证 V35 的两种运行模式与原有功能。
仓库没有本批真实 SEM 原图及对应 TXT，合成图验证不代表真实样品识别准确率已确认。
请结合实际原图复核，严重模糊时尺寸提示不能保证恢复可靠测量。

开发测试（需完整仓库中的 tests/，发布 ZIP 不包含测试数据）：
  MPLBACKEND=Agg python3 -m unittest discover -s tests -v
算法、数据流和兼容说明见同目录 ALGORITHM_V36.md。
