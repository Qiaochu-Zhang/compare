V41 完整使用说明 — slit 四排长度递增，最下排至少 100 nm
======================================================

V41 基于 V40，加入用户确认的 slit 几何规则，默认自动执行，无新增参数或依赖：

  图像从上到下：第 1 排 < 第 2 排 < 第 3 排 < 第 4 排的 Y 向实测长度
  仅第 4 排：Y 向实测长度 >= 100 nm

上三排可以短于 100 nm。行号由图像上下位置确定，不能按长度重排行号。
规则作用于原图实测轮廓和每张图的 PixelSize，亮边与普通内部对比两种定位均适用。
缺失行只比较实际可见对象，不推算缺失尺寸；仍须满足原有四排行级与部分列支持要求。
倒置、长度相等或最下排不足 100 nm 的列不参与统计；补检和长度容差不能绕过此规则。

继续保留 V40 外圈亮边定位/内侧半高测量、V39 零结果恢复、single/compare、六项 nm
指标和单图最多 300 秒硬时限。via 与 trench 的识别测量规则沿用 V40。
本包四文件独立运行，不依赖历史版本；请运行 V41 脚本并换一个新的输出目录。

一、安装与首次运行
------------------
1. 解压 SEM_0830_0831_compare_V41_bundle.zip，进入同名目录。
   应能看到以下四个文件：
     sem_before_after_compare_0831_V41.py
     requirements.txt
     README_V41.txt
     ALGORITHM_V41.md
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
     python sem_before_after_compare_0831_V41.py --input "D:/SEM/samples" --output "D:/SEM/single_V41"
   macOS/Linux 示例：
     python3 sem_before_after_compare_0831_V41.py --input "/data/SEM/samples" --output "/data/SEM/single_V41"
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
  python sem_before_after_compare_0831_V41.py --input "D:/SEM/samples" --output "D:/SEM/single_V41"

显式选择 single，效果相同：
  python sem_before_after_compare_0831_V41.py --mode single --input "D:/SEM/samples" --output "D:/SEM/single_V41_run2"

原有 before/after 对比（不加 --mode 仍可使用）：
  python sem_before_after_compare_0831_V41.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V41"
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
- --output 可指定新的或空的目录；V41 拒绝非空目录，避免旧图或另一模式的结果混入。
  重新运行请换目录名，或自行整理旧输出后使用空目录。
- 不指定时，single 在输入根目录同级创建 SEM_0830_0831_single_V41；
  compare 在 after 根目录同级创建 SEM_0830_0831_compare_V41。

四、测量、统计与作图规则
----------------------
两种模式保留相同的六个尺寸指标，单位为 nm：
  trench：trench_width_nm（上下条平均宽度）、tip_to_tip_y_nm（上下端点 Y 间隙）
  slot：slot_x_width_nm（最下排 X 宽度）、slot_y_length_nm（最下排 Y 长度）
  via：via_x_width_nm（X 宽度）、via_y_height_nm（Y 高度）
对象明细中的源字段分别为 x_width_nm、y_height_nm、tip_gap_y_nm。
slit 和 slot 表示同类结构，输出沿用 slot 名称。

V41 自动执行固定规则：四排按图像位置从上到下编号 1、2、3、4，Y 向实测长度须
L1 < L2 < L3 < L4，且 L4 >= 100 nm。100 nm 只约束最下排，上三排可以短于 100 nm。
每列比较实际可见行，缺失行不推算；全图四个实际行级的长度中位数也必须递增。
违反规则的整列在阵列支持和统计前剔除，不按长度调换行号、不拉长轮廓、不补造对象。
长度沿测量 Y 轴由原图轮廓投影得到，再乘该图的 PixelSize；100 nm 不是 100 像素。
尺寸提示、100% 长度容差、零结果补检均不能放宽此固定规则，无须增加命令参数。


- slot 只统计最下排。四排均有真实轮廓的完整列可以独立保留；部分列必须有全图
  四排行级的实际观测、至少两个受支持列，且至少三排各有两个对象。测量部分列
  默认至少两排行级的边缘支持，须含最下排；没有最下排轮廓时不测量，不提升第三排。
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
- 为避免拉长造成面积/等效直径变化后被另一项默认筛选删除，V41 沿用 V38 的同图尺寸
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
- 每个有效轮廓仍须满足原灰度边缘支持率至少 85%、连续缺边段不超过约 8%。
  另外检查外侧是否有未纳入的前景；“局部支持率 100%”也可能只是内部台阶，不能
  单独证明轮廓完整。V41 沿用 V38，先重找完整定位轮廓，再用同一原图边缘审核。
  少量缺边只允许在真实支持端点间连接，大面积遮挡不会凭模型补出尺寸。
- 如确需恢复圆形度门槛，可显式设 --via-min-circularity 0.75；该设置可能删除椭圆。
  同图偏差筛选显式启用后，从通过形状/提示检查的 via 计算一次中位等效直径 M，
  剔除 |D/M−1| 超限对象；D=2√(面积/π)。不足 3 个有效孔时跳过同图筛选。
  single 每张图独立处理，compare 的 before/after 也各自独立，不按前后差异剔除。
- auto 会试 via40/via60，并复用较可靠的测量与标注。single 每张 via 图片独立选模型；
  compare 由 BEFORE 选并固定给 AFTER。single 不启用历史 AFTER 专用增强。

V40 对常规内部对比定位继续沿用以下完整轮廓复核：
- slit 的首轮候选按 V35 的定位顺序处理；已有通过原图验证的对象不会被增强图案
  替换。补检另跑四排定位，首轮四排行级齐备时，补检对象须落在这些实际行级附近。
- 轮廓验证失败时，在同一格位用原始灰度重新定位，再复核原图边缘、中心偏移和
  图像边界。修复不会直接报告定位掩膜尺寸，也不能凭空白原图产生测量。
- slit 完整性用横向局部背景消除慢变照明；外探距离为 max(8 px,0.5×较小跨度)，
  从边缘之外 2 px 到该距离均须有连续前景。中间有背景间隔的邻位不算目标延伸。
  这是检测证据，不是按该距离延长尺寸。大于该尺度的模糊或极短弱尾部仍需复核。
- 实际弱下半部仍可先重定位完整区域，再执行原图 50% 灰度法线追踪；真正短的
  对象维持原长度。支持率至少 85%、连续缺边约 8% 上限和形状检查继续执行。
- via 优先使用显式尺寸提示；否则从首轮重复近圆对象估计孔径（V35），不足时
  使用首轮整体椭圆（V37），仍不足才从补检候选估计。径向定位使用首轮候选，
  整体轮廓补检继续使用补充候选。可靠的首轮尺度不会被孔间碎片改写。
- via 用孔左/右外侧原图灰度逐行插值背景，以消除弯曲照明；仅将两侧采样范围之内
  的响应当作证据。背景样本不足时退回原局部背景估计。外侧需要连续前景，
  必要时仍进行区域重定位。最终仍由原图灰度边缘测量，不用背景模型代替孔形状。
- 多尺度弱孔搜索、椭圆形貌、格位/最终去重和显式筛选参数保持。V37_* 等继承
  方法/来源字段按实际算法来源保留。
- 重定位窗口保留背景留白，候选需覆盖至少 60% 种子核心，定位面积为原种子的
  0.65–4.5 倍，不能触窗口边缘；无可靠完整轮廓时记录 UNRESOLVED_OUTER_BOUNDARY。
- trench 沿用 V37：上下条各自定位条体与端点，原图平台 50% 交点决定尺寸；
  不优先选择更强的内部明暗台阶，不把全图中位端点复制到所有列。

slit 污染缺失处理详解：
- 四排是实际观测约束。完整列四排齐全时可保留；默认部分列有“最下排 + 任意
  一个以上其他排”的真实轮廓支持即可，但部分列还须满足第四节的重复阵列证据。
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
  python sem_before_after_compare_0831_V41.py --input "D:/SEM/samples" --output "D:/SEM/single_V41_hints" --via-diameter-nm 80 --slit-width-nm 80 --slit-length-nm 210 --trench-width-nm 160 --trench-gap-nm 120

需要主动收紧圆形度/尺寸筛选时（会排除部分椭圆），可用：
  python sem_before_after_compare_0831_V41.py --input "D:/SEM/samples" --output "D:/SEM/single_V41_strict" --via-min-circularity 0.85 --via-max-size-deviation 0.20 --image-timeout-seconds 120

有污染缺失的 slit 默认无需额外参数；要切回每列完整检查：
  python sem_before_after_compare_0831_V41.py --input "D:/SEM/samples" --output "D:/SEM/single_V41_full_columns" --slit-min-column-rows 4

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
  SEM_0830_0831_single_V41_results.xlsx
    工作表：inventory、image_status、all_objects、image_metrics、sample_summary、
            rejected_objects、via_model_choice、errors、settings。
  不生成前后配对表、变化百分比或差值图。

compare 额外输出：
  image_pair_comparison.csv                 图片级 before/after 统计与差值
  overall_comparison.csv                    同条件各区域合并后统计与差值
  object_before_after_comparison.csv        逐对象空间配对结果；不改变实际检测值
  SEM_0830_0831_before_after_V41_results.xlsx
    含通用工作表，以及 image_pair_compare、overall_compare、object_compare。
    前后差值是 after − before；change_percent 是差值 / before × 100%，
    缺失或 before 接近 0 时为空。逐对象匹配只是位置对应推断，需结合标注复核。

V41 长度规则字段：
- settings.json 的 v41_slit_length_rule 保存 minimum_scope=bottom_row_only、
  bottom_min_y_length_nm=100、图像行顺序和严格递增长度规则。
- 对象明细的 v41_slit_length_rule_passed、v41_slit_y_lengths_nm_by_row、
  v41_slit_bottom_min_y_length_nm、v41_slit_y_length_order、v41_slit_minimum_scope
  记录本列实际观测长度和审核结果，缺失行不填数。
- image_status 的 measure_v41_slit_length_* / measure_v41_slit_global_length_order_*
  记录规则启用、因下限/倒置/全图行序被剔除的列数；
  measure_v41_slit_row_median_y_lengths_nm 记录四排行级实际中位数。
- rejected_objects.csv 的 slit_bottom_y_length_below_100nm / slit_y_length_not_increasing /
  slit_column_length_rule_failed / slit_row_y_lengths_not_increasing 区分具体原因。
- 可靠亮圈附近的亮边材料条带不作为另一条较长 slit；明细仍以实际孔内边缘为准。
  measure_v41_slit_rim_band_rejected_count / slit_candidate_is_bright_rim_band 记录此类剔除。
- 最下排实际结果可以大于 100 nm，也可以因无法满足审核而缺失；不会固定输出 100 nm。

继承的 V40 亮边字段与复核：
- all_object_measurements.csv 的 v40_bright_rim_used 表示是否用亮边内侧测量；
  edge_method=V40_bright_rim_inner_gray50_normal_tracking，原图支持字段 edge_* 同时保留。
- v40_rim_boundary_definition 指明测量口径是内部灰度与亮峰之间的内侧 50% 交点；
  v40_rim_prominence_gray 是亮峰超过内/外两侧参考的显著度，
  v40_rim_peak_offset_px 是内侧交点到峰值的距离，不能当作完整亮边厚度。
- slit_locator_path / via_locator_path 可为 V40_bright_rim；亮边测量的
  slit_boundary_rule / via_boundary_rule 为 V40_original_bright_rim_inner_slope。
- 亮边测量的 confidence_flag 和 measurement_confidence_flag 均为 BRIGHT_RIM_REVIEW，
  优先于继承的零结果复核标签；auto_fallback_* / v39_zero_recovery_* 仍记录是否经历补检。
- image_status 的 v40_bright_rim_measured_count 是本图最终亮边测量数；
  measure_v40_rim_candidate_count / examined_count 是定位候选/审核数，不是测量数。
  slit 另有 measure_v40_slit_rim_observed_count / added_count / grid_locked /
  grid_rejected_count；via 另有 measure_v40_via_rim_measured_count。
- settings 保存 v40_bright_rim_enabled / rule / limits / review_flag。
  slot_outside_bright_rim_grid 表示对象与充分观测的亮圈阵列冲突；此时直接剔除，不平移轮廓。

继承的 V39 零结果补检字段：
- image_status 的 v39_zero_recovery_attempted / used 表示是否尝试/采用继承的零结果补检。
- v39_zero_recovery_attempts / candidate_counts 记录尝试方法及各次候选数量；
  via/slit 的数量为已通过原图初审的定位候选，trench 的数量为有效配对数。
- v39_zero_recovery_method 为 V39_multiscale_structure 或 V39_opposite_polarity；
  v39_feature_polarity 为常规内部对比搜索使用的 dark / bright；亮边始终寻找原图亮峰，另记 v40_rim_polarity=bright。
- 对象明细记录 v39_zero_recovery_used / method、v39_feature_polarity；采用新增补检时
  confidence_flag=ZERO_RESULT_RECOVERY_REVIEW，selection_mode 为实际方法名。
  原 measurement_confidence_flag 保留 RELAXED_FALLBACK_REVIEW 兼容值；新增字段
  进一步区分继承的搜索恢复，不表示继续放宽原图边缘或形貌要求。
- auto_fallback_level_used=3 表示原极性多尺度/结构补检，4 表示相反极性补检；
  状态沿用 OK_AUTO_FALLBACK_L3 / L4。原有 L1/L2 仍表示既有增强补检。
- measure_v39_via_sparse_bootstrap 表示最终完整保留了四个经原图、形貌和用户范围
  审核的 2×2 孔。少一个孔或其中一个最终未通过审核，该少量孔阵列仍拒测。
- settings 的 v39_zero_recovery_*、v39_zero_candidate_limits、v39_sparse_via_rule
  保存启用状态、触发规则、数量/尺度上限和四孔支撑条件。
- 补检成功的标注有 V41 ZERO RECOVERY (...): REVIEW 横幅；图像显示原始 SEM。

字段约定：
- condition 始终保存目录原始编号，single 中即样品编号；stage 为 single、before 或 after。
- pair_key 在 single 中仅为编号/Region/图案的分组标识，不意味着建立了前后配对。
- selected_v25_pattern_key、内部 _v32_* 等名称为继承兼容字段/函数，仍实际运行 V41。
- circularity、via_equivalent_diameter_nm、via_image_median_diameter_nm、
  via_size_deviation_fraction、via_size_reference_count、via_size_filter_applied 用于复核 via 筛选。
- input_width_tolerance / input_length_tolerance 是最终生效值；input_size_tolerance
  是是否显式指定旧统一参数的记录。edge_* 用于复核原图边缘证据。
- via 新增 via_convexity、via_ellipse_residual_p95、via_edge_roughness、via_shape_rule；
  via_radial_cv 在 V41 中为椭圆归一化后的半径变异系数，不与旧版原始半径 CV 直接比较。
- slit 新增 slot_observed_rows / slot_missing_rows、slot_is_partial_column、
  slot_supported_column_count。slot_complete_column_count 表示实际四排齐全的列数，
  可以为 0 而仍有可测的部分列。object_id / slot_column_id_v23 可有间隔，表示阵列列号。
- image_status 中 measure_v36_slot_* 记录全图支持情况、各排行数、缺失格位、缺失最下排、
  部分列测量数量；measure_v36_observed_component_count 记录 via 整体轮廓补充数量。
  部分列阵列不足时查看 slot_insufficient_array_support；四排齐全的列可单独保留。
- 完整边缘新字段：boundary_relocalized（是否使用区域重定位路径）、
  boundary_recovery_attempted、boundary_outer_contact_fraction（外侧未解释前景比例）、
  boundary_seed_outer_contact_fraction、boundary_search_roi_px（复核窗口）、
  boundary_seed_area_ratio、boundary_region_threshold_gray、boundary_recovery_method。
  非重定位对象的部分字段为空；重定位标记不等于测量准确率证明，应核对原图。
- image_status 的 measure_v37_boundary_relocalized_count / measure_v37_local_contrast_candidates
  记录相关计数；trench 的 measure_v37_trench_* 记录候选数、噪声和实际有效对数。
- settings.json 的 v38_slit_localization_rule / v38_slit_boundary_rule 保存 slit 新规则；
  v37_boundary_completeness_rule / v37_region_recovery_limits / v37_local_candidate_rule /
  v37_trench_rule 保留继承的规则和区域限制。
- 常规 slit_locator_path 标记 V35_primary（沿用首轮定位顺序）或 V37_supplement（补检），
  不是重新导入旧版本测量结果。slit_boundary_rule 标记 V38 局部连续外侧前景规则。
- image_status 的 measure_v38_slit_primary_accepted_count / supplement_added_count /
  local_repaired_count / primary_row_lock / complete_column_support_only 记录逐对象审核阶段
  的首轮、补充、修复数量及行级锁定、仅保留完整列状态；各项均带 measure_v38_slit_ 前缀。
  最终统计数量以 n_measurements 和对象明细为准，审核后仍可能因阵列/尺寸被剔除。
- 常规 via_boundary_rule 记录继承 V38 R2 的局部背景/连续性复核；via_locator_path 区分
  V35_primary_radial 与 V37_observed_component；via_scale_locator 区分 user_size_hint、
  V35_primary_near_round、V37_primary_whole_component、V37_supplemental_scale；
  从重复亮圈估计尺度时为 V40_bright_rim_scale，并记录 v40_via_rim_scale_used。
- image_status 的 measure_v38_via_primary_scale_locked / primary_scale_candidate_count /
  primary_radial_count（均带 measure_v38_via_ 前缀）记录首轮尺度是否锁定、支撑对象数
  以及首轮径向定位数。后者不是最终测量数；以 n_measurements 和明细为准。
- settings 的 v38_via_scale_rule / v38_via_boundary_rule 保存继承的融合规则；script_version
  为 2026-10-16-V41-SLIT-LENGTH-ORDER，source_algorithm 为
  V41_slit_length_order_on_V40_bright_rim。
- settings.json 还保存 via_shape_limits、via_shape_rule、slot_array_support_rule、
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
   V41 原图边缘不足/外轮廓仍截断时，可在 rejected_objects.csv 查看
   original_image_edge_not_supported 及 edge_status=UNRESOLVED_OUTER_BOUNDARY。
   形貌/缺失相关原因：via_irregular_boundary、via_crosses_lattice_cell、
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

V41 测试实际导入本版脚本，覆盖空间行序、递增/相等/倒置、仅最下排的 100 nm 下限、
每图 PixelSize、缺失行、普通暗槽/亮边槽、尺寸提示不可放宽规则，以及实际 single/compare、
CLI/worker、统计、Excel、作图、auto 和硬超时。
全仓库回归 369 项全部通过，1348.641 秒；其中 33 项直接加载 V41。
完整运行环境及 ZIP 一致性结果见仓库 md/2026-10-16_V41_slit-length-order.md。
V40 及更早版本保留历史行为；旧的从上到下变短的 slit 样例在 V41 中属于拒测反例。
仓库没有本批真实 SEM 原图及对应 TXT，合成验证不代表实际召回率或绝对尺寸准确度。
Linux/Python 3.12 为本次运行环境，Windows/macOS 原生运行未验证。

开发测试（需完整仓库中的 tests/，发布 ZIP 不包含测试数据）：
  MPLBACKEND=Agg python3 -m unittest discover -s tests -v
算法、数据流和兼容说明见同目录 ALGORITHM_V41.md。
