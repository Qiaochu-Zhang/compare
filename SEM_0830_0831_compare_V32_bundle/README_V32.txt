V32 — 默认 via60，按图中孔径搜索，放宽相对 CD 筛选
=================================================

运行
----
Python 3.9 或更新版本，安装：
  pip install numpy pandas scipy opencv-python matplotlib openpyxl

  python sem_before_after_compare_0831_V32.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V32"

默认就是 --via-pattern via60，不需要额外设置；before/after 都使用 via60。
仅在确实需要时使用 --via-pattern via40 或 --via-pattern auto。
auto 在 BEFORE 上比较两个模型，并把选定型号固定给对应 AFTER。
脚本为独立文件，不需要 V31、V27、V30 或 Yakun 源码。

本次修正
--------
1. V31 默认 auto，可能选中 via40；这会改变尺度、平滑和筛选，并非仅名称不同。
   V32 默认固定 via60，文件名以 40 开头也不会覆盖这个设置。
2. via60 的 60 nm 只表示标称图案型号，不再作为 CD 的固定上限依据。
   从图中至少 4 个完整、近圆且尺寸相近的候选估计搜索孔径；before/after 分别估计。
   CD 仍由原图灰度边缘、PixelSize 和最终坐标轴计算，不强制等于 60 nm，也不放大已有结果。
   候选不足时回到标称搜索尺度，并在结果中明确记录。
3. via 的跨度筛选改为图中搜索孔径的 0.50–2.00 倍，补检为 0.45–2.10 / 0.40–2.20 倍。
   同时调整候选尺寸、等效直径和中间轮廓检查，避免只放宽最终输出却被上游过滤。
   这些是筛选范围，不是保证可识别的 CD 范围；邻孔间距、有效边缘和形状仍会限制检测。
4. slot210 没有 210 nm 的 CD 硬上限；相对平均面积范围从 0.1–10 倍放宽为 0.05–20 倍。
   低灰度/背景校正候选面积上限从图像面积的 3.5% 放宽为 10%。
5. trench160 没有 160 nm 的 CD 硬上限；图内宽度/间隙相对中位数容差放宽为 85%/80%，
   MAD 系数从 5 调到 7，补检同步放宽，避免补检反而用更严格的旧参数。
6. 保留 V31 原图边缘审核、形状和阵列检查、slot 完整四行/仅统计最下排、单图 300 秒时限。

输入与配对
----------
自动发现 after 下的纯数字子文件夹，按相同名称匹配 before，保留前导零。
支持 TIF/TIFF/PNG/JPG/JPEG；对应文件夹图片数量须相同。
按文件名自然排序，每 3 张依次为 trench、slot、via，从 Region 2 开始。
每张图片必须有同名 TXT，包含 PixelSize=<数值>，单位 nm/pixel。
处理前请检查 inventory_mapping.csv 是否符合实际拍摄顺序。

输出与复核
----------
  SEM_0830_0831_before_after_V32_results.xlsx
  all_object_measurements.csv / image_metric_summary.csv
  image_pair_comparison.csv / overall_comparison.csv / object_before_after_comparison.csv
  inventory_mapping.csv / image_status.csv / rejected_objects.csv
  via_model_choice.csv / processing_errors.csv / settings.json
  annotated/、plots/

重点查看：
  selected_pattern_key / source_pattern_key：实际模型，默认应为 via60。
  via_nominal_diameter_nm：型号标称值 60，不是实测 CD。
  via_search_diameter_nm / via_search_scale_source：图中估计的搜索孔径及来源。
  x_width_nm / y_height_nm：实测 X/Y CD。
  edge_*：原图边缘审核信息。
  settings.json 的 cd_filter_policy：本次运行的相对尺寸策略。

历史 selected_v25_pattern_key 字段保留兼容，含义同 selected_pattern_key。
V27/V24 等内部字段来自继承算法，不表示实际运行了另一个版本。
via_model_choice.csv 在默认模式只记录 via60；不会自动试 via40。

默认单图上限 300 秒；--image-timeout-seconds 120 可缩短，上限不能超过 300。
超时跳过当前图并继续，未完成标注不会发布。

验证与限制
----------
新增 65–240 nm 合成 via、较大像素孔径、真实独立 worker、超出 160/210 nm 的 trench/slot 测试。
沿用 V31 的原图边缘精度、四行约束、空白拒测、超时和输出检查。
仓库没有这批真实 SEM 原图及 PixelSize TXT，因此尚未验证这批图片的实际准确率。
请用新输出目录重跑，并复核标注轮廓；不要只把旧结果中的 via40 改名为 via60。

算法细节见 ALGORITHM_V32.md。
