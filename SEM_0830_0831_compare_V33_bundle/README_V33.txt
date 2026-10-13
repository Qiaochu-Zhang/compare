V33 — via 圆形度 / 同图尺寸筛选，slit / via / trench 尺寸辅助识别
================================================================

安装与运行
----------
Python 3.9 或更新版本，安装本目录 requirements.txt：
  pip install -r requirements.txt

  python sem_before_after_compare_0831_V33.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V33"

脚本独立运行，不需要 V32 或其他旧版本文件。建议使用新的输出目录。
默认 via60，沿用四排 slot / 最下排测量、原图边缘复核和单图 300 秒上限。

1. 仅针对 via 的质量筛选
------------------------
  --via-min-circularity 0.75
    最终原图轮廓的圆形度下限。圆形度 = 4π × 面积 / 周长²，范围 0–1，
    越接近 1 越圆；低于阈值的 via 不进入有效标注和统计。
    默认 0.75；例如 0.85 更严格，0.65 更宽松。设 0 关闭这项圆形度限制，
    原图证据、完整性、solidity、长宽比等其他形状检查仍保留。
    此参数不会筛掉长条 slit 或 trench。补检也不能降低这个用户阈值。

  --via-max-size-deviation 0.30
    等效直径 D = 2√(面积/π)。先取得同一张图内通过边缘、圆形度及所填尺寸
    范围检查的 via，再一次性计算 D 的中位数 M；|D/M - 1| 超过阈值则剔除。
    默认 0.30，即允许相对中位直径 ±30%；设 0 关闭，取值须小于 1。
    至少需要 3 个有效孔；不足时跳过同图尺寸筛选，并在结果中记录。
    BEFORE 和 AFTER 各自计算，绝不按前后尺寸差异删除真实变化。

示例（仅调 via 筛选）：
  python sem_before_after_compare_0831_V33.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V33" --via-min-circularity 0.85 --via-max-size-deviation 0.20

2. 可选的大致尺寸，单位均为 nm
------------------------------
  --via-diameter-nm        预计 via 直径
  --slit-width-nm          预计 slit 宽度（--slot-width-nm 是同义写法）
  --slit-length-nm         预计最下排 slit 长度（--slot-length-nm 是同义写法）
  --trench-width-nm        预计 trench 宽度
  --trench-gap-nm          预计上下 trench 端点间隙，不是横向阵列间距
  --size-tolerance 0.50    所填尺寸的相对容差，默认 ±50%，须大于 0 且小于 1

只填知道的项目，其余省略。尺寸必须是有限正数，不知道时不要填 0。
程序用每张图自己的 PixelSize 将 nm 转成像素：
  - via 直径决定圆孔搜索尺度和候选模板尺度范围。
  - slit 宽度决定去背景窗口、候选尺寸范围与缺失位置搜索；长度辅助平滑。
    四排 slit 的长度不同，长度范围只检查最下排测量对象，其他三排独立识别。
  - trench 宽度决定垂直暗条响应、峰间隔和边缘搜索窗口；间隙用于选择成对端点。
  - 最终测量还须落在已提供尺寸的容差范围内；没有填写的尺寸维度不加此限制。

例如预计 via 为 80 nm、容差 0.50，则最终等效直径允许 40–120 nm。
如果真实前后变化超出这个范围，应放宽容差、调整尺寸或省略该尺寸提示。
以下只是参数写法示例，数值不是本批图片的推荐值，请替换为自己的预计尺寸：
  python sem_before_after_compare_0831_V33.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V33" --via-diameter-nm 80 --slit-width-nm 80 --slit-length-nm 210 --trench-width-nm 160 --trench-gap-nm 120 --size-tolerance 0.50

尺寸提示用于搜索和筛选，不会把输出强制改成输入尺寸，也不会补画缺乏边缘证据的对象。
所有参数统一用于本次运行中的所有条件及前后图片；条件尺寸差别很大时，可分批运行。

3. 输入与输出
--------------
与 V32 一致：自动发现 AFTER 下的纯数字子目录，与 BEFORE 同名目录匹配，保留前导零。
支持 TIF/TIFF/PNG/JPG/JPEG，每张图片须有同名 TXT，内容包含 PixelSize=<数值>，单位 nm/pixel。
对应目录图片数量相同，按文件名自然排序；每 3 张依次为 trench、slot、via，Region 从 2 开始。
先检查 inventory_mapping.csv 是否符合实际拍摄顺序。

主要输出：
  SEM_0830_0831_before_after_V33_results.xlsx
  all_object_measurements.csv / image_metric_summary.csv
  image_pair_comparison.csv / overall_comparison.csv / object_before_after_comparison.csv
  rejected_objects.csv / image_status.csv / processing_errors.csv / settings.json
  via_model_choice.csv / inventory_mapping.csv / annotated/ / plots/

复核字段：
  circularity                           最终 via 圆形度
  via_equivalent_diameter_nm             最终 via 等效直径
  via_image_median_diameter_nm            同图参考中位直径
  via_size_deviation_fraction            相对中位直径偏差
  via_size_reference_count               参考孔数
  via_size_filter_applied                是否实际执行同图尺寸筛选
  input_*                               此次输入的筛选和尺寸参数
  via_search_scale_source                user_size_hint / candidate_pool_median / nominal_fallback
  edge_*                                原图边缘复核结果

rejected_objects.csv 中给出剔除原因及相关数值，例如：
  via_circularity_below_threshold / via_size_outlier / via_diameter_outside_hint_range
  slit_size_outside_hint_range / trench_size_outside_hint_range
settings.json 的 detection_options 保存完整参数。
被剔除的对象不进入有效统计；若启用拒绝标注，它们只作为拒绝对象显示。

4. 验证范围
------------
测试覆盖参数校验、真实子进程传参、同图尺寸异常剔除、补检不放宽用户阈值、
模糊 slit 尺寸辅助搜索、trench 宽度/间隙、空白拒测、标注及统计一致性。
继承测试继续检查原灰度精度、四排完整性、独立进程超时和处理后续图片。
仓库没有本批真实 SEM 原图及对应 PixelSize TXT；合成图通过不代表本批图片已完成准确率验证。
严重模糊、无可用边缘或尺寸提示偏差过大的图片仍可能被拒绝测量，应复核标注。

详见 ALGORITHM_V33.md。查看所有参数：
  python sem_before_after_compare_0831_V33.py --help
