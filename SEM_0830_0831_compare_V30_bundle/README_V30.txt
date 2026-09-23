SEM Before/After Compare V30
===========================

主程序：sem_before_after_compare_0831_V30.py
单文件独立运行，不需要 V29 文件或其他检测脚本。

V30 修复内容
------------
1. Slot 新增多尺度局部背景校正候选搜索，减少亮度渐变、对比度较低时的漏检。
2. 在阵列位置附近，用 ROI 边缘估计背景平面，再从背景校正图提取完整轮廓。
   旧方法可能反复找不到结构，也可能只检测到结构的一部分；V30 也会对已有
   候选重新检查局部轮廓，避免直接把截断的局部候选作为宽度。
3. 有足够候选支持的行，使用实际候选中心校正行位置，减少最下排预测偏移。
4. 保留四排、完整列、面积、形状和边界完整性检查；只测完整列中的最下排。
   无法通过检查的结构仍不强制输出测量值。
5. 保留上一轮四排拟合性能修复，并显示实际行位置、行距、补检尝试数、
   新增结构数、每排数量，以及局部对比不足/触及ROI边界/偏离位置/形状等诊断。
   attempts=90 表示尝试了 90 次，并不是已经找到了 90 个 slot。
6. 每张图片默认最多计算 300 秒（5 分钟）。使用独立子进程运行，超时后终止
   该图、丢弃未完成结果，记录 SKIPPED_TIMEOUT / IMAGE_TIMEOUT 并继续后续图片。
   限时涵盖该图的读取、分割、所有补检、自动 fallback、标注，以及 BEFORE via
   的 via40/via60 两次模型试算；两种模型不会各自重置为 5 分钟。
   via 复用胜出试算结果，不再重复测量 BEFORE。
   若 auto 模式的 BEFORE via 超时或失败、无法确定型号，对应 AFTER 会记为
   SKIPPED_VIA_MODEL_UNAVAILABLE，继续下一组；可在已知型号时使用 --via-pattern。
7. image_status.csv 和 processing_errors.csv 在处理期间持续更新。
   超时不记为尺寸 0，缺失一侧的图片不会生成有效前后尺寸差。

Trench 和 Via 沿用内嵌 V25 检测/测量逻辑。Slot 的检测和局部轮廓提取升级为 V30，
测量仍使用完整四排列的最下排轮廓在拟合坐标轴上的 X/Y 跨度。
source_algorithm 对 slot 标记为 V30_slot_recovery，其余为 V25_auto_fallback_only。

输入及配对规则
--------------
沿用 V29：只扫描 AFTER 根目录下纯数字名称的直接子文件夹，按编号排序，
在 BEFORE 找完全同名的文件夹。保留前导零，如 007 只匹配 007。
支持 tif/tiff/png/jpg/jpeg（大小写不限、可混用），每个 condition 的前后数量
必须相同，数量不限。每边分别按文件名自然顺序排列，按位置前后配对。
每张图片必须有同名 TXT，包含 PixelSize=<数字>（nm/pixel）。
图案按 trench、slot、via 循环；第 1-3 张为 Region 2，第 4-6 张为 Region 3，
继续递增。缺少文件夹、数量不同或 PixelSize 有问题时，输出诊断。

安装与运行
----------
pip install numpy pandas scipy opencv-python matplotlib openpyxl

python sem_before_after_compare_0831_V30.py

指定目录：
python sem_before_after_compare_0831_V30.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V30"

已知 via 型号时：
python sem_before_after_compare_0831_V30.py --via-pattern via40
python sem_before_after_compare_0831_V30.py --via-pattern via60

将单图时限缩短至 120 秒（只允许大于 0、且不超过 300 秒）：
python sem_before_after_compare_0831_V30.py --image-timeout-seconds 120

默认输入：
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\830SEM_before_treat
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\0831SEM_10-14_topview_after

默认输出：AFTER 根目录同级的 SEM_0830_0831_compare_V30。
请先停止旧进程，再运行 V30 脚本；运行中的 V29 不会自动升级。

主要输出
--------
SEM_0830_0831_before_after_V30_results.xlsx
inventory_mapping.csv
image_status.csv
all_object_measurements.csv
image_metric_summary.csv
image_pair_comparison.csv
overall_comparison.csv
object_before_after_comparison.csv
rejected_objects.csv
via_v25_model_choice.csv
processing_errors.csv
settings.json
annotated/ （via auto 试算图位于对应图片目录下的 _via_trials/）
plots/

验证范围
--------
回归测试包含亮度梯度、低对比、最下排比上三排更弱、行偏移、亮结构极性、
纯背景不误检，以及真实子进程强制超时、超时后继续与 via 共用单图任务。
测试时用较短时限验证超时行为，无须每次等待完整 5 分钟。
尚未获得用户出现漏检的原始 SEM 图，不能据合成图测试保证该原图一定检出。
如仍无结果，请结合标注图、slot grid 行位置及 diagnostics 检查失败阶段。
