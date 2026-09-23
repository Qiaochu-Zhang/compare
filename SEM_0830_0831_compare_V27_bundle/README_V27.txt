SEM Before/After Compare V27 - V25 only
========================================

主程序：
  sem_before_after_compare_0831_V27_V25only.py

特点
----
1. 单文件独立运行；V25算法已经完整包含在本文件中。
2. 不加载 sem_cd_measure_200k_batch_V1_6.py。
3. 只处理两个根目录下的 10、11、12、13、14 子文件夹。
4. 每个目标子文件夹的 tif/tiff 数量不限，只要求 before/after 对应子文件夹内数量一致。
   每边分别按文件名自然数字顺序排列（例如 1、2、10），按相同位置一一配对。
   图案排列沿用原规则，按 trench、slot、via 循环，Region 从 2 开始递增：
     1: Region 2 trench
     2: Region 2 slot
     3: Region 2 via
     4: Region 3 trench
     5: Region 3 slot
     6: Region 3 via
     7: Region 4 trench
     8: Region 4 slot
     9: Region 4 via
    10: Region 5 trench
    11: Region 5 slot
    12: Region 5 via
    ……后续依此类推；最后一组不足 3 张也可处理，不要求总数为 3 的倍数。
   不同 condition 可有不同照片数量，例如 10 文件夹为 8 张、11 文件夹为 12 张。
   对应两边数量不一致或缺少文件夹时，跳过该 condition 的前后两组并记录错误，
   其他数量一致的 condition 继续处理；两边都为空时无照片可处理。
   数量校验仅统计当前子文件夹内的 tif/tiff（扩展名不区分大小写），不计 TXT 等文件。
5. 每张 TIF 独立读取同名 TXT 中 PixelSize=<number>。
6. before/after 按 Condition + 自然排序位置配对，不要求文件名相同；
   Region + Pattern 由排序位置生成。两边照片应保持相同拍摄顺序及图案排列。
7. trench/slot/via 的识别、边缘、QC、点阵/列规则与 zero-result fallback 全部使用 V25。
8. 当前 1_xxxx 文件名没有 40/60 信息，因此 via 默认：
     - BEFORE 上分别运行 V25 via40 与 via60；
     - 按 V25 测量有效性/严格检测/点阵与边缘质量选一个；
     - 对应 AFTER 固定使用同一个 via key。
   如果明确知道类型，可强制：--via-pattern via40 或 --via-pattern via60。

默认输入
--------
Before:
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\830SEM_before_treat

After:
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\0831SEM_10-14_topview_after

默认输出
--------
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\SEM_0830_0831_compare_V27

安装依赖
--------
pip install numpy pandas scipy opencv-python matplotlib openpyxl

运行
----
python sem_before_after_compare_0831_V27_V25only.py

指定前后及输出目录：
python sem_before_after_compare_0831_V27_V25only.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_output"

强制 via40：
python sem_before_after_compare_0831_V27_V25only.py --via-pattern via40

强制 via60：
python sem_before_after_compare_0831_V27_V25only.py --via-pattern via60

主要输出
--------
SEM_0830_0831_before_after_V27_results.xlsx
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
annotated/
plots/

指标
----
Trench:
  trench_width_nm = V25 上下两个 edge-detected rectangle 宽度的算术平均
  tip_to_tip_y_nm = V25 定义的 Y 轴 tip-to-tip 距离

Slot:
  slot_x_width_nm = V25 完整四排列模型中最下排 slot 的 X 宽度
  slot_y_length_nm = 同一最下排 slot 的 Y 长度

Via:
  via_x_width_nm = V25 via X span
  via_y_height_nm = V25 via Y span

image_pair_comparison.csv 是最直接查看每张对应图 before/after 平均尺寸变化的表。
