SEM Before/After Compare V29 - V25 only
========================================

主程序：
  sem_before_after_compare_0831_V29_V25only.py

特点
----
1. 单文件独立运行；V25算法已经完整包含在本文件中。
2. 不加载 sem_cd_measure_200k_batch_V1_6.py。
3. 自动扫描 after 根目录下所有以纯数字命名的直接子文件夹，按编号数值排序处理。
   文件夹名称由数字 0-9 组成，例如 0、1、2、007、25、100；不要求编号连续。
   以 after 中的文件夹名单为准，在 before 中查找完全同名的文件夹。
   before 独有的文件夹、非数字名称文件夹、数字名称的普通文件均不参与处理，
   也不会递归扫描更深层目录。
   保留前导零，例如 after/007 只匹配 before/007，不会自动匹配 before/7。
   before 缺少对应文件夹时记录错误并跳过该编号，继续处理其他有效编号。
   after 不存在、不是目录或没有任何纯数字子文件夹时，给出明确错误提示。
   例如 after 中只有 2、7、25 时，只比较 before/after 的这三个编号。
4. 支持 tif、tiff、png、jpg、jpeg，扩展名不区分大小写，可混合使用。
   每个目标子文件夹的图片数量不限，只要求 before/after 对应子文件夹内数量一致。
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
   不同 condition 可有不同照片数量，例如 2 文件夹为 8 张、25 文件夹为 12 张。
   对应两边数量不一致或缺少文件夹时，跳过该 condition 的前后两组并记录错误，
   其他数量一致的 condition 继续处理；两边都为空时无照片可处理。
   数量校验仅统计当前子文件夹内的上述图片格式，不计 TXT 等文件或嵌套文件夹。
   前后对应图片的格式可以不同，例如 before 为 PNG、after 为 JPG。
5. 每张图片独立读取同名 TXT 中 PixelSize=<number>。
   例如 image_1.png / image_1.jpg / image_1.tif 均对应 image_1.txt。
6. before/after 按 Condition + 自然排序位置配对，不要求文件名相同；
   Region + Pattern 由排序位置生成。两边照片应保持相同拍摄顺序及图案排列。
   结果表和跨文件夹汇总图也按编号数值排序，例如 2、7、25。
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
C:\Users\z00027644\Documents\倾斜刻蚀\SEM\SEM_0830_0831_compare_V29

安装依赖
--------
pip install numpy pandas scipy opencv-python matplotlib openpyxl

运行
----
python sem_before_after_compare_0831_V29_V25only.py

指定前后及输出目录：
python sem_before_after_compare_0831_V29_V25only.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_output"

强制 via40：
python sem_before_after_compare_0831_V29_V25only.py --via-pattern via40

强制 via60：
python sem_before_after_compare_0831_V29_V25only.py --via-pattern via60

主要输出
--------
SEM_0830_0831_before_after_V29_results.xlsx
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

Slot 长时间停在 fallback mask 日志（V29 性能修复）
------------------------------------------------
例如：fallback mask enh_adaptive: raw=13, added=13, pool=13
表示该阈值掩膜提取了 13 个候选、加入 13 个、当前候选池共有 13 个；
这不是最下排测量完成，也不是报错。后面仍需要拟合完整四排、补齐缺失结构，
通过完整四排列检查后，才统计最下排。

旧代码在形状有效候选少于 8 个或初次拟合失败时，还会追加最多 180 个
低灰度定位候选。旧四排拟合逐个枚举行距、起点并遍历全部候选，最坏计算量
接近候选数的四次方，而且这部分原先没有进度输出，可能看起来一直卡住。
单靠最后一行日志无法判定某张原图实际卡在哪一步。

本修复保留所有原有行距假设、评分、平分时的选择顺序及测量/QC标准：
1. 用行区间支持数的评分上界排除不可能胜出的模型，并分批计算距离。
2. 四排拟合开始时显示 candidates/pitches，长任务约每 5 秒显示进度和耗时。
3. 显示低灰度定位阶段、实际新增定位候选数和四排拟合候选数。
4. 低灰度阈值图逐张处理，不再同时保留 18 张全图掩膜；达到候选上限即停止。

本机合成候选测试：45 个候选的行拟合从约 31 秒降到 0.23 秒；
180 个随机候选约 13 秒。以上仅为行拟合耗时，原图整套处理耗时取决于
分辨率、噪声、缺失结构数、局部补检及 zero-result fallback 是否触发。

使用方式：停止旧进程，替换为本包的同名 Python 文件，按原命令重新运行。
运行中的旧进程不会自动加载修改。如果仍慢，保留新增日志中最后一个阶段、
candidates/pitches 数量与耗时，用于继续定位。
