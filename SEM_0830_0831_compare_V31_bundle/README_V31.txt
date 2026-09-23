V31 — V27 定位/输出 + Yakun 灰度边缘 + V30 有限补检
===================================================

运行
----
Python 3.9 或更新版本，安装：
  pip install numpy pandas scipy opencv-python matplotlib openpyxl

  python sem_before_after_compare_0831_V31.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V31"

脚本为独立文件，无需同时安装或运行 Yakun、V27、V30。
如果明确知道 via 类型，可加 --via-pattern via40 或 --via-pattern via60。

输入与 V30 相同
---------------
自动读取 after 下的纯数字子文件夹，按相同名称匹配 before，保留前导零。
支持 TIF/TIFF/PNG/JPG/JPEG；对应文件夹的图片数量须相同，不再固定为 9 张。
图片按文件名自然排序，每 3 张依次为 trench、slot、via，从 Region 2 开始。
每张图片必须有同名 TXT，其中包含 PixelSize=<数值>，单位 nm/pixel。
处理前请检查 inventory_mapping.csv 是否符合实际拍摄顺序。

V31 的主要改变
---------------
1. 以用户上传的 V27 为定位/阵列几何基准，恢复四行 slot 的定位方法。
   已通过原图边缘检查的列不再被 V30 补检轮廓替换。
2. 借鉴并内置 Yakun 的带方向梯度、两侧灰度参考、50% 亚像素交点和闭环追踪。
   slot/via 沿实际轮廓法线找边，不把目标强行拟合成矩形或椭圆。
   局部背景平面校正减少亮度渐变影响；增强图仅辅助定位，最终轮廓检查使用原图。
   8 位原图保留原灰度，高位深图仅作浮点线性换算，不裁剪两端灰度、不再次量化。
   短缺口连接相邻的有效边缘点，避免原候选上的假尖角撑大尺寸。
3. V30 的背景校正补检仅用于缺列/缺点或未通过原图边缘检查的候选。
   补检得到的轮廓仍须经过原图边缘、形状、位置检查。
4. slot 仍要求同一列的 4 行均有效，仅统计最下排 X/Y 尺寸。
   via 保留阵列/尺寸/形状检查。trench 保留 V27 的上下矩形宽度平均值、轴向 tip gap 定义。
5. 单张图片最长 300 秒，包含启动、补检和自动 via40/via60 两次尝试。
   到时终止该图片的独立进程，记录 SKIPPED_TIMEOUT / IMAGE_TIMEOUT，继续后续图片。
   未完成的临时标注不会发布；进程终止清理可能另用几秒。
   运行较久时每 15 秒显示已用时间。--image-timeout-seconds 120 可缩短上限，不能超过 300。
   auto 模式下若 BEFORE via 超时，配对 AFTER 记为 SKIPPED_VIA_MODEL_UNAVAILABLE，
   后续区域继续处理；这样不会用不同的 via 型号做前后比较。

输出参照 V27
-------------
  SEM_0830_0831_before_after_V31_results.xlsx
  all_object_measurements.csv
  image_metric_summary.csv
  image_pair_comparison.csv
  overall_comparison.csv
  object_before_after_comparison.csv
  inventory_mapping.csv
  image_status.csv
  rejected_objects.csv
  via_v25_model_choice.csv（保留原文件名，实际比较 V31 的两个 via 模型）
  processing_errors.csv
  settings.json
  annotated/、plots/

Excel 的表页和主要尺寸字段保持 V27 形式。
slot 上三排为青色轮廓和列.行号，最下排绿色 S 编号及 X×Y；via 保留 V 编号与底排颜色；
trench 保留矩形、轴点、宽度与 tip gap 箭头。最终标注尺寸来自最终测量轮廓。

新增的 edge_* 字段记录边缘支持率、局部对比度、位移和测量方法。
v27_x_width_nm / v27_y_height_nm / v27_tip_gap_y_nm 是细化前候选的参考值，
不是另一次完整 V27 运行；v31_locator 表示使用的候选来源。
边缘检查不通过的对象记录在 rejected_objects.csv，不进入统计；
同一图片的早期候选被拒绝后可能被补检救回，因此拒测明细不等于最终失败对象数。

验证范围
---------
已加入有已知 50% 灰度边界的合成图测试、V27 输出兼容检查，以及真实子进程超时测试。
椭圆、圆角矩形及旋转目标用于检查形状和尺寸；低对比、背景渐变、弱最下排用于补检验证。
这些是数值与流程验证。仓库未提供用户报告误测的原始 SEM 图片及 PixelSize TXT，
因此不能据此保证真实 SEM 的绝对测量精度；50% 灰度边界也不等同于材料边界的物理标定。
复核时请同时查看原图、annotated 标注和 edge_* 信息。

算法细节和来源校验值见 ALGORITHM_V31.md。
