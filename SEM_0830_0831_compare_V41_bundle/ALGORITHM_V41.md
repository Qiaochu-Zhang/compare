# V41：slit 四排向下增长，最下排至少 100 nm

V41 基于 V40。用户确认 **100 nm 下限只作用于最下排**，上三排允许短于 100 nm。新规则默认启用，普通内部对比与亮圈定位、正常测量和所有补检路径统一执行。无需新增参数或依赖。

## 四排的实际顺序

行号仍按原图位置自上而下为 1、2、3、4；不按长度重排行号，也不倒置图像来满足先验。Y 长度沿既有测量 Y 轴，由原图轮廓投影跨度乘每张图的 PixelSize 得到，单位 nm。

完整列要求：`L1 < L2 < L3 < L4` 且 `L4 >= 100 nm`。

部分列比较按实际行号排序的可见对象，例如只有第 2、4 排时要求 `L2 < L4` 且 `L4 >= 100 nm`。缺失行没有插值长度；缺最下排不能提升第 3 排来测量。无最下排的辅助列没有 100 nm 检查，也不产生最下排结果。

四排全图实际长度中位数也须从上到下递增，避免不同部分列拼出与四排趋势冲突的图案。只对已经实际观测到的四排行级执行此检查，不补造不存在的中位数。

## 审核位置与处理

- `_v41_slot_length_audit` 检查每列原图细化后的轮廓。所有可见长度必须有限且为正；较低行的实测长度严格大于较高行，同时实际中心必须符合图像上下顺序。最下排单独执行 100 nm 下限。
- 比较只使用 `1e-7 nm` 的浮点运算量级误差，不加入额外尺寸容差、固定增长比例或尺寸补偿。
- `_v41_filter_slot_length_columns` 在计算最终阵列支持前剔除违规整列，记录具体对象与原因；其他正确列可保留。剔除后重新拟合测量轴，并继续审核剩余列直到稳定，保证尺寸与最终轴一致。
- 充分观测的四排行级中位数不递增时，拒绝这组列。完整列/部分列支持门槛随后照常执行；因长度失败被剔除的列不贡献阵列支撑。
- 仅保留完整列的回退发生后，再按实际保留列重新拟合轴并审核长度；不会通过丢弃部分列而绕过新下限。
- 最下排最终 nm 长度同时接受显式 `--slit-length-nm` 和 `--length-tolerance` 检查。宽容差与已知长度提示不能放宽向下递增和 100 nm 固定下限。
- `_v41_exclude_slit_rim_bands` 在已观测到可靠亮圈时，剔除属于同一亮圈的普通候选并重新送入实际内沿。避免最下排孔内长度不足 100 nm 后，相反极性补检把较长的亮边材料条带当作另一个 slit。关联扩展距为 `clip(2×内沿到亮峰距离+3,3,16) px`；内沿仍接受原图和长度审核。
- 输出保持原图实测值，不重排行号、不伸长边缘、不填零、不将预计值或 100 nm 直接作为结果。

## 保留的定位与测量

V40 的亮圈内孔定位、原图法线“内侧上升—亮峰—外侧回落”和内侧 50% 交点继续工作；有效边缘支持至少 85%、连续缺边不超过 8%。传统暗槽/弱槽仍走 V35/V37/V38 融合定位与原图复核。

slit 继续要求四个实际行级。完整四排一列可保留；部分列仍须全图四排实际观测、至少两列且三排有重复支持，每个测量列包含最下排及用户要求的实际行数。只统计最下排，辅助轮廓不等于统计对象。

via 继续使用普通径向、整体轮廓与亮边内沿互补、原图形貌、实际二维阵列及格位去重；trench 继续按原图条体与端点测量。两者不执行 slit 的 100 nm 或向下增长规则。

正常及 L1/L2 无结果时，继承的 V39 L3 原极性多尺度和 L4 相反极性流程仍运行；这些流程调用同一最终 slit 测量入口，因此不绕过 V41 规则。所有尝试及 auto 两模型仍共享单图最多 300 秒硬时限。

## 输出与审计

- 入口：`sem_before_after_compare_0831_V41.py`、`main_v41`、`V41DetectionOptions`；目录与工作簿名称带 V41。
- `script_version=2026-10-16-V41-SLIT-LENGTH-ORDER`，`source_algorithm=V41_slit_length_order_on_V40_bright_rim`。
- settings：`v41_slit_length_rule` 记录下限 100 nm、`minimum_scope=bottom_row_only`、图像行顺序、严格递增和缺失行规则。
- 对象：`v41_slit_length_rule_passed`、`v41_slit_bottom_min_y_length_nm`、`v41_slit_y_lengths_nm_by_row`、`v41_slit_y_length_order`、`v41_slit_minimum_scope`；长度 JSON 仅含实际观测行，绝不为缺失行补值。
- image_status：`measure_v41_slit_length_rule_enabled`、`measure_v41_slit_length_rejected_columns`、`measure_v41_slit_length_floor_rejected_columns`、`measure_v41_slit_length_order_rejected_columns`、`measure_v41_slit_global_length_order_rejected_columns`、`measure_v41_slit_row_median_y_lengths_nm`。
- `measure_v41_slit_rim_band_rejected_count` 记录被可靠内沿替代的亮边材料候选数，原因码为 `slit_candidate_is_bright_rim_band`。
- 剔除原因：`slit_bottom_y_length_below_100nm`、`slit_y_length_invalid`、`slit_y_length_not_increasing`、`slit_column_length_rule_failed`、`slit_row_y_lengths_not_increasing`。
- `slot_column_id_v41` 是长度审核时的候选列编号；最终 `slot_column_id_v23` / object_id 继续表示实际保留阵列的列号。
- 原 `v40_*` 亮边与 `v39_*` 零结果审计字段保留实际来源。亮边测量继续带 `BRIGHT_RIM_REVIEW`。

single / compare、数字样品目录与自然排序、Region 从 2 开始、六项 nm 指标、对象加权统计和均值 ± 样本 SD 继续保留。缺失结果不填零，n=1 标准差为空。

## 验证与限制

新版本正向合成 slit 自上而下为 21/31/43/55 px，PixelSize=2 时约为 42/62/86/110 nm，明确覆盖“上三排短于 100 nm，最下排符合下限”。V40 的 55/43/31/21 px 递减场景在 V41 中作为拒测反例，历史源码和测试不改。

`tests/test_v41_slit_rules.py` 直接加载 V41，覆盖精确下限、PixelSize 换算、相等/倒置/局部逆序、空间行序、部分列与全图趋势、原图暗槽/亮边、缺失、尺寸参数不可放宽、独立 worker、全部零结果补检拒测及真实 compare 差值。

`tests/test_v41_workflows.py` 实际验证 V41 的 single/compare、自然排序、统计/作图/Excel、auto、真实 CLI/worker 和超时；真实 CLI 的 slit 输入与预期尺寸改为递增长度，via/trench 仍验证原有尺寸。

全仓库回归 **369 项全部通过，1348.641 秒**，其中 33 项直接加载 V41。运行环境和 ZIP 一致性结果见 `md/2026-10-16_V41_slit-length-order.md`。本批真实 SEM 原图/TXT 未提供，真实召回率、误检率和物理边界标定尚未验证；Windows/macOS 原生运行也未验证。缺失或被污染的最下排、原图边缘不足、错误 PixelSize、无法确认四排或长度不符合用户给定规则时，仍会拒测。
