# SEM before/after comparison

当前版本：[V32](SEM_0830_0831_compare_V32_bundle/README_V32.txt)。

- 默认固定 **via60**，标称 60 nm 与实际孔径分开，按图中候选估计搜索尺度。
- 取消依赖 40/60 nm 标称值的固定 CD 截断，放宽 slot/trench 的相对尺寸筛选。
- 保留 V31 原图边缘审核，弱边缘补检后再次审核，不强行拟合矩形或椭圆。
- 单张图片超过 5 分钟自动跳过，记录原因并继续处理。
- 保留完整四排列检查，统计最下排 slot 的 X/Y 尺寸。
- 延续 V29 的数字文件夹发现、多格式图片与自然排序前后配对。

下载 [V32 压缩包](SEM_0830_0831_compare_V32_bundle.zip)，解压后运行：

```bash
pip install numpy pandas scipy opencv-python matplotlib openpyxl
python sem_before_after_compare_0831_V32.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V32"
```

默认 `--via-pattern via60`，不再自动试 via40。确有需要时可显式选 `auto` 或 `via40`。
用新目录重跑；旧的 via40 结果不能通过改名变成 via60 测量。
默认单图时限 300 秒，可用 `--image-timeout-seconds 120` 缩短。
输出中的 `selected_pattern_key` 表示实际型号，`via_search_diameter_nm` 表示搜索孔径，
`x_width_nm/y_height_nm` 为实测 CD。历史 V31 及此前版本、上传的 V27/Yakun 包保留。

详见 [V32 算法与验证说明](SEM_0830_0831_compare_V32_bundle/ALGORITHM_V32.md)。
已增加大于 60 nm 的合成 via 和较大 slot/trench 验证，尚无本批真实 SEM 原图可复核。

开发验证：

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```
