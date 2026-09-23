# SEM before/after comparison

当前版本：[V30](SEM_0830_0831_compare_V30_bundle/README_V30.txt)。

- 改进 slot 在亮度不均、低对比图中的定位与完整轮廓提取。
- 单张图片超过 5 分钟自动跳过，记录原因并继续处理。
- 保留完整四排列检查，统计最下排 slot 的 X/Y 尺寸。
- 延续 V29 的数字文件夹发现、多格式图片与自然排序前后配对。

下载 [V30 压缩包](SEM_0830_0831_compare_V30_bundle.zip)，解压后运行：

```bash
pip install numpy pandas scipy opencv-python matplotlib openpyxl
python sem_before_after_compare_0831_V30.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V30"
```

默认单图时限 300 秒，可用 `--image-timeout-seconds 120` 缩短。
已有 V29 文件保留，便于复查历史结果。

开发验证：

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```
