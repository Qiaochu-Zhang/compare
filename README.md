# SEM before/after comparison

当前版本：[V31](SEM_0830_0831_compare_V31_bundle/README_V31.txt)。

- 结合 V27 的定位/输出、Yakun 的原图灰度边缘、V30 的缺点恢复。
- 保留有原图支持的形状；弱边缘补检后再次审核，不强行拟合矩形或椭圆。
- 单张图片超过 5 分钟自动跳过，记录原因并继续处理。
- 保留完整四排列检查，统计最下排 slot 的 X/Y 尺寸。
- 延续 V29 的数字文件夹发现、多格式图片与自然排序前后配对。

下载 [V31 压缩包](SEM_0830_0831_compare_V31_bundle.zip)，解压后运行：

```bash
pip install numpy pandas scipy opencv-python matplotlib openpyxl
python sem_before_after_compare_0831_V31.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V31"
```

默认单图时限 300 秒，可用 `--image-timeout-seconds 120` 缩短。
输出表格及标注参考 V27，新增边缘支持率和候选尺寸，便于核查误测。
已有版本及用户上传的 V27、Yakun 源码包保留。

详见 [算法与验证说明](SEM_0830_0831_compare_V31_bundle/ALGORITHM_V31.md)。目前验证基于合成图和流程测试，尚未复核用户报告误测的原始 SEM 图片。

开发验证：

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```
