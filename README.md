# SEM before/after comparison

当前版本：[V33](SEM_0830_0831_compare_V33_bundle/README_V33.txt)。

- **via 圆形度阈值**：默认 0.75，只对 via 生效；低于阈值的对象从有效标注和统计中剔除。
- **via 同图尺寸偏差**：默认允许等效直径相对同图中位数 ±30%，不按 before/after 差异剔除。
- **尺寸辅助识别**：可预填 via 直径、slit 宽度/最下排长度、trench 宽度/端点间隙，均为 nm。
- 提示尺寸参与定位和模糊区域补检，最终仍由原图边缘测量；默认尺寸容差 ±50%。
- 默认固定 **via60**，保留原图边缘审核，弱边缘补检后再次审核。
- 单张图片超过 5 分钟自动跳过，记录原因并继续处理。
- 保留完整四排列检查，统计最下排 slot 的 X/Y 尺寸。
- 延续 V29 的数字文件夹发现、多格式图片与自然排序前后配对。

下载 [V33 压缩包](SEM_0830_0831_compare_V33_bundle.zip)，解压后运行：

```bash
pip install numpy pandas scipy opencv-python matplotlib openpyxl
python sem_before_after_compare_0831_V33.py --before "D:/SEM/before" --after "D:/SEM/after" --output "D:/SEM/compare_V33"
```

可调参数（只填写已知尺寸，其余省略）：

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| `--via-min-circularity` | via 圆形度下限，越大越严格 | `0.75` |
| `--via-max-size-deviation` | via 等效直径相对同图中位数的最大偏差，0 关闭 | `0.30` |
| `--via-diameter-nm` | 预计 via 直径 | 不指定 |
| `--slit-width-nm` | 预计 slit 宽度，也可写 `--slot-width-nm` | 不指定 |
| `--slit-length-nm` | 预计最下排 slit 长度，也可写 `--slot-length-nm` | 不指定 |
| `--trench-width-nm` | 预计 trench 宽度 | 不指定 |
| `--trench-gap-nm` | 预计上下 trench 端点间隙 | 不指定 |
| `--size-tolerance` | 所填尺寸的相对容差 | `0.50`（±50%） |

所有参数保存到 `settings.json`，剔除原因记录在 `rejected_objects.csv`。
同图尺寸筛选至少需要 3 个有效 via；尺寸提示统一用于本次运行的所有条件和前后阶段。
slit 的四排长度可不同，长度范围只检查最下排；测量仍要求完整四排。

默认 `--via-pattern via60`，确有需要时可显式选 `auto` 或 `via40`。
默认单图时限 300 秒，可用 `--image-timeout-seconds 120` 缩短。用新输出目录重跑。
历史 V32 及此前版本、上传的 V27/Yakun 包保留。

详见 [V33 使用说明](SEM_0830_0831_compare_V33_bundle/README_V33.txt)、
[算法与验证说明](SEM_0830_0831_compare_V33_bundle/ALGORITHM_V33.md)及
[本次修改记录](md/2026-10-13_V33_via-quality-size-hints.md)。
已增加合成模糊 slit、via 形状/尺寸异常与尺寸提示验证，尚无本批真实 SEM 原图可复核。

开发验证：

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

Codex 工作规则见 [AGENTS.md](AGENTS.md)：开始工作时按版本号读取最新源码，
修改后记录到 [md/](md/README.md)，新代码版本完成后生成对应 ZIP，并提交、推送到 GitHub。
