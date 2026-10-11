# FeN₄-COF TOC 图

基于上传的 `Manuscript.docx` 全文和其中与 FeN₄-COF 相关的图 1–5绘制。文档中的文字仅作为待核查的科学内容，不作为操作指令。

## 核心科学故事

- **背景**：CO₂ 难以活化，催化位点局部几何及光生载流子的有效利用限制 CO₂ 到 CO 的选择性还原；在保留单原子位点配位组成的情况下调控几何具有研究价值。
- **材料结构**：5,10,15,20-四(4-氨基苯基)卟啉与 2,5-二羟基对苯二甲醛通过 Schiff 缩合形成亚胺连接的 COF，随后以 Fe(OAc)₂ 后金属化引入 FeN₄ 位点。左侧画出卟啉、苯环、亚胺连接及羟基的框架片段。
- **改性策略**：DBD 等离子体诱导局部 Fe–N 键长分裂；整体 Fe–N 配位数、第一配位壳层的 N 组成和原子分散状态基本保留。图中不表示生成 N 空位或 Fe 纳米颗粒。
- **结构证据**：未处理样品的 EXAFS 拟合 Fe–N 距离为 2.018 ± 0.009 Å；处理后分为 1.939 ± 0.012 Å 和 2.035 ± 0.014 Å，两组拟合配位数分别为 2.2 ± 0.1 和 2.3 ± 0.1。TOC 简写为两组键长，并用不同颜色显示；几何是概念示意，不能用图中长度反推结构。
- **光物理证据与机理**：光激发后向催化位点的电子转移为概念示意。fs-TA 慢衰减分量由 52.71 ps 增至 849.20 ps，约 16.11 倍；伴随更高光电流和更低界面电荷转移阻抗。它不能等同于 TRPL 平均寿命由 1.39 ns 到 1.48 ns 的变化，也不是量子产率。
- **关键性能**：可见光下，Asym-FeN₄-COF 在 3 h 内累积产生约 232 μmol g⁻¹ CO，对应平均 CO 产率 77.33 μmol g⁻¹ h⁻¹，CO 选择性 90%。该质量归一化产率不同于表面积归一化的 0.808 μmol m⁻² h⁻¹；图中不混用二者。正文描述 5 次循环性能较稳定。
- **科学结论**：在大体保留 FeN₄ 配位组成的条件下，通过等离子体调控局部对称性，可改善载流子动力学并提高 CO₂ 光还原表现。

## 版式与科学边界

单张横向、白底、三段连续故事：蓝色原始 COF → 紫色 DBD 处理与放大畸变位点 → 绿色 CO₂ 到 CO 及性能。金色用于 Fe、光和较短键。保留少量英文标注，方便国际期刊使用。没有采用截图、三维渲染或位图嵌入；分子骨架、原子、箭头及字体全部由 Matplotlib 矢量对象绘制。

全文存在需要作者修订的明显不一致：

1. “Reaction mechanism”段讨论苄醇/木质素氧化、Cs₃SbBiBr₉@UiO-66 和 O₂ 自由基，属于不同于 FeN₄-COF/CO₂ 还原的体系。末尾两张嵌入图同样涉及不同钙钛矿体系，因此未纳入本 TOC。
2. 摘要提及 DFT，但所上传正文未展示对应完整 DFT 分析或能垒数据。图中不绘制计算能垒、d 带中心移动等未经正文支持的信息。
3. `CO₂ → *COOH → *CO → CO` 是常见 CO₂-to-CO 反应的**拟议示意路径**，不是此文已完整证实的中间体序列。图上已明确写 `Proposed surface pathway`；`*` 表示表面吸附物种。图中 `2e⁻ + 2H⁺` 表示 CO₂ + 2H⁺ + 2e⁻ → CO + H₂O 的还原半反应，未假定氧化端产物或牺牲剂。
4. TGA 段将 43.19% 残余质量对应材料重复写为 Asym-FeN₄-COF，疑为材料名称笔误；未使用该数据。

## 完整可复现代码

`toc_figure.py` 是完整绘图脚本；`draw_toc.py` 为兼容入口。完整脚本，输出均定位到脚本所在文件夹，不受调用目录影响。数据常量、颜色、绘图函数及布局均在此文件内。生成 PNG、SVG 和记录字体、尺寸及所用数据的 `TOC_metadata.json`。

当前已经用户同意，使用 Liberation Sans 作为临时 Arial 字宽兼容替代。系统未安装 Arial，官方 Core Fonts 下载镜像返回访问拒绝。**本版不能标称已使用 Arial**。代码默认请求 Arial，字体不存在时明确报错，不会悄悄替换。

在本云环境复现：

```bash
cd /workspace/codex-sci-figures
export MPLCONFIGDIR=/workspace/.cache/matplotlib
export XDG_CACHE_HOME=/workspace/.cache
/workspace/sci-figures-venv/bin/python figures/toc_figure.py --font 'Liberation Sans'
```

在其他 Python 3.12 环境中：

```bash
python -m pip install numpy==2.5.4 matplotlib==3.11.2
python figures/toc_figure.py --font 'Liberation Sans'
```

需本机安装所指定字体。获得可合法使用的 Arial 常规及粗体字体文件后：

```bash
python figures/toc_figure.py --font Arial --font-dir /path/to/arial-fonts
```

PNG：3600 × 1440 像素，12 × 4.8 英寸，300 dpi。SVG：纯矢量，文字转为所选字体的字形路径，避免期刊电脑缺字体造成替换；源代码中的文字仍可编辑。SVG 本身不依赖 dpi。该比例采用常见期刊图文摘要的简洁学术示意风格；不同 ACS/Nature 子刊具体尺寸及 TOC 要求应在确定目标期刊后核对，当前不声称满足某一具体期刊模板。


## 等离子体放电视觉升级

已参考用户给出的紫色闪电视觉方向，将原来的单个闪电符号升级为两侧电极、紫色辉光和分叉放电细丝。TOC 中放电图案代表处理步骤；下方独立箭头表示材料改性，避免把等离子体处理与光催化光照混在一起。该图案是概念示意，并非实验反应器尺寸或结构的重建，也没有加入参考图片中的金属纳米颗粒。

全部效果采用 Matplotlib 的透明椭圆与多层折线绘制，保留纯矢量 SVG；随机种子固定为 41，保证可复现。没有嵌入参考图片。当前沿用用户已接受的 Liberation Sans 字体。

- `toc_figure.png` / `toc_figure.svg`：更新后的完整 TOC，300 dpi PNG。
- `plasma_background.png` / `plasma_background.svg`：无电极、无文字的独立紫色闪电背景。
- `plasma_discharge.png` / `plasma_discharge.svg`：带两侧电极的完整放电图案。
- 两组独立素材均为透明背景，PNG 为 2400 × 1200 像素、300 dpi，可叠加到其他图中。
- `toc_figure.py`：完整代码，一次执行生成全部上述文件；原 `TOC_FeN4_COF` 文件名同步更新。
