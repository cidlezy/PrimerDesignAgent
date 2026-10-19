# Primer Design Agent

本地网页工具：上传基因和载体序列文件，选择双酶切或同源重组克隆方式，自动设计扩增引物，并生成 SnapGene 兼容的 GenBank 文件和引物订购表。

## 功能特性

### 序列解析

- 支持 SnapGene `.dna`、GenBank `.gb` / `.gbk`、FASTA `.fa` / `.fasta` 格式
- 自动读取文件中的 feature 注释，可作为插入片段的选择来源
- 自动检测 SnapGene 文件中已保存的引物（含序列和 Tm）

### 双酶切克隆

- 内置 24 种常用限制酶，可自由选择 Forward / Reverse 两端的酶
- 自动过滤**插入片段内部存在酶切位点**的酶（避免切到目的基因）
- 自动过滤**载体上不存在酶切位点**的酶（避免无效酶）
- 自动添加 NEB 推荐的**保护碱基**序列（如 EcoRI 前加 GCG，NotI 前加 GCGACC）
- 上下游酶选择联动：选定 Forward 酶后 Reverse 下拉框自动排除重复

### 同源重组（Gibson / HiFi）

- 指定载体上的插入位置和同源臂长度（默认 20 bp）
- 自动计算 Forward / Reverse 引物的同源臂序列

### 引物设计引擎

- 使用 Primer3 最近邻模型计算 Tm（SantaLucia 参数，50 mM Na+，无 Mg2+）
- 在给定长度范围内穷举所有退火长度，按"离目标 Tm 的距离 + 两条引物 Tm 差"打分
- 硬性约束：平均 Tm 必须在目标 ±(max_dtm × 2) 范围内，超出直接拒绝
- 3' 端优先 GC clamp
- 计算 hairpin Tm、自身 dimer ΔG、异源 dimer ΔG

### PCR 扩增建议

根据插入片段实际 GC 含量、两条引物 Tm、二级结构，自动生成建议：

- 高 GC（≥70%）：Touchdown PCR 方案 + DMSO / 甜菜碱 + 高 GC 聚合酶
- GC 偏高（65-70%）：Touchdown 或 DMSO
- 低 GC（≤35%）：较低退火温度，不需 Touchdown
- 两条引物 Tm 差 >5 °C：按较低 Tm 设退火，提示补偿方案
- 长片段（>3 kb）：按 30 s/kb 算延伸时间，推荐聚合酶
- hairpin / dimer 强：提示换备选引物对或加 DMSO

### 输出

- SnapGene 兼容 GenBank 文件（含 `primer_bind` 特征，SnapGene 打开时会把引物识别为真引物）
- Tab 分隔的引物订购表（名称、序列、长度、Tm、GC）
- 备选引物对列表（最多 5 组）

## 快速开始

### 安装依赖

```bash
pip install -r requirements.txt
```

注意：`sgffp` 需要从 GitHub 安装，`requirements.txt` 已写好 URL。如果安装失败，单独运行：

```bash
pip install --no-build-isolation 'sgffp @ git+https://github.com/merv1n34k/sgffp.git'
```

### 启动

```bash
python -m uvicorn app:app --host 127.0.0.1 --port 8300
```

打开浏览器访问 http://127.0.0.1:8300

## 使用方法

### 第一步：上传文件

点击"基因序列"和"载体序列"分别上传文件，然后点**解析文件**。工具会：

1. 显示两个文件的名称、长度、拓扑（环状/线性）
2. 扫描 24 种常用限制酶，自动过滤不可用的酶
3. 在黄色提示框中说明哪些酶被隐藏、为什么
4. 显示基因文件中已存在的引物（如果有）

### 第二步：设置参数

**克隆方式**：选"双酶切"或"同源重组"

**双酶切参数**：
- Forward / Reverse 端酶：下拉框只显示可用酶，标签里带保护碱基序列
- 选了 Forward 酶后，Reverse 下拉框自动排除它

**同源重组参数**：
- 插入位置：载体上 1-based 位置
- 同源臂长度：默认 20 bp

**插入片段**：从基因文件里选哪一段做克隆。如果基因文件只包含目的基因本身，选"整条基因序列"即可。如果文件是"载体+基因"合成的质粒，可以从 feature 列表里选特定区域。

**PCR 参数**：
| 参数 | 默认值 | 说明 |
|---|---|---|
| 目标 Tm | 60 °C | 两条引物平均 Tm 的目标值 |
| 最大 Tm 差 | 2 °C | Forward 和 Reverse Tm 允许的最大差值 |
| 最短引物 | 18 bp | 退火段最短长度 |
| 最长引物 | 36 bp | 退火段最长长度 |

### 第三步：查看结果

结果页面显示：

1. **Forward / Reverse 引物**：保护碱基橙色、退火段绿色；Tm、GC 含量、退火段长度、全长
2. **Tm 差、产物长度、Tm 标准**（本地 Primer3 估算）
3. **二级结构检查**：hairpin Tm、dimer ΔG
4. **备选引物对**：最多 5 组按打分排序
5. **PCR 扩增建议**：根据 GC、Tm、二级结构自动生成
6. **输出文件**：可下载的 GenBank 和订购表

## 技术说明

### Tm 计算标准

使用 Primer3 最近邻模型（`PRIMER_TM_FORMULA=1`, SantaLucia 1998 参数）：

| 条件 | 值 |
|---|---|
| 单价阳离子 Na+ | 50 mM |
| 二价阳离子 Mg2+ | 0 mM |
| dNTP | 0 mM |
| DNA 浓度 | 250 nM |

注意：这是**本地估算值**，与 SnapGene 内部算法的整数 Tm 可能有 ±2 °C 偏差。实际 PCR 时建议用 Primer3 的 Tm 减 3 °C 作为起始退火温度，再通过梯度 PCR 微调。

### 保护碱基表

保护碱基（也叫夹头碱基 / clamp bases）数据基于 NEB 的限制酶 DNA 末端切割效率表和常见克隆实验实践整理。代码内嵌于 `core.py` 的 `PROTECTIVE_SEQUENCES` 字典。如需扩展酶列表，直接在字典里加条目即可。

### GenBank 输出

输出格式参考了 [snapgene-bridge](https://github.com/Ryuwang3/snapgene-bridge) 项目的 `snapgene_genbank.py`（Apache-2.0）。GenBank 文件里带 `JOURNAL   Exported ... from SnapGene` 行时，SnapGene 会把 `primer_bind` 特征导入为**真引物**并计算 Tm。没有这一行的话，`primer_bind` 只会作为普通注释。

## 项目结构

```
PrimerDesignAgent/
├── app.py               # FastAPI 后端（路由、API 接口）
├── core.py              # 核心设计引擎（酶选择、引物设计、PCR 建议）
├── genbank_out.py       # SnapGene 兼容 GenBank 生成器
├── templates/
│   └── index.html       # 单页前端（原生 HTML/JS，无框架）
├── requirements.txt     # Python 依赖
└── README.md
```

## 已知限制

- Tm 是本地 Primer3 估算，不是 SnapGene 实际计算值
- 环状载体的同源臂跨原点处理尚未完善（默认把载体当线性处理）
- 无法直接生成 SnapGene 二进制 `.dna` 文件；输出的 GenBank 在 SnapGene 里打开后需另存为 `.dna` 保留引物信息
- 不支持 Golden Gate、Gateway、TOPO 等其他克隆方式
- 不检查引物 3' 端与载体其他区域的非特异性结合（没有全基因组比对）

## License

MIT
