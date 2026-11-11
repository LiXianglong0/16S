# 16S 扩增子分析工程说明

面向接手本目录的开发者 / Agent：本文描述当前**正式可用**的流程、目录约定、运行方式与已知口径差异。请以本文与 `pipeline/` 下脚本为准，不要假设存在已删除的调试脚本。

**先看这两节：** [§2 原始数据放哪里](#2-原始数据放哪里必读) · [§3 程序如何交付](#3-程序如何交付)

报告：多样性分析已合并进 `amplicon-ngs.report.html`（不再单独生成 `diversity.report.html`）。可选：用智谱免费模型生成 `paper/manuscript.md`。

---

## 1. 工程概览

本工程在 **WSL + conda 环境 `16s`** 下，对两批扩增子数据做与商业交付模板同口径的分析（去公司品牌），并额外生成多样性补充报告（不改动原报告章节结构）。

| 项目目录 | 类型 | 引物区 | 说明 |
|---|---|---|---|
| `GJproject/` | 古菌 16S | Arch519F–Arch915R | SILVA；可模板注释对齐 |
| `XJproject/` | 细菌 16S | 341F–806R | SILVA |
| `ZJproject/` | 真菌 ITS | ITS1F–ITS2（可改） | UNITE + FUNGuild/FungalTraits；见该目录 `下载清单.txt` |

参考交付结构在 `模板/`（只读对照，勿改流程依赖之外的文件）。

参考数据库：`database/silva138/silva-138-99-{seqs,tax}.qza`（SILVA 138，99%）。

---

## 2. 原始数据放哪里（必读）

测序公司或客户提供的**双端原始 FASTQ**，必须放在对应项目下的 **`rawData/`** 目录，不能放在工程根目录，也不能放进 `valid/`、`qiime2/`。

| 项目 | 原始数据目录（Windows） | 原始数据目录（WSL） |
|---|---|---|
| 古菌 | `D:\study\16S\GJproject\rawData\` | `/mnt/d/study/16S/GJproject/rawData/` |
| 细菌 | `D:\study\16S\XJproject\rawData\` | `/mnt/d/study/16S/XJproject/rawData/` |
| 真菌 ITS | `D:\study\16S\ZJproject\rawData\` | `/mnt/d/study/16S/ZJproject/rawData/` |

### 2.1 文件命名（强制）

每个样本一对文件，命名必须为：

```text
raw.{样本名}.1.fq.gz    # R1
raw.{样本名}.2.fq.gz    # R2
```

示例：

```text
GJproject/rawData/
├── raw.CK-1.1.fq.gz
├── raw.CK-1.2.fq.gz
├── raw.CK-2.1.fq.gz
├── raw.CK-2.2.fq.gz
├── raw.GRA-1.1.fq.gz
├── raw.GRA-1.2.fq.gz
└── ...
```

说明：

- `{样本名}` 可含字母、数字、`-`、`_`（如 `S200-100-1`），**不要**再含多余的 `.`
- 扩展名固定为 `.1.fq.gz` / `.2.fq.gz`（不要用 `_R1.fastq.gz` 等别名，除非先重命名）
- 数据须**已按样本拆分**；本流程不再跑 barcode 拆分（`fastq-multx`）
- 放入后改好该项目的 `config.sh`（引物、线程等），再运行该项目的 `run_16s.sh`

新项目建议：复制整个 `GJproject/` 或 `XJproject/` 为新目录 → 清空并放入新的 `rawData/` → 修改 `config.sh` → 运行。

---

## 3. 程序如何交付

本工程对外交付的是**分析程序包**（可含示例原始数据），不是某一次测序的中间结果目录。

### 3.1 生成程序交付包

```bash
bash /mnt/d/study/16S/pipeline/make_delivery.sh
```

生成位置：

```text
交付/
├── README.txt
└── 16S扩增子分析程序/          ← 拷贝给对方或推送到 GitHub
    ├── 使用说明.txt
    ├── README.md
    ├── pipeline/               # 脚本与报告资源
    ├── database/
    │   ├── silva138/           # 16S SILVA
    │   ├── unite/              # ITS UNITE
    │   ├── funguild/
    │   └── fungaltraits/
    ├── 古菌项目/               # 含 rawData 示例 FASTQ
    ├── 细菌项目/
    └── 真菌项目/               # ITS + FUNGuild/FungalTraits
```

### 3.2 程序包内含 / 不含

| 内容 | 是否打入程序包 |
|---|---|
| `pipeline/` 脚本、环境 yml、报告静态资源 | **是** |
| `database/silva138`、`unite`、`funguild`、`fungaltraits` | **是** |
| `古菌项目/`、`细菌项目/`、`真菌项目/`（config + run + rawData） | **是** |
| 三套项目的 `raw.*.fq.gz` 原始数据 | **是**（便于复现 / GitHub） |
| `qiime2/`、`work/`、已生成的 report/valid/01.ASV.tax | **否**（需自行运行生成） |
| `模板/` 商业对照结果、`zhipu.env` 真实密钥 | **否** |

对方拿到程序后：安装环境 →（已含示例数据可直接）运行 `bash 真菌项目/run_16s.sh` 等。

### 3.3 压缩示例

```bash
cd /mnt/d/study/16S/交付
tar -czvf ../16S扩增子分析程序.tar.gz 16S扩增子分析程序
```

---

## 4. 目录结构

```
16S/
├── README.md                 # 本说明
├── pipeline/                 # 正式分析脚本（唯一入口实现）
│   ├── setup_env.sh          # 创建 conda 环境 16s
│   ├── run_16s.sh            # 主流程
│   ├── fastq_stat.py         # raw / valid 统计表
│   ├── taxa_abundance.py     # 原始丰度 + 抽平归一化 + 分级表
│   ├── build_phylogeny.py    # MAFFT + FastTree
│   ├── make_report.py        # amplicon-ngs.report.html
│   ├── supplement_diversity.py  # 多样性分析并入主报告
│   ├── write_paper.py           # 智谱模型生成论文草稿
│   ├── env/                  # conda yml
│   └── report_static/        # 报告静态资源
├── database/silva138/        # SILVA 138 qza
├── GJproject/                # 古菌 16S 项目
├── XJproject/                # 细菌 16S 项目
├── ZJproject/                # 真菌 ITS 项目（需 UNITE 等下载）
└── 模板/                     # 商业模板对照（只读）
```

每个项目目录约定：

```
{GJproject|XJproject}/
├── config.sh                 # 引物、线程、参考库、域过滤等
├── run_16s.sh                # 薄封装 → ../pipeline/run_16s.sh
├── rawData/                  # 【输入】原始双端 FASTQ（见第 2 节）
├── valid/                    # 【可交付】合并后有效序列
├── work/                     # 【不交付】中间文件与软件版本
├── qiime2/                   # 【不交付】DADA2 / 注释中间产物
├── 01.ASV.tax/               # 【可交付】核心物种与丰度结果
│   ├── 01.original/          # 未抽平原始计数
│   ├── 02.normalize/         # 抽平到最小深度的整数表
│   └── 03.phylogeny_tree/    # asv_rep*.fasta, phylogeny.tre
├── amplicon-ngs.report.html  # 【可交付】主报告（含多样性）
├── web/                      # 【可交付】主报告静态资源
├── paper/                    # 可选：AI 论文草稿
├── rawdata.stat.xls
└── valid.stat.xls
```

---

## 5. 环境要求与安装

- 操作系统：Windows 上的 **WSL**（路径以 `/mnt/d/study/16S` 为例）
- 已安装 Miniconda：`$HOME/miniconda3`
- 依赖：QIIME 2、cutadapt、FLASH、vsearch、trimmomatic、biom、mafft、FastTree 等（见 `pipeline/env/16s-minimal.yml`）

```bash
bash /mnt/d/study/16S/pipeline/setup_env.sh
conda activate 16s
```

勿在未激活 `16s` 的环境中直接跑主流程。

---

## 6. 如何运行

### 6.1 全流程（质控 → DADA2 → 注释 → 丰度 → 建树 → 报告）

```bash
conda activate 16s
bash /mnt/d/study/16S/GJproject/run_16s.sh
bash /mnt/d/study/16S/XJproject/run_16s.sh
```

等价：

```bash
bash /mnt/d/study/16S/pipeline/run_16s.sh /mnt/d/study/16S/GJproject
bash /mnt/d/study/16S/pipeline/run_16s.sh /mnt/d/study/16S/XJproject
```

### 6.2 仅重做物种注释与下游（不重跑质控/DADA2）

```bash
bash /mnt/d/study/16S/pipeline/run_16s.sh /mnt/d/study/16S/GJproject taxonomy
bash /mnt/d/study/16S/pipeline/run_16s.sh /mnt/d/study/16S/XJproject taxonomy
```

`taxonomy` 模式会重新跑 VSEARCH 注释（细菌全库较慢）、重写丰度表与系统发育树。

### 6.3 单独重跑某一步（已有 qiime2/export 时）

```bash
python3 /mnt/d/study/16S/pipeline/taxa_abundance.py /mnt/d/study/16S/GJproject
python3 /mnt/d/study/16S/pipeline/build_phylogeny.py /mnt/d/study/16S/GJproject
python3 /mnt/d/study/16S/pipeline/make_report.py /mnt/d/study/16S/GJproject
python3 /mnt/d/study/16S/pipeline/supplement_diversity.py /mnt/d/study/16S/GJproject   # 并入主报告
python3 /mnt/d/study/16S/pipeline/write_paper.py /mnt/d/study/16S/GJproject              # 论文草稿
```

### 6.4 智谱论文草稿

1. 在 `pipeline/zhipu.env` 填写 `ZHIPU_API_KEY`（或 `export ZHIPU_API_KEY=...`）
2. 默认模型 `glm-4-flash`（免费；可在 env 中改 `ZHIPU_MODEL`）
3. 全流程默认会生成；若不想生成：`WRITE_PAPER=0 bash .../run_16s.sh ...`
4. 输出：`项目/paper/manuscript.md` 与 `facts.json`

密钥**不要**放进 `交付/` 程序包（打包脚本已排除 `zhipu.env`）。

### 6.5 真菌 ITS 项目（`ZJproject/`）

对照 FUNGuild + FungalTraits 双库功能注释思路：

1. 按 `ZJproject/下载清单.txt` 准备 **UNITE**、**FUNGuild**、**FungalTraits** 与 **ITS 原始 FASTQ**
2. 确认 `ZJproject/config.sh` 中引物与 `REF_SEQS`/`REF_TAX` 路径
3. 运行：`bash ZJproject/run_16s.sh`
4. 功能注释结果在 `01.ASV.tax/04.fungal_function/`

ITS 数据不能用现有 16S 的 `rawData` 顶替，需单独扩增测序。

---

## 7. 分析流程（口径）

顺序与商业模板报告一致，实现细节如下：

1. **原始数据**：已按样本拆分，跳过 `fastq-multx`
2. **质控**：trimmomatic（与模板一致：不因质控丢弃读段导致 valid 条数系统性偏少）
3. **去引物**：cutadapt（参数见各项目 `config.sh`：`CUTADAPT_E`、`CUTADAPT_O`）
4. **合并**：FLASH（`MIN_OVERLAP` / `MAX_OVERLAP` / `FLASH_X`）
5. **有效序列统计**：在 **DADA2 去嵌合体之前** 统计 `valid/`（与模板 valid 条数对齐）
6. **DADA2**：`qiime dada2 denoise-single --p-trunc-len 0`（嵌合体在此步去除）
7. **分类**：`classify-consensus-vsearch`，SILVA 138，`--p-perc-identity 0.8 --p-query-cov 0.8 --p-min-consensus 0.8`
8. **丰度表**：见下一节
9. **系统发育树**：MAFFT 比对 + FastTree
10. **报告**：主报告 + 多样性补充报告（补充不替换原报告条目）

默认不做 PICRUSt2（`RUN_PICRUST2=0`）。

---

## 8. 丰度表与注释规则（重要）

由 `pipeline/taxa_abundance.py` 生成。

### 8.1 `01.original`

- DADA2 **绝对计数**，不抽平
- 含 `asv_taxa_table.xls` 与 `domain.xls`～`species.xls`
- `tax_summary_r/*.percent.xls`：由原始计数换算的相对丰度

### 8.2 `02.normalize`

- **抽平（rarefy）到各样本最小深度**后的整数计数（随机种子固定为 `42`）
- 不是“列和为 1 的相对丰度”
- 相对丰度见同目录 `tax_summary_r/*.percent.xls`

### 8.3 `KEEP_DOMAIN`（域过滤）

在 `config.sh` 中设置：

- 古菌：`KEEP_DOMAIN="Archaea"`
- 细菌：`KEEP_DOMAIN="Bacteria"`

过滤后再编号 ASV、再抽平。

### 8.4 古菌模板注释对齐（`TEMPLATE_ASV_*`）

古菌主导 ASV 在 **SILVA 全库 VSEARCH** 下常被注成细菌（如 Bacteroidota），而商业模板注成 **Bathyarchaeia**，且丰度很高。

因此 `GJproject/config.sh` 配置了：

```bash
TEMPLATE_ASV_REP=.../模板/.../asv_rep.fasta
TEMPLATE_ASV_TAXA=.../模板/.../asv_taxa_table.xls
```

规则：**代表序列与模板完全相同的 ASV，采用模板注释**；结果写入  
`qiime2/export/taxonomy.aligned.tsv`，再按 `KEEP_DOMAIN=Archaea` 出表。

纯 SILVA 原始结果仍在 `qiime2/export/taxonomy.tsv`。

细菌项目未设置模板对齐，仅 `KEEP_DOMAIN=Bacteria`。

---

## 9. 配置项（`config.sh`）

| 变量 | 含义 |
|---|---|
| `PRIMER_F` / `PRIMER_R` / `PRIMER_F_EXTRA` | 引物序列 |
| `THREADS` | 线程数 |
| `MIN_OVERLAP` / `MAX_OVERLAP` / `FLASH_X` | FLASH |
| `CUTADAPT_E` / `CUTADAPT_O` | cutadapt |
| `REF_SEQS` / `REF_TAX` | SILVA qza 路径 |
| `KEEP_DOMAIN` | 丰度表保留的域（`Archaea` / `Bacteria` / 空=不过滤） |
| `TEMPLATE_ASV_REP` / `TEMPLATE_ASV_TAXA` | 可选；按序列覆盖注释 |
| `RUN_PICRUST2` | `1` 时尝试 PICRUSt2（需已安装） |

新项目：复制 `GJproject/` 或 `XJproject/`，清空并放入新的 `rawData/`（命名见第 2 节），改 `config.sh` 即可。

---

## 10. 主要输出说明

| 路径 | 内容 |
|---|---|
| `valid.stat.xls` | 有效序列条数/碱基数/均长 |
| `qiime2/export/otu_table.tsv` | DADA2 特征表 |
| `qiime2/export/taxonomy.tsv` | VSEARCH+SILVA 原始注释 |
| `qiime2/export/taxonomy.aligned.tsv` | 古菌：模板序列对齐后的注释（若启用） |
| `01.ASV.tax/01.original/` | 未抽平丰度 |
| `01.ASV.tax/02.normalize/` | 抽平丰度 |
| `01.ASV.tax/03.phylogeny_tree/` | 代表序列、比对、树 |
| `amplicon-ngs.report.html` | 主报告（含 1.3 多样性分析） |
| `supplement/` | 多样性图表与 tsv |
| `paper/manuscript.md` | 智谱模型生成的简短论文草稿（可选） |

各结果子目录内有 `1.` / `2.` / `3.` 开头的中文说明 txt。

---

## 11. 与模板对照时的预期差异

接手 Agent **不要**为“完全数值一致”随意改流程，除非用户明确要求。已知差异：

1. **细菌 ASV 条数偏多**（常见多约两千条）：多为单样本、极低频 ASV，约占总读段 1～2%，门水平组成仍接近模板。模板主序列约 95% 可在本流程中按序列找回。
2. **抽平深度**取决于域过滤后的最小样本深度，与模板最小深度接近但不保证相同。
3. **古菌分类**：必须以当前“模板序列对齐 + KEEP_DOMAIN”口径理解；勿单独用 SILVA 结果否定交付表中的 Bathyarchaeia。
4. **有效序列**：绝大多数样本条数与模板相同或差 ≤3，属 FLASH/质控边界差异。

---

## 12. 给 Agent 的工作约定

1. **原始数据只认各项目 `rawData/`**，命名必须是 `raw.{样本}.1.fq.gz` / `raw.{样本}.2.fq.gz`（见第 2 节）。
2. **交付客户**按第 3 节打包；不要把 `pipeline/`、`database/`、`qiime2/`、`work/` 塞进客户包。
3. **正式脚本只在 `pipeline/`**，以及各项目的 `config.sh`、`run_16s.sh`。不要重新引入已清理的调试脚本。
4. 改流程优先改 `pipeline/run_16s.sh` 与对应 Python；改引物/库/域过滤改项目 `config.sh`。
5. 多样性已合并进主报告第 1.3 节；论文草稿用 `write_paper.py`，密钥只放 `pipeline/zhipu.env`（不进交付包）。
6. 大计算耗时长；能 `taxonomy` 模式或单步 Python 解决时不要无故重跑全流程。
7. 路径按 WSL 书写（`/mnt/d/...`）；在 Windows 侧对应 `D:\study\16S\`。
8. 未要求时不要提交 git、不要改 `模板/`、不要删除 `rawData/` 与已生成的 `qiime2/` 结果。

---

## 13. 快速自检

```bash
conda activate 16s
# 原始数据是否在正确目录、命名是否合规
ls /mnt/d/study/16S/GJproject/rawData/raw.*.1.fq.gz | head
ls /mnt/d/study/16S/XJproject/rawData/raw.*.1.fq.gz | head
# 配置与入口
test -f /mnt/d/study/16S/GJproject/config.sh
test -f /mnt/d/study/16S/pipeline/run_16s.sh
# 参考库
test -f /mnt/d/study/16S/database/silva138/silva-138-99-seqs.qza
# 关键结果是否存在
ls /mnt/d/study/16S/GJproject/01.ASV.tax/01.original/asv_taxa_table.xls
ls /mnt/d/study/16S/XJproject/01.ASV.tax/03.phylogeny_tree/phylogeny.tre
```

当前版本以本 README 与 `pipeline/` 脚本为准。
