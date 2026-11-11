#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""根据有效序列统计表生成与模板相同版式的 HTML 报告，不包含测序公司信息。"""

import os
import shutil
import sys

PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_WEB = os.path.join(PIPELINE_DIR, "report_static", "web")


def load_config(path):
    """读取项目 config.sh 中的 KEY=VALUE。"""
    cfg = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip().strip('"').strip("'")
    return cfg


def read_valid_table(path):
    """读取 valid.stat.xls，返回不含表头的数据行。"""
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise SystemExit(f"统计表为空: {path}")
        for line in handle:
            line = line.rstrip("\n")
            if line.strip():
                rows.append(line)
    if not rows:
        raise SystemExit(f"统计表没有样本: {path}")
    return rows


def read_versions(path):
    """读取流程实际记录的软件版本。"""
    if not os.path.isfile(path):
        return []
    items = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if "\t" in line:
                items.append(line.split("\t", 1))
    return items


def copy_static(project_dir):
    """复制报告页面依赖的样式和脚本，不包含公司标识图片。"""
    dest = os.path.join(project_dir, "web")
    if os.path.isdir(dest):
        shutil.rmtree(dest)
    shutil.copytree(STATIC_WEB, dest)


def write_table_json(project_dir, rows):
    """写出 report.js 读取的表格数据，行之间用 \\n，列之间用 \\t。"""
    escaped = []
    for row in rows:
        cell = row.replace("\\", "\\\\").replace('"', '\\"').replace("\t", "\\t")
        escaped.append(cell)
    data = "\\n".join(escaped)
    text = (
        "let tableJsonNoSearch = []\n\n"
        "let tableJson = [{\"id\":\"data_stat\",\"info\":{\"header\":["
        "{\"title\":\"Sample\",\"searchType\":\"text\"},"
        "{\"title\":\"Sequences\",\"searchType\":\"num\"},"
        "{\"title\":\"Bases(bp)\",\"searchType\":\"num\"},"
        "{\"title\":\"AverageLength(bp)\",\"searchType\":\"num\"}"
        "],\"data\":\"" + data + "\"}}]\n"
    )
    with open(os.path.join(project_dir, "web", "table.json.js"), "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def existing_outputs(project_dir):
    """列出已经生成的结果文件，供报告正文引用。"""
    names = [
        "rawdata.stat.xls",
        "valid.stat.xls",
        os.path.join("qiime2", "export", "otu_table.tsv"),
        os.path.join("qiime2", "export", "dna-sequences.fasta"),
        os.path.join("qiime2", "export", "taxonomy.tsv"),
        os.path.join("qiime2", "denoising-stats.qza"),
        os.path.join("qiime2", "picrust2"),
    ]
    found = []
    for name in names:
        if os.path.exists(os.path.join(project_dir, name)):
            found.append(name.replace("\\", "/"))
    stats_dir = os.path.join(project_dir, "qiime2", "export-stats")
    if os.path.isdir(stats_dir):
        for name in sorted(os.listdir(stats_dir)):
            found.append("qiime2/export-stats/" + name)
    return found


def build_html(project_dir, cfg, rows, versions):
    """组装报告 HTML。实验步骤文字保留，页眉和公司联系方式去掉。"""
    organism = cfg.get("ORGANISM", "")
    region = cfg.get("REGION", "")
    primer_f_name = cfg.get("PRIMER_F_NAME", "")
    primer_r_name = cfg.get("PRIMER_R_NAME", "")
    primer_f = cfg.get("PRIMER_F", "")
    primer_f_extra = cfg.get("PRIMER_F_EXTRA", "")
    primer_r = cfg.get("PRIMER_R", "")
    threads = cfg.get("THREADS", "8")
    extra = f"；兼容正向引物 {primer_f_extra}" if primer_f_extra else ""
    project_info = (
        f"本项目分析对象为{organism}，扩增区域为 {region}，"
        f"正向引物 {primer_f_name}（{primer_f}{extra}），"
        f"反向引物 {primer_r_name}（{primer_r}）。"
        f"原始数据已按样本拆分，样本数 {len(rows)}。"
        "每个样本包含双端文件 raw.{样本名}.1.fq.gz 与 raw.{样本名}.2.fq.gz。"
    )
    root = project_dir.replace("\\", "/")
    commands = f"""#1) 导入数据-单端：
qiime tools import --input-path {root}/qiime2/manifest.tsv --output-path {root}/qiime2/data.qza --type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2

# manifest 格式
sample-id	absolute-filepath
样本名	{root}/valid/trim.样本名.fq.gz

#2) DADA2 滤除有噪声的序列，校正测序错误，去除嵌合序列和单体，再去冗余
qiime dada2 denoise-single --i-demultiplexed-seqs {root}/qiime2/data.qza --p-n-threads {threads} --p-trim-left 0 --p-trunc-len 0 --o-table {root}/qiime2/table.qza --o-representative-sequences {root}/qiime2/rep-seqs.qza --o-denoising-stats {root}/qiime2/denoising-stats.qza --o-base-transition-stats {root}/qiime2/base-transition-stats.qza

#3) 导出 biom
qiime tools export --input-path {root}/qiime2/table.qza --output-path {root}/qiime2/export

#4) biom 转为表格
biom convert -i {root}/qiime2/export/feature-table.biom -o {root}/qiime2/export/otu_table.tsv --to-tsv"""

    outputs = existing_outputs(project_dir)
    output_html = "".join(f"<p>{item}</p>" for item in outputs) if outputs else "<p>结果文件将在流程结束后出现在项目目录中。</p>"

    version_rows = "".join(f"<tr><td>{name}</td><td>{ver}</td></tr>" for name, ver in versions)
    if not version_rows:
        version_rows = "<tr><td colspan=\"2\">运行流程后会记录实际版本</td></tr>"

    return f"""<!DOCTYPE html>
<html lang="zh-CN">

<head>
    <link rel="icon" href="data:,">
    <meta charset="UTF-8">
    <title>微生物群落多样性测序</title>
    <link type="text/css" rel="stylesheet" href="web/js/bootstrap/css/bootstrap.css">
    <link type="text/css" rel="stylesheet" href="web/js/bootstrap-table-master/dist/bootstrap-table.css">
    <link type="text/css" rel="stylesheet" href="web/css/report.css">
    <link type="text/css" rel="stylesheet" href="web/js/bootstrap-advanced-sortable/css/bootstrap-advanced-sortable.css">
    <link rel="stylesheet" href="web/css/font-awesome-4.7.0/css/font-awesome.min.css"/>
    <script src="web/js/jquery-3.1.0.min.js" type="text/javascript"></script>
    <script type="text/javascript" src="web/js/bootstrap/js/bootstrap.js"></script>
    <script type="text/javascript" src="web/js/bootstrap-table-master/dist/bootstrap-table.js"></script>
    <script type="text/javascript" src="web/js/bootstrap-advanced-sortable/js/bootstrap-advanced-sortable.js"></script>
</head>

<body>

<div class="leftsidebar_box" id="myScrollspy">
    <div class="left-logo-box">
        <div style="color:#1297AD;font-weight:bold;font-size:22px;padding:12px 0;">16S</div>
    </div>
    <div class="bs-docs-sidebar">
        <ul class="nav bs-docs-sidenav">
            <li><a href="#chapter0" style="color: white;">背景介绍 <i class="fa fa-caret-down"></i></a>
                <ul class="nav-stacked">
                    <li><a href="#chapter0_1">项目概况</a></li>
                </ul>
            </li>
            <li><a href="#chapter1" style="color: white;">工作流程 <i class="fa fa-caret-down"></i></a>
                <ul class="nav-stacked">
                    <li><a href="#chapter1_1">1.1 NGS测序实验流程</a></li>
                    <li><a href="#chapter1_1_1">&nbsp;&nbsp;1.1.1 基因组DNA抽提</a></li>
                    <li><a href="#chapter1_1_2">&nbsp;&nbsp;1.1.2 PCR扩增</a></li>
                    <li><a href="#chapter1_1_3">&nbsp;&nbsp;1.1.3 荧光定量</a></li>
                    <li><a href="#chapter1_1_4">&nbsp;&nbsp;1.1.4 二代测序文库构建</a></li>
                    <li><a href="#chapter1_1_5">&nbsp;&nbsp;1.1.5 NGS测序</a></li>
                    <li><a href="#chapter1_2">1.2 生物信息分析流程</a></li>
                    <li><a href="#chapter1_2_1">&nbsp;&nbsp;1.2.1 数据拆分</a></li>
                    <li><a href="#chapter1_2_2">&nbsp;&nbsp;1.2.2 数据统计</a></li>
                    <li><a href="#chapter1_2_3">&nbsp;&nbsp;1.2.3 qiime2导入valid数据</a></li>
                </ul>
            </li>
            <li><a href="#chapter2">附录 <i class="fa fa-caret-down"></i></a>
                <ul class="nav-stacked">
                    <li><a href="#chapter2_1">2.1 文件解压缩方法</a></li>
                    <li><a href="#chapter2_2">2.2 文件打开或浏览方法</a></li>
                    <li><a href="#chapter2_3">2.3 分析软件</a></li>
                </ul>
            </li>
        </ul>
    </div>
</div>

<div id="content">
    <div id="head">
        <div id="title">
            <h2>微生物群落多样性测序</h2>
        </div>
    </div>
    <div id="floatbar">
        <a id="top" href="#main"><img src="web/image/Top.jpg" /></a>
    </div>
    <div id="main">

        <h2 id="chapter0" class="main-title">背景介绍</h2>
        <h3 id="chapter0_1" class="sub-title">项目概况</h3>
        <p>微生物是一群以分解代谢为主的生物类群，其生物活性十分丰富。根据微生物的生物特性，从物种、生理代谢类群及遗传背景几个方面探讨微生物多样性的问题，使人们更加的了解微生物的多样性，更加全面的认识微生物的多样性，为解决人类目前面临的诸多问题提供更多的途径。基因组时代的到来，为研究微生物群落结构提供了新的技术平台，将一个全新的微生物世界展现在人们面前，它充分显示了微生物的生物多样性，并揭示了生物进化的动力。其中，16S rRNA基因测序已经成为当前研究原核微生物群落组成及其分布的重要手段。</p>
        <p>{project_info}</p>

        <h2 id="chapter1" class="main-title">工作流程</h2>
        <h3 id="chapter1_1" class="sub-title">1.1&nbsp;NGS测序实验流程</h3>
        <br>
        <div align="center">
            <img src="web/image/illumina_lab_flow.jpg" width="100%" height="100%">
        </div>
        <br>
        <h4 id="chapter1_1_1" class="third-title">1.1.1&nbsp;基因组DNA抽提</h4>
        <p>完成基因组DNA抽提后，利用1%琼脂糖凝胶电泳检测抽提的基因组DNA。</p>
        <br>
        <h4 id="chapter1_1_2" class="third-title">1.1.2&nbsp;PCR扩增</h4>
        <p>按指定测序区域，合成带有barcode的特异引物。</p>
        <p>为保证后续数据分析的准确性及可靠性，需满足两个条件，1) 尽可能使用低循环数扩增；2) 保证每个样本扩增的循环数一致。随机选取具有代表性的样本进行预实验，确保在最低循环数中使绝大多数样本能够扩增出浓度合适的产物。</p>
        <p>PCR 采用TransGen AP221-02：TransStart Fastpfu DNA Polymerase；</p>
        <p>PCR仪：ABI GeneAmp® 9700型；</p>
        <p>全部样本按照正式实验条件进行，每个样本3个重复，将同一样本的PCR产物混合后用2%琼脂糖凝胶电泳检测，使用AxyPrepDNA凝胶回收试剂盒（AXYGEN公司）切胶回收PCR产物，Tris_HCl洗脱； 2%琼脂糖电泳检测。</p>
        <p>PCR 采用TransStart Fastpfu DNA Polymerase，20μl反应体系：</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5×FastPfu Buffer .....................4 μl</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.5 mM dNTPs.........................2 μl</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Forward Primer(5 μM).......... 0.8μl</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Reverse Primer(5 μM) ......... 0.8μl</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FastPfu Polymerase..............0.4 μl</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Template DNA ......................10ng</p>
        <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;补ddH2O至.........................20 μl</p>
        <br>
        <p>PCR反应参数</p>
        <p>&nbsp;a. 1cycles× (5 minutes at 95°C)；</p>
        <p>&nbsp;b. cycles× (30 seconds at 95°C；30 seconds at 55°C；45 seconds at 72°C)；</p>
        <p>&nbsp;c. 10 minutes at 72°C，10°C until halted by user。</p>
        <br>
        <h4 id="chapter1_1_3" class="third-title">1.1.3&nbsp;荧光定量</h4>
        <p>参照电泳初步定量结果，将PCR产物用QuantiFluor™ -ST蓝色荧光定量系统（Promega公司）进行检测定量，之后按照每个样本的测序量要求，进行相应比例的混合。</p>
        <br>
        <h4 id="chapter1_1_4" class="third-title">1.1.4&nbsp;NGS文库构建</h4>
        <p>&nbsp;1) 连接“Y”字形接头；</p>
        <p>&nbsp;2) 使用磁珠筛选去除接头自连片段；</p>
        <p>&nbsp;3) 利用PCR扩增进行文库模板的富集；</p>
        <p>&nbsp;4) 氢氧化钠变性，产生单链DNA片段。</p>
        <br>
        <h4 id="chapter1_1_5" class="third-title">1.1.5&nbsp;NGS测序</h4>
        <p>&nbsp;1) DNA片段的一端与引物碱基互补，固定在芯片上；</p>
        <p>&nbsp;2) 另一端随机与附近的另外一个引物互补，也被固定住，形成“桥 (bridge)”；</p>
        <p>&nbsp;3) PCR扩增，产生DNA簇；</p>
        <p>&nbsp;4) DNA扩增子线性化成为单链；</p>
        <p>&nbsp;5) 加入改造过的DNA聚合酶和带有4种荧光标记的dNTP，每次循环只合成一个碱基；</p>
        <p>&nbsp;6) 用激光扫描反应板表面，读取每条模板序列第一轮反应所聚合上去的核苷酸种类；</p>
        <p>&nbsp;7) 将“荧光基团”和“终止基团”化学切割，恢复3'端粘性，继续聚合第二个核苷酸；</p>
        <p>&nbsp;8) 统计每轮收集到的荧光信号结果，获知模板DNA片段的序列；</p>
        <br>
        <h3 id="chapter1_2" class="sub-title">1.2 生物信息分析流程</h3>
        <h4 id="chapter1_2_1" class="sub-title">1.2.1 数据拆分</h4>
        <p>第一步：拆 barcode。原流程软件为 fastq-multx。本项目收到的原始数据已经按样本拆成双端文件，因此不再重复拆分。</p>
        <p>第二步：质控。使用软件：trimmomatic。参数：LEADING:3 TRAILING:3 MINLEN:50。不使用滑窗截短，避免把仍可合并的读段丢掉。</p>
        <p>第三步：匹配引物。使用软件：cutadapt。正向引物 {primer_f_name}，反向引物 {primer_r_name}，错配比例 {cfg.get("CUTADAPT_E", "0.45")}，最短匹配 {cfg.get("CUTADAPT_O", "8")} bp，未匹配到引物的序列丢弃。</p>
        <p>第四步：合并 R1、R2。使用软件：flash。最小重叠 {cfg.get("MIN_OVERLAP", "8")} bp，最大重叠 {cfg.get("MAX_OVERLAP", "300")} bp，最大错配密度 {cfg.get("FLASH_X", "0.6")}。</p>
        <p>第五步：去除嵌合体。原流程软件为 identify_chimeric_seqs.py。模板给出的 valid 序列条数与原始读段对数相同，这一步没有从有效序列中扣掉嵌合体。本流程与该口径一致，嵌合体在 DADA2 denoise-single 中去除。</p>
        <br>
        <h4 id="chapter1_2_2" class="sub-title">1.2.2 数据统计</h4>
        <p>各样本序列统计见下表：</p>
        <b class="tuming">有效序列统计</b>
        <table id="data_stat" style="width: 1000px;"></table>
        <p>Sample：样本名称；</p>
        <p>Sequences：总序列条数；</p>
        <p>Bases(bp)：样本优化序列碱基数；</p>
        <p>AverageLength(bp)：优化数据的序列平均长度。</p>

        <h4 id="chapter1_2_3" class="sub-title">1.2.3 qiime2导入valid数据</h4>
        <p>DADA2 滤除有噪声的序列，校正不确定序列中的错误，去除嵌合序列，去除单体（singletons，出现频率仅有一次的序列），然后对这些序列进行去冗余，得到 ASV 表。</p>
        <pre style="white-space:pre-wrap;background:#f7f7f7;padding:12px;">{commands}</pre>
        <p>本项目结果文件：</p>
        {output_html}

        <h2 id="chapter2" class="main-title">附录</h2>
        <h3 id="chapter2_1" class="sub-title">2.1&nbsp;文件解压缩方法</h3>
        <p>结果中的压缩包为 Linux 下的文件，压缩包使用 tar -zcvf 命令压缩。不同系统的解压方法如下：</p>
        <p>Unix/Linux/Mac用户：使用 tar -zxvf *.tar.gz 命令</p>
        <p>Windows用户：使用 WinRAR 软件解压缩</p>
        <br>
        <h3 id="chapter2_2" class="sub-title">2.2&nbsp;文件打开或浏览方法</h3>
        <p>如果在本附录中无特殊说明，结果文件为 Linux 系统下文本文件。Unix/Linux 用户可以使用 more 或 less 命令查看文本文件内容。对于 Windows 用户，一般文本文件可以使用写字板或 Excel 打开。推荐使用文本编辑器。当文件比较大时，打开文件可能导致 Windows 系统卡顿，建议使用性能较好的计算机或者使用更适合处理大量数据的 Unix/Linux 系统打开。</p>
        <p>数据中可能包含部分图像文件，一般图像文件后缀名为 .png、.jpg、.gif 等。对于图像文件，Windows 用户可以使用图片浏览器打开，Linux/Unix 用户使用 display 命令打开。</p>
        <p>后缀名为 svg 的文件为文本格式描述的图像文件，Windows 用户需要安装 Adobe Illustrator 软件打开。Linux/Unix 用户可以使用 rsvg-view 命令查看。结果中如提供 pdf 格式的矢量图，可利用 Adobe Illustrator 对该格式图片进行编辑。</p>
        <p>Linux 下的表格均为制表符（Tab）分割的文本。为了便于阅读，建议使用 Excel 或 OpenOffice 等办公软件用表格形式打开，打开时请用 Tab 分割方式。</p>
        <br>
        <h3 id="chapter2_3" class="sub-title">2.3&nbsp;分析软件</h3>
        <p>原分析记录中的软件和参数：</p>
        <table class="table table-bordered" style="width: 800px;">
            <tr><th>分析内容</th><th>软件</th><th>版本</th><th>参数</th></tr>
            <tr><td>生成ASV表</td><td>DADA2</td><td>v1.18.0</td><td>default</td></tr>
            <tr><td>分类学分析</td><td>uclust algorithm</td><td>v1.2.22q</td><td>置信度阈值为0.8</td></tr>
            <tr><td>PICRUSt2功能预测分析</td><td>qiime2</td><td>2020.11</td><td>default</td></tr>
        </table>
        <p>本流程实际使用的软件版本：</p>
        <table class="table table-bordered" style="width: 800px;">
            <tr><th>软件</th><th>版本</th></tr>
            {version_rows}
        </table>
        <p>分类学分析在配置参考序列和注释表后执行。开源环境使用 QIIME 2 的 classify-consensus-vsearch，序列相似度阈值 0.8，共识阈值 0.8，对应原记录中的 uclust 置信度 0.8。PICRUSt2 在配置文件中将 RUN_PICRUST2 设为 1 且本机已安装该软件时执行。</p>
        <br>
    </div>
</div>

<script type="text/javascript" src="web/table.json.js"></script>
<script type="text/javascript" src="web/js/report.js"></script>
</body>
</html>
"""


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 make_report.py 项目目录")
    project_dir = os.path.abspath(sys.argv[1])
    valid_stat = os.path.join(project_dir, "valid.stat.xls")
    if not os.path.isfile(valid_stat):
        raise SystemExit(f"还没有有效序列统计表: {valid_stat}")
    cfg = load_config(os.path.join(project_dir, "config.sh"))
    rows = read_valid_table(valid_stat)
    versions = read_versions(os.path.join(project_dir, "work", "software_versions.txt"))
    copy_static(project_dir)
    write_table_json(project_dir, rows)
    html = build_html(project_dir, cfg, rows, versions)
    out = os.path.join(project_dir, "amplicon-ngs.report.html")
    with open(out, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(html)
    print(out)


if __name__ == "__main__":
    main()
