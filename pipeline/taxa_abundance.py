#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 ASV 表和 SILVA 注释写成两份丰度表。

01.original：DADA2 原始计数（可按 KEEP_DOMAIN 只保留目标域）。
02.normalize：抽平到各样本最小深度后的整数计数（与模板口径一致）。
相对丰度百分比另外写在 tax_summary_r/*.percent.xls。

若 config 中设置了 TEMPLATE_ASV_REP + TEMPLATE_ASV_TAXA，
则对序列完全相同的 ASV 采用模板注释（用于古菌引物产物与商业模板口径对齐）。
"""

import os
import random
import sys
from collections import OrderedDict

RANKS = [
    ("domain", "d__"),
    ("phylum", "p__"),
    ("class", "c__"),
    ("order", "o__"),
    ("family", "f__"),
    ("genus", "g__"),
    ("species", "s__"),
]

RAREFY_SEED = 42


def load_config(project):
    cfg = {}
    path = os.path.join(project, "config.sh")
    if not os.path.isfile(path):
        return cfg
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip().strip('"').strip("'")
    return cfg


def load_fasta(path):
    seqs = {}
    name = None
    parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).replace("-", "").upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts).replace("-", "").upper()
    return seqs


def read_otu(path):
    samples = []
    counts = OrderedDict()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("# "):
                continue
            parts = line.split("\t")
            if parts[0] in ("#OTU ID", "OTU ID"):
                samples = parts[1:]
                continue
            counts[parts[0]] = [int(round(float(x))) for x in parts[1:]]
    if not samples or not counts:
        raise SystemExit(f"特征表为空: {path}")
    return samples, counts


def read_taxonomy(path):
    taxa = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] in ("Feature ID", "FeatureID"):
                continue
            taxon = parts[1].strip() if len(parts) > 1 else "Unassigned"
            taxa[parts[0]] = taxon.replace("; ", ";")
    return taxa


def overlay_template_taxonomy(taxa, our_fasta, tmpl_fasta, tmpl_taxa_xls):
    """序列完全相同的 ASV，用模板注释覆盖。返回覆盖条数。"""
    ours = load_fasta(our_fasta)
    tmpls = load_fasta(tmpl_fasta)
    by_seq = {seq: fid for fid, seq in ours.items()}
    tmpl_tax = {}
    with open(tmpl_taxa_xls, "r", encoding="utf-8") as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            tmpl_tax[parts[0]] = parts[-1].replace("; ", ";")
    n = 0
    for asv, seq in tmpls.items():
        if asv not in tmpl_tax:
            continue
        feature = by_seq.get(seq)
        if feature is None:
            continue
        taxa[feature] = tmpl_tax[asv]
        n += 1
    return n


def write_taxonomy(path, taxa, feature_order):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("Feature ID\tTaxon\tConfidence\n")
        for feature in feature_order:
            handle.write(f"{feature}\t{taxa.get(feature, 'Unassigned')}\t\n")
        for feature, taxon in taxa.items():
            if feature not in feature_order:
                handle.write(f"{feature}\t{taxon}\t\n")


def keep_domain(taxon, domain):
    """按域/界过滤。细菌古菌用 d__；真菌 UNITE 常用 k__Fungi。"""
    if not domain:
        return True
    compact = (taxon or "").replace("; ", ";")
    candidates = [f"d__{domain}", f"k__{domain}"]
    for prefix in candidates:
        if compact.startswith(prefix + ";") or compact == prefix:
            return True
    return False


def ranks_for_taxa(taxa_values):
    """16S 常用 d__ 开头；ITS/UNITE 常用 k__ 开头。"""
    for taxon in taxa_values:
        t = (taxon or "").replace("; ", ";")
        if t.startswith("k__"):
            return [
                ("kingdom", "k__"),
                ("phylum", "p__"),
                ("class", "c__"),
                ("order", "o__"),
                ("family", "f__"),
                ("genus", "g__"),
                ("species", "s__"),
            ]
    return RANKS


def rank_label(taxon, level, ranks=None):
    ranks = ranks or RANKS
    if not taxon or taxon == "Unassigned":
        return "Unassigned"
    parts = [p for p in taxon.split(";") if p]
    if level >= len(parts):
        return "Unclassified"
    name = parts[level]
    prefix = ranks[level][1]
    if name in (prefix, prefix.rstrip("_"), "uncultured", "norank", "unidentified"):
        return "Unclassified"
    return name


def sample_depths(counts, n):
    depth = [0] * n
    for values in counts.values():
        for i, value in enumerate(values):
            depth[i] += value
    return depth


def rarefy_table(counts, n, depth, seed=RAREFY_SEED):
    rng = random.Random(seed)
    features = list(counts.keys())
    out = OrderedDict((f, [0] * n) for f in features)
    for i in range(n):
        pool = []
        for feature in features:
            c = counts[feature][i]
            if c:
                pool.extend([feature] * c)
        if len(pool) < depth:
            raise SystemExit(f"样本列 {i} 深度 {len(pool)} < 抽平深度 {depth}")
        for feature in rng.sample(pool, depth):
            out[feature][i] += 1
    return OrderedDict((f, v) for f, v in out.items() if sum(v) > 0)


def write_table(path, samples, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(["ASV ID"] + samples + ["taxonomy"]) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def write_rank(path, samples, label_counts, as_float):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ordered = sorted(label_counts.items(), key=lambda item: -sum(item[1]))
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(["ASV ID"] + samples) + "\n")
        for label, values in ordered:
            if as_float:
                cells = [f"{v:.6f}" for v in values]
            else:
                cells = [str(int(v)) for v in values]
            handle.write("\t".join([label] + cells) + "\n")


def aggregate_ranks(items, n, ranks=None):
    ranks = ranks or RANKS
    maps = [OrderedDict() for _ in ranks]
    for _total, _feature, values, taxon in items:
        for level, abs_map in enumerate(maps):
            label = rank_label(taxon, level, ranks)
            if label not in abs_map:
                abs_map[label] = [0] * n
            for i in range(n):
                abs_map[label][i] += values[i]
    return maps


def to_percent(label_counts, n, depth):
    out = OrderedDict()
    for label, values in label_counts.items():
        out[label] = [
            (values[i] / depth[i] if depth[i] else 0.0) for i in range(n)
        ]
    return out


def build(project):
    project = os.path.abspath(project)
    cfg = load_config(project)
    domain = cfg.get("KEEP_DOMAIN", "").strip()
    otu_path = os.path.join(project, "qiime2", "export", "otu_table.tsv")
    tax_path = os.path.join(project, "qiime2", "export", "taxonomy.tsv")
    fasta_path = os.path.join(project, "qiime2", "export", "dna-sequences.fasta")
    if not os.path.isfile(otu_path) or not os.path.isfile(tax_path):
        raise SystemExit("缺少 otu_table.tsv 或 taxonomy.tsv")
    samples, counts = read_otu(otu_path)
    taxa = read_taxonomy(tax_path)
    n = len(samples)

    tmpl_fa = cfg.get("TEMPLATE_ASV_REP", "").strip()
    tmpl_tax = cfg.get("TEMPLATE_ASV_TAXA", "").strip()
    if tmpl_fa and tmpl_tax:
        if not os.path.isfile(tmpl_fa) or not os.path.isfile(tmpl_tax):
            raise SystemExit(f"模板注释文件不存在: {tmpl_fa} / {tmpl_tax}")
        if not os.path.isfile(fasta_path):
            raise SystemExit("缺少 dna-sequences.fasta，无法按序列对齐模板注释")
        n_over = overlay_template_taxonomy(taxa, fasta_path, tmpl_fa, tmpl_tax)
        # 写出对齐后的注释，供建树与复查
        aligned_tax = os.path.join(project, "qiime2", "export", "taxonomy.aligned.tsv")
        write_taxonomy(aligned_tax, taxa, list(counts.keys()))
        print(f"按模板序列对齐注释 {n_over} 条 → {aligned_tax}", flush=True)

    if domain:
        counts = OrderedDict(
            (f, v)
            for f, v in counts.items()
            if keep_domain(taxa.get(f, "Unassigned"), domain)
        )
        if not counts:
            raise SystemExit(f"按 KEEP_DOMAIN={domain} 过滤后无 ASV")
        print(f"KEEP_DOMAIN={domain}，保留 {len(counts)} 条 ASV", flush=True)

    depth = sample_depths(counts, n)
    min_depth = min(depth) if depth else 0
    if min_depth <= 0:
        raise SystemExit("存在深度为 0 的样本，无法抽平")

    rarefied = rarefy_table(counts, n, min_depth)
    rare_depth = [min_depth] * n

    items = []
    for feature, values in counts.items():
        items.append((sum(values), feature, values, taxa.get(feature, "Unassigned")))
    items.sort(key=lambda item: (-item[0], item[1]))

    name_of = {}
    original_rows = []
    for index, (_total, feature, values, taxon) in enumerate(items, start=1):
        name = f"ASV{index}"
        name_of[feature] = name
        original_rows.append([name] + [str(v) for v in values] + [taxon])

    normalize_rows = []
    for _total, feature, _values, taxon in items:
        name = name_of[feature]
        rvals = rarefied.get(feature, [0] * n)
        if sum(rvals) == 0:
            continue
        normalize_rows.append([name] + [str(v) for v in rvals] + [taxon])

    root = os.path.join(project, "01.ASV.tax")
    write_table(os.path.join(root, "01.original", "asv_taxa_table.xls"), samples, original_rows)
    write_table(os.path.join(root, "02.normalize", "asv_taxa_table.xls"), samples, normalize_rows)

    ranks = ranks_for_taxa([t for *_rest, t in items])
    orig_maps = aggregate_ranks(items, n, ranks)
    rare_list = [
        (sum(rarefied[f]), f, rarefied[f], taxa.get(f, "Unassigned"))
        for f in rarefied
    ]
    rare_maps = aggregate_ranks(rare_list, n, ranks)

    for level, (rank_name, _prefix) in enumerate(ranks):
        write_rank(
            os.path.join(root, "01.original", f"{rank_name}.xls"),
            samples,
            orig_maps[level],
            False,
        )
        write_rank(
            os.path.join(root, "02.normalize", f"{rank_name}.xls"),
            samples,
            rare_maps[level],
            False,
        )
        pct = to_percent(rare_maps[level], n, rare_depth)
        write_rank(
            os.path.join(root, "02.normalize", "tax_summary_r", f"{rank_name}.percent.xls"),
            samples,
            pct,
            True,
        )
        pct0 = to_percent(orig_maps[level], n, depth)
        write_rank(
            os.path.join(root, "01.original", "tax_summary_r", f"{rank_name}.percent.xls"),
            samples,
            pct0,
            True,
        )

    write_folder_notes(
        root, domain, min_depth, len(items), len(normalize_rows), bool(tmpl_fa and tmpl_tax)
    )
    print(f"抽平深度={min_depth}，原始 ASV={len(items)}，抽平后 ASV={len(normalize_rows)}")
    print(root)


def write_folder_notes(root, domain, min_depth, n_orig, n_rare, used_template):
    domain_line = (
        f"- 已按 KEEP_DOMAIN={domain} 过滤，只保留该域 ASV。\n"
        if domain
        else "- 未做域过滤，保留全部注释结果。\n"
    )
    tmpl_line = (
        "- 与模板代表序列完全相同的 ASV，注释采用模板口径（用于组成图对齐）。\n"
        if used_template
        else ""
    )
    notes = {
        ("01.original", "1.原始丰度说明.txt"): (
            "本文件夹：01.original（原始丰度，不做抽平）\n\n"
            "文件说明\n"
            "- asv_taxa_table.xls：ASV×样本原始计数，末列为物种注释。\n"
            "- domain.xls～species.xls：各分类水平汇总后的原始计数。\n"
            "- tax_summary_r/*.percent.xls：由原始计数换算的相对丰度。\n\n"
            "数值为 DADA2 绝对计数。基础注释为 QIIME 2 VSEARCH + SILVA 138。\n"
            f"{domain_line}{tmpl_line}"
            f"本表共 {n_orig} 条 ASV。\n"
        ),
        ("02.normalize", "2.归一化丰度说明.txt"): (
            "本文件夹：02.normalize（抽平归一化，与模板口径一致）\n\n"
            "文件说明\n"
            "- asv_taxa_table.xls：ASV×样本抽平后整数计数，同一样本之和等于抽平深度。\n"
            "- domain.xls～species.xls：各分类水平汇总后的抽平计数。\n"
            "- tax_summary_r/*.percent.xls：由抽平计数换算的相对丰度（列和为 1）。\n\n"
            f"归一化方式：无放回随机抽平到最小样本深度 {min_depth}（随机种子 {RAREFY_SEED}）。\n"
            f"{domain_line}{tmpl_line}"
            f"抽平后保留 {n_rare} 条 ASV。\n"
        ),
        ("03.phylogeny_tree", "3.系统发育树说明.txt"): (
            "本文件夹：03.phylogeny_tree（系统发育树）\n\n"
            "由 build_phylogeny.py 生成：asv_rep.fasta、asv_rep_aligned.fasta、phylogeny.tre。\n"
            "ASV 编号与 01.original/asv_taxa_table.xls 一致。\n"
            "用于 Faith's PD、UniFrac 等依赖系统发育的分析。\n"
        ),
    }
    for (subdir, name), text in notes.items():
        folder = os.path.join(root, subdir)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, name), "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 taxa_abundance.py 项目目录")
    build(sys.argv[1])


if __name__ == "__main__":
    main()
